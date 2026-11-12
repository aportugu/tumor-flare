# Generic Firth/bootstrap and survival models on normalized controlled inputs.
tf_ash_covariates <- c("age10","ecog2","female","prior_lines","hrca","true_emd")
tf_ash_specs <- function() {
  data.frame(field=c(tf_ash_covariates,"hgb","albumin","log2_ldh","log2_ferritin","log2_crp"),
    label=c("Age, per 10 years","ECOG 2–5 vs 0–1","Female vs male","Prior treatment lines, per line",
      "At least one HRCA vs none","True EMD vs none/paraskeletal only","Hemoglobin, per g/dL",
      "Albumin, per g/dL","LDH, per doubling","Ferritin, per doubling","CRP, per doubling"),
    family=c(rep("Baseline",6),rep("Laboratory",5)),
    primary_question=c(rep(FALSE,5),TRUE,rep(FALSE,5)),stringsAsFactors=FALSE)
}

tf_ash_model_data <- function(cohort) {
  b <- tf_baseline_data(cohort); o <- tf_outcome_data(cohort)
  ecog <- suppressWarnings(as.numeric(b$ecog))
  d <- data.frame(course_id=cohort$course_id,patient_cluster_id=cohort$patient_cluster_id,
    site_key=cohort$site_key,flare=as.integer(cohort$analysis_flare),
    age10=ifelse(is.finite(b$age)&b$age>0,b$age/10,NA_real_),
    ecog2=ifelse(ecog %in% 0:5,as.integer(ecog>=2),NA_integer_),
    female=ifelse(b$sex=="Female",1L,ifelse(b$sex=="Male",0L,NA_integer_)),
    prior_lines=ifelse(is.finite(b$prior_lots)&b$prior_lots>=0,b$prior_lots,NA_real_),
    hrca=as.integer(b$hrca),true_emd=ifelse(b$emd=="True EMD with/without paraskeletal",1L,
      ifelse(b$emd %in% c("None","Paraskeletal only"),0L,NA_integer_)),
    hgb=b$hgb,albumin=b$albumin,stringsAsFactors=FALSE)
  for(field in c("ldh","ferritin","crp")) {
    z<-b[[field]]; d[[paste0("log2_",field)]]<-ifelse(is.finite(z)&z>0,log2(ifelse(z>0,z,NA_real_)),NA_real_)
  }
  known <- cohort$analysis_response %in% tf_response_levels[1:7]
  for(threshold in c("pr","vgpr","cr")) {
    levels <- switch(threshold,pr=tf_response_levels[1:4],vgpr=tf_response_levels[1:3],cr=tf_response_levels[1:2])
    d[[paste0("response_",threshold)]]<-ifelse(known,as.integer(cohort$analysis_response %in% levels),NA_integer_)
  }
  for(field in c("crs_1","icans_1","toci_treatment","toci_prophylaxis","any_admission","infection_30d","severe_infection_30d")) d[[field]]<-as.integer(o[[field]])
  d$hospital_days_admitted<-ifelse(!is.na(o$hospital_days)&o$hospital_days>0,o$hospital_days,NA_real_)
  if(anyDuplicated(d$course_id)) stop("ASH model input must have unique course IDs.")
  if(anyNA(d$patient_cluster_id)||any(!nzchar(d$patient_cluster_id))) stop("Documented cluster IDs are required.")
  d
}

tf_ash_missingness <- function(d,fields,scenario) {
  dplyr::bind_rows(lapply(fields,function(field) dplyr::bind_rows(lapply(c(0,1),function(group) {
    z<-d[d$flare==group,,drop=FALSE]
    data.frame(scenario=scenario,field=field,group=if(group)"Flare candidate" else "No included flare",
      total=nrow(z),evaluable=sum(is.finite(z[[field]])),unavailable=sum(!is.finite(z[[field]])))
  }))))
}

tf_ash_firth_fit <- function(d,response,predictors,term) {
  warnings<-character(); failure<-NULL
  if(!nrow(d)||length(unique(d[[response]]))<2L||length(unique(d[[term]]))<2L)
    return(list(coefficient=NA_real_,diagnostic="No outcome or predictor variation",converged=FALSE))
  formula<-stats::reformulate(predictors,response=response)
  fit<-tryCatch(withCallingHandlers(logistf::logistf(formula,data=d,pl=FALSE,
    control=logistf::logistf.control(maxit=200,lconv=1e-7,gconv=1e-7,xconv=1e-7)),
    warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),
    error=function(e){failure<<-conditionMessage(e);NULL})
  ok<-!is.null(fit)&&all(is.finite(fit$coefficients))&&length(fit$conv)==3L&&
    all(is.finite(fit$conv))&&all(abs(fit$conv)<=1.1e-7)&&
    !any(grepl("converg|singular|infinite|maximum number",warnings,ignore.case=TRUE))
  coef<-if(ok)unname(fit$coefficients[term]) else NA_real_
  list(coefficient=coef,converged=ok&&is.finite(coef),
    diagnostic=if(ok)"Firth penalized likelihood converged" else paste(c(failure,warnings,"Convergence criteria not met"),collapse="; "))
}

tf_ash_cluster_indices <- function(cluster_ids) {
  blocks<-split(seq_along(cluster_ids),cluster_ids)
  unlist(blocks[sample.int(length(blocks),length(blocks),replace=TRUE)],use.names=FALSE)
}

tf_ash_bootstrap_inference <- function(coefficient,draws,B,minimum=ceiling(.9*B)) {
  valid<-draws[is.finite(draws)]; se<-if(length(valid)>1L)stats::sd(valid) else NA_real_
  ok<-length(valid)>=minimum&&is.finite(coefficient)&&is.finite(se)&&se>0
  ci<-if(ok)exp(coefficient+c(-1,1)*stats::qnorm(.975)*se) else c(NA_real_,NA_real_)
  if(any(!is.finite(ci)))ok<-FALSE
  list(successes=length(valid),se=if(ok)se else NA_real_,estimate=if(is.finite(coefficient))exp(coefficient) else NA_real_,
    conf_low=if(ok)ci[1] else NA_real_,conf_high=if(ok)ci[2] else NA_real_,
    p_value=if(ok)2*stats::pnorm(-abs(coefficient/se)) else NA_real_,inference_available=ok)
}

tf_ash_firth_cluster <- function(d,response,predictors,term,label,family,scenario,
    model_id,B=1000L,seed=20261106L,cache_dir=NULL,progress=TRUE) {
  fields<-unique(c(response,predictors)); use<-stats::complete.cases(d[fields])&
    apply(d[fields],1,function(x)all(is.finite(x)))
  z<-d[use,,drop=FALSE]
  row<-data.frame(model_id=model_id,scenario=scenario,family=family,label=label,response=response,
    term=term,adjustment=paste(setdiff(predictors,term),collapse="; "),total_courses=nrow(d),
    n_courses=nrow(z),n_patient_clusters=length(unique(z$patient_cluster_id)),
    candidates=sum(z$flare),events=sum(z[[response]]),candidate_events=sum(z[[response]][z$flare==1]),
    excluded_missing=nrow(d)-nrow(z),estimate=NA_real_,conf_low=NA_real_,conf_high=NA_real_,p_value=NA_real_,
    bootstrap_successes=0L,bootstrap_requested=B,bootstrap_seed=seed,bootstrap_se=NA_real_,
    model_status="Not estimable",diagnostic="",stringsAsFactors=FALSE)
  # Cache is keyed to the exact coded analytic inputs, parameters and implementation.
  keyfile<-tempfile();on.exit(unlink(keyfile),add=TRUE)
  saveRDS(list(z=z[unique(c("course_id","patient_cluster_id","flare",fields))],response=response,
    predictors=predictors,term=term,B=B,seed=seed,
    implementation=lapply(list(tf_ash_firth_fit,tf_ash_cluster_indices,tf_ash_bootstrap_inference),function(f)deparse(body(f))),
    logistf_version=as.character(utils::packageVersion("logistf"))),keyfile)
  key<-unname(tools::md5sum(keyfile)); cache<-if(!is.null(cache_dir))file.path(cache_dir,paste0(model_id,"_",key,".rds")) else NULL
  if(!is.null(cache)&&file.exists(cache))return(readRDS(cache))
  point<-tf_ash_firth_fit(z,response,predictors,term);draws<-rep(NA_real_,B)
  failure_counts<-integer()
  if(point$converged&&length(unique(z$patient_cluster_id))>=2L) {
    set.seed(seed)
    for(k in seq_len(B)) {
      zz<-z[tf_ash_cluster_indices(z$patient_cluster_id),,drop=FALSE]
      fit<-tf_ash_firth_fit(zz,response,predictors,term)
      if(fit$converged)draws[k]<-fit$coefficient else {
        reason<-fit$diagnostic
        if(!reason %in% names(failure_counts))failure_counts[reason]<-0L
        failure_counts[reason]<-failure_counts[reason]+1L
      }
      if(progress&&k%%250L==0L)message(model_id,": ",k,"/",B," cluster resamples")
    }
  }
  inf<-tf_ash_bootstrap_inference(point$coefficient,draws,B)
  for(field in c("estimate","conf_low","conf_high","p_value"))row[[field]]<-inf[[field]]
  row$bootstrap_successes<-inf$successes;row$bootstrap_se<-inf$se
  row$model_status<-if(inf$inference_available)"Estimated (exploratory)" else if(point$converged)"Point estimate only; bootstrap inference suppressed" else "Not estimable"
  row$diagnostic<-paste(point$diagnostic,if(!inf$inference_available)paste0("; fewer than ",ceiling(.9*B)," converged resamples or invalid variance") else "; patient-cluster bootstrap Wald inference",sep="")
  diagnostics<-data.frame(model_id=model_id,scenario=scenario,family=family,requested=B,successful=inf$successes,
    seed=seed,minimum_successes=ceiling(.9*B),point_converged=point$converged,
    failure_summary=paste(names(failure_counts),failure_counts,sep=" [",collapse="; "),stringsAsFactors=FALSE)
  membership<-data.frame(model_id=model_id,scenario=scenario,course_id=d$course_id,included=use,
    exclusion_reason=vapply(seq_len(nrow(d)),function(i)if(use[i])"" else paste(fields[!is.finite(as.numeric(d[i,fields]))],collapse="; "),character(1)))
  out<-list(table=row,diagnostics=diagnostics,membership=membership,draws=draws)
  if(!is.null(cache)){dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE);saveRDS(out,cache)}
  out
}

tf_ash_fit_cox <- function(intervals,covariates,scenario,model_label,total_courses) {
  use<-if(length(covariates))stats::complete.cases(intervals[covariates]) else rep(TRUE,nrow(intervals))
  z<-intervals[use,,drop=FALSE]
  base<-edc_fit_exposure_cox(z,"flare_td",counting_process=TRUE)
  row<-base$table;row$scenario<-scenario;row$model<-model_label;row$total_model_eligible_courses<-total_courses
  row$excluded_missing<-total_courses-length(unique(z$course_id));row$adjustment<-paste(covariates,collapse="; ")
  row$ph_global_p<-NA_real_;row$ph_flare_p<-NA_real_
  if(!length(covariates))return(list(table=row,fit=base$fit))
  row[c("estimate","conf_low","conf_high","p_value")]<-NA_real_;row$model_status<-"Not estimable"
  if(!nrow(z)||row$exposed_events==0L||row$unexposed_events==0L)return(list(table=row,fit=NULL))
  warnings<-character();failure<-NULL
  formula<-stats::as.formula(paste("survival::Surv(tstart,tstop,event) ~ flare_td +",paste(covariates,collapse=" + ")))
  fit<-tryCatch(withCallingHandlers(survival::coxph(formula,data=z,cluster=patient_cluster_id,
    ties="efron",robust=TRUE,x=TRUE,model=TRUE),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),
    error=function(e){failure<<-conditionMessage(e);NULL})
  ok<-!is.null(fit)&&all(is.finite(stats::coef(fit)))&&!any(grepl("converg|infinite|singular|iteration",warnings,ignore.case=TRUE))
  if(ok) {
    beta<-unname(stats::coef(fit)["flare_td"]);se<-sqrt(stats::vcov(fit)["flare_td","flare_td"])
    ci<-exp(beta+c(-1,1)*stats::qnorm(.975)*se)
    ok<-all(is.finite(c(beta,se,ci)))&&se>0
    if(ok){row$estimate<-exp(beta);row$conf_low<-ci[1];row$conf_high<-ci[2];row$p_value<-2*stats::pnorm(-abs(beta/se));row$model_status<-"Estimated (exploratory)"}
    ph<-tryCatch(survival::cox.zph(fit),error=function(e)NULL)
    if(!is.null(ph)){row$ph_global_p<-ph$table["GLOBAL","p"];row$ph_flare_p<-ph$table["flare_td","p"]}
  }
  row$diagnostic<-if(ok)paste("Six baseline covariates; Efron ties; patient-cluster robust SE",paste(warnings,collapse="; "),sep="; ") else paste(c(failure,warnings,"Invalid adjusted fit"),collapse="; ")
  list(table=row,fit=if(ok)fit else NULL)
}

tf_survival_scenarios <- function(analysis,cohorts=NULL) {
  if(is.null(cohorts))cohorts<-analysis$cohorts
  historical<-analysis$cohorts$primary
  pos<-match(historical$course_id,analysis$original$course_id)
  for(field in c("pd_date","death_date","last_contact_date"))historical[[field]]<-analysis$original[[field]][pos]
  historical$pfs_excluded<-FALSE
  cohorts$original_master_endpoints<-historical
  rules<-analysis$survival_scenario_rules
  if(!is.null(rules)&&nrow(rules))for(scenario in unique(rules$scenario)) {
    d<-analysis$cohorts$primary
    rr<-rules[rules$scenario==scenario,,drop=FALSE]
    for(i in seq_len(nrow(rr))) {
      field<-rr$field[i];pos<-match(rr$course_id[i],d$course_id)
      if(is.na(pos)||!field %in% c("onset_date","pd_date","death_date","last_contact_date","pfs_excluded"))
        stop("Invalid controlled survival sensitivity rule")
      value<-if(rr$action[i]=="invalidate")NA else rr$value[i]
      if(inherits(d[[field]],"Date"))value<-as.Date(value)else if(is.logical(d[[field]]))value<-as.logical(value)
      d[[field]][pos]<-value
    }
    cohorts[[scenario]]<-d
  }
  cohorts
}

tf_ash_survival <- function(analysis,cohorts,model_data) {
  tables<-members<-landmarks<-list()
  cohorts<-tf_survival_scenarios(analysis,cohorts)
  model_data<-lapply(cohorts,tf_ash_model_data)
  for(scenario in names(cohorts))for(endpoint in c("PFS","OS")) {
    cohort<-cohorts[[scenario]];d<-model_data[[scenario]]
    exclude<-cohort$course_id[cohort$pfs_excluded]
    input<-cohort[c("course_id","patient_cluster_id","start_date","pd_date","death_date","last_contact_date","analysis_flare","onset_date")]
    run<-edc_time_dependent_cox(input,endpoint,exclude)
    intervals<-merge(run$intervals,d[c("course_id",tf_ash_covariates)],by="course_id",sort=FALSE)
    complete<-stats::complete.cases(d[tf_ash_covariates]);cc_ids<-d$course_id[complete]
    n<-length(unique(intervals$course_id))
    for(model in c("Unadjusted: all eligible","Unadjusted: adjusted subset","Adjusted: six baseline variables")) {
      z<-if(model=="Unadjusted: all eligible")intervals else intervals[intervals$course_id %in% cc_ids,,drop=FALSE]
      fit<-tf_ash_fit_cox(z,if(startsWith(model,"Adjusted"))tf_ash_covariates else character(),scenario,model,n)
      fit$table$endpoint<-endpoint
      fit$table$scenario_interpretation<-if(scenario=="original_master_endpoints")
        "Historical audit of superseded master values; not a valid clinical endpoint sensitivity or canonical inference"else
        "Exploratory clinical analysis using authoritative reconciled endpoints"
      tables[[length(tables)+1L]]<-fit$table
    }
    members[[length(members)+1L]]<-data.frame(scenario=scenario,endpoint=endpoint,course_id=d$course_id,
      endpoint_and_timing_eligible=d$course_id %in% intervals$course_id,
      complete_baseline_covariates=complete,
      adjusted_included=d$course_id %in% intervals$course_id&complete)
    lm<-edc_landmark_analysis(input,endpoint,pfs_exclude_course_ids=exclude)
    lm$table$scenario<-scenario
    lm$table$scenario_interpretation<-if(scenario=="original_master_endpoints")
      "Historical audit of superseded master values; not a valid clinical endpoint sensitivity or canonical inference"else
      "Exploratory clinical analysis using authoritative reconciled endpoints"
    landmarks[[length(landmarks)+1L]]<-lm$table
  }
  list(table=dplyr::bind_rows(tables),membership=dplyr::bind_rows(members),landmark_table=dplyr::bind_rows(landmarks))
}

tf_build_ash_extension <- function(analysis,output_dir,B=1000L,seed=20261106L,progress=TRUE) {
  if(!requireNamespace("logistf",quietly=TRUE))stop("The requested logistf package is unavailable in the project R library.")
  cohorts<-list(primary=analysis$cohorts$primary,strict=analysis$cohorts$strict,
    without_dominant_site=analysis$cohorts$primary[!analysis$cohorts$primary$site_key %in% analysis$excluded_site_keys,,drop=FALSE])
  data<-lapply(cohorts,tf_ash_model_data); specs<-tf_ash_specs()
  models<-list();miss<-list();clinical<-list();los<-list()
  cache_dir<-file.path(output_dir,"build","model_cache")
  add_model<-function(d,response,predictors,term,label,family,scenario,id) {
    models[[length(models)+1L]]<<-tf_ash_firth_cluster(d,response,predictors,term,label,family,scenario,id,
      B=B,seed=analysis$model_seeds$seed[match(id,analysis$model_seeds$model_id)],cache_dir=cache_dir,progress=progress)
  }
  for(scenario in names(cohorts)) {
    d<-data[[scenario]]
    miss[[length(miss)+1L]]<-tf_ash_missingness(d,specs$field,scenario)
    for(k in seq_len(nrow(specs))) {
      sp<-specs[k,];add_model(d,"flare",sp$field,sp$field,sp$label,sp$family,scenario,paste(scenario,sp$field,sep="_"))
    }
    for(threshold in c("pr","vgpr","cr")) {
      label<-switch(threshold,pr="Response ≥PR",vgpr="Response ≥VGPR",cr="Response CR/sCR")
      response<-paste0("response_",threshold)
      miss[[length(miss)+1L]]<-tf_ash_missingness(d,response,scenario)
      for(adjusted in c(FALSE,TRUE))add_model(d,response,c("flare",if(adjusted)tf_ash_covariates),"flare",
        paste(label,if(adjusted)"(adjusted)" else "(unadjusted)"),"Response",scenario,
        paste(scenario,response,if(adjusted)"adjusted" else "unadjusted",sep="_"))
    }
    fields<-c("crs_1","icans_1","toci_treatment","toci_prophylaxis","any_admission","infection_30d","severe_infection_30d")
    labels<-c("Any CRS","Any ICANS","Treatment tocilizumab","Prophylactic tocilizumab","Any admission","Infection within 30 days","Severe infection within 30 days")
    # These are post-treatment course-level comparisons, never baseline predictors.
    for(k in seq_along(fields)) {
      field<-fields[k];z<-d[!is.na(d[[field]]),,drop=FALSE]
      tab<-table(factor(z$flare,levels=0:1),factor(z[[field]],levels=0:1))
      test<-if(nrow(z)&&all(rowSums(tab)>0))stats::fisher.test(tab) else NULL
      clinical[[length(clinical)+1L]]<-data.frame(scenario=scenario,field=field,label=labels[k],
        flare_yes=tab[2,2],flare_evaluable=sum(tab[2,]),flare_total=sum(d$flare),
        comparison_yes=tab[1,2],comparison_evaluable=sum(tab[1,]),comparison_total=sum(d$flare==0),
        fisher_p=if(!is.null(test))test$p.value else NA_real_,
        interpretation="Exploratory post-treatment association; Fisher test assumes independent courses")
    }
    z<-d[!is.na(d$hospital_days_admitted),,drop=FALSE]
    # Match the abstract's admitted-course comparison. This course-level rank test
    # is exploratory and does not provide patient-cluster-adjusted inference.
    los_test<-if(all(c(0,1) %in% z$flare))tryCatch(stats::wilcox.test(
      hospital_days_admitted ~ flare,data=z,exact=FALSE),error=function(e)NULL) else NULL
    for(group in c(0,1)) {
      days<-z$hospital_days_admitted[z$flare==group]
      los[[length(los)+1L]]<-data.frame(scenario=scenario,group=if(group)"Flare candidate" else "No included flare",
        admitted_evaluable=length(days),total_courses=sum(d$flare==group),
        median=if(length(days))median(days) else NA_real_,q1=if(length(days))unname(quantile(days,.25)) else NA_real_,
        q3=if(length(days))unname(quantile(days,.75)) else NA_real_,
        wilcoxon_p=if(!is.null(los_test))los_test$p.value else NA_real_,
        test_status=if(!is.null(los_test))"Exploratory; assumes independent courses" else "Not estimable: both groups required",
        definition="Numeric total hospital days among recorded admissions; zero-day outpatient courses excluded; no flare-attributable duration inferred")
    }
  }
  table<-dplyr::bind_rows(lapply(models,`[[`,"table"))
  table$sparse_note<-vapply(seq_len(nrow(table)),function(i) {
    r<-table[i,]
    if(r$response=="flare"&&r$candidates<5L)return(paste("Only",r$candidates,"flare candidates; sparse exploratory estimate"))
    if(r$family=="Response") {
      counts<-c(r$candidate_events,r$candidates-r$candidate_events)
      if(any(counts<5L))return(paste0("Sparse candidate response cells: ",counts[1]," above and ",counts[2]," below threshold"))
    }
    ""
  },character(1))
  table$lab_bh_p<-NA_real_
  for(scenario in names(cohorts)) {
    hit<-table$scenario==scenario&table$family=="Laboratory"
    table$lab_bh_p[hit]<-stats::p.adjust(table$p_value[hit],method="BH",n=5L)
  }
  survival<-tf_ash_survival(analysis,cohorts,data)
  result<-list(models=table,bootstrap_diagnostics=dplyr::bind_rows(lapply(models,`[[`,"diagnostics")),
    model_membership=dplyr::bind_rows(lapply(models,`[[`,"membership")),missingness=dplyr::bind_rows(miss),
    survival=survival$table,survival_membership=survival$membership,landmark_sensitivity=survival$landmark_table,
    clinical=dplyr::bind_rows(clinical),hospitalization=dplyr::bind_rows(los),
    specs=specs,parameters=list(bootstrap_replicates=B,minimum_successes=ceiling(.9*B),seed=seed,
      logistf_version=as.character(utils::packageVersion("logistf")),adjustment_variables=tf_ash_covariates))
  for(name in names(result)[vapply(result,is.data.frame,logical(1))])utils::write.csv(result[[name]],
    file.path(output_dir,paste0("ash_",name,".csv")),row.names=FALSE,na="")
  jsonlite::write_json(result$parameters,file.path(output_dir,"ash_model_parameters.json"),auto_unbox=TRUE,pretty=TRUE)
  result
}
