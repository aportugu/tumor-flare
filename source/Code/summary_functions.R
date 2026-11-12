# Aggregate manuscript display exports and restricted row-level analytical audit.
tf_write <- function(d,output_dir,name)utils::write.csv(d,file.path(output_dir,paste0(name,".csv")),row.names=FALSE,na="")

tf_counts <- function(cases,fields,labels,sections,groups=TRUE) {
  populations<-list(All=cases)
  if(groups) {
    populations[["Definite/Probable"]]<-cases[cases$investigator_confidence %in% c("Definite","Probable"),,drop=FALSE]
    populations[["Possible"]]<-cases[cases$investigator_confidence %in% "Possible",,drop=FALSE]
  }
  dplyr::bind_rows(lapply(names(populations),function(group)dplyr::bind_rows(lapply(seq_along(fields),function(i) {
    z<-populations[[group]];values<-tf_category(z[[fields[i]]]);t<-table(values)
    data.frame(group=group,section=sections[i],panel=sections[i],field=fields[i],label=labels[i],
      category=names(t),n=as.integer(t),denominator=nrow(z),percent=100*as.integer(t)/nrow(z),stringsAsFactors=FALSE)
  }))))
}

tf_export_descriptive <- function(a,out) {
  x<-a$master;c<-a$cases;p<-a$cohorts$primary;s<-a$cohorts$strict
  tf_write(x,out,"canonical_courses_restricted")
  tf_write(c,out,"canonical_cases_restricted")
  tf_write(a$provenance,out,"field_provenance_restricted")
  tf_write(a$register,out,"reconciliation_register_restricted")
  flow<-data.frame(stage=c("Master courses","Supported agents","First eligible BCMA / all GPRC5D courses",
    "Nonbridging eligible courses","Excluded: flare status unavailable","Primary known-status courses",
    "Primary suspected flare","Primary comparison","Definite/Probable-only courses","Definite/Probable-only flare",
    "Clinical returned series","Clinical bridging cases","Clinical later BCMA cases"),
    courses=c(nrow(x),sum(x$product_key %in% tf_supported_products),sum(x$product_key %in% tf_supported_products&x$keep_first_bcma),
    nrow(a$eligible),sum(is.na(a$eligible$analysis_flare)),nrow(p),sum(p$analysis_flare),sum(!p$analysis_flare),nrow(s),sum(s$analysis_flare),
    nrow(c),sum(c$bridging_any,na.rm=TRUE),sum(!c$keep_first_bcma&!c$bridging_any,na.rm=TRUE)))
  tf_write(flow,out,"flow")
  coverage<-dplyr::bind_rows(lapply(unique(x$site),function(site) {
    z<-x[x$site==site,,drop=FALSE];e<-a$eligible[a$eligible$site==site,,drop=FALSE];d<-c[c$site==site,,drop=FALSE]
    data.frame(site=site,master_courses=nrow(z),eligible_courses=nrow(e),unknown_flare=sum(is.na(e$analysis_flare)),
      known_flare_status=sum(!is.na(e$analysis_flare)),original_flags=sum(z$original_flare,na.rm=TRUE),
      returned_records=nrow(d),returned_candidates=sum(d$analysis_flare,na.rm=TRUE))
  }))
  tf_write(coverage,out,"coverage")
  incidence<-dplyr::bind_rows(lapply(names(a$cohorts),function(scenario)tf_incidence(a$cohorts[[scenario]],scenario)))
  incidence<-incidence[incidence$product=="Overall",,drop=FALSE];tf_write(incidence,out,"incidence")
  specs<-c(product="Agent",age="Age (years)",sex="Sex",race="Race",ecog="ECOG",emd="Baseline disease location",
    prior_lots="Prior treatment lines",prior_asct="Prior ASCT",prior_cart="Prior CAR-T",prior_bcma="Prior BCMA therapy",
    prior_gprc5d="Prior GPRC5D therapy",high_marrow="High marrow burden",oligosecretory="Oligo/nonsecretory",
    hrca="At least one HRCA",triple_refractory="Triple refractory",penta_refractory="Penta refractory",
    hgb="Hemoglobin (g/dL)",ldh="LDH (U/L)",albumin="Albumin (g/dL)",igg="IgG (mg/dL)",ferritin="Ferritin (ng/mL)",crp="CRP (mg/L)")
  baseline<-dplyr::bind_rows(lapply(names(a$cohorts),function(scenario)dplyr::mutate(tf_summary_fields(a$cohorts[[scenario]],specs),scenario=scenario)))
  # Categorical known selections use the variable's available denominator.
  # Missing selections show their fraction of the complete group separately.
  baseline$evaluable<-baseline$n
  for(scenario in names(a$cohorts))for(measure in unname(specs)) {
    field<-names(specs)[match(measure,specs)]
    if(is.numeric(a$cohorts[[scenario]][[field]]))next
    all_categories<-unique(baseline$category[baseline$scenario==scenario&baseline$measure==measure])
    for(group in c("All","Flare candidate","No included flare")) {
      d<-a$cohorts[[scenario]]
      if(group!="All")d<-d[d$analysis_flare==(group=="Flare candidate"),,drop=FALSE]
      values<-tf_category(d[[field]])
      available<-sum(!values %in% c("Blank","Unknown","Not applicable","NE/Unknown","Not reassessed"))
      rows<-which(baseline$scenario==scenario&baseline$measure==measure&baseline$group==group)
      absent<-setdiff(all_categories,baseline$category[rows])
      if(length(absent)) {
        extra<-baseline[rep(rows[1],length(absent)),,drop=FALSE]
        extra$category<-absent;extra$n<-0L;extra$summary<-NA_character_
        baseline<-rbind(baseline,extra)
        rows<-which(baseline$scenario==scenario&baseline$measure==measure&baseline$group==group)
      }
      baseline$evaluable[rows]<-available
      missing<-baseline$category[rows] %in% c("Blank","Unknown","Not applicable","NE/Unknown","Not reassessed")
      denominators<-ifelse(missing,nrow(d),available)
      baseline$percent[rows]<-ifelse(denominators>0,100*baseline$n[rows]/denominators,NA_real_)
    }
  }
  tf_write(baseline,out,"table1_baseline")
  fields<-c("management_steroids","management_radiation","management_surgery","management_toci_anakinra",
    "management_nonopioid","management_opioids","bsab_action","resolution_status","local_lesion_outcome",
    "first_postflare_response","subsequent_response","best_index_response","improved_while_bsab_continued")
  labels<-c("Corticosteroids","Radiation","Surgery/procedure","Tocilizumab/anakinra","Non-opioid analgesia","Opioid analgesia",
    "BsAb action","Symptom-resolution status","Local-lesion outcome","First postflare systemic assessment",
    "Later supported systemic response","Best documented index-BsAb response","Improved without permanent BsAb discontinuation")
  sections<-c(rep("Interventions",6),"BsAb treatment",rep("Subsequent course",6))
  tf_write(tf_counts(c,fields,labels,sections),out,"table2_clinical_course")
  fields<-c("pain_flare","palpable_enlargement","imaging_worsening","new_lesion","neurologic_cns","biochemical_or_organ",
    "organ_dysfunction","flare_hospitalization","neurologic_deficit","pathologic_fracture","cord_compression","icu_due_to_flare")
  labels<-c("Pain","Palpable enlargement","Imaging worsening","New lesion","Neurological manifestation","Biochemical/organ manifestation",
    "Organ dysfunction","Primary flare admission","Neurological deficit","Pathological fracture","Cord compression","Primary flare ICU")
  tf_write(tf_counts(c,fields,labels,c(rep("Presentation",6),rep("Clinical impact",6)),FALSE),out,"presentation")
  fields<-c("focality","known_or_new_site","skeletal_flare_type","true_emd_classification","paraskeletal_involvement","spinal_involvement","cns_involvement")
  labels<-c("Extent","Known/new disease site","Skeletal pattern","True EMD","Paraskeletal involvement","Spinal involvement","CNS involvement")
  tf_write(tf_counts(c,fields,labels,rep("Anatomy",length(fields)),FALSE),out,"anatomy")
  anatomy_plot<-c
  anatomy_plot$true_emd_any<-ifelse(c$true_emd_classification %in% c("Visceral EMD (organ involvement)","Non-visceral EMD (no organ involvement)","Both visceral and non-visceral EMD"),"Yes",ifelse(c$true_emd_classification %in% "No true EMD","No","Unknown"))
  tf_write(tf_counts(anatomy_plot,c("paraskeletal_involvement","true_emd_any","spinal_involvement","cns_involvement"),
    c("Paraskeletal involvement","True EMD","Spinal involvement","CNS involvement"),rep("Anatomy",4),FALSE),out,"anatomy_binary")
  start<-x$start_date[match(c$course_id,x$course_id)]
  timeline<-data.frame(case_id=c$case_id,course_id=c$course_id,product=c$bsab_product,confidence=c$investigator_confidence,
    bridging=c$bridging_any,comparative_included=c$comparative_included,
    onset_day=as.numeric(c$flare_onset_date-start),peak_day=as.numeric(c$flare_peak_date-start),
    resolution_day=as.numeric(c$flare_resolution_date-start),onset_lower_day=as.numeric(c$onset_lower_date-start),
    onset_upper_day=as.numeric(c$onset_upper_date-start),onset_estimated=c$onset_estimated,peak_estimated=c$peak_estimated,
    resolution_estimated=c$resolution_estimated,symptom_status=c$resolution_status,local_outcome=c$local_lesion_outcome)
  timeline$case_label<-timeline$case_id
  timeline$onset_days<-timeline$onset_day;timeline$peak_days<-timeline$peak_day;timeline$resolution_days<-timeline$resolution_day
  timeline$onset_lower<-timeline$onset_lower_day;timeline$onset_upper<-timeline$onset_upper_day
  tf_write(timeline,out,"case_timeline")
  days<-timeline$onset_day[!is.na(timeline$onset_day)]
  timing<-data.frame(metric=c("Clinical series","Point onset available","Recorded peak available","Resolution date available","Median onset days","Minimum onset days","Maximum onset days","Onset by day 7","Onset by day 14","Onset by day 30"),
    value=c(nrow(c),length(days),sum(!is.na(timeline$peak_day)),sum(!is.na(timeline$resolution_day)),median(days),min(days),max(days),sum(days<=7),sum(days<=14),sum(days<=30)),
    denominator=c(nrow(c),nrow(c),nrow(c),nrow(c),length(days),length(days),length(days),length(days),length(days),length(days)))
  tf_write(timing,out,"timing_summary")
  summary<-timeline[,setdiff(names(timeline),"course_id"),drop=FALSE]
  for(field in c("pain_flare","palpable_enlargement","imaging_worsening","new_lesion","neurologic_cns","biochemical_or_organ",
    "focality","known_or_new_site","skeletal_flare_type","paraskeletal_involvement","true_emd_classification","cns_involvement",
    "spinal_involvement","organ_dysfunction","flare_hospitalization","icu_due_to_flare","pathologic_fracture","cord_compression",
    "overlap_crs","overlap_icans","overlap_infection","alternative_fracture","alternative_toxicity","other_competing_explanation",
    "bsab_action","management_steroids","management_radiation","management_surgery","management_toci_anakinra",
    "first_postflare_response","subsequent_response","best_index_response","true_pd_within_14d","alternative_diagnosis"))summary[[field]]<-c[[field]]
  for(field in intersect(c("first_response_on_index","subsequent_response_on_index"),names(c)))summary[[field]]<-c[[field]]
  tf_write(summary,out,"case_summary")
  changes<-a$register[!startsWith(a$register$field,"case."),,drop=FALSE]
  change_summary<-as.data.frame(table(changes$field));names(change_summary)<-c("field","change_events")
  tf_write(change_summary,out,"reconciliation_summary")
}

tf_export_models <- function(a,m,out) {
  mapping<-c(models="models",bootstrap_diagnostics="bootstrap_diagnostics",missingness="missingness",
    survival="survival_models",clinical="clinical_comparisons",hospitalization="hospitalization",landmark_sensitivity="landmark_models")
  for(field in names(mapping))tf_write(m[[field]],out,mapping[[field]])
  tf_write(m$model_membership,out,"model_membership_restricted")
  tf_write(m$survival_membership,out,"survival_membership_restricted")
  thresholds<-dplyr::bind_rows(lapply(names(a$cohorts),function(scenario) {
    d<-a$cohorts[[scenario]]
    dplyr::bind_rows(lapply(c(FALSE,TRUE),function(group) {
      z<-d[d$analysis_flare==group,,drop=FALSE];known<-z$analysis_response %in% tf_response_levels[1:7]
      dplyr::bind_rows(lapply(c("pr","vgpr","cr"),function(threshold) {
        levels<-switch(threshold,pr=tf_response_levels[1:4],vgpr=tf_response_levels[1:3],cr=tf_response_levels[1:2])
        n<-sum(known);responders<-sum(z$analysis_response %in% levels);ci<-if(n)binom.test(responders,n)$conf.int else c(NA,NA)
        data.frame(scenario=scenario,group=if(group)"Flare candidate"else"No included flare",threshold=threshold,responders=responders,
          evaluable=n,total_courses=nrow(z),unavailable=nrow(z)-n,percent=if(n)100*responders/n else NA,
          lower_95_percent=100*ci[1],upper_95_percent=100*ci[2])
      }))
    }))
  }))
  tf_write(thresholds,out,"response_thresholds")
  distribution<-dplyr::bind_rows(lapply(names(a$cohorts),function(scenario) {
    d<-a$cohorts[[scenario]]
    dplyr::bind_rows(lapply(c(FALSE,TRUE),function(group) {
      z<-d[d$analysis_flare==group,,drop=FALSE];known<-z$analysis_response %in% tf_response_levels[1:7]
      categories<-c(tf_response_levels[1:7],"Unavailable")
      n<-vapply(categories,function(category)if(category=="Unavailable")sum(!known)else sum(z$analysis_response==category,na.rm=TRUE),integer(1))
      evaluable<-sum(known)
      data.frame(scenario=scenario,group=if(group)"Flare candidate"else"No included flare",category=categories,n=n,
        evaluable=evaluable,total_courses=nrow(z),unavailable=sum(!known),
        percent=ifelse(categories=="Unavailable",100*n/nrow(z),if(evaluable)100*n/evaluable else NA_real_))
    }))
  }))
  tf_write(distribution,out,"response_distribution")
  km<-risk<-landmark_rows<-endpoint_audit<-list()
  cohorts<-tf_survival_scenarios(a,c(a$cohorts,list(without_dominant_site=
    a$cohorts$primary[!a$cohorts$primary$site_key %in% a$excluded_site_keys,,drop=FALSE])))
  for(scenario in names(cohorts))for(endpoint in c("PFS","OS")) {
    d<-cohorts[[scenario]][c("course_id","patient_cluster_id","start_date","pd_date","death_date","last_contact_date","analysis_flare","onset_date")]
    exclusions<-cohorts[[scenario]]$course_id[cohorts[[scenario]]$pfs_excluded]
    run<-edc_landmark_analysis(d,endpoint,pfs_exclude_course_ids=exclusions)
    run$data$scenario<-scenario;landmark_rows[[length(landmark_rows)+1L]]<-run$data
    er<-edc_build_endpoint(d,endpoint,pfs_exclude_course_ids=exclusions)$audit
    er$scenario<-scenario;endpoint_audit[[length(endpoint_audit)+1L]]<-er
    if(scenario=="primary")for(day in c(30,60))for(group in c(0,1)) {
      z<-run$data[run$data$landmark_day==day&run$data$flare_by_landmark==group,,drop=FALSE]
      if(!nrow(z))next
      fit<-survival::survfit(survival::Surv(time_after_landmark_months,status)~1,data=z,conf.type="log-log")
      km[[length(km)+1L]]<-data.frame(endpoint=endpoint,landmark_day=day,group=group,
        time=c(0,fit$time),survival=c(1,fit$surv),lower=c(1,fit$lower),upper=c(1,fit$upper))
      times<-c(0,6,12,24,36,48)
      risk[[length(risk)+1L]]<-data.frame(endpoint=endpoint,landmark_day=day,group=group,time=times,
        n_risk=vapply(times,function(t)sum(z$time_after_landmark_months>=t),integer(1)))
    }
  }
  tf_write(dplyr::bind_rows(km),out,"landmark_km_steps")
  tf_write(dplyr::bind_rows(risk),out,"landmark_risk_counts")
  tf_write(dplyr::bind_rows(landmark_rows),out,"landmark_data_restricted")
  tf_write(dplyr::bind_rows(endpoint_audit),out,"endpoint_audit_restricted")
  final<-a$provenance[a$provenance$is_final_record,,drop=FALSE]
  endpoint_provenance<-list()
  for(endpoint in c("PFS","OS")) {
    d<-edc_build_endpoint(a$master[c("course_id","patient_cluster_id","start_date","pd_date","death_date","last_contact_date","analysis_flare","onset_date")],
      endpoint,pfs_exclude_course_ids=a$master$course_id[a$master$pfs_excluded])$audit
    inputs<-c("start_date",if(endpoint=="PFS")"pd_date","death_date","last_contact_date",if(endpoint=="PFS")"pfs_excluded")
    by_course<-split(final[final$field %in% inputs,,drop=FALSE],final$course_id[final$field %in% inputs])
    refs<-vapply(d$course_id,function(id) {
      r<-by_course[[id]]
      if(is.null(r))return("")
      paste(paste0(r$field,": ",r$source_file," [",r$source_sheet,", row ",r$source_row,"]",
        ifelse(nzchar(r$clarification_pointer),paste0("; ",r$clarification_pointer),""),
        ifelse(nzchar(r$controlled_rule_pointer),paste0("; ",r$controlled_rule_pointer),"")),collapse=" | ")
    },character(1))
    for(field in c("event_date","observed_date","status","followup_days","endpoint_eligible","endpoint_exclusion_reason")) {
      values<-d[[field]]
      if(field %in% c("observed_date","status","followup_days"))values[!d$endpoint_eligible]<-NA
      endpoint_provenance[[length(endpoint_provenance)+1L]]<-data.frame(course_id=d$course_id,field=paste(endpoint,field,sep="."),
        source="Derived from canonical endpoint fields",source_field=paste(inputs,collapse="; "),source_value=tf_value_text(values),
        source_state=ifelse(!d$endpoint_eligible,"Analytically excluded; no inferred censoring",tf_source_state(values)),
        original_master_value="",canonical_value=tf_value_text(values),changed=FALSE,
        reason=ifelse(d$endpoint_eligible,"Canonical event/censor derivation after authoritative source reconciliation",d$endpoint_exclusion_reason),
        source_file="",source_sheet="",source_row="",clarification_pointer="",controlled_rule_pointer="",supporting_references=refs)
    }
  }
  records<-dplyr::bind_rows(a$provenance,dplyr::bind_rows(endpoint_provenance))
  records$record_order<-seq_len(nrow(records));records$is_final_record<-!duplicated(paste(records$course_id,records$field),fromLast=TRUE)
  tf_write(records,out,"field_provenance_restricted")
  attribution_path<-file.path(a$project_root,"response_attribution_review_restricted.csv")
  if(file.exists(attribution_path))invisible(file.copy(attribution_path,file.path(out,"response_attribution_review_restricted.csv"),overwrite=TRUE))
}

tf_write_analysis_summary <- function(a,m,out) {
  p<-a$cohorts$primary;c<-a$cases
  inc<-tf_incidence(p,"primary");inc<-inc[inc$product=="Overall",,drop=FALSE]
  response<-utils::read.csv(file.path(out,"response_thresholds.csv"),stringsAsFactors=FALSE)
  old<-a$original$analysis_response[match(c$course_id,a$original$course_id)]
  new<-c$best_index_response
  survival<-m$survival[m$survival$scenario=="primary"&startsWith(m$survival$model,"Adjusted"),,drop=FALSE]
  emd<-m$models[m$models$scenario=="primary"&m$models$model_id=="primary_true_emd",,drop=FALSE]
  result<-list(counts=list(master=nrow(a$master),eligible=nrow(a$eligible),unknown_flare=sum(is.na(a$eligible$analysis_flare)),
    comparative=nrow(p),flare=sum(p$analysis_flare),comparison=sum(!p$analysis_flare),clinical=nrow(c),sites=length(unique(c$site)),
    strict_courses=nrow(a$cohorts$strict),strict_flare=sum(a$cohorts$strict$analysis_flare)),
    confidence=as.list(table(c$investigator_confidence)),products=as.list(table(c$bsab_product)),incidence=inc,
    timing=utils::read.csv(file.path(out,"timing_summary.csv")),response_thresholds=response[response$scenario=="primary",],
    response_changes=list(all_returned=sum(tf_value_text(old)!=tf_value_text(new)),
      comparative=sum(tf_value_text(old[c$comparative_included])!=tf_value_text(new[c$comparative_included])),
      canonical_unavailable=sum(p$analysis_flare&!p$analysis_response %in% tf_response_levels[1:7])),
    emd=emd,adjusted_survival=survival,bootstrap=list(models=nrow(m$models),requested=m$parameters$bootstrap_replicates,
      minimum_successes=m$parameters$minimum_successes,success_range=range(m$bootstrap_diagnostics$successful)),
    pfs_excluded_returned=sum(a$master$pfs_excluded&a$master$edc_received),
    pfs_excluded_comparative=sum(p$pfs_excluded),source_policy="Returned EDC/site clarification is authoritative; master supplies separate uncollected, noncontradicted concepts",
    survival_sensitivities=setdiff(unique(m$survival$scenario),"original_master_endpoints"),
    historical_endpoint_audit="original_master_endpoints is an audit of superseded values; not a clinically valid sensitivity or canonical inference",
    response_attribution=list(reviewed_courses=nrow(c),reviewed_assessments=2L*nrow(c),
      policy="Only accepted index-course categories enter response thresholds. Local-only, unconfirmed, post-salvage and post-CAR-T observations do not supply a category."))
  jsonlite::write_json(result,file.path(out,"analysis_summary.json"),auto_unbox=TRUE,dataframe="rows",pretty=TRUE,na="null")
}
