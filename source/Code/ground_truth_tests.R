# Generic regression scenarios for precedence, scope, endpoints, and resampling.
tf_ground_truth_tests <- function(master,cases,rules,analysis) {
  results<-list()
  check<-function(name,value)results[[length(results)+1L]]<<-data.frame(test=name,status=if(isTRUE(value))"PASS"else"FAIL")
  x<-analysis$master;c<-analysis$cases;p<-analysis$cohorts$primary;s<-analysis$cohorts$strict
  check("Every authoritative return has a canonical record",nrow(c)==nrow(cases)&&setequal(c$course_id,cases$course_id))
  check("Unique normalized course and case identities",!anyDuplicated(x$course_id)&&!anyDuplicated(c$case_id))
  check("Primary classified eligibility",all(p$eligibility_verified&p$keep_first_bcma&!p$bridging_any)&&!anyNA(p$analysis_flare))
  check("Strict excludes Possible without making controls",!any(s$investigator_confidence %in% "Possible")&&sum(!s$analysis_flare)==sum(!p$analysis_flare))
  j<-match(cases$course_id,x$course_id)
  check("Returned canonical start uses EDC including missing",identical(x$start_date[j],cases$step_up_day1))
  check("Returned product uses EDC",identical(x$product[j],as.character(cases$bsab_product)))
  expected<-vapply(seq_len(nrow(cases)),function(i)tf_best_index_response(cases$first_postflare_response[i],cases$subsequent_response[i],cases$first_response_on_index[i],cases$subsequent_response_on_index[i]),character(1))
  check("All returned responses use site assessments",identical(x$analysis_response[j],expected))
  check("Original master preserved",identical(master$analysis_response,analysis$original$analysis_response[match(master$course_id,analysis$original$course_id)]))
  check("Final provenance has one authoritative record per field",!anyDuplicated(paste(analysis$provenance$course_id[analysis$provenance$is_final_record],
    analysis$provenance$field[analysis$provenance$is_final_record])))
  check("Final first-course provenance reflects recomputed eligibility",all(x$course_id %in%
    analysis$provenance$course_id[analysis$provenance$field=="keep_first_bcma"&analysis$provenance$is_final_record&
    analysis$provenance$source=="Derived from canonical source fields"]))
  check("Final infection-window provenance records canonical derivation",all(x$course_id %in%
    analysis$provenance$course_id[analysis$provenance$field=="infection_30d"&analysis$provenance$is_final_record&
    analysis$provenance$source=="Derived from canonical source fields"]))
  if(all(c("source_file","source_sheet","source_row") %in% names(cases))) {
    r<-analysis$provenance[analysis$provenance$source=="EDC"&analysis$provenance$field=="analysis_response",,drop=FALSE]
    q<-match(r$course_id,cases$course_id)
    check("EDC response provenance preserves file sheet and row",identical(r$source_file,tf_value_text(cases$source_file[q]))&&
      identical(r$source_sheet,tf_value_text(cases$source_sheet[q]))&&identical(r$source_row,tf_value_text(cases$source_row[q])))
  }
  check("Known response does not upgrade from master",tf_best_index_response("PR","Unknown")=="PR")
  check("Later formal response supersedes first apparent PD",tf_best_index_response("PD","PR")=="PR")
  check("Salvage response excluded",tf_best_index_response("MR","CR",TRUE,FALSE)=="MR")
  check("Unknown response stays unavailable",tf_best_index_response("Unknown","Not reassessed")=="NE/Unknown")
  fixture<-cases[1,,drop=FALSE];mm<-master
  k<-match(fixture$course_id,mm$course_id)
  mm$analysis_response[k]<-"CR";fixture$first_postflare_response<-"Unknown";fixture$subsequent_response<-NA_character_
  z<-tf_reconcile(mm,fixture,trace=FALSE)
  check("Unknown/blank EDC never backfilled from master response",z$master$analysis_response[k]=="NE/Unknown")
  fixture$step_up_day1<-as.Date(NA_character_);z<-tf_reconcile(mm,fixture,trace=FALSE)
  check("Missing EDC start never backfilled",is.na(z$master$start_date[k]))
  fixture<-cases[1,,drop=FALSE];fixture$investigator_confidence<-NA_character_;z<-tf_reconcile(mm,fixture,trace=FALSE)
  check("Missing returned classification never inherits original flag",is.na(z$master$analysis_flare[k])&&!fixture$course_id %in% z$cohorts$primary$course_id)
  fixture$investigator_confidence<-"Not tumor flare";mm$original_flare[k]<-TRUE;mm$analysis_flare[k]<-TRUE;z<-tf_reconcile(mm,fixture,trace=FALSE)
  check("Site negative classification overrides master positive",identical(z$master$analysis_flare[k],FALSE))
  fixture<-cases[1,,drop=FALSE];fixture$flare_hospitalization<-"Yes";mm$any_admission[k]<-FALSE;z<-tf_reconcile(mm,fixture,trace=FALSE)
  check("Positive scoped admission establishes course admission",isTRUE(z$master$any_admission[k]))
  mm$hospital_days[k]<-0;z<-tf_reconcile(mm,fixture,trace=FALSE)
  check("Contradicted zero total admission days becomes unavailable",is.na(z$master$hospital_days[k])&&isTRUE(z$master$any_admission[k]))
  fixture$flare_hospitalization<-"No";mm$any_admission[k]<-TRUE;z<-tf_reconcile(mm,fixture,trace=FALSE)
  check("Negative scoped admission preserves other course admissions",isTRUE(z$master$any_admission[k]))
  fixture<-cases[1,,drop=FALSE];fixture$overlap_crs<-"Yes";fixture$overlap_icans<-"Yes"
  mm$crs_grade[k]<-mm$icans_grade[k]<-0;mm$crs_1[k]<-mm$icans_1[k]<-FALSE
  mm$crs_2[k]<-mm$crs_3[k]<-mm$icans_2[k]<-mm$icans_3[k]<-FALSE
  z<-tf_reconcile(mm,fixture,trace=FALSE)
  check("Positive CRS with zero master grade preserves unknown severity",isTRUE(z$master$crs_1[k])&&is.na(z$master$crs_grade[k])&&is.na(z$master$crs_2[k])&&is.na(z$master$crs_3[k]))
  check("Positive ICANS with zero master grade preserves unknown severity",isTRUE(z$master$icans_1[k])&&is.na(z$master$icans_grade[k])&&is.na(z$master$icans_2[k])&&is.na(z$master$icans_3[k]))
  fixture<-cases[1,,drop=FALSE];fixture$overlap_infection<-"Yes";mm$infection_flag[k]<-FALSE;mm$infection_count[k]<-0
  for(v in paste0("infection_date_",1:7))mm[[v]][k]<-as.Date(NA_character_)
  z<-tf_reconcile(mm,fixture,trace=FALSE)
  check("Scoped infection invalidates zero count and negative windows",isTRUE(z$master$infection_flag[k])&&is.na(z$master$infection_count[k])&&is.na(z$master$infection_30d[k])&&is.na(z$master$severe_infection_30d[k]))
  fixture<-cases[1,,drop=FALSE];fixture$first_postflare_response<-"MR";fixture$subsequent_response<-"CR";fixture$subsequent_response_on_index<-FALSE
  correction<-data.frame(course_id=fixture$course_id,field="subsequent_response_on_index",value="TRUE",action="set",reason="Fixture index attribution clarified")
  z<-tf_reconcile(master,fixture,correction,trace=FALSE)
  check("Case attribution clarification precedes best-response derivation",z$master$analysis_response[k]=="CR")
  fixture<-cases[1,,drop=FALSE];fixture$step_up_day1<-as.Date("2000-01-01");fixture$flare_onset_date<-as.Date("2000-01-03")
  fixture$first_postflare_response<-"PD";fixture$first_postflare_response_date<-as.Date("2000-01-08");fixture$true_pd_within_14d<-"Yes"
  z<-tf_reconcile(master,fixture,trace=FALSE)
  check("Site-confirmed true early PD supplies endpoint date",z$master$pd_date[k]==as.Date("2000-01-08"))
  fixture$true_pd_within_14d<-"No";z<-tf_reconcile(master,fixture,trace=FALSE)
  check("Apparent assessment PD not promoted to true endpoint",identical(z$master$pd_date[k],master$pd_date[k]))
  correction<-data.frame(course_id=fixture$course_id,field=c("pd_date","pfs_excluded"),value=c(NA,"TRUE"),action=c("invalidate","set"),reason="Fixture unresolved endpoint")
  z<-tf_reconcile(master,fixture,correction,trace=FALSE)
  check("Contradicted progression remains unavailable",is.na(z$master$pd_date[k])&&isTRUE(z$master$pfs_excluded[k]))
  d<-z$master[k,c("course_id","patient_cluster_id","start_date","pd_date","death_date","last_contact_date","analysis_flare","onset_date")]
  ep<-edc_build_endpoint(d,"PFS",pfs_exclude_course_ids=fixture$course_id)
  check("Invalidated PFS not converted to censoring",nrow(ep$data)==0L&&nrow(ep$exclusions)==1L)
  fixture$course_id<-"NEW_CASE";fixture$case_id<-"NEW_RETURN";z<-tf_reconcile(master,fixture,trace=FALSE)
  check("EDC-only return retained descriptively",fixture$course_id %in% z$cases$course_id)
  check("Unverified EDC-only denominator not invented",!fixture$course_id %in% z$cohorts$primary$course_id)
  duplicate_error<-tryCatch({tf_reconcile(master,rbind(cases,cases[1,]));FALSE},error=function(e)TRUE)
  check("Duplicate authoritative return rejected",duplicate_error)
  set.seed(99);idx<-tf_ash_cluster_indices(c("a","a","b","c","c"));multiplicity<-tabulate(idx,nbins=5)
  check("Bootstrap samples all courses within patient blocks",multiplicity[1]==multiplicity[2]&&multiplicity[4]==multiplicity[5])
  set.seed(5);draws<-rnorm(1000)
  low<-tf_ash_bootstrap_inference(.5,c(draws[1:899],rep(NA_real_,101)),1000)
  boundary<-tf_ash_bootstrap_inference(.5,c(draws[1:900],rep(NA_real_,100)),1000)
  check("Bootstrap inference suppressed below minimum",!low$inference_available&&is.na(low$conf_low))
  check("Bootstrap inference available at minimum",boundary$inference_available&&is.finite(boundary$conf_low))
  fit<-tf_ash_firth_fit(data.frame(y=rep(0:1,each=20),v=rep(0:1,each=20)),"y","v","v")
  check("Firth handles complete separation",fit$converged&&is.finite(fit$coefficient))
  fit<-tf_ash_firth_fit(data.frame(y=rep(0,40),v=rep(0:1,each=20)),"y","v","v")
  check("No outcome variation reported nonestimable",!fit$converged)
  do.call(rbind,results)
}
