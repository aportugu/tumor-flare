# Canonical fieldwise reconciliation. Controlled inputs carry all case-specific
# corrections; this module contains no case identities or clinical narratives.
tf_supported_products <- c("teclistamab","talquetamab","elranatamab","linvoseltamab")

tf_read_controlled_csv <- function(path) {
  d <- utils::read.csv(path,check.names=FALSE,stringsAsFactors=FALSE,na.strings="")
  dates <- intersect(c("start_date","full_dose_date","pd_date","death_date","last_contact_date",
    "step_up_day1","flare_onset_date","flare_peak_date","last_flare_assessment","flare_resolution_date",
    "last_bsab_dose_before_onset","flare_admission_date","bsab_resume_date","local_outcome_date",
    "first_postflare_response_date","subsequent_response_date","alternative_assessment_date",
    "onset_lower_date","onset_upper_date",paste0("infection_date_",1:7)),names(d))
  for(field in dates) d[[field]] <- edc_survival_date(d[[field]],field)
  d
}

tf_source_state <- function(x) {
  v <- as.character(x)
  ifelse(is.na(v)|!nzchar(trimws(v)),"Blank",ifelse(tolower(trimws(v)) %in%
    c("unknown","ne","ne/unknown","not reassessed","not applicable"),v,"Recorded"))
}

tf_value_text <- function(x) {
  out <- as.character(x); out[is.na(out)] <- ""; out
}

tf_best_index_response <- function(first,later,first_on_index=TRUE,later_on_index=TRUE) {
  levels <- tf_response_levels[1:7]
  candidates <- c(if(isTRUE(first_on_index))as.character(first),if(isTRUE(later_on_index))as.character(later))
  supported <- candidates[candidates %in% levels]
  if(!length(supported)) return("NE/Unknown")
  levels[min(match(supported,levels))]
}

tf_recalculate_infection_windows <- function(x) {
  if(!all(paste0("infection_date_",1:7) %in% names(x)))return(x)
  for(i in seq_len(nrow(x))) {
    dates <- as.Date(vapply(1:7,function(k)as.numeric(x[[paste0("infection_date_",k)]][i]),numeric(1)),origin="1970-01-01")
    severe <- vapply(1:7,function(k)isTRUE(x[[paste0("infection_severe_",k)]][i]),logical(1))
    known <- vapply(1:7,function(k)isTRUE(x[[paste0("infection_severity_known_",k)]][i]),logical(1))
    within <- !is.na(dates)&!is.na(x$start_date[i])&dates>=x$start_date[i]&dates<=x$start_date[i]+30
    flag<-x$infection_flag[i];count<-x$infection_count[i]
    documented<-any(!is.na(dates))||isTRUE(flag)||(!is.na(count)&&count>0)
    explicit_no<-(identical(flag,FALSE)||(!is.na(count)&&count==0))&&!documented
    complete<-explicit_no||(!is.na(count)&&count>=0&&count<=7&&sum(!is.na(dates))==count)
    x$infection_30d[i]<-if(any(within))TRUE else if(explicit_no||(complete&&!is.na(x$start_date[i])))FALSE else NA
    x$severe_infection_30d[i]<-if(any(within&severe))TRUE else if(explicit_no||(complete&&!is.na(x$start_date[i])&&all(known[within])))FALSE else NA
  }
  x
}

tf_reconcile <- function(master,cases,corrections=NULL,master_sources=NULL,trace=TRUE) {
  if(anyDuplicated(master$course_id)||anyDuplicated(cases$course_id))stop("Course identities must be unique")
  master$eligibility_verified<-TRUE
  unmatched<-which(!cases$course_id %in% master$course_id)
  if(length(unmatched))for(i in unmatched) {
    row<-master[1,,drop=FALSE]
    for(field in names(row))row[[field]][1]<-NA
    row$course_id<-cases$course_id[i]
    row$patient_cluster_id<-paste0("unlinked_course:",cases$course_id[i])
    row$site<-cases$site[i];row$site_key<-cases$site_key[i]
    row$product<-cases$bsab_product[i];row$product_key<-tolower(cases$bsab_product[i])
    row$original_response<-row$analysis_response<-"NE/Unknown"
    row$pfs_excluded<-FALSE
    row$eligibility_verified<-FALSE
    master<-rbind(master,row)
  }
  original <- master
  x <- master
  x$edc_received <- x$course_id %in% cases$course_id
  x$investigator_confidence <- NA_character_
  x$onset_date <- as.Date(NA_character_)
  x$pfs_excluded <- as.logical(x$pfs_excluded)
  if(anyNA(x$pfs_excluded))x$pfs_excluded[is.na(x$pfs_excluded)]<-FALSE
  audit<-list();provenance<-list()
  add <- function(ids,field,old,value,source,source_field,reason,reference=NULL) {
    if(!trace)return(invisible(NULL))
    supporting<-if(source=="Master")master_sources else cases
    pos<-if(!is.null(supporting))match(ids,supporting$course_id)else rep(NA_integer_,length(ids))
    get_ref<-function(name)if(!is.null(supporting)&&name %in% names(supporting))
      tf_value_text(supporting[[name]][pos])else rep("",length(ids))
    rows<-data.frame(course_id=ids,field=field,source=source,source_field=source_field,
      source_value=tf_value_text(value),source_state=tf_source_state(value),
      original_master_value=tf_value_text(old),canonical_value=tf_value_text(value),
      changed=tf_value_text(old)!=tf_value_text(value),reason=reason,
      source_file=get_ref("source_file"),source_sheet=get_ref("source_sheet"),
      source_row=get_ref("source_row"),clarification_pointer=get_ref("clarification_pointer"),controlled_rule_pointer="",
      supporting_references="",stringsAsFactors=FALSE)
    if(source=="EDC"&&field %in% c("analysis_response","case.first_response_on_index","case.subsequent_response_on_index"))
      rows$controlled_rule_pointer<-get_ref("response_attribution_pointer")
    if(source=="EDC"&&startsWith(field,"case.")&&!is.null(corrections)&&nrow(corrections)) {
      rr<-corrections[corrections$field==sub("^case\\.","",field),,drop=FALSE]
      rp<-match(rows$course_id,rr$course_id)
      selected<-!is.na(rp)
      for(name in intersect(c("source_file","source_sheet","source_row","clarification_pointer","controlled_rule_pointer"),names(rr)))
        rows[[name]][selected]<-tf_value_text(rr[[name]][rp[selected]])
    }
    if(!is.null(reference))for(name in intersect(c("source_file","source_sheet","source_row",
      "clarification_pointer","controlled_rule_pointer","supporting_references"),names(reference)))
      rows[[name]]<-tf_value_text(reference[[name]])
    invalidated<-is.na(value)&grepl("unavailable|contradict|invalid|unresolved",reason,ignore.case=TRUE)
    rows$source_state[invalidated]<-"Unavailable after source reconciliation"
    provenance[[length(provenance)+1L]]<<-rows
    audit[[length(audit)+1L]]<<-rows[rows$changed,,drop=FALSE]
  }
  master_fields<-setdiff(names(x),c("course_id","patient_cluster_id","site_key","edc_received"))
  for(field in master_fields) add(x$course_id,field,x[[field]],x[[field]],"Master",field,"Master supplies a separate concept that site collection does not reassess")
  case_corrections<-integer()
  if(!is.null(corrections)&&nrow(corrections))for(i in seq_len(nrow(corrections))) {
    r<-corrections[i,,drop=FALSE];field<-as.character(r$field)
    if(!field %in% names(cases)||field %in% c("course_id","case_id","source_file","source_sheet","source_row"))next
    k<-match(r$course_id,cases$course_id)
    if(is.na(k))stop("Site case clarification has no authoritative return")
    old<-cases[[field]][k]
    value<-if(r$action=="invalidate")NA else as.character(r$value)
    if(inherits(old,"Date"))value<-as.Date(value)else if(is.logical(old))value<-as.logical(value)else if(is.numeric(old))value<-as.numeric(value)
    cases[[field]][k]<-value
    add(r$course_id,paste0("case.",field),old,value,"Verified site clarification",field,as.character(r$reason),r)
    case_corrections<-c(case_corrections,i)
  }
  for(i in seq_len(nrow(cases))) {
    j<-match(cases$course_id[i],x$course_id)
    confidence<-as.character(cases$investigator_confidence[i])
    flare<-if(confidence %in% c("Definite","Probable","Possible"))TRUE else if(!is.na(confidence)&&confidence=="Not tumor flare")FALSE else NA
    mappings<-list(start_date=cases$step_up_day1[i],product=as.character(cases$bsab_product[i]),
      analysis_flare=flare,onset_date=cases$flare_onset_date[i],investigator_confidence=confidence,
      analysis_response=tf_best_index_response(cases$first_postflare_response[i],cases$subsequent_response[i],
        cases$first_response_on_index[i],cases$subsequent_response_on_index[i]))
    source_fields<-c(start_date="step_up_day1",product="bsab_product",analysis_flare="investigator_confidence",
      onset_date="flare_onset_date",investigator_confidence="investigator_confidence",
      analysis_response="first_postflare_response; subsequent_response; first_response_on_index; subsequent_response_on_index")
    for(field in names(mappings)) {
      old<-x[[field]][j];x[[field]][j]<-mappings[[field]]
      add(x$course_id[j],field,old,mappings[[field]],"EDC",source_fields[[field]],
        "Authoritative returned value, including explicit uncertainty or blank; no master backfill")
    }
    for(field in setdiff(names(cases),c("course_id","case_id","source_file","source_sheet","source_row"))) {
      add(x$course_id[j],paste0("case.",field),NA_character_,cases[[field]][i],"EDC",field,"Scoped index-event field stays separate from whole-course measures")
    }
    # A positive scoped event establishes occurrence somewhere in the course.
    # A negative scoped event cannot negate other events elsewhere in the course.
    implications<-c(overlap_crs="crs_1",overlap_icans="icans_1",flare_hospitalization="any_admission",
      management_steroids="steroids",overlap_infection="any_infection")
    for(field in names(implications)) if(isTRUE(cases[[field]][i]=="Yes")) {
      target<-implications[[field]];old<-x[[target]][j];x[[target]][j]<-TRUE
      add(x$course_id[j],target,old,TRUE,"EDC",field,"Positive index event establishes corresponding whole-course occurrence")
    }
    if(isTRUE(cases$overlap_crs[i]=="Yes")&&isTRUE(x$crs_grade[j]==0)) {
      add(x$course_id[j],"crs_grade",x$crs_grade[j],NA_real_,"EDC","overlap_crs","Occurrence confirmed; maximum grade is not supplied by this field")
      x$crs_grade[j]<-NA_real_
      for(target in c("crs_2","crs_3")) {
        add(x$course_id[j],target,x[[target]][j],NA,"EDC","overlap_crs","Occurrence confirmed; higher severity unavailable")
        x[[target]][j]<-NA
      }
    }
    if(isTRUE(cases$overlap_icans[i]=="Yes")&&isTRUE(x$icans_grade[j]==0)) {
      add(x$course_id[j],"icans_grade",x$icans_grade[j],NA_real_,"EDC","overlap_icans","Occurrence confirmed; maximum grade is not supplied by this field")
      x$icans_grade[j]<-NA_real_
      for(target in c("icans_2","icans_3")) {
        add(x$course_id[j],target,x[[target]][j],NA,"EDC","overlap_icans","Occurrence confirmed; higher severity unavailable")
        x[[target]][j]<-NA
      }
    }
    if(isTRUE(cases$flare_hospitalization[i]=="Yes")&&isTRUE(x$hospital_days[j]==0)) {
      add(x$course_id[j],"hospital_days",x$hospital_days[j],NA_real_,"EDC","flare_hospitalization","Admission confirmed; whole-course total days unavailable rather than replaced by flare-only days")
      x$hospital_days[j]<-NA_real_
    }
    if(isTRUE(cases$overlap_infection[i]=="Yes")) {
      add(x$course_id[j],"infection_flag",x$infection_flag[j],TRUE,"EDC","overlap_infection","Infection occurrence confirmed; scoped timing and severity do not establish all-course counts")
      x$infection_flag[j]<-TRUE
      if(isTRUE(x$infection_count[j]==0)) {
        add(x$course_id[j],"infection_count",x$infection_count[j],NA_real_,"EDC","overlap_infection","Contradicted zero count becomes unavailable")
        x$infection_count[j]<-NA_real_
      }
      if(isTRUE(x$severe_infection[j]==FALSE)) {
        add(x$course_id[j],"severe_infection",x$severe_infection[j],NA,"EDC","overlap_infection","Severity unavailable for the newly confirmed infection")
        x$severe_infection[j]<-NA
      }
    }
    # The early true-PD field plus an explicitly dated PD assessment in that
    # window establishes a true-progression date; apparent PD alone does not.
    date<-cases$first_postflare_response_date[i]
    onset<-cases$flare_onset_date[i]
    if(isTRUE(cases$true_pd_within_14d[i]=="Yes")&&isTRUE(cases$first_postflare_response[i]=="PD")&&
       !is.na(date)&&!is.na(onset)&&date>=onset&&date<=onset+14) {
      old<-x$pd_date[j];x$pd_date[j]<-date
      add(x$course_id[j],"pd_date",old,date,"EDC","true_pd_within_14d; first_postflare_response_date",
        "Site confirms true early progression and supplies the corresponding assessment date")
    }
  }
  if(!is.null(corrections)&&nrow(corrections))for(i in setdiff(seq_len(nrow(corrections)),case_corrections)) {
    r<-corrections[i,,drop=FALSE];j<-match(r$course_id,x$course_id)
    if(is.na(j))stop("Correction refers to an unknown normalized course")
    field<-as.character(r$field)
    in_cases<-field %in% names(cases)&&!field %in% names(x)
    if(in_cases) {
      k<-match(r$course_id,cases$course_id)
      if(is.na(k))stop("Case correction has no authoritative returned case")
      old<-cases[[field]][k]
    } else {
      if(!field %in% names(x))stop("Unknown correction field: ",field)
      old<-x[[field]][j]
    }
    value<-if(r$action=="invalidate")NA else as.character(r$value)
    if(inherits(old,"Date"))value<-as.Date(value) else if(is.logical(old))value<-as.logical(value) else if(is.numeric(old))value<-as.numeric(value)
    if(in_cases)cases[[field]][k]<-value else x[[field]][j]<-value
    add(x$course_id[j],if(in_cases)paste0("case.",field)else field,old,value,"Verified site clarification",field,as.character(r$reason),r)
  }
  derive<-function(field,value,inputs,reason) {
    if(!trace){x[[field]]<<-value;return(invisible(NULL))}
    old<-if(field %in% names(x))x[[field]]else rep(NA_character_,nrow(x))
    previous<-dplyr::bind_rows(provenance)
    previous<-previous[!duplicated(paste(previous$course_id,previous$field),fromLast=TRUE),,drop=FALSE]
    previous<-previous[previous$field %in% inputs,,drop=FALSE]
    by_course<-split(previous,previous$course_id)
    references<-vapply(x$course_id,function(id) {
      r<-by_course[[id]]
      if(is.null(r))return("")
      text<-paste0(r$field,": ",r$source_file," [",r$source_sheet,", row ",r$source_row,"]",
        ifelse(nzchar(r$clarification_pointer),paste0("; ",r$clarification_pointer),""),
        ifelse(nzchar(r$controlled_rule_pointer),paste0("; ",r$controlled_rule_pointer),""))
      paste(unique(text),collapse=" | ")
    },character(1))
    x[[field]]<<-value
    add(x$course_id,field,old,value,"Derived from canonical source fields",paste(inputs,collapse="; "),reason,
      list(source_file="",source_sheet="",source_row="",supporting_references=references))
  }
  derive("product_key",tolower(as.character(x$product)),"product","Normalized authoritative course product")
  bcma<-x$product_key %in% setdiff(tf_supported_products,"talquetamab")
  rank_group<-ifelse(bcma,x$patient_cluster_id,paste0("nonbcma:",x$course_id))
  idx<-order(rank_group,is.na(x$start_date),x$start_date,is.na(x$full_dose_date),x$full_dose_date,x$course_id)
  first<-rep(TRUE,nrow(x));first[idx]<-!duplicated(rank_group[idx])
  derive("keep_first_bcma",first,c("product","start_date","full_dose_date"),
    "First BCMA course is ordered within the supplied patient cluster by authoritative step-up date, then full-dose date and course ID; GPRC5D courses remain eligible")
  infection<-tf_recalculate_infection_windows(x)
  infection_inputs<-c("start_date","infection_flag","infection_count",paste0("infection_date_",1:7),
    paste0("infection_severe_",1:7),paste0("infection_severity_known_",1:7))
  for(field in c("any_infection","severe_infection"))derive(field,x[[field]],
    c(infection_inputs,paste0("case.overlap_infection")),
    "Final course-level occurrence/severity after scoped positive-event reconciliation; unknown new-event severity stays unavailable")
  for(field in c("infection_30d","severe_infection_30d"))derive(field,infection[[field]],infection_inputs,
    "Recomputed after source precedence; positive dated occurrence establishes the window, and incomplete contradicted counts never establish a negative window")
  cases$analysis_flare<-x$analysis_flare[match(cases$course_id,x$course_id)]
  cases$best_index_response<-x$analysis_response[match(cases$course_id,x$course_id)]
  cases$bridging_any<-x$bridging_any[match(cases$course_id,x$course_id)]
  cases$keep_first_bcma<-x$keep_first_bcma[match(cases$course_id,x$course_id)]
  membership<-x$product_key %in% tf_supported_products&x$keep_first_bcma&!x$bridging_any&x$eligibility_verified
  membership[is.na(membership)]<-FALSE
  derive("comparative_eligible",membership,c("product_key","keep_first_bcma","bridging_any","eligibility_verified"),
    "Supported product, first BCMA/all GPRC5D, nonbridging and verified linkage eligibility")
  eligible<-x[membership,,drop=FALSE]
  primary<-eligible[!is.na(eligible$analysis_flare),,drop=FALSE]
  strict<-primary[is.na(primary$investigator_confidence)|primary$investigator_confidence!="Possible",,drop=FALSE]
  cases$comparative_included<-cases$course_id %in% primary$course_id
  records<-if(trace)dplyr::bind_rows(provenance)else data.frame(course_id=character(),field=character(),changed=logical())
  records$record_order<-seq_len(nrow(records))
  records$is_final_record<-!duplicated(paste(records$course_id,records$field),fromLast=TRUE)
  changes<-records[records$changed,,drop=FALSE]
  list(master=x,cases=cases,eligible=eligible,cohorts=list(primary=primary,strict=strict),
    provenance=records,register=changes,original=original)
}
