# Portable analytical entry point. All source-specific rules are controlled data.
# Usage: Rscript run_analysis.R CONTROLLED_DATA_DIRECTORY OUTPUT_DIRECTORY [B]
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<2L)stop("Expected controlled data and output directories")
data_dir<-normalizePath(args[1],mustWork=TRUE)
output_dir<-args[2];dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
output_dir<-normalizePath(output_dir,mustWork=TRUE)
B<-if(length(args)>=3L)as.integer(args[3])else 1000L
if(!is.finite(B)||B<1L)stop("Bootstrap draw count must be positive")
RNGkind("Mersenne-Twister","Inversion","Rejection")
library_path<-Sys.getenv("TF_R_LIBRARY","")
if(nzchar(library_path)).libPaths(c(normalizePath(library_path,mustWork=TRUE),.libPaths()))
script_arg<-grep("^--file=",commandArgs(),value=TRUE)
code_dir<-dirname(normalizePath(sub("^--file=","",script_arg[1]),mustWork=TRUE))
for(file in c("analysis_functions.R","survival_functions.R","reconciliation_functions.R","model_functions.R","summary_functions.R"))source(file.path(code_dir,file))
for(package in c("dplyr","tidyr","tibble","survival","logistf","jsonlite"))if(!requireNamespace(package,quietly=TRUE))stop("Missing dependency: ",package)
master<-tf_read_controlled_csv(file.path(data_dir,"master_courses.csv"))
cases<-tf_read_controlled_csv(file.path(data_dir,"edc_cases.csv"))
rules<-utils::read.csv(file.path(data_dir,"verified_corrections.csv"),check.names=FALSE,stringsAsFactors=FALSE,na.strings="")
source_path<-file.path(data_dir,"master_source_records.csv")
master_sources<-if(file.exists(source_path))utils::read.csv(source_path,stringsAsFactors=FALSE,na.strings="")else NULL
analysis<-tf_reconcile(master,cases,rules,master_sources)
scenario_path<-file.path(data_dir,"survival_scenario_rules.csv")
analysis$survival_scenario_rules<-if(file.exists(scenario_path))utils::read.csv(scenario_path,
  check.names=FALSE,stringsAsFactors=FALSE,na.strings="")else NULL
analysis$project_root<-data_dir
counts<-sort(table(analysis$cases$site_key),decreasing=TRUE)
analysis$excluded_site_keys<-names(counts)[1]
specs<-tf_ash_specs()
ids<-character()
for(scenario in c("primary","strict","without_dominant_site")) {
  ids<-c(ids,paste(scenario,specs$field,sep="_"))
  for(threshold in c("pr","vgpr","cr"))ids<-c(ids,paste(scenario,paste0("response_",threshold),c("unadjusted","adjusted"),sep="_"))
}
seed_path<-file.path(data_dir,"model_seed_manifest.csv")
if(!file.exists(seed_path))seed_path<-file.path(code_dir,"model_seed_manifest.csv")
if(file.exists(seed_path)) {
  seeds<-utils::read.csv(seed_path,stringsAsFactors=FALSE)
  if(!setequal(ids,seeds$model_id)||anyDuplicated(seeds$model_id)||anyNA(seeds$seed))stop("Model seed manifest does not match requested models")
} else stop("A fixed model-ID seed manifest is required")
analysis$model_seeds<-seeds
tf_export_descriptive(analysis,output_dir)
source(file.path(code_dir,"ground_truth_tests.R"))
checks<-tf_ground_truth_tests(master,cases,rules,analysis)
utils::write.csv(checks,file.path(output_dir,"ground_truth_verification.csv"),row.names=FALSE)
if(any(checks$status!="PASS"))stop("Ground-truth reconciliation validation failed")
cat("Reconciliation passed",nrow(checks),"checks. Starting",length(ids),"Firth models with",B,"patient-cluster draws each.\n")
models<-tf_build_ash_extension(analysis,output_dir,B=B,seed=20261106L,progress=TRUE)
tf_export_models(analysis,models,output_dir)
tf_write_analysis_summary(analysis,models,output_dir)
utils::write.csv(seeds,file.path(output_dir,"model_seed_manifest.csv"),row.names=FALSE)
writeLines(c(paste("Bootstrap draws:",B),paste("RNG:",paste(RNGkind(),collapse="; ")),
  "Canonical returned-case fields use authoritative EDC/site data, including unavailable values.",
  "No date, response, event, or covariate imputation.","",capture.output(sessionInfo())),
  file.path(output_dir,"session_info.txt"))
cat("Canonical analysis complete.\n")
