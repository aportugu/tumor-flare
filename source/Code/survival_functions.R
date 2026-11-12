# Survival analyses for the returned tumor-flare EDC extension.
# These functions operate on an analysis copy and never alter source workbooks.
# Required inputs are one row per course and these columns:
# course_id, patient_cluster_id, start_date, pd_date, death_date,
# last_contact_date, analysis_flare (logical), and onset_date.
# Dates must be Date, POSIXt, or ISO date strings. Exposure is an investigator-
# classified candidate flag, not a claim that the event is independently proven.

edc_survival_date <- function(x, column) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.logical(x) && all(is.na(x))) return(as.Date(rep(NA_character_, length(x))))
  if (is.numeric(x)) {
    stop(column, " must be converted to Date before survival analysis; numeric date origins are not inferred.", call. = FALSE)
  }
  raw <- trimws(as.character(x))
  raw[is.na(x) | raw == ""] <- NA_character_
  out <- suppressWarnings(as.Date(raw, format = "%Y-%m-%d"))
  if (any(!is.na(raw) & (is.na(out) | !grepl("^\\d{4}-\\d{2}-\\d{2}$", raw)))) {
    stop(column, " contains non-ISO date text; resolve it in the import layer.", call. = FALSE)
  }
  out
}

edc_validate_survival_data <- function(data) {
  required <- c("course_id", "patient_cluster_id", "start_date", "pd_date",
                "death_date", "last_contact_date", "analysis_flare", "onset_date")
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Missing survival columns: ", paste(missing, collapse = ", "), call. = FALSE)
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  data$course_id <- trimws(as.character(data$course_id))
  if (any(is.na(data$course_id) | data$course_id == "") || anyDuplicated(data$course_id)) {
    stop("Survival input must contain one row per nonmissing, unique course_id.", call. = FALSE)
  }
  if (!is.logical(data$analysis_flare)) stop("analysis_flare must be logical; unknown classifications must remain NA.", call. = FALSE)
  for (column in c("start_date", "pd_date", "death_date", "last_contact_date", "onset_date")) {
    data[[column]] <- edc_survival_date(data[[column]], column)
  }
  data$patient_cluster_id <- trimws(as.character(data$patient_cluster_id))
  missing_cluster <- is.na(data$patient_cluster_id) | data$patient_cluster_id == ""
  prior_fallback <- if ("cluster_fallback_to_course" %in% names(data)) data$cluster_fallback_to_course else rep(FALSE, nrow(data))
  data$cluster_fallback_to_course <- prior_fallback | missing_cluster
  data$patient_cluster_id[missing_cluster] <- paste0("unlinked_course:", data$course_id[missing_cluster])
  data
}

edc_apply_onset_overrides <- function(data, onset_overrides = NULL) {
  data <- edc_validate_survival_data(data)
  if (!"onset_date_before_override" %in% names(data)) data$onset_date_before_override <- data$onset_date
  if (!"onset_override_applied" %in% names(data)) data$onset_override_applied <- rep(FALSE, nrow(data))
  if (is.null(onset_overrides) || nrow(onset_overrides) == 0L) return(data)
  if (!all(c("course_id", "onset_date") %in% names(onset_overrides))) {
    stop("onset_overrides requires course_id and onset_date columns.", call. = FALSE)
  }
  overrides <- as.data.frame(onset_overrides, stringsAsFactors = FALSE)
  overrides$course_id <- as.character(overrides$course_id)
  if (anyDuplicated(overrides$course_id) || any(is.na(overrides$course_id))) {
    stop("onset_overrides must have unique, nonmissing course IDs.", call. = FALSE)
  }
  if (any(!overrides$course_id %in% data$course_id)) stop("An onset override does not match an input course.", call. = FALSE)
  overrides$onset_date <- edc_survival_date(overrides$onset_date, "onset_overrides$onset_date")
  idx <- match(overrides$course_id, data$course_id)
  if (any(is.na(data$analysis_flare[idx]) | !data$analysis_flare[idx])) {
    stop("Onset overrides can be applied only to classified flare candidates.", call. = FALSE)
  }
  data$onset_date[idx] <- overrides$onset_date
  data$onset_override_applied[idx] <- TRUE
  data
}

edc_min_date <- function(x, y) {
  out <- as.Date(rep(NA_character_, length(x)))
  both <- !is.na(x) & !is.na(y)
  out[both] <- pmin(x[both], y[both])
  out[!is.na(x) & is.na(y)] <- x[!is.na(x) & is.na(y)]
  out[is.na(x) & !is.na(y)] <- y[is.na(x) & !is.na(y)]
  out
}

edc_build_endpoint <- function(data, endpoint = c("PFS", "OS"),
                               pfs_exclude_course_ids = character(), onset_overrides = NULL) {
  endpoint <- match.arg(endpoint)
  data <- edc_apply_onset_overrides(data, onset_overrides)
  data$endpoint <- rep(endpoint, nrow(data))
  data$event_date <- if (endpoint == "PFS") edc_min_date(data$pd_date, data$death_date) else data$death_date
  data$event_after_last_contact <- !is.na(data$event_date) & !is.na(data$last_contact_date) &
    data$event_date > data$last_contact_date
  data$status <- as.integer(!is.na(data$event_date) &
    (is.na(data$last_contact_date) | data$event_date <= data$last_contact_date))
  data$observed_date <- data$last_contact_date
  data$observed_date[data$status == 1L] <- data$event_date[data$status == 1L]
  data$followup_days <- as.numeric(data$observed_date - data$start_date)
  data$time_months <- data$followup_days / 30.4375
  data$endpoint_exclusion_reason <- rep(NA_character_, nrow(data))
  reason <- function(condition, text) {
    use <- condition & is.na(data$endpoint_exclusion_reason)
    data$endpoint_exclusion_reason[use] <<- text
  }
  reason(endpoint == "PFS" & data$course_id %in% pfs_exclude_course_ids,
         "Unresolved true-progression date; excluded from primary updated PFS without an inferred event or censor date")
  reason(is.na(data$start_date), "Missing treatment-start date")
  reason(is.na(data$observed_date), "No valid endpoint observation date")
  reason(!is.na(data$followup_days) & data$followup_days < 0, "Endpoint observation precedes treatment start")
  reason(!is.na(data$followup_days) & data$followup_days == 0, "Zero follow-up; no positive survival interval")
  data$endpoint_eligible <- is.na(data$endpoint_exclusion_reason)
  audit_columns <- c("course_id", "patient_cluster_id", "endpoint", "start_date", "pd_date", "death_date",
                     "last_contact_date", "event_date", "observed_date", "status", "followup_days",
                     "endpoint_exclusion_reason")
  exclusions <- data[!data$endpoint_eligible, audit_columns, drop = FALSE]
  notes <- data[data$event_after_last_contact | data$cluster_fallback_to_course,
                c("course_id", "endpoint", "event_after_last_contact", "cluster_fallback_to_course"), drop = FALSE]
  list(data = data[data$endpoint_eligible, , drop = FALSE], exclusions = exclusions,
       audit = data, notes = notes)
}

edc_prepare_survival_timing <- function(endpoint_data) {
  data <- endpoint_data
  data$onset_day <- as.numeric(data$onset_date - data$start_date)
  data$timing_exclusion_reason <- rep(NA_character_, nrow(data))
  unknown <- is.na(data$analysis_flare)
  missing_onset <- !unknown & data$analysis_flare & is.na(data$onset_day)
  before_start <- !unknown & data$analysis_flare & !is.na(data$onset_day) & data$onset_day < 0
  data$timing_exclusion_reason[unknown] <- "Unknown flare classification; not assigned to the no-flare group"
  data$timing_exclusion_reason[missing_onset] <- "Flare candidate has no usable onset date"
  data$timing_exclusion_reason[before_start] <- "Flare onset precedes treatment start"
  data$timing_eligible <- is.na(data$timing_exclusion_reason)
  exclusions <- data[!data$timing_eligible,
    c("course_id", "patient_cluster_id", "endpoint", "analysis_flare", "onset_date", "onset_day", "timing_exclusion_reason"), drop = FALSE]
  retained <- data[data$timing_eligible, , drop = FALSE]
  retained$onset_after_observation <- retained$analysis_flare & retained$onset_day > retained$followup_days
  retained$endpoint_on_onset_day <- retained$analysis_flare & retained$status == 1L & retained$onset_day == retained$followup_days
  retained$onset_at_observation_boundary <- retained$analysis_flare & retained$onset_day == retained$followup_days
  list(data = retained, exclusions = exclusions)
}

edc_fit_exposure_cox <- function(data, exposure_column, counting_process = FALSE) {
  exposure <- data[[exposure_column]]
  event <- if (counting_process) data$event else data$status
  exposed_events <- sum(event == 1L & exposure == 1L)
  unexposed_events <- sum(event == 1L & exposure == 0L)
  rows <- data.frame(
    n_courses = length(unique(data$course_id)), n_patient_clusters = length(unique(data$patient_cluster_id)),
    exposed_courses = length(unique(data$course_id[exposure == 1L])),
    unexposed_courses = length(unique(data$course_id)) - length(unique(data$course_id[exposure == 1L])),
    courses_with_unexposed_time = length(unique(data$course_id[exposure == 0L])),
    events = sum(event == 1L), exposed_events = exposed_events, unexposed_events = unexposed_events,
    estimate = NA_real_, conf_low = NA_real_, conf_high = NA_real_, p_value = NA_real_,
    model_status = "Not estimable", diagnostic = NA_character_, sparse_note = NA_character_,
    stringsAsFactors = FALSE
  )
  if (exposed_events < 5L) {
    rows$sparse_note <- paste0("Only ", exposed_events, " endpoint event(s) in exposed time/group; estimate is sparse and exploratory")
  }
  reason <- if (nrow(data) == 0L) {
    "No eligible courses"
  } else if (length(unique(exposure)) < 2L) {
    "No exposure variation"
  } else if (length(unique(data$patient_cluster_id)) < 2L) {
    "Fewer than two independent patient clusters"
  } else if (exposed_events == 0L || unexposed_events == 0L) {
    "No endpoint events in one exposure state/group; ordinary Cox estimate is not reported"
  } else NA_character_
  if (!is.na(reason)) {
    rows$diagnostic <- reason
    return(list(table = rows, fit = NULL, warnings = character()))
  }
  rhs <- exposure_column
  lhs <- if (counting_process) "survival::Surv(tstart, tstop, event)" else "survival::Surv(followup_after_landmark_days, status)"
  model_formula <- stats::as.formula(paste(lhs, "~", rhs))
  warnings <- character()
  failure <- NULL
  fit <- tryCatch(
    withCallingHandlers(
      survival::coxph(model_formula, data = data, cluster = patient_cluster_id,
                      ties = "efron", robust = TRUE, x = TRUE, model = TRUE),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) { failure <<- conditionMessage(e); NULL }
  )
  bad_warning <- any(grepl("converg|infinite|singular|iteration|ran out", warnings, ignore.case = TRUE))
  if (is.null(fit) || bad_warning) {
    rows$diagnostic <- paste(c(failure, warnings), collapse = "; ")
    return(list(table = rows, fit = NULL, warnings = unique(c(failure, warnings))))
  }
  log_hr <- unname(stats::coef(fit)[1L])
  variance <- unname(stats::vcov(fit)[1L, 1L])
  se <- if (is.finite(variance) && variance > 0) sqrt(variance) else NA_real_
  estimate <- exp(log_hr)
  ci <- exp(log_hr + c(-1, 1) * stats::qnorm(0.975) * se)
  if (any(!is.finite(c(log_hr, variance, se, estimate, ci))) || variance <= 0) {
    rows$diagnostic <- "Nonfinite coefficient/confidence interval or nonpositive robust variance"
    return(list(table = rows, fit = NULL, warnings = unique(warnings)))
  }
  rows$estimate <- estimate
  rows$conf_low <- ci[1L]
  rows$conf_high <- ci[2L]
  rows$p_value <- 2 * stats::pnorm(-abs(log_hr / se))
  rows$model_status <- "Estimated (exploratory)"
  rows$diagnostic <- if (length(warnings)) paste(unique(warnings), collapse = "; ") else "Exposure-only Cox; Efron ties; patient-cluster robust standard errors"
  list(table = rows, fit = fit, warnings = unique(warnings))
}

edc_time_dependent_cox <- function(data, endpoint = c("PFS", "OS"),
                                   pfs_exclude_course_ids = character(), onset_overrides = NULL) {
  endpoint <- match.arg(endpoint)
  endpoint_result <- edc_build_endpoint(data, endpoint, pfs_exclude_course_ids, onset_overrides)
  timing <- edc_prepare_survival_timing(endpoint_result$data)
  eligible <- timing$data
  if (nrow(eligible)) {
    eligible$event <- eligible$status
    intervals <- survival::tmerge(data1 = eligible, data2 = eligible, id = course_id,
      tstop = followup_days, event = event(followup_days, status))
    updates <- eligible[eligible$analysis_flare, c("course_id", "onset_day"), drop = FALSE]
    if (nrow(updates)) {
      intervals <- survival::tmerge(intervals, updates, id = course_id, flare_td = tdc(onset_day))
    } else intervals$flare_td <- 0L
    # tmerge uses (tstart, tstop]. An event at the onset boundary remains in
    # the preceding interval. No unobserved time-of-day or epsilon is added.
    tmerge_counts <- attr(intervals, "tcount")
    intervals <- as.data.frame(intervals)
    intervals$flare_td <- as.integer(intervals$flare_td)
  } else {
    intervals <- eligible
    intervals$tstart <- numeric()
    intervals$tstop <- numeric()
    intervals$event <- integer()
    intervals$flare_td <- integer()
    tmerge_counts <- NULL
  }
  fit <- edc_fit_exposure_cox(intervals, "flare_td", counting_process = TRUE)
  fit$table$endpoint <- endpoint
  fit$table$analysis_type <- "Time-dependent Cox"
  fit$table$landmark_day <- NA_real_
  fit$table$exposed_person_days <- sum((intervals$tstop - intervals$tstart)[intervals$flare_td == 1L])
  fit$table$unexposed_person_days <- sum((intervals$tstop - intervals$tstart)[intervals$flare_td == 0L])
  list(table = fit$table, intervals = intervals, fit = fit$fit, warnings = fit$warnings,
       endpoint_data = endpoint_result$data, endpoint_exclusions = endpoint_result$exclusions,
       endpoint_audit = endpoint_result$audit, endpoint_notes = endpoint_result$notes,
       timing_exclusions = timing$exclusions, timing_data = eligible, tmerge_counts = tmerge_counts)
}

edc_landmark_analysis <- function(data, endpoint = c("PFS", "OS"), landmark_days = c(30, 60),
                                  pfs_exclude_course_ids = character(), onset_overrides = NULL) {
  endpoint <- match.arg(endpoint)
  if (!length(landmark_days) || any(!is.finite(landmark_days) | landmark_days < 0) || anyDuplicated(landmark_days)) {
    stop("landmark_days must be unique, nonnegative finite days.", call. = FALSE)
  }
  endpoint_result <- edc_build_endpoint(data, endpoint, pfs_exclude_course_ids, onset_overrides)
  timing <- edc_prepare_survival_timing(endpoint_result$data)
  results <- lapply(landmark_days, function(day) {
    current <- timing$data
    current$landmark_day <- rep(day, nrow(current))
    current$landmark_exclusion_reason <- rep(NA_character_, nrow(current))
    event_before <- current$status == 1L & current$followup_days <= day
    inadequate <- current$followup_days <= day & !event_before
    current$landmark_exclusion_reason[event_before] <- "Endpoint at or before landmark"
    current$landmark_exclusion_reason[inadequate] <- "Follow-up does not extend beyond landmark"
    excluded <- current[!is.na(current$landmark_exclusion_reason),
      c("course_id", "patient_cluster_id", "endpoint", "landmark_day", "status", "followup_days", "landmark_exclusion_reason"), drop = FALSE]
    current <- current[is.na(current$landmark_exclusion_reason), , drop = FALSE]
    current$flare_by_landmark <- as.integer(current$analysis_flare & current$onset_day <= day)
    current$followup_after_landmark_days <- current$followup_days - day
    current$time_after_landmark_months <- current$followup_after_landmark_days / 30.4375
    fit <- edc_fit_exposure_cox(current, "flare_by_landmark")
    fit$table$endpoint <- endpoint
    fit$table$analysis_type <- "Landmark Cox"
    fit$table$landmark_day <- day
    list(table = fit$table, data = current, exclusions = excluded, fit = fit$fit, warnings = fit$warnings)
  })
  names(results) <- paste0("day", landmark_days)
  list(table = dplyr::bind_rows(lapply(results, `[[`, "table")),
       data = dplyr::bind_rows(lapply(results, `[[`, "data")),
       exclusions = dplyr::bind_rows(lapply(results, `[[`, "exclusions")),
       results = results, endpoint_exclusions = endpoint_result$exclusions,
       timing_exclusions = timing$exclusions)
}

edc_landmark_km_plot <- function(landmark_data, endpoint = NULL, landmark_day = NULL,
                                title = NULL, show_confidence = FALSE) {
  required <- c("endpoint", "landmark_day", "time_after_landmark_months", "status", "flare_by_landmark")
  if (!all(required %in% names(landmark_data))) stop("Pass landmark_data returned by edc_landmark_analysis or edc_run_survival_analyses.", call. = FALSE)
  data <- landmark_data
  if (!is.null(endpoint)) data <- data[data$endpoint == endpoint, , drop = FALSE]
  if (!is.null(landmark_day)) data <- data[data$landmark_day == landmark_day, , drop = FALSE]
  if (nrow(data) && (length(unique(data$endpoint)) != 1L || length(unique(data$landmark_day)) != 1L)) {
    stop("Select exactly one endpoint and landmark_day for a Kaplan-Meier plot.", call. = FALSE)
  }
  if (nrow(data) == 0L) {
    return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = "No eligible landmark follow-up") +
      ggplot2::theme_void() + ggplot2::labs(title = title))
  }
  endpoint <- data$endpoint[1L]
  day <- data$landmark_day[1L]
  labels <- c(paste0("No flare by day ", day), paste0("Flare by day ", day))
  data$landmark_group <- factor(data$flare_by_landmark, levels = c(0, 1), labels = labels)
  curves <- lapply(levels(droplevels(data$landmark_group)), function(group) {
    current <- data[data$landmark_group == group, , drop = FALSE]
    fit <- survival::survfit(survival::Surv(time_after_landmark_months, status) ~ 1, data = current,
      cluster = patient_cluster_id, robust = TRUE)
    data.frame(time = c(0, fit$time), survival = c(1, fit$surv),
               lower = c(1, fit$lower), upper = c(1, fit$upper),
               group = paste0(group, " (n=", nrow(current), ")"), stringsAsFactors = FALSE)
  })
  curve <- dplyr::bind_rows(curves)
  palette <- stats::setNames(c("#1b4965", "#b23a48")[match(levels(droplevels(data$landmark_group)), labels)], unique(curve$group))
  plot <- ggplot2::ggplot(curve, ggplot2::aes(time, survival, color = group))
  if (show_confidence) {
    # Duplicate each event boundary to shade a right-continuous step band,
    # rather than interpolating uncertainty between observed event times.
    band <- dplyr::bind_rows(lapply(split(curve, curve$group), function(z) {
      n <- nrow(z)
      idx <- rep(seq_len(n), each = 2L)[seq_len(2L*n-1L)]
      stepped <- z[idx, , drop = FALSE]
      stepped$time <- c(z$time[1L], rep(z$time[-1L], each = 2L))
      stepped
    }))
    plot <- plot + ggplot2::geom_ribbon(data = band,
      ggplot2::aes(x = time, ymin = lower, ymax = upper, fill = group, group = group),
      inherit.aes = FALSE, color = NA, alpha = 0.12, show.legend = FALSE, na.rm = TRUE) +
      ggplot2::scale_fill_manual(values = palette)
  }
  plot <- plot +
    ggplot2::geom_step(linewidth = 0.8) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1)) +
    ggplot2::labs(title = if (is.null(title)) paste0(endpoint, " after the day-", day, " landmark") else title,
                  x = paste0("Months after day ", day), y = paste0(endpoint, " probability"), color = NULL,
                  caption = "Exploratory fixed-landmark comparison; later flares remain in the no-flare-by-landmark group.") +
    ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(legend.position = "bottom")
  if (show_confidence) {
    plot <- plot + ggplot2::geom_step(ggplot2::aes(y = lower), linewidth = 0.35, linetype = "dotted", show.legend = FALSE, na.rm = TRUE) +
      ggplot2::geom_step(ggplot2::aes(y = upper), linewidth = 0.35, linetype = "dotted", show.legend = FALSE, na.rm = TRUE)
  }
  plot
}

edc_run_survival_analyses <- function(data, pfs_exclude_course_ids = character(),
                                      landmark_days = c(30, 60), onset_overrides = NULL,
                                      scenario = "Primary") {
  data <- edc_apply_onset_overrides(data, onset_overrides)
  # The override is applied once, retaining its provenance in every returned row.
  td <- lapply(c("PFS", "OS"), function(endpoint) edc_time_dependent_cox(data, endpoint, pfs_exclude_course_ids))
  lm <- lapply(c("PFS", "OS"), function(endpoint) edc_landmark_analysis(data, endpoint, landmark_days, pfs_exclude_course_ids))
  names(td) <- names(lm) <- c("PFS", "OS")
  combine <- function(results, field) {
    out <- dplyr::bind_rows(lapply(results, `[[`, field))
    out$scenario <- rep(scenario, nrow(out))
    out
  }
  warning_rows <- list()
  for (endpoint in names(td)) {
    if (length(td[[endpoint]]$warnings)) warning_rows[[length(warning_rows) + 1L]] <-
      data.frame(endpoint = endpoint, analysis_type = "Time-dependent Cox", landmark_day = NA_real_, warning = td[[endpoint]]$warnings)
    for (result in lm[[endpoint]]$results) if (length(result$warnings)) warning_rows[[length(warning_rows) + 1L]] <-
      data.frame(endpoint = endpoint, analysis_type = "Landmark Cox", landmark_day = result$table$landmark_day, warning = result$warnings)
  }
  warning_table <- if (length(warning_rows)) dplyr::bind_rows(warning_rows) else
    data.frame(endpoint = character(), analysis_type = character(), landmark_day = numeric(), warning = character())
  warning_table$scenario <- rep(scenario, nrow(warning_table))
  list(
    time_dependent_table = combine(td, "table"), landmark_table = combine(lm, "table"),
    intervals = combine(td, "intervals"), endpoint_data = combine(td, "endpoint_data"),
    endpoint_audit = combine(td, "endpoint_audit"), endpoint_exclusions = combine(td, "endpoint_exclusions"),
    endpoint_notes = combine(td, "endpoint_notes"), timing_exclusions = combine(td, "timing_exclusions"),
    timing_data = combine(td, "timing_data"), landmark_data = combine(lm, "data"),
    landmark_exclusions = combine(lm, "exclusions"), warnings = warning_table,
    fits = list(time_dependent = lapply(td, `[[`, "fit"),
                landmark = lapply(lm, function(x) lapply(x$results, `[[`, "fit"))),
    details = list(time_dependent = td, landmark = lm)
  )
}
