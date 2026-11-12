# Generic normalized-data analysis helpers.
tf_response_levels <- c("sCR","CR","VGPR","PR","MR","SD","PD","NE/Unknown","Blank")
tf_product_labels <- c(teclistamab="Teclistamab",talquetamab="Talquetamab",elranatamab="Elranatamab",linvoseltamab="Linvoseltamab")
tf_norm <- function(x) tolower(trimws(gsub("\\s+", " ", as.character(x))))

tf_number <- function(x) suppressWarnings(as.numeric(as.character(x)))

tf_binary <- function(x) {
  z <- tf_norm(x)
  out <- rep(NA, length(z)); out[z %in% c("1", "yes", "y", "true")] <- TRUE
  out[z %in% c("0", "no", "n", "false")] <- FALSE
  out
}

tf_category <- function(x) {
  z <- if(is.logical(x)) ifelse(is.na(x),NA_character_,ifelse(x,"Yes","No")) else as.character(x)
  z[is.na(z) | trimws(z) == ""] <- "Blank"; z
}

tf_response <- function(x) {
  z <- tf_norm(x); out <- rep("NE/Unknown", length(z))
  out[is.na(z) | z == ""] <- "Blank"
  # Keep a leading formal category in an annotation without upgrading it from
  # a deeper category mentioned later. Qualified/ambiguous claims stay unknown.
  for (category in c("PD", "SD", "MR", "PR", "CR", "VGPR", "sCR")) {
    pattern <- if (category == "VGPR") "^(vgpr|vpgr)$" else paste0("^", tolower(category), "$")
    out[!is.na(z) & grepl(pattern, z)] <- category
  }
  annotated <- !is.na(z) & grepl("^(scr|cr|vgpr|vpgr|pr|mr|sd|pd)\\s*(\\(|-)",z) & !grepl("compatible|unconfirmed|possible|unclear",z)
  lead <- toupper(sub("^([a-z]+).*$","\\1",z[annotated]))
  lead[lead=="VPGR"] <- "VGPR"; lead[lead=="SCR"] <- "sCR"
  out[annotated] <- lead
  out
}

tf_incidence <- function(data, scenario, group = "product") {
  groups <- c("Overall", unname(tf_product_labels))
  dplyr::bind_rows(lapply(groups, function(g) {
    z <- if (g == "Overall") data else data[data$product == g, , drop = FALSE]
    n <- nrow(z); events <- sum(z$analysis_flare, na.rm = TRUE)
    ci <- if (n) stats::binom.test(events, n)$conf.int else c(NA_real_, NA_real_)
    tibble::tibble(scenario = scenario, product = g, courses = n, flare_candidates = events,
      incidence_percent = if (n) 100 * events/n else NA_real_, lower_95_percent = 100 * ci[1], upper_95_percent = 100 * ci[2])
  }))
}

tf_summary_fields <- function(data, specs, groups = c("All", "Flare candidate", "No included flare")) {
  result <- list()
  for (group in groups) {
    z <- if (group == "All") data else data[data$analysis_flare == (group == "Flare candidate"), , drop = FALSE]
    for (field in names(specs)) {
      value <- z[[field]]; label <- specs[[field]]
      if (is.numeric(value) && !inherits(value, "Date")) {
        known <- value[!is.na(value) & is.finite(value)]
        summary <- if (length(known)) sprintf("%.1f (%.1f, %.1f)", stats::median(known), stats::quantile(known, .25), stats::quantile(known, .75)) else "Unavailable"
        result[[length(result)+1L]] <- tibble::tibble(group = group, measure = label, category = "Median (IQR)", n = length(known), denominator = nrow(z), percent = NA_real_, summary = summary)
        result[[length(result)+1L]] <- tibble::tibble(group = group, measure = label, category = "Blank", n = sum(is.na(value)), denominator = nrow(z), percent = if(nrow(z)) 100*sum(is.na(value))/nrow(z) else NA_real_, summary = NA_character_)
      } else {
        counts <- table(tf_category(value))
        result[[length(result)+1L]] <- tibble::tibble(group = group, measure = label, category = names(counts), n = as.integer(counts), denominator = nrow(z), percent = if(nrow(z)) 100*as.integer(counts)/nrow(z) else NA_real_, summary = NA_character_)
      }
    }
  }
  dplyr::bind_rows(result)
}

tf_baseline_data <- function(x) x
tf_outcome_data <- function(x) x
