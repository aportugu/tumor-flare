#!/usr/bin/env Rscript
# Publication figures from the canonical, reconciled tumor-flare analysis.
# This file uses relative paths and contains no patient identifiers or dates.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
})
args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1L) args[1] else "analysis"
output_dir <- if (length(args) >= 2L) args[2] else "publication"
framework_only <- "--framework-only" %in% args
dir.create(file.path(output_dir, "Figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "Display_data"), recursive = TRUE, showWarnings = FALSE)
font_family <- "Helvetica"
colours <- c("Yes" = "#1B5678", "No" = "#D6DCE0", "Unknown" = "#BB9153", "Not applicable" = "#A9A6B8")
theme_paper <- function(size = 10) {
  theme_classic(base_size = size, base_family = font_family) +
    theme(plot.title = element_text(face = "bold", size = size + 1, hjust = 0),
          plot.subtitle = element_text(size = size - 1, lineheight = 1.1),
          axis.title = element_text(size = size), axis.text = element_text(colour = "#222222"),
          legend.position = "bottom", legend.title = element_blank(),
          plot.margin = margin(8, 10, 8, 8))
}
read_input <- function(name, required = character()) {
  path <- file.path(input_dir, paste0(name, ".csv"))
  if (!file.exists(path)) stop("Missing canonical output: ", path)
  x <- read_csv(path, show_col_types = FALSE, progress = FALSE)
  x <- x %>% mutate(across(where(is.character), ~gsub("\u2265", ">=", gsub("\u2013|\u2014|\u2212", "-", .x))))
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(name, " lacks required columns: ", paste(missing, collapse = ", "))
  if (name == "case_coded") {
    contract <- jsonlite::fromJSON(file.path(input_dir, "public_case_vocabulary.json"))
    if (!setequal(names(x), contract$fields)) stop("Public case fields differ from closed contract")
    for (field in names(contract$finite_category_vocabulary)) {
      if (!all(as.character(x[[field]]) %in% contract$finite_category_vocabulary[[field]])) stop("Case value outside finite vocabulary: ", field)
    }
    for (field in contract$numeric_elapsed_fields) {
      value <- suppressWarnings(as.numeric(x[[field]]))
      if (any(!is.na(x[[field]]) & (!is.finite(value) | value < 0 | value > 36600))) stop("Invalid numeric elapsed value: ", field)
    }
    if (anyDuplicated(x$case_label) || !all(grepl("^TF[0-9]{2}$", x$case_label))) stop("Invalid neutral case labels")
  }
  x
}
save_display_data <- function(x, name) {
  stopifnot(!any(grepl("(^|_)(course_id|patient_id|mrn|source_reference|narrative|calendar_date|start_date)$", names(x))))
  write_csv(x, file.path(output_dir, "Display_data", paste0(name, ".csv")), na = "")
}
save_figure <- function(p, name, width, height) {
  pdf_path <- file.path(output_dir, "Figures", paste0(name, ".pdf"))
  png_path <- file.path(output_dir, "Figures", paste0(name, ".png"))
  ggsave(pdf_path, p, device = grDevices::pdf, width = width, height = height,
         units = "in", family = font_family, useDingbats = FALSE, bg = "white")
  ggsave(png_path, p, device = ragg::agg_png, width = width, height = height,
         units = "in", dpi = 320, bg = "white")
}
panel_title <- function(letter, title) paste0(letter, "  ", title)
box_layer <- function(x, y, w, h, label, fill = "#EEF2F4", size = 3.4) {
  list(annotate("rect", xmin = x - w/2, xmax = x + w/2, ymin = y - h/2,
                ymax = y + h/2, fill = fill, colour = "#34424B", linewidth = 0.45),
       annotate("text", x = x, y = y, label = label, family = font_family,
                size = size, lineheight = 1.1))
}
arrow_layer <- function(x, y, xend, yend, colour = "#34424B") {
  annotate("segment", x = x, y = y, xend = xend, yend = yend, colour = colour,
           linewidth = .55, arrow = grid::arrow(length = grid::unit(2.5, "mm")))
}

# DISPLAY Figure5 START
recognition_framework <- function() {
  ggplot() +
    box_layer(4.7, 9.2, 8.0, 1.05,
              "New focal pain, swelling, or apparent lesion worsening\nafter bispecific antibody initiation") +
    arrow_layer(4.7, 8.67, 4.7, 8.18) +
    box_layer(4.7, 7.55, 8.0, 1.25,
              "Assess immediately for urgent complications\nNeurological compromise, cord compression, fracture,\norgan dysfunction, or hemodynamic instability", "#FAECEA") +
    arrow_layer(4.7, 6.92, 4.7, 6.42) +
    box_layer(4.7, 5.73, 8.0, 1.35,
              "Evaluate competing explanations\nTrue progression, infection, CRS or ICANS, structural injury,\nand other treatment-related or unrelated causes") +
    arrow_layer(4.7, 5.05, 4.7, 4.52) +
    box_layer(4.7, 3.80, 8.0, 1.4,
              "Integrate evidence and individualize treatment\nClinical course, biochemical trajectory, imaging, and biopsy when feasible\nAddress complications and symptoms; reassess the risks and benefits\nof continuing, interrupting, or changing therapy") +
    arrow_layer(4.7, 3.09, 4.7, 2.55) +
    box_layer(4.7, 1.78, 8.0, 1.5,
              "Reassess longitudinally\nDocument symptom recovery, local lesion outcome, and systemic response\nseparately; record treatment changes and unresolved attribution\nRevise the working diagnosis as evidence develops", "#EAF2F5") +
    annotate("text", x = 4.7, y = .38,
             label = "Proposed clinical framework; diagnostic criteria and treatment decisions require prospective validation.",
             family = font_family, size = 3.1, colour = "#333333") +
    coord_cartesian(xlim = c(.2, 9.2), ylim = c(0, 10), expand = FALSE) +
    theme_void(base_family = font_family) + theme(plot.margin = margin(12, 10, 10, 10))
}
save_figure(recognition_framework(), "Figure_5_Recognition_framework", 8.5, 9.0)
save_display_data(data.frame(step = 1:5,
  action = c("Recognize new focal symptoms or apparent lesion worsening", "Assess urgent complications",
             "Evaluate competing explanations", "Integrate evidence and individualize treatment",
             "Reassess longitudinally"), validation = "Proposed framework"), "Figure_5")
# DISPLAY Figure5 END
if (framework_only) quit(save = "no", status = 0L)

# DISPLAY Figure1 START
flow_panel <- function(x, title) {
  stopifnot(all(c("label", "n") %in% names(x)))
  x <- x %>% mutate(y = rev(seq_len(n())), text = paste0(label, "\n", format(n, big.mark = ",")))
  p <- ggplot() + theme_void(base_family = font_family) +
    ggtitle(title) + theme(plot.title = element_text(face = "bold", size = 11, hjust = 0),
                          plot.margin = margin(8, 8, 8, 8))
  for (i in seq_len(nrow(x))) {
    p <- p + box_layer(0, x$y[i], 4.6, .72, x$text[i], size = 3.1)
    if (i < nrow(x)) p <- p + arrow_layer(0, x$y[i] - .36, 0, x$y[i+1] + .38)
  }
  p + coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(.45, max(x$y) + .6), expand = FALSE)
}
comparative_flow_panel <- function(x) {
  if (nrow(x) < 6L) stop("Comparative flow requires four sequential stages and two terminal groups")
  stages <- x[seq_len(nrow(x) - 2L), ]
  groups <- x[(nrow(x)-1L):nrow(x), ]
  stages <- stages %>% mutate(y = rev(seq_len(n())) + 1.5)
  p <- ggplot() + theme_void(base_family = font_family) +
    ggtitle(panel_title("A", "Comparative cohort")) +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0), plot.margin = margin(8, 8, 8, 8))
  for (i in seq_len(nrow(stages))) {
    p <- p + box_layer(0, stages$y[i], 4.6, .72,
                      paste0(stages$label[i], "\n", format(stages$n[i], big.mark = ",")), size = 3.1)
    if (i < nrow(stages)) p <- p + arrow_layer(0, stages$y[i] - .36, 0, stages$y[i+1] + .38)
  }
  for (i in 1:2) {
    xx <- if (i == 1L) -1.3 else 1.3
    p <- p + arrow_layer(0, min(stages$y) - .36, xx, 1.42) +
      box_layer(xx, 1.0, 2.35, .82, paste0(groups$label[i], "\n", format(groups$n[i], big.mark = ",")), size = 2.9)
  }
  p + coord_cartesian(xlim = c(-2.6, 2.6), ylim = c(.4, max(stages$y) + .6), expand = FALSE)
}
case_series_flow_panel <- function(x) {
  if (nrow(x) != 7L) stop("Case-series flow requires reviewed total, three confidence groups, and three eligibility groups")
  confidence_text <- paste0(x$label[2:4], ": ", format(x$n[2:4], big.mark = ","), collapse = "\n")
  eligibility_text <- paste0(x$label[5:7], ": ", format(x$n[5:7], big.mark = ","), collapse = "\n")
  ggplot() +
    box_layer(0, 5.8, 4.6, .85, paste0(x$label[1], "\n", format(x$n[1], big.mark = ",")), size = 3.2) +
    arrow_layer(0, 5.37, 0, 4.92) +
    box_layer(0, 4.2, 4.6, 1.45, paste0("Investigator confidence\n", confidence_text), size = 3.2) +
    arrow_layer(0, 3.47, 0, 2.98) +
    box_layer(0, 2.2, 4.6, 1.55, paste0("Comparative eligibility\n", eligibility_text), size = 3.2) +
    annotate("text", x = 0, y = .85, label = "All site-reviewed courses contribute to the descriptive case series.",
             family = font_family, size = 2.65) +
    coord_cartesian(xlim = c(-2.6, 2.6), ylim = c(.4, 6.6), expand = FALSE) +
    theme_void(base_family = font_family) + ggtitle(panel_title("B", "Site-reviewed case series")) +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0), plot.margin = margin(8, 8, 8, 8))
}
flow <- read_input("flow")
incidence <- read_input("incidence", c("scenario", "product", "courses", "flare_candidates", "incidence_percent", "lower_95_percent", "upper_95_percent"))
incidence_label <- function(scenario, title) {
  x <- incidence %>% filter(.data$scenario == .env$scenario, product == "Overall")
  if (nrow(x) != 1L) stop("Expected one overall incidence row for ", scenario)
  paste0(title, ": ", x$flare_candidates, "/", format(x$courses, big.mark = ","),
         " (", sprintf("%.2f", x$incidence_percent), "%; 95% CI ", sprintf("%.2f", x$lower_95_percent), "-", sprintf("%.2f", x$upper_95_percent), "%)")
}
selection_notes <- ""
if (all(c("stage", "courses") %in% names(flow))) {
  count_stage <- function(stage) {
    value <- flow$courses[flow$stage == stage]
    if (length(value) != 1L) stop("Missing or duplicated cohort stage: ", stage)
    value
  }
  cases_flow <- read_input("case_coded", c("confidence"))
  confidence_counts <- table(factor(cases_flow$confidence, levels = c("Definite", "Probable", "Possible")))
  unknown <- count_stage("Excluded: flare status unavailable")
  selection_notes <- paste0("Selection exclusions: ",
    count_stage("Supported agents") - count_stage("First eligible BCMA / all GPRC5D courses"), " later BCMA courses; ",
    count_stage("First eligible BCMA / all GPRC5D courses") - count_stage("Nonbridging eligible courses"), " bridging courses; ",
    unknown, " with unavailable flare status.")
  flow <- bind_rows(
    tibble(track = "comparative", label = c("Master bispecific courses",
      "First eligible BCMA / all GPRC5D courses", "Nonbridging eligible courses",
      paste0("Known flare status\n(", format(unknown, big.mark = ","), " unavailable excluded)"),
      "Suspected flare", "No included flare"),
      n = c(count_stage("Master courses"), count_stage("First eligible BCMA / all GPRC5D courses"),
      count_stage("Nonbridging eligible courses"), count_stage("Primary known-status courses"),
      count_stage("Primary suspected flare"), count_stage("Primary comparison"))),
    tibble(track = "case_series", label = c("Site-reviewed suspected flare", "Definite", "Probable", "Possible",
      "Primary cohort", "Bridging excluded", "Later BCMA excluded"),
      n = c(count_stage("Clinical returned series"), as.numeric(confidence_counts),
      count_stage("Primary suspected flare"), count_stage("Clinical bridging cases"), count_stage("Clinical later BCMA cases"))))
}
if (!all(c("track", "label", "n") %in% names(flow))) stop("flow needs canonical stage/courses or adapted track/label/n")
fa <- flow %>% filter(track == "comparative")
fb <- flow %>% filter(track == "case_series")
stopifnot(nrow(fa) > 0, nrow(fb) > 0)
save_figure((comparative_flow_panel(fa) + case_series_flow_panel(fb)) +
  plot_annotation(caption = paste(incidence_label("primary", "Primary documented proportion"),
                                  incidence_label("strict", "Definite/Probable-only sensitivity"), selection_notes, sep = "\n"),
    theme = theme(plot.caption = element_text(family = font_family, size = 9, hjust = 0, lineheight = 1.2))),
            "Figure_1_Cohort_flow", 10, 7.8)
save_display_data(bind_rows(flow, incidence %>% mutate(track = "incidence", label = scenario, n = flare_candidates)), "Figure_1")
# DISPLAY Figure1 END

# DISPLAY Figure2 START
timeline <- read_input("case_coded", c("case_label", "onset_days", "peak_days", "resolution_days", "confidence", "bridging"))
# Calendar dates, original identifiers, and free text never enter figure exports.
timeline <- timeline %>% select(any_of(c("case_label", "product", "confidence", "bridging", "eligibility",
      "onset_days", "onset_lower", "onset_upper", "peak_days", "resolution_days",
  "onset_estimated", "peak_estimated", "resolution_estimated")))
for (f in c("onset_estimated", "peak_estimated", "resolution_estimated")) if (!f %in% names(timeline)) timeline[[f]] <- FALSE
for (f in c("onset_lower", "onset_upper")) if (!f %in% names(timeline)) timeline[[f]] <- NA_real_
timeline <- timeline %>% arrange(is.na(onset_days), onset_days, case_label) %>% mutate(row = rev(row_number()))
onset <- timeline %>% filter(is.finite(onset_days))
onset_ecdf <- onset %>% count(onset_days) %>% arrange(onset_days) %>% mutate(cumulative = cumsum(n)/sum(n))
origin_ecdf <- bind_rows(tibble(onset_days = 0, n = 0L, cumulative = 0), onset_ecdf)
ecdf <- ggplot(origin_ecdf, aes(onset_days, cumulative)) + geom_step(colour = "#1B5678", linewidth = .75) +
  geom_point(data = onset_ecdf, size = 1.4, colour = "#1B5678") +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, .02))) +
  scale_x_continuous(breaks = c(0, 7, 14, 30, 60, 90), limits = c(0, 100)) +
  labs(title = panel_title("A", "Recorded onset"), subtitle = paste(nrow(onset), "cases with point onset timing"),
       x = "Days after first step-up dose", y = "Cumulative proportion") + theme_paper()
markers <- timeline %>% pivot_longer(c(onset_days, peak_days, resolution_days), names_to = "event", values_to = "day") %>%
  filter(is.finite(day)) %>% mutate(event = recode(event, onset_days = "Onset", peak_days = "Recorded peak", resolution_days = "Symptom resolution"),
    estimated = case_when(event == "Onset" ~ onset_estimated, event == "Recorded peak" ~ peak_estimated, TRUE ~ resolution_estimated),
    draw_day = pmin(day, 100))
segments <- timeline %>% filter(is.finite(onset_days), is.finite(resolution_days))
intervals <- timeline %>% filter(is.na(onset_days), is.finite(onset_lower), is.finite(onset_upper))
tp <- ggplot() +
  geom_segment(data = segments, aes(x = onset_days, xend = pmin(resolution_days, 100), y = row, yend = row), colour = "#BCC5CB", linewidth = .45) +
  geom_segment(data = intervals, aes(x = onset_lower, xend = onset_upper, y = row, yend = row), colour = "#1B5678", linewidth = 1.2) +
  geom_point(data = markers, aes(draw_day, row, shape = event, colour = event, fill = ifelse(estimated, "Estimated", "Recorded")), size = 2.35, stroke = .65) +
  geom_text(data = markers %>% filter(day > 100), aes(x = 101, y = row), label = ">", size = 2.8, family = font_family) +
  geom_text(data = timeline, aes(x = 116, y = row, label = confidence), size = 2.55, family = font_family) +
  geom_text(data = timeline, aes(x = 139, y = row, label = ifelse(bridging, "Yes", "No")), size = 2.55, family = font_family) +
  annotate("text", x = 116, y = max(timeline$row) + 1, label = "Confidence", family = font_family, fontface = "bold", size = 2.8) +
  annotate("text", x = 139, y = max(timeline$row) + 1, label = "Bridging", family = font_family, fontface = "bold", size = 2.8) +
  geom_vline(xintercept = 104, colour = "#D3D8DB", linewidth = .4) +
  scale_shape_manual(values = c("Onset" = 21, "Recorded peak" = 24, "Symptom resolution" = 22)) +
  scale_colour_manual(values = c("Onset" = "#1B5678", "Recorded peak" = "#9B6A2D", "Symptom resolution" = "#547761")) +
  scale_fill_manual(values = c("Recorded" = "#687E8B", "Estimated" = "white")) +
  scale_y_continuous(breaks = timeline$row, labels = timeline$case_label, expand = expansion(add = c(.5, 1.6))) +
  scale_x_continuous(breaks = c(0, 7, 14, 30, 60, 90, 100)) +
  coord_cartesian(xlim = c(0, 147), clip = "off") +
  labs(title = panel_title("B", "Recorded chronology"), x = "Days after first step-up dose", y = NULL,
       subtitle = "Open markers indicate estimated dates; horizontal onset marks indicate intervals") +
  theme_paper(9) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(size = 7.2)) +
  guides(fill = guide_legend(order = 2, override.aes = list(shape = 21, colour = "#34424B", fill = c("white", "#687E8B"), size = 3)),
         shape = guide_legend(order = 1), colour = guide_legend(order = 1))
save_figure(ecdf / tp + plot_layout(heights = c(1, 3.3)), "Figure_2_Onset_and_chronology", 10.5, 11)
save_display_data(timeline %>% select(-row), "Figure_2")
# DISPLAY Figure2 END

# DISPLAY Figure3 START
status_bars <- function(x, title) {
  stopifnot(all(c("label", "category", "n", "denominator") %in% names(x)))
  x <- x %>% mutate(category = ifelse(is.na(category) | category == "", "Unknown", category))
  if (!all(x$category %in% names(colours))) stop("Unexpected clinical status category")
  labels <- unique(x$label)
  x <- x %>% mutate(label = factor(label, levels = rev(labels)), category = factor(category, levels = names(colours)), proportion = n/denominator)
  yes <- x %>% filter(category == "Yes")
  ggplot(x, aes(proportion, label, fill = category)) + geom_col(width = .64, position = position_stack(reverse = TRUE), colour = "white", linewidth = .25, show.legend = TRUE) +
    geom_text(data = yes, aes(x = 1.02, y = label, label = paste0(n, "/", denominator)), inherit.aes = FALSE,
              hjust = 0, size = 3.0, family = font_family) +
    scale_fill_manual(values = colours, drop = FALSE) +
    scale_x_continuous(labels = scales::label_percent(accuracy = 1), breaks = c(0, .25, .5, .75, 1), expand = expansion(mult = c(0, 0))) +
    coord_cartesian(xlim = c(0, 1.14), clip = "off") +
    labs(title = title, x = "Proportion of site-reviewed cases", y = NULL) + theme_paper(9.5) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
}
presentation <- read_input("presentation", c("field", "label", "panel", "category", "n", "denominator"))
anatomy <- read_input("anatomy_binary", c("field", "label", "category", "n", "denominator"))
manifest <- presentation %>% filter(field %in% c("pain_flare", "palpable_enlargement", "imaging_worsening", "new_lesion", "neurologic_cns", "biochemical_or_organ"))
impact <- presentation %>% filter(field %in% c("organ_dysfunction", "flare_hospitalization", "neurologic_deficit", "pathologic_fracture", "cord_compression", "icu_due_to_flare"))
site <- anatomy %>% filter(field %in% c("true_emd_any", "paraskeletal_involvement", "spinal_involvement", "cns_involvement"))
save_figure((status_bars(manifest, panel_title("A", "Clinical manifestations")) /
             status_bars(impact, panel_title("B", "Clinical impact")) /
             status_bars(site, panel_title("C", "Index anatomical involvement"))) + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
            "Figure_3_Clinical_spectrum", 9.5, 10)
save_display_data(bind_rows(manifest, impact, site), "Figure_3")
# DISPLAY Figure3 END

# DISPLAY Figure4 START
forest <- function(x, title, metric) {
  stopifnot(all(c("label", "estimate", "conf_low", "conf_high", "n_courses", "candidates") %in% names(x)))
  x <- x %>% mutate(row = rev(row_number()))
  rng <- range(c(x$conf_low, x$conf_high, 1), na.rm = TRUE)
  rng <- c(max(.01, rng[1] * .8), rng[2] * 1.2)
  labels <- ggplot(x, aes(y = row)) + geom_text(aes(x = 0, label = label), hjust = 0, family = font_family, size = 3.0) +
    coord_cartesian(xlim = c(0, 1), ylim = c(.35, nrow(x) + .65), clip = "off") + theme_void() + theme(plot.margin = margin(4, 2, 4, 2))
  values <- x %>% mutate(value = ifelse(is.finite(estimate), sprintf("%.2f (%.2f-%.2f)", estimate, conf_low, conf_high), "Not estimated"),
    population = paste0("N=", format(n_courses, big.mark = ","), "; flare=", candidates))
  right <- ggplot(values, aes(y = row)) + geom_text(aes(x = 0, label = value), hjust = 0, family = font_family, size = 3.0) +
    geom_text(aes(x = .60, label = population), hjust = 0, family = font_family, size = 2.75) +
    coord_cartesian(xlim = c(0, 1.35), ylim = c(.35, nrow(x) + .65), clip = "off") + theme_void() + theme(plot.margin = margin(4, 2, 4, 2))
  center <- ggplot(x, aes(estimate, row)) + geom_vline(xintercept = 1, linetype = "dashed", colour = "#89949C", linewidth = .4) +
    geom_errorbar(aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = .12, linewidth = .55, na.rm = TRUE) +
    geom_point(size = 2.1, colour = "#1B5678", na.rm = TRUE) + scale_x_log10(limits = rng) +
    coord_cartesian(ylim = c(.35, nrow(x) + .65)) + labs(x = paste0(metric, " (95% CI)"), y = NULL) + theme_paper(8.5) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          plot.margin = margin(4, 4, 4, 4))
  (labels + center + right + plot_layout(widths = c(1.05, 1.0, 1.55))) +
    plot_annotation(title = title, theme = theme(plot.title = element_text(family = font_family, face = "bold", size = 11)))
}
models <- read_input("models", c("scenario", "family", "label", "estimate", "conf_low", "conf_high", "n_courses", "candidates"))
survival <- read_input("survival_models", c("scenario", "model", "endpoint", "estimate", "conf_low", "conf_high", "n_courses", "exposed_courses"))
baseline <- models %>% filter(scenario == "primary", family == "Baseline")
responses <- models %>% filter(scenario == "primary", family == "Response", !is.na(adjustment), adjustment != "")
surv_primary <- survival %>% filter(scenario == "primary", grepl("Adjusted", model)) %>%
  mutate(label = endpoint, candidates = exposed_courses)
if (nrow(baseline) != 6L || nrow(responses) != 3L || nrow(surv_primary) != 2L) stop("Figure4 model selection is not 6 baseline/3 response/2 survival models")
save_figure(wrap_plots(list(
              wrap_elements(full = forest(baseline, panel_title("A", "Baseline associations with documented suspected flare"), "OR")),
              wrap_elements(full = forest(responses, panel_title("B", "Adjusted best-response associations"), "OR")),
              wrap_elements(full = forest(surv_primary, panel_title("C", "Adjusted time-dependent survival associations"), "HR"))),
              ncol = 1, heights = c(6, 3.3, 2.6)), "Figure_4_Association_models", 12, 9)
save_display_data(bind_rows(baseline, responses, surv_primary) %>% select(any_of(c("scenario", "family", "model", "label", "endpoint", "estimate", "conf_low", "conf_high", "p_value", "n_courses", "n_patient_clusters", "candidates", "events", "adjustment", "bootstrap_successes", "bootstrap_requested", "model_status"))), "Figure_4")
# DISPLAY Figure4 END

# DISPLAY FigureS1 START
km <- read_input("landmark_km_steps", c("endpoint", "landmark_day", "group", "time", "survival", "lower", "upper"))
risk <- read_input("landmark_risk_counts", c("endpoint", "landmark_day", "group", "time", "n_risk"))
group_label <- function(x) ifelse(as.character(x) == "0", "No flare by landmark", ifelse(as.character(x) == "1", "Flare by landmark", as.character(x)))
km <- km %>% mutate(group = group_label(group))
risk <- risk %>% mutate(group = group_label(group))
group_colours <- c("No flare by landmark" = "#65747E", "Flare by landmark" = "#1B5678")
landmark_panel <- function(endpoint_value, day_value, letter) {
  x <- km %>% filter(endpoint == endpoint_value, landmark_day == day_value)
  r <- risk %>% filter(endpoint == endpoint_value, landmark_day == day_value)
  stopifnot(nrow(x) > 0, nrow(r) > 0)
  xmax <- max(r$time)
  x <- x %>% arrange(group, time)
  step_band <- x %>% group_by(group) %>% mutate(next_time = lead(time)) %>% ungroup()
  step_band <- bind_rows(step_band %>% mutate(order = 1), step_band %>% filter(is.finite(next_time)) %>% mutate(time = next_time, order = 0)) %>% arrange(group, time, order)
  # Canonical landmark KM and risk-count time fields are elapsed months.
  p <- ggplot(x, aes(time, survival, colour = group, fill = group)) +
    geom_ribbon(data = step_band, aes(ymin = lower, ymax = upper), alpha = .12, colour = NA) + geom_step(linewidth = .65) +
    scale_colour_manual(values = group_colours) + scale_fill_manual(values = group_colours) +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1)) +
    scale_x_continuous(breaks = sort(unique(r$time)), labels = scales::label_number(accuracy = 1)) + coord_cartesian(xlim = c(0, xmax)) +
    labs(title = panel_title(letter, paste0(endpoint_value, " at day ", day_value)), x = "Months after landmark", y = paste0(endpoint_value, " probability")) +
    theme_paper(9) + theme(legend.position = "bottom")
  groups <- c("Flare by landmark", "No flare by landmark")
  r <- r %>% mutate(row = match(group, rev(groups)))
  rt <- ggplot(r, aes(time, row, label = n_risk, colour = group)) + geom_text(family = font_family, size = 2.7) +
    scale_colour_manual(values = group_colours) + scale_x_continuous(breaks = sort(unique(r$time))) +
    scale_y_continuous(breaks = 1:2, labels = rev(groups)) + coord_cartesian(xlim = c(0, xmax), ylim = c(.4, 2.6), clip = "off") +
    labs(x = NULL, y = NULL, title = "Number at risk") + theme_void(base_family = font_family) +
    theme(axis.text.y = element_text(size = 7, hjust = 1), plot.title = element_text(size = 8.5, hjust = 0),
          plot.margin = margin(3, 10, 3, 8), legend.position = "none")
  p / rt + plot_layout(heights = c(4.2, 1))
}
save_figure(wrap_plots(list(
              wrap_elements(full = landmark_panel("PFS", 30, "A")),
              wrap_elements(full = landmark_panel("OS", 30, "B")),
              wrap_elements(full = landmark_panel("PFS", 60, "C")),
              wrap_elements(full = landmark_panel("OS", 60, "D"))), ncol = 2, nrow = 2),
            "Figure_S1_Landmark_survival", 12, 9.5)
save_display_data(km, "Figure_S1_KM")
save_display_data(risk, "Figure_S1_Risk")
# DISPLAY FigureS1 END
cat("Publication figures exported from", normalizePath(input_dir), "\n")
