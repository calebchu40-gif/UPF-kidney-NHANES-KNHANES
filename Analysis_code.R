upf_source <- list()

upf_source[["NHANES/02_代码/正式五结局/00_config.R"]] <- r"--------------------(# NHANES UPF five-outcome formal continuous-analysis configuration.

PHASE1_VERSION <- "2026-09-07.cholesterol-design-revision"
PROJECT_ROOT_DEFAULT <- "/storage/home/tmu2301/nhanes分析/18_formal_cholesterol_design_reanalysis_20260907/NHANES"
PROJECT_ROOT <- Sys.getenv("NHANES_PHASE1_ROOT", unset = PROJECT_ROOT_DEFAULT)
SOURCE_REV22_ROOT <- Sys.getenv(
  "NHANES_REV22_SOURCE_ROOT",
  unset = "/storage/home/tmu2301/nhanes分析/UPF_全肾结局_REV22_两日膳食_最终分析_20260813"
)
RAW_ROOT <- Sys.getenv(
  "NHANES_RAW_ROOT",
  unset = "/storage/home/tmu2301/nhanes分析/数据/nhanes数据"
)
LEGACY_2007_2020_ROOT <- Sys.getenv(
  "NHANES_LEGACY_ROOT",
  unset = "/storage/home/tmu2301/nhanes分析/UPF_全肾结局_窗口扩展_20260806_2007_2020"
)

DIRS <- list(
  plan = file.path(PROJECT_ROOT, "00_计划与官方依据"),
  config = file.path(PROJECT_ROOT, "01_数据配置"),
  code = file.path(PROJECT_ROOT, "02_代码", "正式五结局"),
  derived = file.path(PROJECT_ROOT, "03_派生数据", "正式五结局_v1"),
  results = file.path(PROJECT_ROOT, "04_结果", "正式五结局_v1"),
  tables = file.path(PROJECT_ROOT, "05_表格", "正式五结局_v1"),
  figures = file.path(PROJECT_ROOT, "06_图形", "正式五结局_v1"),
  logs = file.path(PROJECT_ROOT, "07_运行记录", "正式五结局_v1"),
  tests = file.path(PROJECT_ROOT, "08_测试", "正式五结局_v1"),
  archive = file.path(PROJECT_ROOT, "09_归档", "正式五结局_v1")
)

FORMAL_AUTH_TOKEN <- "I_AUTHORIZE_NHANES_PHASE1_FORMAL_RUN"
formal_analysis_authorized <- function() {
  identical(Sys.getenv("NHANES_PHASE1_FORMAL_AUTH", unset = ""), FORMAL_AUTH_TOKEN) &&
    file.exists(file.path(DIRS$plan, "FORMAL_ANALYSIS_AUTHORIZED.txt"))
}
require_formal_authorization <- function(stage) {
  if (!formal_analysis_authorized()) {
    stop(
      "Phase-1 formal-analysis lock is active for stage '", stage,
      "'. The authorization token and marker are both required.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)
PHASE1_SEED <- 20260824L
MI_M <- 50L
MI_MAXIT <- 60L
MAX_PARALLEL_WORKERS <- 14L
DEFAULT_N_CORES <- 14L

SAFE_MAX_USER_PROCESSES <- 28L
SAFE_WARN_RSS_GIB <- 250
SAFE_STOP_RSS_GIB <- 300
SAFE_MIN_AVAILABLE_MEMORY_GIB <- 400
SAFE_START_AVAILABLE_MEMORY_GIB <- 600

for (.thread_env in c(
  "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
  "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"
)) {
  do.call(Sys.setenv, stats::setNames(list("1"), .thread_env))
}
rm(.thread_env)

n_cores_resolve <- function() {
  x <- suppressWarnings(as.integer(Sys.getenv(
    "NHANES_PHASE1_WORKERS", unset = as.character(DEFAULT_N_CORES))))
  if (!is.finite(x) || x < 1L) x <- DEFAULT_N_CORES
  as.integer(min(x, MAX_PARALLEL_WORKERS))
}

CYCLES <- data.frame(
  cycle = c("E", "F", "G", "H", "I", "J", "P"),
  label = c("2007-2008", "2009-2010", "2011-2012", "2013-2014",
            "2015-2016", "2017-2018", "2017-March 2020"),
  suffix = c("e", "f", "g", "h", "i", "j", "p"),
  duration_years = c(2, 2, 2, 2, 2, 0, 3.2),
  include = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
  stringsAsFactors = FALSE
)
ANALYTIC_CYCLES <- CYCLES$cycle[CYCLES$include]
TOTAL_SURVEY_YEARS <- sum(CYCLES$duration_years[CYCLES$include])

cycle_dir <- c(
  E = "2007-2008", F = "2009-2010", G = "2011-2012",
  H = "2013-2014", I = "2015-2016", J = "2017-2018",
  P = "2019-2020"
)

module_file <- function(cycle, module, stem) {
  stopifnot(cycle %in% ANALYTIC_CYCLES)
  if (cycle == "P") {
    candidates <- c(
      file.path(RAW_ROOT, cycle_dir[[cycle]], module,
                paste0("p_", tolower(stem), ".xpt")),
      file.path(RAW_ROOT, cycle_dir[[cycle]], module,
                paste0("P_", toupper(stem), ".XPT"))
    )
  } else {
    candidates <- c(file.path(
      RAW_ROOT, cycle_dir[[cycle]], module,
      paste0(stem, "_", tolower(cycle), ".xpt")
    ))
  }
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) normalizePath(hit[[1]], winslash = "/", mustWork = TRUE) else candidates[[1]]
}

diet_file <- function(cycle, day, kind = c("iff", "tot")) {
  stopifnot(cycle %in% ANALYTIC_CYCLES, day %in% 1:2)
  kind <- match.arg(kind)
  if (cycle == "P") {
    stem <- sprintf("P_DR%d%s.XPT", day, toupper(kind))
    candidates <- c(
      file.path(Sys.getenv("NHANES_OFFICIAL_XPT_DIR"), stem),
      file.path(DIRS$config, "official_xpt", stem),
      file.path(DIRS$config, "official_xpt", sub("\\.XPT$", ".xpt", stem)),
      file.path(RAW_ROOT, cycle_dir[[cycle]], "Dietary", stem),
      file.path(RAW_ROOT, cycle_dir[[cycle]], "Dietary",
                sub("\\.XPT$", ".xpt", stem)),
      file.path(RAW_ROOT, cycle_dir[[cycle]], "Dietary", tolower(stem)),
      file.path(LEGACY_2007_2020_ROOT, "01_数据", "官方XPT_2017_2020", stem)
    )
  } else {
    stem <- sprintf("dr%d%s_%s.xpt", day, kind, tolower(cycle))
    candidates <- c(file.path(RAW_ROOT, cycle_dir[[cycle]], "Dietary", stem))
  }
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) normalizePath(hit[[1]], winslash = "/", mustWork = TRUE) else candidates[[1]]
}

NOVA_MANIFEST_DIR <- file.path(DIRS$config, "nova_manifest")

EXPOSURES <- list(
  G_R0 = list(variable = "upf_gram_ratio_nonwater", label = "Two-day non-water gram share")
)
EXPOSURE_UNIT_PP <- 10
UNCLASSIFIED_GRAM_SHARE_LIMIT <- 0.05

OUTCOMES <- list(
  CKD = list(variable = "ckd", family = "quasibinomial", domain = "kidney_observed",
             effect_scale = "odds_ratio"),
  eGFR = list(variable = "egfr_2021", family = "gaussian", domain = "egfr_observed",
              effect_scale = "identity"),
  UACR = list(variable = "log_uacr", family = "gaussian", domain = "uacr_observed",
              effect_scale = "log_response_percent"),
  DKD = list(variable = "dkd", family = "quasibinomial", domain = "dkd_domain",
             effect_scale = "odds_ratio"),
  kidney_stones = list(variable = "kidney_stones", family = "quasibinomial", domain = "stone_observed",
                       effect_scale = "odds_ratio")
)

AGE_FORMS <- list(linear = c("age10"))

MODEL_COVARIATES <- list(
  M1 = character(),
  M2 = c(
    "age10", "sex", "race_ethnicity", "education", "cycle_factor",
    "mean_energy_kcal"
  ),
  M3 = c(
    "age10", "sex", "race_ethnicity", "education", "cycle_factor",
    "mean_energy_kcal", "smoking", "alcohol", "physical_activity", "bmi"
  ),
  M4 = c(
    "age10", "sex", "race_ethnicity", "education", "cycle_factor",
    "mean_energy_kcal", "smoking", "alcohol", "physical_activity", "bmi",
    "hypertension", "diabetes", "high_cholesterol"
  )
)
MODEL_EXTENSIONS <- MODEL_COVARIATES
PIR_MISSING_LIMIT <- 0.30
PIR_DECISION_FILE <- file.path(DIRS$derived, "pir_decision.json")

PHYSICAL_ACTIVITY_VARS <- c(
  "PAQ605", "PAQ610", "PAD615", "PAQ620", "PAQ625", "PAD630",
  "PAQ635", "PAQ640", "PAD645", "PAQ650", "PAQ655", "PAD660",
  "PAQ665", "PAQ670", "PAD675"
)
PHYSICAL_ACTIVITY_CUTS <- c(0, 600)

required_packages <- c(
  "haven", "data.table", "survey", "mice", "mitools", "digest", "jsonlite"
)

MICE_BASE_METHODS <- c(
  bmi = "pmm", mean_sbp = "pmm", mean_dbp = "pmm",
  total_cholesterol = "pmm", hba1c = "pmm",
  education = "polyreg", smoking = "polyreg", alcohol = "polyreg",
  physical_activity = "polyreg",
  diabetes_diagnosed = "logreg", insulin_current = "logreg",
  oral_med_current = "logreg", hypertension = "logreg",
  high_cholesterol = "logreg"
)

PHASE1_VARIABLE_MAP_FILES <- c(
  file.path(DIRS$config, "variable_map.csv"),
  file.path(DIRS$config, "phase1_variable_extensions.csv")
)

invisible(TRUE)
)--------------------"

upf_source[["NHANES/02_代码/正式五结局/00_functions.R"]] <- r"--------------------(`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

.this_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "NHANES/02_代码/正式五结局/00_functions.R"), error = function(e) NULL) %||% "00_functions.R"
source(file.path(dirname(normalizePath(.this_file, mustWork = FALSE)), "00_config.R"), local = FALSE)

assert_columns <- function(x, cols, label = deparse(substitute(x))) {
  missing <- setdiff(cols, names(x))
  if (length(missing)) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

rbind_fill <- function(xs) {
  cols <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs <- lapply(xs, function(x) {
    for (nm in setdiff(cols, names(x))) x[[nm]] <- NA
    x[cols]
  })
  do.call(rbind, xs)
}

yes_no_na <- function(x) {
  ifelse(x == 1, 1L, ifelse(x == 2, 0L, NA_integer_))
}

or_composite <- function(...) {
  z <- cbind(...)
  positive <- apply(z == 1, 1, any, na.rm = TRUE)
  all_negative <- rowSums(!is.na(z)) == ncol(z) & apply(z == 0, 1, all)
  ifelse(positive, 1L, ifelse(all_negative, 0L, NA_integer_))
}

derive_ckd_three_state <- function(egfr, acr) {
  positive <- (!is.na(egfr) & egfr < 60) | (!is.na(acr) & acr >= 30)
  negative <- !is.na(egfr) & egfr >= 60 & !is.na(acr) & acr < 30
  ifelse(positive, 1L, ifelse(negative, 0L, NA_integer_))
}

calc_egfr_2021 <- function(creatinine_mg_dl, age_years, female) {
  invalid <- is.na(creatinine_mg_dl) | creatinine_mg_dl <= 0 |
    is.na(age_years) | age_years < 18 | is.na(female)
  k <- ifelse(female, 0.7, 0.9)
  alpha <- ifelse(female, -0.241, -0.302)
  out <- 142 * pmin(creatinine_mg_dl / k, 1)^alpha *
    pmax(creatinine_mg_dl / k, 1)^(-1.200) *
    0.9938^age_years * ifelse(female, 1.012, 1)
  out[invalid] <- NA_real_
  out
}

calibrate_creatinine_dxc <- function(x, cycle) {
  ifelse(cycle == "P" & is.finite(x), 1.051 * x - 0.06945, x)
}

calibrate_uric_acid_dxc <- function(x, cycle) {
  ifelse(cycle == "P" & is.finite(x), 0.9323 * x + 0.2326, x)
}

derive_diabetes_components <- function(diagnosis, insulin, oral_med, hba1c) {
  dx <- ifelse(diagnosis == 1, 1L,
               ifelse(diagnosis %in% c(2, 3), 0L, NA_integer_))
  ins <- yes_no_na(insulin)
  oral <- yes_no_na(oral_med)
  ins[dx == 0L & is.na(ins)] <- 0L
  oral[dx == 0L & is.na(oral)] <- 0L
  lab <- ifelse(is.na(hba1c), NA_integer_, as.integer(hba1c >= 6.5))
  list(
    diagnosed = dx,
    insulin = ins,
    oral_med = oral,
    diabetes = or_composite(dx, ins, oral, lab)
  )
}

derive_hypertension_components <- function(mean_sbp, mean_dbp, diagnosis, medication) {
  measured_positive <- (!is.na(mean_sbp) & mean_sbp >= 140) |
    (!is.na(mean_dbp) & mean_dbp >= 90)
  measured_negative <- !is.na(mean_sbp) & mean_sbp < 140 &
    !is.na(mean_dbp) & mean_dbp < 90
  measured <- ifelse(measured_positive, 1L,
                     ifelse(measured_negative, 0L, NA_integer_))
  dx <- yes_no_na(diagnosis)
  med <- yes_no_na(medication)
  med[dx == 0L & is.na(med)] <- 0L
  list(
    diagnosed = dx,
    medication = med,
    hypertension = or_composite(measured, dx, med)
  )
}

derive_cholesterol_components <- function(total_cholesterol, diagnosis,
                                          medication_recommended,
                                          medication_current) {
  measured <- ifelse(
    is.na(total_cholesterol), NA_integer_, as.integer(total_cholesterol >= 200)
  )
  dx <- yes_no_na(diagnosis)
  recommended <- yes_no_na(medication_recommended)
  current <- yes_no_na(medication_current)
  med <- ifelse(recommended == 0L, 0L,
                ifelse(recommended == 1L, current, NA_integer_))
  list(
    diagnosed = dx,
    medication = med,
    high_cholesterol = or_composite(measured, dx, med)
  )
}

row_mean_valid <- function(x) {
  ans <- rowMeans(as.matrix(x), na.rm = TRUE)
  ans[rowSums(!is.na(x)) == 0] <- NA_real_
  ans
}

nhanes_bp_mean_auscultatory <- function(x, diastolic = FALSE) {
  x <- as.matrix(x)
  apply(x, 1, function(v) {
    observed <- which(is.finite(v))
    if (!length(observed)) return(NA_real_)
    keep <- if (length(observed) == 1L) observed else observed[-1L]
    z <- as.numeric(v[keep])
    if (diastolic) z <- z[z != 0]
    if (!length(z)) NA_real_ else mean(z)
  })
}

nhanes_bp_mean_oscillometric <- function(x, diastolic = FALSE) {
  x <- as.matrix(x)
  apply(x, 1, function(v) {
    z <- as.numeric(v[seq_len(min(3L, length(v)))])
    z <- z[is.finite(z)]
    if (diastolic) z <- z[z != 0]
    if (!length(z)) NA_real_ else mean(z)
  })
}

derive_smoking <- function(ever_100, current) {
  factor(
    ifelse(ever_100 == 2, "Never",
           ifelse(ever_100 == 1 & current == 3, "Former",
                  ifelse(ever_100 == 1 & current %in% c(1, 2), "Current", NA_character_))),
    levels = c("Never", "Former", "Current")
  )
}

derive_alcohol_weekly <- function(cycle, sex, alq101 = NA, alq110 = NA,
                                  alq120q = NA, alq120u = NA,
                                  alq111 = NA, alq121 = NA, alq130 = NA) {
  drinks_week <- rep(NA_real_, length(cycle))
  regular <- cycle %in% c("E", "F", "G", "H", "I")
  q <- as.numeric(alq120q)
  u <- as.numeric(alq120u)
  multiplier <- ifelse(u == 1, 1, ifelse(u == 2, 1 / 4.345,
                 ifelse(u == 3, 1 / 52.1775, NA_real_)))
  valid_q <- is.finite(q) & q >= 0 & !q %in% c(777, 999)
  frequency_week <- q * multiplier
  frequency_week[!valid_q] <- NA_real_
  frequency_week[regular & valid_q & q == 0] <- 0
  p_mid <- c(`0` = 0, `1` = 7, `2` = 6, `3` = 3.5, `4` = 2,
             `5` = 1, `6` = 2.5 / 4.345, `7` = 1 / 4.345,
             `8` = 9 / 52.1775, `9` = 4.5 / 52.1775, `10` = 1.5 / 52.1775)
  p <- cycle == "P" & as.character(alq121) %in% names(p_mid)
  frequency_week[p] <- unname(p_mid[as.character(alq121[p])])
  amount <- as.numeric(alq130)
  valid_amount <- is.finite(amount) & amount > 0 & !amount %in% c(777, 999)
  drinks_week[is.finite(frequency_week) & frequency_week == 0] <- 0
  ok <- is.finite(frequency_week) & valid_amount
  drinks_week[ok] <- frequency_week[ok] * amount[ok]
  drinks_week[regular & alq101 == 2] <- 0
  drinks_week[regular & is.na(alq101) & alq110 == 2] <- 0
  drinks_week[cycle == "P" & alq111 == 2] <- 0
  cut <- ifelse(sex == "Female", 7, ifelse(sex == "Male", 14, NA_real_))
  factor(
    ifelse(is.na(drinks_week) | is.na(cut), NA_character_,
           ifelse(drinks_week < 1, "None_or_very_low",
                  ifelse(drinks_week <= cut, "Low_to_moderate", "Higher"))),
    levels = c("None_or_very_low", "Low_to_moderate", "Higher")
  )
}

gpaq_component <- function(yes_no, days, minutes, met) {
  ans <- rep(NA_real_, length(yes_no))
  ans[yes_no == 2] <- 0
  valid <- yes_no == 1 & is.finite(days) & days >= 0 & days <= 7 &
    is.finite(minutes) & minutes >= 0 & minutes <= 1440
  ans[valid] <- met * days[valid] * minutes[valid]
  ans
}

derive_gpaq_met_category <- function(d) {
  components <- cbind(
    gpaq_component(d$PAQ605, d$PAQ610, d$PAD615, 8),
    gpaq_component(d$PAQ620, d$PAQ625, d$PAD630, 4),
    gpaq_component(d$PAQ635, d$PAQ640, d$PAD645, 4),
    gpaq_component(d$PAQ650, d$PAQ655, d$PAD660, 8),
    gpaq_component(d$PAQ665, d$PAQ670, d$PAD675, 4)
  )
  total <- rowSums(components, na.rm = FALSE)
  category <- factor(
    ifelse(is.na(total), NA_character_,
           ifelse(total == 0, "Inactive",
                  ifelse(total < 600, "Insufficient", "Sufficient"))),
    levels = c("Inactive", "Insufficient", "Sufficient")
  )
  list(total_met_min_week = total, category = category)
}

parse_nova_values <- function(x, label = "NOVA manifest") {
  out <- suppressWarnings(as.integer(gsub("[^1-4]", "", as.character(x))))
  if (any(!is.na(x) & !out %in% 1:4)) stop("Invalid NOVA value in ", label, call. = FALSE)
  out
}

derive_analysis_weight <- function(cycle, wtdr2d, wtdr2dpp) {
  out <- rep(NA_real_, length(cycle))
  regular <- cycle %in% c("E", "F", "G", "H", "I")
  out[regular] <- wtdr2d[regular] * 2 / TOTAL_SURVEY_YEARS
  out[cycle == "P"] <- wtdr2dpp[cycle == "P"] * 3.2 / TOTAL_SURVEY_YEARS
  out
}

make_full_design <- function(d) {
  assert_columns(d, c("SDMVPSU", "SDMVSTRA", "analysis_weight", "cycle"))
  keep <- complete.cases(d[c("SDMVPSU", "SDMVSTRA", "analysis_weight")]) &
    is.finite(d$analysis_weight) & d$analysis_weight > 0
  if (!any(keep)) stop("No positive-weight records for survey design", call. = FALSE)
  dd <- d[keep, , drop = FALSE]
  dd$.strata_cycle <- interaction(dd$cycle, dd$SDMVSTRA, drop = TRUE)
  dd$.psu_cycle <- interaction(dd$cycle, dd$SDMVPSU, drop = TRUE)
  survey::svydesign(
    ids = ~.psu_cycle, strata = ~.strata_cycle,
    weights = ~analysis_weight, nest = TRUE, data = dd
  )
}

phase1_model_covariates <- function(model, age_form, outcome, include_pir) {
  stopifnot(model %in% names(MODEL_COVARIATES), age_form %in% names(AGE_FORMS),
            outcome %in% names(OUTCOMES))
  if (isTRUE(include_pir)) stop("PIR is excluded from the formal five-outcome model", call. = FALSE)
  z <- MODEL_COVARIATES[[model]]
  if (outcome == "DKD") z <- setdiff(z, "diabetes")
  unique(z)
}

phase1_formula <- function(exposure, model, age_form, outcome, include_pir) {
  stopifnot(exposure %in% names(EXPOSURES))
  exposure_term <- sprintf("I(%s/%s)", EXPOSURES[[exposure]]$variable, EXPOSURE_UNIT_PP)
  stats::reformulate(
    c(exposure_term, phase1_model_covariates(model, age_form, outcome, include_pir)),
    response = OUTCOMES[[outcome]]$variable
  )
}

outcome_family <- function(outcome) {
  fam <- OUTCOMES[[outcome]]$family
  if (fam == "quasibinomial") stats::quasibinomial() else stats::gaussian()
}

outcome_domain_design <- function(full_design, outcome, exposure) {
  domain <- full_design$variables[[OUTCOMES[[outcome]]$domain]] %in% TRUE &
    is.finite(full_design$variables[[EXPOSURES[[exposure]]$variable]])
  full_design[domain, ]
}

pool_svy_fits <- function(fits) {
  if (length(fits) != MI_M || any(!vapply(fits, inherits, logical(1), "svyglm"))) {
    stop("Pooling requires ", MI_M, "/", MI_M, " successful svyglm fits", call. = FALSE)
  }
  term_sets <- lapply(fits, function(x) names(stats::coef(x)))
  if (!all(vapply(term_sets[-1], identical, logical(1), term_sets[[1]]))) {
    stop("Coefficient sets differ across imputations", call. = FALSE)
  }
  df_complete <- min(vapply(fits, stats::df.residual, numeric(1)))
  mitools::MIcombine(
    lapply(fits, stats::coef), lapply(fits, stats::vcov),
    df.complete = df_complete
  )
}

pool_svy_summaries <- function(summaries) {
  valid <- function(x) {
    is.list(x) && is.numeric(x$coefficients) && is.matrix(x$variance) &&
      is.numeric(x$df_complete) && length(x$df_complete) == 1L
  }
  if (length(summaries) != MI_M || any(!vapply(summaries, valid, logical(1)))) {
    stop("Pooling requires ", MI_M, "/", MI_M,
         " successful numerical survey-model summaries", call. = FALSE)
  }
  term_sets <- lapply(summaries, function(x) names(x$coefficients))
  if (!all(vapply(term_sets[-1], identical, logical(1), term_sets[[1]]))) {
    stop("Coefficient sets differ across imputations", call. = FALSE)
  }
  matrix_names_match <- vapply(summaries, function(x) {
    identical(rownames(x$variance), names(x$coefficients)) &&
      identical(colnames(x$variance), names(x$coefficients))
  }, logical(1))
  if (any(!matrix_names_match)) {
    stop("Covariance-matrix names do not match coefficient names", call. = FALSE)
  }
  if (any(!vapply(summaries, function(x) {
    all(is.finite(x$coefficients)) && all(is.finite(x$variance)) &&
      is.finite(x$df_complete) && x$df_complete > 0
  }, logical(1)))) {
    stop("Non-finite model summary supplied for pooling", call. = FALSE)
  }
  mitools::MIcombine(
    lapply(summaries, `[[`, "coefficients"),
    lapply(summaries, `[[`, "variance"),
    df.complete = min(vapply(summaries, `[[`, numeric(1), "df_complete"))
  )
}

user_resource_snapshot <- function() {
  user <- Sys.getenv("USER", unset = "")
  if (!nzchar(user)) stop("USER is unavailable for resource audit", call. = FALSE)
  lines <- system2("ps", c("-u", user, "-o", "pid=,rss="), stdout = TRUE)
  fields <- strsplit(trimws(lines[nzchar(trimws(lines))]), "[[:space:]]+")
  rss_kib <- vapply(fields, function(x) as.numeric(x[[2]]), numeric(1))
  meminfo <- readLines("/proc/meminfo", warn = FALSE)
  available_line <- grep("^MemAvailable:", meminfo, value = TRUE)
  available_kib <- if (length(available_line)) {
    as.numeric(sub("^MemAvailable:[[:space:]]+([0-9]+).*", "\\1", available_line[[1]]))
  } else {
    NA_real_
  }
  list(
    process_count = length(fields),
    rss_gib = sum(rss_kib) / 1024^2,
    available_memory_gib = available_kib / 1024^2
  )
}

resource_preflight <- function(workers = n_cores_resolve(), stage = "model") {
  snapshot <- user_resource_snapshot()
  projected_processes <- snapshot$process_count + workers
  if (projected_processes > SAFE_MAX_USER_PROCESSES) {
    stop(stage, " preflight rejected: current user processes=",
         snapshot$process_count, ", workers=", workers,
         ", projected=", projected_processes, ", safe maximum=",
         SAFE_MAX_USER_PROCESSES, call. = FALSE)
  }
  if (snapshot$rss_gib >= SAFE_STOP_RSS_GIB) {
    stop(stage, " preflight rejected: current summed RSS=",
         sprintf("%.2f", snapshot$rss_gib), " GiB", call. = FALSE)
  }
  if (is.finite(snapshot$available_memory_gib) &&
      snapshot$available_memory_gib < SAFE_MIN_AVAILABLE_MEMORY_GIB) {
    stop(stage, " preflight rejected: system MemAvailable=",
         sprintf("%.2f", snapshot$available_memory_gib), " GiB", call. = FALSE)
  }
  if (snapshot$rss_gib >= SAFE_WARN_RSS_GIB ||
      (is.finite(snapshot$available_memory_gib) &&
       snapshot$available_memory_gib < SAFE_START_AVAILABLE_MEMORY_GIB)) {
    warning(stage, " resource preflight is in warning range", call. = FALSE)
  }
  message(sprintf(
    "RESOURCE_PREFLIGHT_OK stage=%s workers=%d processes_now=%d projected=%d rss_gib=%.2f mem_available_gib=%.2f",
    stage, workers, snapshot$process_count, projected_processes,
    snapshot$rss_gib, snapshot$available_memory_gib
  ))
  invisible(snapshot)
}

pooled_term_table <- function(pooled, term, effect_scale) {
  stopifnot(effect_scale %in% c("odds_ratio", "identity", "log_response_percent"))
  cf <- pooled$coefficients
  idx <- match(term, names(cf))
  if (is.na(idx)) stop("Exposure term absent from pooled model: ", term, call. = FALSE)
  estimate <- unname(cf[idx])
  se <- unname(sqrt(diag(pooled$variance))[idx])
  df <- pooled$df[idx] %||% Inf
  if (!length(df) || is.na(df) || !is.finite(df) || df <= 0) df <- Inf
  critical <- if (is.finite(df)) stats::qt(0.975, df) else stats::qnorm(0.975)
  low <- unname(estimate - critical * se)
  high <- unname(estimate + critical * se)
  statistic <- unname(estimate / se)
  p <- if (is.finite(df)) 2 * stats::pt(abs(statistic), df, lower.tail = FALSE) else
    2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
  if (effect_scale == "odds_ratio") {
    return(c(
      estimate = exp(estimate), standard_error_link = se,
      ci_low = exp(low), ci_high = exp(high), p_value = p,
      estimate_link = estimate, ci_low_link = low, ci_high_link = high
    ))
  }
  if (effect_scale == "log_response_percent") {
    return(c(
      estimate = 100 * (exp(estimate) - 1), standard_error_link = se,
      ci_low = 100 * (exp(low) - 1), ci_high = 100 * (exp(high) - 1), p_value = p,
      estimate_link = estimate, ci_low_link = low, ci_high_link = high
    ))
  }
  c(
    estimate = estimate, standard_error_link = se,
    ci_low = low, ci_high = high, p_value = p,
    estimate_link = estimate, ci_low_link = low, ci_high_link = high
  )
}

sha256_file <- function(path) {
  digest::digest(path, file = TRUE, algo = "sha256", serialize = FALSE)
}

write_csv_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  utils::write.csv(x, tmp, row.names = FALSE, na = "")
  if (!file.rename(tmp, path)) stop("Atomic rename failed for ", path, call. = FALSE)
  invisible(path)
}

log_stage <- function(stage, message) {
  dir.create(DIRS$logs, recursive = TRUE, showWarnings = FALSE)
  line <- sprintf("%s\t%s\t%s\n", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), stage, message)
  cat(line, file = file.path(DIRS$logs, "phase1_stage.log"), append = TRUE)
  message("[", stage, "] ", message)
}

invisible(TRUE)
)--------------------"

upf_source[["NHANES/02_代码/正式五结局/01_prepare_inputs.R"]] <- r"--------------------(# Build the selected two-day non-water gram-share UPF exposure. Additional
.script_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "NHANES/02_代码/正式五结局/01_prepare_inputs.R"), error = function(e) NULL)
if (is.null(.script_file)) {
  .file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .script_file <- if (length(.file_arg)) sub("^--file=", "", .file_arg[[1]]) else
    "01_prepare_inputs.R"
}
source(file.path(dirname(normalizePath(.script_file, mustWork = FALSE)), "00_functions.R"))

read_xpt_required <- function(path, columns) {
  if (!file.exists(path)) stop("Missing official input: ", path, call. = FALSE)
  x <- haven::read_xpt(path)
  assert_columns(x, columns, basename(path))
  as.data.frame(x[columns])
}

read_manifest <- function(cycle) {
  path <- file.path(NOVA_MANIFEST_DIR, paste0("manifest_", cycle, ".csv"))
  if (!file.exists(path)) stop("Missing NOVA manifest: ", path, call. = FALSE)
  m <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  assert_columns(m, c("code", "desc", "nova", "rule"), basename(path))
  m$code <- as.character(m$code)
  m$nova <- parse_nova_values(m$nova, path)
  if (anyDuplicated(m$code)) stop("Duplicate food code in ", path, call. = FALSE)
  m
}

map_for_cycle <- function(map, cycle) {
  assert_columns(map, c("cycle", "food_code"), "cycle-specific code map")
  map$cycle <- as.character(map$cycle)
  map$food_code <- as.character(map$food_code)
  exact <- map[map$cycle == cycle, , drop = FALSE]
  global <- map[map$cycle %in% c("ALL", "all", "*"), , drop = FALSE]
  global <- global[!global$food_code %in% exact$food_code, , drop = FALSE]
  out <- rbind(exact, global)
  if (anyDuplicated(out$food_code)) stop("Ambiguous mapping for cycle ", cycle, call. = FALSE)
  out
}

diet_roster_day <- function(cycle, day) {
  status <- paste0("DR", day, "DRSTZ")
  kcal <- paste0("DR", day, "TKCAL")
  cols <- c("SEQN", status, kcal, if (day == 2) "DRDINT")
  x <- read_xpt_required(diet_file(cycle, day, "tot"), cols)
  names(x) <- c("SEQN", "recall_status", "total_energy_kcal",
                if (day == 2) "number_of_recalls")
  if (anyDuplicated(x$SEQN)) stop("Duplicate SEQN in ", cycle, " Day ", day, call. = FALSE)
  x
}

collapse_nova_day <- function(cycle, day, roster, water_codes) {
  prefix <- paste0("DR", day, "I")
  columns <- c(
    SEQN = "SEQN", food_code = paste0(prefix, "FDCD"),
    grams = paste0(prefix, "GRMS"), energy = paste0(prefix, "KCAL")
  )
  path <- diet_file(cycle, day, "iff")
  x <- read_xpt_required(path, unname(columns))
  names(x) <- names(columns)
  x$food_code <- ifelse(is.na(x$food_code), NA_character_, sprintf("%.0f", x$food_code))
  x$grams <- as.numeric(x$grams)
  x$energy <- as.numeric(x$energy)
  if (any(x$grams < 0, na.rm = TRUE) || any(x$energy < 0, na.rm = TRUE)) {
    stop("Negative gram or energy value in ", basename(path), call. = FALSE)
  }
  x$grams[!is.finite(x$grams)] <- 0
  x$energy[!is.finite(x$energy)] <- 0

  manifest <- read_manifest(cycle)
  x$nova <- manifest$nova[match(x$food_code, manifest$code)]
  x$is_water <- x$food_code %in% as.character(water_codes)

  base <- roster
  for (k in 1:4) {
    base[[paste0("nova", k, "_g")]] <- 0
    base[[paste0("nova", k, "_kcal")]] <- 0
  }
  base$water_g <- 0
  base$unclassified_g <- 0
  base$unclassified_kcal <- 0
  base$iff_record_count <- 0L

  sum_by_person <- function(value) {
    z <- tapply(value, x$SEQN, sum, na.rm = TRUE)
    ans <- unname(z[match(base$SEQN, names(z))])
    replace(ans, is.na(ans), 0)
  }

  if (nrow(x)) {
    for (k in 1:4) {
      base[[paste0("nova", k, "_g")]] <-
        sum_by_person(ifelse(!x$is_water & x$nova == k, x$grams, 0))
      base[[paste0("nova", k, "_kcal")]] <-
        sum_by_person(ifelse(x$nova == k, x$energy, 0))
    }
    base$water_g <- sum_by_person(ifelse(x$is_water, x$grams, 0))
    base$unclassified_g <- sum_by_person(ifelse(is.na(x$nova) & !x$is_water, x$grams, 0))
    base$unclassified_kcal <- sum_by_person(ifelse(is.na(x$nova), x$energy, 0))
    base$iff_record_count <- sum_by_person(rep(1L, nrow(x)))
  }

  base$classified_nonwater_g <- rowSums(base[paste0("nova", 1:4, "_g")])
  base$total_nonwater_with_unclassified_g <- base$classified_nonwater_g + base$unclassified_g
  base$no_iff_records <- base$iff_record_count == 0
  base$only_plain_water <- base$water_g > 0 & base$total_nonwater_with_unclassified_g == 0
  base$zero_total_energy <- base$recall_status == 1 &
    is.finite(base$total_energy_kcal) & base$total_energy_kcal == 0
  base
}

combine_exposure_days <- function(cycle, d1, d2) {
  one <- d1
  two <- d2
  names(one)[names(one) != "SEQN"] <- paste0(names(one)[names(one) != "SEQN"], "_d1")
  names(two)[names(two) != "SEQN"] <- paste0(names(two)[names(two) != "SEQN"], "_d2")
  z <- merge(one, two, by = "SEQN", all = TRUE, sort = FALSE)
  z$cycle <- cycle

  z$diet_two_day_reliable <- z$recall_status_d1 == 1 &
    z$recall_status_d2 == 1 & z$number_of_recalls_d2 == 2
  z$total_energy_d1 <- as.numeric(z$total_energy_kcal_d1)
  z$total_energy_d2 <- as.numeric(z$total_energy_kcal_d2)
  z$mean_energy_kcal <- (z$total_energy_d1 + z$total_energy_d2) / 2

  z$nova4_energy_d1 <- z$nova4_kcal_d1
  z$nova4_energy_d2 <- z$nova4_kcal_d2
  z$nova4_gram_2d <- z$nova4_g_d1 + z$nova4_g_d2
  classified_cols <- c(paste0("nova", 1:4, "_g_d1"), paste0("nova", 1:4, "_g_d2"))
  z$classified_nonwater_gram_2d <- rowSums(z[classified_cols])
  z$unclassified_gram_2d <- z$unclassified_g_d1 + z$unclassified_g_d2
  z$total_nonwater_with_unclassified_gram_2d <-
    z$classified_nonwater_gram_2d + z$unclassified_gram_2d
  z$unclassified_gram_share <- ifelse(
    z$total_nonwater_with_unclassified_gram_2d > 0,
    z$unclassified_gram_2d / z$total_nonwater_with_unclassified_gram_2d,
    NA_real_
  )
  z$classification_quality_ok <- z$diet_two_day_reliable &
    is.finite(z$unclassified_gram_share) &
    z$unclassified_gram_share <= UNCLASSIFIED_GRAM_SHARE_LIMIT

  total_energy_2d <- z$total_energy_d1 + z$total_energy_d2
  z$upf_energy_ratio_2d <- ifelse(
    z$classification_quality_ok & is.finite(total_energy_2d) & total_energy_2d > 0,
    100 * (z$nova4_energy_d1 + z$nova4_energy_d2) / total_energy_2d,
    NA_real_
  )
  z$upf_energy_mean_daily <- ifelse(
    z$classification_quality_ok & z$total_energy_d1 > 0 & z$total_energy_d2 > 0,
    50 * (z$nova4_energy_d1 / z$total_energy_d1 +
          z$nova4_energy_d2 / z$total_energy_d2),
    NA_real_
  )
  z$upf_gram_ratio_nonwater <- ifelse(
    z$classification_quality_ok & z$classified_nonwater_gram_2d > 0,
    100 * z$nova4_gram_2d / z$classified_nonwater_gram_2d,
    NA_real_
  )

  z$all_day_fasting_d1 <- z$diet_two_day_reliable & z$zero_total_energy_d1 &
    z$no_iff_records_d1 & !z$only_plain_water_d1
  z$all_day_fasting_d2 <- z$diet_two_day_reliable & z$zero_total_energy_d2 &
    z$no_iff_records_d2 & !z$only_plain_water_d2
  z$e_m_not_computable <- z$diet_two_day_reliable &
    (!is.finite(z$total_energy_d1) | !is.finite(z$total_energy_d2) |
       z$total_energy_d1 <= 0 | z$total_energy_d2 <= 0)
  z
}

exposure_qa <- function(x) {
  do.call(rbind, lapply(names(EXPOSURES), function(id) {
    v <- EXPOSURES[[id]]$variable
    z <- x[[v]]
    data.frame(
      exposure = id,
      variable = v,
      n_total = nrow(x),
      n_computable = sum(is.finite(z)),
      missing_fraction = mean(!is.finite(z)),
      min = if (any(is.finite(z))) min(z, na.rm = TRUE) else NA_real_,
      median = if (any(is.finite(z))) stats::median(z, na.rm = TRUE) else NA_real_,
      mean = if (any(is.finite(z))) mean(z, na.rm = TRUE) else NA_real_,
      max = if (any(is.finite(z))) max(z, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}

stage_prepare_inputs <- function() {
  require_formal_authorization("prepare_inputs")
  log_stage("prepare_inputs", "started")
  water_path <- file.path(DIRS$config, "plain_water_codes.csv")
  if (!file.exists(water_path)) stop("Missing plain-water code map", call. = FALSE)
  water_map <- utils::read.csv(water_path, stringsAsFactors = FALSE)
  pieces <- lapply(ANALYTIC_CYCLES, function(cycle) {
    water <- map_for_cycle(water_map, cycle)$food_code
    d1 <- collapse_nova_day(cycle, 1, diet_roster_day(cycle, 1), water)
    d2 <- collapse_nova_day(cycle, 2, diet_roster_day(cycle, 2), water)
    combine_exposure_days(cycle, d1, d2)
  })
  exposure <- rbind_fill(pieces)
  if (any(vapply(EXPOSURES, function(e) {
    z <- exposure[[e$variable]]
    any(is.finite(z) & (z < -1e-8 | z > 100 + 1e-8))
  }, logical(1)))) stop("One or more exposure values fall outside 0-100", call. = FALSE)

  dir.create(DIRS$derived, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(DIRS$derived, "phase1_exposures.rds")
  saveRDS(exposure, path, compress = "xz")
  write_csv_atomic(exposure_qa(exposure), file.path(DIRS$tables, "exposure_build_QA.csv"))
  log_stage("prepare_inputs", paste("completed", nrow(exposure), "participant-cycle rows"))
  path
}

if (sys.nframe() == 0L) stage_prepare_inputs()
)--------------------"

upf_source[["NHANES/02_代码/正式五结局/02_build_dataset.R"]] <- r"--------------------(# NHANES phase-1 participant-level analytic frame.

.script_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "NHANES/02_代码/正式五结局/02_build_dataset.R"), error = function(e) NULL)
if (is.null(.script_file)) {
  .file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .script_file <- if (length(.file_arg)) sub("^--file=", "", .file_arg[[1]]) else
    "02_build_dataset.R"
}
source(file.path(dirname(normalizePath(.script_file, mustWork = FALSE)), "00_functions.R"))

split_source_variables <- function(x) {
  unique(trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE), use.names = FALSE)))
}

read_phase1_variable_contract <- function() {
  pieces <- lapply(PHASE1_VARIABLE_MAP_FILES, function(path) {
    if (!file.exists(path)) stop("Missing variable contract: ", path, call. = FALSE)
    z <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    assert_columns(z, c("cycle", "module", "stem", "source_variables"), basename(path))
    z[c("cycle", "module", "stem", "source_variables")]
  })
  z <- unique(rbind_fill(pieces))
  z <- z[z$cycle %in% ANALYTIC_CYCLES, , drop = FALSE]
  if (!nrow(z)) stop("Variable contract contains no analytic-cycle rows", call. = FALSE)
  z
}

person_module_path <- function(cycle, module, stem) {
  if (module == "Dietary" && stem %in% c("dr1tot", "dr2tot")) {
    day <- as.integer(substr(stem, 3, 3))
    return(diet_file(cycle, day, "tot"))
  }
  module_file(cycle, module, stem)
}

read_keep <- function(path, columns) {
  if (!file.exists(path)) stop("Missing source module: ", path, call. = FALSE)
  x <- haven::read_xpt(path)
  missing <- setdiff(columns, names(x))
  if (length(missing)) {
    stop("Configured variable(s) absent from ", path, ": ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  as.data.frame(x[columns])
}

left_merge_unique <- function(x, y) {
  key <- if ("cycle" %in% names(y)) c("SEQN", "cycle") else "SEQN"
  if (anyDuplicated(y[key])) stop("Duplicate participant merge key", call. = FALSE)
  overlap <- setdiff(intersect(names(x), names(y)), key)
  if (length(overlap)) {
    stop("Unexpected overlapping fields: ", paste(overlap, collapse = ", "), call. = FALSE)
  }
  merge(x, y, by = key, all.x = TRUE, sort = FALSE)
}

read_cycle_person <- function(cycle, contract = read_phase1_variable_contract()) {
  vm <- contract[contract$cycle == cycle &
                   !contract$stem %in% c("dr1iff", "dr2iff"), , drop = FALSE]
  keys <- unique(vm[c("module", "stem")])
  priority <- match(paste(keys$module, keys$stem, sep = "::"),
                    c("Demographics::demo", "Dietary::dr1tot", "Dietary::dr2tot"))
  keys <- keys[order(ifelse(is.na(priority), 99L, priority), keys$module, keys$stem), ]
  pieces <- lapply(seq_len(nrow(keys)), function(i) {
    rows <- vm$module == keys$module[[i]] & vm$stem == keys$stem[[i]]
    columns <- unique(c("SEQN", split_source_variables(vm$source_variables[rows])))
    read_keep(person_module_path(cycle, keys$module[[i]], keys$stem[[i]]), columns)
  })
  demo_index <- which(keys$module == "Demographics" & keys$stem == "demo")
  if (length(demo_index) != 1L) stop("Exactly one demographics file is required", call. = FALSE)
  d <- pieces[[demo_index]]
  if (anyDuplicated(d$SEQN)) stop("Duplicate SEQN in demographics cycle ", cycle, call. = FALSE)
  d$cycle <- cycle
  for (i in setdiff(seq_along(pieces), demo_index)) d <- left_merge_unique(d, pieces[[i]])
  d
}

ensure_raw_columns <- function(d) {
  needed <- c(
    "RIAGENDR", "RIDAGEYR", "RIDRETH1", "DMDEDUC2", "INDFMPIR", "RIDEXPRG",
    "SDMVPSU", "SDMVSTRA", "WTDR2D", "WTDR2DPP", "DR1DRSTZ", "DR2DRSTZ",
    "DRDINT", "DR1TKCAL", "DR2TKCAL", "DIQ010", "DIQ050", "DIQ070", "DID070",
    "SMQ020", "SMQ040", "BPQ020", "BPQ050A", "BPQ080", "BPQ090D", "BPQ100D",
    "KIQ026", "LBXGH", "LBXTC", "LBXSCR", "LBXSUA", "URXUMA", "URXUCR",
    "URDACT", "BMXBMI", "ALQ101", "ALQ110", "ALQ120Q", "ALQ120U", "ALQ111",
    "ALQ121", "ALQ130", PHYSICAL_ACTIVITY_VARS,
    paste0("BPXSY", 1:4), paste0("BPXDI", 1:4),
    paste0("BPXOSY", 1:3), paste0("BPXODI", 1:3)
  )
  for (nm in setdiff(needed, names(d))) d[[nm]] <- NA_real_
  d
}

derive_nonpregnant <- function(sex, age, pregnancy) {
  reproductive_age <- sex == "Female" & age >= 20 & age <= 44
  ifelse(sex == "Male" | (sex == "Female" & !reproductive_age), 1L,
         ifelse(reproductive_age & pregnancy == 2, 1L,
                ifelse(reproductive_age & pregnancy == 1, 0L, NA_integer_)))
}

standardize_cycle <- function(d) {
  d <- ensure_raw_columns(d)
  d$age <- as.numeric(d$RIDAGEYR)
  d$age10 <- (d$age - 50) / 10
  d$sex <- factor(d$RIAGENDR, levels = c(1, 2), labels = c("Male", "Female"))
  d$race_ethnicity <- factor(
    d$RIDRETH1, levels = 1:5,
    labels = c("Mexican American", "Other Hispanic", "Non-Hispanic White",
               "Non-Hispanic Black", "Other/Multiracial")
  )
  d$education <- factor(
    ifelse(d$DMDEDUC2 %in% 1:2, "Below_high_school",
           ifelse(d$DMDEDUC2 == 3, "High_school_GED",
                  ifelse(d$DMDEDUC2 %in% 4:5, "Above_high_school", NA_character_))),
    levels = c("Below_high_school", "High_school_GED", "Above_high_school")
  )
  d$pir <- as.numeric(d$INDFMPIR)
  d$bmi <- as.numeric(d$BMXBMI)
  d$cycle_factor <- factor(d$cycle, levels = ANALYTIC_CYCLES)
  d$smoking <- derive_smoking(d$SMQ020, d$SMQ040)
  d$alcohol <- derive_alcohol_weekly(
    d$cycle, d$sex, d$ALQ101, d$ALQ110, d$ALQ120Q, d$ALQ120U,
    d$ALQ111, d$ALQ121, d$ALQ130
  )
  pa <- derive_gpaq_met_category(d)
  d$physical_activity_met_min_week <- pa$total_met_min_week
  d$physical_activity <- pa$category
  d$mean_sbp <- ifelse(
    d$cycle == "P",
    nhanes_bp_mean_oscillometric(d[paste0("BPXOSY", 1:3)]),
    nhanes_bp_mean_auscultatory(d[paste0("BPXSY", 1:4)])
  )
  d$mean_dbp <- ifelse(
    d$cycle == "P",
    nhanes_bp_mean_oscillometric(d[paste0("BPXODI", 1:3)], diastolic = TRUE),
    nhanes_bp_mean_auscultatory(d[paste0("BPXDI", 1:4)], diastolic = TRUE)
  )
  d$total_cholesterol <- as.numeric(d$LBXTC)
  d$hba1c <- as.numeric(d$LBXGH)
  oral <- ifelse(d$cycle == "E" & !is.na(d$DID070), d$DID070, d$DIQ070)
  dm <- derive_diabetes_components(d$DIQ010, d$DIQ050, oral, d$hba1c)
  d$diabetes_diagnosed <- dm$diagnosed
  d$insulin_current <- dm$insulin
  d$oral_med_current <- dm$oral_med
  d$diabetes <- dm$diabetes
  d$diabetes_original <- d$diabetes
  ht <- derive_hypertension_components(d$mean_sbp, d$mean_dbp, d$BPQ020, d$BPQ050A)
  d$hypertension_diagnosed <- ht$diagnosed
  d$hypertension_med <- ht$medication
  d$hypertension <- ht$hypertension
  ch <- derive_cholesterol_components(
    d$total_cholesterol, d$BPQ080, d$BPQ090D, d$BPQ100D
  )
  d$cholesterol_diagnosed <- ch$diagnosed
  d$cholesterol_med <- ch$medication
  d$high_cholesterol <- ch$high_cholesterol
  d$nonpregnant <- derive_nonpregnant(d$sex, d$age, d$RIDEXPRG)
  d$diet_two_day_reliable_tot <- as.integer(
    d$DR1DRSTZ == 1 & d$DR2DRSTZ == 1 & d$DRDINT == 2
  )
  d$mean_energy_kcal <- rowMeans(cbind(d$DR1TKCAL, d$DR2TKCAL), na.rm = FALSE)
  d$analysis_weight <- derive_analysis_weight(d$cycle, d$WTDR2D, d$WTDR2DPP)
  d
}

derive_phase1_outcomes <- function(d) {
  d$creatinine_dxc <- calibrate_creatinine_dxc(as.numeric(d$LBXSCR), d$cycle)
  d$uric_acid_dxc <- calibrate_uric_acid_dxc(as.numeric(d$LBXSUA), d$cycle)
  d$egfr_2021 <- calc_egfr_2021(d$creatinine_dxc, d$age, d$sex == "Female")
  released_acr <- as.numeric(d$URDACT)
  calculated_acr <- ifelse(
    is.finite(d$URXUMA) & is.finite(d$URXUCR) & d$URXUCR > 0,
    100 * as.numeric(d$URXUMA) / as.numeric(d$URXUCR), NA_real_
  )
  d$acr_mg_g <- ifelse(d$cycle == "E" & is.na(released_acr), calculated_acr, released_acr)
  d$log_uacr <- ifelse(
    is.finite(d$acr_mg_g) & d$acr_mg_g > 0,
    log(d$acr_mg_g), NA_real_
  )
  d$albuminuria <- ifelse(is.na(d$acr_mg_g), NA_integer_, as.integer(d$acr_mg_g >= 30))
  d$ckd <- derive_ckd_three_state(d$egfr_2021, d$acr_mg_g)
  d$dkd <- ifelse(d$diabetes_original == 1L, d$ckd, NA_integer_)
  d$hyperuricemia <- ifelse(
    is.na(d$uric_acid_dxc) | is.na(d$sex), NA_integer_,
    as.integer(ifelse(d$sex == "Male", d$uric_acid_dxc > 7.0,
                      d$uric_acid_dxc > 5.7))
  )
  d$kidney_stones <- yes_no_na(d$KIQ026)
  d
}

attach_phase1_exposures <- function(d, exposure) {
  assert_columns(
    exposure,
    c("SEQN", "cycle", "diet_two_day_reliable", "classification_quality_ok",
      vapply(EXPOSURES, `[[`, character(1), "variable")),
    "phase1_exposures"
  )
  keep <- unique(c(
    "SEQN", "cycle", "diet_two_day_reliable", "classification_quality_ok",
    "unclassified_gram_share", "total_energy_d1", "total_energy_d2",
    "all_day_fasting_d1", "all_day_fasting_d2", "e_m_not_computable",
    vapply(EXPOSURES, `[[`, character(1), "variable")
  ))
  d <- left_merge_unique(d, exposure[intersect(keep, names(exposure))])
  compare <- !is.na(d$diet_two_day_reliable) & !is.na(d$diet_two_day_reliable_tot)
  if (any(d$diet_two_day_reliable[compare] != as.logical(d$diet_two_day_reliable_tot[compare]))) {
    stop("Dietary reliability mismatch between exposure and participant stages", call. = FALSE)
  }
  d
}

create_domains <- function(d) {
  d$design_eligible <- complete.cases(d[c("SDMVPSU", "SDMVSTRA", "analysis_weight")]) &
    is.finite(d$analysis_weight) & d$analysis_weight > 0
  d$base_domain <- (d$age >= 20) %in% TRUE & (d$nonpregnant == 1L) %in% TRUE &
    d$diet_two_day_reliable %in% TRUE & d$classification_quality_ok %in% TRUE
  d$acr_observed <- d$base_domain & !is.na(d$albuminuria)
  d$uacr_observed <- d$base_domain & is.finite(d$log_uacr)
  d$egfr_observed <- d$base_domain & !is.na(d$egfr_2021)
  d$kidney_observed <- d$base_domain & !is.na(d$ckd)
  d$dkd_domain <- d$base_domain & d$diabetes_original %in% 1L & !is.na(d$dkd)
  d$uric_acid_observed <- d$base_domain & !is.na(d$uric_acid_dxc)
  d$stone_observed <- d$base_domain & !is.na(d$kidney_stones)
  d
}

weighted_missing_fraction <- function(design, domain, variable) {
  keep <- design$variables[[domain]] %in% TRUE
  if (!any(keep)) return(NA_real_)
  z <- design[keep, ]
  z$variables$.missing_indicator <- as.numeric(is.na(z$variables[[variable]]))
  unname(stats::coef(survey::svymean(~.missing_indicator, z, na.rm = TRUE))[[1]])
}

audit_pir <- function(d) {
  full <- make_full_design(d)
  domain_map <- c(base = "base_domain", vapply(OUTCOMES, `[[`, character(1), "domain"))
  domain_map <- domain_map[!duplicated(domain_map)]
  rows <- lapply(names(domain_map), function(label) {
    domain <- unname(domain_map[[label]])
    idx <- d$design_eligible & d[[domain]] %in% TRUE
    data.frame(
      scope = label, domain = domain, cycle = "ALL",
      n_domain = sum(idx), n_pir_missing = sum(idx & is.na(d$pir)),
      unweighted_missing_fraction = if (sum(idx)) mean(is.na(d$pir[idx])) else NA_real_,
      weighted_missing_fraction = weighted_missing_fraction(full, domain, "pir"),
      stringsAsFactors = FALSE
    )
  })
  cycle_rows <- lapply(ANALYTIC_CYCLES, function(cycle) {
    idx <- d$design_eligible & d$base_domain %in% TRUE & d$cycle == cycle
    data.frame(
      scope = "base_by_cycle", domain = "base_domain", cycle = cycle,
      n_domain = sum(idx), n_pir_missing = sum(idx & is.na(d$pir)),
      unweighted_missing_fraction = if (sum(idx)) mean(is.na(d$pir[idx])) else NA_real_,
      weighted_missing_fraction = NA_real_, stringsAsFactors = FALSE
    )
  })
  audit <- do.call(rbind, c(rows, cycle_rows))
  overall <- audit$unweighted_missing_fraction[audit$scope == "base"][[1]]
  include_pir <- FALSE
  decision <- list(
    version = PHASE1_VERSION,
    decision_rule = "PIR excluded from the frozen phase-1 v3 model specification",
    threshold = PIR_MISSING_LIMIT,
    observed_fraction = overall,
    include_pir = include_pir,
    decision = "omit_from_all_phase1_models",
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  dir.create(dirname(PIR_DECISION_FILE), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(
    decision, PIR_DECISION_FILE, pretty = TRUE, auto_unbox = TRUE, digits = 16
  )
  write_csv_atomic(audit, file.path(DIRS$tables, "pir_missingness_audit.csv"))
  list(audit = audit, decision = decision)
}

stage_build_dataset <- function() {
  require_formal_authorization("build_dataset")
  log_stage("build_dataset", "started")
  contract <- read_phase1_variable_contract()
  pieces <- lapply(ANALYTIC_CYCLES, function(cycle) {
    log_stage("build_dataset", paste("reading cycle", cycle))
    standardize_cycle(read_cycle_person(cycle, contract))
  })
  d <- rbind_fill(pieces)
  if (anyDuplicated(d[c("SEQN", "cycle")])) stop("Duplicate participant-cycle key", call. = FALSE)
  d <- derive_phase1_outcomes(d)
  exposure_path <- file.path(DIRS$derived, "phase1_exposures.rds")
  if (!file.exists(exposure_path)) stop("Run 01_prepare_inputs.R first", call. = FALSE)
  d <- attach_phase1_exposures(d, readRDS(exposure_path))
  d <- create_domains(d)
  pir <- audit_pir(d)
  out <- file.path(DIRS$derived, "analysis_frame_pre_mice.rds")
  saveRDS(d, out, compress = "xz")
  write_csv_atomic(
    data.frame(
      cycle = ANALYTIC_CYCLES,
      n_design_universe = vapply(ANALYTIC_CYCLES, function(x) sum(d$cycle == x & d$design_eligible), numeric(1)),
      n_base_domain = vapply(ANALYTIC_CYCLES, function(x) sum(d$cycle == x & d$base_domain), numeric(1)),
      stringsAsFactors = FALSE
    ),
    file.path(DIRS$tables, "analysis_frame_counts_by_cycle.csv")
  )
  log_stage("build_dataset", paste("completed", nrow(d), "rows; PIR", pir$decision$decision))
  out
}

if (sys.nframe() == 0L) stage_build_dataset()
)--------------------"

upf_source[["NHANES/02_代码/正式五结局/03_imputation.R"]] <- r"--------------------(# Outcome-specific multiple imputation for the formal five-outcome analysis.

.script_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "NHANES/02_代码/正式五结局/03_imputation.R"), error = function(e) NULL)
if (is.null(.script_file)) {
  .file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .script_file <- if (length(.file_arg)) sub("^--file=", "", .file_arg[[1]]) else
    "03_imputation.R"
}
source(file.path(dirname(normalizePath(.script_file, mustWork = FALSE)), "00_functions.R"))

read_pir_decision <- function() {
  if (!file.exists(PIR_DECISION_FILE)) stop("PIR decision is missing; run stage 02", call. = FALSE)
  jsonlite::read_json(PIR_DECISION_FILE, simplifyVector = TRUE)
}

mi_target_variables <- function(include_pir) {
  if (isTRUE(include_pir)) stop("PIR is excluded from formal five-outcome analysis", call. = FALSE)
  unique(c(
    "bmi", "mean_energy_kcal", "education", "smoking", "alcohol",
    "physical_activity", "mean_sbp", "mean_dbp", "total_cholesterol", "hba1c",
    "diabetes_diagnosed", "insulin_current", "oral_med_current",
    "hypertension", "high_cholesterol"
  ))
}

mi_fixed_candidates <- function(outcome) {
  unique(c(
    "SEQN", "age10", "sex", "race_ethnicity", "cycle_factor",
    "log_analysis_weight", "mi_stratum_code", "mi_psu_code",
    OUTCOMES[[outcome]]$variable,
    vapply(EXPOSURES, `[[`, character(1), "variable")
  ))
}

prepare_mi_subframe <- function(d, outcome, include_pir) {
  domain <- OUTCOMES[[outcome]]$domain
  row_index <- which(d$design_eligible & d[[domain]] %in% TRUE)
  if (!length(row_index)) stop("No records in outcome domain: ", outcome, call. = FALSE)
  z <- d[row_index, , drop = FALSE]
  z$log_analysis_weight <- log(z$analysis_weight)
  z$mi_stratum_code <- as.integer(factor(
    interaction(z$cycle, z$SDMVSTRA, drop = TRUE)
  ))
  z$mi_psu_code <- as.integer(factor(
    interaction(z$cycle, z$SDMVSTRA, z$SDMVPSU, drop = TRUE)
  ))
  binary_targets <- c(
    "diabetes_diagnosed", "insulin_current", "oral_med_current",
    "hypertension", "high_cholesterol"
  )
  for (nm in intersect(binary_targets, names(z))) {
    z[[nm]] <- factor(z[[nm]], levels = c(0, 1))
  }
  targets <- mi_target_variables(include_pir)
  fixed <- mi_fixed_candidates(outcome)
  keep <- unique(c(fixed, targets))
  assert_columns(z, keep, paste0("MI source for ", outcome))
  z <- z[keep]
  list(data = z, row_index = row_index, targets = targets, fixed = fixed)
}

mice_convergence_tables <- function(mids, outcome) {
  conv <- mice::convergence(mids)
  conv$outcome <- outcome
  final <- conv[conv$.it == max(conv$.it, na.rm = TRUE), , drop = FALSE]
  data.frame(
    outcome = outcome,
    maxit = mids$iteration,
    n_chains = mids$m,
    n_variables_checked = nrow(final),
    max_psrf_final = if (any(is.finite(final$psrf))) max(final$psrf, na.rm = TRUE) else NA_real_,
    n_psrf_over_1_10 = sum(final$psrf > 1.10, na.rm = TRUE),
    max_abs_autocorrelation_final = if (any(is.finite(final$ac)))
      max(abs(final$ac), na.rm = TRUE) else NA_real_,
    n_abs_autocorrelation_over_0_50 = sum(abs(final$ac) > 0.50, na.rm = TRUE),
    stringsAsFactors = FALSE
  ) -> summary
  list(detail = conv, summary = summary)
}

build_mice_spec <- function(z, targets, fixed) {
  method <- rep("", ncol(z)); names(method) <- names(z)
  for (nm in intersect(targets, names(MICE_BASE_METHODS))) {
    if (anyNA(z[[nm]]) && length(unique(z[[nm]][!is.na(z[[nm]])])) > 1L) {
      method[[nm]] <- MICE_BASE_METHODS[[nm]]
    }
  }
  predictor_matrix <- matrix(0L, nrow = ncol(z), ncol = ncol(z),
                             dimnames = list(names(z), names(z)))
  imputed_targets <- names(method)[nzchar(method)]
  complete_fixed <- fixed[vapply(z[fixed], function(x) !anyNA(x), logical(1))]
  predictors <- unique(c(imputed_targets, complete_fixed))
  predictors <- setdiff(predictors, c("SEQN"))
  for (target in imputed_targets) {
    predictor_matrix[target, setdiff(predictors, target)] <- 1L
  }
  list(method = method, predictor_matrix = predictor_matrix,
       imputed_targets = imputed_targets, complete_fixed = complete_fixed)
}

imputation_audit <- function(z, outcome, targets, spec) {
  data.frame(
    outcome = outcome,
    variable = targets,
    n = nrow(z),
    n_missing = vapply(z[targets], function(x) sum(is.na(x)), numeric(1)),
    missing_fraction = vapply(z[targets], function(x) mean(is.na(x)), numeric(1)),
    mice_method = unname(spec$method[targets]),
    imputed = unname(nzchar(spec$method[targets])),
    stringsAsFactors = FALSE
  )
}

run_outcome_mice <- function(d, outcome, include_pir) {
  prepared <- prepare_mi_subframe(d, outcome, include_pir)
  z <- prepared$data
  spec <- build_mice_spec(z, prepared$targets, prepared$fixed)
  if (!length(spec$imputed_targets)) stop("No MICE targets have missing values for ", outcome, call. = FALSE)
  set.seed(PHASE1_SEED + match(outcome, names(OUTCOMES)))
  mids <- mice::futuremice(
    z, m = MI_M, maxit = MI_MAXIT,
    method = spec$method, predictorMatrix = spec$predictor_matrix,
    seed = NA,
    parallelseed = PHASE1_SEED + 1000L + match(outcome, names(OUTCOMES)),
    n.core = min(n_cores_resolve(), MI_M), future.plan = "multicore",
    packages = "mice",
    remove.constant = TRUE, remove.collinear = TRUE
  )
  if (mids$m != MI_M) stop("MICE did not create exactly ", MI_M, " imputations", call. = FALSE)
  object <- list(
    version = PHASE1_VERSION, outcome = outcome,
    outcome_variable = OUTCOMES[[outcome]]$variable,
    domain = OUTCOMES[[outcome]]$domain,
    include_pir = include_pir, m = mids$m, maxit = MI_MAXIT,
    original_row_index = prepared$row_index,
    imputed_targets = spec$imputed_targets,
    complete_fixed_predictors = spec$complete_fixed,
    predictor_matrix = spec$predictor_matrix,
    mids = mids
  )
  list(
    object = object,
    audit = imputation_audit(z, outcome, prepared$targets, spec),
    events = if (is.null(mids$loggedEvents)) data.frame() else mids$loggedEvents,
    convergence = mice_convergence_tables(mids, outcome)
  )
}

stage_imputation <- function() {
  require_formal_authorization("imputation")
  log_stage("imputation", "started")
  source_path <- file.path(DIRS$derived, "analysis_frame_pre_mice.rds")
  if (!file.exists(source_path)) stop("Run 02_build_dataset.R first", call. = FALSE)
  d <- readRDS(source_path)
  pir <- read_pir_decision()
  include_pir <- isTRUE(pir$include_pir)
  all_audits <- list(); all_events <- list()
  for (outcome in names(OUTCOMES)) {
    log_stage("imputation", paste("running", outcome, "m=", MI_M, "maxit=", MI_MAXIT))
    ans <- run_outcome_mice(d, outcome, include_pir)
    path <- file.path(DIRS$derived, paste0("mi_", outcome, ".rds"))
    saveRDS(ans$object, path, compress = "xz")
    all_audits[[outcome]] <- ans$audit
    if (nrow(ans$events)) {
      ans$events$outcome <- outcome
      all_events[[outcome]] <- ans$events
    }
    write_csv_atomic(
      as.data.frame(ans$object$predictor_matrix),
      file.path(DIRS$tests, paste0("mice_predictor_matrix_", outcome, ".csv"))
    )
    write_csv_atomic(ans$convergence$detail,
                     file.path(DIRS$tests, paste0("mice_convergence_", outcome, ".csv")))
    write_csv_atomic(ans$convergence$summary,
                     file.path(DIRS$tests, paste0("mice_convergence_summary_", outcome, ".csv")))
  }
  write_csv_atomic(do.call(rbind, all_audits), file.path(DIRS$tables, "mice_missingness_and_methods.csv"))
  events <- if (length(all_events)) rbind_fill(all_events) else
    data.frame(out = integer(), ins = integer(), dep = character(), meth = character(),
               outcome = character())
  write_csv_atomic(events, file.path(DIRS$tests, "mice_logged_events.csv"))
  log_stage("imputation", paste("completed", length(OUTCOMES), "outcome-specific objects"))
  invisible(file.path(DIRS$derived, paste0("mi_", names(OUTCOMES), ".rds")))
}

if (sys.nframe() == 0L) stage_imputation()
)--------------------"

upf_source[["NHANES/02_代码/正式五结局/04_models.R"]] <- r"--------------------(# Formal five-outcome continuous-exposure survey models across 50 imputations.

.script_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "NHANES/02_代码/正式五结局/04_models.R"), error = function(e) NULL)
if (is.null(.script_file)) {
  .file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .script_file <- if (length(.file_arg)) sub("^--file=", "", .file_arg[[1]]) else
    "04_models.R"
}
source(file.path(dirname(normalizePath(.script_file, mustWork = FALSE)), "00_functions.R"))

read_pir_inclusion <- function() {
  FALSE
}

restore_completed_subframe <- function(full, mi_object, imputation_number) {
  completed <- mice::complete(mi_object$mids, action = imputation_number)
  idx <- mi_object$original_row_index
  if (nrow(completed) != length(idx)) stop("Completed MI row count mismatch", call. = FALSE)
  binary_targets <- c(
    "diabetes_diagnosed", "insulin_current", "oral_med_current",
    "hypertension", "high_cholesterol"
  )
  for (nm in intersect(names(completed), names(full))) {
    value <- completed[[nm]]
    if (nm %in% binary_targets && is.factor(value)) value <- as.integer(as.character(value))
    full[[nm]][idx] <- value
  }
  dm_lab <- ifelse(is.na(full$hba1c), NA_integer_, as.integer(full$hba1c >= 6.5))
  full$diabetes <- or_composite(
    full$diabetes_diagnosed, full$insulin_current,
    full$oral_med_current, dm_lab
  )
  full$age10 <- (full$age - 50) / 10
  full
}

make_outcome_grid <- function(outcome, d) {
  natural <- expand.grid(
    outcome = outcome, exposure = names(EXPOSURES), model = names(MODEL_EXTENSIONS),
    age_form = names(AGE_FORMS), stringsAsFactors = FALSE
  )
  natural$sample_rule <- "natural"
  natural$spec_id <- sprintf(
    "%s__%s__%s__%s__%s", natural$outcome, natural$exposure,
    natural$model, natural$age_form, natural$sample_rule
  )
  natural
}

spec_domain <- function(variables, spec) {
  domain <- variables[[OUTCOMES[[spec$outcome]]$domain]] %in% TRUE
  exposure_var <- EXPOSURES[[spec$exposure]]$variable
  domain <- domain & is.finite(variables[[exposure_var]])
  domain
}

fit_one_imputation <- function(imputation_number, d, mi_object, grid, include_pir) {
  completed <- restore_completed_subframe(d, mi_object, imputation_number)
  full_design <- make_full_design(completed)
  summaries <- vector("list", nrow(grid)); names(summaries) <- grid$spec_id
  audit <- vector("list", nrow(grid)); names(audit) <- grid$spec_id
  for (j in seq_len(nrow(grid))) {
    spec <- grid[j, , drop = FALSE]
    form <- phase1_formula(spec$exposure, spec$model, spec$age_form,
                           spec$outcome, include_pir)
    domain <- spec_domain(full_design$variables, spec)
    model_vars <- all.vars(form)
    incomplete <- model_vars[vapply(full_design$variables[model_vars], function(x) {
      any(is.na(x[domain]))
    }, logical(1))]
    if (length(incomplete)) {
      stop("Post-MICE missing model fields in ", spec$spec_id, ": ",
           paste(incomplete, collapse = ", "), call. = FALSE)
    }
    if (!any(domain)) stop("Empty analysis domain: ", spec$spec_id, call. = FALSE)
    design <- full_design[domain, ]
    fit <- survey::svyglm(
      form, design = design, family = outcome_family(spec$outcome), na.action = na.fail
    )
    coefficients <- stats::coef(fit)
    variance <- stats::vcov(fit)
    if (!identical(names(coefficients), rownames(variance)) ||
        !identical(names(coefficients), colnames(variance))) {
      stop("Coefficient/covariance mismatch inside worker: ", spec$spec_id, call. = FALSE)
    }
    summaries[[spec$spec_id]] <- list(
      coefficients = coefficients,
      variance = variance,
      df_complete = stats::df.residual(fit)
    )
    y <- design$variables[[OUTCOMES[[spec$outcome]]$variable]]
    sampling_weights <- as.numeric(stats::weights(design, "sampling"))
    audit[[spec$spec_id]] <- c(
      n_unweighted = nrow(design$variables),
      events = if (OUTCOMES[[spec$outcome]]$family == "quasibinomial") sum(y == 1) else NA_real_,
      weighted_population = sum(sampling_weights),
      effective_sample_size = sum(sampling_weights)^2 / sum(sampling_weights^2),
      survey_design_df = survey::degf(design),
      complete_data_residual_df = stats::df.residual(fit)
    )
    rm(fit, coefficients, variance, design)
  }
  list(summaries = summaries, audit = audit)
}

pool_outcome_grid <- function(imputation_results, grid, include_pir) {
  rows <- lapply(seq_len(nrow(grid)), function(j) {
    spec <- grid[j, , drop = FALSE]
    summaries <- lapply(imputation_results, function(x) x$summaries[[spec$spec_id]])
    pooled <- pool_svy_summaries(summaries)
    exposure_term <- sprintf(
      "I(%s/%s)", EXPOSURES[[spec$exposure]]$variable, EXPOSURE_UNIT_PP
    )
    effect_scale <- OUTCOMES[[spec$outcome]]$effect_scale
    estimate <- pooled_term_table(pooled, exposure_term, effect_scale)
    audit <- imputation_results[[1]]$audit[[spec$spec_id]]
    exposure_index <- match(exposure_term, names(pooled$coefficients))
    link_coefficients <- vapply(
      summaries, function(x) unname(x$coefficients[[exposure_term]]), numeric(1)
    )
    between_variance <- stats::var(link_coefficients)
    mc_error <- sqrt(between_variance / length(link_coefficients))
    total_se <- estimate[["standard_error_link"]]
    data.frame(
      version = PHASE1_VERSION,
      spec_id = spec$spec_id, outcome = spec$outcome, exposure = spec$exposure,
      exposure_variable = EXPOSURES[[spec$exposure]]$variable,
      exposure_increment = EXPOSURE_UNIT_PP,
      model = spec$model, age_form = spec$age_form, sample_rule = spec$sample_rule,
      include_pir = include_pir,
      effect_measure = unname(c(
        odds_ratio = "OR_per_10_percentage_points",
        identity = "beta_per_10_percentage_points",
        log_response_percent = "percent_change_per_10_percentage_points"
      )[[effect_scale]]),
      estimate = estimate[["estimate"]],
      standard_error = estimate[["standard_error_link"]],
      ci_low = estimate[["ci_low"]], ci_high = estimate[["ci_high"]],
      p_value = estimate[["p_value"]],
      estimate_link = estimate[["estimate_link"]],
      ci_low_link = estimate[["ci_low_link"]],
      ci_high_link = estimate[["ci_high_link"]],
      n_unweighted = audit[["n_unweighted"]], events = audit[["events"]],
      weighted_population = audit[["weighted_population"]],
      effective_sample_size = audit[["effective_sample_size"]],
      survey_design_df = audit[["survey_design_df"]],
      complete_data_residual_df = audit[["complete_data_residual_df"]],
      pooled_df = unname(pooled$df[[exposure_index]]),
      fmi_percent = unname(pooled$missinfo[[exposure_index]]),
      between_imputation_variance_link_scale = between_variance,
      monte_carlo_error_link_scale = mc_error,
      monte_carlo_error_over_total_se = mc_error / total_se,
      successful_imputations = length(summaries),
      fit_status = paste0("success_", length(summaries), "_of_", MI_M),
      formula = paste(deparse(phase1_formula(
        spec$exposure, spec$model, spec$age_form, spec$outcome, include_pir
      )), collapse = " "),
      covariates = paste(phase1_model_covariates(
        spec$model, spec$age_form, spec$outcome, include_pir
      ), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

stage_models <- function() {
  require_formal_authorization("models")
  log_stage("models", "started")
  source_path <- file.path(DIRS$derived, "analysis_frame_pre_mice.rds")
  if (!file.exists(source_path)) stop("Pre-MICE frame is missing", call. = FALSE)
  d <- readRDS(source_path)
  include_pir <- read_pir_inclusion()
  all_results <- list()
  for (outcome in names(OUTCOMES)) {
    mi_path <- file.path(DIRS$derived, paste0("mi_", outcome, ".rds"))
    if (!file.exists(mi_path)) stop("Missing outcome MI object: ", mi_path, call. = FALSE)
    mi_object <- readRDS(mi_path)
    if (mi_object$m != MI_M) stop("MI object is not m=", MI_M, ": ", outcome, call. = FALSE)
    grid <- make_outcome_grid(outcome, d)
    log_stage("models", paste("fitting", nrow(grid), "specifications for", outcome,
                               "across", MI_M, "imputations"))
    workers <- n_cores_resolve()
    resource_preflight(workers, paste0("models_", outcome))
    imputation_results <- parallel::mclapply(
      seq_len(MI_M), fit_one_imputation,
      d = d, mi_object = mi_object, grid = grid, include_pir = include_pir,
      mc.cores = workers, mc.preschedule = TRUE
    )
    failed <- vapply(imputation_results, inherits, logical(1), "try-error")
    if (any(failed)) stop("One or more parallel model workers failed for ", outcome, call. = FALSE)
    result <- pool_outcome_grid(imputation_results, grid, include_pir)
    write_csv_atomic(result, file.path(DIRS$results, paste0("phase1_models_", outcome, ".csv")))
    saveRDS(result, file.path(DIRS$results, paste0("phase1_models_", outcome, ".rds")),
            compress = "xz")
    all_results[[outcome]] <- result
    rm(imputation_results, mi_object, grid)
    invisible(gc(full = TRUE))
  }
  results <- do.call(rbind, all_results)
  rownames(results) <- NULL
  write_csv_atomic(results, file.path(DIRS$results, "phase1_models_all_outcomes.csv"))
  saveRDS(results, file.path(DIRS$results, "phase1_models_all_outcomes.rds"), compress = "xz")
  natural_n <- sum(results$sample_rule == "natural")
  expected_n <- length(OUTCOMES) * length(EXPOSURES) *
    length(MODEL_COVARIATES) * length(AGE_FORMS)
  if (natural_n != expected_n) {
    stop("Natural model-grid count is not ", expected_n, ": ", natural_n, call. = FALSE)
  }
  log_stage("models", paste("completed", nrow(results), "pooled specifications"))
  invisible(results)
}

if (sys.nframe() == 0L) stage_models()
)--------------------"

upf_source[["NHANES/02_代码/正式五结局/06_postbuild_descriptives.R"]] <- r"--------------------(# Prespecified exposure/sample descriptives before any model-result inspection.

.script_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "NHANES/02_代码/正式五结局/06_postbuild_descriptives.R"), error = function(e) NULL)
if (is.null(.script_file)) {
  .file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .script_file <- if (length(.file_arg)) sub("^--file=", "", .file_arg[[1]]) else
    "06_postbuild_descriptives.R"
}
source(file.path(dirname(normalizePath(.script_file, mustWork = FALSE)), "00_functions.R"))

stage_postbuild_descriptives <- function() {
  require_formal_authorization("postbuild_descriptives")
  path <- file.path(DIRS$derived, "analysis_frame_pre_mice.rds")
  if (!file.exists(path)) stop("Run build stage first", call. = FALSE)
  d <- readRDS(path)
  full <- make_full_design(d)
  rows <- lapply(names(EXPOSURES), function(id) {
    v <- EXPOSURES[[id]]$variable
    natural <- full$variables$base_domain %in% TRUE & is.finite(full$variables[[v]])
    design <- full[natural, ]
    z <- design$variables[[v]]
    wm <- survey::svymean(stats::as.formula(paste0("~", v)), design, na.rm = TRUE)
    data.frame(
      exposure = id, variable = v,
      n_base = sum(full$variables$base_domain %in% TRUE),
      n_natural = length(z), n_missing_from_base = sum(full$variables$base_domain %in% TRUE) - length(z),
      min = min(z), q25 = unname(stats::quantile(z, .25)), median = stats::median(z),
      mean = mean(z), q75 = unname(stats::quantile(z, .75)), max = max(z), sd = stats::sd(z),
      weighted_mean = unname(stats::coef(wm)[[1]]), weighted_se = unname(survey::SE(wm)[[1]]),
      stringsAsFactors = FALSE
    )
  })
  write_csv_atomic(do.call(rbind, rows), file.path(DIRS$tables, "phase1_exposure_descriptives.csv"))

  pairs <- if (length(EXPOSURES) >= 2L) {
    utils::combn(names(EXPOSURES), 2, simplify = FALSE)
  } else {
    list()
  }
  correlations <- lapply(pairs, function(ids) {
    vars <- vapply(EXPOSURES[ids], `[[`, character(1), "variable")
    keep <- full$variables$base_domain %in% TRUE &
      is.finite(full$variables[[vars[[1]]]]) & is.finite(full$variables[[vars[[2]]]])
    z <- full$variables[keep, vars, drop = FALSE]
    design <- full[keep, ]
    vv <- stats::cov2cor(as.matrix(survey::svyvar(stats::reformulate(vars), design, na.rm = TRUE)))
    data.frame(
      exposure_1 = ids[[1]], exposure_2 = ids[[2]], n_common = nrow(z),
      pearson = stats::cor(z[[1]], z[[2]], method = "pearson"),
      spearman = stats::cor(z[[1]], z[[2]], method = "spearman"),
      survey_weighted_pearson = vv[1, 2], stringsAsFactors = FALSE
    )
  })
  correlation_table <- if (length(correlations)) do.call(rbind, correlations) else
    data.frame(
      exposure_1 = character(), exposure_2 = character(), n_common = integer(),
      pearson = numeric(), spearman = numeric(), survey_weighted_pearson = numeric()
    )
  write_csv_atomic(correlation_table,
                   file.path(DIRS$tables, "phase1_exposure_correlations.csv"))

  sample_rows <- list(); k <- 0L
  for (outcome in names(OUTCOMES)) for (id in names(EXPOSURES)) {
    k <- k + 1L
    domain <- OUTCOMES[[outcome]]$domain
    v <- EXPOSURES[[id]]$variable
    idx <- full$variables[[domain]] %in% TRUE & is.finite(full$variables[[v]])
    y <- full$variables[[OUTCOMES[[outcome]]$variable]][idx]
    sample_rows[[k]] <- data.frame(
      outcome = outcome, exposure = id, n_natural = sum(idx),
      events = if (OUTCOMES[[outcome]]$family == "quasibinomial") sum(y == 1) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  write_csv_atomic(do.call(rbind, sample_rows),
                   file.path(DIRS$tables, "phase1_natural_sample_counts.csv"))
  log_stage("postbuild_descriptives", "completed")
  invisible(TRUE)
}

if (sys.nframe() == 0L) stage_postbuild_descriptives()
)--------------------"

upf_source[["NHANES/02_代码/正式五结局/23_extension_common.R"]] <- r"--------------------(# Shared definitions for the prespecified first- and second-tier extensions.

.extension_common_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "NHANES/02_代码/正式五结局/23_extension_common.R"), error = function(e) NULL)
if (is.null(.extension_common_file)) .extension_common_file <- "23_extension_common.R"
EXT_CODE_DIR <- dirname(normalizePath(.extension_common_file, mustWork = FALSE))
source(file.path(EXT_CODE_DIR, "00_functions.R"))
source(file.path(EXT_CODE_DIR, "04_models.R"))

stopifnot(
  identical(PHASE1_VERSION, "2026-09-07.cholesterol-design-revision"),
  identical(names(EXPOSURES), "G_R0"),
  identical(EXPOSURE_UNIT_PP, 10)
)

EXTENSION_LABEL <- "后续第一二级分析"
EXT_DIRS <- list(
  derived = file.path(DIRS$derived, EXTENSION_LABEL),
  results = file.path(DIRS$results, EXTENSION_LABEL),
  tables = file.path(DIRS$tables, EXTENSION_LABEL),
  logs = file.path(DIRS$logs, EXTENSION_LABEL),
  tests = file.path(DIRS$tests, EXTENSION_LABEL)
)
invisible(lapply(EXT_DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

EXTENSION_VERSION <- "2026-08-25.tier1-tier2-v1"
EXPOSURE_VAR <- "upf_gram_ratio_nonwater"
RCS_PROBS <- c(0.05, 0.35, 0.65, 0.95)
RCS_REFERENCE_PROB <- 0.50
TRIM_PROBS <- c(0.01, 0.99)

effect_scale_for_outcome <- function(outcome) {
  OUTCOMES[[outcome]]$effect_scale
}

weighted_quantiles <- function(design, variable, probs) {
  f <- stats::reformulate(variable)
  z <- survey::svyquantile(f, design, quantiles = probs, ci = FALSE, na.rm = TRUE)
  out <- as.numeric(stats::coef(z))
  names(out) <- sprintf("p%02d", round(100 * probs))
  out
}

build_extension_spec <- function(d) {
  assert_columns(d, c(
    "base_domain", EXPOSURE_VAR, "analysis_weight", "SDMVSTRA", "SDMVPSU",
    "total_energy_d1", "total_energy_d2", "sex"
  ))
  full <- make_full_design(d)
  base <- full$variables$base_domain %in% TRUE &
    is.finite(full$variables[[EXPOSURE_VAR]])
  design <- full[base, ]
  all_probs <- sort(unique(c(
    TRIM_PROBS, RCS_PROBS, RCS_REFERENCE_PROB, c(0.25, 0.50, 0.75)
  )))
  q <- weighted_quantiles(design, EXPOSURE_VAR, all_probs)
  names(q) <- sprintf("p%02d", round(100 * all_probs))
  cuts <- unname(q[c("p25", "p50", "p75")])
  if (any(!is.finite(cuts)) || any(diff(cuts) <= 0)) {
    stop("Invalid weighted quartile cutpoints", call. = FALSE)
  }
  component <- cut(
    design$variables[[EXPOSURE_VAR]],
    breaks = c(-Inf, cuts, Inf), labels = paste0("Q", 1:4),
    include.lowest = TRUE, right = TRUE
  )
  design$variables$.extension_quartile <- factor(component, levels = paste0("Q", 1:4))
  medians <- vapply(paste0("Q", 1:4), function(level) {
    sub <- design[design$variables$.extension_quartile == level, ]
    weighted_quantiles(sub, EXPOSURE_VAR, 0.50)[[1]]
  }, numeric(1))
  list(
    version = EXTENSION_VERSION,
    n_base = nrow(design$variables),
    percentiles = q,
    quartile_cutpoints = cuts,
    quartile_medians = medians,
    rcs_knots = unname(q[c("p05", "p35", "p65", "p95")]),
    rcs_reference = unname(q[["p50"]]),
    trim_limits = unname(q[c("p01", "p99")]),
    energy_rule = list(
      male_low = 800, male_high = 4200,
      female_low = 500, female_high = 3500,
      application = "exclude_if_either_recall_day_outside_sex_specific_range"
    )
  )
}

extension_spec_table <- function(spec) {
  data.frame(
    item = c(
      "base_n", "quartile_p25", "quartile_p50", "quartile_p75",
      "quartile_Q1_median", "quartile_Q2_median", "quartile_Q3_median",
      "quartile_Q4_median", "rcs_knot_p05", "rcs_knot_p35",
      "rcs_knot_p65", "rcs_knot_p95", "rcs_reference_p50",
      "trim_p01", "trim_p99", "male_energy_range", "female_energy_range"
    ),
    value = c(
      spec$n_base, spec$quartile_cutpoints,
      unname(spec$quartile_medians), spec$rcs_knots, spec$rcs_reference,
      spec$trim_limits, "800-4200 kcal/day on both days",
      "500-3500 kcal/day on both days"
    ),
    stringsAsFactors = FALSE
  )
}

add_extension_exposure_fields <- function(d, spec) {
  x <- d[[EXPOSURE_VAR]]
  d$upf_quartile <- cut(
    x, breaks = c(-Inf, spec$quartile_cutpoints, Inf),
    labels = paste0("Q", 1:4), include.lowest = TRUE, right = TRUE
  )
  d$upf_quartile <- factor(d$upf_quartile, levels = paste0("Q", 1:4))
  d$upf_quartile_median10 <- ifelse(
    is.na(d$upf_quartile), NA_real_,
    unname(spec$quartile_medians[as.character(d$upf_quartile)]) / 10
  )
  d$g_r0_trimmed_domain <- d$base_domain %in% TRUE & is.finite(x) &
    x >= spec$trim_limits[[1]] & x <= spec$trim_limits[[2]]
  valid_energy <- is.finite(d$total_energy_d1) & is.finite(d$total_energy_d2)
  male <- as.character(d$sex) == "Male"
  female <- as.character(d$sex) == "Female"
  plausible <- valid_energy & (
    (male & d$total_energy_d1 >= spec$energy_rule$male_low &
       d$total_energy_d1 <= spec$energy_rule$male_high &
       d$total_energy_d2 >= spec$energy_rule$male_low &
       d$total_energy_d2 <= spec$energy_rule$male_high) |
      (female & d$total_energy_d1 >= spec$energy_rule$female_low &
         d$total_energy_d1 <= spec$energy_rule$female_high &
         d$total_energy_d2 >= spec$energy_rule$female_low &
         d$total_energy_d2 <= spec$energy_rule$female_high)
  )
  d$plausible_energy_both_days <- d$base_domain %in% TRUE & plausible
  d
}

rcs_basis <- function(x, knots) {
  knots <- as.numeric(knots) / 10
  x <- as.numeric(x) / 10
  stopifnot(length(knots) == 4L, all(diff(knots) > 0))
  pos3 <- function(z) pmax(z, 0)^3
  k1 <- knots[[1]]; k3 <- knots[[3]]; k4 <- knots[[4]]
  denom <- (k4 - k1)^2
  nonlinear <- vapply(seq_len(2L), function(j) {
    kj <- knots[[j]]
    (pos3(x - kj) - pos3(x - k3) * (k4 - kj) / (k4 - k3) +
       pos3(x - k4) * (k3 - kj) / (k4 - k3)) / denom
  }, numeric(length(x)))
  nonlinear <- matrix(nonlinear, nrow = length(x), ncol = 2L)
  out <- cbind(rcs_linear = x, rcs_nonlin1 = nonlinear[, 1], rcs_nonlin2 = nonlinear[, 2])
  rownames(out) <- NULL
  out
}

add_rcs_fields <- function(d, spec) {
  b <- rcs_basis(d[[EXPOSURE_VAR]], spec$rcs_knots)
  d$rcs_linear <- b[, "rcs_linear"]
  d$rcs_nonlin1 <- b[, "rcs_nonlin1"]
  d$rcs_nonlin2 <- b[, "rcs_nonlin2"]
  d
}

d1_joint_test <- function(summaries, terms) {
  stopifnot(length(summaries) == MI_M, length(terms) >= 1L)
  qhat <- do.call(cbind, lapply(summaries, function(s) {
    unname(s$coefficients[terms])
  }))
  uhat <- simplify2array(lapply(summaries, function(s) {
    s$variance[terms, terms, drop = FALSE]
  }))
  k <- nrow(qhat); m <- ncol(qhat)
  qbar <- rowMeans(qhat)
  ubar <- apply(uhat, c(1, 2), mean)
  b <- stats::cov(t(qhat))
  r <- (1 + 1 / m) * sum(diag(b %*% solve(ubar))) / k
  tval <- k * (m - 1)
  ttilde <- (1 + r) * ubar
  fvalue <- as.numeric(t(qbar) %*% solve(ttilde) %*% qbar / k)
  dfcom <- min(vapply(summaries, `[[`, numeric(1), "df_complete"))
  a <- r * tval / (tval - 2)
  vstar <- ((dfcom + 1) / (dfcom + 3)) * dfcom
  c0 <- 1 / (tval - 4)
  c1 <- vstar - 2 * (1 + a)
  c2 <- vstar - 4 * (1 + a)
  z <- 1 / c2 + c0 * (a^2 * c1 / ((1 + a)^2 * c2)) +
    c0 * (8 * a^2 * c1 / ((1 + a) * c2^2) + 4 * a^2 / ((1 + a) * c2)) +
    c0 * (4 * a^2 / (c2 * c1) + 16 * a^2 * c1 / c2^3) +
    c0 * (8 * a^2 / c2^2)
  df2 <- 4 + 1 / z
  data.frame(
    f_value = fvalue, df1 = k, df2 = df2,
    p_value = stats::pf(fvalue, k, df2, lower.tail = FALSE),
    relative_increase_variance = r,
    complete_data_df = dfcom,
    stringsAsFactors = FALSE
  )
}

transform_link_contrast <- function(link, se, df, effect_scale) {
  if (!is.finite(df) || df <= 0) df <- Inf
  critical <- if (is.finite(df)) stats::qt(0.975, df) else stats::qnorm(0.975)
  low <- link - critical * se
  high <- link + critical * se
  if (effect_scale == "odds_ratio") return(c(exp(link), exp(low), exp(high)))
  if (effect_scale == "log_response_percent") {
    return(100 * c(exp(link) - 1, exp(low) - 1, exp(high) - 1))
  }
  c(link, low, high)
}

invisible(TRUE)
)--------------------"

upf_source[["NHANES/02_代码/正式五结局/24_quartile_trend_models.R"]] <- r"--------------------(# Survey-weighted G-R0 quartiles and median-based trend tests for all five

.script_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "NHANES/02_代码/正式五结局/24_quartile_trend_models.R"), error = function(e) NULL)
if (is.null(.script_file)) {
  .arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .script_file <- if (length(.arg)) sub("^--file=", "", .arg[[1]]) else
    "24_quartile_trend_models.R"
}
code_dir <- dirname(normalizePath(.script_file, mustWork = FALSE))
source(file.path(code_dir, "23_extension_common.R"))
require_formal_authorization("quartile_trend_models")

d <- readRDS(file.path(DIRS$derived, "analysis_frame_pre_mice.rds"))
spec <- build_extension_spec(d)
d <- add_extension_exposure_fields(d, spec)
saveRDS(spec, file.path(EXT_DIRS$derived, "extension_analysis_spec.rds"), compress = "xz")
write_csv_atomic(
  extension_spec_table(spec),
  file.path(EXT_DIRS$tables, "extension_analysis_spec.csv")
)

quartile_formula <- function(outcome, model, type = c("categorical", "trend")) {
  type <- match.arg(type)
  exposure_term <- if (type == "categorical") "upf_quartile" else
    "upf_quartile_median10"
  stats::reformulate(
    c(exposure_term, phase1_model_covariates(model, "linear", outcome, FALSE)),
    response = OUTCOMES[[outcome]]$variable
  )
}

fit_quartile_one_imputation <- function(imputation_number, d, mi_object, outcome) {
  completed <- restore_completed_subframe(d, mi_object, imputation_number)
  full <- make_full_design(completed)
  domain <- full$variables[[OUTCOMES[[outcome]]$domain]] %in% TRUE &
    is.finite(full$variables[[EXPOSURE_VAR]])
  out <- list(categorical = list(), trend = list(), audit = list())
  for (model in names(MODEL_COVARIATES)) {
    for (type in c("categorical", "trend")) {
      form <- quartile_formula(outcome, model, type)
      model_vars <- all.vars(form)
      incomplete <- model_vars[vapply(full$variables[model_vars], function(x) {
        any(is.na(x[domain]))
      }, logical(1))]
      if (length(incomplete)) {
        stop(
          "Post-MICE missing fields in ", outcome, " ", model, " ", type,
          ": ", paste(incomplete, collapse = ", "), call. = FALSE
        )
      }
      design <- full[domain, ]
      fit <- survey::svyglm(
        form, design = design, family = outcome_family(outcome), na.action = na.fail
      )
      out[[type]][[model]] <- list(
        coefficients = stats::coef(fit), variance = stats::vcov(fit),
        df_complete = stats::df.residual(fit)
      )
      if (type == "categorical") {
        w <- as.numeric(stats::weights(design, "sampling"))
        y <- design$variables[[OUTCOMES[[outcome]]$variable]]
        out$audit[[model]] <- c(
          n_unweighted = nrow(design$variables),
          events = if (OUTCOMES[[outcome]]$family == "quasibinomial") sum(y == 1L) else NA_real_,
          weighted_population = sum(w),
          effective_sample_size = sum(w)^2 / sum(w^2),
          survey_design_df = survey::degf(design),
          complete_data_residual_df = stats::df.residual(fit)
        )
      }
      rm(fit, design)
    }
  }
  out
}

pool_quartile_outcome <- function(imputation_results, outcome) {
  rows <- list(); index <- 0L
  effect_scale <- effect_scale_for_outcome(outcome)
  for (model in names(MODEL_COVARIATES)) {
    cat_summaries <- lapply(imputation_results, function(x) x$categorical[[model]])
    trend_summaries <- lapply(imputation_results, function(x) x$trend[[model]])
    cat_pool <- pool_svy_summaries(cat_summaries)
    trend_pool <- pool_svy_summaries(trend_summaries)
    audit <- imputation_results[[1]]$audit[[model]]
    for (quartile in paste0("Q", 2:4)) {
      index <- index + 1L
      term <- paste0("upf_quartile", quartile)
      est <- pooled_term_table(cat_pool, term, effect_scale)
      rows[[index]] <- data.frame(
        version = EXTENSION_VERSION, outcome = outcome, model = model,
        model_role = unname(c(
          M1 = "unadjusted", M2 = "background_adjusted",
          M3 = "primary", M4 = "clinical_extension"
        )[model]),
        analysis = "weighted_quartile", parameter = quartile, reference = "Q1",
        effect_measure = unname(c(
          odds_ratio = "OR_vs_Q1", identity = "beta_vs_Q1",
          log_response_percent = "percent_change_vs_Q1"
        )[effect_scale]),
        estimate = est[["estimate"]], standard_error_link = est[["standard_error_link"]],
        ci_low = est[["ci_low"]], ci_high = est[["ci_high"]],
        p_value = est[["p_value"]], estimate_link = est[["estimate_link"]],
        n_unweighted = audit[["n_unweighted"]], events = audit[["events"]],
        weighted_population = audit[["weighted_population"]],
        effective_sample_size = audit[["effective_sample_size"]],
        pooled_df = unname(cat_pool$df[match(term, names(cat_pool$coefficients))]),
        successful_imputations = length(cat_summaries),
        formula = paste(deparse(quartile_formula(outcome, model, "categorical")), collapse = " "),
        covariates = paste(phase1_model_covariates(model, "linear", outcome, FALSE), collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
    index <- index + 1L
    trend_term <- "upf_quartile_median10"
    trend_est <- pooled_term_table(trend_pool, trend_term, effect_scale)
    rows[[index]] <- data.frame(
      version = EXTENSION_VERSION, outcome = outcome, model = model,
      model_role = unname(c(
        M1 = "unadjusted", M2 = "background_adjusted",
        M3 = "primary", M4 = "clinical_extension"
      )[model]),
      analysis = "weighted_quartile_trend", parameter = "weighted_median_per_10pp",
      reference = NA_character_,
      effect_measure = unname(c(
        odds_ratio = "OR_per_10pp_quartile_median", identity = "beta_per_10pp_quartile_median",
        log_response_percent = "percent_change_per_10pp_quartile_median"
      )[effect_scale]),
      estimate = trend_est[["estimate"]],
      standard_error_link = trend_est[["standard_error_link"]],
      ci_low = trend_est[["ci_low"]], ci_high = trend_est[["ci_high"]],
      p_value = trend_est[["p_value"]], estimate_link = trend_est[["estimate_link"]],
      n_unweighted = audit[["n_unweighted"]], events = audit[["events"]],
      weighted_population = audit[["weighted_population"]],
      effective_sample_size = audit[["effective_sample_size"]],
      pooled_df = unname(trend_pool$df[match(trend_term, names(trend_pool$coefficients))]),
      successful_imputations = length(trend_summaries),
      formula = paste(deparse(quartile_formula(outcome, model, "trend")), collapse = " "),
      covariates = paste(phase1_model_covariates(model, "linear", outcome, FALSE), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

all_results <- list()
workers <- n_cores_resolve()
for (outcome in names(OUTCOMES)) {
  message("QUARTILE_OUTCOME_START ", outcome)
  resource_preflight(workers, paste0("quartile_", outcome))
  mi_object <- readRDS(file.path(DIRS$derived, paste0("mi_", outcome, ".rds")))
  stopifnot(mi_object$m == MI_M)
  imputation_results <- parallel::mclapply(
    seq_len(MI_M), fit_quartile_one_imputation,
    d = d, mi_object = mi_object, outcome = outcome,
    mc.cores = workers, mc.preschedule = TRUE
  )
  failed <- vapply(imputation_results, inherits, logical(1), "try-error")
  if (any(failed)) stop("Quartile workers failed for ", outcome, call. = FALSE)
  all_results[[outcome]] <- pool_quartile_outcome(imputation_results, outcome)
  rm(mi_object, imputation_results)
  invisible(gc(full = TRUE))
  message("QUARTILE_OUTCOME_DONE ", outcome)
}

results <- do.call(rbind, all_results)
rownames(results) <- NULL
results <- results[order(
  match(results$outcome, names(OUTCOMES)),
  match(results$model, names(MODEL_COVARIATES)),
  match(results$analysis, c("weighted_quartile", "weighted_quartile_trend")),
  match(results$parameter, c("Q2", "Q3", "Q4", "weighted_median_per_10pp"))
), ]
stopifnot(
  nrow(results) == 80L,
  all(results$successful_imputations == MI_M),
  all(is.finite(results$estimate)),
  all(is.finite(results$p_value))
)
write_csv_atomic(results, file.path(EXT_DIRS$results, "weighted_quartile_trend_M1_M4.csv"))
saveRDS(results, file.path(EXT_DIRS$results, "weighted_quartile_trend_M1_M4.rds"), compress = "xz")

primary <- results[results$model == "M3", ]
write_csv_atomic(primary, file.path(EXT_DIRS$tables, "weighted_quartile_trend_M3_review.csv"))

integrity <- data.frame(
  check = c(
    "version", "base_n", "result_rows", "expected_rows",
    "all_50_imputations", "one_common_quartile_spec", "all_results_finite"
  ),
  value = c(
    EXTENSION_VERSION, spec$n_base, nrow(results), 80,
    all(results$successful_imputations == MI_M), TRUE,
    all(is.finite(results$estimate) & is.finite(results$p_value))
  ),
  stringsAsFactors = FALSE
)
write_csv_atomic(integrity, file.path(EXT_DIRS$tests, "quartile_trend_integrity.csv"))
message("QUARTILE_TREND_OK rows=", nrow(results), " base_n=", spec$n_base)

)--------------------"

upf_source[["NHANES/02_代码/正式五结局/26_sensitivity_m3_models.R"]] <- r"--------------------(# Prespecified Model-3 sensitivity analyses:

.script_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "NHANES/02_代码/正式五结局/26_sensitivity_m3_models.R"), error = function(e) NULL)
if (is.null(.script_file)) {
  .arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .script_file <- if (length(.arg)) sub("^--file=", "", .arg[[1]]) else
    "26_sensitivity_m3_models.R"
}
code_dir <- dirname(normalizePath(.script_file, mustWork = FALSE))
source(file.path(code_dir, "23_extension_common.R"))
require_formal_authorization("sensitivity_m3_models")

d <- readRDS(file.path(DIRS$derived, "analysis_frame_pre_mice.rds"))
spec <- readRDS(file.path(EXT_DIRS$derived, "extension_analysis_spec.rds"))
stopifnot(isTRUE(all.equal(spec, build_extension_spec(d), tolerance = 1e-12)))
d <- add_extension_exposure_fields(d, spec)

continuous_m3_formula <- function(outcome) {
  stats::reformulate(
    c(
      sprintf("I(%s/%s)", EXPOSURE_VAR, EXPOSURE_UNIT_PP),
      phase1_model_covariates("M3", "linear", outcome, FALSE)
    ),
    response = OUTCOMES[[outcome]]$variable
  )
}

sensitivity_domain <- function(variables, outcome, sensitivity) {
  domain <- variables[[OUTCOMES[[outcome]]$domain]] %in% TRUE &
    is.finite(variables[[EXPOSURE_VAR]])
  if (sensitivity == "exclude_extreme_energy") {
    domain <- domain & variables$plausible_energy_both_days %in% TRUE
  } else if (sensitivity == "trim_G_R0_p01_p99") {
    domain <- domain & variables$g_r0_trimmed_domain %in% TRUE
  } else {
    stop("Unknown imputed sensitivity: ", sensitivity, call. = FALSE)
  }
  domain
}

fit_sensitivity_one_imputation <- function(imputation_number, d, mi_object, outcome) {
  completed <- restore_completed_subframe(d, mi_object, imputation_number)
  full <- make_full_design(completed)
  form <- continuous_m3_formula(outcome)
  out <- list()
  for (sensitivity in c("exclude_extreme_energy", "trim_G_R0_p01_p99")) {
    domain <- sensitivity_domain(full$variables, outcome, sensitivity)
    model_vars <- all.vars(form)
    incomplete <- model_vars[vapply(full$variables[model_vars], function(x) {
      any(is.na(x[domain]))
    }, logical(1))]
    if (length(incomplete)) {
      stop("Post-MICE missing sensitivity fields: ", paste(incomplete, collapse = ", "), call. = FALSE)
    }
    design <- full[domain, ]
    fit <- survey::svyglm(
      form, design = design, family = outcome_family(outcome), na.action = na.fail
    )
    w <- as.numeric(stats::weights(design, "sampling"))
    y <- design$variables[[OUTCOMES[[outcome]]$variable]]
    out[[sensitivity]] <- list(
      coefficients = stats::coef(fit), variance = stats::vcov(fit),
      df_complete = stats::df.residual(fit),
      audit = c(
        n_unweighted = nrow(design$variables),
        events = if (OUTCOMES[[outcome]]$family == "quasibinomial") sum(y == 1L) else NA_real_,
        weighted_population = sum(w),
        effective_sample_size = sum(w)^2 / sum(w^2),
        survey_design_df = survey::degf(design),
        complete_data_residual_df = stats::df.residual(fit)
      )
    )
    rm(design, fit, w, y)
  }
  out
}

result_from_pool <- function(summaries, outcome, sensitivity) {
  pooled <- pool_svy_summaries(summaries)
  term <- sprintf("I(%s/%s)", EXPOSURE_VAR, EXPOSURE_UNIT_PP)
  estimate <- pooled_term_table(pooled, term, effect_scale_for_outcome(outcome))
  audit <- summaries[[1]]$audit
  data.frame(
    version = EXTENSION_VERSION, outcome = outcome, model = "M3",
    sensitivity = sensitivity,
    estimate = estimate[["estimate"]], ci_low = estimate[["ci_low"]],
    ci_high = estimate[["ci_high"]], p_value = estimate[["p_value"]],
    standard_error_link = estimate[["standard_error_link"]],
    estimate_link = estimate[["estimate_link"]],
    n_unweighted = audit[["n_unweighted"]], events = audit[["events"]],
    weighted_population = audit[["weighted_population"]],
    effective_sample_size = audit[["effective_sample_size"]],
    successful_imputations = length(summaries),
    formula = paste(deparse(continuous_m3_formula(outcome)), collapse = " "),
    covariates = paste(phase1_model_covariates("M3", "linear", outcome, FALSE), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

fit_complete_case <- function(d, outcome) {
  full <- make_full_design(d)
  form <- continuous_m3_formula(outcome)
  domain <- full$variables[[OUTCOMES[[outcome]]$domain]] %in% TRUE &
    is.finite(full$variables[[EXPOSURE_VAR]])
  model_vars <- all.vars(form)
  domain <- domain & stats::complete.cases(full$variables[model_vars])
  design <- full[domain, ]
  fit <- survey::svyglm(
    form, design = design, family = outcome_family(outcome), na.action = na.fail
  )
  coefficients <- stats::coef(fit); variance <- stats::vcov(fit)
  term <- sprintf("I(%s/%s)", EXPOSURE_VAR, EXPOSURE_UNIT_PP)
  idx <- match(term, names(coefficients))
  link <- unname(coefficients[[idx]])
  se <- sqrt(unname(variance[idx, idx]))
  df <- stats::df.residual(fit)
  p <- 2 * stats::pt(abs(link / se), df, lower.tail = FALSE)
  transformed <- transform_link_contrast(link, se, df, effect_scale_for_outcome(outcome))
  w <- as.numeric(stats::weights(design, "sampling"))
  y <- design$variables[[OUTCOMES[[outcome]]$variable]]
  data.frame(
    version = EXTENSION_VERSION, outcome = outcome, model = "M3",
    sensitivity = "complete_case",
    estimate = transformed[[1]], ci_low = transformed[[2]], ci_high = transformed[[3]],
    p_value = p, standard_error_link = se, estimate_link = link,
    n_unweighted = nrow(design$variables),
    events = if (OUTCOMES[[outcome]]$family == "quasibinomial") sum(y == 1L) else NA_real_,
    weighted_population = sum(w),
    effective_sample_size = sum(w)^2 / sum(w^2),
    successful_imputations = 0L,
    formula = paste(deparse(form), collapse = " "),
    covariates = paste(phase1_model_covariates("M3", "linear", outcome, FALSE), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

primary_path <- file.path(DIRS$results, "phase1_models_all_outcomes.csv")
primary <- read.csv(primary_path, stringsAsFactors = FALSE)
primary <- primary[primary$model == "M3" & primary$outcome %in% names(OUTCOMES), ]
stopifnot(nrow(primary) == 5L)
primary_rows <- data.frame(
  version = EXTENSION_VERSION, outcome = primary$outcome, model = "M3",
  sensitivity = "primary_MICE",
  estimate = primary$estimate, ci_low = primary$ci_low, ci_high = primary$ci_high,
  p_value = primary$p_value, standard_error_link = primary$standard_error,
  estimate_link = primary$estimate_link, n_unweighted = primary$n_unweighted,
  events = primary$events, weighted_population = primary$weighted_population,
  effective_sample_size = primary$effective_sample_size,
  successful_imputations = primary$successful_imputations,
  formula = primary$formula, covariates = primary$covariates,
  stringsAsFactors = FALSE
)

workers <- n_cores_resolve(); rows <- list()
for (outcome in names(OUTCOMES)) {
  message("SENSITIVITY_OUTCOME_START ", outcome)
  resource_preflight(workers, paste0("sensitivity_", outcome))
  mi_object <- readRDS(file.path(DIRS$derived, paste0("mi_", outcome, ".rds")))
  imputation_results <- parallel::mclapply(
    seq_len(MI_M), fit_sensitivity_one_imputation,
    d = d, mi_object = mi_object, outcome = outcome,
    mc.cores = workers, mc.preschedule = TRUE
  )
  failed <- vapply(imputation_results, inherits, logical(1), "try-error")
  if (any(failed)) stop("Sensitivity workers failed for ", outcome, call. = FALSE)
  for (sensitivity in c("exclude_extreme_energy", "trim_G_R0_p01_p99")) {
    summaries <- lapply(imputation_results, function(x) x[[sensitivity]])
    rows[[paste(outcome, sensitivity)]] <- result_from_pool(summaries, outcome, sensitivity)
  }
  rows[[paste(outcome, "complete_case")]] <- fit_complete_case(d, outcome)
  rm(mi_object, imputation_results)
  invisible(gc(full = TRUE))
  message("SENSITIVITY_OUTCOME_DONE ", outcome)
}

results <- rbind(primary_rows, do.call(rbind, rows))
rownames(results) <- NULL
results <- results[order(
  match(results$outcome, names(OUTCOMES)),
  match(results$sensitivity, c(
    "primary_MICE", "complete_case", "exclude_extreme_energy", "trim_G_R0_p01_p99"
  ))
), ]
stopifnot(
  nrow(results) == 20L,
  all(is.finite(results$estimate)), all(is.finite(results$p_value)),
  all(results$successful_imputations[results$sensitivity != "complete_case"] == MI_M)
)
write_csv_atomic(results, file.path(EXT_DIRS$results, "sensitivity_M3_all_outcomes.csv"))
write_csv_atomic(results, file.path(EXT_DIRS$tables, "sensitivity_M3_review.csv"))

integrity <- data.frame(
  check = c(
    "version", "rows", "expected_rows", "sensitivity_scenarios",
    "all_imputed_runs_50_of_50", "all_values_finite"
  ),
  value = c(
    EXTENSION_VERSION, nrow(results), 20,
    length(unique(results$sensitivity)),
    all(results$successful_imputations[results$sensitivity != "complete_case"] == MI_M),
    all(is.finite(results$estimate) & is.finite(results$p_value))
  ),
  stringsAsFactors = FALSE
)
write_csv_atomic(integrity, file.path(EXT_DIRS$tests, "sensitivity_integrity.csv"))
message("SENSITIVITY_M3_OK rows=", nrow(results))
)--------------------"

upf_source[["NHANES/02_代码/正式五结局/30_three_module_common.R"]] <- r"--------------------(# Shared utilities for the frozen three explanatory modules.

.three_common_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "NHANES/02_代码/正式五结局/30_three_module_common.R"), error = function(e) NULL)
if (is.null(.three_common_file)) .three_common_file <- "30_three_module_common.R"
THREE_CODE_DIR <- dirname(normalizePath(.three_common_file, mustWork = FALSE))
source(file.path(THREE_CODE_DIR, "23_extension_common.R"))

THREE_MODULE_VERSION <- "2026-08-25.three-modules-v1"
THREE_DIRS <- list(
  derived = file.path(DIRS$derived, "三模块探索"),
  results = file.path(DIRS$results, "三模块探索"),
  tables = file.path(DIRS$tables, "三模块探索"),
  logs = file.path(DIRS$logs, "三模块探索"),
  tests = file.path(DIRS$tests, "三模块探索")
)
invisible(lapply(THREE_DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

three_module_frame <- function() {
  path <- file.path(THREE_DIRS$derived, "analysis_frame_three_modules.rds")
  if (!file.exists(path)) stop("Run 29_prepare_three_module_inputs.R first", call. = FALSE)
  d <- readRDS(path)
  stopifnot(nrow(d) == 66148L, !anyDuplicated(d$SEQN))
  d
}

outcome_domain_vector <- function(d, outcome, required = character()) {
  stopifnot(outcome %in% names(OUTCOMES))
  domain <- d[[OUTCOMES[[outcome]]$domain]] %in% TRUE
  if (length(required)) {
    assert_columns(d, required)
    domain <- domain & stats::complete.cases(d[required])
  }
  domain
}

survey_summary <- function(fit, sample_ids = NULL) {
  coefficients <- stats::coef(fit)
  variance <- stats::vcov(fit)
  if (!identical(names(coefficients), rownames(variance)) ||
      !identical(names(coefficients), colnames(variance))) {
    stop("Coefficient/covariance name mismatch", call. = FALSE)
  }
  list(
    coefficients = coefficients,
    variance = variance,
    df_complete = stats::df.residual(fit),
    n_unweighted = if (is.null(sample_ids)) NA_integer_ else length(sample_ids),
    sample_hash = if (is.null(sample_ids)) NA_character_ else
      digest::digest(as.integer(sample_ids), algo = "sha256", serialize = TRUE)
  )
}

scaled_summary <- function(summary, terms, weights, name) {
  stopifnot(length(terms) == length(weights), all(terms %in% names(summary$coefficients)))
  beta <- summary$coefficients[terms]
  variance <- summary$variance[terms, terms, drop = FALSE]
  estimate <- sum(weights * beta)
  var_estimate <- as.numeric(t(weights) %*% variance %*% weights)
  list(
    coefficients = stats::setNames(estimate, name),
    variance = matrix(var_estimate, 1L, 1L, dimnames = list(name, name)),
    df_complete = summary$df_complete
  )
}

pool_contrast <- function(summaries, terms, weights, name, effect_scale) {
  transformed <- lapply(
    summaries, scaled_summary, terms = terms, weights = weights, name = name
  )
  pooled_term_table(pool_svy_summaries(transformed), name, effect_scale)
}

model_covariate_text <- function(outcome) {
  paste(phase1_model_covariates("M3", "linear", outcome, FALSE), collapse = ";")
}

model_family <- function(outcome) outcome_family(outcome)

read_module_outcome_arg <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  outcome <- if (length(args)) args[[1]] else ""
  if (!outcome %in% names(OUTCOMES)) {
    stop("Usage: Rscript SCRIPT.R {", paste(names(OUTCOMES), collapse = "|"), "}",
         call. = FALSE)
  }
  outcome
}

verify_summary_samples <- function(summaries, label) {
  hashes <- vapply(summaries, `[[`, character(1), "sample_hash")
  ns <- vapply(summaries, `[[`, integer(1), "n_unweighted")
  if (anyNA(hashes) || length(unique(hashes)) != 1L || length(unique(ns)) != 1L) {
    stop(label, " samples differ across imputations", call. = FALSE)
  }
  list(n = ns[[1]], hash = hashes[[1]])
}

invisible(TRUE)
)--------------------"

upf_source[["NHANES/02_代码/正式五结局/32_nova_substitution_models.R"]] <- r"--------------------(# Frozen NOVA1-reference statistical substitution model for one outcome.

.script_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "NHANES/02_代码/正式五结局/32_nova_substitution_models.R"), error = function(e) NULL)
if (is.null(.script_file)) {
  .arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .script_file <- if (length(.arg)) sub("^--file=", "", .arg[[1]]) else
    "32_nova_substitution_models.R"
}
code_dir <- dirname(normalizePath(.script_file, mustWork = FALSE))
source(file.path(code_dir, "30_three_module_common.R"))
require_formal_authorization("nova_substitution_models")

outcome <- read_module_outcome_arg()
d <- three_module_frame()
nova_terms <- paste0("I(nova", 2:4, "_share/10)")
form <- stats::reformulate(
  c(nova_terms, "I(total_classified_nonwater_g/1000)",
    phase1_model_covariates("M3", "linear", outcome, FALSE)),
  response = OUTCOMES[[outcome]]$variable
)

fit_one <- function(i, d, mi_object) {
  completed <- restore_completed_subframe(d, mi_object, i)
  full <- make_full_design(completed)
  required <- c(paste0("nova", 1:4, "_share"), "total_classified_nonwater_g")
  domain <- outcome_domain_vector(full$variables, outcome, required)
  design <- full[domain, ]
  fit <- survey::svyglm(form, design = design, family = model_family(outcome),
                        na.action = na.fail)
  survey_summary(fit, design$variables$SEQN)
}

mi_object <- readRDS(file.path(DIRS$derived, paste0("mi_", outcome, ".rds")))
stopifnot(mi_object$m == MI_M, identical(mi_object$outcome, outcome))
workers <- n_cores_resolve()
resource_preflight(workers, paste0("nova_substitution_", outcome))
message("NOVA_SUBSTITUTION_START outcome=", outcome, " workers=", workers)
summaries <- parallel::mclapply(
  seq_len(MI_M), fit_one, d = d, mi_object = mi_object,
  mc.cores = workers, mc.preschedule = TRUE
)
if (any(vapply(summaries, inherits, logical(1), "try-error"))) {
  stop("NOVA substitution worker failure", call. = FALSE)
}
sample <- verify_summary_samples(summaries, "NOVA substitution")
effect_scale <- effect_scale_for_outcome(outcome)
contrasts <- list(
  NOVA4_vs_NOVA1 = list(terms = nova_terms[[3]], weights = 1, role = "primary"),
  NOVA2_vs_NOVA1 = list(terms = nova_terms[[1]], weights = 1, role = "descriptive"),
  NOVA3_vs_NOVA1 = list(terms = nova_terms[[2]], weights = 1, role = "descriptive"),
  NOVA4_vs_NOVA3 = list(terms = nova_terms[c(3, 2)], weights = c(1, -1), role = "supplemental"),
  NOVA4_vs_NOVA2 = list(terms = nova_terms[c(3, 1)], weights = c(1, -1), role = "supplemental")
)
rows <- lapply(names(contrasts), function(name) {
  z <- contrasts[[name]]
  est <- pool_contrast(summaries, z$terms, z$weights, "substitution_per10pp", effect_scale)
  data.frame(
    version = THREE_MODULE_VERSION, outcome = outcome, contrast = name,
    role = z$role, model = "M3", unit = "per 10 percentage point statistical substitution",
    reference_group = if (grepl("NOVA1$", name)) "NOVA1" else
      sub("^.*_vs_", "", name),
    estimate = est[["estimate"]], ci_low = est[["ci_low"]],
    ci_high = est[["ci_high"]], p_value = est[["p_value"]],
    estimate_link = est[["estimate_link"]],
    standard_error_link = est[["standard_error_link"]],
    n_unweighted = sample$n, sample_hash = sample$hash,
    successful_imputations = MI_M,
    formula = paste(deparse(form), collapse = " "),
    covariates = model_covariate_text(outcome), stringsAsFactors = FALSE
  )
})
results <- do.call(rbind, rows)
joint <- d1_joint_test(summaries, nova_terms)
joint <- cbind(
  data.frame(version = THREE_MODULE_VERSION, outcome = outcome,
             test = "NOVA2_NOVA3_NOVA4_vs_NOVA1_joint_D1",
             n_unweighted = sample$n, successful_imputations = MI_M,
             stringsAsFactors = FALSE), joint
)
write_csv_atomic(results, file.path(
  THREE_DIRS$results, paste0("nova_substitution_", outcome, ".csv")
))
write_csv_atomic(joint, file.path(
  THREE_DIRS$results, paste0("nova_substitution_", outcome, "_joint_D1.csv")
))
message("NOVA_SUBSTITUTION_OK outcome=", outcome, " contrasts=5")

)--------------------"

upf_source[["NHANES/02_代码/正式五结局/33_diet_quality_models.R"]] <- r"--------------------(# Same-sample M3 vs M3+HEI-2015/aMED comparison for one outcome.

.script_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "NHANES/02_代码/正式五结局/33_diet_quality_models.R"), error = function(e) NULL)
if (is.null(.script_file)) {
  .arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .script_file <- if (length(.arg)) sub("^--file=", "", .arg[[1]]) else
    "33_diet_quality_models.R"
}
code_dir <- dirname(normalizePath(.script_file, mustWork = FALSE))
source(file.path(code_dir, "30_three_module_common.R"))
require_formal_authorization("diet_quality_models")

outcome <- read_module_outcome_arg()
d <- three_module_frame()
upf_term <- "I(upf_gram_ratio_nonwater/10)"
scores <- c("hei2015")
reference_formula <- stats::reformulate(
  c(upf_term, phase1_model_covariates("M3", "linear", outcome, FALSE)),
  response = OUTCOMES[[outcome]]$variable
)
adjusted_formulas <- stats::setNames(lapply(scores, function(score) {
  stats::reformulate(
    c(upf_term, score, phase1_model_covariates("M3", "linear", outcome, FALSE)),
    response = OUTCOMES[[outcome]]$variable
  )
}), scores)

fit_one <- function(i, d, mi_object) {
  completed <- restore_completed_subframe(d, mi_object, i)
  full <- make_full_design(completed)
  out <- lapply(scores, function(score) {
    domain <- outcome_domain_vector(full$variables, outcome, score)
    design <- full[domain, ]
    ref <- survey::svyglm(reference_formula, design = design,
                          family = model_family(outcome), na.action = na.fail)
    adj <- survey::svyglm(adjusted_formulas[[score]], design = design,
                          family = model_family(outcome), na.action = na.fail)
    list(
      reference = survey_summary(ref, design$variables$SEQN),
      adjusted = survey_summary(adj, design$variables$SEQN)
    )
  })
  names(out) <- scores
  out
}

mi_object <- readRDS(file.path(DIRS$derived, paste0("mi_", outcome, ".rds")))
stopifnot(mi_object$m == MI_M, identical(mi_object$outcome, outcome))
workers <- n_cores_resolve()
resource_preflight(workers, paste0("diet_quality_", outcome))
message("DIET_QUALITY_START outcome=", outcome, " workers=", workers)
fits <- parallel::mclapply(
  seq_len(MI_M), fit_one, d = d, mi_object = mi_object,
  mc.cores = workers, mc.preschedule = TRUE
)
if (any(vapply(fits, inherits, logical(1), "try-error"))) {
  stop("Diet-quality worker failure", call. = FALSE)
}
effect_scale <- effect_scale_for_outcome(outcome)
comparison_rows <- list()
score_rows <- list()
for (score in scores) {
  ref_summaries <- lapply(fits, function(x) x[[score]]$reference)
  adj_summaries <- lapply(fits, function(x) x[[score]]$adjusted)
  ref_sample <- verify_summary_samples(ref_summaries, paste(score, "reference"))
  adj_sample <- verify_summary_samples(adj_summaries, paste(score, "adjusted"))
  if (!identical(ref_sample, adj_sample)) {
    stop(score, " reference and adjusted samples differ", call. = FALSE)
  }
  ref <- pool_contrast(ref_summaries, upf_term, 1, "upf_per10pp", effect_scale)
  adj <- pool_contrast(adj_summaries, upf_term, 1, "upf_per10pp", effect_scale)
  for (variant in c("reference", "adjusted")) {
    est <- if (variant == "reference") ref else adj
    comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(
      version = THREE_MODULE_VERSION, outcome = outcome, score = score,
      variant = variant, model = if (variant == "reference") "M3" else
        paste0("M3_plus_", score),
      unit = "UPF per 10 percentage points of total non-water food grams",
      estimate = est[["estimate"]], ci_low = est[["ci_low"]],
      ci_high = est[["ci_high"]], p_value = est[["p_value"]],
      estimate_link = est[["estimate_link"]],
      standard_error_link = est[["standard_error_link"]],
      link_difference_from_reference = if (variant == "reference") 0 else
        adj[["estimate_link"]] - ref[["estimate_link"]],
      attenuation_percent = if (variant == "reference" ||
                                   abs(ref[["estimate_link"]]) < 1e-12) NA_real_ else
        100 * (ref[["estimate_link"]] - adj[["estimate_link"]]) /
          ref[["estimate_link"]],
      n_unweighted = ref_sample$n, sample_hash = ref_sample$hash,
      successful_imputations = MI_M,
      alcohol_covariate_retained = TRUE,
      formula = paste(deparse(if (variant == "reference") reference_formula else
        adjusted_formulas[[score]]), collapse = " "),
      covariates = model_covariate_text(outcome), stringsAsFactors = FALSE
    )
  }
  score_est <- pool_contrast(
    adj_summaries, score, 1, paste0(score, "_per1point"), effect_scale
  )
  score_rows[[length(score_rows) + 1L]] <- data.frame(
    version = THREE_MODULE_VERSION, outcome = outcome, score = score,
    unit = "per 1 score point", estimate = score_est[["estimate"]],
    ci_low = score_est[["ci_low"]], ci_high = score_est[["ci_high"]],
    p_value = score_est[["p_value"]],
    estimate_link = score_est[["estimate_link"]],
    standard_error_link = score_est[["standard_error_link"]],
    n_unweighted = ref_sample$n, sample_hash = ref_sample$hash,
    successful_imputations = MI_M, stringsAsFactors = FALSE
  )
}
comparison_results <- do.call(rbind, comparison_rows)
score_results <- do.call(rbind, score_rows)
write_csv_atomic(comparison_results, file.path(
  THREE_DIRS$results, paste0("diet_quality_UPF_", outcome, ".csv")
))
write_csv_atomic(score_results, file.path(
  THREE_DIRS$results, paste0("diet_quality_score_coefficients_", outcome, ".csv")
))
message("DIET_QUALITY_OK outcome=", outcome, " paired_models=2")

)--------------------"

upf_source[["NHANES/02_代码/正式五结局/41_subgroup_m50_models.R"]] <- r"--------------------(# NHANES formal m=50 five-outcome subgroup interaction analysis.

.script_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "NHANES/02_代码/正式五结局/41_subgroup_m50_models.R"), error = function(e) NULL)
if (is.null(.script_file)) {
  .file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .script_file <- if (length(.file_arg)) sub("^--file=", "", .file_arg[[1]]) else
    "41_subgroup_m50_models.R"
}
CODE_DIR <- dirname(normalizePath(.script_file, mustWork = FALSE))
source(file.path(CODE_DIR, "23_extension_common.R"))

args <- commandArgs(trailingOnly = TRUE)
RUN_MODE <- if (length(args)) tolower(args[[1]]) else ""
if (!RUN_MODE %in% c("smoke", "formal")) {
  stop("Usage: Rscript 41_subgroup_m50_models.R {smoke|formal}", call. = FALSE)
}
if (RUN_MODE == "formal") require_formal_authorization("subgroup_m50_formal")

SUBGROUP_VERSION <- "2026-08-29.nhanes-m50-five-outcome-subgroup-v1"
SUBGROUP_LEVELS <- list(
  age_group = c("20-44", "45-64", "65+"),
  sex = c("Male", "Female"),
  smoking = c("Never", "Former", "Current"),
  hypertension = c("No", "Yes"),
  diabetes = c("No", "Yes")
)
SUBGROUP_SOURCE_VARIABLE <- c(
  age_group = "age10", sex = "sex", smoking = "smoking",
  hypertension = "hypertension", diabetes = "diabetes"
)
MIN_BINARY_EVENTS <- 5L
MIN_BINARY_NONEVENTS <- 5L

subgroups_for_outcome <- function(outcome) {
  z <- names(SUBGROUP_LEVELS)
  if (identical(outcome, "DKD")) z <- setdiff(z, "diabetes")
  z
}

subgroup_field <- function(subgroup) paste0(".sg_", subgroup)

add_subgroup_fields <- function(d) {
  d[[subgroup_field("age_group")]] <- factor(
    ifelse(d$age < 45, "20-44", ifelse(d$age < 65, "45-64", "65+")),
    levels = SUBGROUP_LEVELS$age_group
  )
  d[[subgroup_field("sex")]] <- factor(
    as.character(d$sex), levels = SUBGROUP_LEVELS$sex
  )
  d[[subgroup_field("smoking")]] <- factor(
    as.character(d$smoking), levels = SUBGROUP_LEVELS$smoking
  )
  d[[subgroup_field("hypertension")]] <- factor(
    ifelse(d$hypertension == 1, "Yes", ifelse(d$hypertension == 0, "No", NA_character_)),
    levels = SUBGROUP_LEVELS$hypertension
  )
  d[[subgroup_field("diabetes")]] <- factor(
    ifelse(d$diabetes == 1, "Yes", ifelse(d$diabetes == 0, "No", NA_character_)),
    levels = SUBGROUP_LEVELS$diabetes
  )
  d
}

subgroup_formula <- function(outcome, subgroup) {
  stopifnot(outcome %in% names(OUTCOMES), subgroup %in% subgroups_for_outcome(outcome))
  exposure_term <- sprintf(
    "I(%s/%s)", EXPOSURES$G_R0$variable, EXPOSURE_UNIT_PP
  )
  covariates <- phase1_model_covariates("M3", "linear", outcome, FALSE)
  if (subgroup != "age_group") covariates <- setdiff(covariates, unname(SUBGROUP_SOURCE_VARIABLE[[subgroup]]))
  stats::reformulate(
    c(paste0(exposure_term, "*", subgroup_field(subgroup)), covariates),
    response = OUTCOMES[[outcome]]$variable
  )
}

atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary, compress = "xz")
  if (!file.rename(temporary, path)) stop("Atomic RDS rename failed: ", path, call. = FALSE)
  invisible(path)
}

task_file_path <- function(task_dir, outcome, imputation, subgroup) {
  file.path(task_dir, sprintf(
    "task_%s_imp%02d_%s.rds", outcome, as.integer(imputation), subgroup
  ))
}

support_table <- function(design, outcome, subgroup, imputation) {
  variable <- subgroup_field(subgroup)
  levels_expected <- SUBGROUP_LEVELS[[subgroup]]
  y <- design$variables[[OUTCOMES[[outcome]]$variable]]
  family_binary <- identical(OUTCOMES[[outcome]]$family, "quasibinomial")
  strata <- as.character(design$variables$.strata_cycle)
  cluster <- interaction(
    design$variables$.strata_cycle, design$variables$.psu_cycle, drop = TRUE
  )
  sampling_weights <- as.numeric(stats::weights(design, "sampling"))
  rows <- lapply(levels_expected, function(level) {
    hit <- !is.na(design$variables[[variable]]) & design$variables[[variable]] == level
    n <- sum(hit)
    events <- if (family_binary) sum(y[hit] == 1, na.rm = TRUE) else NA_integer_
    n_psu <- length(unique(cluster[hit]))
    n_strata <- length(unique(strata[hit]))
    weights_level <- sampling_weights[hit]
    data.frame(
      outcome = outcome, subgroup = subgroup, level = level,
      imputation = as.integer(imputation), n = as.integer(n),
      events = events,
      non_events = if (family_binary) n - events else NA_integer_,
      n_psu = n_psu, n_strata = n_strata, support_df = n_psu - n_strata,
      weighted_population = sum(weights_level),
      effective_sample_size = sum(weights_level)^2 / sum(weights_level^2),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  invalid <- out$n <= 0L | out$n_psu < 2L | out$support_df <= 0L
  if (family_binary) {
    invalid <- invalid | out$events < MIN_BINARY_EVENTS |
      out$non_events < MIN_BINARY_NONEVENTS
  }
  if (any(invalid)) {
    bad <- out[invalid, c(
      "level", "n", "events", "non_events", "n_psu", "n_strata", "support_df"
    ), drop = FALSE]
    stop(
      outcome, "/", subgroup, " support gate failed at imputation ", imputation,
      ": ", paste(capture.output(print(bad, row.names = FALSE)), collapse = " "),
      call. = FALSE
    )
  }
  out
}

survey_summary <- function(fit) {
  coefficients <- stats::coef(fit)
  variance <- stats::vcov(fit)
  if (!identical(names(coefficients), rownames(variance)) ||
      !identical(names(coefficients), colnames(variance)) ||
      any(!is.finite(coefficients)) || any(!is.finite(variance)) ||
      !isTRUE(all.equal(variance, t(variance), tolerance = 1e-10)) ||
      !is.finite(stats::df.residual(fit)) || stats::df.residual(fit) <= 0) {
    stop("Invalid compact survey-model summary", call. = FALSE)
  }
  list(
    coefficients = coefficients,
    variance = variance,
    df_complete = stats::df.residual(fit)
  )
}

contrast_summary <- function(summary, contrast, name = "subgroup_effect") {
  stopifnot(identical(names(contrast), names(summary$coefficients)))
  estimate <- sum(contrast * summary$coefficients)
  variance <- as.numeric(t(contrast) %*% summary$variance %*% contrast)
  list(
    coefficients = stats::setNames(estimate, name),
    variance = matrix(variance, 1L, 1L, dimnames = list(name, name)),
    df_complete = summary$df_complete
  )
}

effect_measure_label <- function(outcome) {
  unname(c(
    odds_ratio = "OR_per_10_percentage_points",
    identity = "beta_per_10_percentage_points",
    log_response_percent = "percent_change_per_10_percentage_points"
  )[[OUTCOMES[[outcome]]$effect_scale]])
}

interaction_term_names <- function(coefficient_names, subgroup) {
  exposure_term <- sprintf(
    "I(%s/%s)", EXPOSURES$G_R0$variable, EXPOSURE_UNIT_PP
  )
  sg <- subgroup_field(subgroup)
  coefficient_names[
    grepl(exposure_term, coefficient_names, fixed = TRUE) &
      grepl(sg, coefficient_names, fixed = TRUE) &
      coefficient_names != exposure_term
  ]
}

interaction_term_for_level <- function(coefficient_names, subgroup, level) {
  candidates <- interaction_term_names(coefficient_names, subgroup)
  hits <- candidates[grepl(paste0(subgroup_field(subgroup), level), candidates, fixed = TRUE)]
  if (length(hits) != 1L) {
    stop("Expected one interaction term for ", subgroup, "=", level,
         "; found ", length(hits), call. = FALSE)
  }
  hits[[1]]
}

subgroup_d1_joint_test <- function(summaries, terms) {
  stopifnot(length(summaries) == MI_M, length(terms) >= 1L)
  qhat <- do.call(cbind, lapply(summaries, function(s) {
    unname(s$coefficients[terms])
  }))
  k <- nrow(qhat)
  m <- ncol(qhat)
  uhat <- array(
    unlist(lapply(summaries, function(s) {
      s$variance[terms, terms, drop = FALSE]
    }), use.names = FALSE),
    dim = c(k, k, m)
  )
  qbar <- rowMeans(qhat)
  ubar <- matrix(apply(uhat, c(1, 2), mean), nrow = k, ncol = k)
  b <- if (k == 1L) matrix(stats::var(as.numeric(qhat)), 1L, 1L) else
    stats::cov(t(qhat))
  r <- (1 + 1 / m) * sum(diag(b %*% solve(ubar))) / k
  tval <- k * (m - 1)
  ttilde <- (1 + r) * ubar
  fvalue <- as.numeric(t(qbar) %*% solve(ttilde) %*% qbar / k)
  dfcom <- min(vapply(summaries, `[[`, numeric(1), "df_complete"))
  a <- r * tval / (tval - 2)
  vstar <- ((dfcom + 1) / (dfcom + 3)) * dfcom
  c0 <- 1 / (tval - 4)
  c1 <- vstar - 2 * (1 + a)
  c2 <- vstar - 4 * (1 + a)
  z <- 1 / c2 + c0 * (a^2 * c1 / ((1 + a)^2 * c2)) +
    c0 * (8 * a^2 * c1 / ((1 + a) * c2^2) + 4 * a^2 / ((1 + a) * c2)) +
    c0 * (4 * a^2 / (c2 * c1) + 16 * a^2 * c1 / c2^3) +
    c0 * (8 * a^2 / c2^2)
  df2 <- 4 + 1 / z
  out <- data.frame(
    f_value = fvalue, df1 = k, df2 = df2,
    p_value = stats::pf(fvalue, k, df2, lower.tail = FALSE),
    relative_increase_variance = r, complete_data_df = dfcom,
    stringsAsFactors = FALSE
  )
  if (any(!is.finite(unlist(out))) || out$df1 <= 0 || out$df2 <= 0 ||
      out$p_value < 0 || out$p_value > 1) {
    stop("D1 joint test returned an invalid result", call. = FALSE)
  }
  out
}

run_id_default <- paste0(
  if (RUN_MODE == "formal") "NHANES_SUBGROUP_M50_FORMAL_" else "NHANES_SUBGROUP_M50_SMOKE_",
  format(Sys.time(), "%Y%m%dT%H%M%S"), "_P", Sys.getpid()
)
RUN_ID <- Sys.getenv("NHANES_SUBGROUP_RUN_ID", unset = run_id_default)
if (!grepl("^[A-Za-z0-9_]+$", RUN_ID)) stop("Unsafe run ID", call. = FALSE)

analysis_label <- "亚组分析"
if (RUN_MODE == "formal") {
  RESULT_DIR <- file.path(DIRS$results, analysis_label, RUN_ID)
  QC_DIR <- file.path(DIRS$tests, analysis_label, RUN_ID)
} else {
  RESULT_DIR <- file.path(DIRS$tests, analysis_label, "SMOKE_ONLY", RUN_ID)
  QC_DIR <- RESULT_DIR
}
LOG_DIR <- file.path(DIRS$logs, analysis_label, RUN_ID)
TASK_DIR <- file.path(RESULT_DIR, "task_files")
invisible(lapply(c(RESULT_DIR, QC_DIR, LOG_DIR, TASK_DIR), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

frame_path <- file.path(DIRS$derived, "analysis_frame_pre_mice.rds")
frame <- readRDS(frame_path)
stopifnot(
  identical(PHASE1_VERSION, "2026-09-07.cholesterol-design-revision"),
  identical(MI_M, 50L), nrow(frame) == 66148L, !anyDuplicated(frame$SEQN)
)
mi_paths <- stats::setNames(
  file.path(DIRS$derived, paste0("mi_", names(OUTCOMES), ".rds")),
  names(OUTCOMES)
)
if (any(!file.exists(mi_paths))) stop("One or more formal MI objects are missing", call. = FALSE)
mi_objects <- lapply(mi_paths, readRDS)
if (any(vapply(mi_objects, function(x) x$m != MI_M, logical(1)))) {
  stop("All formal MI objects must be m=50", call. = FALSE)
}

imputation_numbers <- if (RUN_MODE == "formal") seq_len(MI_M) else 1L
task_grid <- expand.grid(
  outcome = names(OUTCOMES), imputation = imputation_numbers,
  stringsAsFactors = FALSE
)
expected_model_tasks <- sum(vapply(
  names(OUTCOMES), function(x) length(subgroups_for_outcome(x)), integer(1)
)) * length(imputation_numbers)

fit_outcome_imputation <- function(task_index) {
  task <- task_grid[task_index, , drop = FALSE]
  outcome <- task$outcome[[1]]
  imputation <- as.integer(task$imputation[[1]])
  expected_paths <- vapply(subgroups_for_outcome(outcome), function(subgroup) {
    task_file_path(TASK_DIR, outcome, imputation, subgroup)
  }, character(1))
  if (RUN_MODE == "formal" && all(file.exists(expected_paths))) {
    reusable <- vapply(seq_along(expected_paths), function(j) {
      x <- tryCatch(readRDS(expected_paths[[j]]), error = function(e) NULL)
      is.list(x) && identical(x$version, SUBGROUP_VERSION) &&
        identical(x$run_id, RUN_ID) && identical(x$run_mode, "formal") &&
        identical(x$outcome, outcome) && x$imputation == imputation &&
        identical(x$subgroup, subgroups_for_outcome(outcome)[[j]]) &&
        all(is.finite(x$summary$coefficients)) && all(is.finite(x$summary$variance))
    }, logical(1))
    if (all(reusable)) return(expected_paths)
  }
  completed <- restore_completed_subframe(frame, mi_objects[[outcome]], imputation)
  completed <- add_subgroup_fields(completed)
  full_design <- make_full_design(completed)
  output_paths <- character()
  for (subgroup in subgroups_for_outcome(outcome)) {
    form <- subgroup_formula(outcome, subgroup)
    exposure_variable <- EXPOSURES$G_R0$variable
    domain <- full_design$variables[[OUTCOMES[[outcome]]$domain]] %in% TRUE &
      is.finite(full_design$variables[[exposure_variable]])
    model_variables <- all.vars(form)
    missing_fields <- model_variables[vapply(
      full_design$variables[model_variables], function(x) any(is.na(x[domain])), logical(1)
    )]
    if (length(missing_fields)) {
      stop(outcome, "/", subgroup, " post-MICE missing fields: ",
           paste(missing_fields, collapse = ", "), call. = FALSE)
    }
    design <- full_design[domain, ]
    support <- support_table(design, outcome, subgroup, imputation)
    fit <- survey::svyglm(
      form, design = design, family = outcome_family(outcome), na.action = na.fail
    )
    summary <- survey_summary(fit)
    terms <- interaction_term_names(names(summary$coefficients), subgroup)
    if (length(terms) != length(SUBGROUP_LEVELS[[subgroup]]) - 1L) {
      stop(outcome, "/", subgroup, " interaction-term count mismatch", call. = FALSE)
    }
    object <- list(
      version = SUBGROUP_VERSION, run_id = RUN_ID, run_mode = RUN_MODE,
      outcome = outcome, subgroup = subgroup, imputation = imputation,
      required_imputations = MI_M,
      formula = paste(deparse(form), collapse = " "),
      covariates = phase1_model_covariates("M3", "linear", outcome, FALSE),
      levels = SUBGROUP_LEVELS[[subgroup]],
      sample_n = nrow(design$variables),
      sample_hash = digest::digest(
        as.integer(design$variables$SEQN), algo = "sha256", serialize = TRUE
      ),
      summary = summary, support = support
    )
    path <- task_file_path(TASK_DIR, outcome, imputation, subgroup)
    atomic_save_rds(object, path)
    output_paths <- c(output_paths, path)
    rm(fit, design, object)
  }
  rm(completed, full_design)
  invisible(gc(FALSE))
  output_paths
}

workers <- if (RUN_MODE == "formal") n_cores_resolve() else
  min(5L, nrow(task_grid), n_cores_resolve())
resource_preflight(workers, paste0("subgroup_m50_", RUN_MODE))
message(
  "SUBGROUP_RUN_START run_id=", RUN_ID, " mode=", RUN_MODE,
  " workers=", workers, " outcome_imputation_tasks=", nrow(task_grid),
  " model_tasks=", expected_model_tasks
)
task_results <- parallel::mclapply(
  seq_len(nrow(task_grid)), fit_outcome_imputation,
  mc.cores = workers, mc.preschedule = TRUE
)
failed <- vapply(task_results, inherits, logical(1), "try-error")
if (any(failed)) {
  stop("Subgroup worker failure(s): ", paste(task_results[failed], collapse = " | "),
       call. = FALSE)
}
task_files <- sort(list.files(TASK_DIR, pattern = "^task_.*[.]rds$", full.names = TRUE))
if (length(task_files) != expected_model_tasks) {
  stop("Task-file count mismatch: ", length(task_files), " vs ", expected_model_tasks,
       call. = FALSE)
}

input_lineage <- data.frame(
  role = c("analysis_frame", paste0("mi_", names(mi_paths)), "analysis_code", "frozen_plan"),
  path = c(
    frame_path, unname(mi_paths), normalizePath(.script_file, winslash = "/", mustWork = TRUE),
    file.path(DIRS$plan, "NHANES_m50五结局亚组分析冻结计划_20260829.md")
  ),
  stringsAsFactors = FALSE
)
input_lineage$sha256 <- vapply(input_lineage$path, sha256_file, character(1))
write_csv_atomic(input_lineage, file.path(QC_DIR, "input_lineage.csv"))

task_inventory <- data.frame(
  path = task_files, bytes = file.info(task_files)$size,
  sha256 = vapply(task_files, sha256_file, character(1)),
  stringsAsFactors = FALSE
)
write_csv_atomic(task_inventory, file.path(QC_DIR, "task_inventory.csv"))

if (RUN_MODE == "smoke") {
  objects <- lapply(task_files, readRDS)
  if (length(objects) != 24L || any(vapply(objects, function(x) {
    !identical(x$run_mode, "smoke") || x$imputation != 1L ||
      any(!is.finite(x$summary$coefficients)) || any(!is.finite(x$summary$variance))
  }, logical(1)))) stop("Smoke validation failed", call. = FALSE)
  writeLines(
    c(
      paste0("run_id=", RUN_ID), "status=SMOKE_OK", "imputations=1",
      paste0("model_tasks=", length(task_files)), "expected_formal_model_tasks=1200"
    ),
    file.path(RESULT_DIR, "SMOKE_OK.txt")
  )
  message("SUBGROUP_SMOKE_OK run_id=", RUN_ID, " model_tasks=", length(task_files))
  quit(save = "no", status = 0L)
}

effect_rows <- list()
interaction_rows <- list()
support_rows <- list()
effect_index <- 0L
interaction_index <- 0L
support_index <- 0L
exposure_term <- sprintf(
  "I(%s/%s)", EXPOSURES$G_R0$variable, EXPOSURE_UNIT_PP
)

for (outcome in names(OUTCOMES)) {
  for (subgroup in subgroups_for_outcome(outcome)) {
    paths <- vapply(seq_len(MI_M), function(imputation) {
      task_file_path(TASK_DIR, outcome, imputation, subgroup)
    }, character(1))
    objects <- lapply(paths, readRDS)
    if (any(vapply(seq_len(MI_M), function(i) {
      x <- objects[[i]]
      !identical(x$version, SUBGROUP_VERSION) || !identical(x$run_id, RUN_ID) ||
        !identical(x$run_mode, "formal") || !identical(x$outcome, outcome) ||
        !identical(x$subgroup, subgroup) || x$imputation != i
    }, logical(1)))) stop("Task metadata mismatch for ", outcome, "/", subgroup, call. = FALSE)
    summaries <- lapply(objects, `[[`, "summary")
    coefficient_names <- names(summaries[[1]]$coefficients)
    if (length(unique(vapply(objects, `[[`, character(1), "sample_hash"))) != 1L ||
        length(unique(vapply(objects, `[[`, integer(1), "sample_n"))) != 1L) {
      stop("Analysis samples differ across imputations for ", outcome, "/", subgroup,
           call. = FALSE)
    }
    interactions <- interaction_term_names(coefficient_names, subgroup)
    d1 <- subgroup_d1_joint_test(summaries, interactions)
    interaction_index <- interaction_index + 1L
    interaction_rows[[interaction_index]] <- data.frame(
      version = SUBGROUP_VERSION, run_id = RUN_ID, outcome = outcome,
      subgroup = subgroup, interaction_df = length(interactions),
      f_value = d1$f_value, df1 = d1$df1, df2 = d1$df2,
      p_interaction = d1$p_value,
      relative_increase_variance = d1$relative_increase_variance,
      complete_data_df = d1$complete_data_df,
      successful_imputations = length(summaries), required_imputations = MI_M,
      sample_n = objects[[1]]$sample_n, sample_hash = objects[[1]]$sample_hash,
      formula = objects[[1]]$formula, status = "OK", stringsAsFactors = FALSE
    )
    support <- do.call(rbind, lapply(objects, `[[`, "support"))
    support_index <- support_index + 1L
    support_rows[[support_index]] <- support
    for (level_index in seq_along(SUBGROUP_LEVELS[[subgroup]])) {
      level <- SUBGROUP_LEVELS[[subgroup]][[level_index]]
      contrast <- stats::setNames(rep(0, length(coefficient_names)), coefficient_names)
      contrast[[exposure_term]] <- 1
      if (level_index > 1L) {
        contrast[[interaction_term_for_level(coefficient_names, subgroup, level)]] <- 1
      }
      scalar_summaries <- lapply(
        summaries, contrast_summary, contrast = contrast, name = "subgroup_effect"
      )
      pooled <- pool_svy_summaries(scalar_summaries)
      estimate <- pooled_term_table(
        pooled, "subgroup_effect", OUTCOMES[[outcome]]$effect_scale
      )
      support_level <- support[support$level == level, , drop = FALSE]
      effect_index <- effect_index + 1L
      effect_rows[[effect_index]] <- data.frame(
        version = SUBGROUP_VERSION, run_id = RUN_ID, outcome = outcome,
        subgroup = subgroup, level = level,
        effect_measure = effect_measure_label(outcome),
        estimate = estimate[["estimate"]], ci_low = estimate[["ci_low"]],
        ci_high = estimate[["ci_high"]], p_value = estimate[["p_value"]],
        estimate_link = estimate[["estimate_link"]],
        standard_error_link = estimate[["standard_error_link"]],
        ci_low_link = estimate[["ci_low_link"]],
        ci_high_link = estimate[["ci_high_link"]],
        pooled_df = unname(pooled$df[[1]]),
        fmi_percent = unname(pooled$missinfo[[1]]),
        n_min = min(support_level$n), n_max = max(support_level$n),
        events_min = if (all(is.na(support_level$events))) NA_integer_ else
          min(support_level$events, na.rm = TRUE),
        events_max = if (all(is.na(support_level$events))) NA_integer_ else
          max(support_level$events, na.rm = TRUE),
        min_n_psu = min(support_level$n_psu),
        min_support_df = min(support_level$support_df),
        successful_imputations = length(summaries), required_imputations = MI_M,
        sample_n = objects[[1]]$sample_n, sample_hash = objects[[1]]$sample_hash,
        formula = objects[[1]]$formula, status = "OK", stringsAsFactors = FALSE
      )
    }
  }
}

effects <- do.call(rbind, effect_rows)
interactions <- do.call(rbind, interaction_rows)
support <- do.call(rbind, support_rows)
interactions$bh_q_interaction <- stats::p.adjust(interactions$p_interaction, method = "BH")
combined <- merge(
  effects,
  interactions[, c("outcome", "subgroup", "p_interaction", "bh_q_interaction")],
  by = c("outcome", "subgroup"), all.x = TRUE, sort = FALSE
)
combined <- combined[order(
  match(combined$outcome, names(OUTCOMES)),
  match(combined$subgroup, names(SUBGROUP_LEVELS)),
  ave(seq_len(nrow(combined)), combined$outcome, combined$subgroup, FUN = seq_along)
), ]
rownames(combined) <- NULL

expected_effect_rows <- 4L * sum(lengths(SUBGROUP_LEVELS)) +
  sum(lengths(SUBGROUP_LEVELS[names(SUBGROUP_LEVELS) != "diabetes"]))
expected_support_rows <- expected_effect_rows * MI_M
if (nrow(effects) != expected_effect_rows || nrow(interactions) != 24L ||
    nrow(support) != expected_support_rows ||
    any(effects$status != "OK") || any(interactions$status != "OK") ||
    any(!is.finite(effects$estimate)) || any(!is.finite(effects$ci_low)) ||
    any(!is.finite(effects$ci_high)) || any(!is.finite(effects$p_value)) ||
    any(effects$ci_low > effects$estimate | effects$estimate > effects$ci_high) ||
    any(!is.finite(interactions$p_interaction)) ||
    any(!is.finite(interactions$bh_q_interaction)) ||
    any(interactions$p_interaction < 0 | interactions$p_interaction > 1) ||
    any(interactions$bh_q_interaction < 0 | interactions$bh_q_interaction > 1) ||
    any(support$n <= 0 | support$n_psu < 2 | support$support_df <= 0) ||
    any(is.finite(support$events) &
          (support$events < MIN_BINARY_EVENTS |
             support$non_events < MIN_BINARY_NONEVENTS))) {
  stop("Formal subgroup output hard gate failed", call. = FALSE)
}

paths <- c(
  effects = file.path(RESULT_DIR, "subgroup_effects.csv"),
  interactions = file.path(RESULT_DIR, "subgroup_interactions_bh.csv"),
  combined = file.path(RESULT_DIR, "subgroup_effects_with_interactions.csv"),
  support = file.path(QC_DIR, "subgroup_support_all_imputations.csv"),
  lineage = file.path(QC_DIR, "input_lineage.csv"),
  task_inventory = file.path(QC_DIR, "task_inventory.csv")
)
write_csv_atomic(effects, paths[["effects"]])
write_csv_atomic(interactions, paths[["interactions"]])
write_csv_atomic(combined, paths[["combined"]])
write_csv_atomic(support, paths[["support"]])

run_metadata <- data.frame(
  version = SUBGROUP_VERSION, run_id = RUN_ID, mode = RUN_MODE,
  phase1_version = PHASE1_VERSION, m = MI_M, workers = workers,
  outcome_imputation_tasks = nrow(task_grid), model_tasks = length(task_files),
  effect_rows = nrow(effects), interaction_rows = nrow(interactions),
  support_rows = nrow(support), bh_family_size = nrow(interactions),
  exposure_variable = EXPOSURES$G_R0$variable,
  exposure_increment_percentage_points = EXPOSURE_UNIT_PP,
  model = "M3", created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  stringsAsFactors = FALSE
)
write_csv_atomic(run_metadata, file.path(QC_DIR, "run_metadata.csv"))

manifest_files <- c(
  unname(paths), file.path(QC_DIR, "run_metadata.csv")
)
manifest <- data.frame(
  path = manifest_files, bytes = file.info(manifest_files)$size,
  sha256 = vapply(manifest_files, sha256_file, character(1)),
  stringsAsFactors = FALSE
)
write_csv_atomic(manifest, file.path(QC_DIR, "output_manifest.csv"))
writeLines(
  c(
    paste0("run_id=", RUN_ID), "status=FORMAL_OK", paste0("m=", MI_M),
    paste0("workers=", workers), paste0("model_tasks=", length(task_files)),
    paste0("effect_rows=", nrow(effects)),
    paste0("interaction_rows=", nrow(interactions)),
    paste0("bh_family_size=", nrow(interactions))
  ),
  file.path(RESULT_DIR, "FORMAL_OK.txt")
)
latest_dir <- file.path(DIRS$logs, analysis_label)
dir.create(latest_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(RUN_ID, file.path(latest_dir, "LATEST_NHANES_SUBGROUP_M50_FORMAL.txt"))
message(
  "SUBGROUP_FORMAL_OK run_id=", RUN_ID, " model_tasks=", length(task_files),
  " effect_rows=", nrow(effects), " interaction_rows=", nrow(interactions),
  " bh_family_size=", nrow(interactions)
)
)--------------------"

upf_source[["NHANES/02_代码/00_config.R"]] <- r"--------------------(# REV22 frozen configuration: two reliable 24-hour recalls only.

REV22_VERSION <- "2026-08-14.rev22-confirmed-plan"
REV22_ROOT_DEFAULT <- "/storage/home/tmu2301/nhanes分析/UPF_全肾结局_REV22_两日膳食_最终分析_20260813"
REV22_ROOT <- Sys.getenv("NHANES_REV22_ROOT", unset = REV22_ROOT_DEFAULT)
RAW_ROOT <- Sys.getenv("NHANES_RAW_ROOT", unset = "/storage/home/tmu2301/nhanes分析/数据/nhanes数据")
LEGACY_2007_2020_ROOT <- Sys.getenv(
  "NHANES_LEGACY_ROOT",
  unset = "/storage/home/tmu2301/nhanes分析/UPF_全肾结局_窗口扩展_20260806_2007_2020"
)

DIRS <- list(
  plan = file.path(REV22_ROOT, "00_计划文档"),
  config = file.path(REV22_ROOT, "01_数据配置"),
  code = file.path(REV22_ROOT, "02_代码"),
  results = file.path(REV22_ROOT, "03_结果"),
  tables = file.path(REV22_ROOT, "04_表格"),
  figures = file.path(REV22_ROOT, "05_图表"),
  submission = file.path(REV22_ROOT, "06_投稿包"),
  logs = file.path(REV22_ROOT, "07_运行记录"),
  tests = file.path(REV22_ROOT, "08_测试")
)

FORMAL_AUTH_TOKEN <- "I_HAVE_EXPLICIT_USER_AUTHORIZATION"
formal_analysis_authorized <- function() {
  identical(Sys.getenv("NHANES_REV22_FORMAL_AUTH", unset = ""), FORMAL_AUTH_TOKEN) &&
    file.exists(file.path(DIRS$plan, "FORMAL_ANALYSIS_AUTHORIZED.txt"))
}
require_formal_authorization <- function(stage) {
  if (!formal_analysis_authorized()) {
    stop(
      "REV22 formal-analysis lock is active. Refusing stage '", stage,
      "'. Both the authorization environment token and marker file are required.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)
REV22_SEED <- 20260813L
MI_M <- 100L
MI_MAXIT <- 20L
MAX_PARALLEL_WORKERS <- 20L
DEFAULT_N_CORES <- 10L

n_cores_resolve <- function() {
  x <- suppressWarnings(as.integer(Sys.getenv(
    "NHANES_REV22_WORKERS", unset = as.character(DEFAULT_N_CORES))))
  if (!is.finite(x) || x < 1L) x <- DEFAULT_N_CORES
  as.integer(min(x, MAX_PARALLEL_WORKERS))
}

DIET_STATUS_VARS <- c(day1_status = "DR1DRSTZ", day2_status = "DR2DRSTZ",
                      number_of_recalls = "DRDINT")
DIET_WEIGHT_VARS <- c(regular = "WTDR2D", special = "WTDR2DPP")
PHYSICAL_ACTIVITY_VARS <- c("PAQ605", "PAQ620", "PAQ635", "PAQ650", "PAQ665")

CYCLES <- data.frame(
  cycle = c("E", "F", "G", "H", "I", "J", "P"),
  label = c("2007-2008", "2009-2010", "2011-2012", "2013-2014",
            "2015-2016", "2017-2018", "2017-March 2020"),
  suffix = c("e", "f", "g", "h", "i", "j", "p"),
  duration_years = c(2, 2, 2, 2, 2, 0, 3.2),
  include = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
  stringsAsFactors = FALSE
)
ANALYTIC_CYCLES <- CYCLES$cycle[CYCLES$include]

cycle_dir <- c(
  E = "2007-2008", F = "2009-2010", G = "2011-2012",
  H = "2013-2014", I = "2015-2016", J = "2017-2018",
  P = "2019-2020"
)

module_file <- function(cycle, module, stem) {
  stopifnot(cycle %in% ANALYTIC_CYCLES)
  if (cycle == "P") {
    candidates <- c(
      file.path(RAW_ROOT, cycle_dir[[cycle]], module,
                paste0("p_", tolower(stem), ".xpt")),
      file.path(RAW_ROOT, cycle_dir[[cycle]], module,
                paste0("P_", toupper(stem), ".XPT"))
    )
    hit <- candidates[file.exists(candidates)]
    return(if (length(hit)) normalizePath(hit[[1]], winslash = "/", mustWork = TRUE) else
      candidates[[1]])
  }
  file.path(RAW_ROOT, cycle_dir[[cycle]], module,
            paste0(stem, "_", tolower(cycle), ".xpt"))
}

diet_file <- function(cycle, day, kind = c("iff", "tot")) {
  stopifnot(cycle %in% ANALYTIC_CYCLES, day %in% 1:2)
  kind <- match.arg(kind)
  stem <- sprintf("dr%d%s_%s.xpt", day, kind, tolower(cycle))
  if (cycle == "P") {
    stem <- sprintf("P_DR%d%s.XPT", day, toupper(kind))
    mixed_case_stem <- sub("\\.XPT$", ".xpt", stem)
    lower_stem <- tolower(stem)
    candidates <- c(
      file.path(Sys.getenv("NHANES_OFFICIAL_XPT_DIR"), stem),
      file.path(DIRS$config, "official_xpt", stem),
      file.path(DIRS$config, "official_xpt", mixed_case_stem),
      file.path(RAW_ROOT, cycle_dir[[cycle]], "Dietary", stem),
      file.path(RAW_ROOT, cycle_dir[[cycle]], "Dietary", mixed_case_stem),
      file.path(RAW_ROOT, cycle_dir[[cycle]], "Dietary", lower_stem),
      file.path(LEGACY_2007_2020_ROOT, "01_数据", "官方XPT_2017_2020", stem)
    )
  } else {
    candidates <- c(file.path(RAW_ROOT, cycle_dir[[cycle]], "Dietary", stem))
  }
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) normalizePath(hit[[1]], winslash = "/", mustWork = TRUE) else candidates[[1]]
}

NOVA_MANIFEST_DIR <- Sys.getenv(
  "NHANES_NOVA_MANIFEST_DIR",
  unset = file.path(DIRS$config, "nova_manifest")
)
FPED_DIR <- Sys.getenv(
  "NHANES_FPED_DIR",
  unset = file.path(DIRS$config, "FPED")
)

FPED_SUFFIX <- c(E = "0708", F = "0910", G = "1112", H = "1314", I = "1516", P = "1720")
fped_file <- function(cycle, day) {
  stopifnot(cycle %in% ANALYTIC_CYCLES, day %in% 1:2)
  file.path(FPED_DIR, sprintf("fped_dr%dtot_%s.sas7bdat", day, FPED_SUFFIX[[cycle]]))
}

OUTCOMES <- list(
  CKD = list(variable = "ckd", family = "binomial", domain = "kidney_observed"),
  albuminuria = list(variable = "albuminuria", family = "binomial", domain = "acr_observed"),
  eGFR = list(variable = "egfr_2021", family = "gaussian", domain = "egfr_observed"),
  DKD = list(variable = "dkd", family = "binomial", domain = "dkd_domain"),
  hyperuricemia = list(variable = "hyperuricemia", family = "binomial", domain = "uric_acid_observed"),
  uric_acid = list(variable = "uric_acid_dxc", family = "gaussian", domain = "uric_acid_observed"),
  kidney_stones = list(variable = "kidney_stones", family = "binomial", domain = "stone_observed")
)

MODEL_COVARIATES <- list(
  M1 = c("age", "sex"),
  M2 = c("age", "sex", "race_ethnicity", "education", "pir", "smoking",
         "alcohol", "physical_activity", "cycle_factor", "mean_energy_kcal"),
  M3 = c("age", "sex", "race_ethnicity", "education", "pir", "smoking",
         "alcohol", "physical_activity", "cycle_factor", "mean_energy_kcal",
         "hypertension", "diabetes", "total_cholesterol")
)
MODEL_COVARIATES_DKD <- within(MODEL_COVARIATES, M3 <- setdiff(M3, "diabetes"))

RCS_KNOT_PROBS <- c(0.10, 0.50, 0.90)
RCS_DISPLAY_PROBS <- c(0.01, 0.99)
QUARTILE_PROBS <- c(0.25, 0.50, 0.75)
UPF_UNIT_PP <- 10
UNCLASSIFIED_GRAM_SHARE_LIMIT <- 0.05
SII_TRANSFORM <- "standardized natural log"
UPF_CATEGORY_ANALYSIS_ROLE <- "exploratory_secondary"

SUBGROUPS <- list(
  age_group = c("20-44", "45-64", "65+"),
  sex = c("Male", "Female"),
  race_ethnicity = c("Mexican American", "Other Hispanic", "Non-Hispanic White",
                     "Non-Hispanic Black", "Other/Multiracial"),
  smoking = c("Never", "Former", "Current"),
  hypertension = c("No", "Yes"),
  diabetes = c("No", "Yes")
)

required_packages <- c(
  "haven", "data.table", "survey", "mice", "mitools", "splines",
  "ggplot2", "patchwork", "broom", "openxlsx", "digest"
)

FORMAL_INPUT_CONTRACT <- c(
  "plain_water_codes.csv", "ssb_codes.csv", "upf_category_map.csv",
  "upf_category_record_map.csv",
  "laboratory_bridge.csv", "variable_map.csv", "diet_quality_manifest.csv",
  "mice_diagnostic_adjudications.csv", "formal_input_registry.csv"
)

EXPLORATORY_CONFIG_FILES <- c(
  "upf_category_map.csv", "upf_category_record_map.csv"
)

FORMAL_INPUT_REGISTRY_COLUMNS <- c(
  "logical_role", "path", "schema_version", "expected_sha256",
  "configuration_status"
)
SECOND_REVIEW_METADATA_COLUMNS <- c(
  "configuration_status"
)
FORMAL_CONFIG_SCHEMAS <- list(
  plain_water_codes.csv = c("cycle", "food_code", "plain_water_rule"),
  ssb_codes.csv = c("cycle", "food_code", "ssb_family"),
  upf_category_map.csv = c("category", "variable", "priority", "scope_nova"),
  upf_category_record_map.csv = c(
    "cycle", "food_code", "category", "variable", "priority", "assignment_rule"
  ),
  laboratory_bridge.csv = c(
    "analyte", "cycle", "source_variable", "target_variable", "equation"
  ),
  variable_map.csv = c(
    "cycle", "module", "stem", "source_variables", "required_for_rev22", "availability"
  ),
  diet_quality_manifest.csv = c(
    "record_type", "score", "cycle", "day", "relative_path", "sha256"
  ),
  mice_diagnostic_adjudications.csv = c("issue_key")
)

LABORATORY_BRIDGE_EXPECTED <- c(
  serum_creatinine_identity = "source_value",
  serum_uric_acid_identity = "source_value",
  serum_creatinine_P = "1.051 * LBXSCR - 0.06945",
  serum_uric_acid_P = "0.9323 * LBXSUA + 0.2326"
)

MICE_DIAGNOSTIC_THRESHOLDS <- list(
  minimum_iterations = MI_MAXIT,
  tail_iterations = 5L,
  max_abs_autocorrelation = 0.50,
  max_psrf = 1.10,
  max_unadjudicated_logged_events = 0L
)
MICE_DIAGNOSTIC_ADJUDICATION_FILE <- "mice_diagnostic_adjudications.csv"

STAGE_PROGRESS_UNITS <- c(
  stage_prepare_inputs = 3L,
  stage_build_dataset = 2L,
  stage_diet_quality_scores = 2L,
  stage_multiple_imputation = 5L,
  stage_main_models = 2L,
  stage_explanatory_models = 2L,
  stage_sensitivity_models = 1L,
  stage_prepare_output_contracts = 1L,
  stage_outputs = 2L
)

MICE_METHODS <- c(
  pir = "pmm", total_cholesterol = "pmm", hba1c = "pmm",
  mean_sbp = "pmm", mean_dbp = "pmm",
  education = "polyreg", smoking = "polyreg", alcohol = "polyreg",
  physical_activity = "logreg",
  diabetes_diagnosed = "logreg", insulin_current = "logreg",
  oral_med_current = "logreg", hypertension_diagnosed = "logreg",
  hypertension_med = "logreg"
)

MICE_DERIVED_COVARIATES <- c("diabetes", "hypertension")

invisible(TRUE)
)--------------------"

upf_source[["NHANES/02_代码/00_functions.R"]] <- r"--------------------(`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

.this_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "NHANES/02_代码/00_functions.R"), error = function(e) NULL) %||% "00_functions.R"
source(file.path(dirname(normalizePath(.this_file, mustWork = FALSE)), "00_config.R"), local = FALSE)

assert_columns <- function(x, cols, label = deparse(substitute(x))) {
  missing <- setdiff(cols, names(x))
  if (length(missing)) stop(label, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

rbind_fill <- function(xs) {
  cols <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs <- lapply(xs, function(x) {
    for (nm in setdiff(cols, names(x))) x[[nm]] <- NA
    x[cols]
  })
  do.call(rbind, xs)
}

yes_no_na <- function(x) {
  ifelse(x == 1, 1L, ifelse(x == 2, 0L, NA_integer_))
}

diabetes_diagnosis_binary <- function(x) {
  ifelse(x == 1, 1L, ifelse(x %in% c(2, 3), 0L, NA_integer_))
}

or_composite <- function(...) {
  z <- cbind(...)
  positive <- apply(z == 1, 1, any, na.rm = TRUE)
  all_negative <- apply(z == 0, 1, all, na.rm = FALSE) & rowSums(!is.na(z)) == ncol(z)
  ifelse(positive, 1L, ifelse(all_negative, 0L, NA_integer_))
}

derive_diabetes <- function(self_report, insulin, oral_medication, hba1c) {
  components <- list(
    yes_no_na(self_report), yes_no_na(insulin), yes_no_na(oral_medication),
    ifelse(is.na(hba1c), NA_integer_, as.integer(hba1c >= 6.5))
  )
  do.call(or_composite, components)
}

derive_diabetes_nhanes <- function(self_report, insulin, oral_medication, hba1c) {
  diagnosed <- diabetes_diagnosis_binary(self_report)
  insulin_c <- yes_no_na(insulin)
  oral_c <- yes_no_na(oral_medication)
  structural_no <- diagnosed == 0L
  insulin_c[structural_no & is.na(insulin_c)] <- 0L
  oral_c[structural_no & is.na(oral_c)] <- 0L
  lab_c <- ifelse(is.na(hba1c), NA_integer_, as.integer(hba1c >= 6.5))
  or_composite(diagnosed, insulin_c, oral_c, lab_c)
}

derive_hypertension_nhanes <- function(measured_high, diagnosed, medication) {
  dx <- yes_no_na(diagnosed)
  med <- yes_no_na(medication)
  med[dx == 0L & is.na(med)] <- 0L
  or_composite(measured_high, dx, med)
}

calc_egfr_2021 <- function(creatinine_mg_dl, age_years, female) {
  invalid <- is.na(creatinine_mg_dl) | creatinine_mg_dl <= 0 |
    is.na(age_years) | age_years < 18 | is.na(female)
  k <- ifelse(female, 0.7, 0.9)
  alpha <- ifelse(female, -0.241, -0.302)
  out <- 142 * pmin(creatinine_mg_dl / k, 1)^alpha *
    pmax(creatinine_mg_dl / k, 1)^(-1.200) *
    0.9938^age_years * ifelse(female, 1.012, 1)
  out[invalid] <- NA_real_
  out
}

normalize_equation <- function(x) gsub("[[:space:]]+", "", as.character(x))

read_laboratory_bridge <- function(path = file.path(DIRS$config, "laboratory_bridge.csv"),
                                   require_approval = FALSE) {
  if (!file.exists(path)) stop("Missing laboratory bridge: ", path, call. = FALSE)
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  assert_columns(
    x,
    c("analyte", "cycle", "source_variable", "target_variable",
      "conversion_direction", "equation", "configuration_status"),
    basename(path)
  )
  key <- paste(x$analyte, x$cycle, sep = "::")
  if (anyDuplicated(key)) stop("Duplicate analyte/cycle in laboratory bridge", call. = FALSE)
  expected_keys <- as.vector(outer(
    c("serum_creatinine", "serum_uric_acid"), ANALYTIC_CYCLES,
    paste, sep = "::"
  ))
  if (!setequal(key, expected_keys)) {
    stop("Laboratory bridge does not exactly cover both analytes in every analytic cycle",
         call. = FALSE)
  }
  identity <- x$cycle != "P"
  expected <- ifelse(
    identity,
    ifelse(x$analyte == "serum_creatinine",
           LABORATORY_BRIDGE_EXPECTED[["serum_creatinine_identity"]],
           LABORATORY_BRIDGE_EXPECTED[["serum_uric_acid_identity"]]),
    ifelse(x$analyte == "serum_creatinine",
           LABORATORY_BRIDGE_EXPECTED[["serum_creatinine_P"]],
           LABORATORY_BRIDGE_EXPECTED[["serum_uric_acid_P"]])
  )
  if (any(normalize_equation(x$equation) != normalize_equation(expected))) {
    stop("Laboratory bridge equation differs from the frozen implementation contract",
         call. = FALSE)
  }
  if (require_approval) assert_configuration_status(x, basename(path))
  x
}

assert_bridge_bound <- function(bridge, analyte) {
  if (is.null(bridge)) return(invisible(TRUE))
  assert_columns(bridge, c("analyte", "cycle", "equation"), "laboratory bridge")
  if (!all(ANALYTIC_CYCLES %in% bridge$cycle[bridge$analyte == analyte]))
    stop("Laboratory bridge was not validated for ", analyte, call. = FALSE)
  invisible(TRUE)
}

calibrate_serum_creatinine_dxc <- function(x, cycle, bridge = NULL) {
  assert_bridge_bound(bridge, "serum_creatinine")
  out <- x
  i <- cycle == "P" & !is.na(x)
  out[i] <- 1.051 * x[i] - 0.06945
  out
}

calibrate_serum_uric_acid_dxc <- function(x, cycle, bridge = NULL) {
  assert_bridge_bound(bridge, "serum_uric_acid")
  out <- x
  i <- cycle == "P" & !is.na(x)
  out[i] <- 0.9323 * x[i] + 0.2326
  out
}

derive_ckd_three_state <- function(egfr, acr) {
  positive <- (!is.na(egfr) & egfr < 60) | (!is.na(acr) & acr >= 30)
  negative <- !is.na(egfr) & egfr >= 60 & !is.na(acr) & acr < 30
  ifelse(positive, 1L, ifelse(negative, 0L, NA_integer_))
}

derive_binary_four_item <- function(...) {
  z <- cbind(...)
  yes <- apply(z == 1, 1, any, na.rm = TRUE)
  all_no <- rowSums(z == 2, na.rm = TRUE) == ncol(z) & rowSums(!is.na(z)) == ncol(z)
  ifelse(yes, 1L, ifelse(all_no, 0L, NA_integer_))
}

derive_physical_activity <- function(...) {
  z <- cbind(...)
  yes <- rowSums(z == 1, na.rm = TRUE) > 0
  all_no <- rowSums(z == 2, na.rm = TRUE) == ncol(z) & rowSums(!is.na(z)) == ncol(z)
  ifelse(yes, 1L, ifelse(all_no, 0L, NA_integer_))
}

derive_alcohol_weekly <- function(cycle, sex, alq101 = NA, alq110 = NA,
                                  alq120q = NA, alq120u = NA,
                                  alq111 = NA, alq121 = NA, alq130 = NA) {
  n <- length(cycle)
  drinks_week <- rep(NA_real_, n)
  regular <- cycle %in% c("E", "F", "G", "H", "I")
  q <- as.numeric(alq120q)
  u <- as.numeric(alq120u)
  multiplier <- ifelse(u == 1, 1, ifelse(u == 2, 1/4.345,
                 ifelse(u == 3, 1/52.1775, NA_real_)))
  valid_q <- is.finite(q) & q >= 0 & !q %in% c(777, 999)
  frequency_week <- q * multiplier
  frequency_week[!valid_q] <- NA_real_
  frequency_week[regular & valid_q & q == 0] <- 0
  p_mid <- c(`0`=0, `1`=7, `2`=6, `3`=3.5, `4`=2,
             `5`=1, `6`=2.5/4.345, `7`=1/4.345,
             `8`=9/52.1775, `9`=4.5/52.1775, `10`=1.5/52.1775)
  p <- cycle == "P" & as.character(alq121) %in% names(p_mid)
  frequency_week[p] <- unname(p_mid[as.character(alq121[p])])
  amount <- as.numeric(alq130)
  valid_amount <- is.finite(amount) & amount > 0 & !amount %in% c(777, 999)
  drinks_week[is.finite(frequency_week) & frequency_week == 0] <- 0
  drinks_week[is.finite(frequency_week) & valid_amount] <-
    frequency_week[is.finite(frequency_week) & valid_amount] *
    amount[is.finite(frequency_week) & valid_amount]
  drinks_week[regular & alq101 == 2] <- 0
  drinks_week[regular & is.na(alq101) & alq110 == 2] <- 0
  drinks_week[cycle == "P" & alq111 == 2] <- 0
  cut <- ifelse(sex == "Female", 7, ifelse(sex == "Male", 14, NA_real_))
  factor(ifelse(is.na(drinks_week) | is.na(cut), NA_character_,
                ifelse(drinks_week < 1, "None_or_very_low",
                       ifelse(drinks_week <= cut, "Low_to_moderate", "Higher"))),
         levels = c("None_or_very_low", "Low_to_moderate", "Higher"))
}

row_mean_valid <- function(x) {
  ans <- rowMeans(x, na.rm = TRUE)
  ans[rowSums(!is.na(x)) == 0] <- NA_real_
  ans
}

nhanes_bp_mean_auscultatory <- function(x, diastolic = FALSE) {
  x <- as.matrix(x)
  apply(x, 1, function(v) {
    v <- as.numeric(v)
    observed <- which(is.finite(v))
    if (!length(observed)) return(NA_real_)
    keep <- if (length(observed) == 1L) observed else observed[-1L]
    z <- v[keep]
    if (diastolic) z <- z[z != 0]
    if (!length(z)) return(NA_real_)
    mean(z)
  })
}

nhanes_bp_mean_oscillometric <- function(x, diastolic = FALSE) {
  x <- as.matrix(x)
  apply(x, 1, function(v) {
    z <- as.numeric(v[seq_len(min(length(v), 3L))])
    z <- z[is.finite(z)]
    if (!length(z)) return(NA_real_)
    if (diastolic) z <- z[z != 0]
    if (!length(z)) return(NA_real_)
    mean(z)
  })
}

nhanes_bp_mean <- nhanes_bp_mean_auscultatory

derive_sii_zln <- function(platelet, neutrophil, lymphocyte, center = NULL, scale = NULL) {
  raw <- ifelse(
    is.finite(platelet) & platelet > 0 & is.finite(neutrophil) & neutrophil > 0 &
      is.finite(lymphocyte) & lymphocyte > 0,
    platelet * neutrophil / lymphocyte, NA_real_
  )
  ln <- log(raw)
  center <- center %||% mean(ln, na.rm = TRUE)
  scale <- scale %||% stats::sd(ln, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) stop("Invalid ln(SII) scale", call. = FALSE)
  list(raw = raw, ln = ln, z = (ln - center) / scale, center = center, scale = scale)
}

combine_two_day_nova <- function(day1, day2) {
  required <- c("SEQN", "nova1_g", "nova2_g", "nova3_g", "nova4_g",
                "water_g", "unclassified_g", "ssb_nova4_g", "energy_kcal")
  assert_columns(day1, required, "day1")
  assert_columns(day2, required, "day2")
  names(day1)[-1] <- paste0(names(day1)[-1], "_d1")
  names(day2)[-1] <- paste0(names(day2)[-1], "_d2")
  z <- merge(day1, day2, by = "SEQN", all = TRUE, sort = FALSE)
  for (v in setdiff(names(z), "SEQN")) z[[v]][is.na(z[[v]])] <- 0
  sum_vars <- c("nova1_g", "nova2_g", "nova3_g", "nova4_g", "water_g",
                "unclassified_g", "ssb_nova4_g")
  sum_names <- c("nova1_2d", "nova2_2d", "nova3_2d", "nova4_2d", "water_2d",
                 "unclassified_2d", "ssb_nova4_2d")
  for (j in seq_along(sum_vars)) {
    k <- sum_vars[[j]]
    z[[sum_names[[j]]]] <- z[[paste0(k, "_d1")]] + z[[paste0(k, "_d2")]]
  }
  z$total_classified_nonwater_g <- rowSums(z[c("nova1_2d", "nova2_2d", "nova3_2d", "nova4_2d")])
  z$total_nonwater_with_unclassified_g <- z$total_classified_nonwater_g + z$unclassified_2d
  z$unclassified_share <- ifelse(z$total_nonwater_with_unclassified_g > 0,
                                 z$unclassified_2d / z$total_nonwater_with_unclassified_g, NA_real_)
  z$upf_gram_share <- ifelse(
    z$total_classified_nonwater_g > 0 & z$unclassified_share <= UNCLASSIFIED_GRAM_SHARE_LIMIT,
    100 * z$nova4_2d / z$total_classified_nonwater_g, NA_real_
  )
  z$upf_non_ssb_share <- ifelse(
    z$total_classified_nonwater_g - z$ssb_nova4_2d > 0 &
      z$ssb_nova4_2d <= z$nova4_2d & z$unclassified_share <= UNCLASSIFIED_GRAM_SHARE_LIMIT,
    100 * (z$nova4_2d - z$ssb_nova4_2d) /
      (z$total_classified_nonwater_g - z$ssb_nova4_2d), NA_real_
  )
  z
}

weighted_z <- function(x, w, domain) {
  ok <- domain %in% TRUE & is.finite(x) & is.finite(w) & w > 0
  if (sum(ok) < 2L) return(rep(NA_real_, length(x)))
  center <- sum(w[ok] * x[ok]) / sum(w[ok])
  scale <- sqrt(sum(w[ok] * (x[ok] - center)^2) / sum(w[ok]))
  if (!is.finite(scale) || scale <= 0) stop("Invalid weighted standardization scale", call. = FALSE)
  (x - center) / scale
}

assign_frozen_quartile <- function(x, cutpoints) {
  if (length(cutpoints) != 3L || any(!is.finite(cutpoints)) || is.unsorted(cutpoints))
    stop("Invalid frozen quartile cutpoints", call. = FALSE)
  factor(ifelse(is.na(x), NA_character_, ifelse(x <= cutpoints[1], "Q1",
         ifelse(x <= cutpoints[2], "Q2", ifelse(x <= cutpoints[3], "Q3", "Q4")))),
         levels = paste0("Q", 1:4))
}

canonical_path <- function(path) normalizePath(path, winslash = "/", mustWork = FALSE)

sha256_file <- function(path) {
  digest::digest(path, file = TRUE, algo = "sha256", serialize = FALSE)
}

assert_configuration_status <- function(x, label) {
  assert_columns(x, "configuration_status", label)
  if(any(is.na(x$configuration_status) | x$configuration_status != "final"))
    stop("Configuration is not the final study version: ", label, call. = FALSE)
  invisible(TRUE)
}

read_variable_contract <- function(path = file.path(DIRS$config, "variable_map.csv"),
                                   require_approval = FALSE) {
  if (!file.exists(path)) stop("Missing variable map: ", path, call. = FALSE)
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  assert_columns(
    x,
    c("cycle", "module", "stem", "source_variables", "required_for_rev22",
      "availability", "configuration_status", "version"),
    basename(path)
  )
  x <- x[x$cycle %in% ANALYTIC_CYCLES, , drop = FALSE]
  if (!all(ANALYTIC_CYCLES %in% x$cycle))
    stop("Variable map lacks one or more analytic cycles", call. = FALSE)
  if (any(tolower(x$availability) != "present"))
    stop("Required variable-map rows must explicitly be present", call. = FALSE)
  if (any(!tolower(as.character(x$required_for_rev22)) %in% c("true", "1")))
    stop("Variable map contains a non-required row in the formal contract", call. = FALSE)
  if (require_approval) assert_configuration_status(x, basename(path))
  x
}

variable_source_path <- function(cycle, module, stem) {
  if (module == "Dietary" && stem %in% c("dr1iff", "dr1tot", "dr2iff", "dr2tot")) {
    day <- as.integer(substr(stem, 3, 3))
    kind <- substr(stem, 4, 6)
    return(diet_file(cycle, day, kind))
  }
  module_file(cycle, module, stem)
}

validate_variable_source_schemas <- function(variable_map) {
  keys <- unique(variable_map[c("cycle", "module", "stem")])
  for (i in seq_len(nrow(keys))) {
    rows <- variable_map$cycle == keys$cycle[[i]] &
      variable_map$module == keys$module[[i]] & variable_map$stem == keys$stem[[i]]
    required <- unique(trimws(unlist(
      strsplit(as.character(variable_map$source_variables[rows]), ";", fixed = TRUE),
      use.names = FALSE
    )))
    path <- variable_source_path(keys$cycle[[i]], keys$module[[i]], keys$stem[[i]])
    if (!file.exists(path)) stop("Missing configured XPT: ", path, call. = FALSE)
    header <- haven::read_xpt(path, n_max = 0)
    missing <- setdiff(required, names(header))
    if (length(missing)) {
      stop("Configured XPT schema mismatch in ", path, ": ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
  }
  invisible(TRUE)
}

parse_nova_values <- function(x, label = "NOVA manifest") {
  raw <- trimws(as.character(x))
  unclassified <- is.na(x) | !nzchar(raw) | tolower(raw) %in% c("na", "n/a")
  invalid <- !unclassified & !raw %in% as.character(1:4)
  if (any(invalid)) {
    stop("Illegal non-empty NOVA value(s) in ", label, ": ",
         paste(unique(raw[invalid]), collapse = ", "), call. = FALSE)
  }
  out <- rep(NA_integer_, length(raw))
  out[!unclassified] <- as.integer(raw[!unclassified])
  out
}

validate_nova_manifest_schemas <- function() {
  for (cycle in ANALYTIC_CYCLES) {
    path <- file.path(NOVA_MANIFEST_DIR, paste0("manifest_", cycle, ".csv"))
    x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    assert_columns(x, c("code", "desc", "nova", "rule"), basename(path))
    if (anyDuplicated(as.character(x$code)))
      stop("Duplicate food code in ", path, call. = FALSE)
    parse_nova_values(x$nova, path)
  }
  invisible(TRUE)
}

expected_formal_input_paths <- function(variable_map = NULL) {
  if (is.null(variable_map)) variable_map <- read_variable_contract()
  module_keys <- unique(variable_map[c("cycle", "module", "stem")])
  module_paths <- mapply(
    variable_source_path, module_keys$cycle, module_keys$module, module_keys$stem,
    USE.NAMES = FALSE
  )
  config_paths <- file.path(
    DIRS$config,
    setdiff(FORMAL_INPUT_CONTRACT, "formal_input_registry.csv")
  )
  nova_paths <- file.path(NOVA_MANIFEST_DIR, paste0("manifest_", ANALYTIC_CYCLES, ".csv"))
  fped_paths <- unlist(lapply(ANALYTIC_CYCLES, function(cy)
    vapply(1:2, function(day) fped_file(cy, day), character(1))))
  unique(canonical_path(c(config_paths, nova_paths, fped_paths, module_paths)))
}

read_input_registry <- function(path = file.path(DIRS$config, "formal_input_registry.csv"),
                                require_approval = TRUE) {
  if (!file.exists(path)) stop("Missing formal input registry: ", path, call. = FALSE)
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  assert_columns(x, FORMAL_INPUT_REGISTRY_COLUMNS, basename(path))
  nonempty <- function(z) !is.na(z) & nzchar(trimws(as.character(z)))
  if (any(!nonempty(x$logical_role) | !nonempty(x$path) | !nonempty(x$schema_version)))
    stop("Formal input registry has an empty role, path or schema version", call. = FALSE)
  if (anyDuplicated(x$path)) stop("Duplicate path in formal input registry", call. = FALSE)
  if (any(!grepl("^[0-9a-f]{64}$", tolower(x$expected_sha256))))
    stop("Invalid expected SHA-256 in formal input registry", call. = FALSE)
  if (require_approval) assert_configuration_status(x, basename(path))
  x$resolved_path <- canonical_path(ifelse(
    grepl("^(/|[A-Za-z]:[/\\\\])", x$path), x$path, file.path(REV22_ROOT, x$path)
  ))
  x
}

validate_diet_quality_hashes <- function(path = file.path(DIRS$config, "diet_quality_manifest.csv")) {
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  assert_columns(x, c("record_type", "relative_path", "sha256"), basename(path))
  rows <- x$record_type == "data_file"
  for (i in which(rows)) {
    p <- canonical_path(file.path(DIRS$config, x$relative_path[[i]]))
    if (!file.exists(p)) stop("Missing FPED file declared by diet-quality manifest: ", p, call. = FALSE)
    if (!identical(tolower(sha256_file(p)), tolower(x$sha256[[i]])))
      stop("FPED SHA-256 mismatch: ", p, call. = FALSE)
  }
  invisible(TRUE)
}

assert_formal_inputs <- function(return_manifest = FALSE) {
  missing_packages <- required_packages[!vapply(
    required_packages, requireNamespace, logical(1), quietly = TRUE
  )]
  if (length(missing_packages)) {
    stop("Formal input preflight failed: missing R packages: ",
         paste(missing_packages, collapse = ", "), call. = FALSE)
  }
  variable_map <- read_variable_contract(require_approval = TRUE)
  registry <- read_input_registry(require_approval = TRUE)
  expected <- expected_formal_input_paths(variable_map)
  missing_registry <- setdiff(expected, registry$resolved_path)
  extra_registry <- setdiff(registry$resolved_path, expected)
  if (length(missing_registry) || length(extra_registry)) {
    stop(
      "Formal input registry path set differs from the frozen input contract",
      if (length(missing_registry)) paste0("; missing=", paste(missing_registry, collapse = ",")) else "",
      if (length(extra_registry)) paste0("; extra=", paste(extra_registry, collapse = ",")) else "",
      call. = FALSE
    )
  }
  missing_files <- registry$resolved_path[!file.exists(registry$resolved_path)]
  if (length(missing_files))
    stop("Formal input(s) missing: ", paste(missing_files, collapse = ", "), call. = FALSE)
  actual_sha <- vapply(registry$resolved_path, sha256_file, character(1))
  mismatch <- tolower(actual_sha) != tolower(registry$expected_sha256)
  if (any(mismatch))
    stop("Formal input SHA-256 mismatch: ",
         paste(registry$resolved_path[mismatch], collapse = ", "), call. = FALSE)

  review_configs <- file.path(
    DIRS$config,
    setdiff(FORMAL_INPUT_CONTRACT, "formal_input_registry.csv")
  )
  for (path in review_configs) {
    x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    required_schema <- FORMAL_CONFIG_SCHEMAS[[basename(path)]]
    if (is.null(required_schema)) stop("No frozen schema for ", basename(path), call. = FALSE)
    if (basename(path) %in% EXPLORATORY_CONFIG_FILES) {
      assert_columns(x, c(required_schema, "configuration_status"), basename(path))
      if (any(is.na(x$configuration_status) | !nzchar(trimws(as.character(x$configuration_status)))))
        stop("Exploratory configuration has an empty review status: ", basename(path),
             call. = FALSE)
    } else {
      assert_columns(x, c(required_schema, SECOND_REVIEW_METADATA_COLUMNS), basename(path))
      assert_configuration_status(x, basename(path))
    }
  }
  validate_nova_manifest_schemas()
  validate_variable_source_schemas(variable_map)
  read_laboratory_bridge()
  validate_diet_quality_hashes()

  info <- file.info(registry$resolved_path)
  manifest <- data.frame(
    logical_role = registry$logical_role,
    path = registry$resolved_path,
    schema_version = registry$schema_version,
    bytes = unname(info$size),
    mtime = format(info$mtime, "%Y-%m-%d %H:%M:%S %z"),
    expected_sha256 = tolower(registry$expected_sha256),
    actual_sha256 = tolower(actual_sha),
    configuration_status = registry$configuration_status,
    stringsAsFactors = FALSE
  )
  registry_path <- file.path(DIRS$config, "formal_input_registry.csv")
  registry_info <- file.info(registry_path)
  manifest <- rbind(manifest, data.frame(
    logical_role = "input_registry",
    path = canonical_path(registry_path),
    schema_version = "self-describing",
    bytes = unname(registry_info$size),
    mtime = format(registry_info$mtime, "%Y-%m-%d %H:%M:%S %z"),
    expected_sha256 = NA_character_,
    actual_sha256 = sha256_file(registry_path),
    configuration_status = "final",
    stringsAsFactors = FALSE
  ))
  if (return_manifest) manifest else invisible(TRUE)
}

derive_analysis_weight <- function(cycle, wtdr2d, wtdr2dpp) {
  out <- rep(NA_real_, length(cycle))
  regular <- cycle %in% c("E", "F", "G", "H", "I")
  special <- cycle == "P"
  out[regular] <- wtdr2d[regular] * 2 / 13.2
  out[special] <- wtdr2dpp[special] * 3.2 / 13.2
  out
}

base_domain <- function(d) {
  with(d, analysis_weight > 0 & age >= 20 & nonpregnant == 1 &
         diet_two_day_reliable == 1 & is.finite(upf_gram_share) &
         total_classified_nonwater_g > 0)
}

make_full_design <- function(d) {
  assert_columns(d, c("SDMVPSU", "SDMVSTRA", "analysis_weight", "cycle"))
  keep <- complete.cases(d[c("SDMVPSU", "SDMVSTRA", "analysis_weight")]) &
    is.finite(d$analysis_weight) & d$analysis_weight > 0
  if (!any(keep)) stop("No positive-weight records for the full survey design", call. = FALSE)
  dd <- d[keep, , drop = FALSE]
  dd$.strata_cycle <- interaction(dd$cycle, dd$SDMVSTRA, drop = TRUE)
  survey::svydesign(
    ids = ~SDMVPSU, strata = ~.strata_cycle, weights = ~analysis_weight,
    nest = TRUE, data = dd
  )
}

weighted_cutpoints <- function(design, variable, probs) {
  f <- stats::as.formula(paste0("~", variable))
  q <- survey::svyquantile(f, design, quantiles = probs, ci = FALSE, na.rm = TRUE)
  as.numeric(q[[1]])
}

model_formula <- function(outcome, model = "M3", exposure = "I(upf_gram_share/10)") {
  covars <- if (identical(outcome, "DKD")) MODEL_COVARIATES_DKD[[model]] else MODEL_COVARIATES[[model]]
  stats::as.formula(paste(OUTCOMES[[outcome]]$variable, "~", paste(c(exposure, covars), collapse = " + ")))
}

model_covariates <- function(outcome, model = "M3") {
  if (identical(outcome, "DKD")) MODEL_COVARIATES_DKD[[model]] else MODEL_COVARIATES[[model]]
}

outcome_domain_design <- function(full_design, outcome) {
  domain <- full_design$variables[[OUTCOMES[[outcome]]$domain]]
  domain[is.na(domain)] <- FALSE
  full_design[domain, ]
}

pool_svy_fits <- function(fits, df_complete = NULL) {
  if (length(fits) != MI_M || any(!vapply(fits, inherits, logical(1), "svyglm")))
    stop("Formal pooling requires 100/100 successful svyglm fits", call. = FALSE)
  term_sets <- lapply(fits, function(x) names(stats::coef(x)))
  if (!all(vapply(term_sets[-1], identical, logical(1), term_sets[[1]])))
    stop("Coefficient sets differ across imputations", call. = FALSE)
  df_complete <- df_complete %||% min(vapply(fits, stats::df.residual, numeric(1)))
  mitools::MIcombine(
    lapply(fits, stats::coef), lapply(fits, stats::vcov),
    df.complete = df_complete
  )
}

bh_within_family <- function(p) stats::p.adjust(p, method = "BH")

slim_svyglm <- function(fit, modifier = NULL) {
  mf <- stats::model.frame(fit)
  y <- stats::model.response(mf)
  is_binary <- all(stats::na.omit(unique(y)) %in% c(0, 1))
  des <- fit$survey.design
  w <- if (!is.null(des)) stats::weights(des, type = "sampling") else
    stats::weights(fit, type = "prior")
  w <- w[is.finite(w) & w > 0]
  modifier_levels <- NULL
  if (!is.null(modifier) && modifier %in% names(mf) && is.factor(mf[[modifier]]))
    modifier_levels <- levels(mf[[modifier]])
  f <- stats::formula(fit)
  if (!is.null(f) && !is.null(environment(f))) environment(f) <- globalenv()
  structure(
    list(
      coefficients = stats::coef(fit),
      vcov = stats::vcov(fit),
      df.residual = stats::df.residual(fit),
      formula = f,
      modifier_levels = modifier_levels,
      reliability = list(
        unweighted_n = nrow(mf),
        unweighted_events = if (is_binary) sum(y == 1, na.rm = TRUE) else NA_real_,
        design_df = if (!is.null(des)) survey::degf(des) else stats::df.residual(fit),
        effective_n = if (length(w)) sum(w)^2 / sum(w^2) else NA_real_,
        max_to_median_weight_ratio = if (length(w)) max(w) / stats::median(w) else NA_real_,
        sample_hash = digest::digest(sort(as.character(row.names(mf))), algo = "sha256"),
        converged = !identical(fit$converged, FALSE),
        finite_coefficients = all(is.finite(stats::coef(fit)))
      )
    ),
    class = c("svyglm_slim", "svyglm")
  )
}

coef.svyglm_slim <- function(object, ...) object$coefficients
vcov.svyglm_slim <- function(object, ...) object$vcov
df.residual.svyglm_slim <- function(object, ...) object$df.residual
formula.svyglm_slim <- function(x, ...) x$formula

tidy.svyglm_slim <- function(x, conf.int = FALSE, conf.level = 0.95,
                             exponentiate = FALSE, ...) {
  cf <- stats::coef(x)
  se <- sqrt(diag(stats::vcov(x)))
  df <- stats::df.residual(x)
  if (length(df) == 1L) df <- rep(df, length(cf))
  statistic <- cf / se
  p.value <- 2 * stats::pt(-abs(statistic), df)
  out <- data.frame(
    term = names(cf),
    estimate = as.numeric(cf),
    std.error = as.numeric(se),
    statistic = as.numeric(statistic),
    p.value = as.numeric(p.value),
    df.residual = df,
    nobs = x$reliability$unweighted_n,
    stringsAsFactors = FALSE
  )
  if (isTRUE(exponentiate)) {
    out$estimate <- exp(out$estimate)
    out$std.error <- out$estimate * out$std.error
    out$statistic <- out$estimate / out$std.error
  }
  out
}

invisible(TRUE)
)--------------------"

upf_source[["NHANES/02_代码/01_prepare_inputs.R"]] <- r"--------------------(# Two-day dietary exposure preparation. No model fitting occurs here.
.script_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "NHANES/02_代码/01_prepare_inputs.R"), error = function(e) NULL)
if (is.null(.script_file)) .script_file <- "01_prepare_inputs.R"
source(file.path(dirname(normalizePath(.script_file, mustWork = FALSE)), "00_functions.R"))

read_xpt_required <- function(path, columns) {
  if (!file.exists(path)) stop("Missing official input: ", path, call. = FALSE)
  x <- haven::read_xpt(path)
  assert_columns(x, columns, basename(path))
  as.data.frame(x[columns])
}

read_manifest <- function(cycle) {
  path <- file.path(NOVA_MANIFEST_DIR, paste0("manifest_", cycle, ".csv"))
  if (!file.exists(path)) stop("Missing frozen NOVA manifest: ", path, call. = FALSE)
  m <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  assert_columns(m, c("code", "desc", "nova", "rule"), basename(path))
  m$code <- as.character(m$code)
  m$nova <- parse_nova_values(m$nova, path)
  if (anyDuplicated(m$code)) stop("Duplicate food code in ", path, call. = FALSE)
  m
}

detect_iff_columns <- function(day) {
  prefix <- paste0("DR", day, "I")
  c(seqn = "SEQN", food_code = paste0(prefix, "FDCD"),
    grams = paste0(prefix, "GRMS"), energy = paste0(prefix, "KCAL"))
}

map_for_cycle <- function(map, cycle) {
  if (!"cycle" %in% names(map)) return(map)
  exact <- map[map$cycle == cycle, , drop = FALSE]
  global <- map[map$cycle %in% c("ALL", "all", "*"), , drop = FALSE]
  global <- global[!global$food_code %in% exact$food_code, , drop = FALSE]
  out <- rbind(exact, global)
  if (anyDuplicated(out$food_code))
    stop("Ambiguous effective mapping for cycle ", cycle, call. = FALSE)
  out
}

validate_code_map <- function(map, label) {
  assert_columns(map, "food_code", label)
  map$food_code <- as.character(map$food_code)
  if ("cycle" %in% names(map)) map$cycle <- as.character(map$cycle)
  if (any(is.na(map$food_code) | !nzchar(trimws(map$food_code))))
    stop(label, " contains an empty food_code", call. = FALSE)
  key <- if ("cycle" %in% names(map)) paste(map$cycle, map$food_code) else map$food_code
  if (anyDuplicated(key)) stop("Duplicate cycle/food_code in ", label, call. = FALSE)
  if ("cycle" %in% names(map)) {
    missing_cycles <- ANALYTIC_CYCLES[!vapply(ANALYTIC_CYCLES, function(cy)
      any(map$cycle %in% c(cy, "ALL", "all", "*")), logical(1))]
    if (length(missing_cycles)) stop(label, " has no mapping rows for cycle(s): ",
                                     paste(missing_cycles, collapse = ", "), call. = FALSE)
  }
  map
}

validate_record_category_map <- function(map) {
  assert_columns(map, c("cycle", "food_code", "category", "variable", "priority"),
                 "upf_category_map")
  map$cycle <- as.character(map$cycle)
  map$food_code <- as.character(map$food_code)
  map$category <- as.character(map$category)
  map$variable <- as.character(map$variable)
  map$priority <- as.integer(map$priority)
  if (any(is.na(map$food_code) | !nzchar(trimws(map$food_code))) ||
      any(is.na(map$priority)))
    stop("UPF category record map has empty code or priority", call. = FALSE)
  if (anyDuplicated(paste(map$cycle, map$food_code)))
    stop("Each cycle/food_code may have only one frozen UPF category", call. = FALSE)
  schema <- unique(map[c("category", "variable", "priority")])
  if (any(is.na(schema$category) | !nzchar(schema$category) |
          is.na(schema$variable) | !nzchar(schema$variable)) ||
      anyDuplicated(schema$category) || anyDuplicated(schema$variable) ||
      anyDuplicated(schema$priority) ||
      any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", schema$variable)))
    stop("Invalid UPF category schema", call. = FALSE)
  map
}

diet_roster_day <- function(cycle, day) {
  status <- paste0("DR", day, "DRSTZ")
  kcal <- paste0("DR", day, "TKCAL")
  cols <- c("SEQN", status, kcal, if (day == 2) "DRDINT")
  x <- read_xpt_required(diet_file(cycle, day, "tot"), cols)
  names(x) <- c("SEQN", "recall_status", "total_energy_kcal",
                if (day == 2) "number_of_recalls")
  if (anyDuplicated(x$SEQN)) stop("Duplicate SEQN in dietary total file", call. = FALSE)
  x
}

collapse_nova_day <- function(cycle, day, roster, water_codes, ssb_codes,
                              category_map) {
  cm <- detect_iff_columns(day)
  path <- diet_file(cycle, day, "iff")
  x <- haven::read_xpt(path)
  required <- unname(cm)
  assert_columns(x, required, basename(path))
  x <- as.data.frame(x[required])
  names(x) <- names(cm)
  x$food_code <- ifelse(is.na(x$food_code), NA_character_,
                        sprintf("%.0f", x$food_code))
  x$grams <- as.numeric(x$grams)
  x$energy <- as.numeric(x$energy)
  if (any(x$grams < 0, na.rm = TRUE))
    stop("Negative IFF grams in ", basename(path), call. = FALSE)
  if (any(x$energy < 0, na.rm = TRUE))
    stop("Negative IFF energy in ", basename(path), call. = FALSE)
  x$grams[!is.finite(x$grams)] <- 0
  x$energy[!is.finite(x$energy)] <- 0

  manifest <- read_manifest(cycle)
  idx <- match(x$food_code, manifest$code)
  x$nova <- manifest$nova[idx]
  x$is_water <- x$food_code %in% as.character(water_codes)
  x$is_ssb <- x$food_code %in% as.character(ssb_codes)
  if (any(x$is_ssb & x$nova != 4, na.rm = TRUE))
    stop("SSB code is not a strict NOVA4 subset in cycle ", cycle, call. = FALSE)

  cmap <- map_for_cycle(category_map, cycle)
  cidx <- match(x$food_code, cmap$food_code)
  x$category_variable <- cmap$variable[cidx]
  if (any(!is.na(x$category_variable) & x$nova != 4, na.rm = TRUE))
    stop("A frozen UPF category contains a non-NOVA4 code in cycle ", cycle, call. = FALSE)
  x$nova[x$is_water] <- NA_integer_

  base <- merge(roster, data.frame(SEQN = unique(roster$SEQN)), by = "SEQN",
                all.x = TRUE, sort = FALSE)
  base$nova1_g <- base$nova2_g <- base$nova3_g <- base$nova4_g <- 0
  base$water_g <- base$unclassified_g <- base$ssb_nova4_g <- 0
  base$energy_kcal <- 0
  base$iff_record_count <- 0L

  sum_by_person <- function(value) {
    out <- tapply(value, x$seqn, sum, na.rm = TRUE)
    unname(out[match(base$SEQN, names(out))])
  }
  if (nrow(x)) {
    for (k in 1:4) {
      v <- ifelse(!x$is_water & x$nova == k, x$grams, 0)
      base[[paste0("nova", k, "_g")]] <- replace(sum_by_person(v),
                                                   is.na(sum_by_person(v)), 0)
    }
    vectors <- list(
      water_g = ifelse(x$is_water, x$grams, 0),
      unclassified_g = ifelse(is.na(x$nova) & !x$is_water, x$grams, 0),
      ssb_nova4_g = ifelse(x$is_ssb & x$nova == 4, x$grams, 0),
      energy_kcal = x$energy,
      iff_record_count = rep(1L, nrow(x))
    )
    for (nm in names(vectors)) {
      z <- sum_by_person(vectors[[nm]])
      base[[nm]] <- replace(z, is.na(z), 0)
    }
    for (v in unique(category_map$variable)) {
      z <- sum_by_person(ifelse(x$category_variable == v & x$nova == 4, x$grams, 0))
      base[[paste0(v, "__day_g")]] <- replace(z, is.na(z), 0)
    }
  }
  nonwater <- rowSums(base[paste0("nova", 1:4, "_g")]) + base$unclassified_g
  base$no_iff_records <- base$iff_record_count == 0
  base$only_plain_water <- base$water_g > 0 & nonwater == 0
  base$total_energy_missing <- !is.finite(base$total_energy_kcal)
  base$zero_total_energy <- base$recall_status == 1 &
    is.finite(base$total_energy_kcal) & base$total_energy_kcal == 0
  base
}

combine_category_days <- function(z, d1, d2, category_map) {
  schema <- unique(category_map[c("category", "variable", "priority")])
  schema <- schema[order(schema$priority), , drop = FALSE]
  for (v in schema$variable) {
    c1 <- paste0(v, "__day_g")
    a <- d1[[c1]][match(z$SEQN, d1$SEQN)]
    b <- d2[[c1]][match(z$SEQN, d2$SEQN)]
    grams <- replace(a, is.na(a), 0) + replace(b, is.na(b), 0)
    z[[v]] <- grams
    stem <- sub("_g_2d$", "", v)
    z[[paste0(stem, "_share_total")]] <- ifelse(
      z$total_classified_nonwater_g > 0,
      100 * grams / z$total_classified_nonwater_g, NA_real_)
    z[[paste0(stem, "_share_upf")]] <- ifelse(
      z$nova4_2d > 0, 100 * grams / z$nova4_2d, NA_real_)
  }
  z$category_nova4_covered_g <- rowSums(z[schema$variable])
  z$category_nova4_unmapped_g <- pmax(z$nova4_2d - z$category_nova4_covered_g, 0)
  z
}

prepare_two_day_exposure <- function(water_map, ssb_map, category_map) {
  water_map <- validate_code_map(water_map, "plain_water_codes")
  ssb_map <- validate_code_map(ssb_map, "ssb_codes")
  category_map <- validate_record_category_map(category_map)
  pieces <- lapply(ANALYTIC_CYCLES, function(cycle) {
    r1 <- diet_roster_day(cycle, 1)
    r2 <- diet_roster_day(cycle, 2)
    roster <- merge(r1, r2, by = "SEQN", all = TRUE, suffixes = c("_d1", "_d2"), sort = FALSE)
    roster_day1 <- data.frame(SEQN = roster$SEQN, recall_status = roster$recall_status_d1,
                              total_energy_kcal = roster$total_energy_kcal_d1)
    roster_day2 <- data.frame(SEQN = roster$SEQN, recall_status = roster$recall_status_d2,
                              total_energy_kcal = roster$total_energy_kcal_d2,
                              number_of_recalls = roster$number_of_recalls)
    w <- map_for_cycle(water_map, cycle)$food_code
    s <- map_for_cycle(ssb_map, cycle)$food_code
    d1 <- collapse_nova_day(cycle, 1, roster_day1, w, s, category_map)
    d2 <- collapse_nova_day(cycle, 2, roster_day2, w, s, category_map)
    base_cols <- c("SEQN", "nova1_g", "nova2_g", "nova3_g", "nova4_g", "water_g",
                   "unclassified_g", "ssb_nova4_g", "energy_kcal")
    z <- combine_two_day_nova(d1[base_cols], d2[base_cols])
    z <- combine_category_days(z, d1, d2, category_map)
    if (any(z$category_nova4_unmapped_g > 1e-8, na.rm = TRUE))
      stop("Observed NOVA4 grams lack a frozen mutually-exclusive category in cycle ",
           cycle, call. = FALSE)
    z$cycle <- cycle
    z$day1_no_iff_records <- d1$no_iff_records[match(z$SEQN, d1$SEQN)]
    z$day2_no_iff_records <- d2$no_iff_records[match(z$SEQN, d2$SEQN)]
    z$day1_only_plain_water <- d1$only_plain_water[match(z$SEQN, d1$SEQN)]
    z$day2_only_plain_water <- d2$only_plain_water[match(z$SEQN, d2$SEQN)]
    z$day1_zero_total_energy <- d1$zero_total_energy[match(z$SEQN, d1$SEQN)]
    z$day2_zero_total_energy <- d2$zero_total_energy[match(z$SEQN, d2$SEQN)]
    z$diet_two_day_reliable_input <-
      d1$recall_status[match(z$SEQN, d1$SEQN)] == 1 &
      d2$recall_status[match(z$SEQN, d2$SEQN)] == 1 &
      d2$number_of_recalls[match(z$SEQN, d2$SEQN)] == 2
    z$two_day_zero_nova_denominator <- z$total_classified_nonwater_g <= 0
    z$unclassified_over_limit <- is.finite(z$unclassified_share) &
      z$unclassified_share > UNCLASSIFIED_GRAM_SHARE_LIMIT
    z$all_day_fasting_d1 <- z$diet_two_day_reliable_input & z$day1_zero_total_energy &
      z$day1_no_iff_records & !z$day1_only_plain_water
    z$all_day_fasting_d2 <- z$diet_two_day_reliable_input & z$day2_zero_total_energy &
      z$day2_no_iff_records & !z$day2_only_plain_water
    z$only_plain_water_d1 <- z$diet_two_day_reliable_input & z$day1_only_plain_water
    z$only_plain_water_d2 <- z$diet_two_day_reliable_input & z$day2_only_plain_water
    z
  })
  rbind_fill(pieces)
}

stage_prepare_inputs <- function() {
  require_formal_authorization("prepare_inputs")
  assert_formal_inputs()
  water_path <- file.path(DIRS$config, "plain_water_codes.csv")
  ssb_path <- file.path(DIRS$config, "ssb_codes.csv")
  category_path <- file.path(DIRS$config, "upf_category_record_map.csv")
  if (!all(file.exists(c(water_path, ssb_path, category_path))))
    stop("Frozen water, SSB and UPF-category maps are required", call. = FALSE)
  water <- utils::read.csv(water_path, stringsAsFactors = FALSE)
  ssb <- utils::read.csv(ssb_path, stringsAsFactors = FALSE)
  water <- validate_code_map(water, "plain_water_codes")
  ssb <- validate_code_map(ssb, "ssb_codes")
  for (cy in ANALYTIC_CYCLES) {
    overlap <- intersect(map_for_cycle(water, cy)$food_code,
                         map_for_cycle(ssb, cy)$food_code)
    if (length(overlap)) stop("Plain-water and SSB maps overlap in cycle ", cy,
                              call. = FALSE)
  }
  exposure <- prepare_two_day_exposure(
    water, ssb, utils::read.csv(category_path, stringsAsFactors = FALSE))
  path <- file.path(DIRS$config, "two_day_upf_exposure.rds")
  saveRDS(exposure, path, compress = "xz")
  path
}

if (sys.nframe() == 0L) stage_prepare_inputs()
)--------------------"

upf_source[["NHANES/02_代码/07_diet_quality.R"]] <- r"--------------------(# Two-day HEI-2015 and aMED construction from cycle-specific USDA FPED files.
.script_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "NHANES/02_代码/07_diet_quality.R"), error = function(e) NULL)
if (is.null(.script_file)) .script_file <- "07_diet_quality.R"
source(file.path(dirname(normalizePath(.script_file, mustWork = FALSE)), "00_functions.R"))

FPED_REQUIRED <- c(
  "FP_F_TOTAL", "FP_F_CITMLB", "FP_F_OTHER", "FP_V_TOTAL", "FP_V_DRKGR",
  "FP_V_LEGUMES", "FP_V_STARCHY_POTATO", "FP_G_WHOLE", "FP_G_REFINED",
  "FP_D_TOTAL", "FP_PF_MEAT", "FP_PF_CUREDMEAT", "FP_PF_ORGAN",
  "FP_PF_POULT", "FP_PF_SEAFD_HI", "FP_PF_SEAFD_LOW", "FP_PF_EGGS",
  "FP_PF_SOY", "FP_PF_NUTSDS", "FP_PF_LEGUMES", "FP_ADD_SUGARS"
)

HEI2015_COMPONENT_MAX <- c(
  total_fruits = 5, whole_fruits = 5, total_vegetables = 5,
  greens_beans = 5, whole_grains = 10, dairy = 10,
  total_protein = 5, seafood_plant = 5, fatty_acids = 10,
  refined_grains = 10, sodium = 10, added_sugars = 10,
  saturated_fats = 10
)

AMED_COMPONENTS <- c(
  "vegetables", "legumes", "fruits", "nuts", "whole_grains", "fish",
  "mufa_sfa", "red_processed_meat", "alcohol"
)

read_fped_day <- function(cycle, day) {
  path <- fped_file(cycle, day)
  if (!file.exists(path)) stop("Missing USDA FPED: ", path, call. = FALSE)
  x <- as.data.frame(haven::read_sas(path))
  names(x) <- toupper(names(x))
  prefix <- paste0("DR", day, "T_")
  names(x) <- sub(paste0("^", prefix), "FP_", names(x))
  assert_columns(x, c("SEQN", FPED_REQUIRED), basename(path))
  if (anyDuplicated(x$SEQN)) stop("Duplicate SEQN in ", path, call. = FALSE)
  x$cycle <- cycle
  x
}

combine_fped_two_day <- function(d1, d2) {
  assert_columns(d1, c("SEQN", FPED_REQUIRED), "FPED day 1")
  assert_columns(d2, c("SEQN", FPED_REQUIRED), "FPED day 2")
  if (anyDuplicated(d1$SEQN) || anyDuplicated(d2$SEQN))
    stop("FPED day files must contain one row per SEQN", call. = FALSE)
  cols <- intersect(grep("^FP_", names(d1), value = TRUE),
                    grep("^FP_", names(d2), value = TRUE))
  z <- merge(d1[c("SEQN", cols)], d2[c("SEQN", cols)], by = "SEQN", all = TRUE,
             suffixes = c("_d1", "_d2"), sort = FALSE)
  z$fped_day1_present <- z$SEQN %in% d1$SEQN
  z$fped_day2_present <- z$SEQN %in% d2$SEQN
  for (v in cols) {
    a <- as.numeric(z[[paste0(v, "_d1")]])
    b <- as.numeric(z[[paste0(v, "_d2")]])
    z[[paste0(v, "_2d")]] <- ifelse(is.finite(a) & is.finite(b), a + b, NA_real_)
  }
  z$fped_two_day_complete <- z$fped_day1_present & z$fped_day2_present &
    stats::complete.cases(z[paste0(FPED_REQUIRED, "_2d")])
  z
}

score_adequacy <- function(x, maximum_standard, points) {
  out <- points * pmin(pmax(x / maximum_standard, 0), 1)
  out[!is.finite(x)] <- NA_real_
  out
}

score_moderation <- function(x, maximum_standard, zero_standard, points = 10) {
  out <- points * pmin(pmax((zero_standard - x) /
                              (zero_standard - maximum_standard), 0), 1)
  out[!is.finite(x)] <- NA_real_
  out
}

fped_two_day_constituents <- function(fped) {
  req <- paste0(FPED_REQUIRED, "_2d")
  assert_columns(fped, req, "two-day FPED")
  getv <- function(x) as.numeric(fped[[paste0(x, "_2d")]])
  data.frame(
    total_fruits = getv("FP_F_TOTAL"),
    whole_fruits = getv("FP_F_CITMLB") + getv("FP_F_OTHER"),
    total_vegetables = getv("FP_V_TOTAL") + getv("FP_V_LEGUMES"),
    greens_beans = getv("FP_V_DRKGR") + getv("FP_V_LEGUMES"),
    whole_grains = getv("FP_G_WHOLE"),
    refined_grains = getv("FP_G_REFINED"),
    dairy = getv("FP_D_TOTAL"),
    total_protein = getv("FP_PF_MEAT") + getv("FP_PF_CUREDMEAT") +
      getv("FP_PF_ORGAN") + getv("FP_PF_POULT") + getv("FP_PF_SEAFD_HI") +
      getv("FP_PF_SEAFD_LOW") + getv("FP_PF_EGGS") + getv("FP_PF_SOY") +
      getv("FP_PF_NUTSDS") + getv("FP_PF_LEGUMES"),
    seafood_plant = getv("FP_PF_SEAFD_HI") + getv("FP_PF_SEAFD_LOW") +
      getv("FP_PF_SOY") + getv("FP_PF_NUTSDS") + getv("FP_PF_LEGUMES"),
    added_sugars_tsp = getv("FP_ADD_SUGARS"),
    stringsAsFactors = FALSE
  )
}

hei2015_two_day <- function(fped, energy_kcal_2d, sodium_mg_2d,
                            satfat_g_2d, mufa_g_2d, pufa_g_2d) {
  n <- nrow(fped)
  inputs <- list(energy_kcal_2d, sodium_mg_2d, satfat_g_2d, mufa_g_2d, pufa_g_2d)
  if (any(vapply(inputs, length, integer(1)) != n))
    stop("HEI nutrient vectors must match FPED rows", call. = FALSE)
  cst <- fped_two_day_constituents(fped)
  energy <- as.numeric(energy_kcal_2d)
  e1000 <- energy / 1000
  satfat <- as.numeric(satfat_g_2d)
  mufa <- as.numeric(mufa_g_2d)
  pufa <- as.numeric(pufa_g_2d)
  valid <- is.finite(energy) & energy > 0 & is.finite(sodium_mg_2d) &
    is.finite(satfat) & satfat >= 0 & is.finite(mufa) & mufa >= 0 &
    is.finite(pufa) & pufa >= 0 & stats::complete.cases(cst)
  fatty_ratio <- ifelse(satfat > 0, (mufa + pufa) / satfat,
                        ifelse(satfat == 0 & (mufa + pufa) > 0, Inf, NA_real_))
  fatty_score <- 10 * pmin(pmax((fatty_ratio - 1.2) / (2.5 - 1.2), 0), 1)
  fatty_score[is.na(fatty_ratio)] <- NA_real_
  scores <- data.frame(
    hei_total_fruits = score_adequacy(cst$total_fruits / e1000, 0.8, 5),
    hei_whole_fruits = score_adequacy(cst$whole_fruits / e1000, 0.4, 5),
    hei_total_vegetables = score_adequacy(cst$total_vegetables / e1000, 1.1, 5),
    hei_greens_beans = score_adequacy(cst$greens_beans / e1000, 0.2, 5),
    hei_whole_grains = score_adequacy(cst$whole_grains / e1000, 1.5, 10),
    hei_dairy = score_adequacy(cst$dairy / e1000, 1.3, 10),
    hei_total_protein = score_adequacy(cst$total_protein / e1000, 2.5, 5),
    hei_seafood_plant = score_adequacy(cst$seafood_plant / e1000, 0.8, 5),
    hei_fatty_acids = fatty_score,
    hei_refined_grains = score_moderation(cst$refined_grains / e1000, 1.8, 4.3),
    hei_sodium = score_moderation(as.numeric(sodium_mg_2d) / 1000 / e1000, 1.1, 2.0),
    hei_added_sugars = score_moderation(100 * 16 * cst$added_sugars_tsp / energy, 6.5, 26),
    hei_saturated_fats = score_moderation(100 * 9 * satfat / energy, 8, 16),
    stringsAsFactors = FALSE
  )
  scores$hei2015 <- rowSums(scores, na.rm = FALSE)
  scores[!valid, ] <- NA_real_
  if (any(scores$hei2015 < 0 | scores$hei2015 > 100, na.rm = TRUE))
    stop("HEI-2015 score outside 0--100", call. = FALSE)
  scores
}

weighted_median <- function(x, w, domain = rep(TRUE, length(x))) {
  ok <- domain %in% TRUE & is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  ord <- order(x[ok])
  xx <- x[ok][ord]; ww <- w[ok][ord]
  xx[which(cumsum(ww) >= sum(ww) / 2)[1]]
}

amed_constituents_two_day <- function(fped, satfat_g_2d, mufa_g_2d,
                                      alcohol_g_2d) {
  req <- paste0(FPED_REQUIRED, "_2d")
  assert_columns(fped, req, "two-day FPED")
  n <- nrow(fped)
  if (any(c(length(satfat_g_2d), length(mufa_g_2d), length(alcohol_g_2d)) != n))
    stop("aMED nutrient vectors must match FPED rows", call. = FALSE)
  getv <- function(x) as.numeric(fped[[paste0(x, "_2d")]]) / 2
  sat <- as.numeric(satfat_g_2d) / 2
  data.frame(
    vegetables = pmax(getv("FP_V_TOTAL") - getv("FP_V_STARCHY_POTATO"), 0),
    legumes = getv("FP_V_LEGUMES"), fruits = getv("FP_F_TOTAL"),
    nuts = getv("FP_PF_NUTSDS"), whole_grains = getv("FP_G_WHOLE"),
    fish = getv("FP_PF_SEAFD_HI") + getv("FP_PF_SEAFD_LOW"),
    mufa_sfa = ifelse(sat > 0, (as.numeric(mufa_g_2d) / 2) / sat, NA_real_),
    red_processed_meat = getv("FP_PF_MEAT") + getv("FP_PF_CUREDMEAT"),
    alcohol = as.numeric(alcohol_g_2d) / 2,
    stringsAsFactors = FALSE
  )
}

amed_weighted_medians <- function(constituents, sex, weight, domain) {
  assert_columns(constituents, AMED_COMPONENTS, "aMED constituents")
  sex <- as.character(sex)
  out <- do.call(rbind, lapply(c("Male", "Female"), function(s) {
    idx <- domain %in% TRUE & sex == s
    vals <- vapply(setdiff(AMED_COMPONENTS, "alcohol"), function(v)
      weighted_median(constituents[[v]], weight, idx), numeric(1))
    data.frame(sex = s, component = names(vals), median = unname(vals),
               stringsAsFactors = FALSE)
  }))
  if (any(!is.finite(out$median))) stop("Unable to freeze aMED weighted medians", call. = FALSE)
  out
}

amed_score <- function(constituents, sex, weight, domain,
                       medians = amed_weighted_medians(constituents, sex, weight, domain)) {
  assert_columns(constituents, AMED_COMPONENTS, "aMED constituents")
  assert_columns(medians, c("sex", "component", "median"), "aMED medians")
  sex <- as.character(sex)
  n <- nrow(constituents)
  out <- data.frame(matrix(NA_integer_, nrow = n, ncol = length(AMED_COMPONENTS)))
  names(out) <- paste0("amed_", AMED_COMPONENTS)
  beneficial <- setdiff(AMED_COMPONENTS, c("red_processed_meat", "alcohol"))
  for (v in beneficial) {
    cut <- medians$median[match(paste(sex, v), paste(medians$sex, medians$component))]
    out[[paste0("amed_", v)]] <- ifelse(is.finite(constituents[[v]]) & is.finite(cut),
                                         as.integer(constituents[[v]] > cut), NA_integer_)
  }
  cut <- medians$median[match(paste(sex, "red_processed_meat"),
                              paste(medians$sex, medians$component))]
  out$amed_red_processed_meat <- ifelse(
    is.finite(constituents$red_processed_meat) & is.finite(cut),
    as.integer(constituents$red_processed_meat < cut), NA_integer_)
  alc_ok <- is.finite(constituents$alcohol) &
    ((sex == "Female" & constituents$alcohol >= 5 & constituents$alcohol <= 15) |
     (sex == "Male" & constituents$alcohol >= 10 & constituents$alcohol <= 25))
  out$amed_alcohol <- ifelse(is.finite(constituents$alcohol) & sex %in% c("Male", "Female"),
                             as.integer(alc_ok), NA_integer_)
  out$amed <- rowSums(out, na.rm = FALSE)
  out$amed[!(domain %in% TRUE)] <- NA_integer_
  if (any(out$amed < 0 | out$amed > 9, na.rm = TRUE))
    stop("aMED score outside 0--9", call. = FALSE)
  attr(out, "weighted_medians") <- medians
  out
}

build_diet_quality_scores <- function(fped, nutrients, analysis) {
  assert_columns(fped, c("SEQN", paste0(FPED_REQUIRED, "_2d")), "two-day FPED")
  assert_columns(nutrients, c("SEQN", "energy_kcal_2d", "sodium_mg_2d", "satfat_g_2d",
                              "mufa_g_2d", "pufa_g_2d", "alcohol_g_2d"), "nutrient totals")
  assert_columns(analysis, c("SEQN", "sex", "analysis_weight", "base_domain"), "analysis frame")
  if (anyDuplicated(fped$SEQN) || anyDuplicated(nutrients$SEQN) || anyDuplicated(analysis$SEQN))
    stop("Diet-quality inputs must be unique by SEQN", call. = FALSE)
  z <- merge(fped, nutrients, by = "SEQN", all = FALSE, sort = FALSE)
  z <- merge(z, analysis[c("SEQN", "sex", "analysis_weight", "base_domain")],
             by = "SEQN", all = FALSE, sort = FALSE)
  hei <- hei2015_two_day(z, z$energy_kcal_2d, z$sodium_mg_2d, z$satfat_g_2d,
                         z$mufa_g_2d, z$pufa_g_2d)
  ac <- amed_constituents_two_day(z, z$satfat_g_2d, z$mufa_g_2d, z$alcohol_g_2d)
  amed <- amed_score(ac, z$sex, z$analysis_weight, z$base_domain)
  out <- cbind(z[c("SEQN")], hei, amed)
  attr(out, "amed_weighted_medians") <- attr(amed, "weighted_medians")
  out
}

build_pooled_diet_quality_scores <- function(fped_by_cycle, nutrients_by_cycle, analysis) {
  if (!length(fped_by_cycle) || !length(nutrients_by_cycle) ||
      length(fped_by_cycle) != length(nutrients_by_cycle)) {
    stop("Pooled diet-quality inputs must contain matched analytic cycles", call. = FALSE)
  }
  fped <- do.call(rbind, fped_by_cycle)
  nutrients <- do.call(rbind, nutrients_by_cycle)
  build_diet_quality_scores(fped, nutrients, analysis)
}

validate_diet_quality_contract <- function() {
  stopifnot(sum(HEI2015_COMPONENT_MAX) == 100, length(AMED_COMPONENTS) == 9L)
  x <- as.data.frame(setNames(replicate(length(FPED_REQUIRED), c(1, 2), simplify = FALSE),
                              paste0(FPED_REQUIRED, "_2d")))
  protein <- fped_two_day_constituents(x)
  h <- hei2015_two_day(x, c(2000, 2000), c(2200, 4000), c(17.7, 35.6),
                       c(30, 10), c(20, 5))
  stopifnot(nrow(h) == 2L, all(h$hei2015 >= 0 & h$hei2015 <= 100),
            identical(protein$total_protein, c(10, 20)),
            identical(protein$seafood_plant, c(5, 10)),
            all.equal(score_moderation(c(1.1, 2), 1.1, 2), c(10, 0)) == TRUE,
            all.equal(score_adequacy(c(0, 0.8), 0.8, 5), c(0, 5)) == TRUE)
  invisible(TRUE)
}

read_two_day_nutrients <- function(cycle) {
  vars <- c(energy_kcal = "TKCAL", sodium_mg = "TSODI", satfat_g = "TSFAT",
            mufa_g = "TMFAT", pufa_g = "TPFAT", alcohol_g = "TALCO")
  read_day <- function(day) {
    path <- diet_file(cycle, day, "tot")
    if (!file.exists(path)) stop("Missing dietary total file: ", path, call. = FALSE)
    prefix <- paste0("DR", day)
    wanted <- c("SEQN", paste0(prefix, unname(vars)))
    x <- as.data.frame(haven::read_xpt(path))
    assert_columns(x, wanted, basename(path))
    x <- x[wanted]
    names(x) <- c("SEQN", names(vars))
    if (anyDuplicated(x$SEQN)) stop("Duplicate SEQN in ", path, call. = FALSE)
    names(x)[-1] <- paste0(names(x)[-1], "_d", day)
    x
  }
  z <- merge(read_day(1), read_day(2), by = "SEQN", all = FALSE, sort = FALSE)
  out <- data.frame(SEQN = z$SEQN)
  for (nm in names(vars)) {
    a <- as.numeric(z[[paste0(nm, "_d1")]])
    b <- as.numeric(z[[paste0(nm, "_d2")]])
    out[[paste0(nm, "_2d")]] <- ifelse(is.finite(a) & is.finite(b) & a >= 0 & b >= 0,
                                          a + b, NA_real_)
  }
  out
}

stage_diet_quality_scores <- function() {
  require_formal_authorization("diet_quality_scores")
  analysis_path <- file.path(DIRS$config, "analysis_frame_pre_mice.rds")
  if (!file.exists(analysis_path)) stop("Run dataset construction first", call. = FALSE)
  analysis <- readRDS(analysis_path)
  assert_columns(analysis, c("SEQN", "cycle", "sex", "analysis_weight", "base_domain"),
                 "analysis frame")
  if (anyDuplicated(analysis$SEQN)) stop("Analysis frame must be unique by SEQN", call. = FALSE)
  existing <- intersect(c("hei2015", "amed"), names(analysis))
  if (length(existing)) stop("Diet-quality fields already exist; refusing an ambiguous rerun",
                             call. = FALSE)
  cycles <- intersect(ANALYTIC_CYCLES, unique(as.character(analysis$cycle)))
  if (!length(cycles)) stop("No supported analytic cycle in analysis frame", call. = FALSE)
  fped_by_cycle <- lapply(cycles, function(cycle) {
    z <- combine_fped_two_day(read_fped_day(cycle, 1), read_fped_day(cycle, 2))
    z$cycle <- cycle
    z
  })
  nutrients_by_cycle <- lapply(cycles, function(cycle) {
    z <- read_two_day_nutrients(cycle)
    z$cycle <- cycle
    z
  })
  persons <- analysis[as.character(analysis$cycle) %in% cycles,
                      c("SEQN", "sex", "analysis_weight", "base_domain"), drop = FALSE]
  scores <- build_pooled_diet_quality_scores(fped_by_cycle, nutrients_by_cycle, persons)
  if (anyDuplicated(scores$SEQN)) stop("Diet-quality result duplicated SEQN", call. = FALSE)
  score_cols <- setdiff(names(scores), "SEQN")
  idx <- match(analysis$SEQN, scores$SEQN)
  analysis[score_cols] <- lapply(scores[score_cols], function(v) v[idx])
  dir.create(DIRS$config, recursive = TRUE, showWarnings = FALSE)
  score_path <- file.path(DIRS$config, "diet_quality_scores.rds")
  saveRDS(scores, score_path, compress = "xz")
  saveRDS(analysis, analysis_path, compress = "xz")
  c(scores = score_path, analysis_frame = analysis_path)
}

if (sys.nframe() == 0L) {
  stage_diet_quality_scores()
}
)--------------------"

upf_source[["GA/01_代码/10_ga_shared.R"]] <- r"--------------------(suppressPackageStartupMessages({
  library(survey)
  library(mitools)
  library(VGAM)
  library(svyVGAM)
})
options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)

GA_LEVELS <- c("Low", "Moderate", "High", "Very high")
G_LEVELS <- c("G1", "G2", "G3a", "G3b", "G4", "G5")
A_LEVELS <- c("A1", "A2", "A3")
GA_MAP <- c(
  G1_A1 = "Low", G2_A1 = "Low",
  G1_A2 = "Moderate", G2_A2 = "Moderate", G3a_A1 = "Moderate",
  G1_A3 = "High", G2_A3 = "High", G3a_A2 = "High", G3b_A1 = "High",
  G3a_A3 = "Very high", G3b_A2 = "Very high", G3b_A3 = "Very high",
  G4_A1 = "Very high", G4_A2 = "Very high", G4_A3 = "Very high",
  G5_A1 = "Very high", G5_A2 = "Very high", G5_A3 = "Very high"
)

derive_ga_fields <- function(d, egfr = "egfr_2021", uacr = "acr_mg_g") {
  ok <- is.finite(d[[egfr]]) & is.finite(d[[uacr]])
  g <- cut(d[[egfr]], breaks = c(-Inf, 15, 30, 45, 60, 90, Inf), right = FALSE,
           labels = c("G5", "G4", "G3b", "G3a", "G2", "G1"))
  a <- cut(d[[uacr]], breaks = c(-Inf, 30, 300, Inf), right = FALSE,
           labels = c("A1", "A2", "A3"))
  risk <- rep(NA_character_, nrow(d))
  risk[ok] <- unname(GA_MAP[paste(g[ok], a[ok], sep = "_")])
  d$g_category <- factor(as.character(g), levels = G_LEVELS)
  d$a_category <- factor(as.character(a), levels = A_LEVELS)
  d$ga_risk <- ordered(risk, levels = GA_LEVELS)
  d$ga_domain <- ok
  d
}

add_fixed_quartiles <- function(d, base_domain, exposure, weight) {
  keep <- d[[base_domain]] %in% TRUE & is.finite(d[[exposure]]) & is.finite(d[[weight]]) & d[[weight]] > 0
  qfun <- function(p) {
    o <- order(d[[exposure]][keep]); x <- d[[exposure]][keep][o]; w <- d[[weight]][keep][o]
    x[which(cumsum(w) / sum(w) >= p)[1]]
  }
  cuts <- vapply(c(.25, .50, .75), qfun, numeric(1))
  if (any(diff(cuts) <= 0)) stop("Invalid weighted UPF quartile cuts")
  d$upf_quartile <- cut(d[[exposure]], c(-Inf, cuts, Inf), labels = paste0("Q", 1:4), include.lowest = TRUE)
  d$upf_quartile <- factor(d$upf_quartile, levels = paste0("Q", 1:4))
  attr(d, "upf_quartile_cuts") <- cuts
  d
}

weighted_indicator_rows <- function(design, database, group = "Overall") {
  rows <- lapply(GA_LEVELS, function(level) {
    design$variables$.indicator <- as.numeric(as.character(design$variables$ga_risk) == level)
    est <- svymean(~.indicator, design, na.rm = TRUE)
    ci <- confint(est)
    data.frame(database = database, group = group, risk = level,
               n_unweighted = sum(as.character(design$variables$ga_risk) == level, na.rm = TRUE),
               weighted_proportion = unname(coef(est)), se = unname(SE(est)),
               ci_low = ci[1], ci_high = ci[2])
  })
  do.call(rbind, rows)
}

weighted_ga_descriptives <- function(design, database) {
  overall <- weighted_indicator_rows(design, database, "Overall")
  by_q <- do.call(rbind, lapply(paste0("Q", 1:4), function(q) {
    sub <- design[design$variables$upf_quartile == q, ]
    weighted_indicator_rows(sub, database, q)
  }))
  list(risk = rbind(overall, by_q))
}

compact_fit <- function(fit, df_complete) {
  b <- coef(fit); v <- vcov(fit)
  if (!identical(names(b), rownames(v)) || any(!is.finite(b)) || any(!is.finite(v))) {
    stop("Non-finite or unnamed coefficient/covariance output")
  }
  list(coefficients = b, variance = v, df_complete = max(1, as.numeric(df_complete)))
}

pool_scalar <- function(summaries, term, sign_multiplier = 1, label = term, exponentiate = TRUE) {
  transformed <- lapply(summaries, function(s) {
    if (!term %in% names(s$coefficients)) stop("Missing term: ", term)
    b <- sign_multiplier * unname(s$coefficients[term])
    v <- unname(s$variance[term, term, drop = FALSE])
    list(coefficients = setNames(b, label),
         variance = matrix(v, 1, 1, dimnames = list(label, label)),
         df_complete = s$df_complete)
  })
  pooled <- mitools::MIcombine(
    lapply(transformed, `[[`, "coefficients"),
    lapply(transformed, `[[`, "variance"),
    df.complete = min(vapply(transformed, `[[`, numeric(1), "df_complete"))
  )
  b <- unname(pooled$coefficients[[1]]); se <- sqrt(pooled$variance[1, 1]); df <- unname(pooled$df[[1]])
  crit <- if (is.finite(df)) qt(.975, df) else qnorm(.975)
  lo <- b - crit * se; hi <- b + crit * se
  p <- if (is.finite(df)) 2 * pt(abs(b / se), df, lower.tail = FALSE) else 2 * pnorm(abs(b / se), lower.tail = FALSE)
  data.frame(estimate_link = b, standard_error_link = se, df = df, p_value = p,
             estimate = if (exponentiate) exp(b) else b,
             ci_low = if (exponentiate) exp(lo) else lo,
             ci_high = if (exponentiate) exp(hi) else hi,
             successful_imputations = length(summaries))
}

d1_joint_from_contrasts <- function(summaries, terms, contrast_matrix) {
  q <- do.call(cbind, lapply(summaries, function(s) as.numeric(contrast_matrix %*% s$coefficients[terms])))
  u <- simplify2array(lapply(summaries, function(s) contrast_matrix %*% s$variance[terms, terms, drop = FALSE] %*% t(contrast_matrix)))
  k <- nrow(q); m <- ncol(q); qbar <- rowMeans(q); ubar <- apply(u, c(1, 2), mean); b <- cov(t(q))
  r <- (1 + 1 / m) * sum(diag(b %*% solve(ubar))) / k
  f <- as.numeric(t(qbar) %*% solve((1 + r) * ubar) %*% qbar / k)
  tval <- k * (m - 1); dfcom <- min(vapply(summaries, `[[`, numeric(1), "df_complete"))
  a <- r * tval / (tval - 2); vstar <- ((dfcom + 1) / (dfcom + 3)) * dfcom
  c0 <- 1 / (tval - 4); c1 <- vstar - 2 * (1 + a); c2 <- vstar - 4 * (1 + a)
  z <- 1 / c2 + c0 * (a^2 * c1 / ((1 + a)^2 * c2)) +
    c0 * (8 * a^2 * c1 / ((1 + a) * c2^2) + 4 * a^2 / ((1 + a) * c2)) +
    c0 * (4 * a^2 / (c2 * c1) + 16 * a^2 * c1 / c2^3) + c0 * (8 * a^2 / c2^2)
  df2 <- 4 + 1 / z
  data.frame(f_statistic = f, df1 = k, df2 = df2, p_value = pf(f, k, df2, lower.tail = FALSE),
             relative_increase_variance = r, complete_data_df = dfcom)
}

fit_parallel_ordinal <- function(d, covariates, design_fun) {
  form <- reformulate(c("upf_per10", covariates), response = "ga_risk")
  des <- design_fun(d)
  des <- des[des$variables$ga_domain %in% TRUE & is.finite(des$variables$upf_per10), ]
  fit <- suppressWarnings(svyVGAM::svy_vglm(form, VGAM::cumulative(parallel = TRUE, reverse = FALSE), design = des, crit = "coef"))
  dfc <- survey::degf(des) - (length(coef(fit)) - 1)
  compact_fit(fit, dfc)
}

fit_threshold_stack <- function(d, covariates, design_fun) {
  d <- d[d$ga_domain %in% TRUE & is.finite(d$upf_per10), , drop = FALSE]
  nrisk <- as.integer(d$ga_risk)
  thresholds <- c(T1 = 2L, T2 = 3L, T3 = 4L)
  pieces <- lapply(names(thresholds), function(nm) {
    z <- d; z$threshold <- nm; z$.worse <- as.integer(nrisk >= thresholds[[nm]]); z
  })
  stacked <- do.call(rbind, pieces); rownames(stacked) <- NULL
  stacked$threshold <- factor(stacked$threshold, levels = names(thresholds))
  rhs <- c("0 + threshold", paste0("threshold:", c("upf_per10", covariates)))
  form <- as.formula(paste(".worse ~", paste(rhs, collapse = " + ")))
  des <- design_fun(stacked)
  fit <- svyglm(form, design = des, family = quasibinomial(), na.action = na.fail)
  out <- compact_fit(fit, df.residual(fit))
  candidates <- names(out$coefficients)[grepl("upf_per10", names(out$coefficients), fixed = TRUE)]
  term_map <- vapply(names(thresholds), function(nm) {
    hit <- candidates[grepl(nm, candidates, fixed = TRUE)]
    if (length(hit) != 1L) stop("Cannot map threshold exposure term for ", nm)
    hit
  }, character(1))
  out$exposure_terms <- term_map
  out
}

fit_multinomial <- function(d, covariates, design_fun) {
  d$ga_risk <- factor(as.character(d$ga_risk), levels = GA_LEVELS)
  form <- reformulate(c("upf_per10", covariates), response = "ga_risk")
  des <- design_fun(d)
  des <- des[des$variables$ga_domain %in% TRUE & is.finite(des$variables$upf_per10), ]
  fit <- suppressWarnings(svyVGAM::svy_vglm(form, VGAM::multinomial(refLevel = "Low"), design = des, crit = "coef"))
  dfc <- survey::degf(des) - (length(coef(fit)) - 1)
  out <- compact_fit(fit, dfc)
  terms <- names(out$coefficients)[grepl("upf_per10", names(out$coefficients), fixed = TRUE)]
  if (length(terms) != 3L) stop("Expected three multinomial UPF terms")
  out$exposure_terms <- setNames(terms, GA_LEVELS[-1])
  out
}

write_ga_cells <- function(design, database, path) {
  rows <- list(); k <- 0L
  for (q in c("Overall", paste0("Q", 1:4))) {
    desq <- if (q == "Overall") design else design[design$variables$upf_quartile == q, ]
    for (g in G_LEVELS) for (a in A_LEVELS) {
      k <- k + 1L
      desq$variables$.cell <- as.numeric(as.character(desq$variables$g_category) == g & as.character(desq$variables$a_category) == a)
      est <- svymean(~.cell, desq, na.rm = TRUE); ci <- confint(est)
      rows[[k]] <- data.frame(database = database, group = q, g_category = g, a_category = a,
                              n_unweighted = sum(desq$variables$.cell, na.rm = TRUE),
                              weighted_proportion = unname(coef(est)), se = unname(SE(est)),
                              ci_low = ci[1], ci_high = ci[2])
    }
  }
  write.csv(do.call(rbind, rows), path, row.names = FALSE)
}
)--------------------"

upf_source[["GA/01_代码/11_ga_nhanes.R"]] <- r"--------------------(root <- "/storage/home/tmu2301/nhanes分析/18_formal_cholesterol_design_reanalysis_20260907/GA"
.libPaths(c("/storage/home/tmu2301/nhanes分析/11_UPF六结局_GA与跨库整合_20260903/R_library", .libPaths()))
source(file.path(root, "01_代码", "10_ga_shared.R"))
nh_code <- "/storage/home/tmu2301/nhanes分析/18_formal_cholesterol_design_reanalysis_20260907/NHANES/02_代码/正式五结局"
source(file.path(nh_code, "00_functions.R"))
source(file.path(nh_code, "04_models.R"))

out_dir <- file.path(root, "02_结果", "G_A", "NHANES")
qc_dir <- file.path(root, "03_质量控制", "G_A", "NHANES")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

frame_path <- file.path(DIRS$derived, "analysis_frame_pre_mice.rds")
mids_path <- file.path(DIRS$derived, "mi_CKD.rds")
d <- readRDS(frame_path); mi <- readRDS(mids_path)
stopifnot(mi$m == 50L)
d <- derive_ga_fields(d)
d$ga_domain <- d$base_domain %in% TRUE & d$ga_domain & is.finite(d$upf_gram_ratio_nonwater)
d$upf_per10 <- d$upf_gram_ratio_nonwater / 10
d <- add_fixed_quartiles(d, "base_domain", "upf_gram_ratio_nonwater", "analysis_weight")

design_fun <- function(x) make_full_design(x)
full <- design_fun(d); ga_design <- full[full$variables$ga_domain %in% TRUE, ]
desc <- weighted_ga_descriptives(ga_design, "NHANES")
write.csv(desc$risk, file.path(out_dir, "ga_risk_prevalence_overall_and_quartile.csv"), row.names = FALSE)
write_ga_cells(ga_design, "NHANES", file.path(out_dir, "ga_18cell_distribution_overall_and_quartile.csv"))
write.csv(data.frame(statistic = c("Q1_Q2", "Q2_Q3", "Q3_Q4"),
                     upf_percent = attr(d, "upf_quartile_cuts")),
          file.path(qc_dir, "upf_weighted_quartile_cuts.csv"), row.names = FALSE)
cat("NH_GA_DESCRIPTIVES_DONE n=", sum(d$ga_domain), "\n", sep = "")

models <- names(MODEL_COVARIATES)
completed_for <- function(i) {
  z <- restore_completed_subframe(d, mi, i)
  z$ga_risk <- ordered(as.character(z$ga_risk), levels = GA_LEVELS)
  z
}

ordinal_by_model <- parallel::mclapply(models, function(model) {
  covars <- phase1_model_covariates(model, "linear", "CKD", FALSE)
  if (model == "M1") {
    one <- fit_parallel_ordinal(completed_for(1), covars, design_fun)
    out <- replicate(50, one, simplify = FALSE)
  } else {
    out <- lapply(1:50, function(i) fit_parallel_ordinal(completed_for(i), covars, design_fun))
  }
  cat("NH_GA_ORDINAL_MODEL_DONE ", model, "\n", sep = "")
  out
}, mc.cores = n_cores_resolve(), mc.preschedule = TRUE)
names(ordinal_by_model) <- models

ordinal_rows <- do.call(rbind, lapply(models, function(model) {
  ans <- pool_scalar(ordinal_by_model[[model]], "upf_per10", sign_multiplier = -1,
                     label = "log_OR_worse_GA_risk", exponentiate = TRUE)
  cbind(database = "NHANES", model = model, analysis = "proportional_odds",
        contrast = "Higher versus lower KDIGO G-A risk", ans)
}))
write.csv(ordinal_rows, file.path(out_dir, "ga_proportional_odds_M1_M4.csv"), row.names = FALSE)

m3_covars <- phase1_model_covariates("M3", "linear", "CKD", FALSE)
threshold_summaries <- parallel::mclapply(1:50, function(i) {
  fit_threshold_stack(completed_for(i), m3_covars, design_fun)
}, mc.cores = n_cores_resolve(), mc.preschedule = TRUE)
terms <- unname(threshold_summaries[[1]]$exposure_terms[c("T1", "T2", "T3")])
if (!all(vapply(threshold_summaries, function(x) identical(unname(x$exposure_terms[c("T1", "T2", "T3")]), terms), logical(1)))) stop("Threshold term mapping changed across imputations")
L <- rbind(c(-1, 1, 0), c(-1, 0, 1))
po_test <- d1_joint_from_contrasts(threshold_summaries, terms, L)
po_test$database <- "NHANES"; po_test$model <- "M3"; po_test$test <- "Equality of UPF slopes across three cumulative thresholds"
write.csv(po_test, file.path(out_dir, "ga_proportional_odds_assumption_M3.csv"), row.names = FALSE)

threshold_labels <- c(T1 = "Moderate/high/very high vs low", T2 = "High/very high vs low/moderate", T3 = "Very high vs lower risk")
threshold_rows <- do.call(rbind, lapply(names(threshold_labels), function(nm) {
  ans <- pool_scalar(threshold_summaries, threshold_summaries[[1]]$exposure_terms[[nm]], exponentiate = TRUE)
  cbind(database = "NHANES", model = "M3", threshold = nm, contrast = threshold_labels[[nm]], ans)
}))
write.csv(threshold_rows, file.path(out_dir, "ga_cumulative_threshold_effects_M3.csv"), row.names = FALSE)
cat("NH_GA_PO_TEST_DONE p=", format(po_test$p_value, digits = 8), "\n", sep = "")

use_multinomial <- is.finite(po_test$p_value) && po_test$p_value < 0.05
if (use_multinomial) {
  multi_by_model <- parallel::mclapply(models, function(model) {
    covars <- phase1_model_covariates(model, "linear", "CKD", FALSE)
    if (model == "M1") {
      one <- fit_multinomial(completed_for(1), covars, design_fun)
      out <- replicate(50, one, simplify = FALSE)
    } else {
      out <- lapply(1:50, function(i) fit_multinomial(completed_for(i), covars, design_fun))
    }
    cat("NH_GA_MULTINOMIAL_MODEL_DONE ", model, "\n", sep = "")
    out
  }, mc.cores = n_cores_resolve(), mc.preschedule = TRUE)
  names(multi_by_model) <- models
  multi_rows <- do.call(rbind, lapply(models, function(model) {
    do.call(rbind, lapply(GA_LEVELS[-1], function(level) {
      term <- multi_by_model[[model]][[1]]$exposure_terms[[level]]
      ans <- pool_scalar(multi_by_model[[model]], term, exponentiate = TRUE)
      cbind(database = "NHANES", model = model, analysis = "multinomial",
            contrast = paste(level, "vs Low"), ans)
    }))
  }))
  write.csv(multi_rows, file.path(out_dir, "ga_multinomial_M1_M4.csv"), row.names = FALSE)
  selected <- multi_rows[multi_rows$model == "M3", ]
  selected$selection_reason <- "Proportional-odds slope equality test P<0.05"
} else {
  selected <- ordinal_rows[ordinal_rows$model == "M3", ]
  selected$selection_reason <- "Proportional-odds slope equality test P>=0.05"
}
write.csv(selected, file.path(out_dir, "ga_primary_selected_M3.csv"), row.names = FALSE)

qc <- data.frame(
  check = c("m", "ga_unweighted_n", "risk_levels_observed", "ordinal_rows", "po_test_finite", "selected_model", "all_selected_finite"),
  value = c(mi$m, sum(d$ga_domain), length(unique(na.omit(as.character(d$ga_risk)))), nrow(ordinal_rows),
            is.finite(po_test$p_value), if (use_multinomial) "multinomial" else "proportional_odds",
            all(is.finite(selected$estimate) & is.finite(selected$p_value)))
)
write.csv(qc, file.path(qc_dir, "ga_integrity_checks.csv"), row.names = FALSE)
write.csv(data.frame(role = c("analysis_frame", "CKD_mids"), path = c(frame_path, mids_path),
                     md5 = unname(tools::md5sum(c(frame_path, mids_path)))),
          file.path(qc_dir, "input_lineage.csv"), row.names = FALSE)
cat("NH_GA_COMPLETE selected=", if (use_multinomial) "multinomial" else "proportional_odds", "\n", sep = "")
)--------------------"

upf_source[["01_code/modules/RCS/12_recalculate_rcs_v3.R"]] <- r"--------------------(args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !toupper(args[[1]]) %in% c("NHANES", "KNHANES")) {
  stop("Usage: Rscript 12_recalculate_rcs_v3.R NHANES|KNHANES", call. = FALSE)
}
cohort <- toupper(args[[1]])
figure_root <- Sys.getenv(
  "FIGURE_PROJECT_ROOT",
  "/storage/home/tmu2301/nhanes分析/18_formal_cholesterol_design_reanalysis_20260907/RCS"
)
source_dir <- file.path(figure_root, "02_绘图数据_v3重算", "统计重算源")
log_dir <- file.path(figure_root, "08_运行日志")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(log_dir, paste0("12_recalculate_rcs_v3_", tolower(cohort), ".log"))
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)
cat("cohort=", cohort, "\nstarted_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")

write_csv <- function(x, filename) {
  path <- file.path(source_dir, filename)
  readr::write_excel_csv(x, path, na = "")
  path
}

write_manifest <- function(paths, prefix) {
  paths <- normalizePath(paths, mustWork = TRUE)
  info <- file.info(paths)
  manifest <- data.frame(
    cohort = cohort,
    path = paths,
    bytes = info$size,
    modified = format(info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
    sha256 = vapply(paths, digest::digest, character(1), file = TRUE, algo = "sha256", serialize = FALSE),
    stringsAsFactors = FALSE
  )
  write_csv(manifest, paste0(prefix, "_manifest.csv"))
}

rcs3_basis_scaled <- function(x_pp, knots_pp, unit_pp = 10) {
  x <- as.numeric(x_pp) / unit_pp
  knots <- as.numeric(knots_pp) / unit_pp
  if (length(knots) != 3L || any(!is.finite(knots)) || any(diff(knots) <= 0)) {
    stop("RCS requires three finite ordered knots", call. = FALSE)
  }
  k1 <- knots[[1]]; k2 <- knots[[2]]; k3 <- knots[[3]]
  nonlinear <- (
    pmax(x - k1, 0)^3 -
      pmax(x - k2, 0)^3 * (k3 - k1) / (k3 - k2) +
      pmax(x - k3, 0)^3 * (k2 - k1) / (k3 - k2)
  ) / (k3 - k1)^2
  cbind(rcs_linear = x, rcs_nonlin1 = nonlinear)
}

if (cohort == "NHANES") {
  nh_root <- Sys.getenv(
    "NHANES_ROOT",
    "/storage/home/tmu2301/nhanes分析/18_formal_cholesterol_design_reanalysis_20260907/NHANES"
  )
  nh_code <- file.path(nh_root, "02_代码", "正式五结局")
  source(file.path(nh_code, "23_extension_common.R"), local = FALSE)
  require_formal_authorization("v3_rcs_recalculation")

  frame_path <- file.path(DIRS$derived, "analysis_frame_pre_mice.rds")
  d <- readRDS(frame_path)
  full0 <- make_full_design(d)
  base <- full0$variables$base_domain %in% TRUE &
    is.finite(full0$variables[[EXPOSURE_VAR]])
  design0 <- full0[base, ]
  probs <- c(0.01, 0.10, 0.50, 0.90, 0.99)
  quantiles <- weighted_quantiles(design0, EXPOSURE_VAR, probs)
  names(quantiles) <- c("p01", "p10", "p50", "p90", "p99")
  knots_pp <- unname(quantiles[c("p10", "p50", "p90")])
  reference_pp <- 10
  display_pp <- unname(quantiles[c("p01", "p99")])
  grid_pp_v3 <- sort(unique(c(
    seq(display_pp[[1]], display_pp[[2]], length.out = 121L),
    reference_pp, knots_pp
  )))
  basis <- rcs3_basis_scaled(d[[EXPOSURE_VAR]], knots_pp, EXPOSURE_UNIT_PP)
  d$rcs_linear_v3 <- basis[, "rcs_linear"]
  d$rcs_nonlin1_v3 <- basis[, "rcs_nonlin1"]

  rcs_formula_v3 <- function(outcome) {
    stats::reformulate(
      c(
        "rcs_linear_v3", "rcs_nonlin1_v3",
        phase1_model_covariates("M3", "linear", outcome, FALSE)
      ),
      response = OUTCOMES[[outcome]]$variable
    )
  }

  fit_one <- function(imputation_number, outcome, mi_object) {
    completed <- restore_completed_subframe(d, mi_object, imputation_number)
    full <- make_full_design(completed)
    domain <- full$variables[[OUTCOMES[[outcome]]$domain]] %in% TRUE &
      is.finite(full$variables[[EXPOSURE_VAR]])
    form <- rcs_formula_v3(outcome)
    model_vars <- all.vars(form)
    incomplete <- model_vars[vapply(full$variables[model_vars], function(x) {
      any(is.na(x[domain]))
    }, logical(1))]
    if (length(incomplete)) {
      stop("Post-MICE missing V3 RCS fields: ", paste(incomplete, collapse = ", "), call. = FALSE)
    }
    design <- full[domain, ]
    fit <- survey::svyglm(
      form, design = design, family = outcome_family(outcome), na.action = na.fail
    )
    w <- as.numeric(stats::weights(design, "sampling"))
    y <- design$variables[[OUTCOMES[[outcome]]$variable]]
    list(
      coefficients = stats::coef(fit),
      variance = stats::vcov(fit),
      df_complete = stats::df.residual(fit),
      audit = c(
        n_unweighted = nrow(design$variables),
        events = if (OUTCOMES[[outcome]]$family == "quasibinomial") sum(y == 1L) else NA_real_,
        weighted_population = sum(w),
        effective_sample_size = sum(w)^2 / sum(w^2),
        survey_design_df = survey::degf(design),
        complete_data_residual_df = stats::df.residual(fit)
      )
    )
  }

  curve_from_pool_v3 <- function(pooled, outcome) {
    grid <- grid_pp_v3
    basis_grid <- rcs3_basis_scaled(grid, knots_pp, EXPOSURE_UNIT_PP)
    basis_ref <- as.numeric(rcs3_basis_scaled(reference_pp, knots_pp, EXPOSURE_UNIT_PP)[1, ])
    terms <- c("rcs_linear_v3", "rcs_nonlin1_v3")
    beta <- pooled$coefficients[terms]
    variance <- pooled$variance[terms, terms, drop = FALSE]
    df <- min(pooled$df[match(terms, names(pooled$coefficients))], na.rm = TRUE)
    effect_scale <- effect_scale_for_outcome(outcome)
    do.call(rbind, lapply(seq_along(grid), function(i) {
      contrast <- as.numeric(basis_grid[i, ] - basis_ref)
      link <- sum(contrast * beta)
      se <- sqrt(as.numeric(t(contrast) %*% variance %*% contrast))
      transformed <- transform_link_contrast(link, se, df, effect_scale)
      data.frame(
        cohort = "NHANES", version = "v3_rcs_3knot_p10ref_20260828",
        outcome = outcome, model = "M3", exposure_value = grid[[i]],
        reference_value = reference_pp, knot_10 = knots_pp[[1]],
        knot_50 = knots_pp[[2]], knot_90 = knots_pp[[3]],
        display_p01 = display_pp[[1]], display_p99 = display_pp[[2]],
        effect_measure = unname(c(
          odds_ratio = "OR_vs_reference",
          identity = "beta_difference_vs_reference",
          log_response_percent = "percent_difference_vs_reference"
        )[effect_scale]),
        estimate = transformed[[1]], conf_low = transformed[[2]], conf_high = transformed[[3]],
        estimate_link = link, standard_error_link = se,
        stringsAsFactors = FALSE
      )
    }))
  }

  d1_joint_test_v3 <- function(summaries, terms) {
    aliases <- paste0("b", seq_along(terms))
    qhat <- lapply(summaries, function(z) {
      stats::setNames(as.numeric(z$coefficients[terms]), aliases)
    })
    uhat <- lapply(summaries, function(z) {
      v <- z$variance[terms, terms, drop = FALSE]
      dimnames(v) <- list(aliases, aliases)
      v
    })
    df_complete <- min(vapply(summaries, `[[`, numeric(1), "df_complete"))
    test <- mitml::testConstraints(
      qhat = qhat, uhat = uhat, constraints = aliases,
      method = "D1", df.com = df_complete
    )$test[1, ]
    values <- unname(test[c("F.value", "df1", "df2", "P(>F)")])
    if (any(!is.finite(values)) || test[["df1"]] <= 0 || test[["df2"]] <= 0) {
      stop("NHANES V3 D1 test returned invalid results", call. = FALSE)
    }
    list(
      f_value = unname(test[["F.value"]]), df1 = unname(test[["df1"]]),
      df2 = unname(test[["df2"]]), p_value = unname(test[["P(>F)"]])
    )
  }

  test_rows <- list(); curve_rows <- list(); pooled_models <- list()
  workers <- n_cores_resolve()
  for (outcome in names(OUTCOMES)) {
    cat("outcome_start=", outcome, "\n", sep = "")
    mi_object <- readRDS(file.path(DIRS$derived, paste0("mi_", outcome, ".rds")))
    stopifnot(mi_object$m == MI_M)
    summaries <- parallel::mclapply(
      seq_len(MI_M), fit_one, outcome = outcome, mi_object = mi_object,
      mc.cores = workers, mc.preschedule = TRUE
    )
    if (any(vapply(summaries, inherits, logical(1), "try-error"))) {
      stop("NHANES V3 RCS failed for ", outcome, call. = FALSE)
    }
    pooled <- pool_svy_summaries(summaries)
    overall <- d1_joint_test_v3(summaries, c("rcs_linear_v3", "rcs_nonlin1_v3"))
    nonlinear <- d1_joint_test_v3(summaries, "rcs_nonlin1_v3")
    audit <- summaries[[1]]$audit
    test_rows[[outcome]] <- data.frame(
      cohort = "NHANES", version = "v3_rcs_3knot_p10ref_20260828",
      outcome = outcome, model = "M3", knot_probabilities = "0.10;0.50;0.90",
      knots = paste(format(knots_pp, digits = 12), collapse = ";"),
      reference_probability = NA_real_, reference = reference_pp,
      display_probabilities = "0.01;0.99", display_min = display_pp[[1]],
      display_max = display_pp[[2]], p_overall = overall$p_value,
      p_nonlinear = nonlinear$p_value, overall_f = overall$f_value,
      overall_df1 = overall$df1, overall_df2 = overall$df2,
      nonlinear_f = nonlinear$f_value, nonlinear_df1 = nonlinear$df1,
      nonlinear_df2 = nonlinear$df2, n_unweighted = audit[["n_unweighted"]],
      events = audit[["events"]], effective_sample_size = audit[["effective_sample_size"]],
      successful_imputations = length(summaries),
      formula = paste(deparse(rcs_formula_v3(outcome)), collapse = " "),
      covariates = paste(phase1_model_covariates("M3", "linear", outcome, FALSE), collapse = ";"),
      stringsAsFactors = FALSE
    )
    curve_rows[[outcome]] <- curve_from_pool_v3(pooled, outcome)
    pooled_models[[outcome]] <- list(
      coefficients = pooled$coefficients, variance = pooled$variance,
      df = pooled$df, missinfo = pooled$missinfo
    )
    rm(mi_object, summaries, pooled)
    invisible(gc(full = TRUE))
    cat("outcome_done=", outcome, "\n", sep = "")
  }
  tests <- do.call(rbind, test_rows); rownames(tests) <- NULL
  curves <- do.call(rbind, curve_rows); rownames(curves) <- NULL
  stopifnot(
    nrow(tests) == 5L, nrow(curves) == 5L * length(grid_pp_v3),
    all(tests$successful_imputations == 50L),
    all(is.finite(tests$p_overall)), all(is.finite(tests$p_nonlinear)),
    all(is.finite(curves$estimate)), all(curves$conf_low <= curves$estimate),
    all(curves$estimate <= curves$conf_high)
  )
  curve_path <- write_csv(curves, "NHANES_RCS_M3_v3_curves.csv")
  test_path <- write_csv(tests, "NHANES_RCS_M3_v3_tests.csv")
  pooled_path <- file.path(source_dir, "NHANES_RCS_M3_v3_pooled_models.rds")
  saveRDS(pooled_models, pooled_path, compress = "xz")
  lineage <- data.frame(
    cohort = "NHANES", role = c("analysis_frame", paste0("mi_", names(OUTCOMES))),
    path = c(frame_path, file.path(DIRS$derived, paste0("mi_", names(OUTCOMES), ".rds"))),
    stringsAsFactors = FALSE
  )
  lineage$sha256 <- vapply(lineage$path, digest::digest, character(1), file = TRUE, algo = "sha256", serialize = FALSE)
  lineage_path <- write_csv(lineage, "NHANES_RCS_M3_v3_input_lineage.csv")
  write_manifest(c(curve_path, test_path, pooled_path, lineage_path), "NHANES_RCS_M3_v3")
}

if (cohort == "KNHANES") {
  kn_root <- Sys.getenv(
    "KNHANES_ROOT",
    "/project/05_KNHANES_外部验证_2019_2024_20260820"
  )
  kn_code <- file.path(kn_root, "04_代码")
  source(file.path(kn_code, "07_extended_models.R"), local = FALSE)
  require_user_authorization("v3_rcs_recalculation")
  require_nova_gate("v3_rcs_recalculation")
  inputs <- load_model_inputs()
  frame_q <- derive_fixed_exposure_categories(inputs$frame)
  base <- frame_q$base_domain %in% TRUE
  knots_pp <- weighted_quantile(
    frame_q$upf_grams_pct[base], frame_q$analysis_weight[base], c(0.10, 0.50, 0.90)
  )
  display_pp <- weighted_quantile(
    frame_q$upf_grams_pct[base], frame_q$analysis_weight[base], c(0.01, 0.99)
  )
  reference_pp <- 10
  grid_pp_v3 <- sort(unique(c(
    seq(display_pp[[1]], display_pp[[2]], length.out = 121L),
    reference_pp, knots_pp
  )))

  fit_kn_v3 <- function(frame, imp, outcome_name) {
    spec <- OUTCOMES[[outcome_name]]
    covars <- if (outcome_name == "DKD") MODEL_COVARIATES_DKD$M3 else MODEL_COVARIATES$M3
    knots <- knots_pp / UPF_UNIT_PP
    ref_x <- reference_pp / UPF_UNIT_PP
    frame <- derive_fixed_exposure_categories(frame)
    form <- build_model_formula(spec$variable, c("upf_per10", "upf_rcs_nonlinear"), covars)
    fits <- lapply(seq_len(imp$m), function(i) {
      d <- completed_analysis_dataset(frame, imp, i)
      d$upf_rcs_nonlinear <- rcs_nonlinear_basis(d$upf_per10, knots)
      tryCatch(fit_svy_formula(d, form, spec), error = function(e) NULL)
    })
    if (any(vapply(fits, is.null, logical(1)))) {
      stop("KNHANES V3 RCS failed in at least one imputation for ", outcome_name)
    }
    invisible(lapply(fits, validate_compact_fit, family = spec$family,
                     label = paste0(outcome_name, " V3 RCS")))
    term_names <- names(fits[[1]]$coefficients)
    grid_pp <- grid_pp_v3
    curve <- dplyr::bind_rows(lapply(seq_along(grid_pp), function(j) {
      grid_x <- grid_pp[[j]] / UPF_UNIT_PP
      L <- stats::setNames(rep(0, length(term_names)), term_names)
      L[["upf_per10"]] <- grid_x - ref_x
      L[["upf_rcs_nonlinear"]] <-
        rcs_nonlinear_basis(grid_x, knots) - rcs_nonlinear_basis(ref_x, knots)
      z <- pool_scalar_contrast(fits, L)
      is_reference <- abs(grid_pp[[j]] - reference_pp) <=
        .Machine$double.eps * max(1, abs(reference_pp))
      if (is_reference) {
        z$estimate_link <- 0
        z$std_error_link <- 0
        z$conf_low_link <- 0
        z$conf_high_link <- 0
        z$p_value <- 1
      }
      if (spec$family == "binomial") {
        z$estimate <- exp(z$estimate_link)
        z$conf_low <- exp(z$conf_low_link)
        z$conf_high <- exp(z$conf_high_link)
        z$effect_scale <- "OR"
      } else {
        z$estimate <- z$estimate_link
        z$conf_low <- z$conf_low_link
        z$conf_high <- z$conf_high_link
        z$effect_scale <- "difference"
      }
      z$cohort <- "KNHANES"
      z$version <- "v3_rcs_3knot_p10ref_20260828"
      z$outcome <- outcome_name
      z$model <- "M3"
      z$upf_grams_pct <- grid_pp[[j]]
      z$reference_pp <- reference_pp
      z$knot_10 <- knots_pp[[1]]
      z$knot_50 <- knots_pp[[2]]
      z$knot_90 <- knots_pp[[3]]
      z$display_p01 <- display_pp[[1]]
      z$display_p99 <- display_pp[[2]]
      z$status <- "OK"
      z
    }))
    nonlinear <- mi_d1_test(fits, "upf_rcs_nonlinear")
    nonlinear$test <- "nonlinearity"
    overall <- mi_d1_test(fits, c("upf_per10", "upf_rcs_nonlinear"))
    overall$test <- "overall_spline"
    tests <- dplyr::bind_rows(nonlinear, overall)
    tests$cohort <- "KNHANES"
    tests$version <- "v3_rcs_3knot_p10ref_20260828"
    tests$outcome <- outcome_name
    tests$model <- "M3"
    tests$knot_probabilities <- "0.10;0.50;0.90"
    tests$knots <- paste(format(knots_pp, digits = 12), collapse = ";")
    tests$reference_probability <- NA_real_
    tests$reference <- reference_pp
    tests$display_probabilities <- "0.01;0.99"
    tests$display_min <- display_pp[[1]]
    tests$display_max <- display_pp[[2]]
    tests$n_unweighted <- fits[[1]]$n
    tests$events <- fits[[1]]$events
    tests$successful_imputations <- length(fits)
    tests$formula <- paste(deparse(form), collapse = " ")
    tests$covariates <- paste(covars, collapse = ";")
    tests$status <- "OK"
    list(curve = curve, tests = tests)
  }

  outcomes <- names(OUTCOMES)
  workers <- min(4L, N_WORKERS)
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = workers)
  by_outcome <- future.apply::future_lapply(
    outcomes,
    function(outcome_name) {
      configure_worker_runtime()
      local_inputs <- load_model_inputs()
      fit_kn_v3(local_inputs$frame, local_inputs$imp, outcome_name)
    },
    future.seed = KN_SEED,
    future.packages = c("survey", "mitools", "mitml", "mice", "dplyr")
  )
  curves <- dplyr::bind_rows(lapply(by_outcome, `[[`, "curve"))
  tests <- dplyr::bind_rows(lapply(by_outcome, `[[`, "tests"))
  stopifnot(
    nrow(curves) == length(outcomes) * length(grid_pp_v3), nrow(tests) == 8L,
    all(tests$successful_imputations == 50L), all(tests$status == "OK"),
    all(is.finite(tests$p_value)), all(is.finite(curves$estimate)),
    all(is.finite(curves$conf_low)), all(is.finite(curves$conf_high)),
    all(curves$conf_low <= curves$estimate), all(curves$estimate <= curves$conf_high)
  )
  curve_path <- write_csv(curves, "KNHANES_RCS_M3_v3_curves.csv")
  test_path <- write_csv(tests, "KNHANES_RCS_M3_v3_tests.csv")
  lineage <- data.frame(
    cohort = "KNHANES",
    role = c("analysis_frame", "mids", "code_review_marker", "nova_gate_marker", "active_nova_map"),
    path = c(inputs$frame_path, inputs$mice_path, CODE_REVIEW_MARKER, NOVA_GATE_MARKER, ACTIVE_NOVA_MAP),
    stringsAsFactors = FALSE
  )
  lineage$sha256 <- vapply(lineage$path, digest::digest, character(1), file = TRUE, algo = "sha256", serialize = FALSE)
  lineage_path <- write_csv(lineage, "KNHANES_RCS_M3_v3_input_lineage.csv")
  write_manifest(c(curve_path, test_path, lineage_path), "KNHANES_RCS_M3_v3")
}

cat("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")
cat("status=OK\n")
)--------------------"

upf_source[["01_code/modules/SUA/04_nhanes_serum_uric_acid_extended_core.R"]] <- r"--------------------(options(stringsAsFactors = FALSE)

project_root <- "/storage/home/tmu2301/nhanes分析/18_formal_cholesterol_design_reanalysis_20260907/NHANES"
code_dir <- file.path(project_root, "02_代码", "正式五结局")
source(file.path(code_dir, "23_extension_common.R"), local = FALSE)

MI_M <- 50L
PHASE1_VERSION <- "2026-09-02.serum-uric-acid-extension-v1"
EXTENSION_VERSION <- PHASE1_VERSION
OUTCOMES <- list(
  serum_uric_acid = list(
    variable = "uric_acid_dxc", family = "gaussian",
    domain = "uric_acid_observed", effect_scale = "identity"
  ),
  hyperuricemia = list(
    variable = "hyperuricemia", family = "quasibinomial",
    domain = "uric_acid_observed", effect_scale = "odds_ratio"
  )
)
EXPOSURES <- list(
  G_R0 = list(variable = "upf_gram_ratio_nonwater",
              label = "Two-day non-water gram share")
)

extension_root <- Sys.getenv(
  "SERUM_UA_EXTENSION_ROOT",
  unset = "/storage/home/tmu2301/nhanes分析/18_formal_cholesterol_design_reanalysis_20260907/SUA"
)
result_dir <- file.path(extension_root, "03_扩展分析", "NHANES")
qc_dir <- file.path(extension_root, "05_质量控制", "NHANES_extended")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

frame_path <- file.path(project_root, "03_派生数据", "analysis_frame_pre_mice.rds")
mi_path <- file.path(project_root, "03_派生数据", "mi_uric_acid.rds")
main_path <- file.path(extension_root, "02_主分析", "NHANES",
                       "NHANES_serum_uric_acid_M1_M4.csv")
stopifnot(file.exists(frame_path), file.exists(mi_path), file.exists(main_path))
d <- readRDS(frame_path)
mi_object <- readRDS(mi_path)
stopifnot(mi_object$m == MI_M)

spec <- build_extension_spec(d)
d <- add_extension_exposure_fields(d, spec)
write.csv(extension_spec_table(spec), file.path(qc_dir, "extension_analysis_spec.csv"),
          row.names = FALSE)

compact_fit <- function(data, formula, outcome_name, domain_extra = rep(TRUE, nrow(data))) {
  full <- make_full_design(data)
  model_vars <- all.vars(formula)
  domain <- full$variables[[OUTCOMES[[outcome_name]]$domain]] %in% TRUE &
    domain_extra[match(rownames(full$variables), rownames(data))]
  if (anyNA(domain)) stop("Domain alignment failed")
  incomplete <- model_vars[vapply(full$variables[model_vars], function(x) {
    any(is.na(x[domain]))
  }, logical(1))]
  if (length(incomplete)) stop("Missing model fields: ", paste(incomplete, collapse = ","))
  design <- full[domain, ]
  family <- if (OUTCOMES[[outcome_name]]$family == "quasibinomial") {
    stats::quasibinomial()
  } else stats::gaussian()
  fit <- survey::svyglm(formula, design = design, family = family, na.action = na.fail)
  y <- design$variables[[OUTCOMES[[outcome_name]]$variable]]
  list(
    coefficients = stats::coef(fit), variance = stats::vcov(fit),
    df_complete = stats::df.residual(fit),
    n = nrow(design$variables),
    events = if (OUTCOMES[[outcome_name]]$family == "quasibinomial") sum(y == 1L) else NA_real_,
    design_df = survey::degf(design)
  )
}

pool_summary_list <- function(summaries) {
  if (length(summaries) == 1L) {
    s <- summaries[[1]]
    return(list(
      coefficients = s$coefficients, variance = s$variance,
      df = rep(s$df_complete, length(s$coefficients)),
      missinfo = rep(0, length(s$coefficients))
    ))
  }
  pool_svy_summaries(summaries)
}

term_row <- function(summaries, term, outcome, model, module) {
  pooled <- pool_summary_list(summaries)
  z <- pooled_term_table(pooled, term, OUTCOMES[[outcome]]$effect_scale)
  data.frame(
    database = "NHANES", outcome = outcome, model = model, module = module,
    term = term, estimate = z[["estimate"]], ci_low = z[["ci_low"]],
    ci_high = z[["ci_high"]], p_value = z[["p_value"]],
    standard_error_link = z[["standard_error_link"]],
    estimate_link = z[["estimate_link"]], n_unweighted = summaries[[1]]$n,
    events = summaries[[1]]$events, successful_imputations = length(summaries),
    stringsAsFactors = FALSE
  )
}

model_formula <- function(outcome, model, exposure_terms) {
  stats::reformulate(
    c(exposure_terms, phase1_model_covariates(model, "linear", outcome, FALSE)),
    response = OUTCOMES[[outcome]]$variable
  )
}

fit_imputation_bundle <- function(i) {
  x <- restore_completed_subframe(d, mi_object, i)
  x <- add_extension_exposure_fields(x, spec)
  x_rcs <- add_rcs_fields(x, spec)
  out <- list(quartile = list(), trend = list(), sensitivity = list(),
              hyperuricemia = list())
  for (model in names(MODEL_COVARIATES)) {
    fq <- model_formula("serum_uric_acid", model, "upf_quartile")
    ft <- model_formula("serum_uric_acid", model, "upf_quartile_median10")
    out$quartile[[model]] <- compact_fit(x, fq, "serum_uric_acid")
    out$trend[[model]] <- compact_fit(x, ft, "serum_uric_acid")
    fh <- model_formula("hyperuricemia", model,
                        sprintf("I(%s/10)", EXPOSURE_VAR))
    out$hyperuricemia[[model]] <- compact_fit(x, fh, "hyperuricemia")
  }
  fr <- model_formula("serum_uric_acid", "M3",
                      c("rcs_linear", "rcs_nonlin1", "rcs_nonlin2"))
  out$rcs <- compact_fit(x_rcs, fr, "serum_uric_acid")
  fm3 <- model_formula("serum_uric_acid", "M3",
                       sprintf("I(%s/10)", EXPOSURE_VAR))
  out$sensitivity$exclude_extreme_energy <- compact_fit(
    x, fm3, "serum_uric_acid", x$plausible_energy_both_days %in% TRUE
  )
  out$sensitivity$trim_G_R0_p01_p99 <- compact_fit(
    x, fm3, "serum_uric_acid", x$g_r0_trimmed_domain %in% TRUE
  )
  out
}

workers <- n_cores_resolve()
resource_preflight(workers, "nhanes_serum_uric_acid_extended_core")
cat("NHANES_EXTENDED_START workers=", workers, " m=", MI_M, "\n", sep = "")
bundles <- parallel::mclapply(seq_len(MI_M), fit_imputation_bundle,
                              mc.cores = workers, mc.preschedule = TRUE)
if (any(vapply(bundles, inherits, logical(1), "try-error"))) {
  print(bundles[vapply(bundles, inherits, logical(1), "try-error")])
  stop("NHANES extended analysis failed")
}

quartile_rows <- list()
trend_rows <- list()
hyper_rows <- list()
for (model in names(MODEL_COVARIATES)) {
  qsum <- lapply(bundles, function(z) z$quartile[[model]])
  for (term in paste0("upf_quartileQ", 2:4)) {
    quartile_rows[[paste(model, term)]] <- term_row(
      qsum, term, "serum_uric_acid", model, "weighted_quartile"
    )
  }
  tsum <- lapply(bundles, function(z) z$trend[[model]])
  trend_rows[[model]] <- term_row(
    tsum, "upf_quartile_median10", "serum_uric_acid", model, "quartile_trend"
  )
  hsum <- lapply(bundles, function(z) z$hyperuricemia[[model]])
  hyper_rows[[model]] <- term_row(
    hsum, sprintf("I(%s/10)", EXPOSURE_VAR), "hyperuricemia", model,
    "threshold_supplement"
  )
}
quartile <- do.call(rbind, quartile_rows)
trend <- do.call(rbind, trend_rows)
hyperuricemia <- do.call(rbind, hyper_rows)

rcs_summaries <- lapply(bundles, `[[`, "rcs")
rcs_pooled <- pool_svy_summaries(rcs_summaries)
rcs_nonlinearity <- d1_joint_test(rcs_summaries, c("rcs_nonlin1", "rcs_nonlin2"))
rcs_overall <- d1_joint_test(
  rcs_summaries, c("rcs_linear", "rcs_nonlin1", "rcs_nonlin2")
)
rcs_tests <- rbind(
  transform(rcs_nonlinearity, database = "NHANES", outcome = "serum_uric_acid",
            model = "M3", test = "nonlinearity"),
  transform(rcs_overall, database = "NHANES", outcome = "serum_uric_acid",
            model = "M3", test = "overall_spline")
)

pool_contrast <- function(summaries, contrast) {
  q <- lapply(summaries, function(s) sum(contrast * s$coefficients))
  u <- lapply(summaries, function(s) as.numeric(t(contrast) %*% s$variance %*% contrast))
  p <- mitools::MIcombine(q, u, df.complete = min(vapply(
    summaries, `[[`, numeric(1), "df_complete"
  )))
  est <- as.numeric(p$coefficients)
  se <- sqrt(as.numeric(p$variance))
  df <- as.numeric(p$df)[[1]]
  crit <- stats::qt(0.975, df)
  c(estimate = est, ci_low = est - crit * se, ci_high = est + crit * se,
    p_value = 2 * stats::pt(abs(est / se), df, lower.tail = FALSE))
}

grid <- seq(spec$percentiles[["p01"]], spec$percentiles[["p99"]], length.out = 121)
ref_basis <- as.numeric(rcs_basis(spec$rcs_reference, spec$rcs_knots)[1, ])
names(ref_basis) <- c("rcs_linear", "rcs_nonlin1", "rcs_nonlin2")
coef_names <- names(rcs_summaries[[1]]$coefficients)
rcs_curve <- do.call(rbind, lapply(grid, function(value) {
  b <- as.numeric(rcs_basis(value, spec$rcs_knots)[1, ])
  names(b) <- names(ref_basis)
  L <- setNames(rep(0, length(coef_names)), coef_names)
  L[names(b)] <- b - ref_basis
  z <- pool_contrast(rcs_summaries, L)
  data.frame(database = "NHANES", outcome = "serum_uric_acid", model = "M3",
             upf_grams_pct = value, reference_pp = spec$rcs_reference,
             estimate = z[["estimate"]], ci_low = z[["ci_low"]],
             ci_high = z[["ci_high"]], p_value = z[["p_value"]])
}))

primary <- read.csv(main_path)
primary <- primary[primary$model == "M3", ]
primary_sens <- data.frame(
  database = "NHANES", outcome = "serum_uric_acid", model = "M3",
  module = "sensitivity", sensitivity = "primary_MICE",
  estimate = primary$estimate, ci_low = primary$ci_low, ci_high = primary$ci_high,
  p_value = primary$p_value, n_unweighted = primary$n_unweighted,
  successful_imputations = primary$successful_imputations
)

sens_rows <- list(primary_MICE = primary_sens)
for (scenario in c("exclude_extreme_energy", "trim_G_R0_p01_p99")) {
  ss <- lapply(bundles, function(z) z$sensitivity[[scenario]])
  row <- term_row(ss, sprintf("I(%s/10)", EXPOSURE_VAR), "serum_uric_acid",
                  "M3", "sensitivity")
  row$sensitivity <- scenario
  sens_rows[[scenario]] <- row
}

cc_vars <- unique(c(
  OUTCOMES$serum_uric_acid$variable, OUTCOMES$serum_uric_acid$domain,
  EXPOSURE_VAR, phase1_model_covariates("M3", "linear", "serum_uric_acid", FALSE)
))
cc_domain <- d$uric_acid_observed %in% TRUE & stats::complete.cases(d[, cc_vars])
cc_formula <- model_formula("serum_uric_acid", "M3", sprintf("I(%s/10)", EXPOSURE_VAR))
cc_summary <- compact_fit(d, cc_formula, "serum_uric_acid", cc_domain)
cc_row <- term_row(list(cc_summary), sprintf("I(%s/10)", EXPOSURE_VAR),
                   "serum_uric_acid", "M3", "sensitivity")
cc_row$sensitivity <- "complete_case"
sens_rows$complete_case <- cc_row
sensitivity <- dplyr::bind_rows(sens_rows)

base_design <- make_full_design(d)
base_design <- base_design[base_design$variables$uric_acid_observed %in% TRUE, ]
ua_mean <- survey::svymean(~uric_acid_dxc, base_design, na.rm = TRUE)
ua_quantiles <- survey::svyquantile(
  ~uric_acid_dxc, base_design, quantiles = c(0.01, 0.25, 0.50, 0.75, 0.99),
  ci = TRUE, na.rm = TRUE
)
descriptive <- data.frame(
  database = "NHANES", n_unweighted = nrow(base_design$variables),
  weighted_mean_mg_dl = as.numeric(stats::coef(ua_mean)),
  weighted_mean_se = as.numeric(sqrt(diag(stats::vcov(ua_mean)))),
  weighted_hyperuricemia_percent = as.numeric(
    stats::coef(survey::svymean(~hyperuricemia, base_design, na.rm = TRUE))
  ) * 100
)

write.csv(quartile, file.path(result_dir, "weighted_quartile_M1_M4.csv"), row.names = FALSE)
write.csv(trend, file.path(result_dir, "weighted_quartile_trend_M1_M4.csv"), row.names = FALSE)
write.csv(rcs_tests, file.path(result_dir, "RCS_M3_tests.csv"), row.names = FALSE)
write.csv(rcs_curve, file.path(result_dir, "RCS_M3_curve.csv"), row.names = FALSE)
write.csv(sensitivity, file.path(result_dir, "sensitivity_M3.csv"), row.names = FALSE)
write.csv(hyperuricemia, file.path(result_dir, "hyperuricemia_M1_M4_supplement.csv"), row.names = FALSE)
write.csv(descriptive, file.path(result_dir, "serum_uric_acid_weighted_descriptive.csv"), row.names = FALSE)
capture.output(ua_quantiles, file = file.path(qc_dir, "weighted_uric_acid_quantiles.txt"))
writeLines(capture.output(sessionInfo()), file.path(qc_dir, "sessionInfo.txt"))

stopifnot(nrow(quartile) == 12L, nrow(trend) == 4L, nrow(rcs_tests) == 2L,
          nrow(rcs_curve) == 121L, nrow(sensitivity) == 4L,
          nrow(hyperuricemia) == 4L, all(is.finite(rcs_tests$p_value)))
cat("NHANES_EXTENDED_OK\n")
print(trend)
print(rcs_tests[, c("test", "p_value")])
print(sensitivity[, c("sensitivity", "estimate", "ci_low", "ci_high", "p_value",
                      "n_unweighted")])
print(hyperuricemia[, c("model", "estimate", "ci_low", "ci_high", "p_value",
                       "n_unweighted", "events")])
)--------------------"

upf_source[["01_code/modules/SUA/06_nhanes_serum_uric_acid_subgroups.R"]] <- r"--------------------(options(stringsAsFactors = FALSE)

project_root <- "/storage/home/tmu2301/nhanes分析/18_formal_cholesterol_design_reanalysis_20260907/NHANES"
code_dir <- file.path(project_root, "02_代码", "正式五结局")
source(file.path(code_dir, "23_extension_common.R"), local = FALSE)
MI_M <- 50L
PHASE1_VERSION <- "2026-09-02.serum-uric-acid-extension-v1"
OUTCOMES <- list(
  serum_uric_acid = list(
    variable = "uric_acid_dxc", family = "gaussian",
    domain = "uric_acid_observed", effect_scale = "identity"
  )
)
EXPOSURES <- list(G_R0 = list(variable = "upf_gram_ratio_nonwater"))

extension_root <- Sys.getenv(
  "SERUM_UA_EXTENSION_ROOT",
  unset = "/storage/home/tmu2301/nhanes分析/18_formal_cholesterol_design_reanalysis_20260907/SUA"
)
result_dir <- file.path(extension_root, "03_扩展分析", "NHANES")
qc_dir <- file.path(extension_root, "05_质量控制", "NHANES_subgroup")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

d <- readRDS(file.path(project_root, "03_派生数据", "analysis_frame_pre_mice.rds"))
mi <- readRDS(file.path(project_root, "03_派生数据", "mi_uric_acid.rds"))
stopifnot(mi$m == MI_M)

subgroup_levels <- list(
  age_group = c("20-44", "45-64", "65+"),
  sex = c("Male", "Female"),
  smoking = c("Never", "Former", "Current"),
  hypertension = c("No", "Yes"),
  diabetes = c("No", "Yes")
)
subgroup_source <- c(
  age_group = "age10", sex = "sex", smoking = "smoking",
  hypertension = "hypertension", diabetes = "diabetes"
)

add_subgroup <- function(x, subgroup) {
  value <- switch(
    subgroup,
    age_group = ifelse(x$age < 45, "20-44", ifelse(x$age < 65, "45-64", "65+")),
    sex = as.character(x$sex), smoking = as.character(x$smoking),
    hypertension = ifelse(x$hypertension == 1, "Yes",
                          ifelse(x$hypertension == 0, "No", NA_character_)),
    diabetes = ifelse(x$diabetes == 1, "Yes",
                      ifelse(x$diabetes == 0, "No", NA_character_))
  )
  x$.subgroup <- factor(value, levels = subgroup_levels[[subgroup]])
  x
}

fit_one_subgroup_imputation <- function(i, subgroup) {
  x <- restore_completed_subframe(d, mi, i)
  x <- add_subgroup(x, subgroup)
  covars <- setdiff(
    phase1_model_covariates("M3", "linear", "serum_uric_acid", FALSE),
    if (subgroup == "age_group") character() else subgroup_source[[subgroup]]
  )
  exposure_term <- sprintf("I(%s/10)", EXPOSURE_VAR)
  form <- stats::reformulate(
    c(paste0(exposure_term, "*.subgroup"), covars), response = "uric_acid_dxc"
  )
  full <- make_full_design(x)
  domain <- full$variables$uric_acid_observed %in% TRUE
  needed <- all.vars(form)
  if (any(vapply(full$variables[needed], function(v) any(is.na(v[domain])), logical(1)))) {
    stop("Post-MICE missing values in subgroup model")
  }
  design <- full[domain, ]
  fit <- survey::svyglm(form, design = design, family = stats::gaussian(), na.action = na.fail)
  cf <- stats::coef(fit)
  vc <- stats::vcov(fit)
  strata <- as.character(design$variables$.strata_cycle)
  cluster <- interaction(
    design$variables$.strata_cycle, design$variables$.psu_cycle, drop = TRUE
  )
  support <- do.call(rbind, lapply(subgroup_levels[[subgroup]], function(level) {
    hit <- design$variables$.subgroup == level
    data.frame(
      subgroup = subgroup, level = level, imputation = i, n = sum(hit),
      n_psu = length(unique(cluster[hit])), n_strata = length(unique(strata[hit]))
    )
  }))
  if (any(support$n <= 0 | support$n_psu <= support$n_strata)) {
    stop("Subgroup support gate failed: ", subgroup)
  }
  list(coefficients = cf, variance = vc, df_complete = stats::df.residual(fit),
       support = support)
}

pool_scalar <- function(summaries, contrast) {
  q <- lapply(summaries, function(s) sum(contrast * s$coefficients))
  u <- lapply(summaries, function(s) as.numeric(t(contrast) %*% s$variance %*% contrast))
  p <- mitools::MIcombine(q, u, df.complete = min(vapply(
    summaries, `[[`, numeric(1), "df_complete"
  )))
  estimate <- as.numeric(p$coefficients)
  se <- sqrt(as.numeric(p$variance))
  df <- as.numeric(p$df)[[1]]
  critical <- stats::qt(0.975, df)
  data.frame(
    estimate = estimate, standard_error = se,
    ci_low = estimate - critical * se, ci_high = estimate + critical * se,
    p_value = 2 * stats::pt(abs(estimate / se), df, lower.tail = FALSE),
    df = df
  )
}

fit_subgroup <- function(subgroup) {
  summaries <- lapply(seq_len(MI_M), fit_one_subgroup_imputation, subgroup = subgroup)
  terms <- names(summaries[[1]]$coefficients)
  exposure_term <- grep("^I\\(upf_gram_ratio_nonwater/10\\)$", terms, value = TRUE)
  interaction_terms <- grep("I\\(upf_gram_ratio_nonwater/10\\):\\.subgroup|\\.subgroup.*:I\\(upf_gram_ratio_nonwater/10\\)",
                            terms, value = TRUE)
  stopifnot(length(exposure_term) == 1L,
            length(interaction_terms) == length(subgroup_levels[[subgroup]]) - 1L)
  effects <- do.call(rbind, lapply(seq_along(subgroup_levels[[subgroup]]), function(j) {
    L <- setNames(rep(0, length(terms)), terms)
    L[[exposure_term]] <- 1
    if (j > 1L) {
      label <- subgroup_levels[[subgroup]][[j]]
      hit <- interaction_terms[grepl(label, interaction_terms, fixed = TRUE)]
      stopifnot(length(hit) == 1L)
      L[[hit]] <- 1
    }
    z <- pool_scalar(summaries, L)
    transform(z, database = "NHANES", outcome = "serum_uric_acid",
              model = "M3", subgroup = subgroup,
              level = subgroup_levels[[subgroup]][[j]])
  }))
  if (length(interaction_terms) == 1L) {
    L <- setNames(rep(0, length(terms)), terms)
    L[[interaction_terms]] <- 1
    scalar <- pool_scalar(summaries, L)
    interaction <- data.frame(
      f_value = (scalar$estimate / scalar$standard_error)^2,
      df1 = 1, df2 = scalar$df, p_value = scalar$p_value,
      relative_increase_variance = NA_real_,
      complete_data_df = min(vapply(summaries, `[[`, numeric(1), "df_complete"))
    )
  } else {
    interaction <- d1_joint_test(summaries, interaction_terms)
  }
  interaction$database <- "NHANES"
  interaction$outcome <- "serum_uric_acid"
  interaction$model <- "M3"
  interaction$subgroup <- subgroup
  support <- do.call(rbind, lapply(summaries, `[[`, "support"))
  list(effects = effects, interaction = interaction, support = support)
}

workers <- min(n_cores_resolve(), length(subgroup_levels))
resource_preflight(workers, "nhanes_serum_uric_acid_subgroups")
cat("NHANES_SUBGROUP_START workers=", workers, " m=", MI_M, "\n", sep = "")
fits <- parallel::mclapply(names(subgroup_levels), fit_subgroup,
                           mc.cores = workers, mc.preschedule = TRUE)
if (any(vapply(fits, inherits, logical(1), "try-error"))) {
  print(fits[vapply(fits, inherits, logical(1), "try-error")])
  stop("NHANES subgroup analysis failed")
}
effects <- do.call(rbind, lapply(fits, `[[`, "effects"))
interactions <- do.call(rbind, lapply(fits, `[[`, "interaction"))
support <- do.call(rbind, lapply(fits, `[[`, "support"))
interactions$bh_q_interaction_within_outcome <- stats::p.adjust(
  interactions$p_value, method = "BH"
)

write.csv(effects, file.path(result_dir, "subgroup_effects.csv"), row.names = FALSE)
write.csv(interactions, file.path(result_dir, "subgroup_interactions.csv"), row.names = FALSE)
write.csv(support, file.path(qc_dir, "subgroup_support.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(qc_dir, "sessionInfo.txt"))
stopifnot(nrow(interactions) == 5L, nrow(support) == sum(lengths(subgroup_levels)) * MI_M,
          all(is.finite(interactions$p_value)))
cat("NHANES_SUBGROUP_OK\n")
print(interactions[, c("subgroup", "p_value", "bh_q_interaction_within_outcome")])
print(effects[, c("subgroup", "level", "estimate", "ci_low", "ci_high", "p_value")])
)--------------------"

upf_source[["01_code/modules/SUA/08_nhanes_serum_uric_acid_diet_nova.R"]] <- r"--------------------(options(stringsAsFactors = FALSE)

project_root <- "/storage/home/tmu2301/nhanes分析/18_formal_cholesterol_design_reanalysis_20260907/NHANES"
code_dir <- file.path(project_root, "02_代码", "正式五结局")
source(file.path(code_dir, "30_three_module_common.R"), local = FALSE)
MI_M <- 50L
PHASE1_VERSION <- "2026-09-02.serum-uric-acid-extension-v1"
THREE_MODULE_VERSION <- PHASE1_VERSION
OUTCOMES <- list(
  serum_uric_acid = list(
    variable = "uric_acid_dxc", family = "gaussian",
    domain = "uric_acid_observed", effect_scale = "identity"
  )
)
EXPOSURES <- list(G_R0 = list(variable = "upf_gram_ratio_nonwater"))

extension_root <- Sys.getenv(
  "SERUM_UA_EXTENSION_ROOT",
  unset = "/storage/home/tmu2301/nhanes分析/18_formal_cholesterol_design_reanalysis_20260907/SUA"
)
result_dir <- file.path(extension_root, "03_扩展分析", "NHANES")
qc_dir <- file.path(extension_root, "05_质量控制", "NHANES_diet_nova")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

d <- readRDS(file.path(
  project_root, "03_派生数据", "正式五结局_v1", "三模块探索",
  "analysis_frame_three_modules.rds"
))
mi <- readRDS(file.path(project_root, "03_派生数据", "mi_uric_acid.rds"))
stopifnot(nrow(d) == 66148L, mi$m == MI_M,
          all(c("hei2015", "nova1_share", "nova2_share", "nova3_share",
                "nova4_share", "total_classified_nonwater_g") %in% names(d)))

exposure_term <- sprintf("I(%s/10)", EXPOSURE_VAR)
m3_covars <- phase1_model_covariates("M3", "linear", "serum_uric_acid", FALSE)

fit_summary <- function(x, formula, required) {
  full <- make_full_design(x)
  domain <- full$variables$uric_acid_observed %in% TRUE &
    stats::complete.cases(full$variables[, required, drop = FALSE])
  design <- full[domain, ]
  fit <- survey::svyglm(formula, design = design, family = stats::gaussian(),
                        na.action = na.fail)
  survey_summary(fit, design$variables$SEQN)
}

diet_formula_reference <- stats::reformulate(
  c(exposure_term, m3_covars), response = "uric_acid_dxc"
)
diet_formula_adjusted <- stats::reformulate(
  c(exposure_term, "I(hei2015/10)", m3_covars), response = "uric_acid_dxc"
)
nova_terms <- paste0("I(nova", 2:4, "_share/10)")
nova_formula <- stats::reformulate(
  c(nova_terms, "I(total_classified_nonwater_g/1000)", m3_covars),
  response = "uric_acid_dxc"
)

fit_one <- function(i) {
  x <- restore_completed_subframe(d, mi, i)
  diet_required <- unique(c("hei2015", "uric_acid_dxc", "upf_gram_ratio_nonwater",
                            m3_covars))
  nova_required <- unique(c(
    paste0("nova", 1:4, "_share"), "total_classified_nonwater_g",
    "uric_acid_dxc", m3_covars
  ))
  list(
    diet_reference = fit_summary(x, diet_formula_reference, diet_required),
    diet_adjusted = fit_summary(x, diet_formula_adjusted, diet_required),
    nova = fit_summary(x, nova_formula, nova_required)
  )
}

workers <- n_cores_resolve()
resource_preflight(workers, "nhanes_serum_uric_acid_diet_nova")
cat("NHANES_DIET_NOVA_START workers=", workers, " m=", MI_M, "\n", sep = "")
bundles <- parallel::mclapply(seq_len(MI_M), fit_one,
                              mc.cores = workers, mc.preschedule = TRUE)
if (any(vapply(bundles, inherits, logical(1), "try-error"))) {
  print(bundles[vapply(bundles, inherits, logical(1), "try-error")])
  stop("NHANES diet/NOVA analysis failed")
}

reference <- lapply(bundles, `[[`, "diet_reference")
adjusted <- lapply(bundles, `[[`, "diet_adjusted")
nova <- lapply(bundles, `[[`, "nova")
ref_sample <- verify_summary_samples(reference, "HEI reference")
adj_sample <- verify_summary_samples(adjusted, "HEI adjusted")
stopifnot(identical(ref_sample, adj_sample))

diet_sets <- list(M3_same_sample = reference, M3_plus_HEI2015 = adjusted)
diet_rows <- do.call(rbind, lapply(names(diet_sets), function(model_name) {
  summaries <- diet_sets[[model_name]]
  z <- pooled_term_table(pool_svy_summaries(summaries), exposure_term, "identity")
  data.frame(
    database = "NHANES", outcome = "serum_uric_acid",
    model = model_name,
    estimate = z[["estimate"]], ci_low = z[["ci_low"]], ci_high = z[["ci_high"]],
    p_value = z[["p_value"]], standard_error_link = z[["standard_error_link"]],
    n_unweighted = ref_sample$n, sample_hash = ref_sample$hash,
    successful_imputations = MI_M
  )
}))

nova_sample <- verify_summary_samples(nova, "NOVA substitution")
contrasts <- list(
  NOVA4_vs_NOVA1 = list(terms = nova_terms[[3]], weights = 1, role = "primary"),
  NOVA2_vs_NOVA1 = list(terms = nova_terms[[1]], weights = 1, role = "descriptive"),
  NOVA3_vs_NOVA1 = list(terms = nova_terms[[2]], weights = 1, role = "descriptive"),
  NOVA4_vs_NOVA3 = list(terms = nova_terms[c(3, 2)], weights = c(1, -1), role = "supportive"),
  NOVA4_vs_NOVA2 = list(terms = nova_terms[c(3, 1)], weights = c(1, -1), role = "supportive")
)
nova_rows <- do.call(rbind, lapply(names(contrasts), function(name) {
  con <- contrasts[[name]]
  z <- pool_contrast(nova, con$terms, con$weights, name, "identity")
  data.frame(
    database = "NHANES", outcome = "serum_uric_acid", model = "M3",
    contrast = name, role = con$role, estimate = z[["estimate"]],
    ci_low = z[["ci_low"]], ci_high = z[["ci_high"]], p_value = z[["p_value"]],
    standard_error_link = z[["standard_error_link"]],
    n_unweighted = nova_sample$n, sample_hash = nova_sample$hash,
    successful_imputations = MI_M
  )
}))
nova_joint <- d1_joint_test(nova, nova_terms)
nova_joint$database <- "NHANES"
nova_joint$outcome <- "serum_uric_acid"
nova_joint$model <- "M3"
nova_joint$test <- "NOVA2_NOVA3_NOVA4_vs_NOVA1_joint_D1"

write.csv(diet_rows, file.path(result_dir, "HEI2015_adjustment.csv"), row.names = FALSE)
write.csv(nova_rows, file.path(result_dir, "NOVA_substitution_contrasts.csv"), row.names = FALSE)
write.csv(nova_joint, file.path(result_dir, "NOVA_substitution_joint_D1.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(qc_dir, "sessionInfo.txt"))
stopifnot(nrow(diet_rows) == 2L, nrow(nova_rows) == 5L,
          all(is.finite(diet_rows$estimate)), all(is.finite(nova_rows$estimate)))
cat("NHANES_DIET_NOVA_OK\n")
print(diet_rows)
print(nova_rows)
print(nova_joint[, c("test", "p_value")])
)--------------------"

upf_source[["01_code/modules/SUA/final_aligned_rcs.R"]] <- r"--------------------(options(stringsAsFactors = FALSE)

out_root <- Sys.getenv(
  "TABLE_REBUILD_OUT",
  unset = "/storage/home/tmu2301/nhanes分析/18_formal_cholesterol_design_reanalysis_20260907/SUA/final_alignment"
)
rcs_out <- file.path(out_root, "RCS_aligned")
dir.create(rcs_out, recursive = TRUE, showWarnings = FALSE)

nh_root <- "/storage/home/tmu2301/nhanes分析/18_formal_cholesterol_design_reanalysis_20260907/NHANES"
nh_code <- file.path(nh_root, "02_代码", "正式五结局")
MI_M <- 50L
EXTENSION_VERSION <- "2026-09-03.serum-uric-acid-rcs-aligned-v1"
OUTCOMES <- list(
  serum_uric_acid = list(
    variable = "uric_acid_dxc", family = "gaussian",
    domain = "uric_acid_observed", effect_scale = "identity"
  )
)
EXPOSURES <- list(
  G_R0 = list(variable = "upf_gram_ratio_nonwater", label = "Two-day non-water gram share")
)
source(file.path(nh_code, "23_extension_common.R"), local = FALSE)
MI_M <- 50L
OUTCOMES <- list(
  serum_uric_acid = list(
    variable = "uric_acid_dxc", family = "gaussian",
    domain = "uric_acid_observed", effect_scale = "identity"
  )
)
EXPOSURES <- list(
  G_R0 = list(variable = "upf_gram_ratio_nonwater", label = "Two-day non-water gram share")
)

nh_frame <- readRDS(file.path(nh_root, "03_派生数据", "analysis_frame_pre_mice.rds"))
nh_imp <- readRDS(file.path(nh_root, "03_派生数据", "mi_uric_acid.rds"))
stopifnot(nh_imp$m == MI_M)
nh_full <- make_full_design(nh_frame)
nh_base <- nh_full$variables$base_domain %in% TRUE &
  is.finite(nh_full$variables[[EXPOSURE_VAR]])
nh_design <- nh_full[nh_base, ]
nh_knots <- as.numeric(weighted_quantiles(
  nh_design, EXPOSURE_VAR, c(0.10, 0.50, 0.90)
))
nh_display <- as.numeric(weighted_quantiles(
  nh_design, EXPOSURE_VAR, c(0.01, 0.99)
))
stopifnot(length(nh_knots) == 3L, all(diff(nh_knots) > 0), all(diff(nh_display) > 0))

rcs_basis3 <- function(x, knots_pp) {
  knots <- as.numeric(knots_pp) / 10
  x <- as.numeric(x) / 10
  stopifnot(length(knots) == 3L, all(diff(knots) > 0))
  pos3 <- function(z) pmax(z, 0)^3
  k1 <- knots[[1]]; k2 <- knots[[2]]; k3 <- knots[[3]]
  nonlinear <- (
    pos3(x - k1) - pos3(x - k2) * (k3 - k1) / (k3 - k2) +
      pos3(x - k3) * (k2 - k1) / (k3 - k2)
  ) / ((k3 - k1)^2)
  cbind(rcs_linear = x, rcs_nonlin1 = nonlinear)
}

nh_compact_fit <- function(data, formula) {
  full <- make_full_design(data)
  domain <- full$variables$uric_acid_observed %in% TRUE
  model_vars <- all.vars(formula)
  stopifnot(!any(vapply(full$variables[model_vars], function(x) any(is.na(x[domain])), logical(1))))
  design <- full[domain, ]
  fit <- survey::svyglm(
    formula, design = design, family = stats::gaussian(), na.action = na.fail
  )
  list(
    coefficients = stats::coef(fit), variance = stats::vcov(fit),
    df_complete = stats::df.residual(fit), n = nrow(design$variables),
    design_df = survey::degf(design)
  )
}

nh_fit_one <- function(i) {
  x <- restore_completed_subframe(nh_frame, nh_imp, i)
  basis <- rcs_basis3(x[[EXPOSURE_VAR]], nh_knots)
  x$rcs_linear <- basis[, "rcs_linear"]
  x$rcs_nonlin1 <- basis[, "rcs_nonlin1"]
  formula <- stats::reformulate(
    c("rcs_linear", "rcs_nonlin1",
      phase1_model_covariates("M3", "linear", "serum_uric_acid", FALSE)),
    response = "uric_acid_dxc"
  )
  nh_compact_fit(x, formula)
}

nh_workers <- n_cores_resolve()
resource_preflight(nh_workers, "nhanes_serum_uric_acid_rcs_aligned")
nh_fits <- parallel::mclapply(
  seq_len(MI_M), function(i) {
    tryCatch(nh_fit_one(i), error = function(e) {
      structure(conditionMessage(e), class = "try-error")
    })
  }, mc.cores = nh_workers, mc.preschedule = TRUE
)
if (any(vapply(nh_fits, inherits, logical(1), "try-error"))) {
  print(unique(unlist(nh_fits[vapply(nh_fits, inherits, logical(1), "try-error")])))
  stop("NHANES aligned RCS fitting failed", call. = FALSE)
}
stopifnot(length(nh_fits) == MI_M)

nh_pooled <- pool_svy_summaries(nh_fits)
nh_nonlin_term <- pooled_term_table(nh_pooled, "rcs_nonlin1", "identity")
nh_nonlin_df <- as.numeric(nh_pooled$df[match("rcs_nonlin1", names(nh_pooled$coefficients))])
nh_nonlin_se <- as.numeric(nh_nonlin_term[["standard_error_link"]])
nh_nonlin_est <- as.numeric(nh_nonlin_term[["estimate_link"]])
nh_nonlinear <- data.frame(
  f_value = (nh_nonlin_est / nh_nonlin_se)^2,
  df1 = 1,
  df2 = nh_nonlin_df,
  p_value = as.numeric(nh_nonlin_term[["p_value"]]),
  relative_increase_variance = as.numeric(
    nh_pooled$missinfo[match("rcs_nonlin1", names(nh_pooled$coefficients))]
  ),
  complete_data_df = min(vapply(nh_fits, `[[`, numeric(1), "df_complete"))
)
nh_overall <- d1_joint_test(nh_fits, c("rcs_linear", "rcs_nonlin1"))
nh_tests <- rbind(
  transform(nh_nonlinear, database = "NHANES", outcome = "serum_uric_acid",
            model = "M3", test = "nonlinearity"),
  transform(nh_overall, database = "NHANES", outcome = "serum_uric_acid",
            model = "M3", test = "overall_spline")
)

nh_pool_contrast <- function(fits, contrast) {
  q <- lapply(fits, function(s) sum(contrast * s$coefficients))
  u <- lapply(fits, function(s) as.numeric(t(contrast) %*% s$variance %*% contrast))
  pooled <- mitools::MIcombine(
    q, u, df.complete = min(vapply(fits, `[[`, numeric(1), "df_complete"))
  )
  estimate <- as.numeric(pooled$coefficients)
  se <- sqrt(as.numeric(pooled$variance))
  df <- as.numeric(pooled$df)[[1]]
  crit <- stats::qt(0.975, df)
  c(
    estimate = estimate, ci_low = estimate - crit * se,
    ci_high = estimate + crit * se,
    p_value = 2 * stats::pt(abs(estimate / se), df, lower.tail = FALSE)
  )
}

nh_grid <- seq(nh_display[[1]], nh_display[[2]], length.out = 121)
nh_ref <- 10
nh_ref_basis <- as.numeric(rcs_basis3(nh_ref, nh_knots)[1, ])
names(nh_ref_basis) <- c("rcs_linear", "rcs_nonlin1")
nh_coef_names <- names(nh_fits[[1]]$coefficients)
nh_curve <- do.call(rbind, lapply(nh_grid, function(value) {
  basis <- as.numeric(rcs_basis3(value, nh_knots)[1, ])
  names(basis) <- names(nh_ref_basis)
  contrast <- stats::setNames(rep(0, length(nh_coef_names)), nh_coef_names)
  contrast[names(basis)] <- basis - nh_ref_basis
  z <- nh_pool_contrast(nh_fits, contrast)
  data.frame(
    database = "NHANES", outcome = "serum_uric_acid", model = "M3",
    upf_grams_pct = value, reference_pp = nh_ref,
    knot_10 = nh_knots[[1]], knot_50 = nh_knots[[2]], knot_90 = nh_knots[[3]],
    estimate = z[["estimate"]], ci_low = z[["ci_low"]],
    ci_high = z[["ci_high"]], p_value = z[["p_value"]]
  )
}))

write.csv(nh_tests,file.path(rcs_out,"NHANES_serum_uric_acid_RCS_tests.csv"),row.names=FALSE)
write.csv(nh_curve,file.path(rcs_out,"NHANES_serum_uric_acid_RCS_curve.csv"),row.names=FALSE)
saveRDS(nh_fits,file.path(rcs_out,"NHANES_serum_uric_acid_RCS_fit_summaries.rds"))
stopifnot(nh_ref==10,length(nh_knots)==3,all(nh_curve$reference_pp==10))
cat("ALIGNED_SUA_RCS_COMPLETED\n")
)--------------------"

upf_source[["Diagnostics/01_code/00_inputs.R"]] <- r"--------------------(# Read-only loaders for the frozen m=50 analyses. No original files are changed.
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           BLIS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1",
           KNHANES_ANALYSIS_MODE="formal", KNHANES_MI_M="50", KNHANES_MI_MAXIT="20",
           KNHANES_N_CORES="4")
base <- "/storage/home/tmu2301/nhanes分析"
.libPaths(c(file.path(base,"R_library"), .libPaths()))
audit_root <- "/storage/home/tmu2301/nhanes分析/18_formal_cholesterol_design_reanalysis_20260907/Diagnostics"
nh_root <- "/storage/home/tmu2301/nhanes分析/18_formal_cholesterol_design_reanalysis_20260907/NHANES"
kn_root <- file.path(base,"05_KNHANES_外部验证_2019_2024_20260820")
outcomes_order <- c("CKD","eGFR","UACR","serum_uric_acid","DKD","kidney_stones")

initialize_inputs <- function(database) {
  if(database=="NHANES") {
    fp <- file.path(nh_root,"03_派生数据","正式五结局_v1","analysis_frame_pre_mice.rds")
    frame <- readRDS(fp)
    uric_fp <- file.path(nh_root,"03_派生数据","analysis_frame_pre_mice.rds")
    uric <- readRDS(uric_fp)
    frozen <- c("SEQN","base_domain","upf_gram_ratio_nonwater","analysis_weight",
                "SDMVSTRA","SDMVPSU","age","age10","sex","race_ethnicity",
                "education","cycle_factor","mean_energy_kcal","smoking","alcohol",
                "physical_activity","bmi")
    stopifnot(all(vapply(frozen,function(v) identical(frame[[v]],uric[[v]]),logical(1))))
    mi_paths <- setNames(file.path(nh_root,"03_派生数据","正式五结局_v1",
                                  paste0("mi_",names(OUTCOMES),".rds")),names(OUTCOMES))
    mi_paths["serum_uric_acid"] <- file.path(nh_root,"03_派生数据","mi_uric_acid.rds")
    mi <- lapply(mi_paths,readRDS)
    stopifnot(all(vapply(mi,function(x)x$m==50 && x$mids$m==50,logical(1))))
    specs <- OUTCOMES
    specs$serum_uric_acid <- list(variable="uric_acid_dxc",family="gaussian",
                                 domain="uric_acid_observed",effect_scale="identity")
    getter <- function(outcome,i) {
      d <- if(outcome=="serum_uric_acid") uric else frame
      restore_completed_subframe(d,mi[[outcome]],i)
    }
    design_maker <- make_full_design
    covars <- MODEL_COVARIATES$M3
    paths <- c(fp,uric_fp,unname(mi_paths))
    expvar <- "upf_gram_ratio_nonwater"; agevar <- "age10"
  } else {
    stopifnot(ANALYSIS_MODE=="formal",MI_M==50L)
    inputs <- load_model_inputs(); frame <- inputs$frame; mi <- inputs$imp
    old_fp <- file.path(kn_root,"03_派生数据","FORMAL.bak_R9_20260827",
                        "N4_PERSON_FRAME_20260821T121405241_P1029320","analysis_frame_pre_mice.rds")
    old <- readRDS(old_fp); ix <- match(as.character(frame$ID),as.character(old$ID))
    stopifnot(!anyNA(ix),!anyDuplicated(frame$ID),!anyDuplicated(old$ID),mi$m==50L)
    frozen <- c("ID","person_uid","base_domain","upf_grams_pct","analysis_weight",
                "psu_uid","strata_uid","age","sex","education","survey_year",
                "energy_kcal","smoking","alcohol","physical_activity","bmi")
    frozen <- intersect(frozen,intersect(names(frame),names(old)))
    frozen_check <- vapply(frozen,function(v) isTRUE(all.equal(frame[[v]],old[[v]][ix],
                          tolerance=0,check.attributes=TRUE)),logical(1))
    if(!all(frozen_check)) {
      print(lapply(frozen[!frozen_check],function(v) list(variable=v,
                    comparison=all.equal(frame[[v]],old[[v]][ix],tolerance=0,check.attributes=TRUE))))
      stop("Frozen KNHANES fields differ")
    }
    frame$uric_acid_mg_dl <- old$uric_acid_mg_dl[ix]
    frame$uric_acid_observed <- frame$base_domain %in% TRUE & is.finite(frame$uric_acid_mg_dl)
    frame <- derive_fixed_exposure_categories(frame)
    specs <- OUTCOMES; specs$UACR <- specs$log_uacr; specs$log_uacr <- NULL
    specs$serum_uric_acid <- list(variable="uric_acid_mg_dl",family="gaussian",domain="uric_acid_observed")
    getter <- function(outcome,i) completed_analysis_dataset(frame,mi,i)
    design_maker <- make_survey_design; covars <- MODEL_COVARIATES$M3
    paths <- c(inputs$frame_path,inputs$mice_path,old_fp)
    expvar <- "upf_grams_pct"; agevar <- "age"
  }
  specs <- specs[intersect(outcomes_order,names(specs))]
  for(o in names(specs)) specs[[o]]$effect_scale <- if(o=="UACR") "log_response_percent" else
    if(specs[[o]]$family %in% c("binomial","quasibinomial")) "odds_ratio" else "identity"
  list(frame=frame,getter=getter,design=design_maker,specs=specs,covars=covars,
       expvar=expvar,agevar=agevar,paths=paths,m=50L)
}

load_nhanes_moisture <- function(frame) {
  parts <- lapply(ANALYTIC_CYCLES,function(cycle) {
    a <- as.data.frame(haven::read_xpt(diet_file(cycle,1,"tot")))
    b <- as.data.frame(haven::read_xpt(diet_file(cycle,2,"tot")))
    stopifnot(!anyDuplicated(a$SEQN),!anyDuplicated(b$SEQN))
    ac <- c("SEQN","DR1TMOIS","DR1DRSTZ","DR1_320Z")
    bc <- c("SEQN","DR2TMOIS","DR2DRSTZ","DR2_320Z")
    a <- a[intersect(ac,names(a))]; b <- b[intersect(bc,names(b))]
    z <- merge(a,b,by="SEQN",all=TRUE,sort=FALSE); z$cycle <- cycle; z
  })
  z <- do.call(rbind,parts); stopifnot(!anyDuplicated(z$SEQN))
  ix <- match(frame$SEQN,z$SEQN)
  z <- z[ix,]; z$SEQN <- frame$SEQN
  z$moisture_valid <- z$DR1DRSTZ %in% 1 & z$DR2DRSTZ %in% 1 &
    is.finite(z$DR1TMOIS) & is.finite(z$DR2TMOIS)
  z$moisture_kg_day <- ifelse(z$moisture_valid,(z$DR1TMOIS+z$DR2TMOIS)/2000,NA_real_)
  z
}
)--------------------"

upf_source[["Diagnostics/01_code/01_preflight.R"]] <- r"--------------------(args <- commandArgs(TRUE); stopifnot(length(args)==1L,args[1] %in% c("NHANES","KNHANES"))
script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(script),"00_inputs.R"))
database <- args[1]
if(database=="NHANES") source(file.path(nh_root,"02_代码","正式五结局","04_models.R")) else
  source(file.path(kn_root,"04_代码","00_model_functions.R"))
inp <- initialize_inputs(database)
out <- file.path(audit_root,"02_preflight"); dir.create(out,recursive=TRUE,showWarnings=FALSE)
manifest <- data.frame(path=inp$paths,size=file.info(inp$paths)$size,
  sha256=vapply(inp$paths,function(p)digest::digest(file=p,algo="sha256"),character(1)))
write.csv(manifest,file.path(out,paste0(database,"_input_manifest.csv")),row.names=FALSE)
rows <- lapply(names(inp$specs),function(o) {
  d <- inp$getter(o,1); s <- inp$specs[[o]]; hit <- d[[s$domain]] %in% TRUE
  stopifnot(!anyNA(d[hit,c(s$variable,inp$expvar,inp$covars)]))
  data.frame(database=database,outcome=o,n=sum(hit),events=if(s$effect_scale=="odds_ratio")
    sum(d[[s$variable]][hit]==1) else NA_integer_,m=inp$m,age_min=min(d$age[hit]),
    age_max=max(d$age[hit]),formula=paste(deparse(reformulate(c(".upf10",inp$covars),s$variable)),collapse=" "))
})
tab <- do.call(rbind,rows); print(tab,row.names=FALSE)
write.csv(tab,file.path(out,paste0(database,"_domain_audit.csv")),row.names=FALSE)
if(database=="NHANES") {
  z <- load_nhanes_moisture(inp$frame); hit <- inp$frame$stone_observed %in% TRUE
  summary <- data.frame(cycle=ANALYTIC_CYCLES,n=vapply(ANALYTIC_CYCLES,function(c)sum(hit & inp$frame$cycle==c),integer(1)),
    complete_moisture=vapply(ANALYTIC_CYCLES,function(c)sum(hit & inp$frame$cycle==c & z$moisture_valid),integer(1)))
  print(summary); print(quantile(z$moisture_kg_day[hit],na.rm=TRUE));
  write.csv(summary,file.path(out,"NHANES_moisture_availability.csv"),row.names=FALSE)
  saveRDS(z,file.path(out,"NHANES_moisture_linked_server_only.rds"))
}
writeLines(capture.output(sessionInfo()),file.path(out,paste0(database,"_sessionInfo.txt")))
cat(database,"PREFLIGHT_PASSED\n")
)--------------------"

upf_source[["Diagnostics/01_code/02_diagnostic_models.R"]] <- r"--------------------(# Prespecified diagnostic models; never overwrites a published/frozen analysis.
args <- commandArgs(TRUE)
stopifnot(length(args)>=2L,args[1] %in% c("NHANES","KNHANES"),args[2] %in% c("smoke","formal"))
script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(script),"00_inputs.R"))
database <- args[1]; mode <- args[2]
if(database=="NHANES") source(file.path(nh_root,"02_代码","正式五结局","04_models.R")) else
  source(file.path(kn_root,"04_代码","00_model_functions.R"))
inp <- initialize_inputs(database)
output <- file.path(audit_root,"03_models",mode,database)
dir.create(output,recursive=TRUE,showWarnings=FALSE)
workers <- if(mode=="formal") n_cores_resolve() else 1L
indices <- if(mode=="formal") 1:50 else 1:2
source_hash <- digest::digest(list(readLines(script,warn=FALSE),
                   readLines(file.path(dirname(script),"00_inputs.R"),warn=FALSE)),algo="sha256")
source_manifest <- read.csv(file.path(audit_root,"02_preflight",paste0(database,"_input_manifest.csv")))
stopifnot(all(vapply(seq_len(nrow(source_manifest)),function(j)
  identical(digest::digest(file=source_manifest$path[j],algo="sha256"),source_manifest$sha256[j]),logical(1))))
input_hash <- digest::digest(source_manifest,algo="sha256")
if(database=="NHANES") {
  moisture <- readRDS(file.path(audit_root,"02_preflight","NHANES_moisture_linked_server_only.rds"))
  stopifnot(identical(moisture$SEQN,inp$frame$SEQN))
}
weighted_knots <- function(x,w) {
  keep <- is.finite(x)&is.finite(w)&w>0; x<-x[keep];w<-w[keep]
  ix<-order(x);x<-x[ix];w<-w[ix];cdf<-cumsum(w)/sum(w)
  vapply(c(.1,.5,.9),function(p)x[which(cdf>=p)[1]],numeric(1))
}
rcs_nonlinear <- function(x,k) {
  stopifnot(length(k)==3L,all(diff(k)>0))
  (pmax(x-k[1],0)^3-pmax(x-k[2],0)^3*(k[3]-k[1])/(k[3]-k[2])+
     pmax(x-k[3],0)^3*(k[2]-k[1])/(k[3]-k[2]))/(k[3]-k[1])^2/10
}
compact <- function(fit,design,spec,formula) {
  b<-coef(fit);V<-vcov(fit);df<-as.numeric(df.residual(fit))
  stopifnot(all(is.finite(b)),all(is.finite(V)),all(diag(V)>=0),is.finite(df),df>0,
            identical(names(b),rownames(V)),identical(names(b),colnames(V)))
  y<-design$variables[[spec$variable]]
  list(coefficients=b,variance=V,df_complete=df,n=nrow(design$variables),
       events=if(spec$effect_scale=="odds_ratio")sum(y==1) else NA_integer_,
       design_df=survey::degf(design),formula=paste(deparse(formula),collapse=" "))
}
knots_rows <- list()
for(o in names(inp$specs)) {
  spec<-inp$specs[[o]]; d0<-inp$getter(o,1);hit<-d0[[spec$domain]] %in% TRUE
  k<-weighted_knots(d0$age[hit],d0$analysis_weight[hit])
  knots_rows[[o]]<-data.frame(database=database,outcome=o,knot_p10=k[1],knot_p50=k[2],knot_p90=k[3])
  f_main<-reformulate(c(".upf10",inp$covars),spec$variable)
  f_rcs<-reformulate(c(".upf10",inp$covars,".age_nonlinear"),spec$variable)
  f_group_old<-reformulate(c(".upf10*.age_group",setdiff(inp$covars,inp$agevar)),spec$variable)
  f_group_age<-reformulate(c(".upf10*.age_group",inp$covars),spec$variable)
  formulas<-list(M3_replicated=f_main,M3_age_rcs3=f_rcs,
                 age_group_retains_age=f_group_age)
  for(batch_start in seq(1,length(indices),by=workers)) {
    batch<-indices[seq(batch_start,min(batch_start+workers-1L,length(indices)))]
    results<-parallel::mclapply(batch,function(i) {
      path<-file.path(output,sprintf("%s_imp%02d.rds",o,i))
      if(file.exists(path)) {
        cached<-readRDS(path)
        if(identical(cached$source_hash,source_hash)&&identical(cached$input_hash,input_hash)&&
           identical(cached$outcome,o)&&cached$imputation==i&&cached$mode==mode)
          return(list(ok=TRUE,reused=TRUE,path=path))
        stop("Cache lineage mismatch: ",path)
      }
      d<-inp$getter(o,i);d$.upf10<-d[[inp$expvar]]/10
      d$.age_nonlinear<-rcs_nonlinear(d$age,k)
      d$.age_group<-factor(ifelse(d$age<45,"G1",ifelse(d$age<65,"G2","G3")),levels=c("G1","G2","G3"))
      if(o=="kidney_stones")d$.moisture_kg_day<-moisture$moisture_kg_day
      full<-inp$design(d);base_hit<-full$variables[[spec$domain]] %in% TRUE
      fam<-if(spec$effect_scale=="odds_ratio")quasibinomial() else gaussian()
      summaries<-lapply(names(formulas),function(model) {
        valid<-base_hit
        if(model %in% c("M3_moisture_complete","M3_plus_moisture"))
          valid<-valid & is.finite(full$variables$.moisture_kg_day)
        ds<-full[valid,]; f<-formulas[[model]]
        stopifnot(nrow(ds$variables)>0,!anyNA(ds$variables[all.vars(f)]))
        fit<-survey::svyglm(f,design=ds,family=fam,na.action=na.fail)
        compact(fit,ds,spec,f)
      });names(summaries)<-names(formulas)
      support<-table(d$.age_group[d[[spec$domain]] %in% TRUE])
      obj<-list(database=database,outcome=o,imputation=i,mode=mode,m=inp$m,
                source_hash=source_hash,input_hash=input_hash,knots=k,support=support,
                effect_scale=spec$effect_scale,summaries=summaries)
      temporary<-paste0(path,".tmp-",Sys.getpid());saveRDS(obj,temporary,compress="gzip")
      stopifnot(file.rename(temporary,path))
      list(ok=TRUE,reused=FALSE,path=path)
    },mc.cores=workers,mc.preschedule=TRUE)
    if(any(vapply(results,inherits,logical(1),"try-error")))stop("Diagnostic batch failed")
    cat(format(Sys.time()),database,o,"completed",max(batch),"/",length(indices),"imputations\n")
    flush.console()
  }
}
write.csv(do.call(rbind,knots_rows),file.path(output,"fixed_age_knots.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(output,"sessionInfo.txt"))
writeLines(c(paste("source_sha256",source_hash),paste("input_manifest_sha256",input_hash),
             paste("mode",mode),paste("m",length(indices)),paste("workers",workers),"status PASSED"),
           file.path(output,"COMPLETED.txt"))
cat(database,mode,"ALL_DIAGNOSTIC_FITS_PASSED\n")
)--------------------"

upf_source[["Diagnostics/01_code/03_pool_diagnostics.R"]] <- r"--------------------(script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(script),"00_inputs.R"))
dir.create(file.path(audit_root,"04_summary"),recursive=TRUE,showWarnings=FALSE)
pool_contrast<-function(summaries,L,scale) {
  estimates<-lapply(summaries,function(s)setNames(as.numeric(L%*%s$coefficients),"effect"))
  variances<-lapply(summaries,function(s)matrix(as.numeric(L%*%s$variance%*%L),1,1,dimnames=list("effect","effect")))
  p<-mitools::MIcombine(estimates,variances,df.complete=min(vapply(summaries,`[[`,numeric(1),"df_complete")))
  beta<-unname(p$coefficients[1]);se<-sqrt(p$variance[1,1]);df<-p$df[1]
  ci<-beta+c(-1,1)*qt(.975,df)*se
  transform<-switch(scale,odds_ratio=exp,log_response_percent=function(x)100*expm1(x),identity=identity)
  data.frame(estimate=transform(beta),ci_low=transform(ci[1]),ci_high=transform(ci[2]),
             beta_link=beta,se_link=se,df=df,p_value=2*pt(abs(beta/se),df,lower.tail=FALSE),
             n=summaries[[1]]$n,events=summaries[[1]]$events,m=length(summaries))
}
d1<-function(summaries,terms) {
  qhat<-do.call(cbind,lapply(summaries,function(s)unname(s$coefficients[terms])))
  k<-nrow(qhat);m<-ncol(qhat)
  uhat<-array(unlist(lapply(summaries,function(s)s$variance[terms,terms,drop=FALSE]),use.names=FALSE),c(k,k,m))
  qbar<-rowMeans(qhat);ubar<-matrix(apply(uhat,c(1,2),mean),k,k)
  b<-if(k==1L)matrix(var(as.numeric(qhat)),1,1) else cov(t(qhat))
  r<-(1+1/m)*sum(diag(b%*%solve(ubar)))/k;tval<-k*(m-1)
  fvalue<-as.numeric(t(qbar)%*%solve((1+r)*ubar)%*%qbar/k)
  dfcom<-min(vapply(summaries,`[[`,numeric(1),"df_complete"))
  a<-r*tval/(tval-2);vstar<-((dfcom+1)/(dfcom+3))*dfcom
  c0<-1/(tval-4);c1<-vstar-2*(1+a);c2<-vstar-4*(1+a)
  z<-1/c2+c0*(a^2*c1/((1+a)^2*c2))+
    c0*(8*a^2*c1/((1+a)*c2^2)+4*a^2/((1+a)*c2))+
    c0*(4*a^2/(c2*c1)+16*a^2*c1/c2^3)+c0*(8*a^2/c2^2)
  df2<-4+1/z
  data.frame(f_value=fvalue,df1=k,df2=df2,p_interaction=pf(fvalue,k,df2,lower.tail=FALSE))
}
rows<-list();joints<-list();nonlinear<-list();row_i<-0L;joint_i<-0L;non_i<-0L
for(database in c("NHANES")) {
  folder<-file.path(audit_root,"03_models","formal",database)
  stopifnot(file.exists(file.path(folder,"COMPLETED.txt")))
  outcomes<-if(database=="NHANES")outcomes_order else setdiff(outcomes_order,"kidney_stones")
  for(o in outcomes) {
    objects<-lapply(1:50,function(i)readRDS(file.path(folder,sprintf("%s_imp%02d.rds",o,i))))
    stopifnot(all(vapply(objects,function(z)z$outcome==o&&z$database==database&&z$mode=="formal",logical(1))))
    for(model in names(objects[[1]]$summaries)) {
      summaries<-lapply(objects,function(z)z$summaries[[model]])
      terms<-names(summaries[[1]]$coefficients)
      stopifnot(all(vapply(summaries,function(s)identical(names(s$coefficients),terms),logical(1))))
      labels<-if(grepl("^age_group",model))c("20-44","45-64","65+") else "Overall"
      for(level in seq_along(labels)) {
        L<-setNames(rep(0,length(terms)),terms);L[".upf10"]<-1
        if(level>1L) {
          term<-terms[grepl(".upf10",terms,fixed=TRUE)&grepl(paste0(".age_groupG",level),terms,fixed=TRUE)]
          stopifnot(length(term)==1L);L[term]<-1
        }
        row_i<-row_i+1L
        rows[[row_i]]<-cbind(data.frame(database=database,outcome=o,model=model,level=labels[level],
                     effect_scale=objects[[1]]$effect_scale),pool_contrast(summaries,L,objects[[1]]$effect_scale))
      }
      if(grepl("^age_group",model)) {
        terms_int<-terms[grepl(".upf10",terms,fixed=TRUE)&grepl(".age_group",terms,fixed=TRUE)]
        stopifnot(length(terms_int)==2L);joint_i<-joint_i+1L
        joints[[joint_i]]<-cbind(data.frame(database=database,outcome=o,model=model),d1(summaries,terms_int))
      }
      if(model=="M3_age_rcs3") {
        L<-setNames(rep(0,length(terms)),terms);L[".age_nonlinear"]<-1;non_i<-non_i+1L
        nonlinear[[non_i]]<-cbind(data.frame(database=database,outcome=o),pool_contrast(summaries,L,"identity"))
      }
    }
  }
}
out<-file.path(audit_root,"04_summary")
write.csv(do.call(rbind,rows),file.path(out,"diagnostic_effects.csv"),row.names=FALSE)
write.csv(do.call(rbind,joints),file.path(out,"age_interaction_tests.csv"),row.names=FALSE)
write.csv(do.call(rbind,nonlinear),file.path(out,"age_nonlinearity_tests.csv"),row.names=FALSE)
print(subset(do.call(rbind,rows),level=="Overall"),row.names=FALSE)
cat("ALL_DIAGNOSTIC_POOLING_PASSED\n")
)--------------------"

upf_source[["01_code/05_formal_mi_worker.R"]] <- r"--------------------(Sys.setenv(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
suppressPackageStartupMessages({library(mice);library(jsonlite);library(digest)})
root<-'/storage/home/tmu2301/nhanes分析/18_formal_cholesterol_design_reanalysis_20260907'
a<-commandArgs(TRUE);id<-a[1];chain<-as.integer(a[2]);target<-as.integer(a[3])
manifest<-read_json(file.path(root,'03_qc/mi_job_manifest.json'),simplifyVector=TRUE)
spec<-manifest[manifest$id==id,,drop=FALSE];stopifnot(nrow(spec)==1,chain>=1,chain<=50,target == 20L)
tag<-sprintf('%s_C%02d',id,chain);t0<-Sys.time()
for(sub in c('02_mi/chains','02_mi/checkpoints','03_qc/optimizers','05_logs/mi'))dir.create(file.path(root,sub),recursive=TRUE,showWarnings=FALSE)
obj<-readRDS(spec$template);d<-readRDS(spec$frame)
z<-obj$mids$data;pred<-obj$mids$predictorMatrix;method<-obj$mids$method
stopifnot(identical(as.numeric(z$SEQN),as.numeric(d$SEQN[obj$original_row_index])))
z$high_cholesterol<-factor(d$high_cholesterol[obj$original_row_index],levels=c(0,1))
z$mi_psu_code<-factor(z$mi_psu_code);pred[,c('mi_stratum_code','cycle_factor')]<-0L
unscaled<-z;centers<-scales<-c()
for(v in setdiff(names(z)[vapply(z,is.numeric,logical(1))],c('SEQN','mi_stratum_code','mi_psu_code'))) {
 s<-sd(z[[v]],na.rm=TRUE);if(!is.finite(s)||s==0)next
 centers[v]<-mean(z[[v]],na.rm=TRUE);scales[v]<-s;z[[v]]<-(z[[v]]-centers[v])/s
}
fingerprint<-digest(list(z,pred,method,centers,scales,protocol='psu_scaled_exact_pmm_v2',nnet.maxit=1000,nnet.MaxNWts=10000,mice_version=as.character(packageVersion('mice')),nnet_version=as.character(packageVersion('nnet'))),algo='sha256')
.fit_audit<-list()
trace('multinom',where=asNamespace('nnet'),print=FALSE,exit=quote({
 .ret<-returnValue();.GlobalEnv$.fit_audit[[length(.GlobalEnv$.fit_audit)+1L]]<-data.frame(model='multinom',converged=.ret$convergence==0,iterations=NA_integer_,n_predictors=.ret$n[1],n_classes=.ret$n[3],objective=.ret$value)
 if(.ret$convergence!=0)stop('Multinomial optimizer failed; formal candidate stopped')
}))
trace('glm.fit',where=asNamespace('stats'),print=FALSE,exit=quote({
 .ret<-returnValue();.GlobalEnv$.fit_audit[[length(.GlobalEnv$.fit_audit)+1L]]<-data.frame(model='glm.fit',converged=isTRUE(.ret$converged),iterations=.ret$iter,n_predictors=length(.ret$coefficients),n_classes=2,objective=.ret$deviance)
 if(!isTRUE(.ret$converged) || any(!is.finite(.ret$coefficients)))stop('Binary optimizer failed or rank deficient; formal candidate stopped')
}))
state_path<-file.path(root,'05_logs/mi',paste0(tag,'.json'))
write_state<-function(status,it,error=NULL)write_json(list(tag=tag,id=id,chain=chain,status=status,iteration=it,target=target,seconds=as.numeric(difftime(Sys.time(),t0,units='secs')),error=error),state_path,auto_unbox=TRUE,pretty=TRUE)
chain_path<-file.path(root,'02_mi/chains',paste0(tag,'.rds'));imp<-NULL
if(file.exists(chain_path)) {
 previous<-readRDS(chain_path);stopifnot(identical(previous$fingerprint,fingerprint));imp<-previous$mids;rm(previous)
}
restore_object<-function(imp) {
 restored<-imp;restored$data<-unscaled
 for(v in names(scales)) {
  if(!is.null(restored$imp[[v]]) && nrow(restored$imp[[v]])>0) {
   if(identical(unname(method[v]),'pmm')) {
    observed<-!is.na(z[[v]]);donors_scaled<-z[[v]][observed];donors_original<-unscaled[[v]][observed]
    for(k in seq_len(ncol(restored$imp[[v]]))) {
     ix<-match(restored$imp[[v]][[k]],donors_scaled);stopifnot(!anyNA(ix));restored$imp[[v]][[k]]<-donors_original[ix]
    }
   } else restored$imp[[v]]<-restored$imp[[v]]*scales[v]+centers[v]
  }
  if(v %in% dimnames(restored$chainMean)[[1]]) {
   restored$chainMean[v,,]<-restored$chainMean[v,,]*scales[v]+centers[v]
   restored$chainVar[v,,]<-restored$chainVar[v,,]*scales[v]^2
  }
 }
 completed<-complete(restored);scaled_complete<-complete(imp)
 for(v in names(unscaled)) {
  seen<-!is.na(unscaled[[v]]);stopifnot(identical(completed[[v]][seen],unscaled[[v]][seen]))
  if(v %in% names(scales))stopifnot(isTRUE(all.equal(as.numeric(completed[[v]]),as.numeric(scaled_complete[[v]]*scales[v]+centers[v]),tolerance=1e-10)))
 }
 stopifnot(!anyNA(completed[names(method)[nzchar(method)]]))
 ans<-obj;ans$mids<-restored;ans$m<-1L;ans$maxit<-imp$iteration;ans$predictor_matrix<-restored$predictorMatrix
 ans$standardization<-list(center=centers,scale=scales);ans$source_mids_sha256<-digest(file=spec$template,algo='sha256');ans$continuation_rng<-NULL
 ans$revision<-'20260907 cholesterol routing and categorical PSU in MI';ans$input_fingerprint<-fingerprint
 ans
}
tryCatch({
 start<-if(is.null(imp))1L else imp$iteration+1L
 if(start<=target)for(it in seq.int(start,target)) {
  .fit_audit<-list()
  if(is.null(imp))imp<-mice(z,m=1,maxit=1,method=method,predictorMatrix=pred,seed=202609070L+chain,printFlag=FALSE,nnet.MaxNWts=10000,nnet.maxit=1000)
  else imp<-mice.mids(imp,maxit=1,printFlag=FALSE,nnet.MaxNWts=10000,nnet.maxit=1000)
  au<-do.call(rbind,.fit_audit);au$iteration<-it
  ap<-file.path(root,'03_qc/optimizers',paste0(tag,'.csv'));exists<-file.exists(ap)
  write.table(au,ap,sep=',',row.names=FALSE,col.names=!exists,append=exists)
  saveRDS(list(mids=imp,fingerprint=fingerprint),paste0(chain_path,'.tmp'),compress=FALSE);stopifnot(file.rename(paste0(chain_path,'.tmp'),chain_path))
  if(it==1 || it %% 20==0)saveRDS(restore_object(imp),file.path(root,'02_mi/checkpoints',sprintf('%s_I%03d.rds',tag,it)),compress='gzip')
  write_state('running',it)
 }
 stopifnot(imp$iteration==target)
 saveRDS(restore_object(imp),file.path(root,'02_mi/checkpoints',sprintf('%s_I%03d.rds',tag,target)),compress='gzip')
 write_state('completed',imp$iteration)
},error=function(e){
 if(length(.fit_audit))write.csv(do.call(rbind,.fit_audit),file.path(root,'03_qc/optimizers',paste0(tag,'_failed.csv')),row.names=FALSE)
 write_state('failed',if(is.null(imp))0 else imp$iteration,conditionMessage(e));stop(e)
})
)--------------------"

upf_source[["01_code/09_pool_checkpoint.R"]] <- r"--------------------(Sys.setenv(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
suppressPackageStartupMessages({library(mice);library(jsonlite);library(posterior)})
root<-'/storage/home/tmu2301/nhanes分析/18_formal_cholesterol_design_reanalysis_20260907'
a<-commandArgs(TRUE);id<-a[1];iteration<-as.integer(a[2]);nchains<-if(length(a)>=3)as.integer(a[3]) else 50L
manifest<-read_json(file.path(root,'03_qc/mi_job_manifest.json'),simplifyVector=TRUE);spec<-manifest[manifest$id==id,,drop=FALSE];stopifnot(nrow(spec)==1,nchains>=2)
dir.create(file.path(root,'03_qc/mi_diagnostics'),recursive=TRUE,showWarnings=FALSE)
paths<-file.path(root,'02_mi/checkpoints',sprintf('%s_C%02d_I%03d.rds',id,seq_len(nchains),iteration));stopifnot(all(file.exists(paths)))
combined<-NULL;obj<-NULL
for(p in paths){x<-readRDS(p);if(is.null(combined)){obj<-x;combined<-x$mids}else{stopifnot(identical(x$input_fingerprint,obj$input_fingerprint));combined<-ibind(combined,x$mids)}}
stopifnot(combined$m==nchains,combined$iteration==iteration)
obj$mids<-combined;obj$m<-nchains;obj$maxit<-iteration
targets<-names(combined$method)[nzchar(combined$method)]
diagnostics<-list();ranks<-list()
for(param in c('mean','sd')){
 cv<-suppressWarnings(mice::convergence(combined,parameter=param));cv<-cv[cv$vrb %in% targets,];cv$id<-id;cv$parameter<-param;diagnostics[[param]]<-cv
}
for(field in c('chainMean','chainVar'))for(v in targets){
 x<-combined[[field]][v,seq.int(max(1L,floor(iteration/2)+1L),iteration),,drop=FALSE];dim(x)<-c(ceiling(iteration/2),nchains)
 ranks[[paste(field,v)]]<-data.frame(id=id,variable=v,parameter=field,n_missing=combined$nmis[v],iteration=iteration,m=nchains,rhat=suppressWarnings(posterior::rhat(x)),ess_bulk=suppressWarnings(posterior::ess_bulk(x)),ess_tail=suppressWarnings(posterior::ess_tail(x)))
}
write.csv(do.call(rbind,diagnostics),file.path(root,'03_qc/mi_diagnostics',sprintf('%s_I%03d_convergence.csv',id,iteration)),row.names=FALSE)
write.csv(do.call(rbind,ranks),file.path(root,'03_qc/mi_diagnostics',sprintf('%s_I%03d_rank.csv',id,iteration)),row.names=FALSE)
source(file.path(root,'NHANES/02_代码/正式五结局/04_models.R'))
MI_M<-nchains
if(spec$outcome=='serum_uric_acid')OUTCOMES<-list(serum_uric_acid=list(variable='uric_acid_dxc',family='gaussian',domain='uric_acid_observed',effect_scale='identity'))
d<-readRDS(spec$frame);grid<-make_outcome_grid(spec$outcome,d)
if(spec$branch=='day1')grid<-grid[grid$model=='M3',]
fits<-lapply(seq_len(nchains),fit_one_imputation,d=d,mi_object=obj,grid=grid,include_pir=FALSE)
res<-pool_outcome_grid(fits,grid,FALSE);res$iteration<-iteration;res$branch<-spec$branch
out<-file.path(root,'04_results/checkpoints');dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(res,file.path(out,sprintf('%s_I%03d.csv',id,iteration)),row.names=FALSE)
saveRDS(fits,file.path(out,sprintf('%s_I%03d_fits.rds',id,iteration)),compress='gzip')
dir.create(file.path(root,'02_mi/combined'),recursive=TRUE,showWarnings=FALSE)
saveRDS(obj,file.path(root,'02_mi/combined',sprintf('%s_I%03d.rds',id,iteration)),compress='gzip')
cat('CHECKPOINT_POOLED',id,iteration,nchains,'\n');print(res[c('outcome','model','estimate','ci_low','ci_high','monte_carlo_error_over_total_se')])
)--------------------"

upf_source[["01_code/18_table1_cholesterol_summary.R"]] <- r"--------------------(Sys.setenv(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
root <- '/storage/home/tmu2301/nhanes分析/18_formal_cholesterol_design_reanalysis_20260907'
source(file.path(root,'NHANES/02_代码/正式五结局/23_extension_common.R'))
d <- readRDS(file.path(DIRS$derived,'analysis_frame_pre_mice.rds'))
spec <- readRDS(file.path(EXT_DIRS$derived,'extension_analysis_spec.rds'))
d <- add_extension_exposure_fields(d,spec)
stopifnot(all(na.omit(unique(d$high_cholesterol)) %in% c(0,1)))
full <- make_full_design(d)
base <- full$variables$base_domain %in% TRUE & is.finite(full$variables[[EXPOSURE_VAR]])
stopifnot(sum(base)==spec$n_base)
des <- full[base,]
res <- lapply(c('Overall',paste0('Q',1:4)),function(g){
 z <- if(g=='Overall') des else des[des$variables$upf_quartile==g,]
 x <- z$variables$high_cholesterol
 est <- survey::svymean(~high_cholesterol,z,na.rm=TRUE)
 data.frame(group=g,n_total=length(x),n_nonmissing=sum(!is.na(x)),n_missing=sum(is.na(x)),n_positive=sum(x==1,na.rm=TRUE),weighted_prevalence=unname(coef(est)),weighted_se=unname(survey::SE(est)),weighted_percent=100*unname(coef(est)),se_percentage_points=100*unname(survey::SE(est)),missing_unweighted_percent=100*mean(is.na(x)),survey_design_df=survey::degf(z),base_n=spec$n_base,q25=spec$quartile_cutpoints[1],q50=spec$quartile_cutpoints[2],q75=spec$quartile_cutpoints[3],denominator='observed high_cholesterol within base_domain; pre-MI',quartile_rule='saved weighted cuts; right=TRUE; include.lowest=TRUE',weight='analysis_weight',design='cycle x SDMVSTRA; cycle x SDMVPSU; nest=TRUE')
})
res <- do.call(rbind,res)
stopifnot(sum(res$n_total[-1])==res$n_total[1],sum(res$n_missing[-1])==res$n_missing[1])
write.csv(res,file.path(root,'04_results/Table1_high_cholesterol_observed_base_domain.csv'),row.names=FALSE)
print(res[,1:10]);print(spec$quartile_cutpoints)
)--------------------"

upf_source[["01_code/19_table1_F_and_S2_missingness.R"]] <- r"--------------------(Sys.setenv(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
root <- '/storage/home/tmu2301/nhanes分析/18_formal_cholesterol_design_reanalysis_20260907'
source(file.path(root,'NHANES/02_代码/正式五结局/23_extension_common.R'))
d <- readRDS(file.path(DIRS$derived,'analysis_frame_pre_mice.rds'))
spec <- readRDS(file.path(EXT_DIRS$derived,'extension_analysis_spec.rds'))
d <- add_extension_exposure_fields(d,spec);full<-make_full_design(d)
base<-full$variables$base_domain %in% TRUE & is.finite(full$variables[[EXPOSURE_VAR]])
des<-full[base,];observed<-des[!is.na(des$variables$high_cholesterol),]
f<-survey::svychisq(~high_cholesterol+upf_quartile,observed,statistic='F')
res<-data.frame(variable='high_cholesterol',test='Rao-Scott adjusted F',statistic=unname(f$statistic),df_num=unname(f$parameter[1]),df_den=unname(f$parameter[2]),p_value=f$p.value,base_n=nrow(des$variables),observed_n=nrow(observed$variables),missing_n=sum(is.na(des$variables$high_cholesterol)),method=f$method)
write.csv(res,file.path(root,'04_results/Table1_high_cholesterol_RaoScott_F.csv'),row.names=FALSE);print(res)
manifest<-jsonlite::read_json(file.path(root,'03_qc/mi_job_manifest.json'),simplifyVector=TRUE)
rows<-lapply(seq_len(nrow(manifest)),function(i){
 id<-manifest$id[i];obj<-readRDS(file.path(root,'02_mi/combined',paste0(id,'_I020.rds')));mi<-obj$mids
 stopifnot(mi$m==50,mi$iteration==20)
 nm<-mi$nmis;actual<-colSums(is.na(mi$data));stopifnot(identical(as.integer(nm),as.integer(actual[names(nm)])))
 data.frame(id=id,branch=manifest$branch[i],outcome=manifest$outcome[i],variable=names(nm),frame_n=nrow(mi$data),missing_n=as.integer(nm),missing_percent=100*as.integer(nm)/nrow(mi$data),imputation_method=unname(mi$method[names(nm)]))
})
all<-do.call(rbind,rows);write.csv(all,file.path(root,'04_results/S2_PanelA_missingness_all_frames.csv'),row.names=FALSE)
groups<-split(all,interaction(all$branch,all$variable,drop=TRUE))
ranges<-do.call(rbind,lapply(groups,function(z){data.frame(branch=z$branch[1],variable=z$variable[1],frames=nrow(z),frame_n_min=min(z$frame_n),frame_n_max=max(z$frame_n),missing_n_min=min(z$missing_n),missing_n_max=max(z$missing_n),missing_percent_min=min(z$missing_percent),missing_percent_max=max(z$missing_percent),any_imputed=any(nzchar(z$imputation_method)),definition='min/max across outcome-specific pre-imputation frames; percentage computed in each frame before extrema')}))
write.csv(ranges,file.path(root,'04_results/S2_PanelA_missingness_ranges.csv'),row.names=FALSE);print(ranges[ranges$variable=='high_cholesterol',])
)--------------------"

upf_source[["helpers/extend_rcs.R"]] <- r"--------------------(root <- '/storage/home/tmu2301/nhanes分析/18_formal_cholesterol_design_reanalysis_20260907'
out <- file.path(root,'RCS_low_extension_20260911'); dir.create(out,showWarnings=FALSE)
src <- file.path(root,'RCS/02_绘图数据_v3重算/统计重算源')
models <- readRDS(file.path(src,'NHANES_RCS_M3_v3_pooled_models.rds'))
old <- read.csv(file.path(src,'NHANES_RCS_M3_v3_curves.csv'))
basis <- function(x,k) {x<-x/10;k<-k/10; c(x,(pmax(x-k[1],0)^3-pmax(x-k[2],0)^3*(k[3]-k[1])/(k[3]-k[2])+pmax(x-k[3],0)^3*(k[2]-k[1])/(k[3]-k[2]))/(k[3]-k[1])^2)}
rows<-list(); audit<-list()
for (oc in names(models)) {
 d<-old[old$outcome==oc,]; p<-models[[oc]]; terms<-c('rcs_linear_v3','rcs_nonlin1_v3')
 k<-as.numeric(d[1,c('knot_10','knot_50','knot_90')]); df<-min(p$df[match(terms,names(p$coefficients))],na.rm=TRUE)
 pred<-function(x) {L<-basis(x,k)-basis(10,k); b<-sum(L*p$coefficients[terms]); se<-sqrt(as.numeric(t(L)%*%p$variance[terms,terms]%*%L)); v<-c(b,b-qt(.975,df)*se,b+qt(.975,df)*se); if(oc %in% c('CKD','DKD','kidney_stones'))v<-exp(v);if(oc=='UACR')v<-100*expm1(v);v}
 recovered<-t(vapply(d$exposure_value,pred,numeric(3))); err<-max(abs(recovered-as.matrix(d[,c('estimate','conf_low','conf_high')])));stopifnot(err<1e-8)
 grid<-seq(0,min(d$exposure_value),length.out=17)[-17]; z<-t(vapply(grid,pred,numeric(3)))
 rows[[oc]]<-data.frame(x=grid,estimate=z[,1],conf_low=z[,2],conf_high=z[,3],survey='NHANES',outcome=oc)
 audit[[oc]]<-data.frame(outcome=oc,old_points=nrow(d),max_abs_error=err,new_points=length(grid))
}
sdir<-file.path(root,'SUA/final_alignment/RCS_aligned'); fits<-readRDS(file.path(sdir,'NHANES_serum_uric_acid_RCS_fit_summaries.rds')); stopifnot(length(fits)==50)
d<-read.csv(file.path(sdir,'NHANES_serum_uric_acid_RCS_curve.csv')); k<-as.numeric(d[1,c('knot_10','knot_50','knot_90')])
pred<-function(x){L<-setNames(rep(0,length(fits[[1]]$coefficients)),names(fits[[1]]$coefficients));L[c('rcs_linear','rcs_nonlin1')]<-basis(x,k)-basis(10,k);p<-mitools::MIcombine(lapply(fits,function(s)sum(L*s$coefficients)),lapply(fits,function(s)as.numeric(t(L)%*%s$variance%*%L)),df.complete=min(vapply(fits,`[[`,numeric(1),'df_complete'))); b<-as.numeric(p$coefficients);se<-sqrt(as.numeric(p$variance)); c(b,b-qt(.975,p$df)*se,b+qt(.975,p$df)*se)}
recovered<-t(vapply(d$upf_grams_pct,pred,numeric(3)));err<-max(abs(recovered-as.matrix(d[,c('estimate','ci_low','ci_high')])));stopifnot(err<1e-8)
grid<-seq(0,min(d$upf_grams_pct),length.out=17)[-17];z<-t(vapply(grid,pred,numeric(3)))*59.48
rows[['serum_uric_acid']]<-data.frame(x=grid,estimate=z[,1],conf_low=z[,2],conf_high=z[,3],survey='NHANES',outcome='serum_uric_acid')
audit[['serum_uric_acid']]<-data.frame(outcome='serum_uric_acid',old_points=nrow(d),max_abs_error=err,new_points=length(grid))
write.csv(do.call(rbind,rows),file.path(out,'extension_display_values.csv'),row.names=FALSE)
write.csv(do.call(rbind,audit),file.path(out,'reproduction_audit.csv'),row.names=FALSE)
print(do.call(rbind,audit))
)--------------------"

upf_source[["KN/source_snapshot/00_config.R"]] <- r"--------------------(# KNHANES 2019-2024 external-validation configuration.

KN_PROJECT_ID <- "KNHANES_2019_2024_REANALYSIS_20260820"
KN_VERSION <- "2026-08-26.v1.1-four-outcome-aligned"
KN_ROOT_DEFAULT <- "/storage/home/tmu2301/nhanes分析/05_KNHANES_外部验证_2019_2024_20260820"
KN_ROOT <- Sys.getenv("KNHANES_ROOT", unset = KN_ROOT_DEFAULT)

DIRS <- list(
  plan = file.path(KN_ROOT, "00_计划文档"),
  raw = file.path(KN_ROOT, "01_官方原始数据_只读"),
  config = file.path(KN_ROOT, "02_数据配置"),
  derived = file.path(KN_ROOT, "03_派生数据"),
  code = file.path(KN_ROOT, "04_代码"),
  results = file.path(KN_ROOT, "05_结果"),
  tables = file.path(KN_ROOT, "06_表格"),
  figures = file.path(KN_ROOT, "07_图表"),
  cross_database = file.path(KN_ROOT, "08_跨数据库比较"),
  logs = file.path(KN_ROOT, "09_运行记录"),
  qc = file.path(KN_ROOT, "10_质量控制"),
  archive = file.path(KN_ROOT, "99_归档")
)

EXTRACT_ROOT_DEFAULT <- file.path(DIRS$raw, "解压SAS_只读")
EXTRACT_ROOT <- Sys.getenv("KNHANES_EXTRACT_ROOT", unset = EXTRACT_ROOT_DEFAULT)
CODEBOOK_PATHS <- c(
  phase8 = file.path(
    DIRS$plan, "官方文档归档", "2019_2021_官网当前版",
    "제8기(2019-2021) 영양조사 코드자료집.xlsx"
  ),
  year2022 = file.path(
    DIRS$plan, "官方文档归档", "2022_2024_官网当前版",
    "제9기 1차년도(2022) 영양조사 코드자료집.xlsx"
  ),
  year2023 = file.path(
    DIRS$plan, "官方文档归档", "2022_2024_官网当前版",
    "제9기 2차년도(2023) 영양조사 코드자료집.xlsx"
  ),
  year2024 = file.path(
    DIRS$plan, "官方文档归档", "2022_2024_官网当前版",
    "제9기 3차년도(2024) 영양조사 코드자료집.xlsx"
  )
)
CODEBOOK_PATH <- unname(CODEBOOK_PATHS[["phase8"]])
KOREAN_UPF_WORKING_GROUP_PATH <- file.path(
  DIRS$plan, "外部方法学资料",
  "Korean_UPF_Working_Group_public_classification_20260818.xlsx"
)
NHANES_PROJECT_ROOT <- file.path(
  dirname(KN_ROOT), "06_NHANES_UPF七结局_第一阶段探索_20260824"
)
NHANES_RESULT_PATH <- file.path(
  NHANES_PROJECT_ROOT, "04_结果", "正式五结局_v1", "phase1_models_all_outcomes.csv"
)
NHANES_OUTCOME_KEY <- c(
  eGFR = "eGFR", log_uacr = "UACR", CKD = "CKD", DKD = "DKD"
)
USER_GUIDE_PATHS <- c(
  phase8 = file.path(
    DIRS$plan, "官方文档归档", "2019_2021_官网当前版",
    "국민건강영양조사 제8기(2019-2021) 원시자료 이용지침서.pdf"
  ),
  phase9 = file.path(
    DIRS$plan, "官方文档归档", "2022_2024_官网当前版",
    "국민건강영양조사 제9기(2022-2024) 원시자료 이용지침서.pdf"
  )
)
USER_GUIDE_PATH <- unname(USER_GUIDE_PATHS)

YEARS <- 2019:2024
YEAR_FILES <- data.frame(
  year = YEARS,
  all = file.path(EXTRACT_ROOT, sprintf("hn%d_all.sas7bdat", 19:24)),
  rc24 = file.path(EXTRACT_ROOT, sprintf("hn%d_24rc.sas7bdat", 19:24)),
  stringsAsFactors = FALSE
)

options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

KN_SEED <- 20260820L
MI_M <- suppressWarnings(as.integer(Sys.getenv("KNHANES_MI_M", unset = "2")))
MI_MAXIT <- suppressWarnings(as.integer(Sys.getenv("KNHANES_MI_MAXIT", unset = "2")))
if (!is.finite(MI_M) || MI_M < 1L) stop("KNHANES_MI_M must be a positive integer")
if (!is.finite(MI_MAXIT) || MI_MAXIT < 1L) stop("KNHANES_MI_MAXIT must be positive")
ANALYSIS_MODE <- if (identical(MI_M, 2L) && identical(MI_MAXIT, 2L)) {
  "smoke"
} else if (identical(MI_MAXIT, 20L) && (identical(MI_M, 50L) || identical(MI_M, 100L))) {
  "formal"
} else {
  "unsupported"
}
ANALYSIS_SCOPE <- if (identical(ANALYSIS_MODE, "smoke")) {
  "SOFTWARE_SMOKE_ONLY_PROVISIONAL_NOVA_NO_INFERENCE"
} else if (identical(ANALYSIS_MODE, "formal")) {
  "FORMAL_ANALYSIS"
} else {
  "UNSUPPORTED_EXECUTION_SETTINGS"
}
ANALYSIS_SCOPE_DIR <- if (identical(ANALYSIS_MODE, "formal")) {
  "FORMAL"
} else if (identical(ANALYSIS_MODE, "smoke")) {
  "SMOKE_ONLY"
} else {
  "UNSUPPORTED"
}
ANALYSIS_DIRS <- list(
  derived = file.path(DIRS$derived, ANALYSIS_SCOPE_DIR),
  results = file.path(DIRS$results, ANALYSIS_SCOPE_DIR),
  tables = file.path(DIRS$tables, ANALYSIS_SCOPE_DIR),
  figures = file.path(DIRS$figures, ANALYSIS_SCOPE_DIR),
  cross_database = file.path(DIRS$cross_database, ANALYSIS_SCOPE_DIR),
  logs = file.path(DIRS$logs, ANALYSIS_SCOPE_DIR),
  qc = file.path(DIRS$qc, ANALYSIS_SCOPE_DIR)
)
MAX_PARALLEL_WORKERS <- 14L
requested_workers <- suppressWarnings(as.integer(
  Sys.getenv("KNHANES_N_CORES", unset = "14")
))
if (!is.finite(requested_workers) || requested_workers < 1L) requested_workers <- 1L
available_workers <- tryCatch(
  as.integer(parallelly::availableCores()),
  error = function(e) as.integer(parallel::detectCores(logical = TRUE))
)
if (!is.finite(available_workers) || available_workers < 1L) available_workers <- 1L
N_WORKERS <- min(requested_workers, MAX_PARALLEL_WORKERS, available_workers)

UPF_UNIT_PP <- 10
RCS_KNOT_PROBS <- c(0.10, 0.50, 0.90)
RCS_DISPLAY_PROBS <- c(0.01, 0.99)
QUARTILE_PROBS <- c(0.25, 0.50, 0.75)
ENERGY_SENSITIVITY_RANGE <- c(500, 5000)
UNCLASSIFIED_GRAM_SHARE_LIMIT <- 0
NOVA_RECONCILIATION_ABS_TOL_G <- 0.10
NOVA_RECONCILIATION_REL_TOL <- 1e-4
NOVA_RULESET_VERSION <- "KN_KOREAN_NOVA_2019_2024_20260820_v1_context_identity"
LAB_LIMIT_RULE_VERSION <- "KNHANES_VIII_IX_YEAR_SPECIFIC_QUALIFIER_LOD_SQRT2_v1"
LAB_LIMIT_RULES <- data.frame(
  analyte = rep(c("urine_creatinine", "urine_albumin"), each = 6L),
  year = rep(YEARS, times = 2L),
  default_lower = c(
    rep(1.0, 3L), rep(NA_real_, 3L),
    rep(1.0, 3L), rep(3.0, 3L)
  ),
  default_upper = c(
    rep(NA_real_, 6L),
    rep(5000.0, 3L), rep(4400.0, 3L)
  ),
  qualifier_expected = c(
    c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
    rep(TRUE, 6L)
  ),
  stringsAsFactors = FALSE
)
HIGH_CHOLESTEROL_THRESHOLD <- 240
HYPERTENSION_BP_THRESHOLD <- c(sbp = 140, dbp = 90)
DIABETES_HBA1C_THRESHOLD <- 6.5
DIABETES_GLUCOSE_THRESHOLD <- 126
MIN_BINARY_EVENTS <- 5L
MIN_BINARY_NONEVENTS <- 5L

OUTCOMES <- list(
  eGFR = list(variable = "egfr_2021", family = "gaussian", domain = "egfr_observed"),
  log_uacr = list(variable = "log_uacr", family = "gaussian", domain = "uacr_observed"),
  CKD = list(variable = "ckd", family = "binomial", domain = "kidney_observed"),
  DKD = list(variable = "dkd", family = "binomial", domain = "dkd_domain")
)

MODEL_COVARIATES <- list(
  M1 = character(),
  M2 = c("age", "sex", "education", "survey_year", "energy_kcal"),
  M3 = c("age", "sex", "education", "survey_year", "energy_kcal",
         "smoking", "alcohol", "physical_activity", "bmi"),
  M4 = c("age", "sex", "education", "survey_year", "energy_kcal",
         "smoking", "alcohol", "physical_activity", "bmi",
         "hypertension", "diabetes", "high_cholesterol")
)
MODEL_COVARIATES_DKD <- within(MODEL_COVARIATES, M4 <- setdiff(M4, "diabetes"))

IMPUTED_COVARIATES <- c(
  "education", "smoking", "alcohol", "physical_activity",
  "bmi", "hypertension", "diabetes", "high_cholesterol"
)

SUBGROUPS <- list(
  age_group = c("20-44", "45-64", "65+"),
  sex = c("Male", "Female"),
  smoking = c("Never", "Former", "Current"),
  hypertension = c("No", "Yes"),
  diabetes = c("No", "Yes")
)

REQUIRED_PACKAGES <- c(
  "haven", "survey", "mice", "mitools", "dplyr", "tidyr", "purrr", "stringr",
  "readr", "data.table", "splines", "ggplot2", "broom", "openxlsx", "digest",
  "withr", "future", "future.apply", "parallelly", "mitml", "jsonlite", "testthat", "lintr", "renv", "patchwork",
  "cowplot", "forestploter", "gtsummary", "labelled"
)

SMOKE_USER_AUTH_MARKER <- file.path(DIRS$plan, "PREANALYSIS_SMOKE_AUTHORIZED_BY_USER.txt")
FORMAL_USER_AUTH_MARKER <- file.path(DIRS$plan, "FORMAL_ANALYSIS_AUTHORIZED_BY_USER.txt")
CODE_REVIEW_MARKER <- file.path(DIRS$code, "CODE_REVIEW_COMPLETED.txt")
SMOKE_NOVA_GATE_MARKER <- file.path(DIRS$config, "NOVA_GATE_SMOKE_PASSED.txt")
FORMAL_NOVA_GATE_MARKER <- file.path(DIRS$config, "NOVA_GATE_FORMAL_PASSED.txt")
FORMAL_NOVA_MAP <- file.path(DIRS$config, "nova_context_adjudicated.csv")
SMOKE_NOVA_MAP <- file.path(DIRS$config, "nova_context_smoke_only.csv")
NOVA_CONTENT_REVIEW_INPUT_MANIFEST <- file.path(
  DIRS$qc, "nova_content_review_input_manifest.csv"
)
NOVA_CONTENT_REVIEW_SCOPE_MARKER <- file.path(
  DIRS$qc, "NOVA_CONTENT_REVIEW_SCOPE_FROZEN.txt"
)
OUTCOME_DEFINITION_MANIFEST <- file.path(
  DIRS$config, "outcome_definitions_v1.csv"
)
ACTIVE_NOVA_MAP <- if (identical(MI_M, 2L) && identical(MI_MAXIT, 2L)) {
  SMOKE_NOVA_MAP
} else {
  FORMAL_NOVA_MAP
}
NOVA_GATE_MARKER <- if (identical(MI_M, 2L) && identical(MI_MAXIT, 2L)) {
  SMOKE_NOVA_GATE_MARKER
} else {
  FORMAL_NOVA_GATE_MARKER
}

read_gate_values <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  kv_lines <- lines[grepl("^[A-Za-z0-9_]+=", lines)]
  keys <- sub("=.*$", "", kv_lines)
  values <- sub("^[^=]*=", "", kv_lines)
  stats::setNames(values, keys)
}

require_user_authorization <- function(stage) {
  if (identical(MI_M, 2L) && identical(MI_MAXIT, 2L)) {
    marker_path <- SMOKE_USER_AUTH_MARKER
    expected_scope <- "PREANALYSIS_AND_M2_SMOKE_ONLY"
    mode <- "smoke"
  } else if (identical(MI_MAXIT, 20L) && (identical(MI_M, 50L) || identical(MI_M, 100L))) {
    marker_path <- FORMAL_USER_AUTH_MARKER
    expected_scope <- paste0("FORMAL_M", MI_M, "_ANALYSIS")
    mode <- "formal"
    if (!identical(Sys.getenv("KNHANES_ALLOW_FORMAL", unset = "NO"), "YES")) {
      stop("Formal mode additionally requires KNHANES_ALLOW_FORMAL=YES", call. = FALSE)
    }
  } else {
    stop("Execution settings are neither authorized smoke (m=2/maxit=2) nor formal (m=50 or m=100, maxit=20)",
         call. = FALSE)
  }
  if (!file.exists(marker_path)) {
    stop("User-authorization marker is absent for ", mode, "; refusing stage: ",
         stage, call. = FALSE)
  }
  marker <- read_gate_values(marker_path)
  if (!identical(unname(marker[["status"]]), "AUTHORIZED") ||
      !identical(unname(marker[["project_id"]]), KN_PROJECT_ID) ||
      !identical(unname(marker[["scope"]]), expected_scope)) {
    stop("User-authorization marker content does not match this project/mode", call. = FALSE)
  }
  invisible(TRUE)
}

report_field <- function(lines, key, report_path) {
  clean <- trimws(lines)
  hit <- grep(paste0("^", key, "="), clean, value = TRUE)
  if (length(hit) != 1L) {
    stop("Review report must contain exactly one machine-readable ", key,
         " field: ", report_path, call. = FALSE)
  }
  value <- sub(paste0("^", key, "="), "", hit)
  if (!nzchar(value)) stop("Empty review-report field ", key, ": ", report_path)
  value
}

validate_code_review_report <- function(report_path, expected_reviewer, freeze_id,
                                        manifest_sha256) {
  assert_file <- function(path) {
    if (!file.exists(path)) stop("Review report is missing: ", path, call. = FALSE)
  }
  assert_file(report_path)
  lines <- readLines(report_path, warn = FALSE, encoding = "UTF-8")
  expected <- c(
    verdict = "PASSED", project_id = KN_PROJECT_ID,
    reviewer = expected_reviewer, freeze_id = freeze_id,
    snapshot_manifest_sha256 = manifest_sha256
  )
  actual <- vapply(names(expected), function(key) {
    report_field(lines, key, report_path)
  }, character(1))
  if (!identical(unname(actual), unname(expected))) {
    stop("Review report verdict/reviewer/project/freeze lineage is invalid: ",
         report_path, call. = FALSE)
  }
  invisible(TRUE)
}

require_code_review <- function(stage) {
  require_user_authorization(stage)
  if (!file.exists(CODE_REVIEW_MARKER)) {
    stop("Code-review marker is absent; refusing stage: ", stage, call. = FALSE)
  }
  marker <- readLines(CODE_REVIEW_MARKER, warn = FALSE, encoding = "UTF-8")
  if (!any(trimws(marker) == "status=PASSED")) {
    stop("Code-review marker is not PASSED; refusing stage: ", stage, call. = FALSE)
  }
  project_line <- grep("^project_id=", marker, value = TRUE)
  if (length(project_line) != 1L ||
      !identical(sub("^project_id=", "", project_line), KN_PROJECT_ID)) {
    stop("Code-review marker is for a different project", call. = FALSE)
  }
  manifest_line <- grep("^snapshot_manifest=", marker, value = TRUE)
  sha_line <- grep("^snapshot_manifest_sha256=", marker, value = TRUE)
  if (length(manifest_line) != 1L || length(sha_line) != 1L) {
    stop("Code-review marker lacks snapshot lineage", call. = FALSE)
  }
  manifest_path <- sub("^snapshot_manifest=", "", manifest_line)
  expected_sha <- sub("^snapshot_manifest_sha256=", "", sha_line)
  if (!file.exists(manifest_path) ||
      !identical(digest::digest(file = manifest_path, algo = "sha256"), expected_sha)) {
    stop("Code-review snapshot manifest hash mismatch", call. = FALSE)
  }
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE,
                              check.names = FALSE, fileEncoding = "UTF-8-BOM")
  if (!all(c("file", "sha256") %in% names(manifest)) ||
      anyDuplicated(manifest$file)) {
    stop("Code-review snapshot manifest schema is invalid", call. = FALSE)
  }
  if (!"scope" %in% names(manifest)) {
    stop("Code-review manifest lacks scope", call. = FALSE)
  }
  code_rows <- manifest$scope == "active_code"
  reviewed_files <- as.character(manifest$file[code_rows])
  current_files <- list.files(
    DIRS$code,
    pattern = "(\\.(R|py)$)|^(DESCRIPTION|LICENSE|renv\\.lock)$",
    recursive = FALSE, full.names = FALSE
  )
  if (!identical(sort(reviewed_files), sort(current_files))) {
    stop("Code-review manifest does not cover the complete current active code set",
         call. = FALSE)
  }
  current_paths <- file.path(DIRS$code, reviewed_files)
  if (any(!file.exists(current_paths))) stop("A reviewed active-code file is missing", call. = FALSE)
  current_sha <- vapply(current_paths, digest::digest, character(1),
                        file = TRUE, algo = "sha256")
  if (!identical(unname(current_sha), as.character(manifest$sha256[code_rows]))) {
    stop("Current active code differs from the reviewed frozen snapshot", call. = FALSE)
  }
  round_line <- grep("^review_round=", marker, value = TRUE)
  review_round <- suppressWarnings(as.integer(sub("^review_round=", "", round_line)))
  if (length(round_line) != 1L || !is.finite(review_round) || review_round < 1L) {
    stop("Code-review marker has no valid completed review round", call. = FALSE)
  }
  freeze_line <- grep("^freeze_id=", marker, value = TRUE)
  if (length(freeze_line) != 1L || !nzchar(sub("^freeze_id=", "", freeze_line))) {
    stop("Code-review marker has no valid freeze_id", call. = FALSE)
  }
  freeze_id <- sub("^freeze_id=", "", freeze_line)
  report_paths <- character()
  for (reviewer in c("A", "B")) {
    path_line <- grep(paste0("^reviewer_", reviewer, "_report="), marker, value = TRUE)
    report_sha_line <- grep(paste0("^reviewer_", reviewer, "_report_sha256="),
                            marker, value = TRUE)
    if (length(path_line) != 1L || length(report_sha_line) != 1L) {
      stop("Code-review marker lacks Reviewer ", reviewer, " report lineage", call. = FALSE)
    }
    report_path <- sub(paste0("^reviewer_", reviewer, "_report="), "", path_line)
    report_paths <- c(report_paths, report_path)
    report_sha <- sub(paste0("^reviewer_", reviewer, "_report_sha256="), "",
                      report_sha_line)
    if (!file.exists(report_path) ||
        !identical(digest::digest(file = report_path, algo = "sha256"), report_sha)) {
      stop("Reviewer ", reviewer, " report hash mismatch", call. = FALSE)
    }
    validate_code_review_report(
      report_path, paste0("Reviewer_", reviewer), freeze_id, expected_sha
    )
  }
  if (length(unique(normalizePath(report_paths, mustWork = TRUE))) != 2L) {
    stop("Reviewer A and B reports must be distinct files", call. = FALSE)
  }
  invisible(TRUE)
}

require_nova_gate <- function(stage) {
  require_code_review(stage)
  if (!file.exists(NOVA_GATE_MARKER)) {
    stop("NOVA QA gate is absent; refusing stage: ", stage, call. = FALSE)
  }
  marker <- readLines(NOVA_GATE_MARKER, warn = FALSE, encoding = "UTF-8")
  marker_values <- read_gate_values(NOVA_GATE_MARKER)
  expected_status <- if (identical(MI_M, 2L) && identical(MI_MAXIT, 2L)) {
    c("PASSED", "PASSED_FOR_SMOKE_ONLY")
  } else {
    "PASSED"
  }
  if (!(unname(marker_values[["status"]]) %in% expected_status) ||
      !identical(unname(marker_values[["project_id"]]), KN_PROJECT_ID)) {
    stop("NOVA gate is not valid for this project/execution mode: ", stage,
         call. = FALSE)
  }
  map_line <- grep("^final_map=", marker, value = TRUE)
  sha_line <- grep("^final_map_sha256=", marker, value = TRUE)
  if (length(map_line) != 1L || length(sha_line) != 1L) {
    stop("NOVA marker lacks map lineage", call. = FALSE)
  }
  map_path <- sub("^final_map=", "", map_line)
  expected_sha <- sub("^final_map_sha256=", "", sha_line)
  if (!file.exists(map_path) || !identical(digest::digest(file = map_path, algo = "sha256"),
                                            expected_sha)) {
    stop("NOVA map hash no longer matches its gate", call. = FALSE)
  }
  if (!identical(normalizePath(map_path, mustWork = TRUE),
                 normalizePath(ACTIVE_NOVA_MAP, mustWork = TRUE))) {
    stop("NOVA gate points to the wrong map for the current execution mode", call. = FALSE)
  }
  if (identical(unname(marker_values[["status"]]), "PASSED_FOR_SMOKE_ONLY")) {
    if (!identical(MI_M, 2L) || !identical(MI_MAXIT, 2L)) {
      stop("Smoke-only NOVA map cannot authorize formal analysis", call. = FALSE)
    }
    return(invisible(TRUE))
  }
  provenance_line <- grep("^content_review_provenance=", marker, value = TRUE)
  provenance_sha_line <- grep("^content_review_provenance_sha256=", marker, value = TRUE)
  if (length(provenance_line) != 1L || length(provenance_sha_line) != 1L) {
    stop("NOVA marker lacks content-review provenance lineage", call. = FALSE)
  }
  provenance_path <- sub("^content_review_provenance=", "", provenance_line)
  provenance_sha <- sub("^content_review_provenance_sha256=", "", provenance_sha_line)
  if (!file.exists(provenance_path) ||
      !identical(digest::digest(file = provenance_path, algo = "sha256"), provenance_sha)) {
    stop("NOVA content-review provenance hash mismatch", call. = FALSE)
  }
  invisible(TRUE)
}
)--------------------"

upf_source[["KN/source_snapshot/00_functions.R"]] <- r"--------------------(.this_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "KN/source_snapshot/00_functions.R"), error = function(e) NULL)
if (is.null(.this_file)) .this_file <- "00_functions.R"
.code_dir <- dirname(normalizePath(.this_file, mustWork = FALSE))
source(file.path(.code_dir, "00_config.R"), local = FALSE)

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

assert_columns <- function(x, cols, label = deparse(substitute(x))) {
  missing <- setdiff(cols, names(x))
  if (length(missing)) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

assert_files <- function(paths, label = "input files") {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Missing ", label, ": ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

sha256_file <- function(path) digest::digest(file = path, algo = "sha256")

configure_worker_runtime <- function() {
  options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)
  Sys.setenv(
    OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1"
  )
  invisible(list(
    lonely_psu = getOption("survey.lonely.psu"),
    adjust_domain_lonely = getOption("survey.adjust.domain.lonely"),
    library_paths = .libPaths()
  ))
}

read_key_value_marker <- function(path) {
  assert_files(path, "marker")
  lines <- trimws(readLines(path, warn = FALSE, encoding = "UTF-8"))
  lines <- lines[nzchar(lines) & grepl("=", lines, fixed = TRUE)]
  keys <- sub("=.*$", "", lines)
  values <- sub("^[^=]*=", "", lines)
  stats::setNames(values, keys)
}

marker_run_id <- function(path) {
  assert_files(path, "latest-run marker")
  lines <- trimws(readLines(path, warn = FALSE, encoding = "UTF-8"))
  run_line <- grep("^run_id=", lines, value = TRUE)
  if (length(run_line) == 1L) return(sub("^run_id=", "", run_line))
  simple <- lines[nzchar(lines) & !grepl("=", lines, fixed = TRUE)]
  if (length(simple) == 1L) return(simple[[1]])
  stop("Latest-run marker has no unique run_id: ", path, call. = FALSE)
}

commit_latest_marker <- function(run_id, latest_path, required_paths,
                                 manifest_path = NULL, parent_run_ids = character(),
                                 extra_lines = character()) {
  assert_files(required_paths, "stage outputs before marker commit")
  if (is.null(manifest_path)) {
    hits <- required_paths[basename(required_paths) == "output_manifest.csv"]
    if (length(hits) == 1L) manifest_path <- hits[[1]]
  }
  if (!is.null(manifest_path)) assert_files(manifest_path, "stage output manifest")
    lines <- c(
      "status=PASSED", paste0("run_id=", run_id),
      paste0("created_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3%z"))
    )
    if (exists("ANALYSIS_SCOPE", inherits = TRUE)) {
      lines <- c(lines, paste0("analysis_scope=", get("ANALYSIS_SCOPE", inherits = TRUE)))
    }
  if (!is.null(manifest_path)) {
    lines <- c(lines, paste0("output_manifest=", manifest_path),
               paste0("output_manifest_sha256=", sha256_file(manifest_path)))
  }
    if (length(parent_run_ids)) {
      lines <- c(lines, paste0("parent_run_ids=", paste(parent_run_ids, collapse = "|")))
    }
    if (length(extra_lines)) {
      if (any(!grepl("^[A-Za-z0-9_]+=.+$", extra_lines)) ||
          anyDuplicated(sub("=.*$", "", extra_lines))) {
        stop("Latest-marker extra lines must be unique non-empty key=value pairs")
      }
      lines <- c(lines, extra_lines)
    }
  write_lines_atomic(lines, latest_path)
  invisible(latest_path)
}

verify_manifest_entry <- function(path, manifest_path) {
  assert_files(c(path, manifest_path), "manifest verification inputs")
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE,
                              check.names = FALSE, fileEncoding = "UTF-8-BOM")
  assert_columns(manifest, c("path", "sha256"), "output manifest")
  hit <- which(as.character(manifest$path) == path)
  if (length(hit) != 1L || !identical(as.character(manifest$sha256[[hit]]),
                                      sha256_file(path))) {
    stop("Artifact hash is absent from or differs from manifest: ", path)
  }
  invisible(TRUE)
}

timestamp_id <- function(prefix) {
  now <- Sys.time()
  millis <- sprintf("%03d", as.integer((as.numeric(now) %% 1) * 1000))
  paste0(prefix, "_", format(now, "%Y%m%dT%H%M%S"), millis, "_P", Sys.getpid())
}

write_csv_atomic <- function(x, path, na = "") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  readr::write_excel_csv(x, tmp, na = na)
  if (file.exists(path)) file.remove(path)
  if (!file.rename(tmp, path)) stop("Atomic CSV write failed: ", path, call. = FALSE)
  invisible(path)
}

write_lines_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  writeLines(enc2utf8(x), tmp, useBytes = TRUE)
  if (file.exists(path)) file.remove(path)
  if (!file.rename(tmp, path)) stop("Atomic text write failed: ", path, call. = FALSE)
  invisible(path)
}

read_sas_selected <- function(path, columns = NULL) {
  if (!file.exists(path)) stop("Missing SAS file: ", path, call. = FALSE)
  if (is.null(columns)) {
    haven::read_sas(path, encoding = "CP949")
  } else {
    haven::read_sas(path, col_select = tidyselect::any_of(columns), encoding = "CP949")
  }
}

as_numeric_clean <- function(x, invalid = numeric()) {
  z <- suppressWarnings(as.numeric(x))
  z[z %in% invalid] <- NA_real_
  z
}

lab_qualifier_text <- function(x) {
  out <- if (inherits(x, "haven_labelled") || inherits(x, "labelled")) {
    as.character(haven::as_factor(x, levels = "labels"))
  } else {
    as.character(x)
  }
  out[is.na(x)] <- ""
  trimws(out)
}

apply_lab_limits <- function(value, qualifier, lower = NA_real_, upper = NA_real_,
                             lower_divisor = sqrt(2)) {
  raw <- suppressWarnings(as.numeric(value))
  q <- lab_qualifier_text(qualifier)
  below <- nzchar(q) & grepl("^\\s*<", q)
  above <- nzchar(q) & grepl("^\\s*>", q)
  if (any(below) && !is.finite(lower)) stop("Lower qualifier found without LOD")
  if (any(above) && !is.finite(upper)) stop("Upper qualifier found without ULOQ")
  primary <- raw
  primary[below] <- lower / lower_divisor
  primary[above] <- upper
  alternative <- raw
  alternative[below] <- lower / 2
  alternative[above] <- upper
  data.frame(
    value_primary = primary,
    value_alternative = alternative,
    qualifier_text = q,
    below_lod = below,
    above_uloq = above,
    stringsAsFactors = FALSE
  )
}

qualifier_numeric_boundary <- function(x) {
  q <- lab_qualifier_text(x)
  match <- regexpr("[-+]?[0-9]*\\.?[0-9]+", q, perl = TRUE)
  out <- rep(NA_real_, length(q))
  hit <- match > 0L
  out[hit] <- suppressWarnings(as.numeric(regmatches(q, match)))
  out
}

apply_lab_limits_from_qualifier <- function(
    value, qualifier, default_lower = NA_real_, default_upper = NA_real_,
    lower_divisor = sqrt(2)) {
  raw <- suppressWarnings(as.numeric(value))
  q <- lab_qualifier_text(qualifier)
  below <- nzchar(q) & grepl("^\\s*<", q)
  above <- nzchar(q) & grepl("^\\s*>", q)
  unknown <- nzchar(q) & !below & !above
  if (any(unknown)) {
    stop("Unsupported laboratory qualifier text: ",
         paste(unique(q[unknown]), collapse = " | "), call. = FALSE)
  }
  parsed <- qualifier_numeric_boundary(q)
  lower <- rep_len(as.numeric(default_lower), length(raw))
  upper <- rep_len(as.numeric(default_upper), length(raw))
  lower_boundary <- ifelse(is.finite(parsed), parsed, lower)
  upper_boundary <- ifelse(is.finite(parsed), parsed, upper)
  if (any(below & !is.finite(lower_boundary))) {
    stop("Lower qualifier found without a numeric/default detection limit", call. = FALSE)
  }
  if (any(above & !is.finite(upper_boundary))) {
    stop("Upper qualifier found without a numeric/default reporting limit", call. = FALSE)
  }
  primary <- raw
  primary[below] <- lower_boundary[below] / lower_divisor
  primary[above] <- upper_boundary[above]
  alternative <- raw
  alternative[below] <- lower_boundary[below] / 2
  alternative[above] <- upper_boundary[above]
  data.frame(
    value_primary = primary,
    value_alternative = alternative,
    qualifier_text = q,
    qualifier_boundary = parsed,
    boundary_used = ifelse(below, lower_boundary, ifelse(above, upper_boundary, NA_real_)),
    below_lod = below,
    above_uloq = above,
    stringsAsFactors = FALSE
  )
}

label_lookup <- function(x) {
  labs <- attr(x, "labels", exact = TRUE)
  if (is.null(labs) || !length(labs)) {
    return(data.frame(raw_value = character(), value_label = character(),
                      stringsAsFactors = FALSE))
  }
  data.frame(
    raw_value = as.character(unname(labs)),
    value_label = names(labs),
    stringsAsFactors = FALSE
  )
}

calc_egfr_2021 <- function(creatinine_mg_dl, age_years, female) {
  invalid <- is.na(creatinine_mg_dl) | creatinine_mg_dl <= 0 |
    is.na(age_years) | age_years < 18 | is.na(female)
  k <- ifelse(female, 0.7, 0.9)
  alpha <- ifelse(female, -0.241, -0.302)
  out <- 142 * pmin(creatinine_mg_dl / k, 1)^alpha *
    pmax(creatinine_mg_dl / k, 1)^(-1.200) *
    0.9938^age_years * ifelse(female, 1.012, 1)
  out[invalid] <- NA_real_
  out
}

derive_ckd_three_state <- function(egfr, acr) {
  positive <- (!is.na(egfr) & egfr < 60) | (!is.na(acr) & acr >= 30)
  negative <- !is.na(egfr) & egfr >= 60 & !is.na(acr) & acr < 30
  ifelse(positive, 1L, ifelse(negative, 0L, NA_integer_))
}

yes_no_na <- function(x) {
  x <- as.integer(x)
  ifelse(x %in% 1L, 1L, ifelse(x %in% c(0L, 8L), 0L, NA_integer_))
}

medication_yes_no_na <- function(x) {
  x <- as.integer(x)
  ifelse(x %in% 1:4, 1L, ifelse(x %in% c(5L, 8L), 0L, NA_integer_))
}

or_composite <- function(...) {
  parts <- list(...)
  n <- max(vapply(parts, length, integer(1)))
  if (n == 0L) return(integer())
  parts <- lapply(parts, function(p) {
    p <- as.integer(p)
    if (length(p) == 1L && n > 1L) rep(p, n) else p
  })
  any_true <- rep(FALSE, n)
  all_known_false <- rep(TRUE, n)
  for (p in parts) {
    any_true <- any_true | (p %in% 1L)
    all_known_false <- all_known_false & (p %in% 0L)
  }
  out <- rep(NA_integer_, n)
  out[any_true] <- 1L
  out[all_known_false & !any_true] <- 0L
  out
}

derive_hypertension_composite <- function(sbp, dbp, med) {
  measured <- ifelse(
    is.na(sbp) & is.na(dbp), NA_integer_,
    ifelse((!is.na(sbp) & sbp >= HYPERTENSION_BP_THRESHOLD[["sbp"]]) |
             (!is.na(dbp) & dbp >= HYPERTENSION_BP_THRESHOLD[["dbp"]]), 1L, 0L)
  )
  or_composite(measured, medication_yes_no_na(med))
}

derive_diabetes_composite <- function(dx, insulin, oral, glu, hba1c) {
  glucose <- ifelse(is.na(glu), NA_integer_,
                    ifelse(glu >= DIABETES_GLUCOSE_THRESHOLD, 1L, 0L))
  hba1c_lab <- ifelse(is.na(hba1c), NA_integer_,
                      ifelse(hba1c >= DIABETES_HBA1C_THRESHOLD, 1L, 0L))
  or_composite(yes_no_na(dx), yes_no_na(insulin), yes_no_na(oral), glucose, hba1c_lab)
}

derive_high_cholesterol_composite <- function(total_chol, med,
                                              threshold = HIGH_CHOLESTEROL_THRESHOLD) {
  measured <- ifelse(is.na(total_chol), NA_integer_,
                     ifelse(total_chol >= threshold, 1L, 0L))
  or_composite(measured, medication_yes_no_na(med))
}

validate_outcome_definition_manifest <- function(path = OUTCOME_DEFINITION_MANIFEST) {
  assert_files(path, "outcome-definition manifest")
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                       fileEncoding = "UTF-8-BOM")
  required <- c("outcome_id", "label", "type", "source_variables",
                "formula_or_threshold", "missing_logic", "definition_version",
                "reference_url", "reference_note")
  assert_columns(x, required, "outcome-definition manifest")
  expected_outcomes <- names(OUTCOMES)
  if (nrow(x) != length(expected_outcomes) || anyDuplicated(x$outcome_id) ||
      !setequal(x$outcome_id, expected_outcomes) ||
      !all(x$definition_version == "OUTCOME_DEFINITIONS_20260820_v1") ||
      any(!nzchar(x$reference_url)) || any(!nzchar(x$formula_or_threshold)) ||
      !grepl("uACR>0",
             x$formula_or_threshold[x$outcome_id == "log_uacr"], fixed = TRUE) ||
      !grepl("eGFR<60 or uACR>=30",
             x$formula_or_threshold[x$outcome_id == "CKD"], fixed = TRUE)) {
    stop("Outcome-definition manifest does not match the frozen code constants",
         call. = FALSE)
  }
  invisible(x)
}

weighted_quantile <- function(x, w, probs = QUARTILE_PROBS) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(rep(NA_real_, length(probs)))
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  cw <- cumsum(w[ord]) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1))
}

safe_factor <- function(x, levels, labels) {
  factor(ifelse(x %in% levels, x, NA), levels = levels, labels = labels)
}
)--------------------"

upf_source[["KN/source_snapshot/00_dependencies.R"]] <- r"--------------------(.this_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "KN/source_snapshot/00_dependencies.R"), error = function(e) NULL)
if (is.null(.this_file)) .this_file <- "00_dependencies.R"
.code_dir <- dirname(normalizePath(.this_file, mustWork = FALSE))
source(file.path(.code_dir, "00_config.R"), local = FALSE)

r_minor <- strsplit(R.version$minor, ".", fixed = TRUE)[[1]][[1]]
project_library <- file.path(DIRS$code, "R_library", paste(R.version$major, r_minor,
                                                             sep = "."))
if (dir.exists(project_library)) .libPaths(unique(c(project_library, .libPaths())))

REQUIRED_PACKAGES <- setdiff(REQUIRED_PACKAGES,c("lintr","renv"))
missing_packages <- REQUIRED_PACKAGES[
  !vapply(REQUIRED_PACKAGES, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}

suppressPackageStartupMessages({
  library(haven)
  library(survey)
  library(mice)
  library(mitools)
  library(mitml)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(readr)
  library(data.table)
  library(splines)
  library(ggplot2)
  library(broom)
  library(openxlsx)
  library(digest)
  library(withr)
  library(future.apply)
  library(parallelly)
  library(jsonlite)
  library(testthat)
  library(patchwork)
  library(cowplot)
  library(forestploter)
  library(gtsummary)
  library(labelled)
})

loaded_ok <- vapply(
  REQUIRED_PACKAGES,
  function(pkg) paste0("package:", pkg) %in% search() || requireNamespace(pkg, quietly = TRUE),
  logical(1)
)
if (!all(loaded_ok)) {
  stop("Package load verification failed: ",
       paste(REQUIRED_PACKAGES[!loaded_ok], collapse = ", "), call. = FALSE)
}

invisible(data.frame(
  package = REQUIRED_PACKAGES,
  version = vapply(REQUIRED_PACKAGES, function(pkg) as.character(packageVersion(pkg)),
                   character(1)),
  loaded_ok = loaded_ok,
  stringsAsFactors = FALSE
))
)--------------------"

upf_source[["KN/source_snapshot/00_model_functions.R"]] <- r"--------------------(# Shared complex-survey and Rubin-pooling helpers.

.this_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "KN/source_snapshot/00_model_functions.R"), error = function(e) NULL)
if (is.null(.this_file)) .this_file <- "00_model_functions.R"
.code_dir <- dirname(normalizePath(.this_file, mustWork = FALSE))
source(file.path(.code_dir, "00_dependencies.R"), local = FALSE)
source(file.path(.code_dir, "00_functions.R"), local = FALSE)

latest_run_file <- function(marker_name, relative_file) {
  marker <- file.path(ANALYSIS_DIRS$derived, marker_name)
  assert_files(marker, marker_name)
  run_id <- marker_run_id(marker)
  path <- file.path(ANALYSIS_DIRS$derived, run_id, relative_file)
  assert_files(path, relative_file)
  kv <- read_key_value_marker(marker)
  manifest_path <- kv[["output_manifest"]]
  if (!identical(unname(kv[["status"]]), "PASSED") || is.null(manifest_path) ||
      is.null(kv[["output_manifest_sha256"]]) || !file.exists(manifest_path) ||
      !identical(sha256_file(manifest_path), unname(kv[["output_manifest_sha256"]]))) {
    stop("Derived-stage latest marker provenance is invalid: ", marker_name)
  }
  verify_manifest_entry(path, manifest_path)
  path
}

verify_mids_frame_lineage <- function(imp, frame_path) {
  frame_sha <- sha256_file(frame_path)
  stored_sha <- attr(imp, "knhanes_frame_sha256", exact = TRUE)
  if (is.null(stored_sha) || !identical(stored_sha, frame_sha)) {
    stop("MICE object is not bound to the current analysis frame hash")
  }
  stored_scope <- attr(imp, "knhanes_analysis_scope", exact = TRUE)
  if (is.null(stored_scope) || !identical(stored_scope, ANALYSIS_SCOPE)) {
    stop("MICE object belongs to a different smoke/formal analysis scope")
  }
  invisible(frame_sha)
}

load_model_inputs <- function() {
  frame_path <- latest_run_file(
    "LATEST_ANALYSIS_FRAME_RUN_ID.txt", "analysis_frame_pre_mice.rds"
  )
  mice_path <- latest_run_file(
    "LATEST_MICE_RUN_ID.txt", "covariate_imputations.mids.rds"
  )
  imp <- readRDS(mice_path)
  frame_sha <- verify_mids_frame_lineage(imp, frame_path)
  frame <- readRDS(frame_path)
  if (!"analysis_scope" %in% names(frame) ||
      !identical(unique(as.character(frame$analysis_scope)), ANALYSIS_SCOPE)) {
    stop("Analysis frame belongs to a different smoke/formal scope")
  }
  list(frame = frame, imp = imp,
       frame_path = frame_path, mice_path = mice_path,
       frame_sha256 = frame_sha, mice_sha256 = sha256_file(mice_path))
}

derive_fixed_exposure_categories <- function(frame) {
  domain <- frame$base_domain %in% TRUE
  cuts <- weighted_quantile(
    frame$upf_grams_pct[domain], frame$analysis_weight[domain], QUARTILE_PROBS
  )
  if (any(!is.finite(cuts)) || any(diff(cuts) <= 0)) {
    stop("Invalid weighted UPF quartile cut points", call. = FALSE)
  }
  q <- cut(frame$upf_grams_pct, breaks = c(-Inf, cuts, Inf),
           labels = paste0("Q", 1:4), include.lowest = TRUE, right = TRUE)
  medians <- vapply(levels(q), function(level) {
    hit <- domain & q == level
    weighted_quantile(frame$upf_grams_pct[hit], frame$analysis_weight[hit], 0.5)
  }, numeric(1))
  frame$upf_per10 <- frame$upf_grams_pct / UPF_UNIT_PP
  frame$upf_quartile <- q
  frame$upf_q_median_per10 <- medians[as.character(q)] / UPF_UNIT_PP
  attr(frame, "upf_quartile_cuts") <- cuts
  attr(frame, "upf_quartile_medians") <- medians
  frame
}

completed_analysis_dataset <- function(frame, imp, i) {
  if (!"upf_per10" %in% names(frame)) frame <- derive_fixed_exposure_categories(frame)
  covars <- mice::complete(imp, action = i)
  if (!identical(as.character(covars$person_uid), as.character(frame$person_uid))) {
    stop("MICE/frame person order mismatch at imputation ", i, call. = FALSE)
  }
  out <- frame
  for (v in IMPUTED_COVARIATES) out[[v]] <- covars[[v]]
  out$.imp <- i
  out
}

completed_analysis_datasets <- function(frame, imp) {
  frame <- derive_fixed_exposure_categories(frame)
  lapply(seq_len(imp$m), function(i) completed_analysis_dataset(frame, imp, i))
}

make_survey_design <- function(data) {
  survey::svydesign(
    ids = ~psu_uid, strata = ~strata_uid, weights = ~analysis_weight,
    nest = TRUE, data = data
  )
}

outcome_family <- function(spec) {
  if (identical(spec$family, "binomial")) stats::quasibinomial() else stats::gaussian()
}

fit_svy_formula <- function(data, formula, spec) {
  data$.model_domain <- data[[spec$domain]] %in% TRUE
  design <- make_survey_design(data)
  domain_flag <- data$.model_domain
  design_domain <- subset(design, .model_domain)
  fit <- survey::svyglm(formula, design = design_domain,
                        family = outcome_family(spec), na.action = na.omit)
  model_complete <- stats::complete.cases(
    stats::model.frame(formula, data = data, na.action = na.pass)
  )
  actual <- domain_flag & model_complete
  cluster_key <- interaction(data$strata_uid, data$psu_uid, drop = TRUE)
  list(
    coefficients = stats::coef(fit),
    variance = stats::vcov(fit),
    residual_df = as.numeric(fit$df.residual),
    n = sum(actual),
    events = if (identical(spec$family, "binomial")) {
      sum(data[[spec$variable]][actual] == 1, na.rm = TRUE)
    } else NA_integer_,
    design_df = survey::degf(design_domain),
    n_psu = length(unique(cluster_key[actual])),
    n_strata = length(unique(data$strata_uid[actual]))
  )
}

validate_compact_fit <- function(x, family, label = "survey fit") {
  if (is.null(x)) stop(label, " is NULL", call. = FALSE)
  if (!length(x$coefficients) || any(!is.finite(x$coefficients)) ||
      any(!is.finite(x$variance)) || any(diag(x$variance) < 0) ||
      !is.finite(x$residual_df) || x$residual_df <= 0 || x$n <= 0) {
    stop(label, " has non-finite/rank-deficient estimates or invalid df", call. = FALSE)
  }
  if (identical(family, "binomial") &&
      (x$events < MIN_BINARY_EVENTS || (x$n - x$events) < MIN_BINARY_NONEVENTS)) {
    stop(label, " has too few events or non-events", call. = FALSE)
  }
  invisible(TRUE)
}

pool_svy_fits <- function(fit_objects, outcome, model, exposure_form, family) {
  if (!length(fit_objects) || any(vapply(fit_objects, is.null, logical(1)))) {
    return(data.frame(
      outcome = outcome, model = model, exposure_form = exposure_form,
      term = NA_character_, estimate = NA_real_, std_error_link = NA_real_,
      conf_low = NA_real_, conf_high = NA_real_, p_value = NA_real_,
      df = NA_real_, fmi = NA_real_, complete_data_df = NA_real_,
      effect_scale = NA_character_,
      status = "FAILED_IN_AT_LEAST_ONE_IMPUTATION", stringsAsFactors = FALSE
    ))
  }
  invisible(lapply(seq_along(fit_objects), function(i) {
    validate_compact_fit(fit_objects[[i]], family, paste0(outcome, "/", model, "/", i))
  }))
  if (length(unique(vapply(fit_objects, function(x) {
    paste(names(x$coefficients), collapse = "|")
  },
                           character(1)))) != 1L) {
    stop("Coefficient sets differ between imputations for ", outcome, "/", model,
         "/", exposure_form, call. = FALSE)
  }
  coefficients <- lapply(fit_objects, `[[`, "coefficients")
  variances <- lapply(fit_objects, `[[`, "variance")
  df_complete <- min(vapply(fit_objects, `[[`, numeric(1), "residual_df"))
  if (!is.finite(df_complete) || df_complete <= 0) stop("Invalid finite complete-data df")
  if (length(fit_objects) == 1L) {
    pooled <- list(
      coefficients = coefficients[[1]], variance = variances[[1]],
      df = rep(df_complete, length(coefficients[[1]])),
      missinfo = rep(0, length(coefficients[[1]]))
    )
  } else {
    pooled <- mitools::MIcombine(
      results = coefficients, variances = variances, df.complete = df_complete
    )
  }
  beta <- as.numeric(pooled$coefficients)
  names(beta) <- names(pooled$coefficients)
  se <- sqrt(diag(pooled$variance))
  df <- pooled$df
  if (length(df) == 1L) df <- rep(df, length(beta))
  crit <- ifelse(is.finite(df), stats::qt(0.975, df), stats::qnorm(0.975))
  p <- ifelse(
    is.finite(df), 2 * stats::pt(abs(beta / se), df = df, lower.tail = FALSE),
    2 * stats::pnorm(abs(beta / se), lower.tail = FALSE)
  )
  fmi <- pooled$missinfo
  if (is.null(fmi)) fmi <- rep(NA_real_, length(beta))
  if (length(fmi) == 1L) fmi <- rep(fmi, length(beta))
  low <- beta - crit * se
  high <- beta + crit * se
  if (identical(family, "binomial")) {
    beta <- exp(beta); low <- exp(low); high <- exp(high)
    scale <- "OR"
  } else {
    scale <- "beta"
  }
  data.frame(
    outcome = outcome, model = model, exposure_form = exposure_form,
    term = names(pooled$coefficients), estimate = beta, std_error_link = se,
    conf_low = low, conf_high = high, p_value = p, df = df, fmi = fmi,
    complete_data_df = df_complete, effect_scale = scale, status = "OK",
    stringsAsFactors = FALSE
  )
}

build_model_formula <- function(outcome_var, exposure_terms, covariates) {
  stats::reformulate(c(exposure_terms, covariates), response = outcome_var)
}

fit_mi_model_core <- function(dataset_getter, n_datasets, outcome_name, model_name,
                              exposure_form = c("continuous", "quartile", "trend"),
                              complete_case = FALSE) {
  exposure_form <- match.arg(exposure_form)
  spec <- OUTCOMES[[outcome_name]]
  covars <- if (identical(outcome_name, "DKD")) {
    MODEL_COVARIATES_DKD[[model_name]]
  } else {
    MODEL_COVARIATES[[model_name]]
  }
  exposure_terms <- switch(
    exposure_form,
    continuous = "upf_per10",
    quartile = "upf_quartile",
    trend = "upf_q_median_per10"
  )
  form <- build_model_formula(spec$variable, exposure_terms, covars)
  fits <- vector("list", n_datasets)
  errors <- rep(NA_character_, n_datasets)
  for (i in seq_len(n_datasets)) {
    dat <- dataset_getter(i)
    if (complete_case) {
      needed <- unique(c(spec$variable, covars, exposure_terms))
      dat[[spec$domain]] <- dat[[spec$domain]] %in% TRUE &
        stats::complete.cases(dat[, needed, drop = FALSE])
    }
    fits[[i]] <- tryCatch(
      fit_svy_formula(dat, form, spec),
      error = function(e) { errors[[i]] <<- conditionMessage(e); NULL }
    )
  }
  pooled <- pool_svy_fits(fits, outcome_name, model_name, exposure_form, spec$family)
  good <- !vapply(fits, is.null, logical(1))
  metadata <- data.frame(
    outcome = outcome_name, model = model_name, exposure_form = exposure_form,
    imputation = seq_along(fits), success = good, error = errors,
    n = vapply(fits, function(z) if (is.null(z)) NA_integer_ else z$n, integer(1)),
    events = vapply(fits, function(z) if (is.null(z)) NA_integer_ else z$events, integer(1)),
    design_df = vapply(fits, function(z) if (is.null(z)) NA_real_ else z$design_df, numeric(1)),
    residual_df = vapply(fits, function(z) if (is.null(z)) NA_real_ else z$residual_df,
                         numeric(1)),
    n_psu = vapply(fits, function(z) if (is.null(z)) NA_integer_ else z$n_psu, integer(1)),
    n_strata = vapply(fits, function(z) if (is.null(z)) NA_integer_ else z$n_strata, integer(1)),
    stringsAsFactors = FALSE
  )
  list(result = pooled, metadata = metadata, formula = form)
}

fit_mi_model <- function(datasets, outcome_name, model_name,
                         exposure_form = c("continuous", "quartile", "trend"),
                         complete_case = FALSE) {
  fit_mi_model_core(
    dataset_getter = function(i) datasets[[i]], n_datasets = length(datasets),
    outcome_name = outcome_name, model_name = model_name,
    exposure_form = exposure_form, complete_case = complete_case
  )
}

fit_mi_model_stream <- function(frame, imp, outcome_name, model_name,
                                exposure_form = c("continuous", "quartile", "trend"),
                                complete_case = FALSE, transform_data = NULL) {
  frame <- derive_fixed_exposure_categories(frame)
  getter <- function(i) {
    d <- completed_analysis_dataset(frame, imp, i)
    if (!is.null(transform_data)) d <- transform_data(d)
    d
  }
  fit_mi_model_core(
    dataset_getter = getter, n_datasets = imp$m, outcome_name = outcome_name,
    model_name = model_name, exposure_form = exposure_form,
    complete_case = complete_case
  )
}

exposure_rows <- function(result) {
  result[grepl("^upf_per10$|^upf_quartileQ[2-4]$|^upf_q_median_per10$",
               result$term), , drop = FALSE]
}
)--------------------"

upf_source[["KN/source_snapshot/04_build_analysis_frame.R"]] <- r"--------------------(# Build the survey-design universe, person-level NOVA exposure, four outcomes,

.this_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "KN/source_snapshot/04_build_analysis_frame.R"), error = function(e) NULL)
if (is.null(.this_file)) {
  .file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.file_arg)) .this_file <- sub("^--file=", "", .file_arg[[1]])
}
if (is.null(.this_file)) .this_file <- "04_build_analysis_frame.R"
.code_dir <- dirname(normalizePath(.this_file, mustWork = FALSE))
source(file.path(.code_dir, "00_dependencies.R"), local = FALSE)
source(file.path(.code_dir, "00_functions.R"), local = FALSE)

ADJUDICATED_NOVA_MAP <- ACTIVE_NOVA_MAP

PERSON_VARS <- c(
  "ID", "year", "age", "sex", "edu", "kstrata", "psu", "wt_tot",
  "N_EN", "N_INTK", "N_WAT_C", "HE_prg", "N_PRG", "N_BF",
  "BS1_1", "BS3_1", "BD1",
  "BD1_11", "BD2_1", "BD2_14", "pa_aerobic", "HE_crea",
  "HE_Ualb", "HE_Ualb_etc", "HE_Ucrea", "HE_Ucrea_etc",
  "HE_HbA1c", "HE_glu", "DE1_dg",
  "DE1_31", "DE1_32", "HE_sbp", "HE_dbp", "DI1_2",
  "HE_chol", "HE_BMI", "DI2_2"
)

RC_BUILD_VARS <- c(
  "ID", "year", "N_DCODE", "N_DNAME", "N_FCODE", "N_FNAME", "N_ENAME",
  "N_FD_NS", "N_PRODUCT", "N_CNAME", "N_KINDG1", "NF_INTK", "NF_EN"
)

normalize_kn_text <- function(x) {
  x <- enc2utf8(as.character(x))
  x[is.na(x)] <- ""
  trimws(stringr::str_replace_all(x, "[[:space:]]+", " "))
}

normalize_kn_code <- function(x, width = 5L) {
  z <- gsub("[^0-9]", "", as.character(x))
  z[!nzchar(z)] <- NA_character_
  ifelse(is.na(z), NA_character_, stringr::str_pad(z, width, pad = "0"))
}

derive_smoking_knhanes <- function(lifetime, current) {
  out <- rep(NA_character_, length(lifetime))
  out[lifetime %in% c(1, 3)] <- "Never"
  out[lifetime == 2 & current == 3] <- "Former"
  out[lifetime == 2 & current %in% c(1, 2)] <- "Current"
  factor(out, levels = c("Never", "Former", "Current"))
}

derive_alcohol_knhanes <- function(ever, frequency, amount_exact, amount_group, sex) {
  freq_week <- rep(NA_real_, length(ever))
  freq_mid <- c(`1` = 0, `2` = 0.5 / 4.345, `3` = 1 / 4.345,
                `4` = 3 / 4.345, `5` = 2.5, `6` = 5.5)
  hit <- as.character(frequency) %in% names(freq_mid)
  freq_week[hit] <- unname(freq_mid[as.character(frequency[hit])])
  freq_week[ever == 1] <- 0

  amount <- suppressWarnings(as.numeric(amount_exact))
  amount[!is.finite(amount) | amount <= 0 | amount >= 88] <- NA_real_
  amount_mid <- c(`1` = 1.5, `2` = 3.5, `3` = 5.5, `4` = 8, `5` = 10)
  use_group <- is.na(amount) & as.character(amount_group) %in% names(amount_mid)
  amount[use_group] <- unname(amount_mid[as.character(amount_group[use_group])])
  drinks_week <- freq_week * amount
  drinks_week[freq_week == 0] <- 0
  threshold <- ifelse(sex == "Female", 7, ifelse(sex == "Male", 14, NA_real_))
  out <- ifelse(
    !is.finite(drinks_week) | !is.finite(threshold), NA_character_,
    ifelse(drinks_week < 1, "None_or_very_low",
           ifelse(drinks_week <= threshold, "Low_to_moderate", "Higher"))
  )
  factor(out, levels = c("None_or_very_low", "Low_to_moderate", "Higher"))
}

read_adjudicated_nova_map <- function(path = ADJUDICATED_NOVA_MAP) {
  if (!file.exists(path)) stop("Missing adjudicated NOVA map: ", path, call. = FALSE)
  x <- data.table::fread(
    path, encoding = "UTF-8",
    colClasses = c(N_DCODE = "character", N_FCODE = "character",
                   N_KINDG1 = "character")
  )
  required <- c(
    "year", "N_DCODE", "N_DNAME", "N_FCODE", "N_FNAME", "N_ENAME", "N_FD_NS",
    "N_PRODUCT", "N_CNAME", "N_KINDG1", "nova_final", "plain_water_final",
    "is_ssb_final"
  )
  assert_columns(x, required, "adjudicated NOVA map")
  x$N_DCODE <- normalize_kn_text(x$N_DCODE)
  x$N_FCODE <- normalize_kn_code(x$N_FCODE)
  x$N_KINDG1 <- normalize_kn_text(x$N_KINDG1)
  x$nova_final <- suppressWarnings(as.integer(x$nova_final))
  x$plain_water_final <- as.logical(x$plain_water_final)
  x$is_ssb_final <- as.logical(x$is_ssb_final)
  if (any(!x$nova_final %in% 1:4 & !as.logical(x$plain_water_final))) {
    stop("Every non-water NOVA context must have nova_final 1-4", call. = FALSE)
  }
  key <- c("year", "N_DCODE", "N_DNAME", "N_FCODE", "N_FNAME", "N_ENAME",
           "N_FD_NS", "N_PRODUCT", "N_CNAME", "N_KINDG1")
  if (anyDuplicated(x[, ..key])) stop("Duplicate context in adjudicated NOVA map", call. = FALSE)
  x
}

normalize_rc_context <- function(x, year) {
  x <- data.table::as.data.table(x)
  if (!"N_ENAME" %in% names(x)) x[, N_ENAME := ""]
  for (v in c("N_DCODE", "N_DNAME", "N_FCODE", "N_FNAME", "N_ENAME", "N_FD_NS",
              "N_PRODUCT", "N_CNAME", "N_KINDG1")) {
    data.table::set(x, j = v, value = normalize_kn_text(x[[v]]))
  }
  x[, N_FCODE := normalize_kn_code(N_FCODE)]
  x[, year := as.integer(year)]
  x[, NF_INTK := suppressWarnings(as.numeric(NF_INTK))]
  x[, NF_EN := suppressWarnings(as.numeric(NF_EN))]
  x
}

collapse_nova_person_year <- function(year, rc_path, nova_map) {
  message("PERSON_NOVA_BEGIN year=", year)
  rc <- normalize_rc_context(read_sas_selected(rc_path, RC_BUILD_VARS), year)
  map <- data.table::copy(nova_map[nova_map$year == year])
  key <- c("year", "N_DCODE", "N_DNAME", "N_FCODE", "N_FNAME", "N_ENAME",
           "N_FD_NS", "N_PRODUCT", "N_CNAME", "N_KINDG1")
  data.table::setkeyv(rc, key)
  data.table::setkeyv(map, key)
  joined <- map[rc]
  if (nrow(joined) != nrow(rc)) stop("NOVA context join changed row count", call. = FALSE)
  missing_map <- is.na(joined$nova_final) & !joined$plain_water_final %in% TRUE
  if (any(missing_map)) {
    stop("Unmapped 24RC records in year ", year, ": ", sum(missing_map), call. = FALSE)
  }
  joined[, classified_g := ifelse(plain_water_final %in% TRUE, 0, NF_INTK)]
  joined[, upf_g := ifelse(nova_final == 4L & !plain_water_final %in% TRUE, NF_INTK, 0)]
  joined[, ssb_upf_g := ifelse(is_ssb_final %in% TRUE & nova_final == 4L, NF_INTK, 0)]
  for (k in 1:4) {
    joined[, paste0("nova", k, "_g") := ifelse(
      nova_final == k & !plain_water_final %in% TRUE, NF_INTK, 0
    )]
  }
  out <- joined[, .(
    item_rows = .N,
    item_intake_g = sum(NF_INTK, na.rm = TRUE),
    item_energy_kcal = sum(NF_EN, na.rm = TRUE),
    plain_water_g = sum(ifelse(plain_water_final %in% TRUE, NF_INTK, 0), na.rm = TRUE),
    nova1_g = sum(nova1_g, na.rm = TRUE),
    nova2_g = sum(nova2_g, na.rm = TRUE),
    nova3_g = sum(nova3_g, na.rm = TRUE),
    nova4_g = sum(nova4_g, na.rm = TRUE),
    ssb_upf_g = sum(ssb_upf_g, na.rm = TRUE),
    classified_nonwater_g = sum(classified_g, na.rm = TRUE)
  ), by = ID]
  out[, upf_grams_pct := ifelse(classified_nonwater_g > 0,
                                 100 * nova4_g / classified_nonwater_g, NA_real_)]
  out[, non_ssb_denominator_g := classified_nonwater_g - ssb_upf_g]
  out[, upf_non_ssb_pct := ifelse(non_ssb_denominator_g > 0,
                                   100 * (nova4_g - ssb_upf_g) /
                                     non_ssb_denominator_g, NA_real_)]
  out[, year := as.integer(year)]
  message("PERSON_NOVA_END year=", year, " people=", nrow(out))
  as.data.frame(out)
}

read_person_year <- function(year, all_path) {
  x <- as.data.frame(read_sas_selected(all_path, PERSON_VARS))
  for (v in c("HE_Ualb_etc", "HE_Ucrea_etc")) {
    if (!v %in% names(x)) x[[v]] <- NA_character_
  }
  assert_columns(x, c("ID", "year", "age", "sex", "kstrata", "psu", "wt_tot"),
                 paste0(year, " ALL"))
  if (anyDuplicated(x$ID)) stop("Duplicate person ID in year ", year, call. = FALSE)
  x$year <- as.integer(year)
  x
}

standardize_person_frame <- function(d) {
  d$ID <- as.character(d$ID)
  d$person_uid <- paste(d$year, d$ID, sep = "_")
  d$age <- suppressWarnings(as.numeric(d$age))
  d$sex <- factor(d$sex, levels = c(1, 2), labels = c("Male", "Female"))
  d$education <- factor(
    ifelse(d$edu %in% c(1, 2), "Below_high_school",
           ifelse(d$edu == 3, "High_school", ifelse(d$edu == 4,
                                                     "Above_high_school", NA_character_))),
    levels = c("Below_high_school", "High_school", "Above_high_school")
  )
  d$smoking <- derive_smoking_knhanes(d$BS1_1, d$BS3_1)
  d$alcohol <- derive_alcohol_knhanes(d$BD1, d$BD1_11, d$BD2_14, d$BD2_1, d$sex)
  d$physical_activity <- factor(
    ifelse(d$pa_aerobic == 1, "Yes", ifelse(d$pa_aerobic == 0, "No", NA_character_)),
    levels = c("No", "Yes")
  )
  d$bmi <- suppressWarnings(as.numeric(d$HE_BMI))
  d$sbp <- suppressWarnings(as.numeric(d$HE_sbp))
  d$dbp <- suppressWarnings(as.numeric(d$HE_dbp))
  d$hypertension_original <- derive_hypertension_composite(
    d$sbp, d$dbp, d$DI1_2
  )
  d$diabetes_original <- derive_diabetes_composite(
    d$DE1_dg, d$DE1_31, d$DE1_32,
    suppressWarnings(as.numeric(d$HE_glu)), suppressWarnings(as.numeric(d$HE_HbA1c))
  )
  d$high_cholesterol_original <- derive_high_cholesterol_composite(
    suppressWarnings(as.numeric(d$HE_chol)), d$DI2_2
  )
  d$hypertension <- factor(
    ifelse(d$hypertension_original == 1L, "Yes",
           ifelse(d$hypertension_original == 0L, "No", NA_character_)),
    levels = c("No", "Yes")
  )
  d$diabetes <- factor(
    ifelse(d$diabetes_original == 1L, "Yes",
           ifelse(d$diabetes_original == 0L, "No", NA_character_)),
    levels = c("No", "Yes")
  )
  d$high_cholesterol <- factor(
    ifelse(d$high_cholesterol_original == 1L, "Yes",
           ifelse(d$high_cholesterol_original == 0L, "No", NA_character_)),
    levels = c("No", "Yes")
  )
  d$energy_kcal <- suppressWarnings(as.numeric(d$N_EN))
  d$official_intake_g <- suppressWarnings(as.numeric(d$N_INTK))
  pregnant_exam <- ifelse(
    d$HE_prg == 1, TRUE,
    ifelse(d$HE_prg %in% c(0, 8), FALSE, NA)
  )
  pregnant_nutrition_fallback <- ifelse(
    d$N_PRG == 1, TRUE,
    ifelse(d$N_PRG %in% c(2, 8), FALSE, NA)
  )
  d$pregnant <- ifelse(!is.na(pregnant_exam), pregnant_exam,
                       pregnant_nutrition_fallback)
  d$pregnancy_source <- ifelse(
    !is.na(pregnant_exam), "HE_prg",
    ifelse(!is.na(pregnant_nutrition_fallback), "N_PRG_fallback", "unknown")
  )
  d$breastfeeding_nutrition <- ifelse(
    d$year <= 2021 & d$N_PRG == 2, TRUE,
    ifelse(d$year >= 2022 & d$N_BF == 1, TRUE,
           ifelse((d$year <= 2021 & d$N_PRG %in% c(1, 8)) |
                    (d$year >= 2022 & d$N_BF %in% c(2, 8)), FALSE, NA))
  )
  d$survey_year <- factor(d$year, levels = YEARS)
  d$analysis_weight <- suppressWarnings(as.numeric(d$wt_tot)) / length(YEARS)
  d$strata_uid <- factor(as.character(d$kstrata))
  d$psu_uid <- factor(as.character(d$psu))
  d
}

derive_knhanes_outcomes <- function(d) {
  d$creatinine_mg_dl <- suppressWarnings(as.numeric(d$HE_crea))
  limits_for <- function(analyte) {
    rules <- LAB_LIMIT_RULES[LAB_LIMIT_RULES$analyte == analyte, ]
    idx <- match(d$year, rules$year)
    if (anyNA(idx)) stop("No laboratory-limit rule for ", analyte, " and a study year")
    rules[idx, ]
  }
  ualb_limits <- limits_for("urine_albumin")
  ucrea_limits <- limits_for("urine_creatinine")
  ualb <- apply_lab_limits_from_qualifier(
    d$HE_Ualb, d$HE_Ualb_etc, ualb_limits$default_lower, ualb_limits$default_upper
  )
  ucrea <- apply_lab_limits_from_qualifier(
    d$HE_Ucrea, d$HE_Ucrea_etc, ucrea_limits$default_lower, ucrea_limits$default_upper
  )
  d$urine_albumin_ug_ml <- ualb$value_primary
  d$urine_albumin_ug_ml_alt_lod <- ualb$value_alternative
  d$urine_creatinine_mg_dl <- ucrea$value_primary
  d$urine_creatinine_mg_dl_alt_lod <- ucrea$value_alternative
  d$urine_albumin_below_lod <- ualb$below_lod
  d$urine_albumin_above_uloq <- ualb$above_uloq
  d$urine_creatinine_below_lod <- ucrea$below_lod
  d$urine_albumin_qualifier_boundary <- ualb$boundary_used
  d$urine_creatinine_qualifier_boundary <- ucrea$boundary_used
  d$egfr_2021 <- calc_egfr_2021(d$creatinine_mg_dl, d$age, d$sex == "Female")
  d$acr_mg_g <- ifelse(
    is.finite(d$urine_albumin_ug_ml) & is.finite(d$urine_creatinine_mg_dl) &
      d$urine_creatinine_mg_dl > 0,
    100 * d$urine_albumin_ug_ml / d$urine_creatinine_mg_dl,
    NA_real_
  )
  d$log_uacr <- ifelse(is.finite(d$acr_mg_g) & d$acr_mg_g > 0,
                       log(d$acr_mg_g), NA_real_)
  d$ckd <- derive_ckd_three_state(d$egfr_2021, d$acr_mg_g)
  d$acr_mg_g_alt_lod <- ifelse(
    is.finite(d$urine_albumin_ug_ml_alt_lod) &
      is.finite(d$urine_creatinine_mg_dl_alt_lod) &
      d$urine_creatinine_mg_dl_alt_lod > 0,
    100 * d$urine_albumin_ug_ml_alt_lod / d$urine_creatinine_mg_dl_alt_lod,
    NA_real_
  )
  d$log_uacr_alt_lod <- ifelse(
    is.finite(d$acr_mg_g_alt_lod) & d$acr_mg_g_alt_lod > 0,
    log(d$acr_mg_g_alt_lod), NA_real_
  )
  d$ckd_alt_lod <- derive_ckd_three_state(d$egfr_2021, d$acr_mg_g_alt_lod)
  d$dkd <- ifelse(d$diabetes_original == 1L, d$ckd, NA_integer_)
  d$dkd_alt_lod <- ifelse(d$diabetes_original == 1L, d$ckd_alt_lod, NA_integer_)
  d
}

laboratory_limit_summary <- function(d) {
  flags <- c("urine_albumin_below_lod", "urine_albumin_above_uloq",
             "urine_creatinine_below_lod")
  flag_rows <- dplyr::bind_rows(lapply(flags, function(v) {
    data.frame(year = d$year, flag = v, present = d[[v]] %in% TRUE) |>
      dplyr::group_by(year, flag) |>
      dplyr::summarise(n = sum(present), .groups = "drop")
  }))
  impacts <- data.frame(
    outcome = c("CKD", "DKD"),
    changed_primary_vs_lod_half = c(
      sum(d$ckd != d$ckd_alt_lod, na.rm = TRUE),
      sum(d$dkd != d$dkd_alt_lod, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  list(flags = flag_rows, impacts = impacts)
}

make_domains <- function(d) {
  not_pregnant <- d$pregnant %in% FALSE
  d$base_domain <- d$age >= 20 & not_pregnant & !is.na(d$item_rows) &
    is.finite(d$classified_nonwater_g) & d$classified_nonwater_g > 0 &
    is.finite(d$upf_grams_pct) & is.finite(d$energy_kcal)
  d$egfr_observed <- d$base_domain & !is.na(d$egfr_2021)
  d$uacr_observed <- d$base_domain & !is.na(d$log_uacr)
  d$acr_observed <- d$base_domain & !is.na(d$acr_mg_g)
  d$kidney_observed <- d$base_domain & !is.na(d$ckd)
  d$dkd_domain <- d$kidney_observed & d$diabetes_original == 1L
  d
}

reconcile_person_totals <- function(d) {
  z <- d[!is.na(d$item_rows), , drop = FALSE]
  z$intake_diff_g <- z$item_intake_g - z$official_intake_g
  z$energy_diff_kcal <- z$item_energy_kcal - z$energy_kcal
  tol_intake <- pmax(NOVA_RECONCILIATION_ABS_TOL_G,
                     NOVA_RECONCILIATION_REL_TOL * abs(z$official_intake_g))
  tol_energy <- pmax(NOVA_RECONCILIATION_ABS_TOL_G,
                     NOVA_RECONCILIATION_REL_TOL * abs(z$energy_kcal))
  z$intake_pass <- is.finite(z$intake_diff_g) & abs(z$intake_diff_g) <= tol_intake
  z$energy_pass <- is.finite(z$energy_diff_kcal) & abs(z$energy_diff_kcal) <= tol_energy
  if (any(!z$intake_pass) || any(!z$energy_pass)) {
    stop("Person-level 24RC reconciliation failed", call. = FALSE)
  }
  z[c("person_uid", "year", "ID", "item_rows", "official_intake_g", "item_intake_g",
      "intake_diff_g", "energy_kcal", "item_energy_kcal", "energy_diff_kcal",
      "intake_pass", "energy_pass")]
}

outcome_domain_summary <- function(d) {
  full_design <- survey::svydesign(
    ids = ~psu_uid, strata = ~strata_uid, weights = ~analysis_weight,
    nest = TRUE, data = d
  )
  do.call(rbind, lapply(names(OUTCOMES), function(outcome) {
    spec <- OUTCOMES[[outcome]]
    domain <- d[[spec$domain]] %in% TRUE
    y <- d[[spec$variable]]
    design_domain <- full_design
    design_domain$variables$.domain_tmp <- domain
    data.frame(
      outcome = outcome,
      family = spec$family,
      domain = spec$domain,
      n = sum(domain),
      events = if (spec$family == "binomial") sum(y[domain] == 1, na.rm = TRUE) else NA_integer_,
      n_psu = length(unique(d$psu_uid[domain])),
      n_strata = length(unique(d$strata_uid[domain])),
      design_df = survey::degf(subset(design_domain, .domain_tmp)),
      stringsAsFactors = FALSE
    )
  }))
}

stage_build_analysis_frame <- function() {
  require_nova_gate("build_analysis_frame")
  assert_files(c(YEAR_FILES$all, YEAR_FILES$rc24, ADJUDICATED_NOVA_MAP,
                 OUTCOME_DEFINITION_MANIFEST),
               "analysis-frame inputs")
  validate_outcome_definition_manifest()
  run_id <- timestamp_id("N4_PERSON_FRAME")
  derived_run <- file.path(ANALYSIS_DIRS$derived, run_id)
  qc_run <- file.path(ANALYSIS_DIRS$qc, run_id)
  log_run <- file.path(ANALYSIS_DIRS$logs, run_id)
  dir.create(derived_run, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_run, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_run, recursive = TRUE, showWarnings = FALSE)
  sink(file.path(log_run, "build_analysis_frame.log"), split = TRUE)
  on.exit(sink(), add = TRUE)
  cat("run_id=", run_id, "\n", sep = "")

  map <- read_adjudicated_nova_map()
  exposures <- dplyr::bind_rows(lapply(seq_len(nrow(YEAR_FILES)), function(i) {
    collapse_nova_person_year(YEAR_FILES$year[[i]], YEAR_FILES$rc24[[i]], map)
  }))
  people <- dplyr::bind_rows(lapply(seq_len(nrow(YEAR_FILES)), function(i) {
    read_person_year(YEAR_FILES$year[[i]], YEAR_FILES$all[[i]])
  }))
  if (anyDuplicated(paste(people$year, people$ID))) stop("Duplicate year-ID in ALL", call. = FALSE)
  d <- merge(people, exposures, by = c("year", "ID"), all.x = TRUE, sort = FALSE)
  d <- make_domains(derive_knhanes_outcomes(standardize_person_frame(d)))
  d$analysis_scope <- ANALYSIS_SCOPE
  d$nova_map_sha256 <- sha256_file(ADJUDICATED_NOVA_MAP)
  d$outcome_definition_manifest_sha256 <- sha256_file(OUTCOME_DEFINITION_MANIFEST)
  design_universe <- is.finite(d$analysis_weight) & d$analysis_weight > 0 &
    !is.na(d$strata_uid) & !is.na(d$psu_uid)
  d <- d[design_universe, , drop = FALSE]
  if (anyDuplicated(d$person_uid)) stop("Duplicate person_uid after design filter", call. = FALSE)
  if (length(unique(d$strata_uid)) != 54L) stop("Official six-year kstrata must have 54 levels")
  yearly_psu <- tapply(as.character(d$psu_uid), d$year, function(x) length(unique(x)))
  if (!identical(as.integer(yearly_psu), c(192L, 166L, 190L, 184L, 192L, 192L)) ||
      length(unique(d$psu_uid)) != 1116L) {
    stop("Positive-wt_tot design universe must retain 192/166/190/184/192/192 and 1116 PSUs")
  }
  design_check <- survey::svydesign(
    ids = ~psu_uid, strata = ~strata_uid, weights = ~analysis_weight,
    nest = TRUE, data = d
  )
  if (!identical(as.integer(survey::degf(design_check)), 1062L)) {
    stop("Official six-year design degrees of freedom must equal 1062")
  }

  reconciliation <- reconcile_person_totals(d)
  domains <- outcome_domain_summary(d)
  lab_limits <- laboratory_limit_summary(d)
  pregnancy <- dplyr::bind_rows(lapply(YEARS, function(y) {
    z <- d[d$year == y, , drop = FALSE]
    data.frame(
      year = y,
      status = c("pregnant_final", "not_pregnant_final", "pregnancy_unknown",
                 "source_HE_prg", "source_N_PRG_fallback", "breastfeeding_nutrition"),
      n = c(sum(z$pregnant %in% TRUE), sum(z$pregnant %in% FALSE), sum(is.na(z$pregnant)),
            sum(z$pregnancy_source == "HE_prg"),
            sum(z$pregnancy_source == "N_PRG_fallback"),
            sum(z$breastfeeding_nutrition %in% TRUE)),
      stringsAsFactors = FALSE
    )
  }))
  flow <- data.frame(
    step = c("ALL participants", "Positive wt_tot/design universe", "Age >=20",
             "Not pregnant", "Linked day-1 24RC", "Positive classified non-water grams",
             "Base analysis domain"),
    n = c(nrow(people), nrow(d), sum(d$age >= 20, na.rm = TRUE),
          sum(d$age >= 20 & d$pregnant %in% FALSE, na.rm = TRUE),
          sum(d$age >= 20 & d$pregnant %in% FALSE & !is.na(d$item_rows), na.rm = TRUE),
          sum(d$age >= 20 & d$pregnant %in% FALSE & !is.na(d$item_rows) &
                d$classified_nonwater_g > 0, na.rm = TRUE), sum(d$base_domain)),
    stringsAsFactors = FALSE
  )

  frame_path <- file.path(derived_run, "analysis_frame_pre_mice.rds")
  saveRDS(d, frame_path, compress = "xz")
  input_paths <- c(YEAR_FILES$all, YEAR_FILES$rc24, ADJUDICATED_NOVA_MAP,
                   OUTCOME_DEFINITION_MANIFEST,
                   unname(CODEBOOK_PATHS), unname(USER_GUIDE_PATHS),
                   CODE_REVIEW_MARKER, NOVA_GATE_MARKER)
  input_manifest <- data.frame(
    path = input_paths, bytes = file.info(input_paths)$size,
    sha256 = vapply(input_paths, sha256_file, character(1)), stringsAsFactors = FALSE
  )
  input_manifest_path <- file.path(log_run, "input_manifest.csv")
  write_csv_atomic(input_manifest, input_manifest_path)
  write_csv_atomic(reconciliation, file.path(qc_run, "person_reconciliation.csv"))
  write_csv_atomic(domains, file.path(qc_run, "outcome_domain_summary.csv"))
  write_csv_atomic(flow, file.path(qc_run, "sample_flow.csv"))
  write_csv_atomic(lab_limits$flags, file.path(qc_run, "laboratory_limit_counts.csv"))
  write_csv_atomic(lab_limits$impacts, file.path(qc_run, "laboratory_limit_classification_impact.csv"))
  write_csv_atomic(pregnancy, file.path(qc_run, "pregnancy_breastfeeding_counts.csv"))
  manifest_paths <- c(frame_path, input_manifest_path,
                      file.path(qc_run, "person_reconciliation.csv"),
                      file.path(qc_run, "outcome_domain_summary.csv"),
                      file.path(qc_run, "sample_flow.csv"),
                      file.path(qc_run, "laboratory_limit_counts.csv"),
                      file.path(qc_run, "laboratory_limit_classification_impact.csv"),
                      file.path(qc_run, "pregnancy_breastfeeding_counts.csv"))
  manifest <- data.frame(
    path = manifest_paths,
    bytes = file.info(manifest_paths)$size,
    sha256 = vapply(manifest_paths, sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
  write_csv_atomic(manifest, file.path(log_run, "output_manifest.csv"))
  writeLines(capture.output(sessionInfo()), file.path(log_run, "sessionInfo.txt"))
  commit_latest_marker(run_id, file.path(ANALYSIS_DIRS$derived, "LATEST_ANALYSIS_FRAME_RUN_ID.txt"),
                       c(manifest_paths, input_manifest_path,
                         file.path(log_run, "output_manifest.csv")))
  cat("base_domain_n=", sum(d$base_domain), "\n", sep = "")
  cat("run_id=", run_id, "\n", sep = "")
  invisible(list(run_id = run_id, frame = frame_path, domains = domains, flow = flow))
}

if (sys.nframe() == 0L) stage_build_analysis_frame()
)--------------------"

upf_source[["KN/source_snapshot/05_imputation.R"]] <- r"--------------------(# Multiple imputation of harmonized covariates only.

.this_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "KN/source_snapshot/05_imputation.R"), error = function(e) NULL)
if (is.null(.this_file)) {
  .file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.file_arg)) .this_file <- sub("^--file=", "", .file_arg[[1]])
}
if (is.null(.this_file)) .this_file <- "05_imputation.R"
.code_dir <- dirname(normalizePath(.this_file, mustWork = FALSE))
source(file.path(.code_dir, "00_dependencies.R"), local = FALSE)
source(file.path(.code_dir, "00_functions.R"), local = FALSE)

latest_frame_path <- function() {
  marker <- file.path(ANALYSIS_DIRS$derived, "LATEST_ANALYSIS_FRAME_RUN_ID.txt")
  assert_files(marker, "latest analysis-frame marker")
  run_id <- marker_run_id(marker)
  path <- file.path(ANALYSIS_DIRS$derived, run_id, "analysis_frame_pre_mice.rds")
  assert_files(path, "analysis frame")
  kv <- read_key_value_marker(marker)
  if (!identical(unname(kv[["status"]]), "PASSED") ||
      is.null(kv[["output_manifest"]]) || is.null(kv[["output_manifest_sha256"]]) ||
      !file.exists(kv[["output_manifest"]]) ||
      !identical(sha256_file(kv[["output_manifest"]]),
                 unname(kv[["output_manifest_sha256"]]))) {
    stop("Analysis-frame latest marker provenance is invalid", call. = FALSE)
  }
  verify_manifest_entry(path, kv[["output_manifest"]])
  path
}

fill_auxiliary <- function(x, binary = FALSE) {
  missing <- !is.finite(x)
  if (binary) {
    observed <- x[!missing]
    fill <- if (length(observed)) as.numeric(mean(observed) >= 0.5) else 0
  } else {
    observed <- x[!missing]
    fill <- if (length(observed)) stats::median(observed) else 0
  }
  list(value = ifelse(missing, fill, x), missing = as.integer(missing))
}

build_mice_input <- function(frame) {
  assert_columns(
    frame,
    c("person_uid", "base_domain", "age", "sex", "survey_year", "energy_kcal",
      "upf_grams_pct", "analysis_weight", "kstrata", "psu_uid",
      IMPUTED_COVARIATES,
      vapply(OUTCOMES, `[[`, character(1), "variable")),
    "analysis frame"
  )
  x <- frame[c(
    "person_uid", "base_domain", "age", "sex", "survey_year", "energy_kcal",
    "upf_grams_pct", "analysis_weight", "kstrata", "psu_uid", IMPUTED_COVARIATES
  )]
  x$.row_id <- seq_len(nrow(x))
  if (any(!is.finite(x$analysis_weight) | x$analysis_weight <= 0)) {
    stop("MICE input contains invalid survey weights")
  }
  x$log_analysis_weight <- log(x$analysis_weight)
  x$strata_mi <- factor(x$kstrata)
  x$psu_cluster_n <- ave(x$.row_id, x$psu_uid, FUN = length)
  x$psu_mean_log_weight <- ave(x$log_analysis_weight, x$psu_uid, FUN = mean)
  x$psu_mean_age <- ave(x$age, x$psu_uid, FUN = function(z) mean(z, na.rm = TRUE))
  for (nm in names(OUTCOMES)) {
    spec <- OUTCOMES[[nm]]
    raw <- suppressWarnings(as.numeric(frame[[spec$variable]]))
    aux <- fill_auxiliary(raw, binary = identical(spec$family, "binomial"))
    x[[paste0("aux_", nm)]] <- aux$value
    x[[paste0("aux_", nm, "_missing")]] <- aux$missing
  }
  x
}

configure_mice <- function(x) {
  method <- rep("", ncol(x)); names(method) <- names(x)
  method[c("education", "smoking", "alcohol")] <- "polyreg"
  method[c("physical_activity", "hypertension", "diabetes",
           "high_cholesterol")] <- "logreg"
  method["bmi"] <- "pmm"

  predictor <- matrix(0L, nrow = ncol(x), ncol = ncol(x),
                      dimnames = list(names(x), names(x)))
  redundant_lab_missing_flags <- c(
    "aux_CKD_missing", "aux_DKD_missing"
  )
  if (!all(redundant_lab_missing_flags %in% names(x))) {
    stop("MICE input lacks expected laboratory missingness flags")
  }
  stable_predictors <- setdiff(c(
    "age", "sex", "survey_year", "energy_kcal", "upf_grams_pct",
    "log_analysis_weight", "strata_mi", "psu_cluster_n",
    "psu_mean_log_weight", "psu_mean_age",
    grep("^aux_", names(x), value = TRUE), IMPUTED_COVARIATES
  ), redundant_lab_missing_flags)
  for (target in IMPUTED_COVARIATES) {
    predictor[target, setdiff(stable_predictors, target)] <- 1L
  }
  predictor[, c("person_uid", ".row_id", "base_domain", "analysis_weight",
                "kstrata", "psu_uid")] <- 0L

  where <- matrix(FALSE, nrow = nrow(x), ncol = ncol(x),
                  dimnames = list(NULL, names(x)))
  domain <- x$base_domain %in% TRUE
  for (target in IMPUTED_COVARIATES) {
    where[, target] <- domain & is.na(x[[target]])
  }
  list(method = method, predictor = predictor, where = where, ignore = !domain)
}

validate_imputations <- function(imp, base_domain) {
  checks <- lapply(seq_len(imp$m), function(i) {
    d <- mice::complete(imp, action = i)
    missing_n <- vapply(IMPUTED_COVARIATES, function(v) {
      sum(base_domain & is.na(d[[v]]))
    }, integer(1))
    invalid_n <- vapply(IMPUTED_COVARIATES, function(v) {
      original <- imp$data[[v]]
      if (is.factor(original)) {
        sum(base_domain & !is.na(d[[v]]) & !as.character(d[[v]]) %in% levels(original))
      } else if (identical(v, "bmi")) {
        observed <- suppressWarnings(as.numeric(original))
        completed <- suppressWarnings(as.numeric(d[[v]]))
        range <- range(observed, na.rm = TRUE)
        sum(base_domain & is.finite(completed) &
              (completed < range[[1]] | completed > range[[2]]))
      } else 0L
    }, integer(1))
    data.frame(imputation = i, variable = names(missing_n),
               missing_in_base_domain = unname(missing_n),
               invalid_level_or_range = unname(invalid_n), stringsAsFactors = FALSE)
  })
  out <- dplyr::bind_rows(checks)
  if (any(out$missing_in_base_domain > 0L) || any(out$invalid_level_or_range > 0L)) {
    stop("MICE completion or value-range gate failed", call. = FALSE)
  }
  out
}

imputation_hashes <- function(imp) {
  hashes <- vapply(seq_len(imp$m), function(i) {
    d <- mice::complete(imp, action = i)[IMPUTED_COVARIATES]
    digest::digest(d, algo = "sha256", serialize = TRUE)
  }, character(1))
  out <- data.frame(imputation = seq_len(imp$m), sha256 = hashes,
                    duplicated = duplicated(hashes) | duplicated(hashes, fromLast = TRUE),
                    stringsAsFactors = FALSE)
  if (any(out$duplicated)) stop("MICE produced duplicate completed imputations")
  out
}

chain_diagnostics <- function(imp) {
  cm <- imp$chainMean
  cv <- imp$chainVar
  summarize_array <- function(a, metric) {
    if (is.null(a) || !length(a)) return(data.frame())
    variables <- dimnames(a)[[1]] %||% paste0("slice_", seq_len(dim(a)[[1]]))
    dplyr::bind_rows(lapply(seq_along(variables), function(i) {
      z <- a[i, , , drop = TRUE]
      finite <- is.finite(z)
      data.frame(
        metric = metric, variable = variables[[i]],
        finite_fraction = mean(finite),
        minimum = if (any(finite)) min(z[finite]) else NA_real_,
        maximum = if (any(finite)) max(z[finite]) else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  }
  dplyr::bind_rows(summarize_array(cm, "chainMean"), summarize_array(cv, "chainVar"))
}

mice_convergence_diagnostics <- function(imp) {
  if (!is.finite(imp$iteration) || imp$iteration < 3L) {
    z <- data.frame(
      iteration = as.integer(imp$iteration), variable = IMPUTED_COVARIATES,
      ac = NA_real_, psrf = NA_real_,
      diagnostic_note = "Not estimable when maxit < 3; smoke-test placeholder",
      stringsAsFactors = FALSE
    )
    last <- z
    last$soft_flag <- TRUE
    return(list(all = z, last = last))
  }
  z <- tryCatch(as.data.frame(mice::convergence(imp)), error = function(e) {
    stop("mice::convergence failed: ", conditionMessage(e), call. = FALSE)
  })
  if (!nrow(z)) stop("mice::convergence returned no diagnostics", call. = FALSE)
  variable_col <- intersect(c("vrb", "variable"), names(z))
  iteration_col <- intersect(c(".it", "iteration"), names(z))
  if (!length(variable_col) || !length(iteration_col) ||
      !all(c("ac", "psrf") %in% names(z))) {
    stop("Unexpected mice::convergence output schema", call. = FALSE)
  }
  names(z)[names(z) == variable_col[[1]]] <- "variable"
  names(z)[names(z) == iteration_col[[1]]] <- "iteration"
  z$variable <- as.character(z$variable)
  z$iteration <- suppressWarnings(as.integer(z$iteration))
  target <- z$variable %in% IMPUTED_COVARIATES
  missing_metrics <- any(vapply(split(z[target, ], z$variable[target]), function(q) {
    !any(is.finite(q$ac)) || !any(is.finite(q$psrf))
  }, logical(1)))
  if (!all(IMPUTED_COVARIATES %in% z$variable) ||
      any(!is.finite(z$iteration[target])) ||
      missing_metrics) {
    stop("MICE convergence diagnostics do not cover every target with finite values",
         call. = FALSE)
  }
  last <- dplyr::bind_rows(lapply(split(z[target, ], z$variable[target]), function(q) {
    q <- q[which.max(q$iteration), , drop = FALSE]
    q$soft_flag <- !is.finite(q$ac) | !is.finite(q$psrf) |
      abs(q$ac) > 0.10 | q$psrf > 1.10
    q
  }))
  list(all = z, last = last)
}

enforce_mice_convergence_gate <- function(last, analysis_mode = ANALYSIS_MODE) {
  assert_columns(last, c("variable", "ac", "psrf", "soft_flag"),
                 "MICE last-iteration convergence diagnostics")
  target <- last[last$variable %in% IMPUTED_COVARIATES, , drop = FALSE]
  if (nrow(target) != length(IMPUTED_COVARIATES) || anyDuplicated(target$variable) ||
      !setequal(target$variable, IMPUTED_COVARIATES)) {
    stop("MICE convergence gate does not contain exactly one row per target",
         call. = FALSE)
  }
  if (identical(analysis_mode, "formal")) {
    failed <- !is.finite(target$ac) | !is.finite(target$psrf) |
      as.logical(target$soft_flag) | abs(target$ac) > 0.10 | target$psrf > 1.10
    failed[is.na(failed)] <- TRUE
    if (any(failed)) {
      stop("Formal MICE convergence hard gate failed for: ",
           paste(target$variable[failed], collapse = ", "), call. = FALSE)
    }
  } else if (!identical(analysis_mode, "smoke")) {
    stop("MICE convergence gate received unsupported analysis mode", call. = FALSE)
  }
  invisible(TRUE)
}

audit_mice_logged_events <- function(logged, maxit = NULL) {
  empty <- data.frame(
    it = integer(), im = integer(), dep = character(), meth = character(),
    out = character(), worker_id = integer(), local_im = integer(),
    allowed = logical(), classification = character(), stringsAsFactors = FALSE
  )
  if (is.null(logged) || !nrow(logged)) return(empty)
  assert_columns(logged, c("it", "im", "dep", "meth", "out"),
                 "MICE logged events")
  expected_method <- c(
    education = "polyreg", smoking = "polyreg", alcohol = "polyreg",
    physical_activity = "logreg", hypertension = "logreg",
    diabetes = "logreg", high_cholesterol = "logreg", bmi = "pmm"
  )
  dep <- as.character(logged$dep)
  meth <- as.character(logged$meth)
  out <- trimws(as.character(logged$out))
  it <- suppressWarnings(as.integer(logged$it))
  im <- suppressWarnings(as.integer(logged$im))
  worker_id <- if ("worker_id" %in% names(logged)) {
    suppressWarnings(as.integer(logged$worker_id))
  } else rep(NA_integer_, nrow(logged))
  local_im <- if ("local_im" %in% names(logged)) {
    suppressWarnings(as.integer(logged$local_im))
  } else im
  it_upper_bound <- NULL
  if (!is.null(maxit)) {
    maxit <- suppressWarnings(as.integer(maxit))
    if (length(maxit) != 1L || is.na(maxit) || maxit < 1L) {
      stop("MICE logged-event gate received invalid maxit", call. = FALSE)
    }
    it_upper_bound <- maxit
  }
  year_dummy_pattern <- paste0("^survey_year(", paste(YEARS, collapse = "|"), ")$")
  allowed <- !is.na(it) & it >= 1L & !is.na(im) & im >= 1L &
    dep %in% IMPUTED_COVARIATES &
    meth == unname(expected_method[dep]) &
    grepl(year_dummy_pattern, out)
  if (!is.null(it_upper_bound)) {
    allowed <- allowed & it <= it_upper_bound
  }
  allowed[is.na(allowed)] <- FALSE
  data.frame(
    it = it, im = im, dep = dep, meth = meth, out = out,
    worker_id = worker_id, local_im = local_im, allowed = allowed,
    classification = ifelse(
      allowed, "ALLOWED_REDUNDANT_SURVEY_YEAR_DUMMY_AUTO_REMOVAL",
      "BLOCKING_UNRECOGNIZED_MICE_ACTION"
    ), stringsAsFactors = FALSE
  )
}

enforce_mice_logged_event_gate <- function(audit) {
  assert_columns(audit, c("allowed", "classification", "dep", "out"),
                 "MICE logged-event audit")
  if (nrow(audit) && any(!audit$allowed)) {
    bad <- audit[!audit$allowed, , drop = FALSE]
    stop("Unrecognized MICE automatic action(s): ",
         paste(paste0(bad$dep, "->", bad$out), collapse = " | "),
         call. = FALSE)
  }
  invisible(TRUE)
}

empty_mice_worker_events <- function() {
  data.frame(
    it = integer(), im = integer(), dep = character(), meth = character(),
    out = character(), worker_id = integer(), local_im = integer(),
    stringsAsFactors = FALSE
  )
}

normalize_mice_worker_events <- function(logged, worker_id, im_offset,
                                         worker_m) {
  if (is.null(logged) || !nrow(logged)) return(empty_mice_worker_events())
  assert_columns(logged, c("it", "im", "dep", "meth", "out"),
                 paste0("MICE worker ", worker_id, " logged events"))
  local_im <- suppressWarnings(as.integer(logged$im))
  if (any(is.na(local_im)) || any(local_im < 0L) ||
      any(local_im > as.integer(worker_m))) {
    stop("MICE worker logged-event im index is invalid", call. = FALSE)
  }
  global_im <- local_im
  positive <- local_im > 0L
  global_im[positive] <- local_im[positive] + as.integer(im_offset)
  data.frame(
    it = suppressWarnings(as.integer(logged$it)), im = global_im,
    dep = as.character(logged$dep), meth = as.character(logged$meth),
    out = as.character(logged$out), worker_id = as.integer(worker_id),
    local_im = local_im, stringsAsFactors = FALSE
  )
}

combine_mids_preserving_worker_events <- function(imps, requested_m,
                                                  requested_worker_m = NULL) {
  if (!length(imps) || any(!vapply(imps, mice::is.mids, logical(1)))) {
    stop("Every worker result must be a mids object", call. = FALSE)
  }
  worker_m <- vapply(imps, function(z) as.integer(z$m), integer(1))
  if (is.null(requested_worker_m)) requested_worker_m <- worker_m
  requested_worker_m <- suppressWarnings(as.integer(requested_worker_m))
  if (length(requested_worker_m) != length(worker_m) ||
      any(is.na(requested_worker_m)) || any(requested_worker_m < 1L) ||
      sum(requested_worker_m) != as.integer(requested_m) ||
      any(worker_m < 1L) || sum(worker_m) != as.integer(requested_m) ||
      !identical(worker_m, requested_worker_m)) {
    stop("MICE worker results do not cover the requested imputations", call. = FALSE)
  }
  offsets <- c(0L, head(cumsum(worker_m), -1L))
  logs <- lapply(seq_along(imps), function(i) {
    normalize_mice_worker_events(
      imps[[i]]$loggedEvents, worker_id = i, im_offset = offsets[[i]],
      worker_m = worker_m[[i]]
    )
  })
  event_rows <- vapply(logs, nrow, integer(1))
  coverage <- data.frame(
    worker_id = seq_along(imps), requested_m = requested_worker_m,
    returned_m = worker_m, im_offset = offsets,
    first_global_im = offsets + 1L, last_global_im = offsets + worker_m,
    event_rows = event_rows, collection_status = "COLLECTED_BEFORE_IBIND",
    stringsAsFactors = FALSE
  )
  combined <- imps[[1]]
  if (length(imps) > 1L) {
    for (i in 2:length(imps)) combined <- mice::ibind(combined, imps[[i]])
  }
  combined_events <- do.call(rbind, logs)
  rownames(combined_events) <- NULL
  combined$loggedEvents <- if (nrow(combined_events)) combined_events else NULL
  attr(combined, "knhanes_worker_event_coverage") <- coverage
  attr(combined, "knhanes_worker_events_complete") <- TRUE
  combined
}

enforce_mice_worker_event_coverage <- function(imp, requested_m,
                                               requested_workers,
                                               expected_worker_m = NULL) {
  coverage <- attr(imp, "knhanes_worker_event_coverage")
  required <- c(
    "worker_id", "requested_m", "returned_m", "im_offset",
    "first_global_im", "last_global_im", "event_rows", "collection_status"
  )
  if (!isTRUE(attr(imp, "knhanes_worker_events_complete")) ||
      !is.data.frame(coverage)) {
    stop("MICE worker-event completeness attestation is absent", call. = FALSE)
  }
  assert_columns(coverage, required, "MICE worker-event coverage")
  requested_m <- suppressWarnings(as.integer(requested_m))
  requested_workers <- suppressWarnings(as.integer(requested_workers))
  if (length(requested_m) != 1L || is.na(requested_m) || requested_m < 1L ||
      length(requested_workers) != 1L || is.na(requested_workers) ||
      requested_workers < 1L || requested_workers > requested_m) {
    stop("MICE worker-event gate received invalid requested dimensions",
         call. = FALSE)
  }
  expected_workers <- seq_len(requested_workers)
  if (is.null(expected_worker_m)) {
    expected_worker_m <- rep(requested_m %/% requested_workers,
                             requested_workers)
    remainder <- requested_m %% requested_workers
    if (remainder > 0L) expected_worker_m[seq_len(remainder)] <-
      expected_worker_m[seq_len(remainder)] + 1L
  }
  expected_worker_m <- suppressWarnings(as.integer(expected_worker_m))
  expected_offsets <- c(0L, head(cumsum(expected_worker_m), -1L))
  expected_first <- expected_offsets + 1L
  expected_last <- expected_offsets + expected_worker_m
  numeric_fields <- c(
    "worker_id", "requested_m", "returned_m", "im_offset",
    "first_global_im", "last_global_im", "event_rows"
  )
  if (nrow(coverage) != length(expected_workers) ||
      length(expected_worker_m) != length(expected_workers) ||
      any(is.na(expected_worker_m)) || any(expected_worker_m < 1L) ||
      sum(expected_worker_m) != requested_m ||
      any(vapply(coverage[numeric_fields], function(z) {
        any(is.na(suppressWarnings(as.integer(z))))
      }, logical(1))) ||
      !identical(as.integer(coverage$worker_id), expected_workers) ||
      !identical(as.integer(coverage$requested_m), expected_worker_m) ||
      !identical(as.integer(coverage$returned_m), expected_worker_m) ||
      !identical(as.integer(coverage$im_offset), expected_offsets) ||
      !identical(as.integer(coverage$first_global_im), expected_first) ||
      !identical(as.integer(coverage$last_global_im), expected_last) ||
      !identical(as.integer(imp$m), requested_m) ||
      any(coverage$collection_status != "COLLECTED_BEFORE_IBIND")) {
    stop("MICE worker-event coverage is incomplete", call. = FALSE)
  }
  logged <- imp$loggedEvents
  observed <- integer(length(expected_workers))
  if (!is.null(logged) && nrow(logged)) {
    assert_columns(logged, c("worker_id", "local_im"),
                   "combined MICE worker events")
    if (any(!logged$worker_id %in% expected_workers)) {
      stop("Combined MICE events contain an unknown worker", call. = FALSE)
    }
    event_worker <- suppressWarnings(as.integer(logged$worker_id))
    local_im <- suppressWarnings(as.integer(logged$local_im))
    global_im <- suppressWarnings(as.integer(logged$im))
    worker_index <- match(event_worker, expected_workers)
    if (any(is.na(event_worker)) || any(is.na(local_im)) ||
        any(is.na(global_im)) || any(local_im < 0L) ||
        any(local_im > expected_worker_m[worker_index])) {
      stop("Combined MICE events contain invalid local/global im indices",
           call. = FALSE)
    }
    expected_global_im <- ifelse(
      local_im == 0L, 0L, expected_offsets[worker_index] + local_im
    )
    if (!identical(global_im, as.integer(expected_global_im))) {
      stop("Combined MICE global/local im mapping is invalid", call. = FALSE)
    }
    observed <- tabulate(event_worker,
                         nbins = length(expected_workers))
  }
  if (!identical(as.integer(coverage$event_rows), as.integer(observed))) {
    stop("Combined MICE events do not reconcile to worker coverage", call. = FALSE)
  }
  invisible(TRUE)
}

run_auditable_parallel_mice <- function(data, m, maxit, method,
                                        predictorMatrix, where, ignore,
                                        n.core, parallelseed) {
  cores <- min(max(1L, as.integer(n.core)), as.integer(m))
  worker_m <- rep(as.integer(m) %/% cores, cores)
  remainder <- as.integer(m) %% cores
  if (remainder > 0L) worker_m[seq_len(remainder)] <-
    worker_m[seq_len(remainder)] + 1L
  run_worker <- function(i) {
    configure_worker_runtime()
    mice::mice(
      data = data, m = worker_m[[i]], maxit = maxit, method = method,
      predictorMatrix = predictorMatrix, where = where, ignore = ignore,
      printFlag = FALSE, seed = NA
    )
  }
  if (cores == 1L) {
    set.seed(parallelseed)
    imps <- list(run_worker(1L))
  } else {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = cores)
    imps <- future.apply::future_lapply(
      seq_len(cores), run_worker, future.seed = parallelseed,
      future.packages = "mice"
    )
  }
  imp <- combine_mids_preserving_worker_events(
    imps, requested_m = m, requested_worker_m = worker_m
  )
  imp$parallelseed <- parallelseed
  enforce_mice_worker_event_coverage(imp, requested_m = m,
                                     requested_workers = cores,
                                     expected_worker_m = worker_m)
  imp
}

stage_imputation <- function() {
  require_nova_gate("covariate_imputation")
  frame_path <- latest_frame_path()
  frame <- readRDS(frame_path)
  x <- build_mice_input(frame)
  cfg <- configure_mice(x)
  run_id <- timestamp_id(if (identical(ANALYSIS_MODE, "formal")) {
    "N7_MICE_FORMAL"
  } else {
    "N6_MICE_SMOKE"
  })
  derived_run <- file.path(ANALYSIS_DIRS$derived, run_id)
  qc_run <- file.path(ANALYSIS_DIRS$qc, run_id)
  log_run <- file.path(ANALYSIS_DIRS$logs, run_id)
  dir.create(derived_run, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_run, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_run, recursive = TRUE, showWarnings = FALSE)
  sink(file.path(log_run, "mice.log"), split = TRUE)
  on.exit(sink(), add = TRUE)
  cat("run_id=", run_id, "\n", sep = "")
  cat("m=", MI_M, " maxit=", MI_MAXIT, " workers=", N_WORKERS, "\n", sep = "")

  missing_before <- data.frame(
    variable = IMPUTED_COVARIATES,
    n_missing_base = vapply(IMPUTED_COVARIATES, function(v) {
      sum(x$base_domain %in% TRUE & is.na(x[[v]]))
    }, integer(1)),
    pct_missing_base = vapply(IMPUTED_COVARIATES, function(v) {
      mean(is.na(x[[v]][x$base_domain %in% TRUE])) * 100
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
  write_csv_atomic(missing_before, file.path(qc_run, "missingness_before_mice.csv"))
  methods <- data.frame(variable = names(cfg$method), method = unname(cfg$method),
                        stringsAsFactors = FALSE)
  predictor_export <- data.frame(target = rownames(cfg$predictor), cfg$predictor,
                                 check.names = FALSE, stringsAsFactors = FALSE)
  write_csv_atomic(methods, file.path(qc_run, "mice_methods.csv"))
  write_csv_atomic(predictor_export, file.path(qc_run, "mice_predictor_matrix.csv"))
  if (any(!nzchar(cfg$method[IMPUTED_COVARIATES]))) {
    stop("At least one imputation target has an empty method")
  }

  cores <- min(N_WORKERS, MI_M)
  imp <- run_auditable_parallel_mice(
    data = x,
    m = MI_M,
    maxit = MI_MAXIT,
    method = cfg$method,
    predictorMatrix = cfg$predictor,
    where = cfg$where,
    ignore = cfg$ignore,
    n.core = cores,
    parallelseed = KN_SEED
  )
  attr(imp, "knhanes_frame_sha256") <- sha256_file(frame_path)
  attr(imp, "knhanes_frame_path") <- frame_path
  attr(imp, "knhanes_ruleset_version") <- NOVA_RULESET_VERSION
  attr(imp, "knhanes_analysis_scope") <- ANALYSIS_SCOPE
  pending_path <- file.path(derived_run, "covariate_imputations.mids.pending.rds")
  saveRDS(imp, pending_path, compress = "gzip")
  validation <- validate_imputations(imp, x$base_domain %in% TRUE)
  write_csv_atomic(validation, file.path(qc_run, "mice_completion_check.csv"))
  worker_coverage <- attr(imp, "knhanes_worker_event_coverage")
  enforce_mice_worker_event_coverage(
    imp, requested_m = MI_M, requested_workers = cores
  )
  write_csv_atomic(
    worker_coverage, file.path(qc_run, "mice_worker_event_coverage.csv")
  )
  logged <- imp$loggedEvents
  if (is.null(logged)) logged <- data.frame(message = character())
  write_csv_atomic(logged, file.path(qc_run, "mice_logged_events.csv"))
  logged_audit <- audit_mice_logged_events(logged, maxit = MI_MAXIT)
  write_csv_atomic(logged_audit, file.path(qc_run, "mice_logged_event_audit.csv"))
  enforce_mice_logged_event_gate(logged_audit)
  hashes <- imputation_hashes(imp)
  chains <- chain_diagnostics(imp)
  convergence <- mice_convergence_diagnostics(imp)
  target_chains <- chains$variable %in% IMPUTED_COVARIATES
  if (nrow(chains) && any(target_chains & chains$finite_fraction < 1)) {
    stop("MICE chain diagnostics contain non-finite values")
  }
  write_csv_atomic(hashes, file.path(qc_run, "mice_imputation_hashes.csv"))
  write_csv_atomic(chains, file.path(qc_run, "mice_chain_diagnostics.csv"))
  write_csv_atomic(convergence$all, file.path(qc_run, "mice_convergence_all_iterations.csv"))
  write_csv_atomic(convergence$last, file.path(qc_run, "mice_convergence_last_iteration.csv"))
  enforce_mice_convergence_gate(convergence$last, ANALYSIS_MODE)

  imp_path <- file.path(derived_run, "covariate_imputations.mids.rds")
  if (file.exists(imp_path)) file.remove(imp_path)
  if (!file.rename(pending_path, imp_path)) stop("Could not atomically commit mids object")
  manifest_paths <- c(
    frame_path, imp_path, file.path(qc_run, "missingness_before_mice.csv"),
    file.path(qc_run, "mice_completion_check.csv"),
    file.path(qc_run, "mice_worker_event_coverage.csv"),
    file.path(qc_run, "mice_logged_events.csv"),
    file.path(qc_run, "mice_logged_event_audit.csv"),
    file.path(qc_run, "mice_methods.csv"),
    file.path(qc_run, "mice_predictor_matrix.csv"),
    file.path(qc_run, "mice_imputation_hashes.csv"),
    file.path(qc_run, "mice_chain_diagnostics.csv"),
    file.path(qc_run, "mice_convergence_all_iterations.csv"),
    file.path(qc_run, "mice_convergence_last_iteration.csv")
  )
  manifest <- data.frame(
    path = manifest_paths, bytes = file.info(manifest_paths)$size,
    sha256 = vapply(manifest_paths, sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
  write_csv_atomic(manifest, file.path(log_run, "output_manifest.csv"))
  writeLines(capture.output(sessionInfo()), file.path(log_run, "sessionInfo.txt"))
  write_lines_atomic(as.character(imp$parallelseed %||% KN_SEED),
                     file.path(log_run, "parallelseed.txt"))
  commit_latest_marker(run_id, file.path(ANALYSIS_DIRS$derived, "LATEST_MICE_RUN_ID.txt"),
                       c(manifest_paths, file.path(log_run, "output_manifest.csv"),
                         file.path(log_run, "parallelseed.txt")))
  cat("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")
  invisible(list(run_id = run_id, mids = imp_path))
}

if (sys.nframe() == 0L) stage_imputation()
)--------------------"

upf_source[["KN/source_snapshot/06_main_models.R"]] <- r"--------------------(# Fit the prespecified four-outcome M1-M4 continuous, quartile and trend models.

.this_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "KN/source_snapshot/06_main_models.R"), error = function(e) NULL)
if (is.null(.this_file)) {
  .file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.file_arg)) .this_file <- sub("^--file=", "", .file_arg[[1]])
}
if (is.null(.this_file)) .this_file <- "06_main_models.R"
.code_dir <- dirname(normalizePath(.this_file, mustWork = FALSE))
source(file.path(.code_dir, "00_model_functions.R"), local = FALSE)

stage_main_models <- function() {
  require_nova_gate("main_models")
  inputs <- load_model_inputs()
  if (inputs$imp$m != MI_M) {
    stop("Configured MI_M does not match stored mids object", call. = FALSE)
  }
  grid <- expand.grid(
    outcome = names(OUTCOMES), model = names(MODEL_COVARIATES),
    exposure_form = c("continuous", "quartile", "trend"),
    stringsAsFactors = FALSE
  )
  run_id <- timestamp_id(if (identical(ANALYSIS_MODE, "formal")) {
    "N8_MODELS_FORMAL"
  } else {
    "N6_MODELS_SMOKE"
  })
  result_run <- file.path(ANALYSIS_DIRS$results, run_id)
  qc_run <- file.path(ANALYSIS_DIRS$qc, run_id)
  log_run <- file.path(ANALYSIS_DIRS$logs, run_id)
  dir.create(result_run, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_run, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_run, recursive = TRUE, showWarnings = FALSE)
  sink(file.path(log_run, "main_models.log"), split = TRUE)
  on.exit(sink(), add = TRUE)
  cat("run_id=", run_id, " m=", MI_M, " workers=", N_WORKERS, "\n", sep = "")

  workers <- min(N_WORKERS, length(OUTCOMES))
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = workers)
  fitted_by_outcome <- future.apply::future_lapply(
    names(OUTCOMES),
    function(outcome_name) {
      configure_worker_runtime()
      local_inputs <- load_model_inputs()
      outcome_grid <- expand.grid(
        outcome = outcome_name, model = names(MODEL_COVARIATES),
        exposure_form = c("continuous", "quartile", "trend"),
        stringsAsFactors = FALSE
      )
      lapply(seq_len(nrow(outcome_grid)), function(i) {
        fit_mi_model_stream(
          frame = local_inputs$frame, imp = local_inputs$imp,
          outcome_name = outcome_grid$outcome[[i]], model_name = outcome_grid$model[[i]],
          exposure_form = outcome_grid$exposure_form[[i]]
        )
      })
    },
    future.seed = KN_SEED,
    future.packages = c("survey", "mitools", "mice", "dplyr")
  )
  fitted <- unlist(fitted_by_outcome, recursive = FALSE)
  results_all <- dplyr::bind_rows(lapply(fitted, `[[`, "result"))
  results_exposure <- exposure_rows(results_all)
  metadata <- dplyr::bind_rows(lapply(fitted, `[[`, "metadata"))
  formulas <- dplyr::bind_rows(lapply(fitted, function(x) {
    data.frame(
      outcome = x$metadata$outcome[[1]], model = x$metadata$model[[1]],
      exposure_form = x$metadata$exposure_form[[1]],
      formula = paste(deparse(x$formula), collapse = " "), stringsAsFactors = FALSE
    )
  }))
  results_all$analysis_scope <- ANALYSIS_SCOPE
  results_exposure$analysis_scope <- ANALYSIS_SCOPE
  metadata$analysis_scope <- ANALYSIS_SCOPE
  formulas$analysis_scope <- ANALYSIS_SCOPE
  expected_rows <- length(OUTCOMES) * length(MODEL_COVARIATES) *
    (1L + 3L + 1L)
  if (any(!metadata$success) || nrow(results_exposure) != expected_rows ||
      any(results_exposure$status != "OK") ||
      any(!is.finite(unlist(results_exposure[c(
        "estimate", "std_error_link", "conf_low", "conf_high", "p_value", "df"
      )])))) {
    write_csv_atomic(metadata, file.path(qc_run, "model_fit_metadata.csv"))
    stop("At least one main model failed; partial pooled results are not accepted", call. = FALSE)
  }
  expected_successes <- nrow(grid) * MI_M
  if (sum(metadata$success) != expected_successes) {
    stop("Unexpected successful-fit count", call. = FALSE)
  }

  frame_with_q <- derive_fixed_exposure_categories(inputs$frame)
  qcuts <- data.frame(
    statistic = c("Q1_Q2_cut", "Q2_Q3_cut", "Q3_Q4_cut", paste0("median_Q", 1:4)),
    upf_grams_pct = c(attr(frame_with_q, "upf_quartile_cuts"),
                      attr(frame_with_q, "upf_quartile_medians")),
    stringsAsFactors = FALSE
  )
  qcuts$analysis_scope <- ANALYSIS_SCOPE
  result_paths <- c(
    all_terms = file.path(result_run, "main_models_all_terms.csv"),
    exposure = file.path(result_run, "main_models_exposure_terms.csv"),
    metadata = file.path(qc_run, "model_fit_metadata.csv"),
    formulas = file.path(qc_run, "model_formulas.csv"),
    qcuts = file.path(qc_run, "upf_weighted_quartiles.csv"),
    input_lineage = file.path(qc_run, "main_model_input_lineage.csv")
  )
  write_csv_atomic(results_all, result_paths[["all_terms"]])
  write_csv_atomic(results_exposure, result_paths[["exposure"]])
  write_csv_atomic(metadata, result_paths[["metadata"]])
  write_csv_atomic(formulas, result_paths[["formulas"]])
  write_csv_atomic(qcuts, result_paths[["qcuts"]])
  input_lineage <- data.frame(
    role = c("analysis_frame", "mids", "code_review_marker", "nova_gate_marker",
             "active_nova_map"),
    path = c(inputs$frame_path, inputs$mice_path, CODE_REVIEW_MARKER, NOVA_GATE_MARKER,
             ACTIVE_NOVA_MAP),
    sha256 = vapply(c(inputs$frame_path, inputs$mice_path, CODE_REVIEW_MARKER,
                      NOVA_GATE_MARKER, ACTIVE_NOVA_MAP), sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
  write_csv_atomic(input_lineage, result_paths[["input_lineage"]])
  manifest_paths <- unname(result_paths)
  manifest <- data.frame(
    path = manifest_paths, bytes = file.info(manifest_paths)$size,
    sha256 = vapply(manifest_paths, sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
  write_csv_atomic(manifest, file.path(log_run, "output_manifest.csv"))
  writeLines(capture.output(sessionInfo()), file.path(log_run, "sessionInfo.txt"))
  commit_latest_marker(run_id, file.path(ANALYSIS_DIRS$results, "LATEST_MAIN_MODEL_RUN_ID.txt"),
                       c(manifest_paths, file.path(log_run, "output_manifest.csv")))
  cat("successful_fits=", sum(metadata$success), " expected=", expected_successes, "\n", sep = "")
  cat("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")
  invisible(list(run_id = run_id, results = result_paths[["exposure"]]))
}

if (sys.nframe() == 0L) stage_main_models()
)--------------------"

upf_source[["KN/source_snapshot/07_extended_models.R"]] <- r"--------------------(# Prespecified M3 RCS, subgroup interactions and sensitivity analyses.

.this_file <- tryCatch(file.path(Sys.getenv("UPF_RUN_ROOT"), "KN/source_snapshot/07_extended_models.R"), error = function(e) NULL)
if (is.null(.this_file)) {
  .file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(.file_arg)) .this_file <- sub("^--file=", "", .file_arg[[1]])
}
if (is.null(.this_file)) .this_file <- "07_extended_models.R"
.code_dir <- dirname(normalizePath(.this_file, mustWork = FALSE))
source(file.path(.code_dir, "00_model_functions.R"), local = FALSE)

rcs_nonlinear_basis <- function(x, knots) {
  if (length(knots) != 3L || any(diff(knots) <= 0)) stop("RCS requires 3 ordered knots")
  k1 <- knots[[1]]; k2 <- knots[[2]]; k3 <- knots[[3]]
  (pmax(x - k1, 0)^3 - pmax(x - k2, 0)^3 * (k3 - k1) / (k3 - k2) +
     pmax(x - k3, 0)^3 * (k2 - k1) / (k3 - k2)) / (k3 - k1)^2
}

finite_complete_df <- function(fits) {
  out <- min(vapply(fits, `[[`, numeric(1), "residual_df"))
  if (!is.finite(out) || out <= 0) stop("Invalid complete-data survey residual df")
  out
}

mi_pool_components <- function(fits) {
  if (!length(fits) || any(vapply(fits, is.null, logical(1)))) return(NULL)
  df_complete <- finite_complete_df(fits)
  q <- lapply(fits, `[[`, "coefficients")
  u <- lapply(fits, `[[`, "variance")
  if (length(unique(vapply(q, function(z) paste(names(z), collapse = "|"),
                           character(1)))) != 1L) {
    stop("Coefficient sets differ across imputations")
  }
  pooled <- if (length(fits) == 1L) {
    list(coefficients = q[[1]], variance = u[[1]],
         df = rep(df_complete, length(q[[1]])), missinfo = rep(0, length(q[[1]])))
  } else {
    mitools::MIcombine(results = q, variances = u, df.complete = df_complete)
  }
  list(beta = pooled$coefficients, variance = pooled$variance,
       df = pooled$df, fmi = pooled$missinfo, df_complete = df_complete)
}

pool_scalar_contrast <- function(fits, contrast) {
  q <- lapply(fits, function(z) as.numeric(sum(contrast * z$coefficients)))
  u <- lapply(fits, function(z) as.numeric(t(contrast) %*% z$variance %*% contrast))
  df_complete <- finite_complete_df(fits)
  pooled <- if (length(fits) == 1L) {
    list(coefficients = q[[1]], variance = matrix(u[[1]], 1, 1),
         df = df_complete, missinfo = 0)
  } else {
    mitools::MIcombine(results = q, variances = u, df.complete = df_complete)
  }
  estimate <- as.numeric(pooled$coefficients)
  se <- sqrt(as.numeric(pooled$variance))
  df <- as.numeric(pooled$df)[[1]]
  fmi <- as.numeric(pooled$missinfo %||% NA_real_)[[1]]
  crit <- stats::qt(0.975, df)
  data.frame(
    estimate_link = estimate, std_error_link = se,
    conf_low_link = estimate - crit * se, conf_high_link = estimate + crit * se,
    p_value = 2 * stats::pt(abs(estimate / se), df = df, lower.tail = FALSE),
    df = df, fmi = fmi, complete_data_df = df_complete,
    stringsAsFactors = FALSE
  )
}

mi_d1_test <- function(fits, terms) {
  if (!length(terms)) stop("D1 test has no terms")
  aliases <- paste0("b", seq_along(terms))
  qhat <- lapply(fits, function(z) {
    q <- z$coefficients[terms]
    if (any(!is.finite(q))) stop("Non-finite D1 coefficient")
    stats::setNames(as.numeric(q), aliases)
  })
  uhat <- lapply(fits, function(z) {
    v <- z$variance[terms, terms, drop = FALSE]
    dimnames(v) <- list(aliases, aliases)
    v
  })
  df_complete <- finite_complete_df(fits)
  smoke_df_fallback <- identical(MI_M, 2L) && length(fits) == 2L
  test <- mitml::testConstraints(
    qhat = qhat, uhat = uhat, constraints = aliases, method = "D1",
    df.com = if (smoke_df_fallback) NULL else df_complete
  )$test[1, ]
  test_values <- unname(test[c("F.value", "df1", "df2", "P(>F)", "RIV")])
  if (any(!is.finite(test_values)) || test[["df1"]] <= 0 ||
      test[["df2"]] <= 0 || test[["P(>F)"]] < 0 || test[["P(>F)"]] > 1) {
    stop("D1 test returned a non-finite or invalid result")
  }
  data.frame(
    statistic = unname(test[["F.value"]]), df1 = unname(test[["df1"]]),
    df2 = unname(test[["df2"]]), p_value = unname(test[["P(>F)"]]),
    relative_increase_variance = unname(test[["RIV"]]),
    complete_data_df = df_complete,
    d1_df_method = if (smoke_df_fallback) {
      "smoke_m2_unadjusted"
    } else {
      "finite_complete_adjusted"
    },
    stringsAsFactors = FALSE
  )
}

fit_rcs_outcome_stream <- function(frame, imp, outcome_name, knots_pp, display_pp) {
  spec <- OUTCOMES[[outcome_name]]
  covars <- if (outcome_name == "DKD") MODEL_COVARIATES_DKD$M3 else MODEL_COVARIATES$M3
  knots <- knots_pp / UPF_UNIT_PP
  frame <- derive_fixed_exposure_categories(frame)
  form <- build_model_formula(spec$variable, c("upf_per10", "upf_rcs_nonlinear"), covars)
  fits <- lapply(seq_len(imp$m), function(i) {
    d <- completed_analysis_dataset(frame, imp, i)
    d$upf_rcs_nonlinear <- rcs_nonlinear_basis(d$upf_per10, knots)
    tryCatch(fit_svy_formula(d, form, spec), error = function(e) NULL)
  })
  if (any(vapply(fits, is.null, logical(1)))) stop("RCS failed in at least one imputation")
  invisible(lapply(fits, validate_compact_fit, family = spec$family,
                   label = paste0(outcome_name, " RCS")))
  term_names <- names(fits[[1]]$coefficients)
  grid_pp <- seq(display_pp[[1]], display_pp[[2]], length.out = 121)
  ref_x <- knots[[2]]
  curve <- dplyr::bind_rows(lapply(seq_along(grid_pp), function(j) {
    grid_x <- grid_pp[[j]] / UPF_UNIT_PP
    L <- stats::setNames(rep(0, length(term_names)), term_names)
    L[["upf_per10"]] <- grid_x - ref_x
    L[["upf_rcs_nonlinear"]] <-
      rcs_nonlinear_basis(grid_x, knots) - rcs_nonlinear_basis(ref_x, knots)
    z <- pool_scalar_contrast(fits, L)
    if (spec$family == "binomial") {
      z$estimate <- exp(z$estimate_link)
      z$conf_low <- exp(z$conf_low_link)
      z$conf_high <- exp(z$conf_high_link)
      z$effect_scale <- "OR"
    } else {
      z$estimate <- z$estimate_link
      z$conf_low <- z$conf_low_link
      z$conf_high <- z$conf_high_link
      z$effect_scale <- "difference"
    }
    z$outcome <- outcome_name; z$upf_grams_pct <- grid_pp[[j]]
    z$reference_pp <- knots_pp[[2]]; z$knot_10 <- knots_pp[[1]]
    z$knot_50 <- knots_pp[[2]]; z$knot_90 <- knots_pp[[3]]; z$status <- "OK"
    z
  }))
  nonlinear <- mi_d1_test(fits, "upf_rcs_nonlinear")
  nonlinear$test <- "nonlinearity"
  overall <- mi_d1_test(fits, c("upf_per10", "upf_rcs_nonlinear"))
  overall$test <- "overall_spline"
  tests <- dplyr::bind_rows(nonlinear, overall)
  tests$outcome <- outcome_name; tests$status <- "OK"
  metadata <- data.frame(
    outcome = outcome_name, successful_imputations = length(fits),
    required_imputations = imp$m,
    min_design_df = min(vapply(fits, `[[`, numeric(1), "design_df")),
    min_residual_df = finite_complete_df(fits), stringsAsFactors = FALSE
  )
  list(curve = curve, test = tests, metadata = metadata)
}

subgroup_variable <- function(d, subgroup) {
  if (subgroup == "age_group") {
    factor(ifelse(d$age < 45, "20-44", ifelse(d$age < 65, "45-64", "65+")),
           levels = SUBGROUPS$age_group)
  } else {
    factor(d[[subgroup]], levels = SUBGROUPS[[subgroup]])
  }
}

subgroup_support_table <- function(d, form, spec, outcome_name, subgroup,
                                   imputation) {
  model_complete <- stats::complete.cases(
    stats::model.frame(form, data = d, na.action = na.pass)
  )
  domain <- d[[spec$domain]] %in% TRUE & model_complete
  cluster_key <- interaction(d$strata_uid, d$psu_uid, drop = TRUE)
  rows <- dplyr::bind_rows(lapply(SUBGROUPS[[subgroup]], function(level) {
    hit <- domain & d$.subgroup == level
    n <- sum(hit)
    events <- if (identical(spec$family, "binomial")) {
      sum(d[[spec$variable]][hit] == 1, na.rm = TRUE)
    } else NA_integer_
    n_psu <- length(unique(cluster_key[hit]))
    n_strata <- length(unique(d$strata_uid[hit]))
    data.frame(
      outcome = outcome_name, subgroup = subgroup, level = level,
      imputation = imputation, n = n, events = events,
      non_events = if (identical(spec$family, "binomial")) n - events else NA_integer_,
      n_psu = n_psu, n_strata = n_strata,
      support_df = n_psu - n_strata, stringsAsFactors = FALSE
    )
  }))
  invalid <- rows$n <= 0L | rows$n_psu < 2L | rows$support_df <= 0L
  if (identical(spec$family, "binomial")) {
    invalid <- invalid | rows$events < MIN_BINARY_EVENTS |
      rows$non_events < MIN_BINARY_NONEVENTS
  }
  if (any(invalid)) {
    bad <- paste(rows$level[invalid], "n", rows$n[invalid], "events", rows$events[invalid],
                 "non_events", rows$non_events[invalid], "psu", rows$n_psu[invalid],
                 "strata", rows$n_strata[invalid], sep = "=")
    stop(outcome_name, "/", subgroup,
         " subgroup support gate failed at imputation ", imputation, ": ",
         paste(bad, collapse = " | "), call. = FALSE)
  }
  rows
}

fit_subgroup_outcome_stream <- function(frame, imp, outcome_name, subgroup) {
  spec <- OUTCOMES[[outcome_name]]
  covars <- if (outcome_name == "DKD") MODEL_COVARIATES_DKD$M3 else MODEL_COVARIATES$M3
  if (subgroup == "age_group") covars <- setdiff(covars, "age")
  if (subgroup %in% covars) covars <- setdiff(covars, subgroup)
  frame <- derive_fixed_exposure_categories(frame)
  form <- build_model_formula(spec$variable, "upf_per10*.subgroup", covars)
  bundles <- lapply(seq_len(imp$m), function(i) {
    d <- completed_analysis_dataset(frame, imp, i)
    d$.subgroup <- subgroup_variable(d, subgroup)
    support <- subgroup_support_table(d, form, spec, outcome_name, subgroup, i)
    fit <- tryCatch(fit_svy_formula(d, form, spec), error = function(e) {
      stop(outcome_name, "/", subgroup, " fit failed at imputation ", i,
           ": ", conditionMessage(e), call. = FALSE)
    })
    validate_compact_fit(fit, spec$family,
                         paste0(outcome_name, "/", subgroup, "/", i))
    list(fit = fit, support = support)
  })
  fits <- lapply(bundles, `[[`, "fit")
  support <- dplyr::bind_rows(lapply(bundles, `[[`, "support"))
  term_names <- names(fits[[1]]$coefficients)
  levels_sg <- SUBGROUPS[[subgroup]]
  effects <- dplyr::bind_rows(lapply(seq_along(levels_sg), function(j) {
    L <- stats::setNames(rep(0, length(term_names)), term_names)
    L[["upf_per10"]] <- 1
    if (j > 1L) {
      candidates <- c(paste0("upf_per10:.subgroup", levels_sg[[j]]),
                      paste0(".subgroup", levels_sg[[j]], ":upf_per10"))
      hit <- intersect(candidates, term_names)
      if (length(hit) != 1L) stop("Missing subgroup interaction coefficient")
      L[[hit]] <- 1
    }
    z <- pool_scalar_contrast(fits, L)
    if (spec$family == "binomial") {
      z$estimate <- exp(z$estimate_link); z$conf_low <- exp(z$conf_low_link)
      z$conf_high <- exp(z$conf_high_link); z$effect_scale <- "OR"
    } else {
      z$estimate <- z$estimate_link; z$conf_low <- z$conf_low_link
      z$conf_high <- z$conf_high_link; z$effect_scale <- "beta"
    }
    z$outcome <- outcome_name; z$subgroup <- subgroup; z$level <- levels_sg[[j]]
    z$status <- "OK"; z
  }))
  interaction_terms <- grep("upf_per10:|:upf_per10", term_names, value = TRUE)
  interaction <- mi_d1_test(fits, interaction_terms)
  names(interaction)[names(interaction) == "p_value"] <- "p_interaction"
  interaction$outcome <- outcome_name; interaction$subgroup <- subgroup
  interaction$interaction_df <- length(interaction_terms); interaction$status <- "OK"
  list(effects = effects, interaction = interaction, support = support)
}

apply_alternate_lab_outcome <- function(d, outcome_name) {
  if (outcome_name == "log_uacr") d$log_uacr <- d$log_uacr_alt_lod
  if (outcome_name == "CKD") d$ckd <- d$ckd_alt_lod
  if (outcome_name == "DKD") d$dkd <- d$dkd_alt_lod
  d
}

fit_sensitivities_outcome_stream <- function(frame, imp, outcome_name) {
  spec <- OUTCOMES[[outcome_name]]
  original <- derive_fixed_exposure_categories(frame)
  cc <- fit_mi_model(list(original), outcome_name, "M3", "continuous", complete_case = TRUE)
  no_ssb <- fit_mi_model_stream(
    frame, imp, outcome_name, "M3", "continuous", transform_data = function(d) {
      d$upf_per10 <- d$upf_non_ssb_pct / UPF_UNIT_PP
      d[[spec$domain]] <- d[[spec$domain]] %in% TRUE & is.finite(d$upf_per10)
      d
    }
  )
  energy_fit <- fit_mi_model_stream(
    frame, imp, outcome_name, "M3", "continuous", transform_data = function(d) {
      d[[spec$domain]] <- d[[spec$domain]] %in% TRUE & is.finite(d$energy_kcal) &
        d$energy_kcal >= ENERGY_SENSITIVITY_RANGE[[1]] &
        d$energy_kcal <= ENERGY_SENSITIVITY_RANGE[[2]]
      d
    }
  )
  fits <- list(complete_case = cc, remove_SSB_recalculate = no_ssb,
               energy_500_5000 = energy_fit)
  if (outcome_name %in% c("log_uacr", "CKD", "DKD")) {
    fits$lab_LOD_half <- fit_mi_model_stream(
      frame, imp, outcome_name, "M3", "continuous",
      transform_data = function(d) apply_alternate_lab_outcome(d, outcome_name)
    )
  }
  rows <- dplyr::bind_rows(lapply(names(fits), function(nm) {
    transform(exposure_rows(fits[[nm]]$result), sensitivity = nm)
  }))
  metadata <- dplyr::bind_rows(lapply(names(fits), function(nm) {
    transform(fits[[nm]]$metadata, sensitivity = nm)
  }))
  list(results = rows, metadata = metadata)
}

stage_extended_models <- function() {
  require_nova_gate("extended_models")
  inputs <- load_model_inputs()
  if (inputs$imp$m != MI_M) stop("Configured MI_M does not match stored mids object")
  frame_q <- derive_fixed_exposure_categories(inputs$frame)
  base <- frame_q$base_domain %in% TRUE
  knots_pp <- weighted_quantile(frame_q$upf_grams_pct[base],
                                frame_q$analysis_weight[base], RCS_KNOT_PROBS)
  display_pp <- weighted_quantile(frame_q$upf_grams_pct[base],
                                  frame_q$analysis_weight[base], RCS_DISPLAY_PROBS)
  run_id <- timestamp_id(if (identical(ANALYSIS_MODE, "formal")) {
    "N9_EXTENDED_FORMAL"
  } else {
    "N6_EXTENDED_SMOKE"
  })
  result_run <- file.path(ANALYSIS_DIRS$results, run_id)
  qc_run <- file.path(ANALYSIS_DIRS$qc, run_id)
  log_run <- file.path(ANALYSIS_DIRS$logs, run_id)
  dir.create(result_run, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_run, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_run, recursive = TRUE, showWarnings = FALSE)
  sink(file.path(log_run, "extended_models.log"), split = TRUE)
  on.exit(sink(), add = TRUE)
  cat("run_id=", run_id, " m=", MI_M, " workers=", N_WORKERS, "\n", sep = "")

  old_plan <- future::plan(); on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = min(N_WORKERS, length(OUTCOMES)))
  by_outcome <- future.apply::future_lapply(
    names(OUTCOMES),
    function(outcome_name) {
      configure_worker_runtime()
      local_inputs <- load_model_inputs()
      subgroup_names <- names(SUBGROUPS)
      if (outcome_name == "DKD") subgroup_names <- setdiff(subgroup_names, "diabetes")
      list(
        rcs = fit_rcs_outcome_stream(local_inputs$frame, local_inputs$imp,
                                     outcome_name, knots_pp, display_pp),
        subgroups = lapply(subgroup_names, function(sg) {
          fit_subgroup_outcome_stream(local_inputs$frame, local_inputs$imp, outcome_name, sg)
        }),
        sensitivities = fit_sensitivities_outcome_stream(
          local_inputs$frame, local_inputs$imp, outcome_name
        )
      )
    },
    future.seed = KN_SEED,
    future.packages = c("survey", "mitools", "mitml", "mice", "dplyr")
  )
  rcs_curve <- dplyr::bind_rows(lapply(by_outcome, function(x) x$rcs$curve))
  rcs_test <- dplyr::bind_rows(lapply(by_outcome, function(x) x$rcs$test))
  rcs_meta <- dplyr::bind_rows(lapply(by_outcome, function(x) x$rcs$metadata))
  subgroup_effects <- dplyr::bind_rows(lapply(by_outcome, function(x) {
    dplyr::bind_rows(lapply(x$subgroups, `[[`, "effects"))
  }))
  interactions <- dplyr::bind_rows(lapply(by_outcome, function(x) {
    dplyr::bind_rows(lapply(x$subgroups, `[[`, "interaction"))
  }))
  subgroup_support <- dplyr::bind_rows(lapply(by_outcome, function(x) {
    dplyr::bind_rows(lapply(x$subgroups, `[[`, "support"))
  }))
  interactions$p_interaction_fdr <- p.adjust(interactions$p_interaction, method = "BH")
  sensitivity <- dplyr::bind_rows(lapply(by_outcome, function(x) x$sensitivities$results))
  sensitivity_meta <- dplyr::bind_rows(lapply(by_outcome, function(x) {
    x$sensitivities$metadata
  }))
  for (object_name in c("rcs_curve", "rcs_test", "rcs_meta", "subgroup_effects",
                        "interactions", "subgroup_support", "sensitivity",
                        "sensitivity_meta")) {
    object <- get(object_name)
    object$analysis_scope <- ANALYSIS_SCOPE
    assign(object_name, object)
  }
  required_numeric <- c("estimate", "conf_low", "conf_high", "p_value", "df")
  finite_table <- function(x, cols) {
    all(cols %in% names(x)) && all(is.finite(as.matrix(x[, cols, drop = FALSE])))
  }
  interval_ok <- function(x) {
    all(x$conf_low <= x$estimate & x$estimate <= x$conf_high)
  }
  n_lab_sens <- length(intersect(names(OUTCOMES), c("log_uacr", "CKD", "DKD")))
  expected_subgroup_rows <- (length(OUTCOMES) - 1L) * sum(lengths(SUBGROUPS)) +
    sum(lengths(SUBGROUPS[names(SUBGROUPS) != "diabetes"]))
  expected_subgroup_support <- expected_subgroup_rows * MI_M
  expected_interactions <- (length(OUTCOMES) - 1L) * length(SUBGROUPS) +
    length(SUBGROUPS[names(SUBGROUPS) != "diabetes"])
  expected_sensitivity_rows <- length(OUTCOMES) * 3L + n_lab_sens
  expected_sensitivity_meta <- length(OUTCOMES) +
    (2L * length(OUTCOMES) + n_lab_sens) * MI_M
  if (nrow(rcs_curve) != length(OUTCOMES) * 121L ||
      nrow(rcs_test) != length(OUTCOMES) * 2L ||
      nrow(rcs_meta) != length(OUTCOMES) ||
      nrow(subgroup_effects) != expected_subgroup_rows ||
      nrow(subgroup_support) != expected_subgroup_support ||
      nrow(interactions) != expected_interactions ||
      nrow(sensitivity) != expected_sensitivity_rows ||
      nrow(sensitivity_meta) != expected_sensitivity_meta ||
      any(rcs_curve$status != "OK") || any(rcs_test$status != "OK") ||
      any(interactions$status != "OK") || any(subgroup_effects$status != "OK") ||
      any(subgroup_support$n <= 0L | subgroup_support$n_psu < 2L |
            subgroup_support$support_df <= 0L) ||
      any(is.finite(subgroup_support$events) &
            (subgroup_support$events < MIN_BINARY_EVENTS |
               subgroup_support$non_events < MIN_BINARY_NONEVENTS)) ||
      any(sensitivity$status != "OK") || any(!sensitivity_meta$success) ||
      any(rcs_meta$successful_imputations != MI_M) ||
      any(rcs_meta$required_imputations != MI_M) ||
      !finite_table(rcs_curve, c("estimate", "conf_low", "conf_high", "p_value", "df")) ||
      !finite_table(rcs_test, c("statistic", "df1", "df2", "p_value",
                               "complete_data_df")) ||
      !finite_table(subgroup_effects, required_numeric) ||
      !finite_table(interactions, c("statistic", "df1", "df2", "p_interaction",
                                    "p_interaction_fdr", "complete_data_df")) ||
      !finite_table(sensitivity, required_numeric) ||
      !interval_ok(rcs_curve) || !interval_ok(subgroup_effects) ||
      !interval_ok(sensitivity) || any(rcs_test$p_value < 0 | rcs_test$p_value > 1) ||
      any(interactions$p_interaction < 0 | interactions$p_interaction > 1) ||
      any(interactions$p_interaction_fdr < 0 | interactions$p_interaction_fdr > 1)) {
    stop("Extended-model hard gate failed")
  }
  paths <- c(
    rcs_curve = file.path(result_run, "rcs_curves.csv"),
    rcs_test = file.path(result_run, "rcs_tests_D1.csv"),
    subgroup_effects = file.path(result_run, "subgroup_effects.csv"),
    interactions = file.path(result_run, "subgroup_interactions_bh.csv"),
    sensitivity = file.path(result_run, "sensitivity_results.csv"),
    rcs_meta = file.path(qc_run, "rcs_metadata.csv"),
    subgroup_support = file.path(qc_run, "subgroup_support.csv"),
    sensitivity_meta = file.path(qc_run, "sensitivity_fit_metadata.csv"),
    input_lineage = file.path(qc_run, "extended_model_input_lineage.csv")
  )
  write_csv_atomic(rcs_curve, paths[["rcs_curve"]])
  write_csv_atomic(rcs_test, paths[["rcs_test"]])
  write_csv_atomic(subgroup_effects, paths[["subgroup_effects"]])
  write_csv_atomic(interactions, paths[["interactions"]])
  write_csv_atomic(sensitivity, paths[["sensitivity"]])
  write_csv_atomic(rcs_meta, paths[["rcs_meta"]])
  write_csv_atomic(subgroup_support, paths[["subgroup_support"]])
  write_csv_atomic(sensitivity_meta, paths[["sensitivity_meta"]])
  input_lineage <- data.frame(
    role = c("analysis_frame", "mids", "code_review_marker", "nova_gate_marker",
             "active_nova_map"),
    path = c(inputs$frame_path, inputs$mice_path, CODE_REVIEW_MARKER, NOVA_GATE_MARKER,
             ACTIVE_NOVA_MAP),
    sha256 = vapply(c(inputs$frame_path, inputs$mice_path, CODE_REVIEW_MARKER,
                      NOVA_GATE_MARKER, ACTIVE_NOVA_MAP), sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
  write_csv_atomic(input_lineage, paths[["input_lineage"]])
  manifest <- data.frame(
    path = unname(paths), bytes = file.info(unname(paths))$size,
    sha256 = vapply(unname(paths), sha256_file, character(1)), stringsAsFactors = FALSE
  )
  manifest_path <- file.path(log_run, "output_manifest.csv")
  write_csv_atomic(manifest, manifest_path)
  writeLines(capture.output(sessionInfo()), file.path(log_run, "sessionInfo.txt"))
  commit_latest_marker(run_id, file.path(ANALYSIS_DIRS$results, "LATEST_EXTENDED_MODEL_RUN_ID.txt"),
                       c(unname(paths), manifest_path))
  cat("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")
  invisible(list(run_id = run_id, paths = paths))
}

if (sys.nframe() == 0L) stage_extended_models()
)--------------------"

upf_source[["KN/01_code/04_run_models.R"]] <- r"--------------------(Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',KNHANES_N_CORES='14',KNHANES_MI_M='50',KNHANES_MI_MAXIT='20')
base <- '/storage/home/tmu2301/nhanes分析'
root <- file.path(Sys.getenv("UPF_RUN_ROOT"),"KN")
kn <- file.path(Sys.getenv("UPF_RUN_ROOT"),"KNHANES")
.libPaths(c(file.path(base,'11_UPF六结局_GA与跨库整合_20260903/R_library'),.libPaths()))
source(file.path(root,'source_snapshot/07_extended_models.R'))
source(file.path(root,'reuse_sources/ga_shared.R'))
args <- commandArgs(trailingOnly=TRUE)
mode <- 'candidate'
mipath <- file.path(root,'02_derived/candidate_mids.rds')
framepath <- file.path(root,'02_derived/analysis_frame_fasting_corrected.rds')
imp <- readRDS(mipath); frame <- readRDS(framepath)
MI_M <- as.integer(imp$m)
stopifnot(identical(attr(imp,'knhanes_frame_sha256'),sha256_file(framepath)),identical(frame$person_uid,imp$data$person_uid))
if(mode=='candidate') stopifnot(imp$m==50,file.exists(file.path(root,'03_qc/candidate/MICE_COMPLETED.txt')))
output <- file.path(root,'04_results',mode);dir.create(output,recursive=TRUE,showWarnings=FALSE)
tasksdir <- file.path(output,'tasks');dir.create(tasksdir,showWarnings=FALSE)

ufields <- c('HE_Uacid','HE_Uacid_etc','uric_acid_mg_dl','uric_acid_mg_dl_alt_lod','uric_acid_below_lod','uric_acid_above_uloq','uric_acid_qualifier_boundary')
raw <- dplyr::bind_rows(lapply(2019:2024,function(yr) {
 z <- as.data.frame(haven::read_sas(file.path(Sys.getenv('KNHANES_EXTRACT_ROOT'),sprintf('hn%d_all.sas7bdat',yr%%100)),col_select=c('ID','HEI','HE_Uacid')))
 z$year <- yr;z
}))
idx <- match(paste(frame$year,frame$ID),paste(raw$year,raw$ID))
stopifnot(!anyNA(idx),!anyDuplicated(paste(raw$year,raw$ID)),isTRUE(all.equal(as.numeric(frame$HE_Uacid),as.numeric(raw$HE_Uacid[idx]))))
frame$khei <- as.numeric(raw$HEI[idx])
stopifnot(all(is.na(frame$khei) | (frame$khei>=0 & frame$khei<=100)))
frame$uric_acid_observed <- frame$base_domain %in% TRUE & is.finite(frame$uric_acid_mg_dl)
frame$hyperuricemia <- ifelse(frame$uric_acid_observed,as.integer(ifelse(frame$sex=='Male',frame$uric_acid_mg_dl>7,frame$uric_acid_mg_dl>5.7)),NA_integer_)
OUTCOMES$serum_uric_acid <- list(variable='uric_acid_mg_dl',family='gaussian',domain='uric_acid_observed')
OUTCOMES$hyperuricemia <- list(variable='hyperuricemia',family='binomial',domain='uric_acid_observed')
primary_outcomes <- setdiff(names(OUTCOMES),'hyperuricemia')
for(j in 1:4) frame[[paste0('nova',j,'_share_per10')]] <- 10*frame[[paste0('nova',j,'_g')]]/frame$classified_nonwater_g
frame$total_classified_kg <- frame$classified_nonwater_g/1000
keep <- frame$base_domain %in% TRUE
stopifnot(max(abs(rowSums(frame[keep,paste0('nova',1:4,'_share_per10')])-10))<1e-8)
frame <- derive_fixed_exposure_categories(frame)
frame <- derive_ga_fields(frame)
frame$ga_domain <- frame$ga_domain & frame$base_domain %in% TRUE
saveRDS(frame,file.path(output,'analysis_frame_with_auxiliary_fields.rds'))
augmented_sha <- sha256_file(file.path(output,'analysis_frame_with_auxiliary_fields.rds'))
codepaths <- c(file.path(root,'01_code/04_run_models.R'),list.files(file.path(root,'source_snapshot'),pattern='\\.R$',full.names=TRUE),list.files(file.path(root,'reuse_sources'),pattern='\\.R$',full.names=TRUE))
write.csv(data.frame(path=codepaths,sha256=vapply(codepaths,sha256_file,character(1))),file.path(output,'model_code_manifest.csv'),row.names=FALSE)

find_assignment <- function(expr,name) {
 if(is.call(expr) && identical(expr[[1]],as.name('<-')) && identical(expr[[2]],as.name(name))) return(expr)
 if(is.call(expr) && !as.character(expr[[1]])[1] %in% c('if','{')) return(NULL)
 if(is.call(expr) || is.expression(expr)) for(child in as.list(expr)) {
   ans <- find_assignment(child,name); if(!is.null(ans)) return(ans)
 }
 NULL
}
expr <- find_assignment(parse(file.path(root,'reuse_sources/final_rcs_source.R')),'fit_kn_v3')
stopifnot(!is.null(expr));eval(expr)
knots_pp <- weighted_quantile(frame$upf_grams_pct[keep],frame$analysis_weight[keep],c(.1,.5,.9))
display_pp <- weighted_quantile(frame$upf_grams_pct[keep],frame$analysis_weight[keep],c(.01,.99))
reference_pp <- 10
grid_pp_v3 <- sort(unique(c(seq(display_pp[1],display_pp[2],length.out=121),reference_pp,knots_pp)))
sg_body <- deparse(body(fit_subgroup_outcome_stream),width.cutoff=500)
hit <- grepl('subgroup == "age_group"',sg_body,fixed=TRUE)
stopifnot(sum(hit)==1)
line <- which(hit)
if(!grepl('setdiff',sg_body[line],fixed=TRUE)) {
 stopifnot(grepl('setdiff(covars, "age")',sg_body[line+1],fixed=TRUE));sg_body <- sg_body[-c(line,line+1)]
} else sg_body <- sg_body[-line]
body(fit_subgroup_outcome_stream) <- parse(text=paste(sg_body,collapse='\n'))[[1]]
writeLines(deparse(fit_subgroup_outcome_stream),file.path(output,'age_adjusted_subgroup_function.R'))

custom_fit <- function(outcome,terms='upf_per10',covars=NULL,restriction=NULL) {
 spec <- OUTCOMES[[outcome]]
 if(is.null(covars)) covars <- if(outcome=='DKD') MODEL_COVARIATES_DKD$M3 else MODEL_COVARIATES$M3
 form <- build_model_formula(spec$variable,terms,covars)
 fits <- lapply(seq_len(imp$m),function(i) {
   d <- completed_analysis_dataset(frame,imp,i)
   if(!is.null(restriction)) d[[spec$domain]] <- d[[spec$domain]] %in% TRUE & restriction(d)
   fit_svy_formula(d,form,spec)
 })
 pooled <- pool_svy_fits(fits,outcome,'M3','custom',spec$family)
 wanted <- pooled$term %in% terms
 list(result=pooled[wanted,],n=vapply(fits,`[[`,numeric(1),'n'),events=vapply(fits,`[[`,numeric(1),'events'),formula=paste(deparse(form),collapse=' '))
}
run_one <- function(job) {
 kind <- job$kind; oc <- job$outcome
 if(kind=='main') return(fit_mi_model_stream(frame,imp,oc,job$model,job$form))
 if(kind=='diet') {
   cov <- if(oc=='DKD') MODEL_COVARIATES_DKD$M3 else MODEL_COVARIATES$M3
   restrict <- function(d) is.finite(d$khei)
   a <- custom_fit(oc,covars=cov,restriction=restrict)
   b <- custom_fit(oc,covars=c(cov,'khei'),restriction=restrict)
   stopifnot(identical(a$n,b$n),identical(a$events,b$events));return(list(before=a,after=b))
 }
 if(kind=='nova') return(custom_fit(oc,terms=paste0('nova',2:4,'_share_per10'),covars=c(if(oc=='DKD') MODEL_COVARIATES_DKD$M3 else MODEL_COVARIATES$M3,'total_classified_kg')))
 if(kind=='rcs') return(fit_kn_v3(frame,imp,oc))
 if(kind=='subgroup') return(fit_subgroup_outcome_stream(frame,imp,oc,job$subgroup))
 if(kind=='sensitivity') {
   if(job$scenario=='complete_case') return(fit_mi_model(list(frame),oc,'M3','continuous',complete_case=TRUE))
   restrict <- if(job$scenario=='energy') function(d) is.finite(d$energy_kcal)&d$energy_kcal>=500&d$energy_kcal<=5000 else function(d) d$upf_grams_pct>=display_pp[1]&d$upf_grams_pct<=display_pp[2]
   return(custom_fit(oc,restriction=restrict))
 }
 if(kind=='ga') {
   fits <- lapply(seq_len(imp$m),function(i) fit_parallel_ordinal(completed_analysis_dataset(frame,imp,i),MODEL_COVARIATES[[job$model]],make_survey_design))
   return(pool_scalar(fits,'upf_per10',sign_multiplier=-1,label='log_OR_worse_GA_risk',exponentiate=TRUE))
 }
 if(kind=='ga_assumption') {
   fits <- lapply(seq_len(imp$m),function(i) fit_threshold_stack(completed_analysis_dataset(frame,imp,i),MODEL_COVARIATES$M3,make_survey_design))
   terms <- unname(fits[[1]]$exposure_terms[c('T1','T2','T3')])
   test <- d1_joint_from_contrasts(fits,terms,rbind(c(-1,1,0),c(-1,0,1)))
   return(list(test=test,thresholds=lapply(terms,function(t) pool_scalar(fits,t,exponentiate=TRUE))))
 }
 if(kind=='age_diagnostic') {
   spec<-OUTCOMES[[oc]];domain<-frame[[spec$domain]] %in% TRUE
   k<-weighted_quantile(frame$age[domain],frame$analysis_weight[domain],c(.1,.5,.9))
   cov<-if(oc=='DKD') MODEL_COVARIATES_DKD$M3 else MODEL_COVARIATES$M3
   form<-build_model_formula(spec$variable,'upf_per10',c(cov,'.age_nonlinear'))
   fits<-lapply(seq_len(imp$m),function(i) {
     d<-completed_analysis_dataset(frame,imp,i);d$.age_nonlinear<-rcs_nonlinear_basis(d$age,k)/10
     fit_svy_formula(d,form,spec)
   })
   z<-pool_svy_fits(fits,oc,'M3_age_rcs3','continuous',spec$family)
   return(list(result=z[z$term=='upf_per10',],knots=k,formula=paste(deparse(form),collapse=' ')))
 }
 stop('Unknown task')
}
jobs <- list()
add <- function(...) {jobs[[length(jobs)+1L]] <<- list(...)}
for(oc in names(OUTCOMES)) for(model in names(MODEL_COVARIATES)) for(form in c('continuous','quartile','trend')) add(kind='main',outcome=oc,model=model,form=form)
for(oc in primary_outcomes) {
 for(k in c('diet','nova','rcs')) add(kind=k,outcome=oc)
 for(s in c('complete_case','energy','upf_p01_p99')) add(kind='sensitivity',outcome=oc,scenario=s)
 for(s in names(SUBGROUPS)) if(!(oc=='DKD' && s=='diabetes')) add(kind='subgroup',outcome=oc,subgroup=s)
}
for(model in names(MODEL_COVARIATES)) add(kind='ga',outcome='G_A',model=model)
add(kind='ga_assumption',outcome='G_A')
for(oc in primary_outcomes) add(kind='age_diagnostic',outcome=oc)
for(i in seq_along(jobs)) jobs[[i]]$id <- sprintf('%03d_%s_%s',i,jobs[[i]]$kind,jobs[[i]]$outcome)
jsonlite::write_json(jobs,file.path(output,'job_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
workers <- 14L
cat('MODEL_TASKS',length(jobs),'WORKERS',workers,'START',as.character(Sys.time()),'\n')
results <- parallel::mclapply(jobs,function(job) {
 configure_worker_runtime()
 p <- file.path(tasksdir,paste0(job$id,'.rds'))
 if(file.exists(p)) return(list(id=job$id,status='CACHED'))
 tryCatch({
   ans <- run_one(job)
   saveRDS(list(job=job,result=ans,frame_sha256=sha256_file(framepath),augmented_frame_sha256=augmented_sha,mids_sha256=sha256_file(mipath)),p)
   cat('DONE',job$id,as.character(Sys.time()),'\n')
   list(id=job$id,status='OK')
 },error=function(e) {cat('FAILED',job$id,conditionMessage(e),'\n');list(id=job$id,status='ERROR',message=conditionMessage(e))})
},mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=TRUE)
jsonlite::write_json(results,file.path(output,'task_status.json'),auto_unbox=TRUE,pretty=TRUE)
stopifnot(all(vapply(results,function(x) is.list(x)&&x$status %in% c('OK','CACHED'),logical(1))))
writeLines('All requested candidate model tasks completed; requires result and provenance review.',file.path(output,'MODELS_COMPLETED.txt'))
cat('FINISHED',as.character(Sys.time()),'\n')
)--------------------"

upf_source[["KN/01_code/08_summarize_models.R"]] <- r"--------------------(Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
suppressPackageStartupMessages(library(dplyr))
base <- '/storage/home/tmu2301/nhanes分析';root <- file.path(Sys.getenv("UPF_RUN_ROOT"),"KN")
indir <- file.path(root,'04_results/candidate');out <- file.path(indir,'summary');dir.create(out,showWarnings=FALSE)
stopifnot(file.exists(file.path(indir,'MODELS_COMPLETED.txt')))
jobs <- lapply(list.files(file.path(indir,'tasks'),pattern='\\.rds$',full.names=TRUE),readRDS)
actual_augmented_sha <- digest::digest(file=file.path(indir,'analysis_frame_with_auxiliary_fields.rds'),algo='sha256')
actual_mids_sha <- digest::digest(file=file.path(root,'02_derived/candidate_mids.rds'),algo='sha256')
stopifnot(all(vapply(jobs,function(x)identical(x$augmented_frame_sha256,actual_augmented_sha)&&identical(x$mids_sha256,actual_mids_sha),logical(1))))
getkind <- function(k) Filter(function(x) x$job$kind==k,jobs)
write <- function(d,name) {
 d$result_version <- 'KN_NOVA_IDENTITY_CORRECTION_20260913'
 write.csv(d,file.path(out,name),row.names=FALSE,na='');invisible(d)
}
scale_display <- function(d) {
 d$display_estimate <- d$estimate;d$display_low <- d$conf_low;d$display_high <- d$conf_high
 hit <- d$outcome=='log_uacr'
 for(nm in c('display_estimate','display_low','display_high')) d[[nm]][hit] <- 100*expm1(d[[nm]][hit])
 d$display_scale <- ifelse(hit,'percent_difference',d$effect_scale)
 d
}
main <- bind_rows(lapply(getkind('main'),function(x) {
 a <- x$result$result
 a <- a[grepl('^upf_per10$|^upf_quartileQ[2-4]$|^upf_q_median_per10$',a$term),]
 stopifnot(all(x$result$metadata$success),all(a$status=='OK'))
 a$n <- min(x$result$metadata$n);a$events <- if(all(is.na(x$result$metadata$events))) NA_real_ else min(x$result$metadata$events)
 a
}))
stopifnot(nrow(main)==120L,all(is.finite(as.matrix(main[,c('estimate','std_error_link','conf_low','conf_high','p_value')]))))
main <- scale_display(main);write(main,'KNHANES_main_all_exposure_terms.csv')
meta <- bind_rows(lapply(getkind('main'),function(x) x$result$metadata));stopifnot(nrow(meta)==72L*50L)
write(meta,'main_model_fit_metadata.csv')
m3 <- subset(main,model=='M3'&exposure_form=='continuous'&outcome!='hyperuricemia')
stopifnot(nrow(m3)==5);m3$q_BH5 <- p.adjust(m3$p_value,'BH');write(m3,'KNHANES_primary_M3_BH5.csv')
write(subset(main,outcome=='hyperuricemia'&exposure_form=='continuous'),'KNHANES_hyperuricemia_M1_M4.csv')
diet <- bind_rows(lapply(getkind('diet'),function(x) {
 bind_rows(lapply(c('before','after'),function(m) {
   z <- x$result[[m]];a <- z$result;a$comparison <- m;a$n <- min(z$n);a$events <- if(all(is.na(z$events))) NA_real_ else min(z$events);a
 }))
}))
diet$q_BH5 <- NA_real_;i <- diet$comparison=='after';stopifnot(sum(i)==5);diet$q_BH5[i] <- p.adjust(diet$p_value[i],'BH')
write(scale_display(diet),'KNHANES_KHEI_same_sample.csv')
nova <- bind_rows(lapply(getkind('nova'),function(x) {a<-x$result$result;a$n<-min(x$result$n);a$events<-if(all(is.na(x$result$events))) NA_real_ else min(x$result$events);a}))
nova$q_BH5 <- NA_real_;i <- nova$term=='nova4_share_per10';stopifnot(sum(i)==5);nova$q_BH5[i] <- p.adjust(nova$p_value[i],'BH')
write(scale_display(nova),'KNHANES_NOVA_substitution.csv')
sens <- bind_rows(lapply(getkind('sensitivity'),function(x) {
 a<-x$result$result;a<-a[a$term=='upf_per10',];a$scenario<-x$job$scenario
 a$n<-if(!is.null(x$result$metadata)) min(x$result$metadata$n) else min(x$result$n);a
}))
stopifnot(nrow(sens)==15);write(scale_display(sens),'KNHANES_sensitivity.csv')
age <- bind_rows(lapply(getkind('age_diagnostic'),function(x) {z<-x$result$result;z$knots<-paste(x$result$knots,collapse=';');z$formula<-x$result$formula;z}))
stopifnot(nrow(age)==5);write(scale_display(age),'KNHANES_age_function_diagnostic.csv')
curves <- bind_rows(lapply(getkind('rcs'),function(x) x$result$curve))
rcstests <- bind_rows(lapply(getkind('rcs'),function(x) x$result$test))
stopifnot(nrow(rcstests)==10,all(is.finite(rcstests$p_value)),all(curves$reference_pp==10))
write(scale_display(curves),'KNHANES_RCS_curves.csv');write(rcstests,'KNHANES_RCS_tests.csv')
sg <- bind_rows(lapply(getkind('subgroup'),function(x) x$result$effects))
ints <- bind_rows(lapply(getkind('subgroup'),function(x) x$result$interaction))
stopifnot(nrow(ints)==24);ints$q_BH24<-p.adjust(ints$p_interaction,'BH')
write(scale_display(sg),'KNHANES_subgroup_effects.csv');write(ints,'KNHANES_interactions_BH24.csv')
ga <- bind_rows(lapply(getkind('ga'),function(x) {a<-x$result;a$model<-x$job$model;a}))
write(ga,'KNHANES_GA_M1_M4.csv')
gat <- getkind('ga_assumption')[[1]]$result$test
write(gat,'KNHANES_GA_assumption.csv')
)--------------------"

upf_source[["helpers/table1.R"]] <- r"--------------------(Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
suppressPackageStartupMessages({library(survey);library(haven)})
options(survey.lonely.psu='adjust',survey.adjust.domain.lonely=TRUE)
base <- '/storage/home/tmu2301/nhanes分析';root <- file.path(Sys.getenv("UPF_RUN_ROOT"),"KN")
out <- file.path(root,'04_results/Table1_observed_candidate');dir.create(out,recursive=TRUE,showWarnings=FALSE)
template <- read.csv(file.path(root,'source_snapshot/Table1_baseline_six_outcomes.csv'),check.names=FALSE)
nhp <- file.path(Sys.getenv("UPF_RUN_ROOT"),'NHANES/03_派生数据/正式五结局_v1/三模块探索/analysis_frame_three_modules.rds')
nh <- readRDS(nhp)
kn <- readRDS(file.path(root,'04_results/candidate/analysis_frame_with_auxiliary_fields.rds'))
corrected <- readRDS(file.path(root,'02_derived/analysis_frame_fasting_corrected.rds'))
stopifnot(identical(kn$person_uid,corrected$person_uid))
for(v in c('diabetes','diabetes_original','dkd','dkd_domain','dkd_alt_lod')) kn[[v]] <- corrected[[v]]
canonical <- readRDS(file.path(Sys.getenv("UPF_RUN_ROOT"),'NHANES/03_派生数据/正式五结局_v1/analysis_frame_pre_mice.rds'))
stopifnot(identical(nh$SEQN,canonical$SEQN))
frozen <- c('base_domain','upf_gram_ratio_nonwater','analysis_weight','age','sex','bmi','education','diabetes','ckd','dkd')
stopifnot(all(vapply(frozen,function(v) identical(nh[[v]],canonical[[v]]),logical(1))))
for(db in c('nh','kn')) {d<-get(db);v<-if(db=='nh')'mean_energy_kcal' else 'energy_kcal';d[[v]]<-d[[v]]*4.184;v<-if(db=='nh')'uric_acid_dxc' else 'uric_acid_mg_dl';d[[v]]<-d[[v]]*59.48;d$acr_mg_g<-d$acr_mg_g/8.84;assign(db,d)}
normal <- function(x) tolower(gsub('[^a-z0-9]','',tolower(as.character(x))))
fmtp <- function(p) if(!is.finite(p)) '' else if(p<.001) '<0.001' else sprintf('%.3f',p)
details <- list();missing <- list();cutsall <- list();res <- template
for(panel in 0:1) {
 d <- if(panel==0) nh else kn
 exposure <- if(panel==0) 'upf_gram_ratio_nonwater' else 'upf_grams_pct'
 d$.base <- d$base_domain %in% TRUE
 d$.agegroup <- cut(d$age,c(-Inf,45,65,Inf),right=FALSE,labels=c('20–44 years','45–64 years','≥65 years'))
 d$.bmigroup <- cut(d$bmi,c(-Inf,25,30,Inf),right=FALSE,labels=c('<25 kg/m²','25–<30 kg/m²','≥30 kg/m²'))
 d$.lowegfr <- ifelse(is.finite(d$egfr_2021),as.integer(d$egfr_2021<60),NA)
 d$.highuacr <- ifelse(is.finite(d$acr_mg_g),as.integer(d$acr_mg_g>=30/8.84),NA)
 if(panel==1) for(j in 1:3) d[[paste0('nova',j,'_share')]] <- d[[paste0('nova',j,'_share_per10')]]*10
 ok <- is.finite(d$analysis_weight)&d$analysis_weight>0
 if(panel==0) {
   ok <- ok&!is.na(d$SDMVSTRA)&!is.na(d$SDMVPSU)
   d <- d[ok,];d$.strata <- interaction(d$cycle,d$SDMVSTRA,drop=TRUE);d$.psu <- interaction(d$cycle,d$SDMVPSU,drop=TRUE)
 } else {d <- d[ok,];d$.strata <- d$strata_uid;d$.psu <- d$psu_uid}
 design <- svydesign(ids=~.psu,strata=~.strata,weights=~analysis_weight,nest=TRUE,data=d)
 bdes <- subset(design,.base)
 cuts <- as.numeric(coef(svyquantile(reformulate(exposure),bdes,quantiles=c(.25,.5,.75),ci=FALSE,na.rm=TRUE)))
 design$variables$.q <- cut(design$variables[[exposure]],c(-Inf,cuts,Inf),labels=paste0('Q',1:4),include.lowest=TRUE)
 cutsall[[panel+1]] <- data.frame(database=if(panel==0)'NHANES' else 'KNHANES',cut=c('Q1_Q2','Q2_Q3','Q3_Q4'),value=cuts)
 varmap <- c('2'=exposure,'3'='age','4'='.agegroup','5'='sex','6'='race_ethnicity','7'='education','8'='smoking','9'='alcohol','10'='physical_activity','11'='bmi','12'='.bmigroup','13'=if(panel==0)'mean_energy_kcal' else 'energy_kcal','14'='nova1_share','15'='nova2_share','16'='nova3_share','17'='hei2015','18'='khei','19'='hypertension','20'='diabetes','21'='high_cholesterol','22'='ckd','23'='.lowegfr','24'='.highuacr','25'='dkd','26'='kidney_stones','27'='egfr_2021','28'='acr_mg_g','29'=if(panel==0)'uric_acid_dxc' else 'uric_acid_mg_dl')
 for(r in which(template$panel_rank==panel)) {
   rank <- as.character(template$group_rank[r]);level <- template$Level[r]
   if(rank=='1') {res[r,c('Overall','Q1','Q2','Q3','Q4')] <- c('—',sprintf('≤%.2f',cuts[1]),sprintf('%.2f–%.2f',cuts[1],cuts[2]),sprintf('%.2f–%.2f',cuts[2],cuts[3]),sprintf('>%.2f',cuts[3]));next}
   domain <- design$variables$.base
   if(rank=='25') domain <- domain&design$variables$dkd_domain %in% TRUE
   if(rank!='0') {
     var <- unname(varmap[rank]);stopifnot(length(var)==1,!is.na(var),var %in% names(design$variables))
     val <- design$variables[[var]];available <- !is.na(val)
     categorical <- !rank %in% c('2','3','11','13','14','15','16','17','18','27','28','29')
     if(categorical) {
       if(level=='Yes') {
         indicator <- as.numeric(as.character(val) %in% c('Yes','1'))
       } else {
         found <- unique(as.character(val))[normal(unique(as.character(val)))==normal(level)&!is.na(unique(as.character(val)))]
         stopifnot(length(found)==1);indicator <- as.numeric(as.character(val)==found)
       }
       indicator[!available] <- NA;design$variables$.value <- indicator
     } else design$variables$.value <- as.numeric(val)
     design$variables$.original <- val
     md <- data.frame(panel=panel,variable=var,domain_n=sum(domain),observed_n=sum(domain&available),missing_n=sum(domain&!available))
     missing[[length(missing)+1]] <- md
     dd <- design[domain&available,]
     p <- if(categorical) as.numeric(svychisq(~.original+.q,dd,statistic='F')$p.value) else as.numeric(regTermTest(svyglm(.value~.q,dd),~.q)$p)
     res$P_value[r] <- fmtp(p)
   }
   for(group in c('Overall',paste0('Q',1:4))) {
     hit <- domain & if(group=='Overall') TRUE else design$variables$.q %in% group
     if(rank=='0') {res[r,group] <- format(sum(hit),big.mark=',',trim=TRUE);next}
     hit <- hit&available;dd <- design[hit,]
     if(rank=='28') {
       qq <- as.numeric(coef(svyquantile(~.value,dd,quantiles=c(.25,.5,.75),ci=FALSE,na.rm=TRUE)))
       value <- sprintf('%.2f (%.2f–%.2f)',qq[2],qq[1],qq[3]);estimate <- qq[2];se <- NA_real_
     } else {
       z <- svymean(~.value,dd);estimate <- as.numeric(coef(z));se <- as.numeric(SE(z))
       value <- if(categorical) sprintf('%.1f%% (%.1f)',100*estimate,100*se) else sprintf('%.2f (%.2f)',estimate,se)
     }
     res[r,group] <- value
     details[[length(details)+1]] <- data.frame(panel=panel,characteristic=template$Characteristic[r],level=level,variable=var,group=group,n=sum(hit),estimate=estimate,se=se,p=p,source='observed pre-imputation, survey domain analysis')
   }
 }
}
write.csv(res,file.path(out,'Table1_observed_candidate.csv'),row.names=FALSE,fileEncoding='UTF-8')
write.csv(do.call(rbind,details),file.path(out,'Table1_numeric_lineage.csv'),row.names=FALSE)
write.csv(unique(do.call(rbind,missing)),file.path(out,'Table1_missingness_and_denominators.csv'),row.names=FALSE)
write.csv(do.call(rbind,cutsall),file.path(out,'Table1_weighted_quartiles.csv'),row.names=FALSE)
writeLines(c('Candidate Table 1 uses observed pre-imputation data for both descriptive summaries and P tests.','Denominators vary by variable; missingness is documented separately. DKD summaries use the observed diabetes-positive kidney-status domain.','This changes descriptive presentation only. Main regression analyses continue to use 50 multiply imputed datasets.','Final manuscript/table notes must identify this observed-data convention explicitly; do not combine with the old blanket MI note.'),file.path(out,'TABLE1_METHOD_NOTE.txt'))
cat('TABLE1_OBSERVED_CANDIDATE_COMPLETED',nrow(res),'rows\n')
print(warnings())
)--------------------"

upf_source[["KN/reuse_sources/final_rcs_source.R"]] <- r"--------------------(
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !toupper(args[[1]]) %in% c("NHANES", "KNHANES")) {
  stop("Usage: Rscript 12_recalculate_rcs_v3.R NHANES|KNHANES", call. = FALSE)
}
cohort <- toupper(args[[1]])
figure_root <- Sys.getenv(
  "FIGURE_PROJECT_ROOT",
  "/project/08_NHANES_KNHANES_联合论文图表_20260828"
)
source_dir <- file.path(figure_root, "02_绘图数据_v3重算", "统计重算源")
log_dir <- file.path(figure_root, "08_运行日志")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(log_dir, paste0("12_recalculate_rcs_v3_", tolower(cohort), ".log"))
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)
cat("cohort=", cohort, "\nstarted_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")

write_csv <- function(x, filename) {
  path <- file.path(source_dir, filename)
  readr::write_excel_csv(x, path, na = "")
  path
}

write_manifest <- function(paths, prefix) {
  paths <- normalizePath(paths, mustWork = TRUE)
  info <- file.info(paths)
  manifest <- data.frame(
    cohort = cohort,
    path = paths,
    bytes = info$size,
    modified = format(info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
    sha256 = vapply(paths, digest::digest, character(1), file = TRUE, algo = "sha256", serialize = FALSE),
    stringsAsFactors = FALSE
  )
  write_csv(manifest, paste0(prefix, "_manifest.csv"))
}

rcs3_basis_scaled <- function(x_pp, knots_pp, unit_pp = 10) {
  x <- as.numeric(x_pp) / unit_pp
  knots <- as.numeric(knots_pp) / unit_pp
  if (length(knots) != 3L || any(!is.finite(knots)) || any(diff(knots) <= 0)) {
    stop("RCS requires three finite ordered knots", call. = FALSE)
  }
  k1 <- knots[[1]]; k2 <- knots[[2]]; k3 <- knots[[3]]
  nonlinear <- (
    pmax(x - k1, 0)^3 -
      pmax(x - k2, 0)^3 * (k3 - k1) / (k3 - k2) +
      pmax(x - k3, 0)^3 * (k2 - k1) / (k3 - k2)
  ) / (k3 - k1)^2
  cbind(rcs_linear = x, rcs_nonlin1 = nonlinear)
}

if (cohort == "NHANES") {
  nh_root <- Sys.getenv(
    "NHANES_ROOT",
    "/project/06_NHANES_UPF七结局_第一阶段探索_20260824"
  )
  nh_code <- file.path(nh_root, "02_代码", "正式五结局")
  source(file.path(nh_code, "23_extension_common.R"), local = FALSE)
  require_formal_authorization("v3_rcs_recalculation")

  frame_path <- file.path(DIRS$derived, "analysis_frame_pre_mice.rds")
  d <- readRDS(frame_path)
  full0 <- make_full_design(d)
  base <- full0$variables$base_domain %in% TRUE &
    is.finite(full0$variables[[EXPOSURE_VAR]])
  design0 <- full0[base, ]
  probs <- c(0.01, 0.10, 0.50, 0.90, 0.99)
  quantiles <- weighted_quantiles(design0, EXPOSURE_VAR, probs)
  names(quantiles) <- c("p01", "p10", "p50", "p90", "p99")
  knots_pp <- unname(quantiles[c("p10", "p50", "p90")])
  reference_pp <- 10
  display_pp <- unname(quantiles[c("p01", "p99")])
  grid_pp_v3 <- sort(unique(c(
    seq(display_pp[[1]], display_pp[[2]], length.out = 121L),
    reference_pp, knots_pp
  )))
  basis <- rcs3_basis_scaled(d[[EXPOSURE_VAR]], knots_pp, EXPOSURE_UNIT_PP)
  d$rcs_linear_v3 <- basis[, "rcs_linear"]
  d$rcs_nonlin1_v3 <- basis[, "rcs_nonlin1"]

  rcs_formula_v3 <- function(outcome) {
    stats::reformulate(
      c(
        "rcs_linear_v3", "rcs_nonlin1_v3",
        phase1_model_covariates("M3", "linear", outcome, FALSE)
      ),
      response = OUTCOMES[[outcome]]$variable
    )
  }

  fit_one <- function(imputation_number, outcome, mi_object) {
    completed <- restore_completed_subframe(d, mi_object, imputation_number)
    full <- make_full_design(completed)
    domain <- full$variables[[OUTCOMES[[outcome]]$domain]] %in% TRUE &
      is.finite(full$variables[[EXPOSURE_VAR]])
    form <- rcs_formula_v3(outcome)
    model_vars <- all.vars(form)
    incomplete <- model_vars[vapply(full$variables[model_vars], function(x) {
      any(is.na(x[domain]))
    }, logical(1))]
    if (length(incomplete)) {
      stop("Post-MICE missing V3 RCS fields: ", paste(incomplete, collapse = ", "), call. = FALSE)
    }
    design <- full[domain, ]
    fit <- survey::svyglm(
      form, design = design, family = outcome_family(outcome), na.action = na.fail
    )
    w <- as.numeric(stats::weights(design, "sampling"))
    y <- design$variables[[OUTCOMES[[outcome]]$variable]]
    list(
      coefficients = stats::coef(fit),
      variance = stats::vcov(fit),
      df_complete = stats::df.residual(fit),
      audit = c(
        n_unweighted = nrow(design$variables),
        events = if (OUTCOMES[[outcome]]$family == "quasibinomial") sum(y == 1L) else NA_real_,
        weighted_population = sum(w),
        effective_sample_size = sum(w)^2 / sum(w^2),
        survey_design_df = survey::degf(design),
        complete_data_residual_df = stats::df.residual(fit)
      )
    )
  }

  curve_from_pool_v3 <- function(pooled, outcome) {
    grid <- grid_pp_v3
    basis_grid <- rcs3_basis_scaled(grid, knots_pp, EXPOSURE_UNIT_PP)
    basis_ref <- as.numeric(rcs3_basis_scaled(reference_pp, knots_pp, EXPOSURE_UNIT_PP)[1, ])
    terms <- c("rcs_linear_v3", "rcs_nonlin1_v3")
    beta <- pooled$coefficients[terms]
    variance <- pooled$variance[terms, terms, drop = FALSE]
    df <- min(pooled$df[match(terms, names(pooled$coefficients))], na.rm = TRUE)
    effect_scale <- effect_scale_for_outcome(outcome)
    do.call(rbind, lapply(seq_along(grid), function(i) {
      contrast <- as.numeric(basis_grid[i, ] - basis_ref)
      link <- sum(contrast * beta)
      se <- sqrt(as.numeric(t(contrast) %*% variance %*% contrast))
      transformed <- transform_link_contrast(link, se, df, effect_scale)
      data.frame(
        cohort = "NHANES", version = "v3_rcs_3knot_p10ref_20260828",
        outcome = outcome, model = "M3", exposure_value = grid[[i]],
        reference_value = reference_pp, knot_10 = knots_pp[[1]],
        knot_50 = knots_pp[[2]], knot_90 = knots_pp[[3]],
        display_p01 = display_pp[[1]], display_p99 = display_pp[[2]],
        effect_measure = unname(c(
          odds_ratio = "OR_vs_reference",
          identity = "beta_difference_vs_reference",
          log_response_percent = "percent_difference_vs_reference"
        )[effect_scale]),
        estimate = transformed[[1]], conf_low = transformed[[2]], conf_high = transformed[[3]],
        estimate_link = link, standard_error_link = se,
        stringsAsFactors = FALSE
      )
    }))
  }

  d1_joint_test_v3 <- function(summaries, terms) {
    aliases <- paste0("b", seq_along(terms))
    qhat <- lapply(summaries, function(z) {
      stats::setNames(as.numeric(z$coefficients[terms]), aliases)
    })
    uhat <- lapply(summaries, function(z) {
      v <- z$variance[terms, terms, drop = FALSE]
      dimnames(v) <- list(aliases, aliases)
      v
    })
    df_complete <- min(vapply(summaries, `[[`, numeric(1), "df_complete"))
    test <- mitml::testConstraints(
      qhat = qhat, uhat = uhat, constraints = aliases,
      method = "D1", df.com = df_complete
    )$test[1, ]
    values <- unname(test[c("F.value", "df1", "df2", "P(>F)")])
    if (any(!is.finite(values)) || test[["df1"]] <= 0 || test[["df2"]] <= 0) {
      stop("NHANES V3 D1 test returned invalid results", call. = FALSE)
    }
    list(
      f_value = unname(test[["F.value"]]), df1 = unname(test[["df1"]]),
      df2 = unname(test[["df2"]]), p_value = unname(test[["P(>F)"]])
    )
  }

  test_rows <- list(); curve_rows <- list(); pooled_models <- list()
  workers <- min(4L, n_cores_resolve())
  for (outcome in names(OUTCOMES)) {
    cat("outcome_start=", outcome, "\n", sep = "")
    mi_object <- readRDS(file.path(DIRS$derived, paste0("mi_", outcome, ".rds")))
    stopifnot(mi_object$m == MI_M)
    summaries <- parallel::mclapply(
      seq_len(MI_M), fit_one, outcome = outcome, mi_object = mi_object,
      mc.cores = min(4L, workers), mc.preschedule = TRUE
    )
    if (any(vapply(summaries, inherits, logical(1), "try-error"))) {
      stop("NHANES V3 RCS failed for ", outcome, call. = FALSE)
    }
    pooled <- pool_svy_summaries(summaries)
    overall <- d1_joint_test_v3(summaries, c("rcs_linear_v3", "rcs_nonlin1_v3"))
    nonlinear <- d1_joint_test_v3(summaries, "rcs_nonlin1_v3")
    audit <- summaries[[1]]$audit
    test_rows[[outcome]] <- data.frame(
      cohort = "NHANES", version = "v3_rcs_3knot_p10ref_20260828",
      outcome = outcome, model = "M3", knot_probabilities = "0.10;0.50;0.90",
      knots = paste(format(knots_pp, digits = 12), collapse = ";"),
      reference_probability = NA_real_, reference = reference_pp,
      display_probabilities = "0.01;0.99", display_min = display_pp[[1]],
      display_max = display_pp[[2]], p_overall = overall$p_value,
      p_nonlinear = nonlinear$p_value, overall_f = overall$f_value,
      overall_df1 = overall$df1, overall_df2 = overall$df2,
      nonlinear_f = nonlinear$f_value, nonlinear_df1 = nonlinear$df1,
      nonlinear_df2 = nonlinear$df2, n_unweighted = audit[["n_unweighted"]],
      events = audit[["events"]], effective_sample_size = audit[["effective_sample_size"]],
      successful_imputations = length(summaries),
      formula = paste(deparse(rcs_formula_v3(outcome)), collapse = " "),
      covariates = paste(phase1_model_covariates("M3", "linear", outcome, FALSE), collapse = ";"),
      stringsAsFactors = FALSE
    )
    curve_rows[[outcome]] <- curve_from_pool_v3(pooled, outcome)
    pooled_models[[outcome]] <- list(
      coefficients = pooled$coefficients, variance = pooled$variance,
      df = pooled$df, missinfo = pooled$missinfo
    )
    rm(mi_object, summaries, pooled)
    invisible(gc(full = TRUE))
    cat("outcome_done=", outcome, "\n", sep = "")
  }
  tests <- do.call(rbind, test_rows); rownames(tests) <- NULL
  curves <- do.call(rbind, curve_rows); rownames(curves) <- NULL
  stopifnot(
    nrow(tests) == 5L, nrow(curves) == 5L * length(grid_pp_v3),
    all(tests$successful_imputations == 50L),
    all(is.finite(tests$p_overall)), all(is.finite(tests$p_nonlinear)),
    all(is.finite(curves$estimate)), all(curves$conf_low <= curves$estimate),
    all(curves$estimate <= curves$conf_high)
  )
  curve_path <- write_csv(curves, "NHANES_RCS_M3_v3_curves.csv")
  test_path <- write_csv(tests, "NHANES_RCS_M3_v3_tests.csv")
  pooled_path <- file.path(source_dir, "NHANES_RCS_M3_v3_pooled_models.rds")
  saveRDS(pooled_models, pooled_path, compress = "xz")
  lineage <- data.frame(
    cohort = "NHANES", role = c("analysis_frame", paste0("mi_", names(OUTCOMES))),
    path = c(frame_path, file.path(DIRS$derived, paste0("mi_", names(OUTCOMES), ".rds"))),
    stringsAsFactors = FALSE
  )
  lineage$sha256 <- vapply(lineage$path, digest::digest, character(1), file = TRUE, algo = "sha256", serialize = FALSE)
  lineage_path <- write_csv(lineage, "NHANES_RCS_M3_v3_input_lineage.csv")
  write_manifest(c(curve_path, test_path, pooled_path, lineage_path), "NHANES_RCS_M3_v3")
}

if (cohort == "KNHANES") {
  kn_root <- Sys.getenv(
    "KNHANES_ROOT",
    "/project/05_KNHANES_外部验证_2019_2024_20260820"
  )
  kn_code <- file.path(kn_root, "04_代码")
  source(file.path(kn_code, "07_extended_models.R"), local = FALSE)
  require_user_authorization("v3_rcs_recalculation")
  require_nova_gate("v3_rcs_recalculation")
  inputs <- load_model_inputs()
  frame_q <- derive_fixed_exposure_categories(inputs$frame)
  base <- frame_q$base_domain %in% TRUE
  knots_pp <- weighted_quantile(
    frame_q$upf_grams_pct[base], frame_q$analysis_weight[base], c(0.10, 0.50, 0.90)
  )
  display_pp <- weighted_quantile(
    frame_q$upf_grams_pct[base], frame_q$analysis_weight[base], c(0.01, 0.99)
  )
  reference_pp <- 10
  grid_pp_v3 <- sort(unique(c(
    seq(display_pp[[1]], display_pp[[2]], length.out = 121L),
    reference_pp, knots_pp
  )))

  fit_kn_v3 <- function(frame, imp, outcome_name) {
    spec <- OUTCOMES[[outcome_name]]
    covars <- if (outcome_name == "DKD") MODEL_COVARIATES_DKD$M3 else MODEL_COVARIATES$M3
    knots <- knots_pp / UPF_UNIT_PP
    ref_x <- reference_pp / UPF_UNIT_PP
    frame <- derive_fixed_exposure_categories(frame)
    form <- build_model_formula(spec$variable, c("upf_per10", "upf_rcs_nonlinear"), covars)
    fits <- lapply(seq_len(imp$m), function(i) {
      d <- completed_analysis_dataset(frame, imp, i)
      d$upf_rcs_nonlinear <- rcs_nonlinear_basis(d$upf_per10, knots)
      tryCatch(fit_svy_formula(d, form, spec), error = function(e) NULL)
    })
    if (any(vapply(fits, is.null, logical(1)))) {
      stop("KNHANES V3 RCS failed in at least one imputation for ", outcome_name)
    }
    invisible(lapply(fits, validate_compact_fit, family = spec$family,
                     label = paste0(outcome_name, " V3 RCS")))
    term_names <- names(fits[[1]]$coefficients)
    grid_pp <- grid_pp_v3
    curve <- dplyr::bind_rows(lapply(seq_along(grid_pp), function(j) {
      grid_x <- grid_pp[[j]] / UPF_UNIT_PP
      L <- stats::setNames(rep(0, length(term_names)), term_names)
      L[["upf_per10"]] <- grid_x - ref_x
      L[["upf_rcs_nonlinear"]] <-
        rcs_nonlinear_basis(grid_x, knots) - rcs_nonlinear_basis(ref_x, knots)
      z <- pool_scalar_contrast(fits, L)
      is_reference <- abs(grid_pp[[j]] - reference_pp) <=
        .Machine$double.eps * max(1, abs(reference_pp))
      if (is_reference) {
        z$estimate_link <- 0
        z$std_error_link <- 0
        z$conf_low_link <- 0
        z$conf_high_link <- 0
        z$p_value <- 1
      }
      if (spec$family == "binomial") {
        z$estimate <- exp(z$estimate_link)
        z$conf_low <- exp(z$conf_low_link)
        z$conf_high <- exp(z$conf_high_link)
        z$effect_scale <- "OR"
      } else {
        z$estimate <- z$estimate_link
        z$conf_low <- z$conf_low_link
        z$conf_high <- z$conf_high_link
        z$effect_scale <- "difference"
      }
      z$cohort <- "KNHANES"
      z$version <- "v3_rcs_3knot_p10ref_20260828"
      z$outcome <- outcome_name
      z$model <- "M3"
      z$upf_grams_pct <- grid_pp[[j]]
      z$reference_pp <- reference_pp
      z$knot_10 <- knots_pp[[1]]
      z$knot_50 <- knots_pp[[2]]
      z$knot_90 <- knots_pp[[3]]
      z$display_p01 <- display_pp[[1]]
      z$display_p99 <- display_pp[[2]]
      z$status <- "OK"
      z
    }))
    nonlinear <- mi_d1_test(fits, "upf_rcs_nonlinear")
    nonlinear$test <- "nonlinearity"
    overall <- mi_d1_test(fits, c("upf_per10", "upf_rcs_nonlinear"))
    overall$test <- "overall_spline"
    tests <- dplyr::bind_rows(nonlinear, overall)
    tests$cohort <- "KNHANES"
    tests$version <- "v3_rcs_3knot_p10ref_20260828"
    tests$outcome <- outcome_name
    tests$model <- "M3"
    tests$knot_probabilities <- "0.10;0.50;0.90"
    tests$knots <- paste(format(knots_pp, digits = 12), collapse = ";")
    tests$reference_probability <- NA_real_
    tests$reference <- reference_pp
    tests$display_probabilities <- "0.01;0.99"
    tests$display_min <- display_pp[[1]]
    tests$display_max <- display_pp[[2]]
    tests$n_unweighted <- fits[[1]]$n
    tests$events <- fits[[1]]$events
    tests$successful_imputations <- length(fits)
    tests$formula <- paste(deparse(form), collapse = " ")
    tests$covariates <- paste(covars, collapse = ";")
    tests$status <- "OK"
    list(curve = curve, tests = tests)
  }

  outcomes <- names(OUTCOMES)
  workers <- min(4L, N_WORKERS)
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = workers)
  by_outcome <- future.apply::future_lapply(
    outcomes,
    function(outcome_name) {
      configure_worker_runtime()
      local_inputs <- load_model_inputs()
      fit_kn_v3(local_inputs$frame, local_inputs$imp, outcome_name)
    },
    future.seed = KN_SEED,
    future.packages = c("survey", "mitools", "mitml", "mice", "dplyr")
  )
  curves <- dplyr::bind_rows(lapply(by_outcome, `[[`, "curve"))
  tests <- dplyr::bind_rows(lapply(by_outcome, `[[`, "tests"))
  stopifnot(
    nrow(curves) == length(outcomes) * length(grid_pp_v3), nrow(tests) == 8L,
    all(tests$successful_imputations == 50L), all(tests$status == "OK"),
    all(is.finite(tests$p_value)), all(is.finite(curves$estimate)),
    all(is.finite(curves$conf_low)), all(is.finite(curves$conf_high)),
    all(curves$conf_low <= curves$estimate), all(curves$estimate <= curves$conf_high)
  )
  curve_path <- write_csv(curves, "KNHANES_RCS_M3_v3_curves.csv")
  test_path <- write_csv(tests, "KNHANES_RCS_M3_v3_tests.csv")
  lineage <- data.frame(
    cohort = "KNHANES",
    role = c("analysis_frame", "mids", "code_review_marker", "nova_gate_marker", "active_nova_map"),
    path = c(inputs$frame_path, inputs$mice_path, CODE_REVIEW_MARKER, NOVA_GATE_MARKER, ACTIVE_NOVA_MAP),
    stringsAsFactors = FALSE
  )
  lineage$sha256 <- vapply(lineage$path, digest::digest, character(1), file = TRUE, algo = "sha256", serialize = FALSE)
  lineage_path <- write_csv(lineage, "KNHANES_RCS_M3_v3_input_lineage.csv")
  write_manifest(c(curve_path, test_path, lineage_path), "KNHANES_RCS_M3_v3")
}

cat("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")
cat("status=OK\n")
)--------------------"

upf_source[["KN/reuse_sources/ga_shared.R"]] <- r"--------------------(
suppressPackageStartupMessages({
  library(survey)
  library(mitools)
  library(VGAM)
  library(svyVGAM)
})
options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)

GA_LEVELS <- c("Low", "Moderate", "High", "Very high")
G_LEVELS <- c("G1", "G2", "G3a", "G3b", "G4", "G5")
A_LEVELS <- c("A1", "A2", "A3")
GA_MAP <- c(
  G1_A1 = "Low", G2_A1 = "Low",
  G1_A2 = "Moderate", G2_A2 = "Moderate", G3a_A1 = "Moderate",
  G1_A3 = "High", G2_A3 = "High", G3a_A2 = "High", G3b_A1 = "High",
  G3a_A3 = "Very high", G3b_A2 = "Very high", G3b_A3 = "Very high",
  G4_A1 = "Very high", G4_A2 = "Very high", G4_A3 = "Very high",
  G5_A1 = "Very high", G5_A2 = "Very high", G5_A3 = "Very high"
)

derive_ga_fields <- function(d, egfr = "egfr_2021", uacr = "acr_mg_g") {
  ok <- is.finite(d[[egfr]]) & is.finite(d[[uacr]])
  g <- cut(d[[egfr]], breaks = c(-Inf, 15, 30, 45, 60, 90, Inf), right = FALSE,
           labels = c("G5", "G4", "G3b", "G3a", "G2", "G1"))
  a <- cut(d[[uacr]], breaks = c(-Inf, 30, 300, Inf), right = FALSE,
           labels = c("A1", "A2", "A3"))
  risk <- rep(NA_character_, nrow(d))
  risk[ok] <- unname(GA_MAP[paste(g[ok], a[ok], sep = "_")])
  d$g_category <- factor(as.character(g), levels = G_LEVELS)
  d$a_category <- factor(as.character(a), levels = A_LEVELS)
  d$ga_risk <- ordered(risk, levels = GA_LEVELS)
  d$ga_domain <- ok
  d
}

add_fixed_quartiles <- function(d, base_domain, exposure, weight) {
  keep <- d[[base_domain]] %in% TRUE & is.finite(d[[exposure]]) & is.finite(d[[weight]]) & d[[weight]] > 0
  qfun <- function(p) {
    o <- order(d[[exposure]][keep]); x <- d[[exposure]][keep][o]; w <- d[[weight]][keep][o]
    x[which(cumsum(w) / sum(w) >= p)[1]]
  }
  cuts <- vapply(c(.25, .50, .75), qfun, numeric(1))
  if (any(diff(cuts) <= 0)) stop("Invalid weighted UPF quartile cuts")
  d$upf_quartile <- cut(d[[exposure]], c(-Inf, cuts, Inf), labels = paste0("Q", 1:4), include.lowest = TRUE)
  d$upf_quartile <- factor(d$upf_quartile, levels = paste0("Q", 1:4))
  attr(d, "upf_quartile_cuts") <- cuts
  d
}

weighted_indicator_rows <- function(design, database, group = "Overall") {
  rows <- lapply(GA_LEVELS, function(level) {
    design$variables$.indicator <- as.numeric(as.character(design$variables$ga_risk) == level)
    est <- svymean(~.indicator, design, na.rm = TRUE)
    ci <- confint(est)
    data.frame(database = database, group = group, risk = level,
               n_unweighted = sum(as.character(design$variables$ga_risk) == level, na.rm = TRUE),
               weighted_proportion = unname(coef(est)), se = unname(SE(est)),
               ci_low = ci[1], ci_high = ci[2])
  })
  do.call(rbind, rows)
}

weighted_ga_descriptives <- function(design, database) {
  overall <- weighted_indicator_rows(design, database, "Overall")
  by_q <- do.call(rbind, lapply(paste0("Q", 1:4), function(q) {
    sub <- design[design$variables$upf_quartile == q, ]
    weighted_indicator_rows(sub, database, q)
  }))
  list(risk = rbind(overall, by_q))
}

compact_fit <- function(fit, df_complete) {
  b <- coef(fit); v <- vcov(fit)
  if (!identical(names(b), rownames(v)) || any(!is.finite(b)) || any(!is.finite(v))) {
    stop("Non-finite or unnamed coefficient/covariance output")
  }
  list(coefficients = b, variance = v, df_complete = max(1, as.numeric(df_complete)))
}

pool_scalar <- function(summaries, term, sign_multiplier = 1, label = term, exponentiate = TRUE) {
  transformed <- lapply(summaries, function(s) {
    if (!term %in% names(s$coefficients)) stop("Missing term: ", term)
    b <- sign_multiplier * unname(s$coefficients[term])
    v <- unname(s$variance[term, term, drop = FALSE])
    list(coefficients = setNames(b, label),
         variance = matrix(v, 1, 1, dimnames = list(label, label)),
         df_complete = s$df_complete)
  })
  pooled <- mitools::MIcombine(
    lapply(transformed, `[[`, "coefficients"),
    lapply(transformed, `[[`, "variance"),
    df.complete = min(vapply(transformed, `[[`, numeric(1), "df_complete"))
  )
  b <- unname(pooled$coefficients[[1]]); se <- sqrt(pooled$variance[1, 1]); df <- unname(pooled$df[[1]])
  crit <- if (is.finite(df)) qt(.975, df) else qnorm(.975)
  lo <- b - crit * se; hi <- b + crit * se
  p <- if (is.finite(df)) 2 * pt(abs(b / se), df, lower.tail = FALSE) else 2 * pnorm(abs(b / se), lower.tail = FALSE)
  data.frame(estimate_link = b, standard_error_link = se, df = df, p_value = p,
             estimate = if (exponentiate) exp(b) else b,
             ci_low = if (exponentiate) exp(lo) else lo,
             ci_high = if (exponentiate) exp(hi) else hi,
             successful_imputations = length(summaries))
}

d1_joint_from_contrasts <- function(summaries, terms, contrast_matrix) {
  q <- do.call(cbind, lapply(summaries, function(s) as.numeric(contrast_matrix %*% s$coefficients[terms])))
  u <- simplify2array(lapply(summaries, function(s) contrast_matrix %*% s$variance[terms, terms, drop = FALSE] %*% t(contrast_matrix)))
  k <- nrow(q); m <- ncol(q); qbar <- rowMeans(q); ubar <- apply(u, c(1, 2), mean); b <- cov(t(q))
  r <- (1 + 1 / m) * sum(diag(b %*% solve(ubar))) / k
  f <- as.numeric(t(qbar) %*% solve((1 + r) * ubar) %*% qbar / k)
  tval <- k * (m - 1); dfcom <- min(vapply(summaries, `[[`, numeric(1), "df_complete"))
  a <- r * tval / (tval - 2); vstar <- ((dfcom + 1) / (dfcom + 3)) * dfcom
  c0 <- 1 / (tval - 4); c1 <- vstar - 2 * (1 + a); c2 <- vstar - 4 * (1 + a)
  z <- 1 / c2 + c0 * (a^2 * c1 / ((1 + a)^2 * c2)) +
    c0 * (8 * a^2 * c1 / ((1 + a) * c2^2) + 4 * a^2 / ((1 + a) * c2)) +
    c0 * (4 * a^2 / (c2 * c1) + 16 * a^2 * c1 / c2^3) + c0 * (8 * a^2 / c2^2)
  df2 <- 4 + 1 / z
  data.frame(f_statistic = f, df1 = k, df2 = df2, p_value = pf(f, k, df2, lower.tail = FALSE),
             relative_increase_variance = r, complete_data_df = dfcom)
}

fit_parallel_ordinal <- function(d, covariates, design_fun) {
  form <- reformulate(c("upf_per10", covariates), response = "ga_risk")
  des <- design_fun(d)
  des <- des[des$variables$ga_domain %in% TRUE & is.finite(des$variables$upf_per10), ]
  fit <- suppressWarnings(svyVGAM::svy_vglm(form, VGAM::cumulative(parallel = TRUE, reverse = FALSE), design = des, crit = "coef"))
  dfc <- survey::degf(des) - (length(coef(fit)) - 1)
  compact_fit(fit, dfc)
}

fit_threshold_stack <- function(d, covariates, design_fun) {
  d <- d[d$ga_domain %in% TRUE & is.finite(d$upf_per10), , drop = FALSE]
  nrisk <- as.integer(d$ga_risk)
  thresholds <- c(T1 = 2L, T2 = 3L, T3 = 4L)
  pieces <- lapply(names(thresholds), function(nm) {
    z <- d; z$threshold <- nm; z$.worse <- as.integer(nrisk >= thresholds[[nm]]); z
  })
  stacked <- do.call(rbind, pieces); rownames(stacked) <- NULL
  stacked$threshold <- factor(stacked$threshold, levels = names(thresholds))
  rhs <- c("0 + threshold", paste0("threshold:", c("upf_per10", covariates)))
  form <- as.formula(paste(".worse ~", paste(rhs, collapse = " + ")))
  des <- design_fun(stacked)
  fit <- svyglm(form, design = des, family = quasibinomial(), na.action = na.fail)
  out <- compact_fit(fit, df.residual(fit))
  candidates <- names(out$coefficients)[grepl("upf_per10", names(out$coefficients), fixed = TRUE)]
  term_map <- vapply(names(thresholds), function(nm) {
    hit <- candidates[grepl(nm, candidates, fixed = TRUE)]
    if (length(hit) != 1L) stop("Cannot map threshold exposure term for ", nm)
    hit
  }, character(1))
  out$exposure_terms <- term_map
  out
}

fit_multinomial <- function(d, covariates, design_fun) {
  d$ga_risk <- factor(as.character(d$ga_risk), levels = GA_LEVELS)
  form <- reformulate(c("upf_per10", covariates), response = "ga_risk")
  des <- design_fun(d)
  des <- des[des$variables$ga_domain %in% TRUE & is.finite(des$variables$upf_per10), ]
  fit <- suppressWarnings(svyVGAM::svy_vglm(form, VGAM::multinomial(refLevel = "Low"), design = des, crit = "coef"))
  dfc <- survey::degf(des) - (length(coef(fit)) - 1)
  out <- compact_fit(fit, dfc)
  terms <- names(out$coefficients)[grepl("upf_per10", names(out$coefficients), fixed = TRUE)]
  if (length(terms) != 3L) stop("Expected three multinomial UPF terms")
  out$exposure_terms <- setNames(terms, GA_LEVELS[-1])
  out
}

write_ga_cells <- function(design, database, path) {
  rows <- list(); k <- 0L
  for (q in c("Overall", paste0("Q", 1:4))) {
    desq <- if (q == "Overall") design else design[design$variables$upf_quartile == q, ]
    for (g in G_LEVELS) for (a in A_LEVELS) {
      k <- k + 1L
      desq$variables$.cell <- as.numeric(as.character(desq$variables$g_category) == g & as.character(desq$variables$a_category) == a)
      est <- svymean(~.cell, desq, na.rm = TRUE); ci <- confint(est)
      rows[[k]] <- data.frame(database = database, group = q, g_category = g, a_category = a,
                              n_unweighted = sum(desq$variables$.cell, na.rm = TRUE),
                              weighted_proportion = unname(coef(est)), se = unname(SE(est)),
                              ci_low = ci[1], ci_high = ci[2])
    }
  }
  write.csv(do.call(rbind, rows), path, row.names = FALSE)
}
)--------------------"

upf_source[["helpers/cholesterol.R"]] <- r"--------------------(# Candidate implementation for explicit questionnaire routing, not yet adopted.
derive_cholesterol_routed <- function(total_cholesterol, diagnosis,
                                      medication_recommended, medication_current,
                                      ever_checked, cycle) {
  sizes <- lengths(list(total_cholesterol,diagnosis,medication_recommended,
                       medication_current,ever_checked,cycle))
  stopifnot(length(unique(sizes))==1L, all(cycle %in% c('E','F','G','H','I','P')))
  binary <- function(x) ifelse(x %in% 1,1L,ifelse(x %in% 2,0L,NA_integer_))
  measured <- ifelse(is.na(total_cholesterol),NA_integer_,as.integer(total_cholesterol>=200))
  dx <- binary(diagnosis);recommended<-binary(medication_recommended)
  med<-ifelse(recommended==0L,0L,ifelse(recommended==1L,binary(medication_current),NA_integer_))
  early<-cycle %in% c('E','F')
  dx_route<-early & ever_checked %in% 2 & is.na(diagnosis)
  med_route<-((early & diagnosis %in% 2) | ever_checked %in% 2) &
    is.na(medication_recommended) & is.na(medication_current)
  dx[dx_route]<-0L;med[med_route]<-0L
  x<-cbind(measured,dx,med);answer<-rep(NA_integer_,nrow(x))
  answer[rowSums(x==1,na.rm=TRUE)>0]<-1L
  answer[rowSums(is.na(x))==0 & rowSums(x==1,na.rm=TRUE)==0]<-0L
  list(diagnosed=dx,medication=med,high_cholesterol=answer,
       structural_diagnosis_skip=dx_route,structural_medication_skip=med_route)
}

cases <- data.frame(
 case=c('early_dx_no','early_never_checked','later_never_checked','later_checked_med_unknown',
 'early_refused_dx','later_unknown_test','high_lab_skip','known_dx_high',
 'med_no','med_yes','missing_lab_all_other_negative','later_never_checked_unknown_dx'),
 tc=c(rep(180,6),210,180,180,180,NA,180),
 dx=c(2,NA,2,2,7,2,2,1,2,2,2,9),
 rec=c(rep(NA,8),2,1,NA,NA),cur=c(rep(NA,9),1,NA,NA),
 checked=c(1,2,2,1,1,9,1,1,1,1,1,2),
 cycle=c('E','F','G','P','E','P','E','P','G','P','E','P'),
 expected=c(0,0,0,NA,NA,NA,1,1,0,1,NA,NA))
test_result <- with(cases,derive_cholesterol_routed(tc,dx,rec,cur,checked,cycle))
stopifnot(identical(as.integer(test_result$high_cholesterol),as.integer(cases$expected)))
)--------------------"

upf_source[["helpers/day1.R"]] <- r"--------------------(read_day1_weight <- function(cycle) {
  variable <- if (identical(cycle, "P")) "WTDRD1PP" else "WTDRD1"
  z <- read_xpt_required(diet_file(cycle, 1L, "tot"), c("SEQN", variable))
  names(z) <- c("SEQN", "day1_weight_original")
  z$cycle <- cycle
  z
}

derive_day1_exposure <- function(exposure) {
  required <- c(
    "SEQN", "cycle", "recall_status_d1", "total_energy_kcal_d1",
    paste0("nova", 1:4, "_g_d1"), "unclassified_g_d1",
    "total_nonwater_with_unclassified_g_d1", "classified_nonwater_g_d1",
    "no_iff_records_d1", "only_plain_water_d1", "zero_total_energy_d1"
  )
  assert_columns(exposure, required, "phase1_exposures Day-1 fields")
  z <- exposure[required]
  z$diet_two_day_reliable <- z$recall_status_d1 == 1
  z$unclassified_gram_share <- ifelse(
    z$total_nonwater_with_unclassified_g_d1 > 0,
    z$unclassified_g_d1 / z$total_nonwater_with_unclassified_g_d1,
    NA_real_
  )
  z$classification_quality_ok <- z$diet_two_day_reliable &
    is.finite(z$unclassified_gram_share) &
    z$unclassified_gram_share <= UNCLASSIFIED_GRAM_SHARE_LIMIT
  z$upf_gram_ratio_nonwater <- ifelse(
    z$classification_quality_ok & z$classified_nonwater_g_d1 > 0,
    100 * z$nova4_g_d1 / z$classified_nonwater_g_d1, NA_real_
  )
  z$total_energy_d1 <- as.numeric(z$total_energy_kcal_d1)
  z$all_day_fasting_d1 <- z$diet_two_day_reliable & z$zero_total_energy_d1 &
    z$no_iff_records_d1 & !z$only_plain_water_d1
  if (any(is.finite(z$upf_gram_ratio_nonwater) &
          (z$upf_gram_ratio_nonwater < -1e-8 |
             z$upf_gram_ratio_nonwater > 100 + 1e-8))) {
    stop("Day-1 UPF exposure falls outside 0-100")
  }
  z[c("SEQN", "cycle", "diet_two_day_reliable",
      "classification_quality_ok", "unclassified_gram_share",
      "total_energy_d1", "all_day_fasting_d1",
      "upf_gram_ratio_nonwater")]
}

build_day1_analysis_frame <- function() {
  contract <- read_phase1_variable_contract()
  raw <- rbind_fill(lapply(ANALYTIC_CYCLES, function(cycle) {
    standardize_cycle(read_cycle_person(cycle, contract))
  }))
  if (anyDuplicated(raw[c("SEQN", "cycle")])) stop("Duplicate participant-cycle key")
  weights <- rbind_fill(lapply(ANALYTIC_CYCLES, read_day1_weight))
  raw <- left_merge_unique(raw, weights)
  raw$analysis_weight <- derive_analysis_weight(
    raw$cycle, raw$day1_weight_original, raw$day1_weight_original
  )
  raw$mean_energy_kcal <- as.numeric(raw$DR1TKCAL)
  raw$diet_two_day_reliable_tot <- as.integer(raw$DR1DRSTZ == 1)
  raw <- derive_phase1_outcomes(raw)
  exposure_path <- file.path(DIRS$derived, "phase1_exposures.rds")
  if (!file.exists(exposure_path)) stop("Frozen phase1 exposure input is missing")
  day1 <- derive_day1_exposure(readRDS(exposure_path))
  raw <- left_merge_unique(raw, day1)
  raw <- create_domains(raw)
  raw
}

)--------------------"

upf_source[["helpers/prepare.R"]] <- r"--------------------(# Portable reconstruction of observed analytic frames from official survey files.
prepare_nhanes <- function() {
  code <- file.path(work,'NHANES/02_代码/正式五结局')
  source(file.path(code,'01_prepare_inputs.R'),local=FALSE)
  require_formal_authorization <- function(...) invisible(TRUE)
  stage_prepare_inputs()
  source(file.path(code,'02_build_dataset.R'),local=FALSE)
  stage_build_dataset()
  d <- readRDS(file.path(DIRS$derived,'analysis_frame_pre_mice.rds'))
  source(file.path(work,'helpers/cholesterol.R'),local=FALSE)
  checked <- rep(NA_real_,nrow(d))
  for(cy in ANALYTIC_CYCLES) {
    raw <- haven::read_xpt(module_file(cy,'Questionnaire','bpq'))
    ix <- which(d$cycle==cy);j<-match(d$SEQN[ix],raw$SEQN)
    if('BPQ060' %in% names(raw))checked[ix]<-as.numeric(raw$BPQ060[j])
  }
  fix_chol <- function(x) {
    z<-derive_cholesterol_routed(x$total_cholesterol,x$BPQ080,x$BPQ090D,x$BPQ100D,checked[match(x$SEQN,d$SEQN)],as.character(x$cycle))
    x$cholesterol_diagnosed<-z$diagnosed;x$cholesterol_med<-z$medication;x$high_cholesterol<-z$high_cholesterol;x
  }
  d<-fix_chol(d)
  saveRDS(d,file.path(DIRS$derived,'analysis_frame_pre_mice.rds'))
  saveRDS(d,file.path(work,'NHANES/03_派生数据/analysis_frame_pre_mice.rds'))
  source(file.path(work,'helpers/day1.R'),local=FALSE)
  day<-fix_chol(build_day1_analysis_frame())
  dir.create(file.path(work,'Day1/02_derived'),recursive=TRUE,showWarnings=FALSE)
  saveRDS(day,file.path(work,'Day1/02_derived/day1_analysis_frame_pre_mice.rds'))
  prepare_diet_scores(d)
  source(file.path(code,'03_imputation.R'),local=FALSE)
  six<-c('CKD','eGFR','UACR','DKD','kidney_stones','serum_uric_acid')
  OUTCOMES$serum_uric_acid<<-list(variable='uric_acid_dxc',family='gaussian',domain='uric_acid_observed',effect_scale='identity')
  manifest<-list()
  for(oc in six)for(branch in if(oc=='serum_uric_acid')'main' else c('main','day1')) {
    frame<-if(branch=='main')d else day
    EXPOSURES<<-list(G_R0=list(variable='upf_gram_ratio_nonwater'))
    if(branch=='main' && oc!='UACR') {
      EXPOSURES<<-list(E_R0=list(variable='upf_energy_ratio_2d'),G_R0=list(variable='upf_gram_ratio_nonwater'))
      ex<-readRDS(file.path(DIRS$derived,'phase1_exposures.rds'))
      frame$upf_energy_ratio_2d<-ex$upf_energy_ratio_2d[match(frame$SEQN,ex$SEQN)]
    }
    z<-prepare_mi_subframe(frame,oc,FALSE);cfg<-build_mice_spec(z$data,z$targets,z$fixed)
    mids<-mice::mice(z$data,m=1,maxit=0,method=cfg$method,predictorMatrix=cfg$predictor_matrix,printFlag=FALSE)
    obj<-list(mids=mids,original_row_index=z$row_index,outcome=oc,domain=OUTCOMES[[oc]]$domain,outcome_variable=OUTCOMES[[oc]]$variable,include_pir=FALSE)
    id<-paste(branch,oc,sep='_');tp<-file.path(work,'templates',paste0(id,'.rds'));dir.create(dirname(tp),recursive=TRUE,showWarnings=FALSE);saveRDS(obj,tp)
    fp<-if(branch=='main')file.path(work,'NHANES/03_派生数据/正式五结局_v1/analysis_frame_pre_mice.rds') else file.path(work,'Day1/02_derived/day1_analysis_frame_pre_mice.rds')
    dest<-if(branch=='day1')file.path(work,'Day1/02_derived',paste0('day1_mi_',oc,'.rds')) else if(oc=='serum_uric_acid')file.path(work,'NHANES/03_派生数据/mi_uric_acid.rds') else file.path(work,'NHANES/03_派生数据/正式五结局_v1',paste0('mi_',oc,'.rds'))
    manifest[[id]]<-data.frame(id=id,outcome=oc,branch=branch,template=tp,frame=fp,output=dest,n=length(z$row_index))
  }
  jsonlite::write_json(do.call(rbind,manifest),file.path(work,'03_qc/mi_job_manifest.json'),pretty=TRUE)
  cat('NHANES raw frames and observed-data MI templates constructed\n')
}

prepare_diet_scores <- function(formal) {
  source(file.path(work,'NHANES/02_代码/07_diet_quality.R'),local=FALSE)
  cycles<-c('E','F','G','H','I','P')
  fped<-lapply(cycles,function(cy){z<-combine_fped_two_day(read_fped_day(cy,1),read_fped_day(cy,2));z$cycle<-cy;z})
  nutrients<-lapply(cycles,function(cy){z<-read_two_day_nutrients(cy);z$cycle<-cy;z})
  scores<-build_pooled_diet_quality_scores(fped,nutrients,formal[,c('SEQN','sex','analysis_weight','base_domain')])
  ix<-match(formal$SEQN,scores$SEQN)
  for(v in setdiff(names(scores),'SEQN'))formal[[v]]<-scores[[v]][ix]
  e<-readRDS(file.path(work,'NHANES/03_派生数据/正式五结局_v1/phase1_exposures.rds'));j<-match(formal$SEQN,e$SEQN)
  for(k in 1:4) {
    v1<-paste0('nova',k,'_g_d1');v2<-paste0('nova',k,'_g_d2')
    formal[[paste0('nova',k,'_2d')]]<-e[[v1]][j]+e[[v2]][j]
  }
  formal$total_classified_nonwater_g<-rowSums(formal[paste0('nova',1:4,'_2d')])
  for(k in 1:4)formal[[paste0('nova',k,'_share')]]<-100*formal[[paste0('nova',k,'_2d')]]/formal$total_classified_nonwater_g
  for(k in 1:4)formal[[paste0('nova',k,'_share')]][!is.finite(formal$upf_gram_ratio_nonwater)]<-NA_real_
  dest<-file.path(work,'NHANES/03_派生数据/正式五结局_v1/三模块探索');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  saveRDS(formal,file.path(dest,'analysis_frame_three_modules.rds'));saveRDS(scores,file.path(dest,'diet_quality_scores.rds'))
}

prepare_knhanes <- function() {
  source(file.path(work,'KN/source_snapshot/04_build_analysis_frame.R'),local=FALSE)
  map<-read_adjudicated_nova_map(file.path(work,'KNHANES/02_数据配置/nova_context_adjudicated.csv'))
  ex<-dplyr::bind_rows(lapply(seq_len(nrow(YEAR_FILES)),function(i)collapse_nova_person_year(YEAR_FILES$year[i],YEAR_FILES$rc24[i],map)))
  people<-dplyr::bind_rows(lapply(seq_len(nrow(YEAR_FILES)),function(i)read_person_year(YEAR_FILES$year[i],YEAR_FILES$all[i])))
  d<-merge(people,ex,by=c('year','ID'),all.x=TRUE,sort=FALSE)
  d<-make_domains(derive_knhanes_outcomes(standardize_person_frame(d)))
  dir.create(file.path(work,'KN/02_derived'),recursive=TRUE,showWarnings=FALSE)
  saveRDS(d,file.path(work,'KN/02_derived/analysis_frame_before_design_filter.rds'))
  d<-d[is.finite(d$analysis_weight)&d$analysis_weight>0&!is.na(d$strata_uid)&!is.na(d$psu_uid),]
  extra<-dplyr::bind_rows(lapply(seq_len(nrow(YEAR_FILES)),function(i){x<-as.data.frame(haven::read_sas(YEAR_FILES$all[i]));cols<-intersect(c('ID','HE_fst','HEI','HE_Uacid','HE_Uacid_etc'),names(x));x<-x[,cols,drop=FALSE];x$year<-YEAR_FILES$year[i];x}))
  ix<-match(paste(d$year,d$ID),paste(extra$year,extra$ID));stopifnot(!anyNA(ix))
  for(v in setdiff(names(extra),c('ID','year')))d[[v]]<-extra[[v]][ix]
  valid<-is.finite(d$HE_glu)&is.finite(d$HE_fst)&d$HE_fst>=8
  fpg<-rep(NA_integer_,nrow(d));fpg[valid]<-as.integer(d$HE_glu[valid]>=126)
  a1c<-ifelse(is.finite(d$HE_HbA1c),as.integer(d$HE_HbA1c>=6.5),NA_integer_)
  d$diabetes_original<-or_composite(yes_no_na(d$DE1_dg),yes_no_na(d$DE1_31),yes_no_na(d$DE1_32),fpg,a1c)
  d$diabetes<-factor(ifelse(d$diabetes_original==1,'Yes',ifelse(d$diabetes_original==0,'No',NA_character_)),levels=c('No','Yes'))
  d$dkd<-ifelse(d$diabetes_original %in% 1,d$ckd,NA_integer_)
  d$dkd_alt_lod<-ifelse(d$diabetes_original==1L,d$ckd_alt_lod,NA_integer_)
  d$dkd_domain<-d$kidney_observed & d$diabetes_original==1L
  if(!'HE_Uacid_etc' %in% names(d))d$HE_Uacid_etc<-NA_character_
  ua<-apply_lab_limits_from_qualifier(as.numeric(d$HE_Uacid),d$HE_Uacid_etc)
  d$uric_acid_mg_dl<-ua$value_primary;d$uric_acid_mg_dl_alt_lod<-ua$value_alternative
  d$uric_acid_below_lod<-ua$below_lod;d$uric_acid_above_uloq<-ua$above_uloq;d$uric_acid_qualifier_boundary<-ua$qualifier_boundary
  d$khei<-as.numeric(d$HEI)
  fp<-file.path(work,'KN/02_derived/analysis_frame_fasting_corrected.rds');dir.create(dirname(fp),recursive=TRUE,showWarnings=FALSE);saveRDS(d,fp)
  source(file.path(work,'KN/source_snapshot/05_imputation.R'),local=FALSE)
  x<-build_mice_input(d);cfg<-configure_mice(x);saveRDS(list(data=x,config=cfg),file.path(work,'KN/02_derived/mice_input_and_config.rds'))
  cat('KNHANES raw frame and observed-data MI specification constructed\n')
}
)--------------------"

upf_source[["helpers/outputs.R"]] <- r"--------------------(# Figure generation uses R graphics only. Reference mode rebuilds the published
outcome_order <- c('CKD','eGFR','UACR','serum_uric_acid','DKD','kidney_stones')
outcome_names <- c(CKD='Chronic kidney disease',eGFR='Estimated glomerular filtration rate',UACR='Urinary albumin-to-creatinine ratio',serum_uric_acid='Serum uric acid',DKD='Diabetic kidney disease',kidney_stones='Kidney stones')
axis_names <- c(CKD='Odds ratio',eGFR='Difference in eGFR (mL/min/1.73 m²)',UACR='Relative change in UACR (%)',serum_uric_acid='Difference in serum uric acid (µmol/L)',DKD='Odds ratio',kidney_stones='Odds ratio')
is_binary <- function(o)o %in% c('CKD','DKD','kidney_stones')
cols <- c(NHANES='#C93636',KNHANES='#286AB0')
read_csv <- function(p)read.csv(p,check.names=FALSE,stringsAsFactors=FALSE)
bind <- function(...)dplyr::bind_rows(...)
pdf_open <- function(p,w,h)grDevices::cairo_pdf(p,width=w,height=h,family=Sys.getenv('UPF_PLOT_FONT',unset='Arial'),onefile=TRUE)

plot_rcs_panel <- function(d,h,oc,letter) {
  d<-d[d$outcome==oc,];h<-h[h$outcome==oc,]
  yr<-range(d$conf_low,d$conf_high,finite=TRUE);null<-if(is_binary(oc))1 else 0
  par(mar=c(3.1,3.5,2.1,.5),mgp=c(2,.5,0),tcl=-.2,cex=.7)
  plot(NA,xlim=c(0,75),ylim=yr,xlab='UPF intake (%)',ylab=axis_names[[oc]],bty='l')
  if(nrow(h))for(db in names(cols)) {
    hh<-h[h$survey==db,];if(!nrow(hh))next;ht<-hh$weight_percent/100*diff(yr)*1.5
    rect(hh$xmin,yr[1],hh$xmax,yr[1]+ht,col=adjustcolor(cols[[db]],alpha.f=.13),border=NA)
  }
  abline(h=null,lty=2,col='#777777')
  for(db in names(cols)) {
    z<-d[d$survey==db,];z<-z[order(z$x),];if(!nrow(z))next
    polygon(c(z$x,rev(z$x)),c(z$conf_low,rev(z$conf_high)),col=adjustcolor(cols[[db]],alpha.f=.15),border=NA)
    lines(z$x,z$estimate,col=cols[[db]],lwd=1.3)
  }
  mtext(letter,side=3,adj=0,line=.6,font=2,cex=1.1)
  legend('topleft',legend=intersect(names(cols),unique(d$survey)),col=cols[intersect(names(cols),unique(d$survey))],lty=1,bty='n',cex=.65)
}

forest <- function(d,path,order=outcome_order,width=10.8) {
  sensitivity<-basename(path)=='Figure_S1.pdf'
  d<-d[d$outcome %in% order,];d$outcome<-factor(d$outcome,levels=order)
  d<-d[order(d$outcome,d$row_order),]
  blocks<-split(d,d$outcome,drop=TRUE);nr<-sum(vapply(blocks,nrow,integer(1)))
  pdf_open(path,width,max(2.2,nr*.25+length(blocks)*.62+.65))
  on.exit(dev.off())
  par(mar=c(.2,.3,.2,.3),xaxs='i',yaxs='i',family=Sys.getenv('UPF_PLOT_FONT',unset='Arial'))
  total<-nr+length(blocks)*2.2+1.3
  plot.new();plot.window(xlim=c(0,1),ylim=c(total,0))
  text(c(.13,.30,.70,.875),.4,c(if(sensitivity)'Outcome and\nsensitivity analysis' else 'Model','Database','P value','95% CI'),font=2,cex=.85)
  segments(.01,.9,.99,.9);cursor<-1.15
  for(i in seq_along(blocks)) {
    z<-blocks[[i]];oc<-as.character(z$outcome[1]);n<-nrow(z)
    rect(.01,cursor,.99,cursor+.75,col='#EEF0F0',border=NA)
    text(.022,cursor+.37,outcome_names[[oc]],adj=0,font=2,cex=.85)
    start<-cursor+.75
    lim<-if(is_binary(oc))c(.94,1.21) else switch(oc,eGFR=c(-.65,.15),UACR=c(-1.3,3.1),serum_uric_acid=c(-.65,2.25))
    finite<-is.finite(z$low)&is.finite(z$high);lim<-range(lim,z$low[finite],z$high[finite])
    transform<-if(is_binary(oc))log else identity
    fx<-function(v).40+.235*(transform(v)-transform(lim[1]))/diff(transform(lim))
    null<-if(is_binary(oc))1 else 0
    for(j in seq_len(n)) {
      y<-start+j-.5
      if(j%%2==0)rect(.01,start+j-1,.99,start+j,col='#EEF0F0',border=NA)
      model_label<-if(sensitivity)paste(strwrap(z$model[j],width=23),collapse='\n') else z$model[j]
      text(.13,y,model_label,cex=.74);text(.30,y,z$survey[j],cex=.72)
      if(!finite[j]){text(.875,y,'Not assessed',cex=.74);next}
      color<-if(z$low[j]<=null&&z$high[j]>=null)'#999999' else cols[[z$survey[j]]]
      segments(fx(z$low[j]),y,fx(z$high[j]),y,col=color,lwd=1.3)
      shape<-if(sensitivity)if(z$survey[j]=='NHANES')15 else 16 else if(grepl('HEI',z$model[j]))16 else 18
      points(fx(z$estimate[j]),y,pch=shape,col=color,cex=.8)
      text(.70,y,if(z$p[j]<.001)'<0.001' else sprintf('%.3f',z$p[j]),cex=.74)
      fmt<-if(oc %in% c('UACR','serum_uric_acid'))'%.2f (%.2f, %.2f)' else '%.3f (%.3f, %.3f)'
      text(.875,y,sprintf(fmt,z$estimate[j],z$low[j],z$high[j]),cex=.74)
    }
    segments(fx(null),start,fx(null),start+n,col='#999999',lty=2)
    ticks<-if(is_binary(oc))c(.95,1,1.05,1.1,1.15,1.2) else pretty(lim,n=4)
    ticks<-ticks[ticks>=lim[1]&ticks<=lim[2]]
    segments(.40,start+n,.635,start+n)
    segments(fx(ticks),start+n,fx(ticks),start+n+.10)
    text(fx(ticks),start+n+.33,format(ticks,trim=TRUE),cex=.65)
    text(.5175,start+n+.82,axis_names[[oc]],cex=.73)
    cursor<-start+n+1.45
  }
}

collect_new <- function(work) {
  nh<-file.path(work,'NHANES/04_结果/正式五结局_v1');ua<-file.path(work,'SUA/03_扩展分析/NHANES')
  kn<-file.path(work,'KN/04_results/candidate/summary')
  nr<-file.path(work,'RCS/02_绘图数据_v3重算/统计重算源')
  a<-read_csv(file.path(nr,'NHANES_RCS_M3_v3_curves.csv'))
  curves<-data.frame(x=a$exposure_value,estimate=a$estimate,conf_low=a$conf_low,conf_high=a$conf_high,survey='NHANES',outcome=a$outcome)
  a<-read_csv(file.path(work,'SUA/final_alignment/RCS_aligned/NHANES_serum_uric_acid_RCS_curve.csv'))
  curves<-bind(curves,data.frame(x=a$upf_grams_pct,estimate=a$estimate*59.48,conf_low=a$ci_low*59.48,conf_high=a$ci_high*59.48,survey='NHANES',outcome='serum_uric_acid'))
  if(!any(curves$survey=='NHANES' & curves$outcome=='serum_uric_acid' & abs(curves$x-10)<1e-10))curves<-bind(curves,data.frame(x=10,estimate=0,conf_low=0,conf_high=0,survey='NHANES',outcome='serum_uric_acid'))
  curves<-bind(curves,read_csv(file.path(work,'RCS_low_extension_20260911/extension_display_values.csv')))
  a<-read_csv(file.path(kn,'KNHANES_RCS_curves.csv'));a$outcome[a$outcome=='log_uacr']<-'UACR';fac<-ifelse(a$outcome=='serum_uric_acid',59.48,1)
  curves<-bind(curves,data.frame(x=a$upf_grams_pct,estimate=a$display_estimate*fac,conf_low=a$display_low*fac,conf_high=a$display_high*fac,survey='KNHANES',outcome=a$outcome))
  hei<-bind(lapply(setdiff(outcome_order,'serum_uric_acid'),function(o)read_csv(file.path(nh,'三模块探索',paste0('diet_quality_UPF_',o,'.csv')))))
  a<-read_csv(file.path(ua,'HEI2015_adjustment.csv'));a$variant<-ifelse(a$model=='M3_same_sample','reference','adjusted');hei<-bind(hei,a)
  fac<-ifelse(hei$outcome=='serum_uric_acid',59.48,1)
  diet<-data.frame(outcome=hei$outcome,survey='NHANES',model=ifelse(hei$variant=='reference','Primary model','M3 + HEI'),estimate=hei$estimate*fac,low=hei$ci_low*fac,high=hei$ci_high*fac,p=hei$p_value,n=hei$n_unweighted,row_order=ifelse(hei$variant=='reference',1,2))
  a<-read_csv(file.path(kn,'KNHANES_KHEI_same_sample.csv'));a$outcome[a$outcome=='log_uacr']<-'UACR';fac<-ifelse(a$outcome=='serum_uric_acid',59.48,1)
  diet<-bind(diet,data.frame(outcome=a$outcome,survey='KNHANES',model=ifelse(a$comparison=='before','Primary model','M3 + KHEI'),estimate=a$display_estimate*fac,low=a$display_low*fac,high=a$display_high*fac,p=a$p_value,n=a$n,row_order=ifelse(a$comparison=='before',3,4)))
  main<-read_csv(file.path(work,'04_results/NHANES_main_M1_M4_iter20.csv'));main$q_BH<-NA_real_;i<-main$model=='M3';main$q_BH[i]<-p.adjust(main$p_value[i],'BH')
  list(curves=curves,diet=diet,main=main,kn=kn)
}

observed_display_data <- function(work) {
  nh<-readRDS(file.path(work,'NHANES/03_派生数据/正式五结局_v1/analysis_frame_pre_mice.rds'))
  kn<-readRDS(file.path(work,'KN/02_derived/analysis_frame_fasting_corrected.rds'))
  knall<-readRDS(file.path(work,'KN/02_derived/analysis_frame_before_design_filter.rds'))
  source(file.path(work,'GA/01_代码/10_ga_shared.R'),local=TRUE)
  gs<-hs<-fs<-list()
  for(db in c('NHANES','KNHANES')) {
    d<-if(db=='NHANES')nh else kn
    ok<-is.finite(d$analysis_weight)&d$analysis_weight>0
    if(db=='NHANES')ok<-ok&d$design_eligible
    d<-d[ok,];d<-derive_ga_fields(d)
    d$uric_acid_observed<-d$base_domain %in% TRUE & is.finite(d[[if(db=='NHANES')'uric_acid_dxc' else 'uric_acid_mg_dl']])
    ex<-d[[if(db=='NHANES')'upf_gram_ratio_nonwater' else 'upf_grams_pct']]
    dom<-c(CKD='kidney_observed',eGFR='egfr_observed',UACR='uacr_observed',serum_uric_acid='uric_acid_observed',DKD='dkd_domain',kidney_stones='stone_observed')
    for(oc in outcome_order) {
      if(db=='KNHANES'&&oc=='kidney_stones')next
      use<-d$base_domain %in% TRUE & d[[dom[[oc]]]] %in% TRUE & is.finite(ex)
      bin<-cut(pmin(100,pmax(0,ex)),seq(0,100,5),include.lowest=TRUE,right=FALSE,labels=FALSE)
      for(j in 1:20)hs[[length(hs)+1]]<-data.frame(survey=db,outcome=oc,xmin=(j-1)*5,xmax=j*5,weight_percent=100*sum(d$analysis_weight[use & bin %in% j])/sum(d$analysis_weight[use]))
      fs[[length(fs)+1]]<-data.frame(database=db,node_type='outcome',label=outcome_names[[oc]],outcome=oc,remaining=sum(use))
    }
    g<-d[d$base_domain %in% TRUE & d$ga_domain,]
    for(gg in G_LEVELS)for(aa in A_LEVELS) {
      ii<-g$g_category==gg & g$a_category==aa
      gs[[length(gs)+1]]<-data.frame(database=db,group='Overall',g_category=gg,a_category=aa,risk=unname(GA_MAP[paste(gg,aa,sep='_')]),n_unweighted=sum(ii),weighted_proportion=sum(g$analysis_weight[ii])/sum(g$analysis_weight))
    }
    a<-if(db=='NHANES')nh else knall
    x<-a[[if(db=='NHANES')'upf_gram_ratio_nonwater' else 'upf_grams_pct']]
    u<-is.finite(x);adult<-u & a$age>=20
    np<-adult & if(db=='NHANES')!(a$nonpregnant %in% 0) else a$pregnant %in% FALSE
    counts<-c(nrow(a),sum(u,na.rm=TRUE),sum(adult,na.rm=TRUE),sum(np,na.rm=TRUE),sum(d$base_domain %in% TRUE))
    fs[[length(fs)+1]]<-data.frame(database=db,node_type=c('initial',rep('selection',4)),label=c('Initial study sample','Usable NOVA exposure data','Adults aged 20 years or older','After excluding pregnant participants','Common analytic base'),outcome='',remaining=counts)
  }
  fs[[length(fs)+1]]<-data.frame(database='KNHANES',node_type='outcome',label='Kidney stones',outcome='kidney_stones',remaining=NA_real_)
  list(ga=bind(gs),hist=bind(hs),flow=bind(fs))
}

new_sensitivity <- function(work) {
  p<-list.files(file.path(work,'NHANES/04_结果/正式五结局_v1'),pattern='^sensitivity_M3_all_outcomes.csv$',recursive=TRUE,full.names=TRUE);stopifnot(length(p)==1)
  a<-bind(read_csv(p),read_csv(file.path(work,'SUA/03_扩展分析/NHANES/sensitivity_M3.csv')))
  lev<-c('primary_MICE','complete_case','exclude_extreme_energy','trim_G_R0_p01_p99','Day1')
  labs<-c('Primary model (multiple imputation)','Complete-case analysis','Energy-intake restriction','UPF restricted to weighted P1-P99','Day-1 dietary recall and weights')
  day<-read_csv(file.path(work,'Day1/04_results/day1_M3.csv'));day$sensitivity<-'Day1';a<-bind(a,day)
  fac<-ifelse(a$outcome=='serum_uric_acid',59.48,1)
  d<-data.frame(outcome=a$outcome,survey='NHANES',model=labs[match(a$sensitivity,lev)],estimate=a$estimate*fac,low=a$ci_low*fac,high=a$ci_high*fac,p=a$p_value,row_order=2*match(a$sensitivity,lev)-1)
  k<-file.path(work,'KN/04_results/candidate/summary')
  a<-read_csv(file.path(k,'KNHANES_primary_M3_BH5.csv'));a$scenario<-'primary'
  a<-bind(a,read_csv(file.path(k,'KNHANES_sensitivity.csv')));a$outcome[a$outcome=='log_uacr']<-'UACR'
  fac<-ifelse(a$outcome=='serum_uric_acid',59.48,1);ix<-match(a$scenario,c('primary','complete_case','energy','upf_p01_p99'))
  stopifnot(!anyNA(ix),!anyNA(d$model))
  bind(d,data.frame(outcome=a$outcome,survey='KNHANES',model=labs[ix],estimate=a$display_estimate*fac,low=a$display_low*fac,high=a$display_high*fac,p=a$p_value,row_order=ix*2),data.frame(outcome='serum_uric_acid',survey='NHANES',model=labs[5],estimate=NA_real_,low=NA_real_,high=NA_real_,p=NA_real_,row_order=9))
}

cross_survey_tables <- function(work,dest) {
  a<-read_csv(file.path(work,'04_results/NHANES_main_M1_M4_iter20.csv'))
  a<-a[a$outcome!='kidney_stones',c('outcome','model','estimate_link','standard_error')];names(a)[3:4]<-c('beta_nh','se_nh')
  k<-file.path(work,'KN/04_results/candidate/summary')
  b<-read_csv(file.path(k,'KNHANES_main_all_exposure_terms.csv'));b<-b[b$exposure_form=='continuous' & b$outcome!='hyperuricemia',]
  b$outcome[b$outcome=='log_uacr']<-'UACR';b$beta_kn<-b$estimate;is_or<-b$effect_scale=='OR';b$beta_kn[is_or]<-log(b$estimate[is_or]);b$se_kn<-b$std_error_link
  z<-merge(a,b[,c('outcome','model','beta_kn','se_kn')],by=c('outcome','model'));stopifnot(nrow(z)==20)
  z$difference_link<-z$beta_nh-z$beta_kn;z$se_difference<-sqrt(z$se_nh^2+z$se_kn^2)
  z$p_value<-2*pnorm(abs(z$difference_link/z$se_difference),lower.tail=FALSE);z$q_BH5<-ave(z$p_value,z$model,FUN=function(x)p.adjust(x,'BH'))
  z$estimate<-z$difference_link;z$ci_low<-z$difference_link-qnorm(.975)*z$se_difference;z$ci_high<-z$difference_link+qnorm(.975)*z$se_difference
  binary<-z$outcome %in% c('CKD','DKD','UACR')
  for(v in c('estimate','ci_low','ci_high')){z[[v]][binary]<-exp(z[[v]][binary]);z[[v]][z$outcome=='serum_uric_acid']<-z[[v]][z$outcome=='serum_uric_acid']*59.48}
  z$scale<-ifelse(z$outcome %in% c('CKD','DKD'),'ratio of odds ratios',ifelse(z$outcome=='UACR','ratio of multiplicative effects','NHANES minus KNHANES difference'))
  write.csv(z,file.path(dest,'Table_S9_cross_survey_M1_M4.csv'),row.names=FALSE)
  a<-read_csv(file.path(work,'GA/02_结果/G_A/NHANES/ga_proportional_odds_M1_M4.csv'));b<-read_csv(file.path(k,'KNHANES_GA_M1_M4.csv'))
  a<-a[a$model=='M3',];b<-b[b$model=='M3',];delta<-a$estimate_link-b$estimate_link;se<-sqrt(a$standard_error_link^2+b$standard_error_link^2)
  write.csv(data.frame(outcome='G_A',estimate=exp(delta),ci_low=exp(delta-qnorm(.975)*se),ci_high=exp(delta+qnorm(.975)*se),p_value=2*pnorm(abs(delta/se),lower.tail=FALSE),multiplicity='exploratory, outside primary outcome family'),file.path(dest,'Table_S9_GA_exploratory.csv'),row.names=FALSE)
}

generate_outputs <- function(package_root,work) {
  dest<-file.path(work,'regenerated');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
    z<-collect_new(work);curves<-z$curves;diet<-z$diet
    z$main<-z$main[order(match(z$main$outcome,outcome_order),match(z$main$model,c('M1','M2','M3','M4'))),]
    write.csv(z$main,file.path(dest,'Table2_NHANES_M1_M4.csv'),row.names=FALSE)
    cross_survey_tables(work,dest)
    allcsv<-list.files(work,pattern='\\.csv$',recursive=TRUE,full.names=TRUE)
    allcsv<-allcsv[!grepl('/regenerated/|/reference_regenerated/|/inputs/',allcsv)]
    paths<-allcsv[grepl('04_results|04_结果|04_summary|02_结果',allcsv)]
    for(p in paths){rel<-substring(p,nchar(work)+2);q<-file.path(dest,'tables',rel);dir.create(dirname(q),recursive=TRUE,showWarnings=FALSE);file.copy(p,q,overwrite=TRUE)}

    obs<-observed_display_data(work);hist<-obs$hist;ga<-obs$ga;f<-obs$flow
    write.csv(ga,file.path(dest,'Figure3_GA.csv'),row.names=FALSE);write.csv(f,file.path(dest,'Figure2_flow.csv'),row.names=FALSE)
  pdf_open(file.path(dest,'Figure_4.pdf'),11,7.5);par(mfrow=c(2,3))
  for(i in seq_along(outcome_order))plot_rcs_panel(curves,hist,outcome_order[i],LETTERS[i]);dev.off()
  for(i in seq_along(outcome_order)) {o<-outcome_order[i];pdf_open(file.path(dest,paste0('Figure4_',LETTERS[i],'_',o,'.pdf')),4.5,4);plot_rcs_panel(curves,hist,o,LETTERS[i]);dev.off()}
  forest(diet,file.path(dest,'Figure_5.pdf'))
  for(o in outcome_order)forest(diet[diet$outcome==o,],file.path(dest,paste0('Figure5_',o,'.pdf')),order=o)
  write.csv(curves,file.path(dest,'Figure4_plot_data.csv'),row.names=FALSE);write.csv(diet,file.path(dest,'Figure5_plot_data.csv'),row.names=FALSE)
    sens<-new_sensitivity(work);forest(sens,file.path(dest,'Figure_S1.pdf'));write.csv(sens,file.path(dest,'FigureS1_plot_data.csv'),row.names=FALSE)
  ga<-ga[ga$group=='Overall',]
  pdf_open(file.path(dest,'Figure_3.pdf'),9,4.5);par(mfrow=c(1,2),mar=c(3,4,2,1))
  glevels<-c('G1','G2','G3a','G3b','G4','G5');riskcols<-c(Low='#71BE70',Moderate='#F2D666',High='#EDA648',Very_high='#D85454','Very high'='#D85454')
  for(db in names(cols)) {
    plot.new();plot.window(xlim=c(0,3),ylim=c(6,0));z<-ga[ga$database==db,]
    for(i in seq_len(nrow(z))) {g<-match(z$g_category[i],glevels);a<-match(z$a_category[i],c('A1','A2','A3'));co<-riskcols[z$risk[i]];if(is.na(co))co<-'#D85454';rect(a-1,g-1,a,g,col=co,border='white');text(a-.5,g-.5,sprintf('%.2f%%\nn=%s',100*z$weighted_proportion[i],z$n_unweighted[i]),cex=.7)}
    axis(2,at=seq(.5,5.5),labels=glevels,las=1,tick=FALSE);axis(1,at=c(.5,1.5,2.5),labels=c('A1','A2','A3'),tick=FALSE);title(db)
  };dev.off()
  f$label<-gsub('Unusable dietary/NOVA data|Missing NOVA data','Unusable dietary or NOVA classification data',f$label)
  f$label<-gsub('Usable dietary/NOVA data|Participants with available NOVA data','Usable dietary and NOVA data',f$label)
  pdf_open(file.path(dest,'Figure_2.pdf'),8,9)
  par(mar=rep(.5,4));plot.new();plot.window(xlim=c(0,1),ylim=c(1,0))
  for(i in 1:2) {db<-names(cols)[i];z<-f[f$database==db&f$node_type!='outcome',];x<-c(.25,.75)[i];ys<-seq(.07,.65,length.out=nrow(z));
    for(j in seq_len(nrow(z))) {y<-ys[j];rect(x-.22,y-.03,x+.22,y+.03,border=cols[[db]],col='white');text(x,y,paste0(paste(strwrap(z$label[j],width=34),collapse='\n'),'\n(n=',format(z$remaining[j],big.mark=','),')'),cex=.65);if(j<nrow(z))arrows(x,y+.03,x,ys[j+1]-.03,length=.05)}
    z<-f[f$database==db&f$node_type=='outcome',];if(nrow(z))text(x,.84,paste(paste0(z$outcome,': n=',format(z$remaining,big.mark=',')),collapse='\n'),cex=.7)
  };dev.off()
  writeLines('Results generated from the analysis.',file.path(dest,'RUN_MODE.txt'))
}
)--------------------"

embedded_config <- list()

embedded_config[["KNHANES_outcome_definitions.csv"]] <- r"--------------------(outcome_id,label,type,source_variables,formula_or_threshold,missing_logic,definition_version,reference_url,reference_note
eGFR,eGFR 2021 CKD-EPI,continuous,"HE_crea;age;sex","142*min(Scr/k,1)^alpha*max(Scr/k,1)^-1.200*0.9938^age*(1.012 if female); k=0.7 female/0.9 male; alpha=-0.241 female/-0.302 male","NA if Scr missing/nonpositive, age missing/<18, or sex missing",OUTCOME_DEFINITIONS_20260820_v1,https://doi.org/10.1056/NEJMoa2102953,"Inker et al. 2021 race-free creatinine equation"
log_uacr,Natural log of urine albumin-to-creatinine ratio,continuous,"HE_Ualb;HE_Ucrea","uACR=100*urine albumin ug/mL / urine creatinine mg/dL; log_uacr=ln(uACR); restricted to uACR>0","NA if uACR not calculable or uACR<=0",OUTCOME_DEFINITIONS_20260820_v1,https://kdigo.org/wp-content/uploads/2017/02/KDIGO_2012_CKD_GL.pdf,"Continuous uACR on the log scale aligns with the NHANES five-outcome protocol"
CKD,Chronic kidney disease phenotype,binary,eGFR;uACR,"positive if eGFR<60 or uACR>=30; negative only if both observed with eGFR>=60 and uACR<30","three-state: positive if either criterion observed positive; negative only if both observed negative; otherwise NA",OUTCOME_DEFINITIONS_20260820_v1,https://kdigo.org/wp-content/uploads/2017/02/KDIGO_2012_CKD_GL.pdf,"Cross-sectional KNHANES phenotype cannot establish persistence for >3 months"
DKD,Diabetic kidney disease phenotype,binary,diabetes_original;CKD,"CKD phenotype among participants with diabetes","NA outside diabetes or when CKD is indeterminate",OUTCOME_DEFINITIONS_20260820_v1,https://doi.org/10.2337/dc26-s011,"diabetes_original is the composite definition (dx OR insulin OR oral OR HbA1c>=6.5); cross-sectional epidemiologic DKD phenotype, not an etiologic clinical diagnosis"
)--------------------"

embedded_config[["NHANES_config/diet_quality_manifest.csv"]] <- r"--------------------(record_type,score,component,cycle,day,relative_path,sha256,input_mapping,two_day_rule,scoring_rule,source_url,configuration_status,version

data_file,HEI2015_and_aMED,FPED_person_day,E,1,FPED/fped_dr1tot_0708.sas7bdat,497ea6463df07badffee783ddacfde5d7397be0e6e31387bf228cebc37a4c28e,USDA FPED person-day totals,require both days then sum component equivalents before scoring,not_applicable,https://www.ars.usda.gov/northeast-area/beltsville-md-bhnrc/beltsville-human-nutrition-research-center/food-surveys-research-group/docs/fped-databases/,final,2026-08-14.rev22-confirmed-plan

data_file,HEI2015_and_aMED,FPED_person_day,E,2,FPED/fped_dr2tot_0708.sas7bdat,d27d9b58beca469f3a76b5c5aa1a34e7069a4203984b32a39b568a43f040ff28,USDA FPED person-day totals,require both days then sum component equivalents before scoring,not_applicable,https://www.ars.usda.gov/northeast-area/beltsville-md-bhnrc/beltsville-human-nutrition-research-center/food-surveys-research-group/docs/fped-databases/,final,2026-08-14.rev22-confirmed-plan

data_file,HEI2015_and_aMED,FPED_person_day,F,1,FPED/fped_dr1tot_0910.sas7bdat,f9ffae5beaa00b81939ac54f8cc431b3e7a5e9986828ff87482c3ce8c0da91de,USDA FPED person-day totals,require both days then sum component equivalents before scoring,not_applicable,https://www.ars.usda.gov/northeast-area/beltsville-md-bhnrc/beltsville-human-nutrition-research-center/food-surveys-research-group/docs/fped-databases/,final,2026-08-14.rev22-confirmed-plan

data_file,HEI2015_and_aMED,FPED_person_day,F,2,FPED/fped_dr2tot_0910.sas7bdat,ae0f926c4e0078847ec19e82240dcc91e201cd1fac9df4188731cd4360cf045b,USDA FPED person-day totals,require both days then sum component equivalents before scoring,not_applicable,https://www.ars.usda.gov/northeast-area/beltsville-md-bhnrc/beltsville-human-nutrition-research-center/food-surveys-research-group/docs/fped-databases/,final,2026-08-14.rev22-confirmed-plan

data_file,HEI2015_and_aMED,FPED_person_day,G,1,FPED/fped_dr1tot_1112.sas7bdat,fd96b50d15645c0c88992fa5ab90a3c1d62cc04e910ef11943161aec5bd7bd65,USDA FPED person-day totals,require both days then sum component equivalents before scoring,not_applicable,https://www.ars.usda.gov/northeast-area/beltsville-md-bhnrc/beltsville-human-nutrition-research-center/food-surveys-research-group/docs/fped-databases/,final,2026-08-14.rev22-confirmed-plan

data_file,HEI2015_and_aMED,FPED_person_day,G,2,FPED/fped_dr2tot_1112.sas7bdat,e806cd919932c0d43800b1204c98c9c8827aab09fc582c45ed688d3b3a07e7d7,USDA FPED person-day totals,require both days then sum component equivalents before scoring,not_applicable,https://www.ars.usda.gov/northeast-area/beltsville-md-bhnrc/beltsville-human-nutrition-research-center/food-surveys-research-group/docs/fped-databases/,final,2026-08-14.rev22-confirmed-plan

data_file,HEI2015_and_aMED,FPED_person_day,H,1,FPED/fped_dr1tot_1314.sas7bdat,28acabffb029858e0e9f4da2f5a39d11e8fc573f0aa7e1ad29c32643ce7f7f57,USDA FPED person-day totals,require both days then sum component equivalents before scoring,not_applicable,https://www.ars.usda.gov/northeast-area/beltsville-md-bhnrc/beltsville-human-nutrition-research-center/food-surveys-research-group/docs/fped-databases/,final,2026-08-14.rev22-confirmed-plan

data_file,HEI2015_and_aMED,FPED_person_day,H,2,FPED/fped_dr2tot_1314.sas7bdat,3fcd28ecb3975572c5143cbcb544f4a9a06b98fb47e31b81c1991611207d2016,USDA FPED person-day totals,require both days then sum component equivalents before scoring,not_applicable,https://www.ars.usda.gov/northeast-area/beltsville-md-bhnrc/beltsville-human-nutrition-research-center/food-surveys-research-group/docs/fped-databases/,final,2026-08-14.rev22-confirmed-plan

data_file,HEI2015_and_aMED,FPED_person_day,I,1,FPED/fped_dr1tot_1516.sas7bdat,0e96dfe65c4065b6ee4ebd1a8ee85a0b2494eeb8a002c54723491a31309e7d67,USDA FPED person-day totals,require both days then sum component equivalents before scoring,not_applicable,https://www.ars.usda.gov/northeast-area/beltsville-md-bhnrc/beltsville-human-nutrition-research-center/food-surveys-research-group/docs/fped-databases/,final,2026-08-14.rev22-confirmed-plan

data_file,HEI2015_and_aMED,FPED_person_day,I,2,FPED/fped_dr2tot_1516.sas7bdat,4f497160abece25e9a60542623393c3a6643b6cf38fa60110d50dca814e9f801,USDA FPED person-day totals,require both days then sum component equivalents before scoring,not_applicable,https://www.ars.usda.gov/northeast-area/beltsville-md-bhnrc/beltsville-human-nutrition-research-center/food-surveys-research-group/docs/fped-databases/,final,2026-08-14.rev22-confirmed-plan

data_file,HEI2015_and_aMED,FPED_person_day,P,1,FPED/fped_dr1tot_1720.sas7bdat,d09b477609adc1b7aeeb130929396ae80dc3d6224ff9655cfa01076503196ad6,USDA FPED person-day totals,require both days then sum component equivalents before scoring,not_applicable,https://www.ars.usda.gov/northeast-area/beltsville-md-bhnrc/beltsville-human-nutrition-research-center/food-surveys-research-group/docs/fped-databases/,final,2026-08-14.rev22-confirmed-plan

data_file,HEI2015_and_aMED,FPED_person_day,P,2,FPED/fped_dr2tot_1720.sas7bdat,6695f75e78eb38b371fc6ae081d1dc42d7b55ccc414ae112f4c54826c0260b2b,USDA FPED person-day totals,require both days then sum component equivalents before scoring,not_applicable,https://www.ars.usda.gov/northeast-area/beltsville-md-bhnrc/beltsville-human-nutrition-research-center/food-surveys-research-group/docs/fped-databases/,final,2026-08-14.rev22-confirmed-plan

scoring_component,HEI2015,total_fruits,ALL,1+2,,,FP_F_TOTAL,sum both recalls before densities and component scoring,adequacy 5 points; max at 0.8 cup eq/1000 kcal,https://epi.grants.cancer.gov/hei/calculating-hei-scores.html,final,2026-08-14.rev22-confirmed-plan

scoring_component,HEI2015,whole_fruits,ALL,1+2,,,FP_F_CITMLB + FP_F_OTHER,sum both recalls before densities and component scoring,adequacy 5 points; max at 0.4 cup eq/1000 kcal,https://epi.grants.cancer.gov/hei/calculating-hei-scores.html,final,2026-08-14.rev22-confirmed-plan

scoring_component,HEI2015,total_vegetables,ALL,1+2,,,FP_V_TOTAL + FP_V_LEGUMES,sum both recalls before densities and component scoring,adequacy 5 points; max at 1.1 cup eq/1000 kcal,https://epi.grants.cancer.gov/hei/calculating-hei-scores.html,final,2026-08-14.rev22-confirmed-plan

scoring_component,HEI2015,greens_beans,ALL,1+2,,,FP_V_DRKGR + FP_V_LEGUMES,sum both recalls before densities and component scoring,adequacy 5 points; max at 0.2 cup eq/1000 kcal,https://epi.grants.cancer.gov/hei/calculating-hei-scores.html,final,2026-08-14.rev22-confirmed-plan

scoring_component,HEI2015,whole_grains,ALL,1+2,,,FP_G_WHOLE,sum both recalls before densities and component scoring,adequacy 10 points; max at 1.5 oz eq/1000 kcal,https://epi.grants.cancer.gov/hei/calculating-hei-scores.html,final,2026-08-14.rev22-confirmed-plan

scoring_component,HEI2015,dairy,ALL,1+2,,,FP_D_TOTAL,sum both recalls before densities and component scoring,adequacy 10 points; max at 1.3 cup eq/1000 kcal,https://epi.grants.cancer.gov/hei/calculating-hei-scores.html,final,2026-08-14.rev22-confirmed-plan

scoring_component,HEI2015,total_protein,ALL,1+2,,,FP_PF_MEAT + FP_PF_CUREDMEAT + FP_PF_ORGAN + FP_PF_POULT + FP_PF_SEAFD_HI + FP_PF_SEAFD_LOW + FP_PF_EGGS + FP_PF_SOY + FP_PF_NUTSDS + FP_PF_LEGUMES,sum both recalls before densities and component scoring,adequacy 5 points; max at 2.5 oz eq/1000 kcal,https://epi.grants.cancer.gov/hei/calculating-hei-scores.html,final,2026-08-14.rev22-confirmed-plan

scoring_component,HEI2015,seafood_plant,ALL,1+2,,,FP_PF_SEAFD_HI + FP_PF_SEAFD_LOW + FP_PF_SOY + FP_PF_NUTSDS + FP_PF_LEGUMES,sum both recalls before densities and component scoring,adequacy 5 points; max at 0.8 oz eq/1000 kcal,https://epi.grants.cancer.gov/hei/calculating-hei-scores.html,final,2026-08-14.rev22-confirmed-plan

scoring_component,HEI2015,fatty_acids,ALL,1+2,,,(DR1TMFAT+DR2TMFAT+DR1TPFAT+DR2TPFAT)/(DR1TSFAT+DR2TSFAT),sum both recalls before densities and component scoring,adequacy 10 points; ratio 1.2=0 and 2.5=10,https://epi.grants.cancer.gov/hei/calculating-hei-scores.html,final,2026-08-14.rev22-confirmed-plan

scoring_component,HEI2015,refined_grains,ALL,1+2,,,FP_G_REFINED,sum both recalls before densities and component scoring,moderation 10 points; 1.8=10 and 4.3=0 oz eq/1000 kcal,https://epi.grants.cancer.gov/hei/calculating-hei-scores.html,final,2026-08-14.rev22-confirmed-plan

scoring_component,HEI2015,sodium,ALL,1+2,,,DR1TSODI + DR2TSODI,sum both recalls before densities and component scoring,moderation 10 points; 1.1=10 and 2.0=0 g/1000 kcal,https://epi.grants.cancer.gov/hei/calculating-hei-scores.html,final,2026-08-14.rev22-confirmed-plan

scoring_component,HEI2015,added_sugars,ALL,1+2,,,FP_ADD_SUGARS,sum both recalls before densities and component scoring,moderation 10 points; 6.5%=10 and 26%=0 energy,https://epi.grants.cancer.gov/hei/calculating-hei-scores.html,final,2026-08-14.rev22-confirmed-plan

scoring_component,HEI2015,saturated_fats,ALL,1+2,,,DR1TSFAT + DR2TSFAT,sum both recalls before densities and component scoring,moderation 10 points; 8%=10 and 16%=0 energy,https://epi.grants.cancer.gov/hei/calculating-hei-scores.html,final,2026-08-14.rev22-confirmed-plan

scoring_component,aMED_9,vegetables,ALL,1+2,,,"max((FP_V_TOTAL - FP_V_STARCHY_POTATO)/2,0)",sum both recalls then divide by 2 for average daily intake,1 if above sex-specific survey-weighted median,00_计划文档/最终研究计划_冻结版_20260813.md,final,2026-08-14.rev22-confirmed-plan

scoring_component,aMED_9,legumes,ALL,1+2,,,FP_V_LEGUMES/2,sum both recalls then divide by 2 for average daily intake,1 if above sex-specific survey-weighted median,00_计划文档/最终研究计划_冻结版_20260813.md,final,2026-08-14.rev22-confirmed-plan

scoring_component,aMED_9,fruits,ALL,1+2,,,FP_F_TOTAL/2,sum both recalls then divide by 2 for average daily intake,1 if above sex-specific survey-weighted median,00_计划文档/最终研究计划_冻结版_20260813.md,final,2026-08-14.rev22-confirmed-plan

scoring_component,aMED_9,nuts,ALL,1+2,,,FP_PF_NUTSDS/2,sum both recalls then divide by 2 for average daily intake,1 if above sex-specific survey-weighted median,00_计划文档/最终研究计划_冻结版_20260813.md,final,2026-08-14.rev22-confirmed-plan

scoring_component,aMED_9,whole_grains,ALL,1+2,,,FP_G_WHOLE/2,sum both recalls then divide by 2 for average daily intake,1 if above sex-specific survey-weighted median,00_计划文档/最终研究计划_冻结版_20260813.md,final,2026-08-14.rev22-confirmed-plan

scoring_component,aMED_9,fish,ALL,1+2,,,(FP_PF_SEAFD_HI + FP_PF_SEAFD_LOW)/2,sum both recalls then divide by 2 for average daily intake,1 if above sex-specific survey-weighted median,00_计划文档/最终研究计划_冻结版_20260813.md,final,2026-08-14.rev22-confirmed-plan

scoring_component,aMED_9,mufa_sfa,ALL,1+2,,,(two-day MUFA/2)/(two-day SFA/2),sum both recalls then divide by 2 for average daily intake,1 if above sex-specific survey-weighted median,00_计划文档/最终研究计划_冻结版_20260813.md,final,2026-08-14.rev22-confirmed-plan

scoring_component,aMED_9,red_processed_meat,ALL,1+2,,,(FP_PF_MEAT + FP_PF_CUREDMEAT)/2,sum both recalls then divide by 2 for average daily intake,1 if below sex-specific survey-weighted median,00_计划文档/最终研究计划_冻结版_20260813.md,final,2026-08-14.rev22-confirmed-plan

scoring_component,aMED_9,alcohol,ALL,1+2,,,(DR1TALCO + DR2TALCO)/2,sum both recalls then divide by 2 for average daily intake,1 for women 5-15 g/day or men 10-25 g/day,00_计划文档/最终研究计划_冻结版_20260813.md,final,2026-08-14.rev22-confirmed-plan

)--------------------"

embedded_config[["NHANES_config/laboratory_bridge.csv"]] <- r"--------------------(analyte,cycle,source_variable,source_scale,target_variable,target_scale,conversion_direction,equation,source_unit,target_unit,applicable_cycles,official_url,qa_min_exclusive,qa_max,qa_action,configuration_status,version

serum_creatinine,E,LBXSCR,released_value_on_pre_2017_or_DxC_comparable_scale,creatinine_dxc,2015_2016_DxC_660i,identity,source_value,mg/dL,mg/dL,E,https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BIOPRO_J.htm,0,,flag_nonpositive_or_nonfinite;do_not_silently_drop,final,2026-08-14.rev22-confirmed-plan

serum_uric_acid,E,LBXSUA,released_value_on_pre_2017_or_DxC_comparable_scale,uric_acid_dxc,2015_2016_DxC_660i,identity,source_value,mg/dL,mg/dL,E,https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BIOPRO_J.htm,0,,flag_nonpositive_or_nonfinite;do_not_silently_drop,final,2026-08-14.rev22-confirmed-plan

serum_creatinine,F,LBXSCR,released_value_on_pre_2017_or_DxC_comparable_scale,creatinine_dxc,2015_2016_DxC_660i,identity,source_value,mg/dL,mg/dL,F,https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BIOPRO_J.htm,0,,flag_nonpositive_or_nonfinite;do_not_silently_drop,final,2026-08-14.rev22-confirmed-plan

serum_uric_acid,F,LBXSUA,released_value_on_pre_2017_or_DxC_comparable_scale,uric_acid_dxc,2015_2016_DxC_660i,identity,source_value,mg/dL,mg/dL,F,https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BIOPRO_J.htm,0,,flag_nonpositive_or_nonfinite;do_not_silently_drop,final,2026-08-14.rev22-confirmed-plan

serum_creatinine,G,LBXSCR,released_value_on_pre_2017_or_DxC_comparable_scale,creatinine_dxc,2015_2016_DxC_660i,identity,source_value,mg/dL,mg/dL,G,https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BIOPRO_J.htm,0,,flag_nonpositive_or_nonfinite;do_not_silently_drop,final,2026-08-14.rev22-confirmed-plan

serum_uric_acid,G,LBXSUA,released_value_on_pre_2017_or_DxC_comparable_scale,uric_acid_dxc,2015_2016_DxC_660i,identity,source_value,mg/dL,mg/dL,G,https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BIOPRO_J.htm,0,,flag_nonpositive_or_nonfinite;do_not_silently_drop,final,2026-08-14.rev22-confirmed-plan

serum_creatinine,H,LBXSCR,released_value_on_pre_2017_or_DxC_comparable_scale,creatinine_dxc,2015_2016_DxC_660i,identity,source_value,mg/dL,mg/dL,H,https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BIOPRO_J.htm,0,,flag_nonpositive_or_nonfinite;do_not_silently_drop,final,2026-08-14.rev22-confirmed-plan

serum_uric_acid,H,LBXSUA,released_value_on_pre_2017_or_DxC_comparable_scale,uric_acid_dxc,2015_2016_DxC_660i,identity,source_value,mg/dL,mg/dL,H,https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BIOPRO_J.htm,0,,flag_nonpositive_or_nonfinite;do_not_silently_drop,final,2026-08-14.rev22-confirmed-plan

serum_creatinine,I,LBXSCR,released_value_on_pre_2017_or_DxC_comparable_scale,creatinine_dxc,2015_2016_DxC_660i,identity,source_value,mg/dL,mg/dL,I,https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BIOPRO_J.htm,0,,flag_nonpositive_or_nonfinite;do_not_silently_drop,final,2026-08-14.rev22-confirmed-plan

serum_uric_acid,I,LBXSUA,released_value_on_pre_2017_or_DxC_comparable_scale,uric_acid_dxc,2015_2016_DxC_660i,identity,source_value,mg/dL,mg/dL,I,https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BIOPRO_J.htm,0,,flag_nonpositive_or_nonfinite;do_not_silently_drop,final,2026-08-14.rev22-confirmed-plan

serum_creatinine,P,LBXSCR,Cobas_6000,creatinine_dxc,2015_2016_DxC_660i,Cobas_6000_to_DxC_660i,1.051 * LBXSCR - 0.06945,mg/dL,mg/dL,P,https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BIOPRO_J.htm,0,,flag_nonpositive_or_nonfinite;do_not_silently_drop,final,2026-08-14.rev22-confirmed-plan

serum_uric_acid,P,LBXSUA,Cobas_6000,uric_acid_dxc,2015_2016_DxC_660i,Cobas_6000_to_DxC_660i,0.9323 * LBXSUA + 0.2326,mg/dL,mg/dL,P,https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/BIOPRO_J.htm,0,,flag_nonpositive_or_nonfinite;do_not_silently_drop,final,2026-08-14.rev22-confirmed-plan

)--------------------"

embedded_config[["NHANES_config/phase1_variable_extensions.csv"]] <- r"--------------------(cycle,module,stem,source_variables,purpose,configuration_status

E,Examination,bmx,BMXBMI,continuous BMI for C2/C4,approved_schema_verified_20260824

E,Questionnaire,bpq,BPQ080;BPQ090D;BPQ100D,high cholesterol composite with current medication,approved_schema_verified_20260824

E,Questionnaire,paq,PAQ610;PAD615;PAQ625;PAD630;PAQ640;PAD645;PAQ655;PAD660;PAQ670;PAD675,GPAQ MET-min/week,approved_schema_verified_20260824

F,Examination,bmx,BMXBMI,continuous BMI for C2/C4,approved_schema_verified_20260824

F,Questionnaire,bpq,BPQ080;BPQ090D;BPQ100D,high cholesterol composite with current medication,approved_schema_verified_20260824

F,Questionnaire,paq,PAQ610;PAD615;PAQ625;PAD630;PAQ640;PAD645;PAQ655;PAD660;PAQ670;PAD675,GPAQ MET-min/week,approved_schema_verified_20260824

G,Examination,bmx,BMXBMI,continuous BMI for C2/C4,approved_schema_verified_20260824

G,Questionnaire,bpq,BPQ080;BPQ090D;BPQ100D,high cholesterol composite with current medication,approved_schema_verified_20260824

G,Questionnaire,paq,PAQ610;PAD615;PAQ625;PAD630;PAQ640;PAD645;PAQ655;PAD660;PAQ670;PAD675,GPAQ MET-min/week,approved_schema_verified_20260824

H,Examination,bmx,BMXBMI,continuous BMI for C2/C4,approved_schema_verified_20260824

H,Questionnaire,bpq,BPQ080;BPQ090D;BPQ100D,high cholesterol composite with current medication,approved_schema_verified_20260824

H,Questionnaire,paq,PAQ610;PAD615;PAQ625;PAD630;PAQ640;PAD645;PAQ655;PAD660;PAQ670;PAD675,GPAQ MET-min/week,approved_schema_verified_20260824

I,Examination,bmx,BMXBMI,continuous BMI for C2/C4,approved_schema_verified_20260824

I,Questionnaire,bpq,BPQ080;BPQ090D;BPQ100D,high cholesterol composite with current medication,approved_schema_verified_20260824

I,Questionnaire,paq,PAQ610;PAD615;PAQ625;PAD630;PAQ640;PAD645;PAQ655;PAD660;PAQ670;PAD675,GPAQ MET-min/week,approved_schema_verified_20260824

P,Examination,bmx,BMXBMI,continuous BMI for C2/C4,approved_schema_verified_20260824

P,Questionnaire,bpq,BPQ080;BPQ090D;BPQ100D,high cholesterol composite with current medication,approved_schema_verified_20260824

P,Questionnaire,paq,PAQ610;PAD615;PAQ625;PAD630;PAQ640;PAD645;PAQ655;PAD660;PAQ670;PAD675,GPAQ MET-min/week,approved_schema_verified_20260824

)--------------------"

embedded_config[["NHANES_config/plain_water_codes.csv"]] <- r"--------------------(cycle,food_code,description,nova_in_manifest,plain_water_rule,basis,configuration_status,version

E,92410210,"Carbonated water, unsweetened",4,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,unsweetened carbonated water; no flavor stated,final,2026-08-14.rev22-confirmed-plan

E,94000100,"Water, tap",1,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,tap water,final,2026-08-14.rev22-confirmed-plan

E,94100100,"Water, bottled, unsweetened",1,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,plain bottled water by unsweetened unflavored description,final,2026-08-14.rev22-confirmed-plan

F,92410210,"Carbonated water, unsweetened",4,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,unsweetened carbonated water; no flavor stated,final,2026-08-14.rev22-confirmed-plan

F,94000100,"Water, tap",1,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,tap water,final,2026-08-14.rev22-confirmed-plan

F,94100100,"Water, bottled, unsweetened",1,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,plain bottled water by unsweetened unflavored description,final,2026-08-14.rev22-confirmed-plan

F,94300100,"Water, baby, bottled, unsweetened",1,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,bottled unsweetened water; no flavor stated,final,2026-08-14.rev22-confirmed-plan

G,92410210,"Carbonated water, unsweetened",4,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,unsweetened carbonated water; no flavor stated,final,2026-08-14.rev22-confirmed-plan

G,94000100,"Water, tap",1,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,tap water,final,2026-08-14.rev22-confirmed-plan

G,94100100,"Water, bottled, unsweetened",1,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,plain bottled water by unsweetened unflavored description,final,2026-08-14.rev22-confirmed-plan

G,94300100,"Water, baby, bottled, unsweetened",1,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,bottled unsweetened water; no flavor stated,final,2026-08-14.rev22-confirmed-plan

H,92410210,"Carbonated water, unsweetened",4,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,unsweetened carbonated water; no flavor stated,final,2026-08-14.rev22-confirmed-plan

H,94000100,"Water, tap",1,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,tap water,final,2026-08-14.rev22-confirmed-plan

H,94100100,"Water, bottled, unsweetened",1,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,plain bottled water by unsweetened unflavored description,final,2026-08-14.rev22-confirmed-plan

H,94300100,"Water, baby, bottled, unsweetened",1,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,bottled unsweetened water; no flavor stated,final,2026-08-14.rev22-confirmed-plan

I,92410210,"Carbonated water, unsweetened",4,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,unsweetened carbonated water; no flavor stated,final,2026-08-14.rev22-confirmed-plan

I,94000100,"Water, tap",1,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,tap water,final,2026-08-14.rev22-confirmed-plan

I,94100100,"Water, bottled, unsweetened",1,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,plain bottled water by unsweetened unflavored description,final,2026-08-14.rev22-confirmed-plan

I,94300100,"Water, baby, bottled, unsweetened",1,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,bottled unsweetened water; no flavor stated,final,2026-08-14.rev22-confirmed-plan

P,92410210,"Water, carbonated, plain",4,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,explicit plain carbonated water,final,2026-08-14.rev22-confirmed-plan

P,94000010,"Water, NFS",1,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,unqualified drinking water code,final,2026-08-14.rev22-confirmed-plan

P,94000100,"Water, tap",1,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,tap water,final,2026-08-14.rev22-confirmed-plan

P,94100100,"Water, bottled, plain",1,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,explicit plain bottled water,final,2026-08-14.rev22-confirmed-plan

P,94300100,"Water, baby",1,tap_or_plain_bottled_or_unsweetened_unflavored_carbonated,water code with no flavor or sweetener stated,final,2026-08-14.rev22-confirmed-plan

)--------------------"

embedded_config[["NHANES_config/ssb_codes.csv"]] <- r"--------------------(cycle,food_code,description,nova_in_manifest,ssb_family,basis,configuration_status,version

E,92121000,"Coffee, made from powdered instant mix, with whitener and sugar, instant",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92130020,"Coffee, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92301060,"Tea, NS as to type, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92301160,"Tea, NS as to type, decaffeinated, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92302200,"Tea, leaf, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92302600,"Tea, leaf, decaffeinated, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92305040,"Tea, made from powdered instant, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92305050,"Tea, made from powdered instant, decaffeinated, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92306020,"Tea, herbal, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92400000,"Soft drink, NFS",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92410110,"Carbonated water, sweetened",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92410310,"Soft drink, cola-type",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92410315,"Soft drink, cola type, reduced sugar",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92410330,"Soft drink, cola-type, with higher caffeine",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92410340,"Soft drink, cola-type, decaffeinated",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92410360,"Soft drink, pepper-type",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92410390,"Soft drink, pepper-type, decaffeinated",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92410410,Cream soda,4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92410510,"Soft drink, fruit-flavored, caffeine free",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92410550,"Soft drink, fruit flavored, caffeine containing",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92410810,Chocolate-flavored soda,4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92411510,Cola with fruit or vanilla flavor,4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92411520,Cola with chocolate flavor,4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92416010,Mavi drink,4,other_explicit_sugar_beverage,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92417010,"Soft drink, ale type",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92431000,"Carbonated juice drink, NS as to type of juice",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92432000,Carbonated citrus juice drink,4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92433000,Carbonated noncitrus juice drink,4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92510610,Fruit juice drink,4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92510720,"Fruit punch, made with fruit juice and soda",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92510730,"Fruit punch, made with soda, fruit juice, and sherbet or ice cream",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92511010,Fruit flavored drink (formerly lemonade),4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92511250,"Citrus fruit juice drink, containing 40-50% juice",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92530410,"Fruit flavored drink, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92530510,"Cranberry juice drink or cocktail, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92530610,"Fruit juice drink, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92530950,"Vegetable and fruit juice drink, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92531030,"Fruit juice drink, with thiamin (vitamin B1) and high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92541010,"Fruit flavored drink, made from powdered mix",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92542000,"Fruit flavored drink, made from powdered mix,with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92552020,"Fruit juice drink, reduced sugar, with thiamin (vitamin B1) and high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92552030,"Fruit juice drink, reduced sugar, with vitamin E",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92560000,Fruit-flavored thirst quencher beverage,4,sports_thirst_quencher,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92560100,Gatorade Thirst Quencher sports drink,4,sports_thirst_quencher,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92560200,Powerade sports drink,4,sports_thirst_quencher,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92582100,"Fruit juice drink, with high vitamin C, plus added calcium",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92582110,"Fruit juice drink, with thiamin (vitamin B1) and high vitamin C plus calcium",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92582120,"Fruit flavored drink, reduced sugar, with high vitamin C, plus added calcium",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92650000,Red Bull Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92650100,Full Throttle Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92650200,Monster Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92650205,Mountain Dew AMP Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92650700,Rockstar Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92650800,Vault Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92651000,Energy drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

E,92804000,Shirley Temple,4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92121000,"Coffee, made from powdered instant mix, with whitener and sugar, instant",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92130020,"Coffee, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92301060,"Tea, NS as to type, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92301160,"Tea, NS as to type, decaffeinated, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92302200,"Tea, leaf, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92302600,"Tea, leaf, decaffeinated, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92305040,"Tea, made from powdered instant, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92305050,"Tea, made from powdered instant, decaffeinated, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92306020,"Tea, herbal, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92400000,"Soft drink, NFS",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92410110,"Carbonated water, sweetened",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92410310,"Soft drink, cola-type",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92410315,"Soft drink, cola type, reduced sugar",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92410330,"Soft drink, cola-type, with higher caffeine",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92410340,"Soft drink, cola-type, decaffeinated",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92410360,"Soft drink, pepper-type",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92410390,"Soft drink, pepper-type, decaffeinated",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92410410,Cream soda,4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92410510,"Soft drink, fruit-flavored, caffeine free",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92410550,"Soft drink, fruit flavored, caffeine containing",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92410810,Chocolate-flavored soda,4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92411510,Cola with fruit or vanilla flavor,4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92411520,Cola with chocolate flavor,4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92416010,Mavi drink,4,other_explicit_sugar_beverage,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92417010,"Soft drink, ale type",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92431000,"Carbonated juice drink, NS as to type of juice",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92432000,Carbonated citrus juice drink,4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92433000,Carbonated noncitrus juice drink,4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92510610,Fruit juice drink,4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92510720,"Fruit punch, made with fruit juice and soda",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92510730,"Fruit punch, made with soda, fruit juice, and sherbet or ice cream",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92511010,Fruit flavored drink (formerly lemonade),4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92511250,"Citrus fruit juice drink, containing 40-50% juice",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92530410,"Fruit flavored drink, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92530510,"Cranberry juice drink or cocktail, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92530610,"Fruit juice drink, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92530950,"Vegetable and fruit juice drink, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92531030,"Fruit juice drink, with thiamin (vitamin B1) and high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92541010,"Fruit flavored drink, made from powdered mix",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92542000,"Fruit flavored drink, made from powdered mix,with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92552020,"Fruit juice drink, reduced sugar, with thiamin (vitamin B1) and high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92552030,"Fruit juice drink, reduced sugar, with vitamin E",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92560000,Fruit-flavored thirst quencher beverage,4,sports_thirst_quencher,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92560100,Gatorade Thirst Quencher sports drink,4,sports_thirst_quencher,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92560200,Powerade sports drink,4,sports_thirst_quencher,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92582100,"Fruit juice drink, with high vitamin C, plus added calcium",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92582110,"Fruit juice drink, with thiamin (vitamin B1) and high vitamin C plus calcium",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92582120,"Fruit flavored drink, reduced sugar, with high vitamin C, plus added calcium",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92650000,Red Bull Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92650100,Full Throttle Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92650200,Monster Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92650205,Mountain Dew AMP Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92650700,Rockstar Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92650800,Vault Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92651000,Energy drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

F,92804000,Shirley Temple,4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92121000,"Coffee, made from powdered instant mix, with whitener and sugar, instant",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92130020,"Coffee, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92301060,"Tea, NS as to type, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92301160,"Tea, NS as to type, decaffeinated, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92302200,"Tea, leaf, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92302600,"Tea, leaf, decaffeinated, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92305040,"Tea, made from powdered instant, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92305050,"Tea, made from powdered instant, decaffeinated, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92306020,"Tea, herbal, presweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92307500,"Half and Half beverage, half iced tea and half fruit juice drink (lemonade)",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92400000,"Soft drink, NFS",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92410110,"Carbonated water, sweetened",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92410310,"Soft drink, cola-type",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92410315,"Soft drink, cola type, reduced sugar",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92410330,"Soft drink, cola-type, with higher caffeine",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92410340,"Soft drink, cola-type, decaffeinated",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92410360,"Soft drink, pepper-type",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92410390,"Soft drink, pepper-type, decaffeinated",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92410410,Cream soda,4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92410510,"Soft drink, fruit-flavored, caffeine free",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92410550,"Soft drink, fruit flavored, caffeine containing",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92410810,Chocolate-flavored soda,4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92411510,Cola with fruit or vanilla flavor,4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92411520,Cola with chocolate flavor,4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92417010,"Soft drink, ale type",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92431000,"Carbonated juice drink, NS as to type of juice",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92432000,Carbonated citrus juice drink,4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92433000,Carbonated noncitrus juice drink,4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92510610,Fruit juice drink,4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92510720,"Fruit punch, made with fruit juice and soda",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92510730,"Fruit punch, made with soda, fruit juice, and sherbet or ice cream",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92511010,Fruit flavored drink (formerly lemonade),4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92511250,"Citrus fruit juice drink, containing 40-50% juice",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92530410,"Fruit flavored drink, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92530510,"Cranberry juice drink or cocktail, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92530610,"Fruit juice drink, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92530950,"Vegetable and fruit juice drink, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92531030,"Fruit juice drink, with thiamin (vitamin B1) and high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92541010,"Fruit flavored drink, made from powdered mix",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92542000,"Fruit flavored drink, made from powdered mix,with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92552020,"Fruit juice drink, reduced sugar, with thiamin (vitamin B1) and high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92552030,"Fruit juice drink, reduced sugar, with vitamin E",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92582100,"Fruit juice drink, with high vitamin C, plus added calcium",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92582110,"Fruit juice drink, with thiamin (vitamin B1) and high vitamin C plus calcium",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,92804000,Shirley Temple,4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,95310200,Full Throttle Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,95310400,Monster Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,95310500,Mountain Dew AMP Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,95310550,No Fear Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,95310555,No Fear Motherload Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,95310560,NOS Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,95310600,Red Bull Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,95310700,Rockstar Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,95310800,Vault Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,95311000,Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,95320200,Gatorade Thirst Quencher sports drink,4,sports_thirst_quencher,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,95320500,Powerade sports drink,4,sports_thirst_quencher,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

G,95321000,Fruit-flavored thirst quencher beverage,4,sports_thirst_quencher,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92121010,"Coffee, instant, pre-sweetened with sugar, reconstituted",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92130020,"Coffee, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92130021,"Coffee, decaffeinated, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92305040,"Tea, iced, instant, black, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92305050,"Tea, iced, instant, black, decaffeinated, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92305910,"Tea, iced, instant, green, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92307500,Iced Tea / Lemonade juice drink,4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92308000,"Tea, iced, brewed, black, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92308030,"Tea, iced, brewed, black, decaffeinated, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92308500,"Tea, iced, brewed, green, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92308530,"Tea, iced, brewed, green, decaffeinated, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92400000,"Soft drink, NFS",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92410110,"Carbonated water, sweetened",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92410310,"Soft drink, cola",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92410315,"Soft drink, cola, reduced sugar",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92410340,"Soft drink, cola, decaffeinated",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92410360,"Soft drink, pepper type",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92410390,"Soft drink, pepper type, decaffeinated",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92410410,"Soft drink, cream soda",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92410510,"Soft drink, fruit flavored, caffeine free",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92410550,"Soft drink, fruit flavored, caffeine containing",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92410610,"Soft drink, ginger ale",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92410710,"Soft drink, root beer",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92410810,"Soft drink, chocolate flavored",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92411510,"Soft drink, cola, fruit or vanilla flavored",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92411520,"Soft drink, cola, chocolate flavored",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92432000,"Fruit juice drink, citrus, carbonated",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92433000,"Fruit juice drink, noncitrus, carbonated",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92510610,Fruit juice drink,4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92510720,"Fruit punch, made with fruit juice and soda",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92510730,"Fruit punch, made with soda, fruit juice, and sherbet or ice cream",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92510955,"Lemonade, fruit juice drink",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92510960,"Lemonade, fruit flavored drink",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92511015,Fruit flavored drink,4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92511250,"Fruit juice beverage, 40-50% juice, citrus",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92530410,"Fruit flavored drink, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92530510,"Cranberry juice drink, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92530610,"Fruit juice drink, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92530950,"Vegetable and fruit juice drink, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92541010,"Fruit flavored drink, powdered, reconstituted",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92542000,"Fruit flavored drink, with high vitamin C, powdered, reconstituted",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92552030,"Capri Sun, fruit juice drink",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92582100,"Fruit juice drink, with high vitamin C, plus added calcium",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,92804000,Shirley Temple,4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,95310200,Full Throttle Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,95310400,Monster Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,95310500,Mountain Dew AMP Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,95310550,No Fear Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,95310555,No Fear Motherload Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,95310560,NOS Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,95310600,Red Bull Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,95310700,Rockstar Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,95310800,Vault Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,95311000,Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,95320200,Gatorade G sports drink,4,sports_thirst_quencher,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,95320500,Powerade sports drink,4,sports_thirst_quencher,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

H,95321000,"Sports drink, NFS",4,sports_thirst_quencher,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92121010,"Coffee, instant, pre-sweetened with sugar, reconstituted",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92130020,"Coffee, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92130021,"Coffee, decaffeinated, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92305040,"Tea, iced, instant, black, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92305050,"Tea, iced, instant, black, decaffeinated, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92305910,"Tea, iced, instant, green, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92307500,Iced Tea / Lemonade juice drink,4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92308000,"Tea, iced, brewed, black, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92308030,"Tea, iced, brewed, black, decaffeinated, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92308500,"Tea, iced, brewed, green, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92308530,"Tea, iced, brewed, green, decaffeinated, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92400000,"Soft drink, NFS",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92410110,"Carbonated water, sweetened",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92410310,"Soft drink, cola",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92410315,"Soft drink, cola, reduced sugar",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92410340,"Soft drink, cola, decaffeinated",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92410360,"Soft drink, pepper type",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92410390,"Soft drink, pepper type, decaffeinated",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92410410,"Soft drink, cream soda",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92410510,"Soft drink, fruit flavored, caffeine free",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92410550,"Soft drink, fruit flavored, caffeine containing",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92410610,"Soft drink, ginger ale",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92410710,"Soft drink, root beer",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92410810,"Soft drink, chocolate flavored",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92411510,"Soft drink, cola, fruit or vanilla flavored",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92411520,"Soft drink, cola, chocolate flavored",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92432000,"Fruit juice drink, citrus, carbonated",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92433000,"Fruit juice drink, noncitrus, carbonated",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92510610,Fruit juice drink,4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92510720,"Fruit punch, made with fruit juice and soda",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92510730,"Fruit punch, made with soda, fruit juice, and sherbet or ice cream",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92510955,"Lemonade, fruit juice drink",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92510960,"Lemonade, fruit flavored drink",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92511015,Fruit flavored drink,4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92511250,"Fruit juice beverage, 40-50% juice, citrus",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92530410,"Fruit flavored drink, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92530510,"Cranberry juice drink, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92530610,"Fruit juice drink, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92530950,"Vegetable and fruit juice drink, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92531030,Fruit juice drink (Sunny D),4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92541010,"Fruit flavored drink, powdered, reconstituted",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92542000,"Fruit flavored drink, with high vitamin C, powdered, reconstituted",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92552020,"Fruit juice drink, reduced sugar (Sunny D)",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92552030,Fruit juice drink (Capri Sun),4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92582100,"Fruit juice drink, with high vitamin C, plus added calcium",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92582110,"Fruit juice drink, added calcium (Sunny D)",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,92804000,Shirley Temple,4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,95310200,Energy drink (Full Throttle),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,95310400,Energy drink (Monster),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,95310500,Energy drink (Mountain Dew AMP),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,95310550,Energy drink (No Fear),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,95310555,Energy drink (No Fear Motherload),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,95310560,Energy drink (NOS),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,95310600,Energy drink (Red Bull),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,95310700,Energy drink (Rockstar),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,95310750,Energy drink (SoBe Energize Energy Juice Drink),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,95310800,Energy drink (Vault),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,95311000,Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,95312560,Energy drink (Ocean Spray Cran-Energy Juice Drink),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,95320200,Sports drink (Gatorade G),4,sports_thirst_quencher,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,95320500,Sports drink (Powerade),4,sports_thirst_quencher,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

I,95321000,"Sports drink, NFS",4,sports_thirst_quencher,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92121010,"Coffee, instant, pre-sweetened with sugar, reconstituted",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92130020,"Coffee, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92130021,"Coffee, decaffeinated, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92305040,"Tea, iced, instant, black, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92305050,"Tea, iced, instant, black, decaffeinated, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92305910,"Tea, iced, instant, green, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92307500,Iced Tea / Lemonade juice drink,4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92308000,"Tea, iced, brewed, black, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92308030,"Tea, iced, brewed, black, decaffeinated, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92308500,"Tea, iced, brewed, green, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92308530,"Tea, iced, brewed, green, decaffeinated, pre-sweetened with sugar",4,presweetened_coffee_tea,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92400000,"Soft drink, NFS",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92410310,"Soft drink, cola",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92410315,"Soft drink, cola, reduced sugar",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92410340,"Soft drink, cola, decaffeinated",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92410360,"Soft drink, pepper type",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92410390,"Soft drink, pepper type, decaffeinated",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92410410,"Soft drink, cream soda",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92410510,"Soft drink, fruit flavored, caffeine free",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92410550,"Soft drink, fruit flavored, caffeine containing",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92410610,"Soft drink, ginger ale",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92410710,"Soft drink, root beer",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92410810,"Soft drink, chocolate flavored",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92411510,"Soft drink, cola, fruit or vanilla flavored",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92411520,"Soft drink, cola, chocolate flavored",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92432000,"Fruit juice drink, citrus, carbonated",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92433000,"Fruit juice drink, noncitrus, carbonated",4,regular_carbonated_soft_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92510610,Fruit juice drink,4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92510720,"Fruit punch, made with fruit juice and soda",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92510730,"Fruit punch, made with soda, fruit juice, and sherbet or ice cream",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92510955,"Lemonade, fruit juice drink",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92510960,"Lemonade, fruit flavored drink",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92511015,Fruit flavored drink,4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92511250,"Fruit juice beverage, 40-50% juice, citrus",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92530410,"Fruit flavored drink, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92530510,"Cranberry juice drink, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92530610,"Fruit juice drink, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92530950,"Vegetable and fruit juice drink, with high vitamin C",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92531030,Fruit juice drink (Sunny D),4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92541010,"Fruit flavored drink, powdered, reconstituted",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92542000,"Fruit flavored drink, with high vitamin C, powdered, reconstituted",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92552020,"Fruit juice drink, reduced sugar (Sunny D)",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92552030,Fruit juice drink (Capri Sun),4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92582100,"Fruit juice drink, with high vitamin C, plus added calcium",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92582110,"Fruit juice drink, added calcium (Sunny D)",4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,92804000,Shirley Temple,4,fruit_drink_lemonade,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,95310200,Energy drink (Full Throttle),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,95310400,Energy drink (Monster),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,95310500,Energy drink (Mountain Dew AMP),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,95310550,Energy drink (No Fear),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,95310555,Energy drink (No Fear Motherload),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,95310560,Energy drink (NOS),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,95310600,Energy drink (Red Bull),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,95310700,Energy drink (Rockstar),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,95310750,Energy drink (SoBe Energize Energy Juice Drink),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,95310800,Energy drink (Vault),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,95311000,Energy Drink,4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,95312560,Energy drink (Ocean Spray Cran-Energy Juice Drink),4,energy_drink,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,95320200,Sports drink (Gatorade G),4,sports_thirst_quencher,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,95320500,Sports drink (Powerade),4,sports_thirst_quencher,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

P,95321000,"Sports drink, NFS",4,sports_thirst_quencher,NOVA4 plus regular sugar-bearing beverage description; diet/low-calorie/light/zero/sugar-free/unsweetened/nonreconstituted/alcohol/milk excluded,final,2026-08-14.rev22-confirmed-plan

)--------------------"

embedded_config[["NHANES_config/variable_map.csv"]] <- r"--------------------(cycle,cycle_label,module,stem,concept,source_variables,target_variable,unit_or_codes,required_for_rev22,availability,transform_rule,basis,configuration_status,version

E,2007-2008,Demographics,demo,participant_id,SEQN,SEQN,identifier,true,present,none,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Demographics,demo,sex,RIAGENDR,sex,codes 1 male 2 female,true,present,factor labels,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Demographics,demo,age,RIDAGEYR,age,years,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Demographics,demo,race_ethnicity,RIDRETH1,race_ethnicity,NHANES five-category codes,true,present,factor 1:5,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Demographics,demo,education,DMDEDUC2,education,NHANES adult education codes,true,present,1-2 below HS;3 HS/GED;4-5 above HS,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Demographics,demo,pir,INDFMPIR,pir,ratio,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Demographics,demo,pregnancy,RIDEXPRG,ridexprg,NHANES pregnancy status codes,true,present,derive nonpregnant domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Demographics,demo,survey_psu,SDMVPSU,SDMVPSU,masked variance unit,true,present,none,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Demographics,demo,survey_strata,SDMVSTRA,SDMVSTRA,masked variance stratum,true,present,none,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr1tot,day1_recall_status,DR1DRSTZ,DR1DRSTZ,1 reliable,true,present,require 1 in domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr2tot,day2_recall_status,DR2DRSTZ,DR2DRSTZ,1 reliable,true,present,require 1 in domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr2tot,number_of_recalls,DRDINT,DRDINT,2 means two recalls,true,present,require 2 in domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr1tot,day1_energy,DR1TKCAL,day1_energy_kcal,kcal/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr2tot,day2_energy,DR2TKCAL,day2_energy_kcal,kcal/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr1tot,day1_sodium,DR1TSODI,day1_sodium_mg,mg/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr2tot,day2_sodium,DR2TSODI,day2_sodium_mg,mg/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr1tot,day1_saturated_fat,DR1TSFAT,day1_satfat_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr2tot,day2_saturated_fat,DR2TSFAT,day2_satfat_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr1tot,day1_monounsaturated_fat,DR1TMFAT,day1_mufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr2tot,day2_monounsaturated_fat,DR2TMFAT,day2_mufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr1tot,day1_polyunsaturated_fat,DR1TPFAT,day1_pufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr2tot,day2_polyunsaturated_fat,DR2TPFAT,day2_pufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr1tot,day1_alcohol_nutrient,DR1TALCO,day1_alcohol_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr2tot,day2_alcohol_nutrient,DR2TALCO,day2_alcohol_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Questionnaire,bpq,hypertension_diagnosed,BPQ020,bpq020,NHANES yes/no codes,true,present,yes_no_na,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Questionnaire,bpq,hypertension_medication,BPQ050A,bpq050a,NHANES yes/no codes,true,present,yes_no_na,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Questionnaire,diq,diabetes_diagnosed,DIQ010,diq010,NHANES yes/no/borderline codes,true,present,derive diabetes,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Questionnaire,diq,insulin_current,DIQ050,diq050,NHANES yes/no codes,true,present,derive diabetes,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Questionnaire,smq,smoked_100_cigarettes,SMQ020,smq020,NHANES yes/no codes,true,present,derive smoking,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Questionnaire,smq,current_smoking,SMQ040,smq040,every day/some days/not at all,true,present,derive smoking,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Questionnaire,kiq_u,kidney_stone_history,KIQ026,kiq026,NHANES yes/no codes,true,present,yes_no_na,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Laboratory,ghb,hba1c,LBXGH,hba1c,percent,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Laboratory,tchol,total_cholesterol,LBXTC,total_cholesterol,mg/dL,true,present,numeric; fixed LBXTC,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Laboratory,cbc,platelet_count,LBXPLTSI,platelet,10^3 cells/uL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Laboratory,cbc,absolute_neutrophils,LBDNENO,neutrophil,10^3 cells/uL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Laboratory,cbc,absolute_lymphocytes,LBDLYMNO,lymphocyte,10^3 cells/uL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Laboratory,biopro,serum_creatinine,LBXSCR,creatinine_raw,mg/dL,true,present,laboratory_bridge,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Laboratory,biopro,serum_uric_acid,LBXSUA,uric_acid_raw,mg/dL,true,present,laboratory_bridge,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Laboratory,alb_cr,urine_albumin,URXUMA,urine_albumin_ug_ml,ug/mL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Laboratory,alb_cr,urine_creatinine,URXUCR,urine_creatinine_mg_dl,mg/dL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr1iff,food_code_day1,DR1IFDCD,food_code_day1,FNDDS code,true,present,character key to cycle manifest,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr1iff,food_grams_day1,DR1IGRMS,food_grams_day1,g,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr1iff,food_energy_day1,DR1IKCAL,food_energy_day1,kcal,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr1iff,food_source_day1,DR1FS,food_source_day1,NHANES source code,true,present,retain for adjudication,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr1iff,combination_type_day1,DR1CCMTX,combination_type_day1,NHANES combination type,true,present,retain for adjudication,08_测试/schema_inventory.csv and CDC IFF codebook,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr2iff,food_code_day2,DR2IFDCD,food_code_day2,FNDDS code,true,present,character key to cycle manifest,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr2iff,food_grams_day2,DR2IGRMS,food_grams_day2,g,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr2iff,food_energy_day2,DR2IKCAL,food_energy_day2,kcal,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr2iff,food_source_day2,DR2FS,food_source_day2,NHANES source code,true,present,retain for adjudication,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr2iff,combination_type_day2,DR2CCMTX,combination_type_day2,NHANES combination type,true,present,retain for adjudication,08_测试/schema_inventory.csv and CDC IFF codebook,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Dietary,dr2tot,diet_weight,WTDR2D,analysis_weight,regular two-day dietary weight,true,present,divide by 13.2 pooled years,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Questionnaire,diq,oral_hypoglycemic_medication,DID070,diq070,NHANES yes/no codes,true,present,derive diabetes,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Questionnaire,alq,alcohol_frequency,ALQ101;ALQ110;ALQ120Q;ALQ120U;ALQ130,alcohol,cycle-specific questionnaire,true,present,derive weekly categories using frozen structural-skip map,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Questionnaire,paq,physical_activity,PAQ605;PAQ620;PAQ635;PAQ650;PAQ665,physical_activity,five NHANES yes/no items,true,present,any valid yes active; all valid no inactive; otherwise missing,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Examination,bpx,blood_pressure,BPXSY1;BPXSY2;BPXSY3;BPXSY4;BPXDI1;BPXDI2;BPXDI3;BPXDI4,mean_sbp;mean_dbp,mmHg,true,present,mean of available valid readings,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

E,2007-2008,Laboratory,alb_cr,albumin_creatinine_ratio,URXUMA;URXUCR,acr_mg_g,mg/g,true,present,derive 100 * URXUMA / URXUCR; first random urine,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Demographics,demo,participant_id,SEQN,SEQN,identifier,true,present,none,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Demographics,demo,sex,RIAGENDR,sex,codes 1 male 2 female,true,present,factor labels,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Demographics,demo,age,RIDAGEYR,age,years,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Demographics,demo,race_ethnicity,RIDRETH1,race_ethnicity,NHANES five-category codes,true,present,factor 1:5,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Demographics,demo,education,DMDEDUC2,education,NHANES adult education codes,true,present,1-2 below HS;3 HS/GED;4-5 above HS,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Demographics,demo,pir,INDFMPIR,pir,ratio,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Demographics,demo,pregnancy,RIDEXPRG,ridexprg,NHANES pregnancy status codes,true,present,derive nonpregnant domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Demographics,demo,survey_psu,SDMVPSU,SDMVPSU,masked variance unit,true,present,none,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Demographics,demo,survey_strata,SDMVSTRA,SDMVSTRA,masked variance stratum,true,present,none,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr1tot,day1_recall_status,DR1DRSTZ,DR1DRSTZ,1 reliable,true,present,require 1 in domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr2tot,day2_recall_status,DR2DRSTZ,DR2DRSTZ,1 reliable,true,present,require 1 in domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr2tot,number_of_recalls,DRDINT,DRDINT,2 means two recalls,true,present,require 2 in domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr1tot,day1_energy,DR1TKCAL,day1_energy_kcal,kcal/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr2tot,day2_energy,DR2TKCAL,day2_energy_kcal,kcal/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr1tot,day1_sodium,DR1TSODI,day1_sodium_mg,mg/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr2tot,day2_sodium,DR2TSODI,day2_sodium_mg,mg/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr1tot,day1_saturated_fat,DR1TSFAT,day1_satfat_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr2tot,day2_saturated_fat,DR2TSFAT,day2_satfat_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr1tot,day1_monounsaturated_fat,DR1TMFAT,day1_mufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr2tot,day2_monounsaturated_fat,DR2TMFAT,day2_mufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr1tot,day1_polyunsaturated_fat,DR1TPFAT,day1_pufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr2tot,day2_polyunsaturated_fat,DR2TPFAT,day2_pufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr1tot,day1_alcohol_nutrient,DR1TALCO,day1_alcohol_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr2tot,day2_alcohol_nutrient,DR2TALCO,day2_alcohol_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Questionnaire,bpq,hypertension_diagnosed,BPQ020,bpq020,NHANES yes/no codes,true,present,yes_no_na,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Questionnaire,bpq,hypertension_medication,BPQ050A,bpq050a,NHANES yes/no codes,true,present,yes_no_na,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Questionnaire,diq,diabetes_diagnosed,DIQ010,diq010,NHANES yes/no/borderline codes,true,present,derive diabetes,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Questionnaire,diq,insulin_current,DIQ050,diq050,NHANES yes/no codes,true,present,derive diabetes,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Questionnaire,smq,smoked_100_cigarettes,SMQ020,smq020,NHANES yes/no codes,true,present,derive smoking,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Questionnaire,smq,current_smoking,SMQ040,smq040,every day/some days/not at all,true,present,derive smoking,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Questionnaire,kiq_u,kidney_stone_history,KIQ026,kiq026,NHANES yes/no codes,true,present,yes_no_na,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Laboratory,ghb,hba1c,LBXGH,hba1c,percent,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Laboratory,tchol,total_cholesterol,LBXTC,total_cholesterol,mg/dL,true,present,numeric; fixed LBXTC,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Laboratory,cbc,platelet_count,LBXPLTSI,platelet,10^3 cells/uL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Laboratory,cbc,absolute_neutrophils,LBDNENO,neutrophil,10^3 cells/uL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Laboratory,cbc,absolute_lymphocytes,LBDLYMNO,lymphocyte,10^3 cells/uL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Laboratory,biopro,serum_creatinine,LBXSCR,creatinine_raw,mg/dL,true,present,laboratory_bridge,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Laboratory,biopro,serum_uric_acid,LBXSUA,uric_acid_raw,mg/dL,true,present,laboratory_bridge,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Laboratory,alb_cr,urine_albumin,URXUMA,urine_albumin_ug_ml,ug/mL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Laboratory,alb_cr,urine_creatinine,URXUCR,urine_creatinine_mg_dl,mg/dL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr1iff,food_code_day1,DR1IFDCD,food_code_day1,FNDDS code,true,present,character key to cycle manifest,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr1iff,food_grams_day1,DR1IGRMS,food_grams_day1,g,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr1iff,food_energy_day1,DR1IKCAL,food_energy_day1,kcal,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr1iff,food_source_day1,DR1FS,food_source_day1,NHANES source code,true,present,retain for adjudication,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr1iff,combination_type_day1,DR1CCMTX,combination_type_day1,NHANES combination type,true,present,retain for adjudication,08_测试/schema_inventory.csv and CDC IFF codebook,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr2iff,food_code_day2,DR2IFDCD,food_code_day2,FNDDS code,true,present,character key to cycle manifest,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr2iff,food_grams_day2,DR2IGRMS,food_grams_day2,g,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr2iff,food_energy_day2,DR2IKCAL,food_energy_day2,kcal,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr2iff,food_source_day2,DR2FS,food_source_day2,NHANES source code,true,present,retain for adjudication,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr2iff,combination_type_day2,DR2CCMTX,combination_type_day2,NHANES combination type,true,present,retain for adjudication,08_测试/schema_inventory.csv and CDC IFF codebook,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Dietary,dr2tot,diet_weight,WTDR2D,analysis_weight,regular two-day dietary weight,true,present,divide by 13.2 pooled years,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Questionnaire,diq,oral_hypoglycemic_medication,DIQ070,diq070,NHANES yes/no codes,true,present,derive diabetes,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Questionnaire,alq,alcohol_frequency,ALQ101;ALQ110;ALQ120Q;ALQ120U;ALQ130,alcohol,cycle-specific questionnaire,true,present,derive weekly categories using frozen structural-skip map,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Questionnaire,paq,physical_activity,PAQ605;PAQ620;PAQ635;PAQ650;PAQ665,physical_activity,five NHANES yes/no items,true,present,any valid yes active; all valid no inactive; otherwise missing,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Examination,bpx,blood_pressure,BPXSY1;BPXSY2;BPXSY3;BPXSY4;BPXDI1;BPXDI2;BPXDI3;BPXDI4,mean_sbp;mean_dbp,mmHg,true,present,mean of available valid readings,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

F,2009-2010,Laboratory,alb_cr,albumin_creatinine_ratio,URDACT,acr_mg_g,mg/g,true,present,use released URDACT; first random urine,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Demographics,demo,participant_id,SEQN,SEQN,identifier,true,present,none,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Demographics,demo,sex,RIAGENDR,sex,codes 1 male 2 female,true,present,factor labels,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Demographics,demo,age,RIDAGEYR,age,years,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Demographics,demo,race_ethnicity,RIDRETH1,race_ethnicity,NHANES five-category codes,true,present,factor 1:5,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Demographics,demo,education,DMDEDUC2,education,NHANES adult education codes,true,present,1-2 below HS;3 HS/GED;4-5 above HS,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Demographics,demo,pir,INDFMPIR,pir,ratio,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Demographics,demo,pregnancy,RIDEXPRG,ridexprg,NHANES pregnancy status codes,true,present,derive nonpregnant domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Demographics,demo,survey_psu,SDMVPSU,SDMVPSU,masked variance unit,true,present,none,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Demographics,demo,survey_strata,SDMVSTRA,SDMVSTRA,masked variance stratum,true,present,none,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr1tot,day1_recall_status,DR1DRSTZ,DR1DRSTZ,1 reliable,true,present,require 1 in domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr2tot,day2_recall_status,DR2DRSTZ,DR2DRSTZ,1 reliable,true,present,require 1 in domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr2tot,number_of_recalls,DRDINT,DRDINT,2 means two recalls,true,present,require 2 in domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr1tot,day1_energy,DR1TKCAL,day1_energy_kcal,kcal/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr2tot,day2_energy,DR2TKCAL,day2_energy_kcal,kcal/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr1tot,day1_sodium,DR1TSODI,day1_sodium_mg,mg/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr2tot,day2_sodium,DR2TSODI,day2_sodium_mg,mg/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr1tot,day1_saturated_fat,DR1TSFAT,day1_satfat_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr2tot,day2_saturated_fat,DR2TSFAT,day2_satfat_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr1tot,day1_monounsaturated_fat,DR1TMFAT,day1_mufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr2tot,day2_monounsaturated_fat,DR2TMFAT,day2_mufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr1tot,day1_polyunsaturated_fat,DR1TPFAT,day1_pufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr2tot,day2_polyunsaturated_fat,DR2TPFAT,day2_pufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr1tot,day1_alcohol_nutrient,DR1TALCO,day1_alcohol_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr2tot,day2_alcohol_nutrient,DR2TALCO,day2_alcohol_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Questionnaire,bpq,hypertension_diagnosed,BPQ020,bpq020,NHANES yes/no codes,true,present,yes_no_na,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Questionnaire,bpq,hypertension_medication,BPQ050A,bpq050a,NHANES yes/no codes,true,present,yes_no_na,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Questionnaire,diq,diabetes_diagnosed,DIQ010,diq010,NHANES yes/no/borderline codes,true,present,derive diabetes,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Questionnaire,diq,insulin_current,DIQ050,diq050,NHANES yes/no codes,true,present,derive diabetes,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Questionnaire,smq,smoked_100_cigarettes,SMQ020,smq020,NHANES yes/no codes,true,present,derive smoking,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Questionnaire,smq,current_smoking,SMQ040,smq040,every day/some days/not at all,true,present,derive smoking,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Questionnaire,kiq_u,kidney_stone_history,KIQ026,kiq026,NHANES yes/no codes,true,present,yes_no_na,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Laboratory,ghb,hba1c,LBXGH,hba1c,percent,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Laboratory,tchol,total_cholesterol,LBXTC,total_cholesterol,mg/dL,true,present,numeric; fixed LBXTC,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Laboratory,cbc,platelet_count,LBXPLTSI,platelet,10^3 cells/uL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Laboratory,cbc,absolute_neutrophils,LBDNENO,neutrophil,10^3 cells/uL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Laboratory,cbc,absolute_lymphocytes,LBDLYMNO,lymphocyte,10^3 cells/uL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Laboratory,biopro,serum_creatinine,LBXSCR,creatinine_raw,mg/dL,true,present,laboratory_bridge,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Laboratory,biopro,serum_uric_acid,LBXSUA,uric_acid_raw,mg/dL,true,present,laboratory_bridge,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Laboratory,alb_cr,urine_albumin,URXUMA,urine_albumin_ug_ml,ug/mL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Laboratory,alb_cr,urine_creatinine,URXUCR,urine_creatinine_mg_dl,mg/dL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr1iff,food_code_day1,DR1IFDCD,food_code_day1,FNDDS code,true,present,character key to cycle manifest,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr1iff,food_grams_day1,DR1IGRMS,food_grams_day1,g,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr1iff,food_energy_day1,DR1IKCAL,food_energy_day1,kcal,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr1iff,food_source_day1,DR1FS,food_source_day1,NHANES source code,true,present,retain for adjudication,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr1iff,combination_type_day1,DR1CCMTX,combination_type_day1,NHANES combination type,true,present,retain for adjudication,08_测试/schema_inventory.csv and CDC IFF codebook,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr2iff,food_code_day2,DR2IFDCD,food_code_day2,FNDDS code,true,present,character key to cycle manifest,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr2iff,food_grams_day2,DR2IGRMS,food_grams_day2,g,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr2iff,food_energy_day2,DR2IKCAL,food_energy_day2,kcal,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr2iff,food_source_day2,DR2FS,food_source_day2,NHANES source code,true,present,retain for adjudication,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr2iff,combination_type_day2,DR2CCMTX,combination_type_day2,NHANES combination type,true,present,retain for adjudication,08_测试/schema_inventory.csv and CDC IFF codebook,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Dietary,dr2tot,diet_weight,WTDR2D,analysis_weight,regular two-day dietary weight,true,present,divide by 13.2 pooled years,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Questionnaire,diq,oral_hypoglycemic_medication,DIQ070,diq070,NHANES yes/no codes,true,present,derive diabetes,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Questionnaire,alq,alcohol_frequency,ALQ101;ALQ110;ALQ120Q;ALQ120U;ALQ130,alcohol,cycle-specific questionnaire,true,present,derive weekly categories using frozen structural-skip map,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Questionnaire,paq,physical_activity,PAQ605;PAQ620;PAQ635;PAQ650;PAQ665,physical_activity,five NHANES yes/no items,true,present,any valid yes active; all valid no inactive; otherwise missing,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Examination,bpx,blood_pressure,BPXSY1;BPXSY2;BPXSY3;BPXSY4;BPXDI1;BPXDI2;BPXDI3;BPXDI4,mean_sbp;mean_dbp,mmHg,true,present,mean of available valid readings,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

G,2011-2012,Laboratory,alb_cr,albumin_creatinine_ratio,URDACT,acr_mg_g,mg/g,true,present,use released URDACT; first random urine,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Demographics,demo,participant_id,SEQN,SEQN,identifier,true,present,none,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Demographics,demo,sex,RIAGENDR,sex,codes 1 male 2 female,true,present,factor labels,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Demographics,demo,age,RIDAGEYR,age,years,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Demographics,demo,race_ethnicity,RIDRETH1,race_ethnicity,NHANES five-category codes,true,present,factor 1:5,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Demographics,demo,education,DMDEDUC2,education,NHANES adult education codes,true,present,1-2 below HS;3 HS/GED;4-5 above HS,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Demographics,demo,pir,INDFMPIR,pir,ratio,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Demographics,demo,pregnancy,RIDEXPRG,ridexprg,NHANES pregnancy status codes,true,present,derive nonpregnant domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Demographics,demo,survey_psu,SDMVPSU,SDMVPSU,masked variance unit,true,present,none,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Demographics,demo,survey_strata,SDMVSTRA,SDMVSTRA,masked variance stratum,true,present,none,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr1tot,day1_recall_status,DR1DRSTZ,DR1DRSTZ,1 reliable,true,present,require 1 in domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr2tot,day2_recall_status,DR2DRSTZ,DR2DRSTZ,1 reliable,true,present,require 1 in domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr2tot,number_of_recalls,DRDINT,DRDINT,2 means two recalls,true,present,require 2 in domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr1tot,day1_energy,DR1TKCAL,day1_energy_kcal,kcal/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr2tot,day2_energy,DR2TKCAL,day2_energy_kcal,kcal/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr1tot,day1_sodium,DR1TSODI,day1_sodium_mg,mg/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr2tot,day2_sodium,DR2TSODI,day2_sodium_mg,mg/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr1tot,day1_saturated_fat,DR1TSFAT,day1_satfat_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr2tot,day2_saturated_fat,DR2TSFAT,day2_satfat_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr1tot,day1_monounsaturated_fat,DR1TMFAT,day1_mufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr2tot,day2_monounsaturated_fat,DR2TMFAT,day2_mufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr1tot,day1_polyunsaturated_fat,DR1TPFAT,day1_pufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr2tot,day2_polyunsaturated_fat,DR2TPFAT,day2_pufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr1tot,day1_alcohol_nutrient,DR1TALCO,day1_alcohol_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr2tot,day2_alcohol_nutrient,DR2TALCO,day2_alcohol_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Questionnaire,bpq,hypertension_diagnosed,BPQ020,bpq020,NHANES yes/no codes,true,present,yes_no_na,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Questionnaire,bpq,hypertension_medication,BPQ050A,bpq050a,NHANES yes/no codes,true,present,yes_no_na,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Questionnaire,diq,diabetes_diagnosed,DIQ010,diq010,NHANES yes/no/borderline codes,true,present,derive diabetes,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Questionnaire,diq,insulin_current,DIQ050,diq050,NHANES yes/no codes,true,present,derive diabetes,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Questionnaire,smq,smoked_100_cigarettes,SMQ020,smq020,NHANES yes/no codes,true,present,derive smoking,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Questionnaire,smq,current_smoking,SMQ040,smq040,every day/some days/not at all,true,present,derive smoking,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Questionnaire,kiq_u,kidney_stone_history,KIQ026,kiq026,NHANES yes/no codes,true,present,yes_no_na,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Laboratory,ghb,hba1c,LBXGH,hba1c,percent,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Laboratory,tchol,total_cholesterol,LBXTC,total_cholesterol,mg/dL,true,present,numeric; fixed LBXTC,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Laboratory,cbc,platelet_count,LBXPLTSI,platelet,10^3 cells/uL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Laboratory,cbc,absolute_neutrophils,LBDNENO,neutrophil,10^3 cells/uL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Laboratory,cbc,absolute_lymphocytes,LBDLYMNO,lymphocyte,10^3 cells/uL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Laboratory,biopro,serum_creatinine,LBXSCR,creatinine_raw,mg/dL,true,present,laboratory_bridge,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Laboratory,biopro,serum_uric_acid,LBXSUA,uric_acid_raw,mg/dL,true,present,laboratory_bridge,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Laboratory,alb_cr,urine_albumin,URXUMA,urine_albumin_ug_ml,ug/mL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Laboratory,alb_cr,urine_creatinine,URXUCR,urine_creatinine_mg_dl,mg/dL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr1iff,food_code_day1,DR1IFDCD,food_code_day1,FNDDS code,true,present,character key to cycle manifest,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr1iff,food_grams_day1,DR1IGRMS,food_grams_day1,g,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr1iff,food_energy_day1,DR1IKCAL,food_energy_day1,kcal,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr1iff,food_source_day1,DR1FS,food_source_day1,NHANES source code,true,present,retain for adjudication,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr1iff,combination_type_day1,DR1CCMTX,combination_type_day1,NHANES combination type,true,present,retain for adjudication,08_测试/schema_inventory.csv and CDC IFF codebook,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr2iff,food_code_day2,DR2IFDCD,food_code_day2,FNDDS code,true,present,character key to cycle manifest,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr2iff,food_grams_day2,DR2IGRMS,food_grams_day2,g,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr2iff,food_energy_day2,DR2IKCAL,food_energy_day2,kcal,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr2iff,food_source_day2,DR2FS,food_source_day2,NHANES source code,true,present,retain for adjudication,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr2iff,combination_type_day2,DR2CCMTX,combination_type_day2,NHANES combination type,true,present,retain for adjudication,08_测试/schema_inventory.csv and CDC IFF codebook,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Dietary,dr2tot,diet_weight,WTDR2D,analysis_weight,regular two-day dietary weight,true,present,divide by 13.2 pooled years,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Questionnaire,diq,oral_hypoglycemic_medication,DIQ070,diq070,NHANES yes/no codes,true,present,derive diabetes,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Questionnaire,alq,alcohol_frequency,ALQ101;ALQ110;ALQ120Q;ALQ120U;ALQ130,alcohol,cycle-specific questionnaire,true,present,derive weekly categories using frozen structural-skip map,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Questionnaire,paq,physical_activity,PAQ605;PAQ620;PAQ635;PAQ650;PAQ665,physical_activity,five NHANES yes/no items,true,present,any valid yes active; all valid no inactive; otherwise missing,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Examination,bpx,blood_pressure,BPXSY1;BPXSY2;BPXSY3;BPXSY4;BPXDI1;BPXDI2;BPXDI3;BPXDI4,mean_sbp;mean_dbp,mmHg,true,present,mean of available valid readings,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

H,2013-2014,Laboratory,alb_cr,albumin_creatinine_ratio,URDACT,acr_mg_g,mg/g,true,present,use released URDACT; first random urine,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Demographics,demo,participant_id,SEQN,SEQN,identifier,true,present,none,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Demographics,demo,sex,RIAGENDR,sex,codes 1 male 2 female,true,present,factor labels,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Demographics,demo,age,RIDAGEYR,age,years,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Demographics,demo,race_ethnicity,RIDRETH1,race_ethnicity,NHANES five-category codes,true,present,factor 1:5,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Demographics,demo,education,DMDEDUC2,education,NHANES adult education codes,true,present,1-2 below HS;3 HS/GED;4-5 above HS,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Demographics,demo,pir,INDFMPIR,pir,ratio,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Demographics,demo,pregnancy,RIDEXPRG,ridexprg,NHANES pregnancy status codes,true,present,derive nonpregnant domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Demographics,demo,survey_psu,SDMVPSU,SDMVPSU,masked variance unit,true,present,none,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Demographics,demo,survey_strata,SDMVSTRA,SDMVSTRA,masked variance stratum,true,present,none,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr1tot,day1_recall_status,DR1DRSTZ,DR1DRSTZ,1 reliable,true,present,require 1 in domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr2tot,day2_recall_status,DR2DRSTZ,DR2DRSTZ,1 reliable,true,present,require 1 in domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr2tot,number_of_recalls,DRDINT,DRDINT,2 means two recalls,true,present,require 2 in domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr1tot,day1_energy,DR1TKCAL,day1_energy_kcal,kcal/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr2tot,day2_energy,DR2TKCAL,day2_energy_kcal,kcal/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr1tot,day1_sodium,DR1TSODI,day1_sodium_mg,mg/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr2tot,day2_sodium,DR2TSODI,day2_sodium_mg,mg/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr1tot,day1_saturated_fat,DR1TSFAT,day1_satfat_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr2tot,day2_saturated_fat,DR2TSFAT,day2_satfat_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr1tot,day1_monounsaturated_fat,DR1TMFAT,day1_mufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr2tot,day2_monounsaturated_fat,DR2TMFAT,day2_mufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr1tot,day1_polyunsaturated_fat,DR1TPFAT,day1_pufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr2tot,day2_polyunsaturated_fat,DR2TPFAT,day2_pufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr1tot,day1_alcohol_nutrient,DR1TALCO,day1_alcohol_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr2tot,day2_alcohol_nutrient,DR2TALCO,day2_alcohol_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Questionnaire,bpq,hypertension_diagnosed,BPQ020,bpq020,NHANES yes/no codes,true,present,yes_no_na,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Questionnaire,bpq,hypertension_medication,BPQ050A,bpq050a,NHANES yes/no codes,true,present,yes_no_na,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Questionnaire,diq,diabetes_diagnosed,DIQ010,diq010,NHANES yes/no/borderline codes,true,present,derive diabetes,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Questionnaire,diq,insulin_current,DIQ050,diq050,NHANES yes/no codes,true,present,derive diabetes,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Questionnaire,smq,smoked_100_cigarettes,SMQ020,smq020,NHANES yes/no codes,true,present,derive smoking,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Questionnaire,smq,current_smoking,SMQ040,smq040,every day/some days/not at all,true,present,derive smoking,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Questionnaire,kiq_u,kidney_stone_history,KIQ026,kiq026,NHANES yes/no codes,true,present,yes_no_na,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Laboratory,ghb,hba1c,LBXGH,hba1c,percent,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Laboratory,tchol,total_cholesterol,LBXTC,total_cholesterol,mg/dL,true,present,numeric; fixed LBXTC,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Laboratory,cbc,platelet_count,LBXPLTSI,platelet,10^3 cells/uL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Laboratory,cbc,absolute_neutrophils,LBDNENO,neutrophil,10^3 cells/uL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Laboratory,cbc,absolute_lymphocytes,LBDLYMNO,lymphocyte,10^3 cells/uL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Laboratory,biopro,serum_creatinine,LBXSCR,creatinine_raw,mg/dL,true,present,laboratory_bridge,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Laboratory,biopro,serum_uric_acid,LBXSUA,uric_acid_raw,mg/dL,true,present,laboratory_bridge,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Laboratory,alb_cr,urine_albumin,URXUMA,urine_albumin_ug_ml,ug/mL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Laboratory,alb_cr,urine_creatinine,URXUCR,urine_creatinine_mg_dl,mg/dL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr1iff,food_code_day1,DR1IFDCD,food_code_day1,FNDDS code,true,present,character key to cycle manifest,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr1iff,food_grams_day1,DR1IGRMS,food_grams_day1,g,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr1iff,food_energy_day1,DR1IKCAL,food_energy_day1,kcal,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr1iff,food_source_day1,DR1FS,food_source_day1,NHANES source code,true,present,retain for adjudication,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr1iff,combination_type_day1,DR1CCMTX,combination_type_day1,NHANES combination type,true,present,retain for adjudication,08_测试/schema_inventory.csv and CDC IFF codebook,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr2iff,food_code_day2,DR2IFDCD,food_code_day2,FNDDS code,true,present,character key to cycle manifest,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr2iff,food_grams_day2,DR2IGRMS,food_grams_day2,g,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr2iff,food_energy_day2,DR2IKCAL,food_energy_day2,kcal,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr2iff,food_source_day2,DR2FS,food_source_day2,NHANES source code,true,present,retain for adjudication,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr2iff,combination_type_day2,DR2CCMTX,combination_type_day2,NHANES combination type,true,present,retain for adjudication,08_测试/schema_inventory.csv and CDC IFF codebook,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Dietary,dr2tot,diet_weight,WTDR2D,analysis_weight,regular two-day dietary weight,true,present,divide by 13.2 pooled years,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Questionnaire,diq,oral_hypoglycemic_medication,DIQ070,diq070,NHANES yes/no codes,true,present,derive diabetes,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Questionnaire,alq,alcohol_frequency,ALQ101;ALQ110;ALQ120Q;ALQ120U;ALQ130,alcohol,cycle-specific questionnaire,true,present,derive weekly categories using frozen structural-skip map,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Questionnaire,paq,physical_activity,PAQ605;PAQ620;PAQ635;PAQ650;PAQ665,physical_activity,five NHANES yes/no items,true,present,any valid yes active; all valid no inactive; otherwise missing,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Examination,bpx,blood_pressure,BPXSY1;BPXSY2;BPXSY3;BPXSY4;BPXDI1;BPXDI2;BPXDI3;BPXDI4,mean_sbp;mean_dbp,mmHg,true,present,mean of available valid readings,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

I,2015-2016,Laboratory,alb_cr,albumin_creatinine_ratio,URDACT,acr_mg_g,mg/g,true,present,use released URDACT; first random urine,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Demographics,demo,participant_id,SEQN,SEQN,identifier,true,present,none,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Demographics,demo,sex,RIAGENDR,sex,codes 1 male 2 female,true,present,factor labels,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Demographics,demo,age,RIDAGEYR,age,years,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Demographics,demo,race_ethnicity,RIDRETH1,race_ethnicity,NHANES five-category codes,true,present,factor 1:5,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Demographics,demo,education,DMDEDUC2,education,NHANES adult education codes,true,present,1-2 below HS;3 HS/GED;4-5 above HS,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Demographics,demo,pir,INDFMPIR,pir,ratio,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Demographics,demo,pregnancy,RIDEXPRG,ridexprg,NHANES pregnancy status codes,true,present,derive nonpregnant domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Demographics,demo,survey_psu,SDMVPSU,SDMVPSU,masked variance unit,true,present,none,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Demographics,demo,survey_strata,SDMVSTRA,SDMVSTRA,masked variance stratum,true,present,none,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr1tot,day1_recall_status,DR1DRSTZ,DR1DRSTZ,1 reliable,true,present,require 1 in domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr2tot,day2_recall_status,DR2DRSTZ,DR2DRSTZ,1 reliable,true,present,require 1 in domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr2tot,number_of_recalls,DRDINT,DRDINT,2 means two recalls,true,present,require 2 in domain,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr1tot,day1_energy,DR1TKCAL,day1_energy_kcal,kcal/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr2tot,day2_energy,DR2TKCAL,day2_energy_kcal,kcal/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr1tot,day1_sodium,DR1TSODI,day1_sodium_mg,mg/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr2tot,day2_sodium,DR2TSODI,day2_sodium_mg,mg/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr1tot,day1_saturated_fat,DR1TSFAT,day1_satfat_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr2tot,day2_saturated_fat,DR2TSFAT,day2_satfat_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr1tot,day1_monounsaturated_fat,DR1TMFAT,day1_mufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr2tot,day2_monounsaturated_fat,DR2TMFAT,day2_mufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr1tot,day1_polyunsaturated_fat,DR1TPFAT,day1_pufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr2tot,day2_polyunsaturated_fat,DR2TPFAT,day2_pufa_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr1tot,day1_alcohol_nutrient,DR1TALCO,day1_alcohol_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr2tot,day2_alcohol_nutrient,DR2TALCO,day2_alcohol_g,g/day,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Questionnaire,bpq,hypertension_diagnosed,BPQ020,bpq020,NHANES yes/no codes,true,present,yes_no_na,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Questionnaire,bpq,hypertension_medication,BPQ050A,bpq050a,NHANES yes/no codes,true,present,yes_no_na,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Questionnaire,diq,diabetes_diagnosed,DIQ010,diq010,NHANES yes/no/borderline codes,true,present,derive diabetes,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Questionnaire,diq,insulin_current,DIQ050,diq050,NHANES yes/no codes,true,present,derive diabetes,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Questionnaire,smq,smoked_100_cigarettes,SMQ020,smq020,NHANES yes/no codes,true,present,derive smoking,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Questionnaire,smq,current_smoking,SMQ040,smq040,every day/some days/not at all,true,present,derive smoking,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Questionnaire,kiq_u,kidney_stone_history,KIQ026,kiq026,NHANES yes/no codes,true,present,yes_no_na,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Laboratory,ghb,hba1c,LBXGH,hba1c,percent,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Laboratory,tchol,total_cholesterol,LBXTC,total_cholesterol,mg/dL,true,present,numeric; fixed LBXTC,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Laboratory,cbc,platelet_count,LBXPLTSI,platelet,10^3 cells/uL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Laboratory,cbc,absolute_neutrophils,LBDNENO,neutrophil,10^3 cells/uL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Laboratory,cbc,absolute_lymphocytes,LBDLYMNO,lymphocyte,10^3 cells/uL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Laboratory,biopro,serum_creatinine,LBXSCR,creatinine_raw,mg/dL,true,present,laboratory_bridge,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Laboratory,biopro,serum_uric_acid,LBXSUA,uric_acid_raw,mg/dL,true,present,laboratory_bridge,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Laboratory,alb_cr,urine_albumin,URXUMA,urine_albumin_ug_ml,ug/mL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Laboratory,alb_cr,urine_creatinine,URXUCR,urine_creatinine_mg_dl,mg/dL,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr1iff,food_code_day1,DR1IFDCD,food_code_day1,FNDDS code,true,present,character key to cycle manifest,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr1iff,food_grams_day1,DR1IGRMS,food_grams_day1,g,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr1iff,food_energy_day1,DR1IKCAL,food_energy_day1,kcal,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr1iff,food_source_day1,DR1FS,food_source_day1,NHANES source code,true,present,retain for adjudication,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr1iff,combination_type_day1,DR1CCMTX,combination_type_day1,NHANES combination type,true,present,retain for adjudication,08_测试/schema_inventory.csv and CDC IFF codebook,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr2iff,food_code_day2,DR2IFDCD,food_code_day2,FNDDS code,true,present,character key to cycle manifest,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr2iff,food_grams_day2,DR2IGRMS,food_grams_day2,g,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr2iff,food_energy_day2,DR2IKCAL,food_energy_day2,kcal,true,present,numeric,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr2iff,food_source_day2,DR2FS,food_source_day2,NHANES source code,true,present,retain for adjudication,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr2iff,combination_type_day2,DR2CCMTX,combination_type_day2,NHANES combination type,true,present,retain for adjudication,08_测试/schema_inventory.csv and CDC IFF codebook,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Dietary,dr2tot,diet_weight,WTDR2DPP,analysis_weight,prepandemic special two-day dietary weight,true,present,divide by 13.2 pooled years,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Questionnaire,diq,oral_hypoglycemic_medication,DIQ070,diq070,NHANES yes/no codes,true,present,derive diabetes,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Questionnaire,alq,alcohol_frequency,ALQ111;ALQ121;ALQ130,alcohol,cycle-specific questionnaire,true,present,derive weekly categories using frozen structural-skip map,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Questionnaire,paq,physical_activity,PAQ605;PAQ620;PAQ635;PAQ650;PAQ665,physical_activity,five NHANES yes/no items,true,present,any valid yes active; all valid no inactive; otherwise missing,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Examination,bpxo,blood_pressure,BPXOSY1;BPXOSY2;BPXOSY3;BPXODI1;BPXODI2;BPXODI3,mean_sbp;mean_dbp,mmHg,true,present,mean of available valid readings,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

P,2017-March 2020,Laboratory,alb_cr,albumin_creatinine_ratio,URDACT,acr_mg_g,mg/g,true,present,use released URDACT; first random urine,08_测试/schema_inventory.csv and frozen REV22 code contract,final,2026-08-14.rev22-confirmed-plan

)--------------------"

embedded_config[["Table1_layout.csv"]] <- r"--------------------(Panel,Section,Characteristic,Level,Overall,Q1,Q2,Q3,Q4,panel_rank,section_rank,group_rank,level_rank,P_value
A. NHANES,Sample and exposure,Unweighted participants,n,,,,,,0,0,0,0,
A. NHANES,Sample and exposure,"UPF quartile range, %",Range,,,,,,0,0,1,0,
A. NHANES,Sample and exposure,"UPF, % of non-water food weight",Mean (SE),,,,,,0,0,2,0,
A. NHANES,Demographic and socioeconomic characteristics,"Age, years",Mean (SE),,,,,,0,1,3,0,
A. NHANES,Demographic and socioeconomic characteristics,Age group,20–44 years,,,,,,0,1,4,0,
A. NHANES,Demographic and socioeconomic characteristics,Age group,45–64 years,,,,,,0,1,4,1,
A. NHANES,Demographic and socioeconomic characteristics,Age group,≥65 years,,,,,,0,1,4,2,
A. NHANES,Demographic and socioeconomic characteristics,Sex,Female,,,,,,0,1,5,0,
A. NHANES,Demographic and socioeconomic characteristics,Sex,Male,,,,,,0,1,5,1,
A. NHANES,Demographic and socioeconomic characteristics,Race and ethnicity,Mexican American,,,,,,0,1,6,0,
A. NHANES,Demographic and socioeconomic characteristics,Race and ethnicity,Non-Hispanic Black,,,,,,0,1,6,1,
A. NHANES,Demographic and socioeconomic characteristics,Race and ethnicity,Non-Hispanic White,,,,,,0,1,6,2,
A. NHANES,Demographic and socioeconomic characteristics,Race and ethnicity,Other Hispanic,,,,,,0,1,6,3,
A. NHANES,Demographic and socioeconomic characteristics,Race and ethnicity,Other/Multiracial,,,,,,0,1,6,4,
A. NHANES,Demographic and socioeconomic characteristics,Education,Below high school,,,,,,0,1,7,0,
A. NHANES,Demographic and socioeconomic characteristics,Education,High school/GED,,,,,,0,1,7,1,
A. NHANES,Demographic and socioeconomic characteristics,Education,Above high school,,,,,,0,1,7,2,
A. NHANES,Lifestyle and anthropometric characteristics,Smoking status,Never,,,,,,0,2,8,0,
A. NHANES,Lifestyle and anthropometric characteristics,Smoking status,Former,,,,,,0,2,8,1,
A. NHANES,Lifestyle and anthropometric characteristics,Smoking status,Current,,,,,,0,2,8,2,
A. NHANES,Lifestyle and anthropometric characteristics,Alcohol use,None or very low,,,,,,0,2,9,0,
A. NHANES,Lifestyle and anthropometric characteristics,Alcohol use,Low to moderate,,,,,,0,2,9,1,
A. NHANES,Lifestyle and anthropometric characteristics,Alcohol use,Higher,,,,,,0,2,9,2,
A. NHANES,Lifestyle and anthropometric characteristics,Physical activity,Inactive,,,,,,0,2,10,0,
A. NHANES,Lifestyle and anthropometric characteristics,Physical activity,Insufficient,,,,,,0,2,10,1,
A. NHANES,Lifestyle and anthropometric characteristics,Physical activity,Sufficient,,,,,,0,2,10,2,
A. NHANES,Lifestyle and anthropometric characteristics,"Body mass index, kg/m²",Mean (SE),,,,,,0,2,11,0,
A. NHANES,Lifestyle and anthropometric characteristics,Body mass index category,<25 kg/m²,,,,,,0,2,12,0,
A. NHANES,Lifestyle and anthropometric characteristics,Body mass index category,25–<30 kg/m²,,,,,,0,2,12,1,
A. NHANES,Lifestyle and anthropometric characteristics,Body mass index category,≥30 kg/m²,,,,,,0,2,12,2,
A. NHANES,Dietary profile,"Energy intake, kcal/day",Mean (SE),,,,,,0,3,13,0,
A. NHANES,Dietary profile,"NOVA 1 unprocessed/minimally processed foods, %",Mean (SE),,,,,,0,3,14,0,
A. NHANES,Dietary profile,"NOVA 2 processed culinary ingredients, %",Mean (SE),,,,,,0,3,15,0,
A. NHANES,Dietary profile,"NOVA 3 processed foods, %",Mean (SE),,,,,,0,3,16,0,
A. NHANES,Dietary profile,HEI-2015 total score,Mean (SE),,,,,,0,3,17,0,
A. NHANES,Clinical characteristics,Hypertension,Yes,,,,,,0,4,19,0,
A. NHANES,Clinical characteristics,Diabetes,Yes,,,,,,0,4,20,0,
A. NHANES,Clinical characteristics,High cholesterol,Yes,,,,,,0,4,21,0,
A. NHANES,Kidney measures and outcomes,Chronic kidney disease,Yes,,,,,,0,5,22,0,
A. NHANES,Kidney measures and outcomes,Low eGFR (<60 mL/min/1.73 m²),Yes,,,,,,0,5,23,0,
A. NHANES,Kidney measures and outcomes,High UACR (≥30 mg/g),Yes,,,,,,0,5,24,0,
A. NHANES,Kidney measures and outcomes,Diabetic kidney disease,Yes,,,,,,0,5,25,0,
A. NHANES,Kidney measures and outcomes,Kidney stones,Yes,,,,,,0,5,26,0,
A. NHANES,Kidney measures and outcomes,"eGFR, mL/min/1.73 m²",Mean (SE),,,,,,0,5,27,0,
A. NHANES,Kidney measures and outcomes,"UACR, mg/g",Median (IQR),,,,,,0,5,28,0,
A. NHANES,Kidney measures and outcomes,"Serum uric acid, mg/dL",Mean (SE),,,,,,0,5,29,0,
B. KNHANES,Sample and exposure,Unweighted participants,n,,,,,,1,0,0,0,
B. KNHANES,Sample and exposure,"UPF quartile range, %",Range,,,,,,1,0,1,0,
B. KNHANES,Sample and exposure,"UPF, % of non-water food weight",Mean (SE),,,,,,1,0,2,0,
B. KNHANES,Demographic and socioeconomic characteristics,"Age, years",Mean (SE),,,,,,1,1,3,0,
B. KNHANES,Demographic and socioeconomic characteristics,Age group,20–44 years,,,,,,1,1,4,0,
B. KNHANES,Demographic and socioeconomic characteristics,Age group,45–64 years,,,,,,1,1,4,1,
B. KNHANES,Demographic and socioeconomic characteristics,Age group,≥65 years,,,,,,1,1,4,2,
B. KNHANES,Demographic and socioeconomic characteristics,Sex,Female,,,,,,1,1,5,0,
B. KNHANES,Demographic and socioeconomic characteristics,Sex,Male,,,,,,1,1,5,1,
B. KNHANES,Demographic and socioeconomic characteristics,Education,Below high school,,,,,,1,1,7,0,
B. KNHANES,Demographic and socioeconomic characteristics,Education,High school,,,,,,1,1,7,1,
B. KNHANES,Demographic and socioeconomic characteristics,Education,Above high school,,,,,,1,1,7,2,
B. KNHANES,Lifestyle and anthropometric characteristics,Smoking status,Never,,,,,,1,2,8,0,
B. KNHANES,Lifestyle and anthropometric characteristics,Smoking status,Former,,,,,,1,2,8,1,
B. KNHANES,Lifestyle and anthropometric characteristics,Smoking status,Current,,,,,,1,2,8,2,
B. KNHANES,Lifestyle and anthropometric characteristics,Alcohol use,None or very low,,,,,,1,2,9,0,
B. KNHANES,Lifestyle and anthropometric characteristics,Alcohol use,Low to moderate,,,,,,1,2,9,1,
B. KNHANES,Lifestyle and anthropometric characteristics,Alcohol use,Higher,,,,,,1,2,9,2,
B. KNHANES,Lifestyle and anthropometric characteristics,Physical activity,Yes,,,,,,1,2,10,0,
B. KNHANES,Lifestyle and anthropometric characteristics,"Body mass index, kg/m²",Mean (SE),,,,,,1,2,11,0,
B. KNHANES,Lifestyle and anthropometric characteristics,Body mass index category,<25 kg/m²,,,,,,1,2,12,0,
B. KNHANES,Lifestyle and anthropometric characteristics,Body mass index category,25–<30 kg/m²,,,,,,1,2,12,1,
B. KNHANES,Lifestyle and anthropometric characteristics,Body mass index category,≥30 kg/m²,,,,,,1,2,12,2,
B. KNHANES,Dietary profile,"Energy intake, kcal/day",Mean (SE),,,,,,1,3,13,0,
B. KNHANES,Dietary profile,"NOVA 1 unprocessed/minimally processed foods, %",Mean (SE),,,,,,1,3,14,0,
B. KNHANES,Dietary profile,"NOVA 2 processed culinary ingredients, %",Mean (SE),,,,,,1,3,15,0,
B. KNHANES,Dietary profile,"NOVA 3 processed foods, %",Mean (SE),,,,,,1,3,16,0,
B. KNHANES,Dietary profile,KHEI total score,Mean (SE),,,,,,1,3,18,0,
B. KNHANES,Clinical characteristics,Hypertension,Yes,,,,,,1,4,19,0,
B. KNHANES,Clinical characteristics,Diabetes,Yes,,,,,,1,4,20,0,
B. KNHANES,Clinical characteristics,High cholesterol,Yes,,,,,,1,4,21,0,
B. KNHANES,Kidney measures and outcomes,Chronic kidney disease,Yes,,,,,,1,5,22,0,
B. KNHANES,Kidney measures and outcomes,Low eGFR (<60 mL/min/1.73 m²),Yes,,,,,,1,5,23,0,
B. KNHANES,Kidney measures and outcomes,High UACR (≥30 mg/g),Yes,,,,,,1,5,24,0,
B. KNHANES,Kidney measures and outcomes,Diabetic kidney disease,Yes,,,,,,1,5,25,0,
B. KNHANES,Kidney measures and outcomes,"eGFR, mL/min/1.73 m²",Mean (SE),,,,,,1,5,27,0,
B. KNHANES,Kidney measures and outcomes,"UACR, mg/g",Median (IQR),,,,,,1,5,28,0,
B. KNHANES,Kidney measures and outcomes,"Serum uric acid, mg/dL",Mean (SE),,,,,,1,5,29,0,
)--------------------"

embedded_config[["download_sources.csv"]] <- r"--------------------(file,url
fped_dr1tot_1720.sas7bdat,https://www.ars.usda.gov/ARSUserFiles/80400530/apps/FPED_DR1TOT_1720_sas.EXE
fped_dr2tot_1720.sas7bdat,https://www.ars.usda.gov/ARSUserFiles/80400530/apps/FPED_DR2TOT_1720_sas.EXE
fped_dr1tot_1516.sas7bdat,https://www.ars.usda.gov/ARSUserFiles/80400530/apps/FPED_DR1TOT_1516_sas.exe
fped_dr2tot_1516.sas7bdat,https://www.ars.usda.gov/ARSUserFiles/80400530/apps/FPED_DR2TOT_1516_sas.exe
fped_dr1tot_1314.sas7bdat,https://www.ars.usda.gov/ARSUserFiles/80400530/apps/FPED_DR1TOT_1314_sas.exe
fped_dr2tot_1314.sas7bdat,https://www.ars.usda.gov/ARSUserFiles/80400530/apps/FPED_DR2TOT_1314_sas.exe
fped_dr1tot_1112.sas7bdat,https://www.ars.usda.gov/ARSUserFiles/80400530/apps/FPED_DR1TOT_1112_sas.exe
fped_dr2tot_1112.sas7bdat,https://www.ars.usda.gov/ARSUserFiles/80400530/apps/FPED_DR2TOT_1112_sas.exe
fped_dr1tot_0910.sas7bdat,https://www.ars.usda.gov/ARSUserFiles/80400530/apps/FPED_DR1TOT_0910_sas.exe
fped_dr2tot_0910.sas7bdat,https://www.ars.usda.gov/ARSUserFiles/80400530/apps/FPED_DR2TOT_0910_sas.exe
fped_dr1tot_0708.sas7bdat,https://www.ars.usda.gov/ARSUserFiles/80400530/apps/FPED_DR1TOT_0708_sas.exe
fped_dr2tot_0708.sas7bdat,https://www.ars.usda.gov/ARSUserFiles/80400530/apps/FPED_DR2TOT_0708_sas.exe
)--------------------"

embedded_config[["required_packages.csv"]] <- r"--------------------(Package,Version
VGAM,1.1-14
broom,1.0.11
cowplot,1.1.3
data.table,1.17.8
digest,0.6.39
dplyr,1.1.4
forestploter,1.1.4
future,1.68.0
future.apply,1.20.0
ggplot2,4.0.2
grDevices,4.5.3
gtsummary,2.5.1
haven,2.5.5
jsonlite,2.0.0
labelled,2.16.0
mice,3.19.0
mitml,0.4-5
mitools,2.4
nnet,7.3-20
openxlsx,4.2.8.1
parallel,4.5.3
parallelly,1.45.1
patchwork,1.3.2
posterior,1.7.0
purrr,1.2.0
readr,2.1.6
splines,4.5.3
stats,4.5.3
stringr,1.6.0
survey,4.5
svyVGAM,1.3
testthat,3.3.1
tidyr,1.3.1
tidyselect,1.2.1
tools,4.5.3
utils,4.5.3
withr,3.0.2
)--------------------"

embedded_config[["required_raw_files.csv"]] <- r"--------------------("input_root","relative_path","bytes","sha256"
"nhanes","2007-2008/Demographics/demo_e.xpt",3498080,"02b8e97703d04bc99769c751881bece51ac9ba68cf18a2992a45d3bdfb856147"
"nhanes","2007-2008/Dietary/dr1tot_e.xpt",12831440,"c177e58fe5d06e8a5ecc83bbc45cf3e4b280d3efe20de173c6d32bb2a3fac0ac"
"nhanes","2007-2008/Dietary/dr2tot_e.xpt",6494400,"9980449a908bb2f21a9bc9f68117662bbf109f87bb36acc76480afd880e10ac9"
"nhanes","2007-2008/Questionnaire/bpq_e.xpt",998480,"b802b5855d48fc230c734083ea1e4d4a9d84e1b36bf5e78fcd725c17d6ce77df"
"nhanes","2007-2008/Questionnaire/diq_e.xpt",2789600,"fb4d698fbec26bc8a09b2cd0f3019f815114742dd2c37ccac1aaf5e0123394c5"
"nhanes","2007-2008/Questionnaire/smq_e.xpt",2449520,"881865ee05f4b2be0c03c761f00ad8e355f2cc5e7d3be41666edef6f3dca5862"
"nhanes","2007-2008/Questionnaire/kiq_u_e.xpt",762640,"efaa9108e033b47754b24d7e4bc5b8b373c6844de3e3a82035ca82320d9f44e8"
"nhanes","2007-2008/Laboratory/ghb_e.xpt",111760,"4c0311851825e9e0352c58af193c58d8f7a3af613aa5f522893e3c8f83e53a53"
"nhanes","2007-2008/Laboratory/tchol_e.xpt",196400,"1e137eb4fcccd0c0e104b4b1849f00b09e8a9077e2143e844a09bb673507fff1"
"nhanes","2007-2008/Laboratory/cbc_e.xpt",1567280,"e83558b948eaf1b002e029d69731c01c81935ebf85013358d59d2aa95ad38a26"
"nhanes","2007-2008/Laboratory/biopro_e.xpt",2053360,"6c2cec1e790c09e0319bf5a576a2ab64861d4270679c8113663f83a6b7087ced"
"nhanes","2007-2008/Laboratory/alb_cr_e.xpt",326720,"6bf8018a465499043d5aeb89ac6d38ab4bc25a126ec663979bbbeffe0afdc2d8"
"nhanes","2007-2008/Questionnaire/alq_e.xpt",412960,"3fca970ac6ed096119d47962f3c72157dde7c11b13a499380e9a1bbb33123789"
"nhanes","2007-2008/Questionnaire/paq_e.xpt",1576000,"cc4fe3815c0e3602f69db33a3709bf06b7bca38c21d6b2f74d1bb5c98a5ff35a"
"nhanes","2007-2008/Examination/bpx_e.xpt",2113200,"5a0303e68e0a4d7ebdbd2bc82f87852c31e44f1733782e107219bed9303baed5"
"nhanes","2009-2010/Demographics/demo_f.xpt",3631600,"78c4d4a6e5eafa9da3bbb24a90d140fef5d8c0cda932a7ff8ee609288d4fae07"
"nhanes","2009-2010/Dietary/dr1tot_f.xpt",13640000,"7ff258b1cb2020509d0d0ad3afb49797ca5196d085bac3e45cc6df149d93091c"
"nhanes","2009-2010/Dietary/dr2tot_f.xpt",6820400,"214759831d18d9bc0f0b2975fed277a7836e4d35774c5dc655643b582a28adbe"
"nhanes","2009-2010/Questionnaire/bpq_f.xpt",1161040,"e592c987f07002eb309a63dc7ebf6a7eaa83d740a4fd6ae69af15da40c65a21a"
"nhanes","2009-2010/Questionnaire/diq_f.xpt",1620960,"7e90371cbc98e77a9684eff32ffa0e4ee3e8fc04d97581a8773e19cb9fe901fa"
"nhanes","2009-2010/Questionnaire/smq_f.xpt",2580560,"ee85c82ea55046eb7bf491088b5bb9efecbea8aaf66ae8c89123bad23a5f1901"
"nhanes","2009-2010/Questionnaire/kiq_u_f.xpt",798880,"63b648bc50a6c904ddee61efd6246af8703f2ffff7fa4c9f302017c1e517ac46"
"nhanes","2009-2010/Laboratory/ghb_f.xpt",118960,"98505d8aa707bd06b8b1d488d1ecdf13388bffea6ec1b89929726e7d2b32703b"
"nhanes","2009-2010/Laboratory/tchol_f.xpt",207440,"87a237c6bcf8fe623db0ace18c5629bc6e8fddd98d4202a5a2d449958d13c4db"
"nhanes","2009-2010/Laboratory/cbc_f.xpt",1656000,"47e7f693c7899a2d487502ec4a23accadf73001292dc566e8832364e48edd5f7"
"nhanes","2009-2010/Laboratory/biopro_f.xpt",2187200,"d138350f6a9318bb6a411b67f511e8eded677baa6ca8f4cad70dd86a0f0f428a"
"nhanes","2009-2010/Laboratory/alb_cr_f.xpt",758400,"f902e8892db69054e8898022dc1ced19770cc2fa613bf2655455f6a4f1cfbb67"
"nhanes","2009-2010/Questionnaire/alq_f.xpt",438320,"15a46b2310d5ba90e1488df5f84846d97468247b0aae6608407fa9bf2405d0d1"
"nhanes","2009-2010/Questionnaire/paq_f.xpt",1645280,"02d25e8312a36103ff6ed416dfd597fb109f5a299cc951477337ff8fe694aa84"
"nhanes","2009-2010/Examination/bpx_f.xpt",2219280,"d09af545c732e5fb72c437e870823ec5d9efaed1ef3522210dcc1d8f47b9be57"
"nhanes","2011-2012/Demographics/demo_g.xpt",3753760,"eaf0525d1952626885af3e935415a1f66ad62c18698080e7354789c125af252d"
"nhanes","2011-2012/Dietary/dr1tot_g.xpt",12424880,"d5589de4cd987eef86dd13b3b2bc3db9186ab80a87c96a3968421ddef9235c57"
"nhanes","2011-2012/Dietary/dr2tot_g.xpt",6212880,"fade17bfa233609832197f8dee2dc95d86ba320831034f647aec96487a6316ba"
"nhanes","2011-2012/Questionnaire/bpq_g.xpt",743920,"39c5a0ee1d76e6f3df867c0057986a21f4a028038002c43214c5af8bdec83cae"
"nhanes","2011-2012/Questionnaire/diq_g.xpt",3978560,"a9d5475e0cd66d6a7bc30230345ed0172a9abcbc0a7a4145a35025c827a52a87"
"nhanes","2011-2012/Questionnaire/smq_g.xpt",1946960,"69f2375f036a815082b1e97110e79c505b91e0985c61696f8eaacf94828115dc"
"nhanes","2011-2012/Questionnaire/kiq_u_g.xpt",714640,"eb9b52d4071ee1d16c8c4a4ecd125e8a5e8fae72a889703ecd51160d8960dde9"
"nhanes","2011-2012/Laboratory/ghb_g.xpt",105840,"8ff7cc95461fdfd47d2c9640b944a2bcd45926a3c0e64c98a09a4cef6c0f88d7"
"nhanes","2011-2012/Laboratory/tchol_g.xpt",188960,"08682d7bcee71c1a5c5d977ebfe5f55eda43354a4ff76b07c3fbb487075e0731"
"nhanes","2011-2012/Laboratory/cbc_g.xpt",1508320,"f5fd08c8f7564d73429579620f3dab4aec48cfdb0a908e450c0ea849d80589c0"
"nhanes","2011-2012/Laboratory/biopro_g.xpt",1997040,"8c3a0b1bea79ff79c5600ad6a76c4206a7d689479f94ae81da6360058eda5d05"
"nhanes","2011-2012/Laboratory/alb_cr_g.xpt",377040,"46c5733d0bdb1d4bbe841614a3e6ef6b9dab8501f455af4f841cc43a312e579d"
"nhanes","2011-2012/Questionnaire/alq_g.xpt",451360,"5f06ba3506f9d69cedffa0380ad0cecbc15df2d6d9b7a3644ac3e2678aaac5e5"
"nhanes","2011-2012/Questionnaire/paq_g.xpt",1533680,"6b56874f6c3af191aebcb10316c028bc1ae1d4f89d0aafc00cc080c9f533b32a"
"nhanes","2011-2012/Examination/bpx_g.xpt",2021600,"1dd66110f41465d601cd54f1e3c6010baafb9375f09fbcaa3aef2536876146fd"
"nhanes","2013-2014/Demographics/demo_h.xpt",3833200,"f8f0cbb3085a323d4cde22349b164878fea1e64dbc404e65b5815c7816b547d7"
"nhanes","2013-2014/Dietary/dr1tot_h.xpt",13212960,"4fc45d48cf8a4bdcca615efa895e7e05414ca4d549bd75c8115e06abcc86cdc3"
"nhanes","2013-2014/Dietary/dr2tot_h.xpt",6685520,"12fa3550afbd57ec62d1543c668f1736387f3421fc6668369b5fa4c15a369f91"
"nhanes","2013-2014/Questionnaire/bpq_h.xpt",726720,"aa6ddeab80c73b074ccfd04ad010da404ecb91c048476a194594ce9552d03d49"
"nhanes","2013-2014/Questionnaire/diq_h.xpt",4228960,"c74c7ccef65e6997dfac1db1e73bc3f63dec67c70b2e322ffc14f14ce27429a9"
"nhanes","2013-2014/Questionnaire/smq_h.xpt",2170000,"86d74bb8842ce6d60350944417cbd3e9ee39dd5a9fa1bda33426f276b673861f"
"nhanes","2013-2014/Questionnaire/kiq_u_h.xpt",741440,"a12cf2e0e3fbceb1617bae2c09a57edbd5c9bbdf35e6879fe6a52136ba83bef4"
"nhanes","2013-2014/Laboratory/ghb_h.xpt",112720,"0695894ad55ac96f315a8415401977b0856c402d16762115b534bcd5dfeae89e"
"nhanes","2013-2014/Laboratory/tchol_h.xpt",200240,"0867a9d340dc84eec8dc1f08672e7ffd0a364b57a92533973c460acb4f52d643"
"nhanes","2013-2014/Laboratory/cbc_h.xpt",1586640,"08f6fca4aacebacafb2b311c2edfc6897ad92f2c8f8dab14e18bf41acdc4ddee"
"nhanes","2013-2014/Laboratory/biopro_h.xpt",2127760,"862587c4b96b90cae8a79742a93b7beff6dcf9b8ec10cf0cb245c69314d3b26b"
"nhanes","2013-2014/Laboratory/alb_cr_h.xpt",399600,"65adbd6c87541464e69fa2da54fb2c5e9e0c198fa322943a21ccbcacd8c40670"
"nhanes","2013-2014/Questionnaire/alq_h.xpt",476080,"174e4dd21e05f81df8bc786a6e76253a2cdc07503cd7cde4c7dd1c1b9d3e3098"
"nhanes","2013-2014/Questionnaire/paq_h.xpt",7297920,"fa35acbb489594e8e55350f7402a12b5a607805c28ec16246929d33237bd3269"
"nhanes","2013-2014/Examination/bpx_h.xpt",1809600,"b6641efc30ac9b9de328dfbf4110d0d120439aadca2c2bc89f2e77310d1a0c8e"
"nhanes","2015-2016/Demographics/demo_i.xpt",3756480,"c9297c6c37ae8f78f29be9568fa2a03cf3b112616a39afee04030fc775a66a0d"
"nhanes","2015-2016/Dietary/dr1tot_i.xpt",12851440,"185b8f112cdf22175470fe8cc7d4d20800fe28e8ac37c50dc4bf197d40333b6e"
"nhanes","2015-2016/Dietary/dr2tot_i.xpt",6502560,"0940dad7ba10cb705945665de06329affb3a59495d4947b7249602ae18d92643"
"nhanes","2015-2016/Questionnaire/bpq_i.xpt",559120,"7f7e51cd0497eb648f64eaafaf7ed07825d617c32b2cd6853ed249fce1143a38"
"nhanes","2015-2016/Questionnaire/diq_i.xpt",4144720,"e87587479b29f175b63eee5dd40d582837e3e3fe2665503012085eefdb978e0d"
"nhanes","2015-2016/Questionnaire/smq_i.xpt",2681040,"529ce1c407e99dce937ca75d14c8c630e0c9c30a4a9a851760cd2884bbb28fe0"
"nhanes","2015-2016/Questionnaire/kiq_u_i.xpt",735040,"604bc20f6195c1cf0e7f5a9f9c168915c1a042c3e8628b145057fbb00dec4231"
"nhanes","2015-2016/Laboratory/ghb_i.xpt",108960,"e4bc626cd12f6057c7806aef4f85874c9bf1407a7480d6621a4af142260addd2"
"nhanes","2015-2016/Laboratory/tchol_i.xpt",193760,"e7f51cf7c002db01d023278922d89d72318c679967e72ff259ba532b2d72e54d"
"nhanes","2015-2016/Laboratory/cbc_i.xpt",1543440,"d2605788819fb3507d1a002379d69eba1ed6d3edc2b5413a947f2c2481b21278"
"nhanes","2015-2016/Laboratory/biopro_i.xpt",2056320,"b3f07e9cedb1dd7ddc57ba2146eb258a8813c195ebc5bb858a620286a703c8b0"
"nhanes","2015-2016/Laboratory/alb_cr_i.xpt",552800,"af61c16db8d9e96549bf0fb446049ff0631ae1ccf7e2dcd017d28edba84dbe3d"
"nhanes","2015-2016/Questionnaire/alq_i.xpt",460960,"ed445c40ec3f3f993f78321f86e513340d4f131708c7624d88803bde3909102d"
"nhanes","2015-2016/Questionnaire/paq_i.xpt",6973680,"c508dc563a907b2fcd05048f37008197257a3e1c6eb64563623c60459c1afc9d"
"nhanes","2015-2016/Examination/bpx_i.xpt",1607120,"45d5c9fa101ccc6da82241bd45c99611595e699665eed2cfed232e9070b52cc1"
"nhanes","2019-2020/Demographics/p_demo.xpt",3614720,"2e46c6c26bf77cd8989f64011ace12cbf42c0f3e03414eb59acc5328c8f87913"
"official_xpt","P_DR1TOT.XPT",19243440,"02cf0fcb15cdef9442594747a6699c43a65fd280f0f807e82ca0b888b0d97b5e"
"official_xpt","P_DR2TOT.XPT",9736640,"2013d36582131ffc810bb025e11c8a0c6533b5e59197ad2aa76f4d52db81bb5a"
"nhanes","2019-2020/Questionnaire/p_bpq.xpt",899520,"88af3557dd185cee2533c94c5302b361e89a54eb50754704f35d4e8120868c02"
"nhanes","2019-2020/Questionnaire/p_diq.xpt",3361520,"79153f799fc2171792771c4aa250029c62807f14915f1f81405b03233b1e5ae3"
"nhanes","2019-2020/Questionnaire/p_smq.xpt",1428560,"29b7f6c59ab570c042c866312f0eb4329009e85fadbda6f9427a8bf31ecb3a3f"
"nhanes","2019-2020/Questionnaire/p_kiq_u.xpt",1184720,"e1f6cc9b6b06c826731b453a77c8b44acd53f03e254f1a0690a0c6e3598d2e15"
"nhanes","2019-2020/Laboratory/p_ghb.xpt",167600,"dac9e423f56041c2ec46486cb16be16d3347e0a81a575985882ffad5a60e1195"
"nhanes","2019-2020/Laboratory/p_tchol.xpt",294000,"759cd9c408b5d1b91f6cfb5d5e67920a301f7478c0e90ce7adf4bbd277a45ac8"
"nhanes","2019-2020/Laboratory/p_cbc.xpt",2427760,"80851942250bf23cd483aaf07ccdc31a57346c2d8bd3f1ded988e4a2498181d8"
"nhanes","2019-2020/Laboratory/p_biopro.xpt",3420640,"ffb6146d15a158d9f1a4c22b6e192bd99697e8e33ffe30bffc20d2f39a635083"
"nhanes","2019-2020/Laboratory/p_alb_cr.xpt",835600,"c493f6b74cc0608bb7083b92e074fc3070b9e614fd5e1e34d0ce943e0f6116bf"
"nhanes","2019-2020/Questionnaire/p_alq.xpt",719360,"feb3f25e06c6a3627d5699523d26e2ff82a4767109b15665e9eae93b3b38e722"
"nhanes","2019-2020/Questionnaire/p_paq.xpt",1321440,"d0e13cebf96181949e983ab58d14c58689e801e8711d366d2b5daab3c7b9ee65"
"nhanes","2019-2020/Examination/p_bpxo.xpt",1039840,"5edbcb81e99fe8dad38deb8b642570278b200c525d88ffa29ae3ceb4dd38bde4"
"nhanes","2007-2008/Examination/bmx_e.xpt",1800240,"b67b50dbd12c9510b49b6348b595b9fd073c422c0c6eac182f90ac57b491a97d"
"nhanes","2009-2010/Examination/bmx_f.xpt",1890560,"1fa42e07764037674c3f56e1e76f990e10c32245897597de3d731c57b91d869a"
"nhanes","2011-2012/Examination/bmx_g.xpt",1946720,"4814bfc3047ed400b9d43d285f8c3ea7c940ac6489404a9b699579715d158ec3"
"nhanes","2013-2014/Examination/bmx_h.xpt",2045520,"fd5e9fc6e6aab0a4aee6e699f51497bbc9b62101f7f43aee924c473e38fd9442"
"nhanes","2015-2016/Examination/bmx_i.xpt",1989600,"d31da84e14212b4e58e8340598a5b8e2144fac83333563e966eb1b332e4141d6"
"nhanes","2019-2020/Examination/p_bmx.xpt",2520640,"7038d0da3169420a4a3cbaec09ac586a701b38e1c7b925431b60e13ac24fed66"
"nhanes","2007-2008/Dietary/dr1iff_e.xpt",99090720,"aedfffe65f5a845d69099ad170474bfe022e15eb3baa708e741682bf7fa60550"
"nhanes","2007-2008/Dietary/dr2iff_e.xpt",82524560,"6c605cd47b92e161aa2194448653e0d8beaec0fc97b1c048d76441c320accc0b"
"nhanes","2009-2010/Dietary/dr1iff_f.xpt",102686560,"d2e2cf8195695f8c8a3890b62231fa71b344824fd531a86a483da8e1ccd30034"
"nhanes","2009-2010/Dietary/dr2iff_f.xpt",87828560,"f88569d8c7a147690b5ba36300308db3b041e846d55a37499108aa725252ed80"
"nhanes","2011-2012/Dietary/dr1iff_g.xpt",86034720,"bf40537bd378f95377b7c274198ac0e7b3a1fbf70b2ad9d995ccc72943b2c061"
"nhanes","2011-2012/Dietary/dr2iff_g.xpt",76241360,"676f5d0af9aaa82e0b733ed21b0637864b737ffb4187c3739261199fc032fca8"
"nhanes","2013-2014/Dietary/dr1iff_h.xpt",88309280,"2af08841d72b4b4f03fdeddb3e666e1fee2b86d7b0d99d74b3b6f7125b560af7"
"nhanes","2013-2014/Dietary/dr2iff_h.xpt",75664960,"7aaccb4622a6aaf334655812517f02faa666870b4cd0000efe330b6c8effd9b6"
"nhanes","2015-2016/Dietary/dr1iff_i.xpt",81647760,"36b516fa303611f075e80ce1c44360a8c964261adcc67b7aa3f2df8b937532f5"
"nhanes","2015-2016/Dietary/dr2iff_i.xpt",67669440,"00c11f34fb2ef383d8c974d71c98edfdb9855aafc5aee7fc33ee3cea55b26ef1"
"nhanes","2019-2020/Dietary/P_DR1IFF.xpt",123600000,"264a4dfd66455f106a8ac25a745ebb4e065564b049b45000a5f9c44296d48dec"
"official_xpt","P_DR2IFF.XPT",100473120,"6f5807d841fb9e183ba5ad6a9e0ebcef534f96c4ffb8e59afd01293dc7cf55e0"
"knhanes","hn19_all.sas7bdat",101711872,"ef19bdcbc9ccfdd5cd514a872227ff62e821b43f6da0c7f38fce008e42dcb653"
"knhanes","hn19_24rc.sas7bdat",813027328,"9dec9ad71c41a9db1487bc0115770ad0651f42a63f45ba998b0a5f08837affb6"
"knhanes","hn20_all.sas7bdat",83091456,"3d3516534dc7f9c484a6dead367640d786e5194440485491ca0ba9ec6d9c0023"
"knhanes","hn20_24rc.sas7bdat",602882048,"20cedbfabcf24b9185a3eeeb78f4abc91afcff1307fff0a789c00941343baf4f"
"knhanes","hn21_all.sas7bdat",77856768,"dcd0b06cb433b05ee91cc8c8d72783d7aaef9b6edc60b444bc96b1044edb78b0"
"knhanes","hn21_24rc.sas7bdat",637800448,"4b041f2aa46557125991889ea586bd1e3458eec543a1c28e1c4bb8bd7a002650"
"knhanes","hn22_all.sas7bdat",59244544,"971662027b095f8f5f5108c7a441a4a7e14422e56001ccb5377320d796edd553"
"knhanes","hn22_24rc.sas7bdat",641908736,"ad44462585bc5b74abe7cdbf77ddd6aa2384ce91e0b88846cb13808820ff5667"
"knhanes","hn23_all.sas7bdat",58060800,"62b3a67bd1a86fb459c78b404735a182ec0fe03cd4d35f7420d665b5b1e2741c"
"knhanes","hn23_24rc.sas7bdat",843882496,"b7e691a6d2f28280a19dea9d4039840de7a83fd5097b588652c82ae0570ed08b"
"knhanes","hn24_all.sas7bdat",78962688,"ff74cb84432cb1f10ba63d1a3aba54215ef38e47afef29547f83127aea9fc47f"
"knhanes","hn24_24rc.sas7bdat",842231808,"e707b75dbc19dbc9f2b098b9c3243964fa4bc0f7c5caccf446c704470cba23ea"
)--------------------"

read_config <- function(name, ...) read.csv(text=embedded_config[[name]], ...)

write_config <- function(name, path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  writeChar(embedded_config[[name]],path,eos=NULL,useBytes=TRUE)
}

classification_hashes <- c("classifications/KNHANES_NOVA.csv.gz"="8068aee97e4c82a1f421d182e25f9b7fff6d8ee7f6825b043ed55bdd173b6257","classifications/NHANES_NOVA.csv.gz"="14624c7213b14ba4e9710051e52f799d2de6b378e61dd78bd380b90103250998")

analysis_jobs <- list(
  list(id="five_outcome_quartiles",script="NHANES/02_代码/正式五结局/24_quartile_trend_models.R",args=character(0)),
  list(id="five_outcome_sensitivities",script="NHANES/02_代码/正式五结局/26_sensitivity_m3_models.R",args=character(0)),
  list(id="HEI_CKD",script="NHANES/02_代码/正式五结局/33_diet_quality_models.R",args=c("CKD")),
  list(id="HEI_eGFR",script="NHANES/02_代码/正式五结局/33_diet_quality_models.R",args=c("eGFR")),
  list(id="HEI_UACR",script="NHANES/02_代码/正式五结局/33_diet_quality_models.R",args=c("UACR")),
  list(id="HEI_DKD",script="NHANES/02_代码/正式五结局/33_diet_quality_models.R",args=c("DKD")),
  list(id="HEI_kidney_stones",script="NHANES/02_代码/正式五结局/33_diet_quality_models.R",args=c("kidney_stones")),
  list(id="NOVA_CKD",script="NHANES/02_代码/正式五结局/32_nova_substitution_models.R",args=c("CKD")),
  list(id="NOVA_eGFR",script="NHANES/02_代码/正式五结局/32_nova_substitution_models.R",args=c("eGFR")),
  list(id="NOVA_UACR",script="NHANES/02_代码/正式五结局/32_nova_substitution_models.R",args=c("UACR")),
  list(id="NOVA_DKD",script="NHANES/02_代码/正式五结局/32_nova_substitution_models.R",args=c("DKD")),
  list(id="NOVA_kidney_stones",script="NHANES/02_代码/正式五结局/32_nova_substitution_models.R",args=c("kidney_stones")),
  list(id="five_outcome_subgroups",script="NHANES/02_代码/正式五结局/41_subgroup_m50_models.R",args=c("formal")),
  list(id="five_outcome_final_RCS",script="01_code/modules/RCS/12_recalculate_rcs_v3.R",args=c("NHANES")),
  list(id="SUA_quartiles_sens_hyperuricemia",script="01_code/modules/SUA/04_nhanes_serum_uric_acid_extended_core.R",args=character(0)),
  list(id="SUA_subgroups",script="01_code/modules/SUA/06_nhanes_serum_uric_acid_subgroups.R",args=character(0)),
  list(id="SUA_HEI_NOVA",script="01_code/modules/SUA/08_nhanes_serum_uric_acid_diet_nova.R",args=character(0)),
  list(id="SUA_final_RCS",script="01_code/modules/SUA/final_aligned_rcs.R",args=character(0)),
  list(id="GA_models_assumption",script="GA/01_代码/11_ga_nhanes.R",args=character(0)),
  list(id="age_preflight",script="Diagnostics/01_code/01_preflight.R",args=c("NHANES")),
  list(id="age_models",script="Diagnostics/01_code/02_diagnostic_models.R",args=c("NHANES","formal")),
  list(id="age_pooling",script="Diagnostics/01_code/03_pool_diagnostics.R",args=character(0))
)

args <- commandArgs(TRUE)
arg <- function(key,default='') {z<-grep(paste0('^--',key,'='),args,value=TRUE);if(length(z))sub(paste0('^--',key,'='),'',z[1]) else default}
flag <- function(key)paste0('--',key) %in% args
self <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
package_root <- dirname(normalizePath(self,winslash='/',mustWork=TRUE))
work <- normalizePath(arg('work',file.path(package_root,'run')),winslash='/',mustWork=FALSE)
raw_nh <- arg('raw-nhanes',file.path(package_root,'data/nhanes'))
raw_kn <- arg('raw-knhanes',file.path(package_root,'data/knhanes'))
fped <- arg('fped',file.path(package_root,'data/fped'))
official <- arg('official-xpt',file.path(package_root,'data/official_xpt'))
extra_lib <- arg('r-library',Sys.getenv('UPF_R_LIBRARY'))
if(nzchar(extra_lib)) {
  libs <- strsplit(extra_lib,.Platform$path.sep,fixed=TRUE)[[1]]
  if(any(!dir.exists(libs)))stop('R library directory not found: ',extra_lib)
  .libPaths(unique(c(libs,.libPaths())))
  Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep),R_LIBS_USER=paste(.libPaths(),collapse=.Platform$path.sep))
}
preflight <- function() {
  if(!requireNamespace('digest',quietly=TRUE))stop('Install the R package digest first.')
  valid <- mapply(file_matches,file.path(package_root,names(classification_hashes)),classification_hashes)
  if(any(!valid))stop('Missing or changed classification file: ',paste(names(classification_hashes)[!valid],collapse=', '))
  pk <- read_config('required_packages.csv')
  pk$actual <- vapply(pk$Package,function(p)if(requireNamespace(p,quietly=TRUE))as.character(packageVersion(p)) else 'MISSING',character(1))
  missing <- pk$Package[pk$actual=='MISSING']
  if(length(missing))stop('Missing R packages before analysis: ',paste(missing,collapse=', '))
  same <- mapply(function(a,b)package_version(a)==package_version(b),pk$actual,pk$Version)
  if(any(!same))warning('R versions differ from the reference environment: ',paste(pk$Package[!same],collapse=', '))
  invisible(pk)
}

workers <- as.integer(arg('workers','1'))
stopifnot(is.finite(workers),workers>=1,workers<=14)
rewrite <- function(s) {
  base <- '/storage/home/tmu2301/nhanes分析'
  map <- c('18_formal_cholesterol_design_reanalysis_20260907'='',
           '13_UPF_definition_revision_20260905'='KN',
           '06_NHANES_UPF七结局_第一阶段探索_20260824'='NHANES',
           '05_KNHANES_外部验证_2019_2024_20260820'='KNHANES')
  for(n in names(map))for(prefix in c(base,'/project'))s<-gsub(paste0(prefix,'/',n),file.path(work,map[[n]]),s,fixed=TRUE)
  s<-gsub(base,file.path(work,'dependencies'),s,fixed=TRUE)
  s<-gsub('"KNHANES_N_CORES"','"KNHANES_N_CORES"',s,fixed=TRUE)
  s<-gsub("KNHANES_N_CORES='14'",paste0("KNHANES_N_CORES='",workers,"'"),s,fixed=TRUE)
  s
}
run_script <- function(name,argv=character()) {
  p<-file.path(work,name);stopifnot(file.exists(p))
  logfile<-file.path(work,'logs',paste0(gsub('[^A-Za-z0-9]','_',name),'_',paste(argv,collapse='_'),'.log'))
  status<-system2(file.path(R.home('bin'),'Rscript'),c(shQuote(p),vapply(argv,shQuote,character(1))),stdout=logfile,stderr=logfile)
  if(status!=0)stop('Stage failed: ',name,'; see ',logfile)
}
run_text <- function(name,code) {
  p<-paste0('helpers/',name,'.R')
  writeLines(c('options(warn=1)','work <- Sys.getenv("UPF_RUN_ROOT")',code),file.path(work,p),useBytes=TRUE)
  run_script(p)
}
materialize <- function() {
  for(n in names(upf_source)) {p<-file.path(work,n);dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);writeLines(c('options(warn=1)',rewrite(upf_source[[n]])),p,useBytes=TRUE)}
  for(p in c('logs','03_qc','04_results','KN/03_qc/candidate','NHANES/00_计划与官方依据','NHANES/00_计划文档','NHANES/01_数据配置','NHANES/03_派生数据/正式五结局_v1','NHANES/04_结果/正式五结局_v1','NHANES/05_表格/正式五结局_v1','NHANES/07_运行记录/正式五结局_v1','NHANES/08_测试/正式五结局_v1','KNHANES/02_数据配置'))dir.create(file.path(work,p),recursive=TRUE,showWarnings=FALSE)
  for(name in names(embedded_config)[startsWith(names(embedded_config),'NHANES_config/')]) {
    rel<-substring(name,nchar('NHANES_config/')+1)
    write_config(name,file.path(work,'NHANES/01_数据配置',rel))
  }
  con<-gzfile(file.path(package_root,'classifications/NHANES_NOVA.csv.gz'),'rt')
  nh<-tryCatch(read.csv(con,colClasses='character',check.names=FALSE),finally=close(con))
  stopifnot(setequal(unique(nh$cycle),c('E','F','G','H','I','P')))
  dest<-file.path(work,'NHANES/01_数据配置/nova_manifest');dir.create(dest,showWarnings=FALSE)
  for(cy in unique(nh$cycle))write.csv(nh[nh$cycle==cy,c('code','desc','nova','rule')],
    file.path(dest,paste0('manifest_',cy,'.csv')),row.names=FALSE,na='NA',fileEncoding='UTF-8')
  write_config('KNHANES_outcome_definitions.csv',file.path(work,'KNHANES/02_数据配置/outcome_definitions_v1.csv'))
  write_config('Table1_layout.csv',file.path(work,'KN/source_snapshot/Table1_baseline_six_outcomes.csv'))
  mapgz<-file.path(package_root,'classifications/KNHANES_NOVA.csv.gz')
  target<-file.path(work,'KNHANES/02_数据配置/nova_context_adjudicated.csv')
  {input<-gzfile(mapgz,'rb');output<-file(target,'wb');repeat {b<-readBin(input,'raw',1024*1024);if(!length(b))break;writeBin(b,output)};close(input);close(output)}
  for(p in c('NHANES/00_计划与官方依据','NHANES/00_计划文档'))writeLines('Explicit execution of the public reproduction pipeline.',file.path(work,p,'FORMAL_ANALYSIS_AUTHORIZED.txt'))
  Sys.setenv(UPF_RUN_ROOT=work,NHANES_PHASE1_ROOT=file.path(work,'NHANES'),NHANES_REV22_ROOT=file.path(work,'NHANES'),NHANES_ROOT=file.path(work,'NHANES'),KNHANES_ROOT=file.path(work,'KNHANES'),NHANES_PHASE1_FORMAL_AUTH='I_AUTHORIZE_NHANES_PHASE1_FORMAL_RUN',NHANES_REV22_FORMAL_AUTH='I_HAVE_EXPLICIT_USER_AUTHORIZATION',NHANES_PHASE1_WORKERS=workers,KNHANES_N_CORES=workers,KNHANES_MI_M='50',KNHANES_MI_MAXIT='20',NHANES_SUBGROUP_RUN_ID='FORMAL_CHOL_DESIGN_20260907',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
  if(nzchar(raw_nh))Sys.setenv(NHANES_RAW_ROOT=raw_nh)
  if(nzchar(raw_kn))Sys.setenv(KNHANES_EXTRACT_ROOT=raw_kn)
  if(nzchar(fped))Sys.setenv(NHANES_FPED_DIR=fped)
  if(nzchar(official))Sys.setenv(NHANES_OFFICIAL_XPT_DIR=official)
}

allowed <- c('packages','check','download','check-data','offline','all','prepare-nhanes','prepare-knhanes','models','outputs',
             'work','raw-nhanes','raw-knhanes','fped','official-xpt','workers','r-library')
unknown <- setdiff(sub('=.*','',sub('^--','',args)),allowed)
if(length(unknown))stop('Unknown option: ',paste(unknown,collapse=', '))
if(flag('packages')) {print(read_config('required_packages.csv'),row.names=FALSE);quit(status=0)}
if(flag('check') || !length(args)) {
  bad<-vapply(upf_source,function(s)inherits(try(parse(text=s),silent=TRUE),'try-error'),logical(1))
  if(any(bad))stop('Invalid embedded R source: ',paste(names(upf_source)[bad],collapse=', '))
  cat(length(upf_source),'embedded R modules parsed successfully. No analysis run.\n')
  if(!length(args))cat('Use --all to run; --download to obtain NHANES/FPED; --check-data to verify inputs. See README.md.\n')
  quit(status=0)
}
data_manifest <- function() {
  x <- read_config('required_raw_files.csv',stringsAsFactors=FALSE)
  roots <- c(nhanes=raw_nh,official_xpt=official,knhanes=raw_kn)
  x$path <- file.path(unname(roots[x$input_root]),x$relative_path)
  stem <- toupper(tools::file_path_sans_ext(basename(x$relative_path)))
  year <- ifelse(grepl('^P_',stem),'2017',substr(x$relative_path,1,4))
  x$url <- ifelse(x$input_root=='knhanes','',paste0('https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/',year,'/DataFiles/',stem,'.xpt'))
  x$archive <- FALSE
  f <- read_config('NHANES_config/diet_quality_manifest.csv',stringsAsFactors=FALSE)
  f <- f[f$record_type=='data_file',]
  fn <- basename(f$relative_path)
  links <- read_config('download_sources.csv',stringsAsFactors=FALSE)
  urls <- links$url[match(fn,links$file)]
  stopifnot(!anyNA(urls))
  rbind(x,data.frame(input_root='fped',relative_path=fn,bytes=NA_real_,sha256=f$sha256,
                    path=file.path(fped,fn),url=urls,archive=TRUE))
}

file_matches <- function(path,sha) {
  file.exists(path) && identical(tolower(digest::digest(file=path,algo='sha256')),tolower(sha))
}

fetch_file <- function(url,dest,sha,archive=FALSE) {
  if(file.exists(dest)) {
    if(!file_matches(dest,sha))stop('Data version differs from the manuscript: ',dest,
      '\nKeep the existing file and check the survey release; it has not been overwritten.')
    return(invisible(dest))
  }
  dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
  staging <- tempfile('download_',tmpdir=dirname(dest))
  dir.create(staging)
  on.exit(unlink(staging,recursive=TRUE),add=TRUE)
  downloaded <- file.path(staging,if(archive)'source.exe' else basename(dest))
  old <- options(timeout=max(600,getOption('timeout')))
  on.exit(options(old),add=TRUE)
  message('Downloading ',basename(dest))
  status <- NULL
  for(attempt in 1:3) {
    status <- tryCatch(suppressWarnings(utils::download.file(url,downloaded,mode='wb',method='libcurl',quiet=TRUE)),
                       error=function(e)conditionMessage(e))
    if(is.numeric(status) && status==0L)break
  }
  if(!is.numeric(status) || status!=0L)stop('Download failed after three attempts: ',url,'\n',status)
  candidate <- downloaded
  if(archive) {
    candidate <- file.path(staging,basename(dest))
    seven <- Sys.which(c('7zz','7z')); seven <- seven[nzchar(seven)]
    if(length(seven)) {
      status <- system2(seven[1],c('e','-y',shQuote(paste0('-o',staging)),
                        shQuote(downloaded),shQuote(basename(dest))),stdout=FALSE,stderr=FALSE)
    } else {
      tar <- Sys.which('bsdtar')
      if(!nzchar(tar) && .Platform$OS.type=='windows')tar<-Sys.which('tar')
      if(!nzchar(tar))stop('Install 7zip (7zz/7z) or bsdtar to extract USDA FPED files.')
      status <- system2(tar,c('-xf',shQuote(downloaded),'-C',shQuote(staging),
                             shQuote(basename(dest))),stdout=FALSE,stderr=FALSE)
    }
    if(status!=0L || !file.exists(candidate))stop('Could not extract FPED file: ',basename(dest))
  }
  if(!file_matches(candidate,sha))stop('Downloaded file does not match the manuscript data version: ',url)
  if(!file.rename(candidate,dest))stop('Could not save downloaded data: ',dest)
  invisible(dest)
}

prepare_data <- function(groups,download=TRUE) {
  if(!requireNamespace('digest',quietly=TRUE))stop('Install the R package digest first.')
  x <- data_manifest(); x <- x[x$input_root %in% groups,]
  missing_kn <- x$input_root=='knhanes' & !file.exists(x$path)
  if(any(missing_kn))stop('Obtain KNHANES 2019-2024 ALL and 24RC SAS files from https://knhanes.kdca.go.kr/\n',
    'Place the extracted files in ',raw_kn,' (or use --raw-knhanes=DIR).\nMissing: ',
    paste(basename(x$path[missing_kn]),collapse=', '))
  for(i in seq_len(nrow(x))) {
    if(!file.exists(x$path[i]) && (!download || !nzchar(x$url[i])))
      stop('Missing data file: ',x$path[i])
    fetch_file(x$url[i],x$path[i],x$sha256[i],x$archive[i])
  }
  message(nrow(x),' data files verified.')
  invisible(x)
}

if(flag('download'))prepare_data(c('nhanes','official_xpt','fped'))
if(flag('check-data'))prepare_data(c('nhanes','official_xpt','fped','knhanes'),download=FALSE)
if(!any(vapply(c('all','prepare-nhanes','prepare-knhanes','models','outputs'),flag,logical(1))))quit(status=0)
if(.Platform$OS.type!='unix')stop('Full analysis requires Linux; download and input checks also work on Windows.')
if(flag('all') && dir.exists(work) && length(list.files(work,all.files=TRUE,no..=TRUE)))
  stop('Use an empty output directory for --all; specify --work=NEW_DIRECTORY.')
preflight()
if(flag('all') || flag('prepare-knhanes') || flag('models'))prepare_data('knhanes',download=FALSE)
if(flag('all') || flag('prepare-nhanes') || flag('models'))prepare_data(c('nhanes','official_xpt','fped'),download=!flag('offline'))
materialize()
if(flag('prepare-nhanes') || flag('all')) {
  stopifnot(dir.exists(raw_nh),dir.exists(fped))
  run_text('prepare_nhanes','source(file.path(work,"helpers/prepare.R")); prepare_nhanes()')
}
if(flag('prepare-knhanes') || flag('all')) {
  stopifnot(dir.exists(raw_kn))
  run_text('prepare_knhanes','source(file.path(work,"helpers/prepare.R")); prepare_knhanes()')
}
if(flag('models') || flag('all')) {
  manifest<-jsonlite::read_json(file.path(work,'03_qc/mi_job_manifest.json'),simplifyVector=TRUE)
  for(id in manifest$id) {
    for(chain in 1:50)run_script('01_code/05_formal_mi_worker.R',c(id,as.character(chain),'20'))
    run_script('01_code/09_pool_checkpoint.R',c(id,'20','50'))
  }
  run_text('retain_completed_models',c(
    'manifest<-jsonlite::read_json(file.path(work,"03_qc/mi_job_manifest.json"),simplifyVector=TRUE)',
    'main<-day1<-list()',
    'for(i in seq_len(nrow(manifest))) {s<-manifest[i,];src<-file.path(work,"02_mi/combined",paste0(s$id,"_I020.rds"));obj<-readRDS(src);stopifnot(obj$m==50,obj$mids$iteration==20);dir.create(dirname(s$output),recursive=TRUE,showWarnings=FALSE);file.copy(src,s$output,overwrite=TRUE);z<-read.csv(file.path(work,"04_results/checkpoints",paste0(s$id,"_I020.csv")));stopifnot(all(z$successful_imputations==50),all(is.finite(z$estimate)),all(z$ci_low<=z$estimate & z$estimate<=z$ci_high));if(s$branch=="main")main[[s$id]]<-z else day1[[s$id]]<-z}',
    'main<-do.call(rbind,main);day1<-do.call(rbind,day1);stopifnot(nrow(main)==24,nrow(day1)==5)',
    'write.csv(main,file.path(work,"04_results/NHANES_main_M1_M4_iter20.csv"),row.names=FALSE)',
    'five<-main[main$outcome!="serum_uric_acid",];p<-file.path(work,"NHANES/04_结果/正式五结局_v1");write.csv(five,file.path(p,"phase1_models_all_outcomes.csv"),row.names=FALSE);for(oc in unique(five$outcome))write.csv(five[five$outcome==oc,],file.path(p,paste0("phase1_models_",oc,".csv")),row.names=FALSE)',
    'dir.create(file.path(work,"SUA/02_主分析/NHANES"),recursive=TRUE,showWarnings=FALSE);write.csv(main[main$outcome=="serum_uric_acid",],file.path(work,"SUA/02_主分析/NHANES/NHANES_serum_uric_acid_M1_M4.csv"),row.names=FALSE)',
    'dir.create(file.path(work,"Day1/04_results"),recursive=TRUE,showWarnings=FALSE);write.csv(day1,file.path(work,"Day1/04_results/day1_M3.csv"),row.names=FALSE)'))
  run_text('knhanes_imputation',c(
    'source(file.path(work,"KN/source_snapshot/05_imputation.R"))',
    'fp<-file.path(work,"KN/02_derived/analysis_frame_fasting_corrected.rds");x<-readRDS(file.path(work,"KN/02_derived/mice_input_and_config.rds"));cfg<-x$config',
    'imp<-run_auditable_parallel_mice(x$data,50,20,cfg$method,cfg$predictor,cfg$where,cfg$ignore,n.core=14,parallelseed=KN_SEED)',
    'attr(imp,"knhanes_frame_sha256")<-sha256_file(fp);saveRDS(imp,file.path(work,"KN/02_derived/candidate_mids.rds"));writeLines(capture.output(sessionInfo()),file.path(work,"03_qc/sessionInfo.txt"))',
    'write.csv(chain_diagnostics(imp),file.path(work,"KN/03_qc/candidate/chain_diagnostics.csv"),row.names=FALSE);writeLines("Completed 50 imputations at 20 iterations",file.path(work,"KN/03_qc/candidate/MICE_COMPLETED.txt"))'))
  run_script('KN/01_code/04_run_models.R','candidate')
  for(j in analysis_jobs)run_script(j$script,j$args)
  run_script('KN/01_code/08_summarize_models.R')
  run_script('helpers/table1.R')
  run_script('helpers/extend_rcs.R')
  run_script('01_code/18_table1_cholesterol_summary.R')
  run_script('01_code/19_table1_F_and_S2_missingness.R')
}
if(flag('outputs') || flag('all')) {
  source(file.path(work,'helpers/outputs.R'),local=FALSE)
  generate_outputs(package_root,work)
}
cat('Requested stages finished. Logs and regenerated outputs: ',work,'\n',sep='')
