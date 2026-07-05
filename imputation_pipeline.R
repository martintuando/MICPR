# ==============================================================================
# MICPR: Hybrid Imputation Pipeline for Censored Mortality Data (1-9)
# Based on the methodological framework of Erdman et al. (2021)
# ==============================================================================
library(mice)
library(dplyr)
library(tidyr)

set.seed(42)

# --- 1. Path Configuration ---
INPUT_DIR  <- "C:/Users/Martin/Desktop/Github/MICPR/CDC_data"
OUTPUT_DIR <- "C:/Users/Martin/Desktop/Github/MICPR/Result"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

years_vector <- 2014:2019


# ==============================================================================
# PART A: MENTAL & BEHAVIOURAL DISORDERS
# ==============================================================================
cat("\n=== PROCESSING: Mental & Behavioural Disorders ===\n")

# 1. Load and clean files
mental_list <- list()
for (i in seq_along(years_vector)) {
  file_name <- paste0("Mentalandbehaviouraldisorders3Y", years_vector[i], ".csv")
  file_path <- file.path(INPUT_DIR, file_name)
  
  cat("Loading:", file_name, "\n")
  df <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Force all columns to numeric (converts suppressed text strings to NA)
  df[] <- lapply(df, function(x) as.numeric(as.character(x)))
  df <- df[!is.na(df$`County Code`), ]
  df$Year <- years_vector[i]
  
  mental_list[[i]] <- df
}
mental_long <- bind_rows(mental_list) %>% arrange(`County Code`, Year)

# 2. Rule-Based Substitution (Stage 1)
na_idx_m <- which(is.na(mental_long$Deaths))
if (length(na_idx_m) > 0) {
  for (idx in na_idx_m) {
    c_county <- mental_long$`County Code`[idx]
    c_year   <- mental_long$Year[idx]
    
    prev_y <- mental_long %>% filter(`County Code` == c_county, Year == c_year - 1) %>% pull(Deaths)
    next_y <- mental_long %>% filter(`County Code` == c_county, Year == c_year + 1) %>% pull(Deaths)
    
    prev_y <- ifelse(length(prev_y) == 0, NA, prev_y)
    next_y <- ifelse(length(next_y) == 0, NA, next_y)
    
    if (!is.na(prev_y) && prev_y >= 10 && !is.na(next_y) && next_y >= 10) {
      mental_long$Deaths[idx] <- 10
    } else if ((!is.na(prev_y) && prev_y >= 10 && (is.na(next_y) || next_y == 0)) ||
               (!is.na(next_y) && next_y >= 10 && (is.na(prev_y) || prev_y == 0))) {
      avail_val <- ifelse(!is.na(prev_y) && prev_y >= 10, prev_y, next_y)
      mental_long$Deaths[idx] <- round(avail_val / 2)
    }
  }
}

# 3. MICE on 6Y Aggregate (Stage 2)
cat("Processing 6Y Aggregate for Mental Health...\n")
m6y_path <- file.path(INPUT_DIR, "Mentalandbehaviouraldisorders6Y.csv")
mental_6y <- read.csv(m6y_path, stringsAsFactors = FALSE, check.names = FALSE)

mental_6y$Deaths      <- as.numeric(as.character(mental_6y$Deaths))
mental_6y$Population  <- as.numeric(as.character(mental_6y$Population))
mental_6y$`County Code` <- as.numeric(as.character(mental_6y$`County Code`))
mental_6y             <- mental_6y[!is.na(mental_6y$`County Code`), ]

imp_data_m6y <- mental_6y %>% select(Deaths, Population, `County Code`)
na_idx_m6y   <- which(is.na(imp_data_m6y$Deaths))

imp_data_m6y$Deaths_log <- log(imp_data_m6y$Deaths + 1)
imp_data_m6y$Deaths_log[is.na(imp_data_m6y$Deaths)] <- NA

meth_m6y <- make.method(imp_data_m6y)
meth_m6y["Deaths_log"] <- "pmm"
pred_m6y <- make.predictorMatrix(imp_data_m6y)
pred_m6y["Deaths_log", "Population"]  <- 1
pred_m6y["Deaths_log", "County Code"] <- 0
pred_m6y["Deaths_log", "Deaths"]      <- 0

mice_m6y <- mice(imp_data_m6y, m = 5, maxit = 5, method = meth_m6y, predictorMatrix = pred_m6y, printFlag = FALSE)

for (i in 1:mice_m6y$m) {
  raw_vals <- mice_m6y$imp$Deaths_log[, i]
  back_trans <- exp(raw_vals) - 1
  back_trans[back_trans < 0] <- 0
  finite_v <- back_trans[is.finite(back_trans)]
  
  if (length(finite_v) == 0) {
    scaled <- rep(1, length(raw_vals))
  } else {
    min_v <- min(finite_v); max_v <- max(finite_v)
    scaled <- if (max_v == min_v) rep(1, length(raw_vals)) else 1 + (back_trans - min_v) * (9 - 1) / (max_v - min_v)
  }
  mice_m6y$imp$Deaths_log[, i] <- round(pmax(1, pmin(9, scaled)))
}

comp_m6y <- complete(mice_m6y, 1)
mental_6y_imp <- mental_6y
if (length(na_idx_m6y) > 0) {
  mental_6y_imp$Deaths[na_idx_m6y] <- comp_m6y$Deaths_log[na_idx_m6y]
}

# 4. Temporal Decomposition (Stage 3)
mental_6y_imp_clean <- mental_6y_imp %>% mutate(`County Code` = as.numeric(as.character(`County Code`)))
mental_long <- mental_long %>%
  mutate(`County Code` = as.numeric(as.character(`County Code`))) %>%
  left_join(mental_6y_imp_clean %>% select(`County Code`, Deaths_6yr_imputed = Deaths), by = "County Code")

na_idx_m_final <- which(is.na(mental_long$Deaths))
if (length(na_idx_m_final) > 0) {
  vals_6yr <- mental_long$Deaths_6yr_imputed[na_idx_m_final]
  raw_3yr  <- vals_6yr / 2
  raw_3yr[raw_3yr < 0] <- 0
  mental_long$Deaths[na_idx_m_final] <- round(pmax(1, raw_3yr))
}


# ==============================================================================
# PART B: SUICIDES
# ==============================================================================
cat("\n=== PROCESSING: Suicides ===\n")

# 1. Load and clean annual files
suicide_list <- list()
for (i in seq_along(years_vector)) {
  file_name <- paste0("Suicide3Y", years_vector[i], ".csv")
  file_path <- file.path(INPUT_DIR, file_name)
  
  cat("Loading:", file_name, "\n")
  df <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Force all columns to numeric (converts suppressed text strings to NA)
  df[] <- lapply(df, function(x) as.numeric(as.character(x)))
  df <- df[!is.na(df$`County Code`), ]
  df$Year <- years_vector[i]
  
  suicide_list[[i]] <- df
}
suicide_long <- bind_rows(suicide_list) %>% arrange(`County Code`, Year)

# 2. Rule-Based Substitution (Stage 1)
na_idx_s <- which(is.na(suicide_long$Deaths))
if (length(na_idx_s) > 0) {
  for (idx in na_idx_s) {
    c_county <- suicide_long$`County Code`[idx]
    c_year   <- suicide_long$Year[idx]
    
    prev_y <- suicide_long %>% filter(`County Code` == c_county, Year == c_year - 1) %>% pull(Deaths)
    next_y <- suicide_long %>% filter(`County Code` == c_county, Year == c_year + 1) %>% pull(Deaths)
    
    prev_y <- ifelse(length(prev_y) == 0, NA, prev_y)
    next_y <- ifelse(length(next_y) == 0, NA, next_y)
    
    if (!is.na(prev_y) && prev_y >= 10 && !is.na(next_y) && next_y >= 10) {
      suicide_long$Deaths[idx] <- 10
    } else if ((!is.na(prev_y) && prev_y >= 10 && (is.na(next_y) || next_y == 0)) ||
               (!is.na(next_y) && next_y >= 10 && (is.na(prev_y) || prev_y == 0))) {
      avail_val <- ifelse(!is.na(prev_y) && prev_y >= 10, prev_y, next_y)
      suicide_long$Deaths[idx] <- round(avail_val / 2)
    }
  }
}

# 3. MICE on 6Y Aggregate (Stage 2)
cat("Processing 6Y Aggregate for Suicides...\n")
s6y_path <- file.path(INPUT_DIR, "Suicides6Y.csv")
suicide_6y <- read.csv(s6y_path, stringsAsFactors = FALSE, check.names = FALSE)

suicide_6y$Deaths      <- as.numeric(as.character(suicide_6y$Deaths))
suicide_6y$Population  <- as.numeric(as.character(suicide_6y$Population))
suicide_6y$`County Code` <- as.numeric(as.character(suicide_6y$`County Code`))
suicide_6y             <- suicide_6y[!is.na(suicide_6y$`County Code`), ]

imp_data_s6y <- suicide_6y %>% select(Deaths, Population, `County Code`)
na_idx_s6y   <- which(is.na(imp_data_s6y$Deaths))

imp_data_s6y$Deaths_log <- log(imp_data_s6y$Deaths + 1)
imp_data_s6y$Deaths_log[is.na(imp_data_s6y$Deaths)] <- NA

meth_s6y <- make.method(imp_data_s6y)
meth_s6y["Deaths_log"] <- "pmm"
pred_s6y <- make.predictorMatrix(imp_data_s6y)
pred_s6y["Deaths_log", "Population"]  <- 1
pred_s6y["Deaths_log", "County Code"] <- 0
pred_s6y["Deaths_log", "Deaths"]      <- 0

mice_s6y <- mice(imp_data_s6y, m = 5, maxit = 5, method = meth_s6y, predictorMatrix = pred_s6y, printFlag = FALSE)

for (i in 1:mice_s6y$m) {
  raw_vals <- mice_s6y$imp$Deaths_log[, i]
  back_trans <- exp(raw_vals) - 1
  back_trans[back_trans < 0] <- 0
  finite_v <- back_trans[is.finite(back_trans)]
  
  if (length(finite_v) == 0) {
    scaled <- rep(1, length(raw_vals))
  } else {
    min_v <- min(finite_v); max_v <- max(finite_v)
    scaled <- if (max_v == min_v) rep(1, length(raw_vals)) else 1 + (back_trans - min_v) * (9 - 1) / (max_v - min_v)
  }
  mice_s6y$imp$Deaths_log[, i] <- round(pmax(1, pmin(9, scaled)))
}

comp_s6y <- complete(mice_s6y, 1)
suicide_6y_imp <- suicide_6y
if (length(na_idx_s6y) > 0) {
  suicide_6y_imp$Deaths[na_idx_s6y] <- comp_s6y$Deaths_log[na_idx_s6y]
}

# 4. Temporal Decomposition (Stage 3)
suicide_6y_imp_clean <- suicide_6y_imp %>% mutate(`County Code` = as.numeric(as.character(`County Code`)))
suicide_long <- suicide_long %>%
  mutate(`County Code` = as.numeric(as.character(`County Code`))) %>%
  left_join(suicide_6y_imp_clean %>% select(`County Code`, Deaths_6yr_imputed = Deaths), by = "County Code")

na_idx_s_final <- which(is.na(suicide_long$Deaths))
if (length(na_idx_s_final) > 0) {
  vals_6yr <- suicide_long$Deaths_6yr_imputed[na_idx_s_final]
  raw_3yr  <- vals_6yr / 2
  raw_3yr[raw_3yr < 0] <- 0
  suicide_long$Deaths[na_idx_s_final] <- round(pmax(1, raw_3yr))
}


# ==============================================================================
# PART C: MASTER COMBINATION & WIDE EXPORT
# ==============================================================================
cat("\n=== MERGING AND EXPORTING FINAL ANALYSIS DATASET ===\n")

mental_prep  <- mental_long %>% select(`County Code`, Year, Deaths_Mental = Deaths, Population)
suicide_prep <- suicide_long %>% select(`County Code`, Year, Deaths_Suicide = Deaths)
master_long <- full_join(mental_prep, suicide_prep, by = c("County Code", "Year"))

# Pivot to Wide Format
master_wide <- master_long %>%
  pivot_wider(
    id_cols     = `County Code`,
    names_from  = Year,
    values_from = c(Deaths_Mental, Deaths_Suicide, Population),
    names_sep   = "_"
  )

output_file <- file.path(OUTPUT_DIR, "Master_Imputed_Outcomes_Wide.csv")
write.csv(master_wide, file = output_file, row.names = FALSE)

cat("\n==============================================================================")
cat("\nExecution completed successfully!")
cat("\nSaved Master Wide File to:", output_file)
cat("\nFinal dataset dimensions:", dim(master_wide))
cat("\n==============================================================================\n")
