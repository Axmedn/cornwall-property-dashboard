# ============================================================================
# CORNWALL PROPERTY PRICE ANALYSIS - ULTRA-COMPACT VERSION
# ============================================================================
# 5 Models: GLM (additive, interaction, +year), Spline GLM, Mixed Model
# ============================================================================

library(tidyverse)
library(broom)
library(splines)
library(lme4)
library(patchwork)
# === DATA PREP ===
data_model <- read_csv("data/cornwall_property.csv.gz", show_col_types = FALSE) %>%
  mutate(
    house_size = as.numeric(tfarea),
    rooms = as.numeric(numberrooms),
    energy_rating = as.numeric(CURRENT_ENERGY_EFFICIENCY),
    env_impact = as.numeric(ENVIRONMENT_IMPACT_CURRENT),
    property_type = factor(case_when(
      propertytype == "D" ~ "Detached", propertytype == "S" ~ "Semi-Detached",
      propertytype == "T" ~ "Terraced", propertytype == "F" ~ "Flat",
      TRUE ~ "Other"), levels = c("Detached", "Semi-Detached", "Terraced", "Flat", "Other")),
    build_type = factor(if_else(oldnew == "Y", "New", "Old"), levels = c("Old", "New")),
    towncity = factor(str_to_title(str_trim(towncity))),
    postcode_area = factor(str_sub(str_squish(toupper(postcode)), 1, 3)),
    year = as.numeric(year),
    log_price = log(price)
  ) %>%
  filter(price > 0, house_size > 0, rooms > 0, !is.na(towncity),
         !is.na(year), !is.na(energy_rating), !is.na(env_impact)) %>%
  select(price, log_price, house_size, rooms, energy_rating, env_impact,
         property_type, build_type, towncity, postcode_area, year) %>%
  drop_na()

cat("Loaded:", nrow(data_model), "properties,", min(data_model$year), "-",
    max(data_model$year), "\n\n")

# Variable selection: Keep house_size (r=0.69), energy (0.02); Exclude rooms (r=0.83), env (0.04)
cor_matrix <- cor(data_model %>% select(house_size, rooms, energy_rating, env_impact, log_price))
cat("Correlations with log_price:", paste(round(cor_matrix[,"log_price"], 3), collapse=", "), "\n")
cat("Include: house_size, energy_rating; Exclude: rooms, env_impact\n\n")

# === MODEL FITTING ===
models <- list(
  M1_additive = glm(log_price ~ house_size + energy_rating + property_type +
                      build_type + towncity + year, data = data_model),
  M2_size_type = glm(log_price ~ house_size * property_type + energy_rating +
                       build_type + towncity + year, data = data_model),
  M3_size_year_type = glm(log_price ~ house_size * property_type + energy_rating +
                            build_type + towncity + year * property_type, data = data_model),
  M4_spline_size = glm(log_price ~ ns(house_size, df = 4) * property_type + energy_rating +
                         build_type + towncity + year, data = data_model),
  M5_mixed_town = lmer(log_price ~ house_size * property_type + energy_rating +
                         build_type + year + (1 | towncity),
                       data = data_model, REML = FALSE)
)

# --- AIC + pseudo-R² (GLMs only) -----------------------------------------
model_comparison <- tibble(
  Model = c("M1: Additive", "M2: Size×Type", "M3: Size×Type+Year×Type",
            "M4: Spline(4df)×Type", "M5: Mixed(town RE)"),
  AIC = map_dbl(models, AIC),
  R_squared = c(map_dbl(models[1:4], ~1 - .x$deviance/.x$null.deviance), NA)
)

print(model_comparison)

# --- NEW: test-set RMSE on log-price for all models ----------------------
set.seed(123)
train_idx <- sample(seq_len(nrow(data_model)), 0.8 * nrow(data_model))
train <- data_model[train_idx, ]
test  <- data_model[-train_idx, ]

rmse_log <- function(model, test_data) {
  if (inherits(model, "lmerMod")) {
    preds <- predict(model, newdata = test_data, re.form = NA)
  } else {
    preds <- predict(model, newdata = test_data)
  }
  sqrt(mean((test_data$log_price - preds)^2))
}

rmse_comparison <- tibble(
  Model = names(models),
  Class = map_chr(models, ~class(.x)[1]),
  AIC   = map_dbl(models, AIC),
  RMSE_log = map_dbl(models, ~rmse_log(.x, test))
)

print(rmse_comparison)
cat("\n")

# --- choose final model ---------------------------------------------------
final_model <- models$M2_size_type  # Change to M4_spline_size or M5_mixed_town as needed
cat("Selected:", names(models)[which(map_lgl(models, identical, final_model))], "\n\n")

# Key effects
coefs <- tidy(final_model)
energy_pct <- 100 * (exp(coefs$estimate[coefs$term == "energy_rating"]) - 1)
build_pct <- if("build_typeNew" %in% coefs$term)
  100 * (exp(coefs$estimate[coefs$term == "build_typeNew"]) - 1) else NA

cat("Effects: Energy", sprintf("%+.3f%%", energy_pct))
if(!is.na(build_pct)) cat(", New build", sprintf("%+.1f%%", build_pct))

if(inherits(final_model, "glm")) {
  cat(sprintf("\nR²=%.4f, AIC=%.1f, N=%s\n\n",
              1 - final_model$deviance/final_model$null.deviance, AIC(final_model),
              format(nobs(final_model), big.mark=",")))
} else if(inherits(final_model, "lmerMod")) {
  cat(sprintf("\nAIC=%.1f, Towns=%d, N=%s\n\n", AIC(final_model),
              ngrps(final_model), format(nobs(final_model), big.mark=",")))
}

# === DIAGNOSTICS ===
if(inherits(final_model, "glm")) {
  diag_data <- tibble(fitted = fitted(final_model), std_resid = rstandard(final_model))
  p1 <- ggplot(diag_data, aes(fitted, std_resid)) +
    geom_point(alpha=0.2) + geom_hline(yintercept=0, linetype=2, color="red") +
    geom_smooth(se=F, color="blue") + labs(title="Residuals vs Fitted") + theme_minimal()
  p2 <- ggplot(diag_data, aes(sample=std_resid)) +
    stat_qq(alpha=0.3) + stat_qq_line(color="red") + labs(title="Q-Q Plot") + theme_minimal()
  ggsave("diagnostics.png", p1 + p2, width=10, height=4, dpi=300)
} else if(inherits(final_model, "lmerMod")) {
  diag_data <- tibble(fitted = fitted(final_model), resid = residuals(final_model))
  p1 <- ggplot(diag_data, aes(fitted, resid)) +
    geom_point(alpha=0.2) + geom_hline(yintercept=0, linetype=2, color="red") +
    geom_smooth(se=F, color="blue") + labs(title="LMM: Residuals") + theme_minimal()
  p2 <- ggplot(diag_data, aes(sample=resid)) +
    stat_qq(alpha=0.3) + stat_qq_line(color="red") + labs(title="LMM: Q-Q") + theme_minimal()
  ggsave("diagnostics_lmm.png", p1 + p2, width=10, height=4, dpi=300)
  
  # Random effects
  re_data <- ranef(final_model)$towncity %>%
    as_tibble(rownames="town") %>% rename(effect=`(Intercept)`) %>%
    arrange(effect) %>% slice_head(n=20)
  p_re <- ggplot(re_data, aes(reorder(town, effect), effect)) +
    geom_point(size=3, color="steelblue") + geom_hline(yintercept=0, linetype=2, color="red") +
    coord_flip() + labs(title="Random Effects: Top 20", x="", y="Effect") + theme_minimal()
  ggsave("random_effects.png", p_re, width=8, height=6, dpi=300)
}

# Spline visualization
if(identical(final_model, models$M4_spline_size)) {
  spline_data <- expand_grid(
    house_size = seq(min(data_model$house_size), max(data_model$house_size), length.out=100),
    property_type = "Detached", energy_rating = median(data_model$energy_rating),
    build_type = "Old", towncity = levels(data_model$towncity)[1], year = 2019
  ) %>% mutate(predicted_price = exp(predict(final_model, newdata=.)))
  p_spline <- ggplot(spline_data, aes(house_size, predicted_price)) +
    geom_line(color="steelblue", size=1.2) +
    labs(title="Spline Effect", x="Size (m²)", y="Price (£)") +
    scale_y_continuous(labels=scales::comma) + theme_minimal()
  ggsave("spline_effect.png", p_spline, width=8, height=5, dpi=300)
}

cat("Saved: diagnostics plot(s)\n\n")

# === PREDICTION FUNCTIONS ===

predict_value <- function(house_size, energy_rating, property_type, build_type,
                          towncity, year, conf=0.95, model=final_model) {
  new_data <- tibble(
    house_size, energy_rating, year,
    property_type = factor(property_type, levels=levels(data_model$property_type)),
    build_type = factor(build_type, levels=levels(data_model$build_type)),
    towncity = factor(towncity, levels=levels(data_model$towncity))
  )
  if(inherits(model, "glm")) {
    pred <- predict(model, newdata=new_data, se.fit=TRUE)
    z <- qnorm((1 + conf) / 2)
    list(estimate = exp(pred$fit), lower = exp(pred$fit - z*pred$se.fit),
         upper = exp(pred$fit + z*pred$se.fit))
  } else if(inherits(model, "lmerMod")) {
    list(estimate = exp(predict(model, newdata=new_data, re.form=NA)), lower=NA, upper=NA)
  }
}

predict_trend <- function(house_size, energy_rating, property_type, build_type,
                          towncity, years, model=final_model) {
  trend_data <- tibble(
    house_size, energy_rating,
    property_type = factor(property_type, levels=levels(data_model$property_type)),
    build_type = factor(build_type, levels=levels(data_model$build_type)),
    towncity = factor(towncity, levels=levels(data_model$towncity)),
    year = years
  )
  if(inherits(model, "glm")) {
    preds <- predict(model, newdata=trend_data, se.fit=TRUE)
    trend_data %>% mutate(
      predicted_price = exp(preds$fit),
      lower = exp(preds$fit - 1.96*preds$se.fit),
      upper = exp(preds$fit + 1.96*preds$se.fit)
    )
  } else if(inherits(model, "lmerMod")) {
    trend_data %>% mutate(
      predicted_price = exp(predict(model, newdata=trend_data, re.form=NA)),
      lower = NA, upper = NA
    )
  }
}

predict_improvement <- function(current_size, current_energy, size_add, energy_add,
                                property_type, build_type, towncity, year, model=final_model) {
  current <- predict_value(current_size, current_energy, property_type,
                           build_type, towncity, year, model=model)$estimate
  improved <- predict_value(current_size+size_add, current_energy+energy_add,
                            property_type, build_type, towncity, year, model=model)$estimate
  list(current=current, improved=improved, increase=improved-current,
       pct_increase=100*(improved-current)/current)
}

cat("Functions: predict_value(), predict_trend(), predict_improvement()\n\n")

# === VISUALIZATION ===

plot_trend <- function(trend_df, title="Price Trend") {
  p <- ggplot(trend_df, aes(year, predicted_price)) +
    geom_line(color="steelblue", size=1.2) +
    labs(title=title, x="Year", y="Price (£)") +
    scale_y_continuous(labels=scales::comma) + theme_minimal()
  if(!all(is.na(trend_df$lower)))
    p <- p + geom_ribbon(aes(ymin=lower, ymax=upper), fill="steelblue", alpha=0.2)
  p
}

plot_by_type <- function(data) {
  ggplot(data, aes(property_type, log_price, fill=property_type)) +
    geom_boxplot(alpha=0.7) + labs(title="Price by Type", x="", y="Log Price") +
    theme_minimal() + theme(legend.position="none")
}

plot_locations <- function(data, top_n=15) {
  data %>% group_by(towncity) %>%
    summarise(median_price=median(price), n=n(), .groups="drop") %>%
    filter(n>=10) %>% arrange(desc(median_price)) %>% slice_head(n=top_n) %>%
    ggplot(aes(reorder(towncity, median_price), median_price)) +
    geom_col(fill="steelblue", alpha=0.7) + coord_flip() +
    labs(title=paste("Top", top_n, "Towns"), x="", y="Median Price (£)") +
    scale_y_continuous(labels=scales::comma) + theme_minimal()
}

# === EXTENDED PREDICTIONS (2020-2026) ===

rep_props <- data_model %>% filter(year==max(year)) %>%
  group_by(property_type, towncity) %>%
  summarise(med_size=median(house_size), med_energy=median(energy_rating),
            n=n(), .groups="drop") %>%
  filter(n>=10) %>% slice_max(n, n=5, by=property_type)

extended_preds <- expand_grid(year=2020:2026, rep_props) %>%
  mutate(house_size=med_size, energy_rating=med_energy,
         build_type=factor("Old", levels=levels(data_model$build_type))) %>%
  select(-med_size, -med_energy, -n)

if(inherits(final_model, "glm")) {
  preds <- predict(final_model, newdata=extended_preds, se.fit=TRUE)
  extended_preds <- extended_preds %>% mutate(price_pred=exp(preds$fit))
} else if(inherits(final_model, "lmerMod")) {
  extended_preds <- extended_preds %>%
    mutate(price_pred=exp(predict(final_model, newdata=extended_preds, re.form=NA)))
}

pred_summary <- extended_preds %>% group_by(year, property_type) %>%
  summarise(median_price=median(price_pred), .groups="drop")

historical <- data_model %>% group_by(year, property_type) %>%
  summarise(median_price=median(price), .groups="drop") %>% mutate(type="Historical")

timeline <- bind_rows(
  historical %>% rename(price=median_price) %>% select(year, property_type, price, type),
  pred_summary %>% rename(price=median_price) %>% mutate(type="Predicted")
)

p_timeline <- ggplot(timeline, aes(year, price, color=property_type, linetype=type)) +
  geom_line(size=1) + geom_vline(xintercept=2019.5, linetype="dashed", alpha=0.5) +
  scale_y_continuous(labels=scales::comma) +
  labs(title="Cornwall Trends: 2011-2026", x="Year", y="Median Price (£)",
       color="Type", linetype="") +
  theme_minimal() + theme(legend.position="bottom")

ggsave("timeline_2011_2026.png", p_timeline, width=10, height=6, dpi=300)
cat("Saved: timeline_2011_2026.png\n\n")

# === EXAMPLES ===

cat(strrep("=", 70), "\n", "EXAMPLES\n", strrep("=", 70), "\n\n")

ex_town <- levels(data_model$towncity)[1]

# Example 1: Single prediction
cat("Example 1: 100m² Semi-Detached, Old, Energy=60, 2019\n")
result <- predict_value(100, 60, "Semi-Detached", "Old", ex_town, 2019)
cat(sprintf("  Value: £%s", format(round(result$estimate), big.mark=",")))
if (!is.na(result$lower)) {
  cat(sprintf(" (95%% CI: £%s-£%s)\n\n",
              format(round(result$lower), big.mark = ","),
              format(round(result$upper), big.mark = ",")))
} else {
  cat(" (CI n/a for LMM)\n\n")
}

# Example 2: Trend
cat("Example 2: Trend 2015-2022\n")
trend <- predict_trend(100, 60, "Semi-Detached", "Old", ex_town, 2015:2022)

# FIXED: only round prediction columns, leave year alone
trend_print <- trend %>%
  mutate(across(c(predicted_price, lower, upper),
                ~ if (all(is.na(.))) NA_real_ else round(., 0)))
print(trend_print)
cat("\n")

# Example 3: Improvement
cat("Example 3: +20m² +15 energy to 80m² Terraced\n")
improvement <- predict_improvement(80, 50, 20, 15, "Terraced", "Old", ex_town, 2019)
cat(sprintf("  Current: £%s, Improved: £%s, Increase: £%s (+%.1f%%)\n\n",
            format(round(improvement$current), big.mark=","),
            format(round(improvement$improved), big.mark=","),
            format(round(improvement$increase), big.mark=","),
            improvement$pct_increase))

# Example 4: Model comparison
cat("Example 4: Compare all models\n")
test_property <- list(house_size=100, energy_rating=60, property_type="Semi-Detached",
                      build_type="Old", towncity=ex_town, year=2019)
model_predictions <- tibble(
  Model = names(models),
  Value = map_dbl(models, ~tryCatch(
    do.call(predict_value, c(test_property, list(model=.x)))$estimate,
    error=function(e) NA_real_))
) %>% mutate(Value_fmt=format(round(Value), big.mark=","),
             Diff_M2=Value-Value[Model=="M2_size_type"],
             Diff_pct=100*Diff_M2/Value[Model=="M2_size_type"])

print(model_predictions %>% select(Model, Value_fmt, Diff_pct))
cat("\n")

# Generate plots
ggsave("example_trend.png", plot_trend(trend), width=8, height=5, dpi=300)
ggsave("example_types.png", plot_by_type(data_model), width=8, height=5, dpi=300)
ggsave("example_locations.png", plot_locations(data_model), width=8, height=6, dpi=300)
cat("Saved: example_*.png\n\n")

# === SAVE OUTPUTS ===

saveRDS(models, "all_models.rds")
saveRDS(final_model, "final_model.rds")
saveRDS(data_model, "data_model.rds")
write_csv(model_comparison, "model_comparison.csv")
write_csv(pred_summary, "predictions_2020_2026.csv")
write_csv(model_predictions, "model_predictions.csv")
write_csv(rmse_comparison, "rmse_comparison.csv")  # NEW: save RMSE table

key_coefs <- tidy(final_model, conf.int=TRUE) %>%
  filter(term %in% c("house_size", "energy_rating", "build_typeNew", "year") |
           str_detect(term, "property_type|house_size:property")) %>%
  mutate(effect_pct=100*(exp(estimate)-1))
write_csv(key_coefs, "coefficients.csv")

cat("Saved: all_models.rds, final_model.rds, data_model.rds\n")
cat("       model_comparison.csv, rmse_comparison.csv\n")
cat("       model_predictions.csv, coefficients.csv\n\n")

# === SUMMARY ===

cat(strrep("=", 70), "\n", "COMPLETE\n", strrep("=", 70), "\n\n")
cat("Models:\n")
walk2(names(models), map_dbl(models, AIC),
      ~cat(sprintf("  %s: AIC=%.1f\n", .x, .y)))
cat("\nSelected:", names(models)[which(map_lgl(models, identical, final_model))], "\n")
if(inherits(final_model, "glm")) {
  cat(sprintf("R²=%.4f, AIC=%.1f\n",
              1-final_model$deviance/final_model$null.deviance, AIC(final_model)))
} else if(inherits(final_model, "lmerMod")) {
  cat(sprintf("AIC=%.1f, Towns=%d\n", AIC(final_model), ngrps(final_model)))
}
cat(sprintf("\nEffects: Energy %+.3f%%", energy_pct))
if(!is.na(build_pct)) cat(sprintf(", New build %+.1f%%", build_pct))
cat("\n\nTo switch: final_model <- models$M4_spline_size\n")
cat(strrep("=", 70), "\n")
