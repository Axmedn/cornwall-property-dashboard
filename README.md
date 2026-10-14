# Cornwall Property Valuation Dashboard

An interactive R Shiny dashboard for exploring and estimating Cornwall house prices, built on 62,711 property sales from 2011 to 2019. It is designed for people without a statistics background: every estimate comes with a confidence interval, and the app warns when a sample is small or when an input falls outside the range of the data the model was fitted on.

Individual course project, University of Glasgow (2025).

## What the app does

| Tab | Purpose |
|---|---|
| Explore Relationships | Scatter plots, box plots by property type, and median prices by town or postcode district, with filters for property type and year. Flags any group with fewer than 30 sales. |
| Estimate Property Value | Enter type, town, build, floor area, energy rating and year to get a model-based price with a 95% confidence interval and a ±5-year trend line. |
| Improvement Analysis | Compares the estimated value before and after an extension, an energy upgrade or a change of property type. |
| Model & Statistics | Correlations used for variable selection, the model comparison table, the likelihood-ratio test, coefficient effects as % changes in price, and residual diagnostics. |
| About & Methodology | Data, model choice, training ranges, assumptions and limitations in plain language. |

## Method

1. **Cleaning.** Recoded property type and new/old build, standardised town names, extracted the postcode district, and dropped sales with non-positive price, floor area or room count.
2. **Variable selection.** Floor area and number of rooms are strongly correlated, so rooms was left out. Energy efficiency was kept for testing despite a weak correlation with price.
3. **Models.** Three Gaussian GLMs on log price:
   - M1: additive (floor area, energy rating, property type, build type, town, year)
   - M2: M1 plus a floor area × property type interaction
   - M3: M2 plus a year × property type interaction
4. **Selection.** Compared by AIC and a likelihood-ratio test of M2 against M3. M3 was selected and drives every estimate in the app.
5. **Prediction.** Estimates are back-transformed from the log scale, with 95% confidence intervals for the expected price.

`analysis/model_exploration.R` is the wider exploration behind the app: it also fits a natural-spline model and a mixed model with a random intercept for town, and compares all five on a held-out 20% test set.

## Limitations

- The data stops in 2019, so any year after that is an extrapolation, and the app says so.
- Residual plots for the exploratory fits show heavier tails than a normal distribution, so intervals at the extremes of the price range should be read with care.
- Some towns and postcode districts have few sales, so their estimates are less reliable.
- An energy-upgrade estimate reflects the model's structure, not a measured market effect; energy rating is only weakly related to price.

## Data

HM Land Registry price-paid records linked to Energy Performance Certificate (EPC) data, supplied for the course. The copy in `data/` keeps only the ten columns the code uses. Street addresses have been removed and postcodes cut to the district (for example, `TR1`), so no individual property can be identified.

## How to run

Requires R 4.1 or later.

```r
install.packages(c("shiny", "bslib", "tidyverse", "plotly", "DT", "broom", "mgcv"))
shiny::runApp()   # from the repository root
```

The model exploration script also needs `splines`, `lme4` and `patchwork`. Run it from the repository root.

## Files

```
app.R                          Shiny app (data preparation, models, UI and server)
analysis/model_exploration.R   Five-model comparison, test-set RMSE, prediction helpers
data/cornwall_property.csv.gz  Cleaned extract used by both scripts
```
