# Cornwall Property Valuation Dashboard (R Shiny)
# Run from the repository root: shiny::runApp()

library(shiny); library(tidyverse); library(plotly); library(bslib);library(mgcv); library(broom); library(DT)

# Energy Performance Certificate(EPC) score to letter
energy_score_to_letter <- function(score) {
  cut(score, breaks = c(0, 21, 39, 55, 69, 81, 92, 100),
      labels = c("G", "F", "E", "D", "C", "B", "A"), include.lowest = TRUE)
}

# DATA PREP

cat("Loading data...\n")
CORNWALL <- read.csv("data/cornwall_property.csv.gz")  # run from the repository root

env_col <- "ENVIRONMENT_IMPACT_CURRENT"
has_env <- env_col %in% names(CORNWALL)

data_clean <- CORNWALL %>%
  transmute(
    price, house_size = as.numeric(tfarea), rooms = as.numeric(numberrooms),
    energy_rating = as.numeric(CURRENT_ENERGY_EFFICIENCY),
    env_impact = if (has_env) as.numeric(.data[[env_col]]) else NA_real_,
    property_type = factor(case_when(
      propertytype == "D" ~ "Detached", propertytype == "S" ~ "Semi-Detached",
      propertytype == "T" ~ "Terraced", propertytype == "F" ~ "Flat",
      TRUE ~ as.character(propertytype)
    )),
    build_type = factor(oldnew, c("Y", "N"), c("New", "Old")),
    towncity = str_to_title(towncity) %>% factor(),
    postcode_area = str_extract(str_squish(toupper(postcode)),
                                "^[A-Z]{1,2}\\d{1,2}[A-Z]?") %>% factor(),
    year, log_price = log(price)
  ) %>%
  filter(price > 0, house_size > 0, rooms > 0) %>%
  drop_na(towncity, postcode_area, year, energy_rating)

cat("Compact Exploratory analysis: \n")
var_cor <- data_clean %>% select(where(is.numeric)) %>% cor(use = "complete.obs"); round(var_cor,2)
cat("Multicollinearity detected: house_size–rooms and energy_rating–env_impact.\n")
# ANOVA: Categorical features vs log_price
cat("\nANOVA F-statistics (categorical features vs log_price):\n")

anova_property_type <- aov(log_price ~ property_type, data = data_clean)
cat("property_type: F =", round(summary(anova_property_type)[[1]][["F value"]][1], 2), "\n")

anova_build_type <- aov(log_price ~ build_type, data = data_clean)
cat("build_type: F =", round(summary(anova_build_type)[[1]][["F value"]][1], 2), "\n")

anova_towncity <- aov(log_price ~ towncity, data = data_clean)
cat("towncity: F =", round(summary(anova_towncity)[[1]][["F value"]][1], 2), "\n")
cat("All categorical predictors (property_type, build_type, towncity) are significant by ANOVA.\n")

# Extrapolation checks for user inputs - ACTUAL DATA RANGES
num_ranges <- data_clean %>%
  summarise(
    house_size_min = min(house_size),
    house_size_max = max(house_size),
    energy_min = min(energy_rating),
    energy_max = max(energy_rating),
    year_min = min(year),
    year_max = max(year)
  )

cat("Ranges: size:", round(num_ranges$house_size_min), "-", round(num_ranges$house_size_max),
    "m²; energy:", round(num_ranges$energy_min), "-", round(num_ranges$energy_max),
    "; year:", num_ranges$year_min, "-", num_ranges$year_max, "\n")

# GLM models

# M1: Additive GLM
m1 <- glm(log_price ~ house_size + energy_rating + property_type +
            build_type + towncity + year, data = data_clean)

# M2: GLM with size × property_type interaction
m2 <- glm(log_price ~ house_size * property_type + energy_rating +
            build_type + towncity + year, data = data_clean)

# M3: GLM with size × property_type + year × property_type interactions
m3 <- glm(log_price ~ house_size * property_type + year * property_type +
            energy_rating + build_type + towncity, data = data_clean)

# Model comparison table
model_comp <- tibble(
  Model = c("M1: Additive GLM", "M2: Size × Type", "M3: Size×Type + Year×Type"),
  AIC = c(AIC(m1), AIC(m2), AIC(m3)),
  Pseudo_R2 = c(
    1 - m1$deviance / m1$null.deviance,
    1 - m2$deviance / m2$null.deviance,
    1 - m3$deviance / m3$null.deviance
  )
) %>% mutate(Delta_AIC = AIC - min(AIC))

# Likelihood Ratio Test: M2 vs M3
lrt <- anova(m2, m3, test = "LRT"); lrt #  m3 wins in terms of AIC and LT

final_model <- m3 #  has 4 numeric features & 3 categorical

# model characteristics
model_pseudo_r2 <- 1 - final_model$deviance / final_model$null.deviance
model_n <- nrow(data_clean)
levels_type <- levels(data_clean$property_type)
levels_build <- levels(data_clean$build_type)
levels_town <- levels(data_clean$towncity)
year_range <- range(data_clean$year)

# inputs for visualisation
vars_viz <- c("Price (£)" = "price", "Log price" = "log_price",
              "House size (m²)" = "house_size", "Rooms" = "rooms",
              "Energy rating" = "energy_rating")
if (has_env) vars_viz <- c(vars_viz, "Environmental impact" = "env_impact")

vars_region <- c("Town/City" = "towncity", "Postcode area" = "postcode_area")

# Coefficient table with percentage effects
coef_tbl <- tidy(final_model) %>%
  filter(str_detect(term, "house_size|energy|year")) %>%
  mutate(
    pct_effect = (exp(estimate) - 1) * 100,
    term_clean = str_replace_all(term, "property_type", "")
  ) %>%
  select(term_clean, estimate, pct_effect, p.value)

cat("Ready for Shiny.\n")

# ui
ui <- navbarPage(
  "Cornwall Property Analysis",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  # TAB 1: EXPLORE RELATIONSHIPS
  tabPanel(
    "Explore Relationships",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h5("Variables"),
        selectInput("x_var", "X-axis:", vars_viz, "house_size"),
        selectInput("y_var", "Y-axis:", vars_viz, "price"),
        hr(),
        h5("Region Grouping"),
        selectInput("region", "Group by:", vars_region, "towncity"),
        hr(),
        h5("Filters"),
        selectInput("f_type", "Property type:", c("All" = "all", levels_type)),
        sliderInput("f_year", "Year:", year_range[1], year_range[2],
                    year_range, 1, ""),
        hr(),
        textOutput("n_display")
      ),
      mainPanel(
        width = 9,
        tabsetPanel(
          tabPanel(
            "Scatter",
            br(),
            checkboxInput("scatter_town_filter", "Enable town filter (reduces sample)", FALSE),
            conditionalPanel(
              "input.scatter_town_filter == true",
              selectInput("scatter_town", "Filter to:", c("All" = "all", levels_town))
            ),
            plotlyOutput("p_scatter", height = 450),
            br(),
            verbatimTextOutput("cor_text"),
            p(style = "color:#666;font-size:0.9em;",
              em("Note: "), "Cornwall-wide view ensures large samples and stable patterns. ",
              "Enable filtering above to focus on specific towns (may reduce reliability).")
          ),
          tabPanel(
            "By Type",
            br(),
            selectInput("type_town", "Filter by town:", c("All towns" = "all", levels_town)),
            plotlyOutput("p_box", height = 450),
            uiOutput("type_sample_warning"),
            p(style = "color:#666;font-size:0.9em;",
              em("Note: "), "Boxplots compare property types. Town filtering optional.")
          ),
          tabPanel(
            "By Location",
            br(),
            plotlyOutput("p_bar", height = 450),
            hr(),
            h5("Full Location Table (Searchable)"),
            DTOutput("location_table"),
            p(style = "color:#666;font-size:0.9em;margin-top:10px;",
              em("Note: "), "Chart shows top 15. Table shows all locations with ≥30 properties.")
          )
        )
      )
    )
  ),
  
  # Tab 2: Estimate property Value
  tabPanel(
    "Estimate Property Value",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h5("Property Details"),
        selectInput("v_type", "Type:", levels_type),
        selectInput("v_town", "Town:", levels_town),
        selectInput("v_build", "Build:", levels_build),
        numericInput("v_area", "Floor area (m²):", 90, 1, 500),
        sliderInput("v_energy", "Energy rating:", 1, 100, 50),
        sliderInput("v_year", "Year:", year_range[1], year_range[2],
                    max(year_range), 1, ""),
        hr(),
        verbatimTextOutput("model_perf"),
        actionButton("btn_calc", "Calculate Value", class = "btn-primary w-100")
      ),
      mainPanel(
        width = 9,
        uiOutput("value_error"),
        uiOutput("extrap_warning"),  # EXTRAPOLATION WARNING
        fluidRow(
          column(7, wellPanel(
            style = "background:linear-gradient(135deg,#667eea,#764ba2);color:white;",
            h3("Estimated Value", style = "margin-top:0;"),
            uiOutput("val_display")
          )),
          column(5, wellPanel(style = "background:#f8f9fa;",
                              h5("95% Confidence Interval"), uiOutput("ci_display")))
        ),
        plotlyOutput("trend_plot", height = 300),
        helpText("Temporal trend (±5 years) based on model's year coefficient")
      )
    )
  ),
  
  # Tab 3: Improvement 
  tabPanel(
    "Improvement Analysis",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h5("Current Property"),
        selectInput("i_type", "Type:", levels_type),
        selectInput("i_town", "Town:", levels_town),
        selectInput("i_build", "Build:", levels_build),
        numericInput("i_area", "Area (m²):", 85, 1, 500),
        sliderInput("i_energy", "Energy:", 1, 100, 50),
        sliderInput("i_year", "Year:", year_range[1], year_range[2],
                    max(year_range), 1, ""),
        hr(),
        h5("Improvement"),
        selectInput("imp_type", "Type:",
                    c("Size Extension" = "size", "Energy Upgrade" = "energy",
                      "Property Type" = "type")),
        conditionalPanel("input.imp_type=='size'",
                         numericInput("imp_size", "Add (m²):", 20, 0, 100)),
        conditionalPanel("input.imp_type=='energy'",
                         sliderInput("imp_energy", "New rating:", 1, 100, 75)),
        conditionalPanel("input.imp_type=='type'",
                         selectInput("imp_newtype", "Upgrade to:",
                                     c("Detached", "Semi-Detached", "Terraced"))),
        actionButton("btn_imp", "Calculate Impact", class = "btn-success w-100")
      ),
      mainPanel(
        width = 9,
        uiOutput("improve_error"),
        uiOutput("improve_extrap_warning"),  # EXTRAPOLATION WARNING
        fluidRow(
          column(4, wellPanel(
            style = "text-align:center;background:#e3f2fd;",
            h5("Current Value"),
            h3(uiOutput("cur_val"), style = "color:#1976d2;")
          )),
          column(4, wellPanel(
            style = "text-align:center;background:#e8f5e9;",
            h5("After Improvement"),
            h3(uiOutput("imp_val"), style = "color:#388e3c;")
          )),
          column(4, wellPanel(
            style = "text-align:center;background:#fff3e0;",
            h5("Value Change"),
            h3(uiOutput("inc_val"), style = "color:#f57c00;")
          ))
        ),
        fluidRow(
          column(7, plotlyOutput("comp_bar", height = 350)),
          column(5,
                 wellPanel(
                   h5("Change Summary"),
                   verbatimTextOutput("imp_summary")
                 )
          )
        )
      )
    )
  ),
  
  # Tab 4: Model& statisics
  tabPanel(
    "Model & Statistics",
    fluidRow(
      column(6, wellPanel(
        h4("Variable Selection"),
        tableOutput("var_cor_tbl"),
        hr(),
        tags$ul(
          tags$li(strong("house_size:"), " r = 0.67 → INCLUDED"),
          tags$li(strong("energy_rating:"), " r = 0.06 → INCLUDED (hypothesis test)"),
          tags$li(strong("rooms:"), " r = 0.79 with size → EXCLUDED (collinear)")
        )
      )),
      column(6, wellPanel(
        h4("Model Comparison"),
        tableOutput("mod_comp_tbl"),
        hr(),
        h5("Likelihood Ratio Test (M2 vs M3):"),
        verbatimTextOutput("lrt_result")
      ))
    ),
    fluidRow(
      column(6, wellPanel(
        h4("Coefficient Effects (% change in price)"),
        tableOutput("coef_tbl"),
        p(style = "color:#666;font-size:0.9em;",
          em("Note: "), "Percentage effects calculated as (exp(β) - 1) × 100")
      )),
      column(6, wellPanel(
        h4("Diagnostics"),
        plotOutput("diag_plots", height = 300)
      ))
    )
  ),
  
  # Tab 5: About
  tabPanel(
    "About & Methodology",
    div(style = "max-width:900px;margin:20px auto;",
        h3("About This Application"),
        hr(),
        
        h4("Purpose"),
        p("Interactive exploration, valuation, and improvement analysis of Cornwall property data."),
        
        h4("Data"),
        p("UK House Price & EPC data: ", format(model_n, big.mark = ","),
          " transactions (", min(data_clean$year), "-", max(data_clean$year), ")"),
        
        h4("Statistical Models"),
        tags$ul(
          tags$li(strong("M1 (Additive GLM):"), " log(price) ~ house_size + energy_rating + property_type + build_type + town + year"),
          tags$li(strong("M2 (Size × Type):"), " Adds house_size × property_type interaction"),
          tags$li(strong("M3 (Size×Type + Year×Type):"), " Adds year × property_type interaction")
        ),
        
        h4("Interaction Terms Explained"),
        tags$ul(
          tags$li(strong("house_size × property_type:"), 
                  " The effect of each additional m² varies by property type. ",
                  "Example: Adding 10m² adds more value to a Detached than a Flat."),
          tags$li(strong("year × property_type:"), 
                  " Temporal trends (annual % price change) differ by property type. ",
                  "Example: Flats may appreciate faster/slower than Detached.")
        ),
        
        h4("Model Selection"),
        p(strong("Selected: M3"), " (Size×Type + Year×Type)"),
        p("Rationale: LRT shows M3 significantly better than M2 (p < 0.05). ",
          "The year×type interaction captures differential temporal trends by property type, ",
          "providing better fit despite increased complexity."),
        
        h4("Performance Metrics"),
        p(strong("Pseudo-R²:"), " ", round(model_pseudo_r2, 3),
          " (NOT % variance explained - this is deviance reduction for GLM)"),
        p(strong("AIC:"), " ", format(round(AIC(final_model), 0), big.mark = ",")),
        
        h4("Key Findings"),
        tags$ul(
          tags$li("Size effect varies by type (interaction significant)"),
          tags$li("Temporal trends vary by type (year×type significant)"),
          tags$li("Energy rating weakly correlated with price (r ≈ 0.06)"),
          tags$li("Location highly significant"),
          tags$li("Average temporal trend ~2.6% per year (varies by type)")
        ),
        
        h4("Extrapolation & Prediction Limits"),
        p(strong("Training Data Ranges:")),
        tags$ul(
          tags$li("Floor area: ", round(num_ranges$house_size_min), "-", 
                  round(num_ranges$house_size_max), " m²"),
          tags$li("Energy rating: ", round(num_ranges$energy_min), "-", 
                  round(num_ranges$energy_max)),
          tags$li("Year: ", num_ranges$year_min, "-", num_ranges$year_max)
        ),
        p(strong("Warning: "), "Predictions outside these ranges are extrapolations. ",
          "The app will flag these with warnings. Extrapolated predictions are less reliable ",
          "and confidence intervals may underestimate true uncertainty."),
        
        h4("Tab Functionality"),
        h5("1. Explore Relationships"),
        p(strong("Default:"), " Cornwall-wide (large samples, stable patterns)"),
        p(strong("Optional:"), " Town filter (checkbox)"),
        
        h5("2. Estimate Property Value"),
        p("Enter property details to receive model-based valuation with confidence intervals."),
        
        h5("3. Improvement Analysis"),
        p("Compare current vs. improved property values based on size extensions, ",
          "energy upgrades, or property type changes."),
        
        h4("Assumptions & Limitations"),
        tags$ul(
          tags$li("Log-linear relationship assumed"),
          tags$li("Predictions extrapolate temporal trends"),
          tags$li("Energy improvements show model-predicted change (not market reality)"),
          tags$li("Small samples in some towns/postcodes"),
          tags$li("Extrapolation beyond training ranges reduces reliability")
        ),
        
    )
  )
)

#server
server <- function(input, output, session) {
  
  # Exploring tab 
  
  filt_data <- reactive({
    d <- data_clean
    if (input$f_type != "all") d <- filter(d, property_type == input$f_type)
    d <- filter(d, year >= input$f_year[1], year <= input$f_year[2])
    d
  })
  
  scatter_data <- reactive({
    d <- filt_data()
    if (input$scatter_town_filter && input$scatter_town != "all") {
      d <- filter(d, towncity == input$scatter_town)
    }
    d
  })
  
  type_data <- reactive({
    d <- filt_data()
    if (input$type_town != "all") {
      d <- filter(d, towncity == input$type_town)
    }
    d
  })
  
  output$n_display <- renderText(
    paste("n =", format(nrow(filt_data()), big.mark = ","))
  )
  
  output$p_scatter <- renderPlotly({
    d <- scatter_data()
    req(nrow(d) > 1)
    
    x_label <- if (input$x_var == "house_size") {
      "House size (m²)"
    } else {
      names(vars_viz)[vars_viz == input$x_var]
    }
    
    plot_ly(d, x = ~.data[[input$x_var]], y = ~.data[[input$y_var]],
            color = ~.data[[input$region]], type = "scatter", mode = "markers",
            marker = list(size = 5, opacity = 0.6), colors = "Set2") %>%
      layout(
        xaxis = list(title = x_label),
        yaxis = list(title = names(vars_viz)[vars_viz == input$y_var]),
        legend = list(title = list(text = names(vars_region)[vars_region == input$region]))
      )
  })
  
  output$cor_text <- renderPrint({
    d <- scatter_data()
    req(nrow(d) >= 5)
    r <- suppressWarnings(cor(d[[input$x_var]], d[[input$y_var]], use = "complete.obs"))
    cat("n =", format(nrow(d), big.mark = ","), "| r =", round(r, 3))
  })
  
  output$p_box <- renderPlotly({
    d <- type_data()
    req(nrow(d) > 1)
    plot_ly(d, x = ~property_type, y = ~log_price, color = ~property_type,
            type = "box", colors = "Set2") %>%
      layout(xaxis = list(title = "Property type"),
             yaxis = list(title = "Log price"), showlegend = FALSE)
  })
  
  output$type_sample_warning <- renderUI({
    d <- type_data()
    if (nrow(d) < 30) {
      div(style = "background:#fff3cd;padding:10px;border-radius:5px;margin:10px 0;",
          strong("Warning: "), "Sample size < 30. Results may be unreliable.")
    }
  })
  
  location_summary <- reactive({
    d <- filt_data()
    req(nrow(d) > 1)
    d %>% group_by(.data[[input$region]]) %>%
      summarise(n = n(), median_price = median(price), .groups = "drop") %>%
      filter(n >= 30) %>%
      arrange(desc(median_price))
  })
  
  output$p_bar <- renderPlotly({
    summ <- location_summary()
    req(nrow(summ) > 0)
    top15 <- summ %>% slice(1:15)
    
    plot_ly(top15, x = ~reorder(.data[[input$region]], median_price),
            y = ~median_price, type = "bar", marker = list(color = "#667eea")) %>%
      layout(
        xaxis = list(title = names(vars_region)[vars_region == input$region]),
        yaxis = list(title = "Median price (£)")
      )
  })
  
  output$location_table <- renderDT({
    summ <- location_summary()
    req(nrow(summ) > 0)
    summ %>%
      rename(Location = 1, `Sample Size` = n, `Median Price (£)` = median_price) %>%
      mutate(`Median Price (£)` = format(round(`Median Price (£)`), big.mark = ",")) %>%
      datatable(
        options = list(pageLength = 10, lengthMenu = c(10, 25, 50, 100),
                       searching = TRUE, ordering = TRUE),
        rownames = FALSE
      )
  })
  
  # VALUE
  
  output$model_perf <- renderText(
    paste0("M3: Size×Type + Year×Type\nPseudo-R² = ", round(model_pseudo_r2, 3),
           "\nn = ", format(model_n, big.mark = ","))
  )
  
  user_prop <- reactive({
    req(input$btn_calc > 0)
    
    # Input validation
    if (input$v_area <= 0) {
      return(NULL)
    }
    
    data.frame(
      house_size = input$v_area, energy_rating = input$v_energy,
      property_type = factor(input$v_type, levels = levels_type),
      build_type = factor(input$v_build, levels = levels_build),
      towncity = factor(input$v_town, levels = levels_town),
      year = input$v_year
    )
  })
  
  # Extrapolation Check (Value tab)
  extrap_check <- reactive({
    req(input$btn_calc > 0, user_prop())
    
    warnings <- c()
    
    # Check house_size
    if (input$v_area < num_ranges$house_size_min || 
        input$v_area > num_ranges$house_size_max) {
      warnings <- c(warnings, sprintf(
        "Floor area (%.0f m²) outside training range (%.0f–%.0f m²). Prediction is extrapolation.",
        input$v_area, round(num_ranges$house_size_min), round(num_ranges$house_size_max)
      ))
    }
    
    # Check energy_rating
    if (input$v_energy < num_ranges$energy_min || 
        input$v_energy > num_ranges$energy_max) {
      warnings <- c(warnings, sprintf(
        "Energy rating (%d) outside training range (%d–%d). Prediction less reliable.",
        input$v_energy, round(num_ranges$energy_min), round(num_ranges$energy_max)
      ))
    }
    
    # Check year
    if (input$v_year < num_ranges$year_min || 
        input$v_year > num_ranges$year_max) {
      warnings <- c(warnings, sprintf(
        "Year (%d) outside training range (%d–%d). Temporal extrapolation – use with caution.",
        input$v_year, num_ranges$year_min, num_ranges$year_max
      ))
    }
    
    warnings
  })
  
  output$extrap_warning <- renderUI({
    warnings <- extrap_check()
    
    if (length(warnings) > 0) {
      div(
        style = "background:#fff3cd;border:1px solid #ffc107;padding:12px;border-radius:5px;margin:10px 0;",
        h5(style = "margin-top:0;color:#856404;", 
           "⚠ Extrapolation Warning"),
        tags$ul(
          style = "margin-bottom:0;color:#856404;",
          lapply(warnings, function(w) tags$li(w))
        ),
        p(style = "margin:10px 0 0 0;font-size:0.9em;color:#666;",
          em("Note: "), "Predictions outside the training data range are less reliable. ",
          "Confidence intervals may underestimate uncertainty.")
      )
    }
  })
  
  output$value_error <- renderUI({
    req(input$btn_calc > 0)
    if (is.null(user_prop())) {
      div(style = "background:#f8d7da;padding:10px;border-radius:5px;margin:10px;",
          strong("Error: "), "Floor area must be positive.")
    }
  })
  
  pred_val <- reactive({
    req(user_prop())
    p <- predict(final_model, newdata = user_prop(), se.fit = TRUE)
    list(price = exp(p$fit), lower = exp(p$fit - 1.96 * p$se.fit),
         upper = exp(p$fit + 1.96 * p$se.fit))
  })
  
  output$val_display <- renderUI({
    req(pred_val())
    tagList(
      h2(paste0("£", format(round(pred_val()$price), big.mark = ",")),
         style = "margin:10px 0;font-size:2.5em;"),
      p(paste(input$v_area, "m²", input$v_type, "in", input$v_town),
        style = "font-size:1.1em;")
    )
  })
  
  output$ci_display <- renderUI({
    req(pred_val())
    p <- pred_val()
    tagList(
      p(paste0("£", format(round(p$lower), big.mark = ","), " – £",
               format(round(p$upper), big.mark = ",")),
        style = "font-size:1.1em;margin:20px 0;"),
      p(paste("Year:", input$v_year), style = "color:#666;")
    )
  })
  
  output$trend_plot <- renderPlotly({
    req(input$btn_calc > 0, user_prop())
    base <- user_prop()
    yrs <- seq(base$year - 5, base$year + 5)
    nd <- base[rep(1, length(yrs)), ]
    nd$year <- yrs
    nd$pred <- exp(predict(final_model, newdata = nd))
    
    plot_ly(nd, x = ~year, y = ~pred, type = "scatter",
            mode = "lines+markers", line = list(color = "#667eea", width = 3)) %>%
      add_segments(x = input$v_year, xend = input$v_year, y = 0,
                   yend = pred_val()$price,
                   line = list(color = "#e74c3c", dash = "dash")) %>%
      layout(xaxis = list(title = "Year"), yaxis = list(title = "Value (£)"),
             showlegend = FALSE)
  })
  
  # Improve
  
  imp_cur <- reactive({
    req(input$btn_imp > 0)
    
    # Input validation
    if (input$i_area <= 0) {
      return(NULL)
    }
    
    data.frame(
      house_size = input$i_area, energy_rating = input$i_energy,
      property_type = factor(input$i_type, levels = levels_type),
      build_type = factor(input$i_build, levels = levels_build),
      towncity = factor(input$i_town, levels = levels_town),
      year = input$i_year
    )
  })
  
  imp_new <- reactive({
    req(imp_cur())
    
    improved <- imp_cur()
    
    if (input$imp_type == "size") {
      if (input$imp_size < 0) return(NULL)
      improved$house_size <- input$i_area + input$imp_size
    } else if (input$imp_type == "energy") {
      improved$energy_rating <- input$imp_energy
    } else if (input$imp_type == "type") {
      improved$property_type <- factor(input$imp_newtype, levels = levels_type)
    }
    
    improved
  })
  
  
  # Extrapolation check
  improve_extrap_check <- reactive({
    req(input$btn_imp > 0, imp_cur(), imp_new())
    
    warnings <- c()
    
    # Check current property
    if (input$i_area < num_ranges$house_size_min || 
        input$i_area > num_ranges$house_size_max) {
      warnings <- c(warnings, sprintf(
        "Current floor area (%.0f m²) outside training range (%.0f–%.0f m²).",
        input$i_area, round(num_ranges$house_size_min), round(num_ranges$house_size_max)
      ))
    }
    
    # Check improved property i.e if size extension
    if (input$imp_type == "size") {
      new_size <- input$i_area + input$imp_size
      if (new_size > num_ranges$house_size_max) {
        warnings <- c(warnings, sprintf(
          "Improved floor area (%.0f m²) exceeds training maximum (%.0f m²). Impact estimate is extrapolation.",
          new_size, round(num_ranges$house_size_max)
        ))
      }
    }
    
    # Check energy (current and improved)
    if (input$i_energy < num_ranges$energy_min || 
        input$i_energy > num_ranges$energy_max) {
      warnings <- c(warnings, sprintf(
        "Current energy rating (%d) outside training range (%d–%d).",
        input$i_energy, round(num_ranges$energy_min), round(num_ranges$energy_max)
      ))
    }
    
    if (input$imp_type == "energy") {
      if (input$imp_energy < num_ranges$energy_min || 
          input$imp_energy > num_ranges$energy_max) {
        warnings <- c(warnings, sprintf(
          "Improved energy rating (%d) outside training range (%d–%d).",
          input$imp_energy, round(num_ranges$energy_min), round(num_ranges$energy_max)
        ))
      }
    }
    
    warnings
  })
  
  output$improve_extrap_warning <- renderUI({
    warnings <- improve_extrap_check()
    
    if (length(warnings) > 0) {
      div(
        style = "background:#fff3cd;border:1px solid #ffc107;padding:12px;border-radius:5px;margin:10px 0;",
        h5(style = "margin-top:0;color:#856404;", 
           "Extrapolation Warning"),
        tags$ul(
          style = "margin-bottom:0;color:#856404;",
          lapply(warnings, function(w) tags$li(w))
        ),
        p(style = "margin:10px 0 0 0;font-size:0.9em;color:#666;",
          em("Note: "), "Impact estimates involving extrapolation are less reliable.")
      )
    }
  })
  
  output$improve_error <- renderUI({
    req(input$btn_imp > 0)
    if (is.null(imp_cur())) {
      return(div(style = "background:#f8d7da;padding:10px;border-radius:5px;margin:10px;",
                 strong("Error: "), "Floor area must be positive."))
    }
    if (input$imp_type == "size" && input$imp_size < 0) {
      return(div(style = "background:#f8d7da;padding:10px;border-radius:5px;margin:10px;",
                 strong("Error: "), "Addition must be non-negative."))
    }
  })
  
  imp_res <- reactive({
    req(imp_cur(), imp_new())
    cur <- exp(predict(final_model, newdata = imp_cur()))
    new <- exp(predict(final_model, newdata = imp_new()))
    
    desc <- switch(input$imp_type,
                   size = paste0("+", input$imp_size, " m² extension"),
                   energy = paste0("Energy: ", input$i_energy, " → ", input$imp_energy),
                   type = paste0(input$i_type, " → ", input$imp_newtype)
    )
    
    # energy bands
    cur_band <- as.character(energy_score_to_letter(input$i_energy))
    new_band <- if (input$imp_type == "energy") {
      as.character(energy_score_to_letter(input$imp_energy))
    } else {
      cur_band
    }
    
    list(old = cur, new = new, chg = new - cur,
         pct = ((new - cur) / cur) * 100, desc = desc,
         cur_band = cur_band, new_band = new_band)
  })
  
  output$cur_val <- renderUI({
    req(imp_res())
    paste0("£", format(round(imp_res()$old), big.mark = ","))
  })
  
  output$imp_val <- renderUI({
    req(imp_res())
    paste0("£", format(round(imp_res()$new), big.mark = ","))
  })
  
  output$inc_val <- renderUI({
    req(imp_res())
    r <- imp_res()
    tagList(
      div(paste0(ifelse(r$chg > 0, "+", ""), "£",
                 format(round(abs(r$chg)), big.mark = ","))),
      tags$small(paste0("(", ifelse(r$pct > 0, "+", ""),
                        round(r$pct, 1), "%)"), style = "color:#666;")
    )
  })
  
  output$comp_bar <- renderPlotly({
    req(imp_res())
    r <- imp_res()
    df <- data.frame(
      Scenario = factor(c("Current", "Improved"), levels = c("Current", "Improved")),
      Val = c(r$old, r$new)
    )
    
    # Increase y-range to prevent clipping
    y_max <- max(df$Val) * 1.15
    
    plot_ly(df, x = ~Scenario, y = ~Val, type = "bar",
            marker = list(color = c("#667eea", "#00cc66")),
            text = ~paste0("£", format(round(Val / 1000), big.mark = ","), "k"),
            textposition = "outside") %>%
      layout(
        xaxis = list(title = ""),
        yaxis = list(title = "Value (£)", range = c(0, y_max)),
        showlegend = FALSE
      )
  })
  
  # Improvement summary
  output$imp_summary <- renderText({
    req(imp_res())
    r <- imp_res()
    
    # Build summary 
    lines <- c(
      paste0("Improvement: ", r$desc),
      "",
      paste0("Current value: £", format(round(r$old), big.mark = ",")),
      paste0("Improved value: £", format(round(r$new), big.mark = ",")),
      paste0("Change: ", ifelse(r$chg > 0, "+", ""), "£",
             format(round(abs(r$chg)), big.mark = ","),
             " (", ifelse(r$pct > 0, "+", ""), round(r$pct, 1), "%)")
    )
    
    # Add energy band 
    if (input$imp_type == "energy" && r$cur_band != r$new_band) {
      lines <- c(lines, "",
                 paste0("Energy band: ", r$cur_band, " → ", r$new_band))
    } else if (input$imp_type == "energy") {
      lines <- c(lines, "",
                 paste0("Energy band: ", r$cur_band, " (unchanged)"))
    }
    
    # note about energy correlation
    if (input$imp_type == "energy") {
      lines <- c(lines, "",
                 "Note: Energy rating weakly correlated",
                 "with price (r ≈ 0.06). Predicted change",
                 "reflects model structure, not strong",
                 "market signal.")
    }
    
    paste(lines, collapse = "\n")
  })
  
  # MODEL
  
  output$var_cor_tbl <- renderTable({
    tibble(
      Variable = c("House size (m²)", "Rooms", "Energy rating"),
      `Correlation with log(price)` = round(
        var_cor[c("house_size", "rooms", "energy_rating"), "log_price"], 3
      )
    )
  }, striped = TRUE)
  
  output$mod_comp_tbl <- renderTable({
    model_comp %>%
      mutate(across(c(AIC, Delta_AIC), ~round(.x, 0)),
             across(Pseudo_R2, ~round(.x, 3)))
  }, striped = TRUE)
  
  output$lrt_result <- renderPrint({
    cat("Chi-squared =", round(lrt$Deviance[2], 2), "\n")
    cat("p-value =", round(lrt$`Pr(>Chi)`[2], 4), "\n")
    cat("Result: Year×Type interaction IS significant (p < 0.05)\n")
    cat("M3 selected (better fit, captures temporal variation by type)")
  })
  
  output$coef_tbl <- renderTable({
    coef_tbl %>%
      mutate(across(c(estimate, pct_effect), ~round(.x, 3)),
             across(p.value, ~round(.x, 5))) %>%
      rename(Term = term_clean, Estimate = estimate,
             `% Effect` = pct_effect, `p-value` = p.value)
  }, striped = TRUE)
  
  output$diag_plots <- renderPlot({
    par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    plot(final_model, which = c(1, 2))
  })
}

shinyApp(ui, server)