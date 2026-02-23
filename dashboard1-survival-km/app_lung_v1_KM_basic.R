
#=============================================================
# Program Name : app_lung_v1_km_basic.R
# Description  : Interactive Kaplan-Meier Survival Analysis
#                Dashboard for NCCTG Lung Cancer Trial
#                Features: Age range filter, CI toggle,
#                Risk table, Median survival table,
#                Log-rank test
# Author       : Anbarasu Arumugham
# Date Created : 2026-02-22
# Version      : 1.0
# Input Data   : survival::lung (R built-in dataset)
#=============================================================

library(shiny)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

# Load data
lung_data           <- survival::lung
lung_data$sex_label <- ifelse(lung_data$sex == 1, "Male", "Female")
lung_data$status_km <- ifelse(lung_data$status == 2, 1, 0)

# UI
ui <- fluidPage(

  titlePanel("Clinical Survival Analysis Dashboard"),

  tags$h4(
    "NCCTG Lung Cancer Trial — Kaplan-Meier Survival Analysis",
    style = "color: #2c3e50; margin-bottom: 20px;"
  ),

  sidebarLayout(
    sidebarPanel(

      sliderInput(
        inputId = "age_range",
        label   = "Filter by Age Range:",
        min     = min(lung_data$age, na.rm = TRUE),
        max     = max(lung_data$age, na.rm = TRUE),
        value   = c(min(lung_data$age, na.rm = TRUE),
                    max(lung_data$age, na.rm = TRUE))
      ),

      checkboxInput(
        inputId = "show_ci",
        label   = "Show Confidence Intervals",
        value   = TRUE
      ),

      checkboxInput(
        inputId = "show_risk",
        label   = "Show Risk Table",
        value   = TRUE
      ),

      hr(),
      h4("Dataset Summary"),
      verbatimTextOutput("data_summary")
    ),

    mainPanel(
      plotOutput("km_plot", height = "500px"),
      hr(),
      h4("Median Survival with 95% Confidence Intervals"),
      tableOutput("median_table"),
      hr(),
      h4("Log-Rank Test"),
      verbatimTextOutput("logrank_result")
    )
  )
)

# Server
server <- function(input, output, session) {

  # Reactive filtered dataset
  filtered_data <- reactive({
    lung_data %>%
      filter(
        age >= input$age_range[1],
        age <= input$age_range[2],
        !is.na(sex_label)
      )
  })

  # KM fit - stratified by sex
  km_fit <- reactive({
    df <- filtered_data()
    survfit(Surv(time, status_km) ~ sex_label, data = df)
  })

  # KM Plot
  output$km_plot <- renderPlot({
    ggsurvplot(
      fit               = km_fit(),
      data              = filtered_data(),
      conf.int          = as.logical(input$show_ci),
      risk.table        = as.logical(input$show_risk),
      pval              = TRUE,
      pval.method       = TRUE,
      xlab              = "Time (Days)",
      ylab              = "Survival Probability",
      title             = "Kaplan-Meier Survival Curves by Sex",
      palette           = c("#e74c3c", "#2980b9"),
      ggtheme           = theme_bw(),
      risk.table.height = 0.25,
      surv.median.line  = "hv"
    )
  })

  # Median survival table
  output$median_table <- renderTable({
    s <- summary(km_fit())$table
    data.frame(
      Group                = rownames(s),
      N                    = s[, "records"],
      Events               = s[, "events"],
      Median_Survival_Days = s[, "median"],
      CI_Lower_95          = s[, "0.95LCL"],
      CI_Upper_95          = s[, "0.95UCL"]
    )
  })

  # Log-rank test
  output$logrank_result <- renderPrint({
    df   <- filtered_data()
    test <- survdiff(Surv(time, status_km) ~ sex_label, data = df)
    print(test)
  })

  # Dataset summary
  output$data_summary <- renderPrint({
    df <- filtered_data()
    cat("Total Patients :", nrow(df), "\n")
    cat("Events (Deaths):", sum(df$status_km), "\n")
    cat("Censored       :", sum(df$status_km == 0), "\n")
    cat("Median Follow-up:", median(df$time, na.rm = TRUE), "days\n")
  })
}

shinyApp(ui = ui, server = server)
