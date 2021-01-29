# Packages
require(dplyr)
require(magrittr)
require(readr)
require(tibble)
require(shiny)
require(shinyWidgets)
require(ggplot2)
require(tidybayes)

inverse_logit <- function(logit){
    odds <- exp(logit)
    prob <- odds / (1 + odds)
    return(prob)
}

# Read posteriors
# post_df <- read_rds("Bayesian_Hammer/hammer_posteriors.RDS")
post_df <- read_rds("hammer_posteriors.RDS")

# Prediction function
predict_posterior <- function(pst_df, sx = FALSE, gsurg = FALSE, csurg = FALSE, 
                              hyper = FALSE, apache = FALSE, fluid = FALSE, amb = FALSE,
                              los = FALSE)
{
    # Combine variables
    varlist <- c(sx, gsurg, csurg, hyper, apache, fluid, amb, los)
    
    # Create coefficient list
    coef_names <- names(pst_df)[2:9]
    
    # Select relevant variables
    relevant_variables <- coef_names[varlist]
    
    # Extract from posterior data frame
    pred_df <- pst_df %>% 
        select(all_of(c("b_Intercept",relevant_variables)))
    
    # Create prediction posterior
    predictions <- pred_df %>% rowSums()
    
    # Transform to probabilities and output
    return(inverse_logit(predictions))
}

# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel("Bayesian Estimation of Risk Uncertainty from Hammer model"),
    #---------
    sidebarLayout(sidebarPanel(
                    switchInput(
                        inputId = "male",
                        label = "Sex", 
                        labelWidth = "120px",
                        offLabel = "Female",
                        onLabel = "Male",
                        onStatus = "danger", 
                        offStatus = "success"
                    ),
                    materialSwitch(
                        inputId = "missing.male",
                        label = "", 
                        value = FALSE,
                        status = "warning",
                    ),
                    switchInput(
                        inputId = "general_surgery",
                        label = "General Surgery", 
                        labelWidth = "120px",
                        offLabel = "No",
                        onLabel = "Yes",
                        onStatus = "danger", 
                        offStatus = "success"
                    ),
                    materialSwitch(
                        inputId = "missing.general_surgery",
                        label = "", 
                        value = FALSE,
                        status = "warning",
                    ),
                    switchInput(
                        inputId = "cardiac_surgery",
                        label = "Cardiac Surgery", 
                        labelWidth = "120px",
                        offLabel = "No",
                        onLabel = "Yes",
                        onStatus = "success", 
                        offStatus = "danger"
                    ),
                    materialSwitch(
                        inputId = "missing.cardiac_surgery",
                        label = "", 
                        value = FALSE,
                        status = "warning",
                    ),
                    switchInput(
                        inputId = "hyperglycemia",
                        label = "Hyperglycaemia", 
                        labelWidth = "120px",
                        offLabel = "No",
                        onLabel = "Yes",
                        onStatus = "danger", 
                        offStatus = "success"
                    ),
                    materialSwitch(
                        inputId = "missing.hyperglycemia",
                        label = "", 
                        value = FALSE,
                        status = "warning",
                    ),
                    switchInput(
                        inputId = "high_apache",
                        label = "APACHE-II > 20", 
                        labelWidth = "120px",
                        offLabel = "No",
                        onLabel = "Yes",
                        onStatus = "danger", 
                        offStatus = "success"
                    ),
                    materialSwitch(
                        inputId = "missing.high_apache",
                        label = "", 
                        value = FALSE,
                        status = "warning",
                    ),
                    switchInput(
                        inputId = "fluid_balance_5L",
                        label = "Positive fluid balance > 5L", 
                        labelWidth = "120px",
                        offLabel = "No",
                        onLabel = "Yes",
                        onStatus = "danger", 
                        offStatus = "success"
                    ),
                    materialSwitch(
                        inputId = "missing.fluid_balance_5L",
                        label = "", 
                        value = FALSE,
                        status = "warning",
                    ),
                    switchInput(
                        inputId = "ambulation",
                        label = "Ambulation", 
                        labelWidth = "120px",
                        offLabel = "Yes",
                        onLabel = "No",
                        onStatus = "danger", 
                        offStatus = "success"
                    ),
                    materialSwitch(
                        inputId = "missing.ambulation",
                        label = "", 
                        value = FALSE,
                        status = "warning",
                    ),
                    switchInput(
                        inputId = "los_5",
                        label = "Length of stay > 5 days", 
                        labelWidth = "120px",
                        offLabel = "No",
                        onLabel = "Yes",
                        onStatus = "danger", 
                        offStatus = "success"
                    ),
                    materialSwitch(
                        inputId = "missing.los_5",
                        label = "", 
                        value = FALSE,
                        status = "warning",
                    )
             ),
             #-----------
             mainPanel(plotOutput("distPlot"),
                       sliderInput("CI",
                                   "Highest Density Interval (%):",
                                   min = 5,
                                   max = 100,
                                   value = 95,
                                   step = 5),
                       materialSwitch(
                           inputId = "plotZoom",
                           label = "Zoom plot",
                           status = "primary",
                           right = TRUE
                       )))
)

# Debug
input <- data.frame(
    male = TRUE,
    general_surgery = TRUE,
    cardiac_surgery = FALSE,
    hyperglycemia = FALSE,
    high_apache = TRUE,
    fluid_balance_5L = FALSE,
    ambulation = TRUE,
    los_5 = TRUE
)
input$CI <- 95
input$plotZoom <- TRUE

# Define server logic required to draw a histogram
server <- function(input, output) {
    output$distPlot <- renderPlot({
        
        # Produce prediction posterior
        probs <- predict_posterior(post_df, input$male, input$general_surgery, input$cardiac_surgery,
                          input$hyperglycemia, input$high_apache, input$fluid_balance_5L,
                          !input$ambulation, input$los_5)
        # Credible interval
        if(input$CI != 100)
        {
            CI_width <- input$CI / 100
        }else
        {
            CI_width <- 0.99999999
        }
        risk_CI <-  point_interval(probs * 100, .point = Mode,
                                   .interval = hdci, .width = CI_width)
        # Obtain color
        fill_value <- ifelse(risk_CI$y <= 5, "#4daf4a",
                             ifelse(risk_CI$y <= 15, "#ffff33",
                                    ifelse(risk_CI$y <= 50, "#ff7f00", "#e41a1c")))
        
        # Draw plot
        if(input$plotZoom == TRUE)
        {
            outplot <- data.frame(probs = probs * 100) %>% 
                ggplot(aes(x = probs))+
                geom_density(fill = fill_value, alpha = 0.5)+
                theme_classic()+
                labs(x = "Predicted risk of ICU readmission",
                     y = "Confidence",
                     title = paste("Risk = ",
                                   round(risk_CI$y, 2), "% [",
                                   round(risk_CI$ymin, 2), ", ",
                                   round(risk_CI$ymax, 2), "]", sep = ""))+
                theme(axis.text.y = element_blank(),
                      axis.ticks.y = element_blank())+
                geom_vline(xintercept = risk_CI$y)+
                geom_vline(xintercept = c(risk_CI$ymin,
                                          risk_CI$ymax),
                           linetype = "dashed")+
                theme(plot.title = element_text(size=26))
            
        }else
        {
            outplot <- data.frame(probs = probs * 100) %>% 
                ggplot(aes(x = probs))+
                geom_density(fill = fill_value, alpha = 0.5)+
                theme_classic()+
                labs(x = "Predicted risk of ICU readmission",
                     y = "Confidence",
                     title = paste("Risk = ",
                                   round(risk_CI$y, 2), "% [",
                                   round(risk_CI$ymin, 2), ", ",
                                   round(risk_CI$ymax, 2), "]", sep = ""))+
                theme(axis.text.y = element_blank(),
                      axis.ticks.y = element_blank())+
                geom_vline(xintercept = risk_CI$y)+
                geom_vline(xintercept = c(risk_CI$ymin,
                                          risk_CI$ymax),
                           linetype = "dashed")+
                coord_cartesian(xlim=c(0, 100)) +
                theme(plot.title = element_text(size=26))
        }
        outplot
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
