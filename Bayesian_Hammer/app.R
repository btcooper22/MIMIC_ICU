# Packages
require(dplyr)
require(magrittr)
require(readr)
require(tibble)
require(shiny)
require(shinyWidgets)

# Read posteriors
# post_df <- read_rds("Bayesian_Hammer/hammer_posteriors.RDS")
post_df <- read_rds("hammer_posteriors.RDS")

# Prediction function
predict_posterior <- function(pst_df, sx = FALSE, gsurg = FALSE, csurg = FALSE, 
                              hyper = FALSE, apache = FALSE, fluid = FALSE, amb = FALSE,
                              los = FALSE)
{
    # Debug
    # sx <- male
    # gsurg <- general_surgery
    # csurg <- cardiac_surgery
    # hyper <- hyperglycemia
    # apache <- high_apache
    # fluid <- fluid_balance_5L
    # amb <- ambulation
    # los <- los_5
    
    pst_df <- hammer_pst
    
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
                        onStatus = "danger", 
                        offStatus = "success"
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
                    ),
             ),
             mainPanel(plotOutput("distPlot")))
)

# Define server logic required to draw a histogram
server <- function(input, output) {

}

# Run the application 
shinyApp(ui = ui, server = server)
