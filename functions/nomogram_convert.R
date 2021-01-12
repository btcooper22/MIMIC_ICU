nomogram_convert <- function(input_variable,
                             input_system, output_system)
{
  # Debug
  # input_system <- age_score_system_input
  # output_system <- age_score_system_output
  
  # Create input model
  input_model <- lm(pix ~ input, input_system) %>% 
    coefficients()
  
  # Calculate input pixels
  input_pixels <- input_model[1] + (input_variable * input_model[2])
  
  # Create output model
  output_model <- lm(output ~ pix, output_system) %>% 
    coefficients()
  
  # Calculate output variable
  output_model[1] + (input_pixels * output_model[2])
}