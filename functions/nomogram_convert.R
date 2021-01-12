nomogram_convert <- function(input_variable,
                             input_system, output_system,
                             log = FALSE)
{
  # Debug
  # input_system <- points_system_input
  # output_system <- points_system_output
  
  # Create input model
  input_model <- lm(pix ~ input, input_system) %>% 
    coefficients()
  
  # Calculate input pixels
  input_pixels <- input_model[1] + (input_variable * input_model[2])
  

  if(log == TRUE)
  {
    # Create output model
    output_model <- lm(log(output) ~ pix, output_system) %>% 
      coefficients()
    
    # Calculate output variable
    return(inverse_logit(output_model[1] + (input_pixels * output_model[2])))
  }else
  {
    # Create output model
    output_model <- lm(output ~ pix, output_system) %>% 
      coefficients()
    
    # Calculate output variable
    return(output_model[1] + (input_pixels * output_model[2]))
  }
}