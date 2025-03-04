#' Calibrate
#' 
#' General calibration function for malariasimulation. 
#'
#' @param parameters Malariasimulation model parameter list
#' @param target Values of target variable to calibrate to.
#' @param summary_function A function that take the raw model output as an 
#` argument and produces a vector of the same length as the target variable.
#' @param eq_prevalence Estimate of baseline prevalence (to inform EIR initialisation)
#' @param eq_ft  Estimate of baseline treatment coverage (to inform EIR initialisation)
#' @param human_population Vector of increasing population sizes. 
#' @param eir_limits Limits of EIR
#' @param max_attempts Maximum number of calibration attempts (model runs)
#'
#' @return EIR estimate
#' @export
calibrate <- function(parameters, target, summary_function, eq_prevalence, 
                      eq_ft = 0, human_population = c(1000, 10000, 100000), 
                      eir_limits = c(0, 1500), max_attempts = 10,
                      use_pfpr_6_59_mo = FALSE){
  
  eir <- rep(0, 2)
  objective <- rep(NA, 2)
  parameters$human_population <- human_population[1]
  
  # Estimate starting guess EIR with equilibrium solution
  message("Initialising EIR\n")
  eir[1] <- get_eq_eir(
    target_pfpr = eq_prevalence,
    ft = eq_ft,
    pfpr_6_59_on = use_pfpr_6_59_mo
  )
  if(eir[1] < eir_limits[1]){
    eir[1] <- eir_limits[1] + 1
  }
  if(eir[1] > eir_limits[2]){
    eir[1] <- eir_limits[2] - 1
  }
  
  # Ensure the starting EIR is valid and obtain objective evaluation
  min_eir <- 0
  attempts <- 0
  message("Slice sampling EIR, side 1 \n")
  while(is.na(objective[1])){
    message("Attempt ", attempts + 1, " of ", max_attempts)
    
    objective[1] <- check_fit(
      proposed_eir = eir[1],
      parameters = parameters,
      target = target,
      summary_function = summary_function
    )
    
    print_update(eir, objective)
    attempts <- attempts + 1
    
    # In the case that our starting EIR leads to unwanted elimination:
    ## First approach is to increase the human population if available.
    ## Second approach is to increase EIR and update our lower bound on EIR.
    if(is.na(objective[1])){
      if(parameters$human_population < max(human_population)){
        parameters$human_population <- human_population[which(human_population == parameters$human_population) + 1]
        message("Increasing human population due to elimination")
        message("Running with new population size: ", parameters$human_population, "\n")
      } else {
        parameters$human_population <- human_population[1]
        min_eir <- eir[1]
        eir[1] <- proposal(
          current_eir = eir[1],
          limits = eir_limits,
          direction = "increase",
          step = 2
        )
      }
    }
    if(attempts > max_attempts){
      stop("Failure due to max attempts reached before successful run")
    }
  }
  
  # Estimate the second EIR
  message("Slice sampling EIR, side 2 \n")
  direction <- ifelse(objective[1] > 0, "decrease", "increase")
  if(direction == "decrease"){
    eir_limits[1] <- min_eir
  }
  estimated_eir <- NA
  eir[2] <- eir[1]
  repeat{
    message("Attempt ", attempts + 1, " of ", max_attempts)
    # Propose new EIR
    eir[2] <- proposal(
      current_eir = eir[2],
      limits = eir_limits,
      direction = direction
    )
    # Evaluate objective
    objective[2] <- check_fit(
      proposed_eir = eir[2],
      parameters = parameters,
      target = target,
      summary_function = summary_function
    )
    print_update(eir, objective)
    # Stop or update
    if(!is.na(objective[2])){
      if(good_move(objective, direction)){
        estimated_eir <- linear_interpolate(eir, objective)
        message("Success")
        break
      } else{
        # If we haven't moved far enough, we can update eir1 and fit1  
        eir[1] <- eir[2]
        objective[1] <- objective[2]
      }
    } else {
      if(direction == "decrease"){
        eir_limits[1] <- eir[2]
        direction <- "increase"
      }
    }
    attempts <- attempts + 1
    if(attempts > max_attempts){
      estimated_eir <- mean(eir, na.rm = TRUE)
      message("Terminating as max attempts reached")
      break
    }
  }
  return(estimated_eir)
}
