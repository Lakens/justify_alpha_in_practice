#' Justify your alpha level by minimizing or balancing Type 1 and Type 2 error rates.
#' @param power_function Function that outputs the power, calculated with an analytic function.
#' @param costT1T2 Relative cost of Type 1 errors vs. Type 2 errors.
#' @param priorH1H0 How much more likely a-priori is H1 than H0?
#' @param error Either "minimal" to minimize error rates, or "balance" to balance error rates.
#' @param verbose Print each iteration of the optimization function if TRUE. Defaults to FALSE.
#' @param printplot Print a plot to illustrate the alpha level calculation.
#' @return
#' alpha = alpha or Type 1 error that minimizes or balances combined error rates
#' beta = beta or Type 2 error that minimizes or balances combined error rates
#' errorate = weighted combined error rate
#' objective = value that is the result of the minimization, either 0 (for balance) or the combined weighted error rates
#'
#' @examples
#' ## Optimize power for a independent t-test, smallest effect of interest
#' ## d = 0.5, 100 participants per condition
#' res <- optimal_alpha(power_function = "pwr::pwr.t.test(d = 0.5, n = 100,
#' sig.level = x, type = 'two.sample', alternative = 'two.sided')$power")
#' res$alpha
#' res$beta
#' res$errorate
#' @section References:
#' too be added
#' @importFrom stats optimize
#' @export
#'
optimal_alpha <- function(power_function, costT1T2 = 1, priorH1H0 = 1, error = "minimal", verbose = FALSE, printplot = FALSE) {
  #Define the function to be minimized
  f = function(x, power_function, costT1T2 = 1, priorH1H0 = 1, error = "minimal") {
    y <- 1 - eval(parse(text=paste(power_function)))
    if(verbose == TRUE){
      print(c(x, y, x+y)) #optional: print alpha, beta, and objective
    }
    if(error == "balance"){
      max((costT1T2 * x - priorH1H0 * y)/(priorH1H0 + costT1T2), (priorH1H0 * y - costT1T2 * x)/(priorH1H0 + costT1T2))
    } else if (error == "minimal"){
      (costT1T2 * x + priorH1H0 * y)/(priorH1H0 + costT1T2)
    }
  }
  #Run optimize to find the minimum
  res <- stats::optimize(f,
                         c(0, 1),
                         tol = 0.001,
                         power_function = power_function,
                         costT1T2 = costT1T2,
                         priorH1H0 = priorH1H0,
                         error = error)
  if(error == "balance"){
    beta <- res$minimum - res$objective
  } else if (error == "minimal"){
    beta <- res$objective - res$minimum
  }
  
  
  
  
  #Add plot
  
  alpha_level <- 0
  alpha_list <- numeric(9999)
  beta_list <- numeric(9999)
  w_list <- numeric(9999)
  w_c_list <- numeric(9999)
  for(i in 1:9999) {
    alpha_level <- alpha_level + 0.0001
    alpha_list[i] <- alpha_level
    x <- alpha_level
    beta_list[i] <- 1 - eval(parse(text=paste(power_function)))
    w_list[i] <- (alpha_level + beta_list[i]) / 2
    w_c_list[i] <- (costT1T2 * alpha_level + priorH1H0 * beta_list[i]) / (costT1T2 + priorH1H0)
  }
  
  x <- res$minimum
  
  # Create dataframe for plotting
  plot_data <- data.frame(alpha_list, beta_list, w_list, w_c_list)
  
  w_c_alpha_plot <- ggplot2::ggplot(ggplot2::aes(x=alpha_list, y=w_c_list), data=plot_data,) +
    ggplot2::geom_line(size = 1.3) +
    ggplot2::geom_point(ggplot2::aes(x = res$minimum, y = (costT1T2 * res$minimum + priorH1H0 * (1 - eval(parse(text=paste(power_function))))) / (priorH1H0 + costT1T2)), color="red", size = 3) +
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::scale_x_continuous("alpha", seq(0,1,0.1)) +
    ggplot2::scale_y_continuous("weighted combined error rate", seq(0,1,0.1), limits = c(0,1))
  
  if(printplot == TRUE){
    print(w_c_alpha_plot)
  }
  
  #Store results
  alpha = res$minimum
  beta = 1 - eval(parse(text=paste(power_function)))
  
  invisible(list(alpha = res$minimum,
                 beta = 1 - eval(parse(text=paste(power_function))),
                 errorrate = (costT1T2 * alpha + priorH1H0 * beta) / (costT1T2 + priorH1H0), 
                 objective = res$objective,
                 plot_data = plot_data,
                 plot = w_c_alpha_plot
  )
  )
}


alpha_sample_solve <- function(i, power_function, errorgoal, costT1T2, priorH1H0, error){
  
  res <- optimal_alpha(power_function = paste(stringr::str_replace(power_function, "n = i", paste("n =", i))), costT1T2 = costT1T2, priorH1H0 = priorH1H0, error = error)
  (errorgoal - res$errorrate)^2
}

#' Justify your alpha level by minimizing or balancing Type 1 and Type 2 error rates.
#' @param power_function Function that outputs the power, calculated with an analytic function.
#' @param errorgoal Desired weighted combined error rate
#' @param costT1T2 Relative cost of Type 1 errors vs. Type 2 errors.
#' @param priorH1H0 How much more likely a-priori is H1 than H0?
#' @param error Either "minimal" to minimize error rates, or "balance" to balance error rates.
#' @return
#' alpha = alpha or Type 1 error that minimizes or balances combined error rates
#' beta = beta or Type 2 error that minimizes or balances combined error rates
#' errorrate = weighted combined error rate
#' objective = value that is the result of the minimization, either 0 (for balance) or the combined weighted error rates
#' samplesize = the desired samplesize
#'
#' @examples
#' ## Optimize power for a independent t-test, smallest effect of interest
#' ## d = 0.5, desired weighted combined error rate = 5%
#' res <- SampleSizeError(power_function = "pwr::pwr.t.test(d = 0.5, n = i, sig.level = x,
#' type = 'two.sample', alternative = 'two.sided')$power",errorgoal = 0.05)
#' res$alpha
#' res$beta
#' res$errorrate
#' res$samplesize
#' @section References:
#' too be added
#' @importFrom stats optimize
#' @export
#'
SampleSizeError <- function(power_function, errorgoal = 0.05, costT1T2 = 1, priorH1H0 = 1, error = "minimal") {

  
  samplesize<- optim(20, alpha_sample_solve, lower = 0, upper = Inf, method = "L-BFGS-B",
                  power_function = power_function, errorgoal = errorgoal, costT1T2 = costT1T2, priorH1H0 = priorH1H0, error = error)$par
  samplesize <- ceiling(samplesize)
  res <- optimal_alpha(power_function = paste(stringr::str_replace(power_function, "n = i", paste("n =", samplesize))), costT1T2 = costT1T2, priorH1H0 = priorH1H0, error = error)
  invisible(list(alpha = res$alpha,
                 beta = res$beta,
                 errorrate = (costT1T2 * res$alpha + priorH1H0 * res$beta) / (costT1T2 + priorH1H0), 
                 objective = res$objective,
                 samplesize = samplesize
  )
  )
}



