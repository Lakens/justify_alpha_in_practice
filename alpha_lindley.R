
### t-test ###

###Compute Bayes Factor for t-statistic and degrees of freedom###

bf_t.test <- function(t, n1, n2 = 0, Cauchy = T){
  if (n2 == 0) {
    v <- n1 -1 
    N <- n1
  } else {
    N <- n1*n2/(n1+n2)
    v <- n1 + n2 -2
  }
  if (Cauchy == F) {
    ml1 <- 1/sqrt(1+N)*(1+t^2/((1+N)*v))^(-(v+1)/2)
    ml0 <- (1+t^2/v)^(-(v+1)/2)
    return(ml1/ml0)
  } else {
    return(exp(BayesFactor::ttest.tstat(t, n1, n2)$bf))
  }
}

alpha_t.test_solve <- function(x, n1, n2, evidence, Cauchy){
  (evidence - bf_t.test(x, n1, n2, Cauchy))^2
}

#' Justify your alpha level by avoiding the Lindley paradox or aiming for moderate or strong evidence when using a t-test.
#' @param evidence Desired level of evidence. "Lindley" to avoid the Lindley Paradox, "moderate" to achieve moderate evidence and "strong" to achieve strong evidence. 
#' Users that are more familiar with Bayesian statistics can also directly enter their desired Bayes factor. 
#' @param n1 Sample size in Group 1.
#' @param n2 Sample size in Group 2. Leave blank for a one-sample or paired-sample t.test.
#' @param Cauchy If true the Cauchy prior with scale 0.707 is used if false, the Unit information Bayes factor is used.
#' @return alpha level required for a two-sample t-test.
#' @examples
#' ## Avoid the Lindley paradox for a two sample t-test with 300 participants per condition
#' ttestEvidence("lindley", 300, 300)
#' @section References:
#' too be added
#' @importFrom stats optim pf pt
#' @export
#'
ttestEvidence <- function(evidence, n1, n2 = 0, Cauchy = T) {
  if (evidence == "lindley"){
    evidence = 1
  }
  if (evidence == "moderate"){
    evidence = 3
  }
  if (evidence == "strong"){
    evidence = 10
  }
  
  if (n2 != 0){
    df = n1 + n2 -2
  }
  else  {
    df = n1 -1
  }
  
  crit_t <- optim(1.96, alpha_t.test_solve, lower = 0, upper = Inf, method = "L-BFGS-B",
                  n1 = n1, n2 = n2, evidence = evidence, Cauchy = Cauchy)$par
  alpha <- (1 - pt(crit_t, df))*2
  
  return(alpha)
}

### anova ###

bf_bic <- function(F, df1, df2, repeated=FALSE, report.as="BF10") {
  if (repeated==FALSE) {
    N = df1+df2+1
  } 
  else {
    N = df1+df2
  }
  
  bf = sqrt(N^df1*(1+F*df1/df2)^(-1*N))
  
  if (report.as=="BF01"){
    return(c(B01=bf))
  }
  else {
    return(c(B10=1/bf))
  }
}

alpha_f.test_solve <- function(x, df1, df2, evidence, paired){
  (evidence - bf_bic(x, df1, df2, paired))^2
}

#' Justify your alpha level by avoiding the Lindley paradox or aiming for moderate or strong evidence when using a t-test.
#' @param evidence Desired level of evidence. "Lindley" to avoid the Lindley Paradox, "moderate" to achieve moderate evidence and "strong" to achieve strong evidence. 
#' Users that are more familiar with Bayesian statistics can also directly enter their desired Bayes factor. 
#' @param df1 Numerator degrees of freedom.
#' @param df2 Denominator degrees of freedom.
#' @param paired If true a within subjects design is assumed.
#' @return alpha level required for a two-sample t-test.
#' @examples
#' ## Avoid the Lindley paradox for an anova with 1 numerator and 248 denominator degrees of freedom.
#' ftestEvidence("lindley", 1, 248)
#' @section References:
#' too be added
#' @export
#' @importFrom stats optim pf pt
#'
ftestEvidence <- function(evidence, df1, df2, paired = FALSE){
  
  if (evidence == "lindley"){
    evidence = 1
  }
  if (evidence == "moderate"){
    evidence = 3
  }
  if (evidence == "strong"){
    evidence = 10
  }
  
  crit_f <- optim(5, alpha_f.test_solve, lower = 0, upper = Inf, method = "L-BFGS-B", control = list(maxit = 100000), 
                  df1 = df1, df2 = df2, evidence = evidence, paired = paired)$par
  alpha <- (1 - pf(crit_f, df1, df2))
  
  return(alpha)
}
