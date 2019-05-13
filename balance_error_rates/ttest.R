#This code is designed to calculate optimal alpha levels for t-tests, authored by Joe Mudge (joe.mudge@unb.ca).
#The function used to calculate optimal alphas is: optab(n1=NULL,n2=NULL,d=NULL,T1T2cratio=1,HaHopratio=1,type = c("two.sample", "one.sample", "paired"),tails = c("two.tailed","one.tailed"))
#The arguments 'n1' and 'n2' are the samples sizes of each group
#The argument 'd' is the 'Cohen's d' standardized critical effect size. Cohen's d = difference between group means/pooled within group standard deviation
#The argument 'T1T2cratio' is the cost ratio of Type I errors relative to Type II errors. T1T2cratio is set at 1 as a default, making Type I and Type II errors equally serious.
#The argument 'HaHopratio' is the prior probability of the alternate hypothesis relative to the prior probability of the null hypothesis.  HaHopratio is set at 1 as a default, to not weight alpha and beta by their prior probabilities (assuming they are unknown).
#The argument 'type' is the type of t-test being undertaken and must be "two.sample", "one.sample" or "paired". If ignored, "two.sample" is the default.
#The argument 'tails'is the number of tails being examined and must be either "two.tailed" or "one.tailed".  If ignored, "two.tailed" is the default.
#This code is partially based on code modified from the R package 'pwr'(Champely 2009)

beta.t.test<-function (n1 = NULL, n2 = NULL, d = NULL, sig.level = 0.05, type = c("two.sample", "one.sample", "paired"),
                       tails = c("two.tailed","one.tailed")){

  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > sig.level | sig.level > 1))
    stop(sQuote("sig.level"), " must be numeric in [0, 1]")
  if (!is.null(n1) && n1 < 2)
    stop("number of observations in the first group must be at least 2")
  type  <- match.arg(type)
  tails <- match.arg(tails)
  d     <- abs(d)
  tsample <- switch(type, one.sample = 1, two.sample = 2, paired = 1)
  tside   <- switch(tails, one.tailed = 1, two.tailed = 2)
  if (tside == 1) {
    p.body <- quote({
      nu <- switch(type, one.sample = n1-1, two.sample = n1 + n2 - 2, paired = n1-1)
      pt(qt(sig.level/tside, nu, lower = FALSE), nu, ncp = d * switch(type, one.sample = sqrt(n1), two.sample = (1/sqrt(1/n1 + 1/n2)), paired = sqrt(n1)), lower = FALSE)
    })
  }
  if (tside == 2) {
    p.body <- quote({
      nu <- switch(type, one.sample = n1-1, two.sample = n1 + n2 - 2, paired = n1-1)
      qu <- qt(sig.level/tside, nu, lower = FALSE)
      pt(qu, nu, ncp = d * switch(type, one.sample = sqrt(n1), two.sample = (1/sqrt(1/n1 + 1/n2)), paired = sqrt(n1)), lower = FALSE) +
        pt(-qu, nu, ncp = d * switch(type, one.sample = sqrt(n1), two.sample = (1/sqrt(1/n1 + 1/n2)), paired = sqrt(n1)), lower = TRUE)
    })
  }
  1-eval(p.body)
}

w.average.error<-function (alpha=NULL,
                           n1=NULL,
                           n2=NULL,
                           d=NULL,
                           T1T2cratio=1,
                           HaHopratio=1,
                           type = c("two.sample", "one.sample", "paired"),
                           tails = c("two.tailed","one.tailed")){
  ((alpha*T1T2cratio+HaHopratio*(beta.t.test(n1=n1,n2=n2,d=d,sig.level=alpha,type=type,tails=tails))))/(HaHopratio+T1T2cratio)
}

# Maybe use this instead of w.average.error() if there is no cost
w.average.error.rep <-function (alpha=NULL,
                           n1=NULL,
                           n2=NULL,
                           d=NULL,
                           T1T2cratio=1,
                           HaHopratio=1,
                           priorP  =.005,
                           priorH0.ori = .5,
                           type = c("two.sample", "one.sample", "paired"),
                           tails = c("two.tailed","one.tailed")){
  postH0.ori <- ( 1 + ( (ifelse( priorP < (1/exp(1)), -exp(1)*priorP*log(priorH0.ori), 1) * priorH0.ori) / (1-priorH0.ori) )^-1 )^-1
  postH1.ori <- 1- priorH0.ori

  (alpha*postH0.ori + (beta.t.test(n1=n1,n2=n2,d=d,sig.level=alpha,type=type,tails=tails))*postH1.ori)/2
}

optimize.average.error<-function (f, interval, ..., lower = min(interval), upper = max(interval), maximum = FALSE, tol = .Machine$double.eps^0.25){
  if (maximum) {
    val <- optimize(f=f, interval = interval, ..., lower=lower, upper=upper, maximum=maximum, tol=tol)
    return(list(maximum = val$minimum, objective = val$objective))
  }
  else {
    val <- optimize(f=f, interval = interval, ..., lower=lower, upper=upper, maximum=maximum, tol=tol)
    return(val)
  }
}

min.average.error <-function (n1=NULL,
                              n2=NULL,
                              d=NULL,
                              T1T2cratio=1,
                              HaHopratio=1,
                              type = c("two.sample", "one.sample", "paired"),
                              tails = c("two.tailed","one.tailed")){


  optimize.average.error(w.average.error,
                         c(0,1),
                         tol=0.0000000000001,
                         n1=n1,
                         n2=n2,
                         d=d,
                         T1T2cratio=T1T2cratio,
                         HaHopratio=HaHopratio,
                         type=type,
                         tails=tails)$minimum
  }

optimize.alpha <- function (f, interval,..., lower = min(interval), upper = max(interval), maximum = FALSE, tol = .Machine$double.eps^0.25){
  if (maximum) {
    val <-  optimize(f=f, interval = interval, ..., lower=lower, upper=upper, maximum=maximum, tol=tol)
    return(list(maximum = val, objective = f(val, ...)))
  }
  else {
    val <-  optimize(f=f, interval = interval, ..., lower=lower, upper=upper, maximum=maximum, tol=tol)
    return(val)
  }
}

alpha<-function (n1=NULL,n2=NULL,d=NULL,T1T2cratio=1,HaHopratio=1,type = c("two.sample", "one.sample", "paired"),tails = c("two.tailed","one.tailed")){
  optimize.alpha(w.average.error,
                 c(0,1),
                 tol=0.000000000001,
                 n1=n1,
                 n2=n2,
                 d=d,
                 T1T2cratio=T1T2cratio,
                 HaHopratio=HaHopratio,
                 type=type,
                 tails=tails)$minimum
  }

beta <- function(n1=NULL,
                 n2=NULL,
                 d=NULL,
                 T1T2cratio=1,
                 HaHopratio=1,
                 type = c("two.sample", "one.sample", "paired"),
                 tails = c("two.tailed","one.tailed")){
  ((T1T2cratio+HaHopratio)*min.average.error(n1=n1,n2=n2,d=d,T1T2cratio=T1T2cratio,HaHopratio=HaHopratio,type=type,tails=tails)-T1T2cratio*alpha(n1=n1,n2=n2,d=d,T1T2cratio=T1T2cratio,HaHopratio=HaHopratio,type=type,tails=tails))/HaHopratio
  }



optab<-function (n1=NULL,n2=NULL,d=NULL,T1T2cratio=1,HaHopratio=1,type = c("two.sample", "one.sample", "paired"),tails = c("two.tailed","one.tailed")) {
  list("test type"=match.arg(type),
       "tails"=match.arg(tails),
       "output"=t(data.frame("sample size 1"=n1,
                             "sample size 2"=n2,
                             "Cohen's d effect size"=d,
                             "Type I/II error cost ratio"=T1T2cratio,
                             "Ha/Ho prior probability ratio"=HaHopratio,
                             "average of alpha and beta"=(alpha(n1=n1,n2=n2,d=d,T1T2cratio=T1T2cratio,HaHopratio=HaHopratio,type=type,tails=tails)+HaHopratio*beta(n1=n1,n2=n2,d=d,T1T2cratio=T1T2cratio,HaHopratio=HaHopratio,type=type,tails=tails))/(1+HaHopratio),
                             "cost-weighted average of alpha and beta"=min.average.error(n1=n1,n2=n2,d=d,T1T2cratio=T1T2cratio,HaHopratio=HaHopratio,type=type,tails=tails),
                             "optimal alpha"=alpha(n1=n1,n2=n2,d=d,T1T2cratio=T1T2cratio,HaHopratio=HaHopratio,type=type,tails=tails),
                             "optimal beta"=beta(n1=n1,n2=n2,d=d,T1T2cratio=T1T2cratio,HaHopratio=HaHopratio,type=type,tails=tails),row.names="values"))
  )}


#This code is designed to calculate optimal alpha levels for t-tests, authored by Joe Mudge (joe.mudge@unb.ca).
#The function used to calculate optimal alphas is: optab(n1=NULL,n2=NULL,d=NULL,T1T2cratio=1,HaHopratio=1,type = c("two.sample", "one.sample", "paired"),tails = c("two.tailed","one.tailed"))
#The arguments 'n1' and 'n2' are the samples sizes of each group
#The argument 'd' is the 'Cohen's d' standardized critical effect size. Cohen's d = difference between group means/pooled within group standard deviation
#The argument 'T1T2cratio' is the cost ratio of Type I errors relative to Type II errors. T1T2cratio is set at 1 as a default, making Type I and Type II errors equally serious.
#The argument 'HaHopratio' is the prior probability of the alternate hypothesis relative to the prior probability of the null hypothesis.  HaHopratio is set at 1 as a default, to not weight alpha and beta by their prior probabilities (assuming they are unknown).
#The argument 'type' is the type of t-test being undertaken and must be "two.sample", "one.sample" or "paired". If ignored, "two.sample" is the default.
#The argument 'tails'is the number of tails being examined and must be either "two.tailed" or "one.tailed".  If ignored, "two.tailed" is the default.
#This code is partially based on code modified from the R package 'pwr'(Champely 2009)
