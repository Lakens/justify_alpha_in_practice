library(pwr)
library(TOSTER)
options(scipen = 99)
# Instead of using Mudge et al (2012) functions, 
# we can just directly optimize any existing power function.
# 
# This function calculates the Type 2 error rate (1-power, y) and adds it to the
# Type 1 error rate (alpha, x). 
# 
# It then minimizes the combined error rate. 
# 
# The code below reproduces Mudge et al., code: 
# optab(n1=8, n2=8, d=0.5, T1T2cratio=1, HaHopratio=1, type = "two.sample", tails = "two.tailed")
#
# We create a function for the power analysis where we can input the alpha level (x).
# Then we calculate the Type 2 error (1-power), stored as y
# We add x and y and minimize the combined error.

f <- function(x) {
  y <- 1 - pwr.t.test(d=0.5,
                      n=8, 
                      sig.level = x,
                      type="two.sample",
                      alternative="two.sided")$power
  print(c(x, y, (x+y)/2)) #optional: print alpha, beta, and combined error
  x+y
}
# We store the result (res) of the optimize function.
# We minimize the outcome of function  f, over the input (alpha) range 0 to 1.
# We accept a tollerance of 0.0001 (can be adjusted if results do not converge)
res <- optimize(f, c(0, 1), tol = 0.0001)
# We can request the minimum value:
res$minimum
# [1] 0.2963138
# The 'objective', the combined error rate (x+y) ends up being 0.8004419
res$objective

#increase x a little, total error should increase
y <- 1 - pwr.t.test(d=0.5,
                    n=8, 
                    sig.level = 0.35,
                    type="two.sample",
                    alternative="two.sided")$power
0.35+y
# [1] 0.8027443 instead of 0.8004419

# What would it be if we would use a traditional 0.05 alpha?
y <- 1 - pwr.t.test(d=0.5,
                    n=8, 
                    sig.level = 0.05,
                    type="two.sample",
                    alternative="two.sided")$power
0.05+y
# Combined error rate is [1] 0.8959525  - 0.095 higher.

# What would it be if we would use a 0.005 alpha level?
y <- 1 - pwr.t.test(d=0.5,
                    n=8, 
                    sig.level = 0.005,
                    type="two.sample",
                    alternative="two.sided")$power
0.005+y
# Combined error rate is [1] 0.9785644  - 0.178 higher.

# We can use the same code for any power function.
# For example, for power analysis for equivalence tests
#
f <- function(x) {
  y <- 1 - powerTOSTtwo(alpha = x, 
                        N = 121, 
                        low_eqbound_d=-0.3, 
                        high_eqbound_d=0.3)
  print(c(x, y, x+y)) #optional: print alpha, beta, and combined error
  x+y
}
res <- optimize(f, c(0, 1), tol = 0.0001)
# The alpha level is [1] 0.1943336
res$minimum

# But we might prefer equal error rates, not optimized error rates.
# 
f <- function(x) {
  y <- 1 - pwr.t.test(d=0.5,
                      n=8, 
                      sig.level = x,
                      type="two.sample",
                      alternative="two.sided")$power
  print(c(x, y, x+y)) #optional: print alpha, beta, and combined error
  max(x - y, y - x)
}
# We store the result (res) of the optimize function.
# We minimize the outcome of function  f, over the input (alpha) range 0 to 1.
# We accept a tollerance of 0.0001 (can be adjusted if results do not converge)
res <- optimize(f, c(0, 1), lower = 0, tol = 0.0001)
# We can request the minimum value:
res$minimum
# [1] 0.2963138
# The 'objective', the combined error rate (x+y) ends up being 0.8004419
res$objective



pwr.t.test(d=0.5,
           n=8, 
           sig.level = 0.4045606,
           type="two.sample",
           alternative="two.sided")$power

pwr.t.test(d=0.5,
           n=8, 
           sig.level = 0.05,
           type="two.sample",
           alternative="two.sided")$power


# The cost ratio of Type I errors relative to Type II errors can be changed
# By default Type 1 and Type 2 errors are weighed equally. 
# Multiplying the Type 1 error rate by the relative cost will change the optimium.
# In this case, we weight Type 1 errors 4 times as much as Type 2 errors.
f <- function(x) {
  y <- 1 - pwr.t.test(d=0.5,
                      n=50, 
                      sig.level = x,
                      type="two.sample",
                      alternative="two.sided")$power
  print(c(x, y, x+y)) #optional: print alpha, beta, and combined error
  4*x+y
}
res <- optimize(f, c(0, 1), tol = 0.0000001)
# The Type 1 error rate is now [1] 0.03875288 compared to 0.1301174 with equal weighting
res$minimum
res$objective


#In the example below the prior says H1 is 6 times as likely as H0
f <- function(x) {
  y <- 1 - pwr.t.test(d=0.5,
                      n=50, 
                      sig.level = x,
                      type="two.sample",
                      alternative="two.sided")$power
  print(c(x, y, x+y)) #optional: print alpha, beta, and combined error
  (x+6*y)/7
}
res <- optimize(f, c(0, 1), tol = 0.0000001)
# The Type 1 error rate is now [1] 0.03875288 compared to 0.1301174 with equal weighting
res$minimum
res$objective


#Flexible example where we can set cost and priors.
f <- function(x, costT1T2, prior_H1H0) {
  y <- 1 - pwr.t.test(d=0.5,
                      n=50, 
                      sig.level = x,
                      type="two.sample",
                      alternative="two.sided")$power
  print(c(x, y, x+y)) #optional: print alpha, beta, and combined error
  (costT1T2*x+prior_H1H0*y)/(prior_H1H0+1)
}

res <- optimize(f, c(0, 1), tol = 0.0000001, costT1T2 = 2, prior_H1H0 = 6)



# Flexible example where we can set cost and priors for one-way ANOVA.
f <- function(x, costT1T2, prior_H1H0, error) {
  y <- 1 - pwr.anova.test(f = 0.28,
                          k = 4,
                          n = 20,
                          sig.level = x)$power
  print(c(x, y, x+y)) #optional: print alpha, beta, and combined error
  if(error == "balance"){
    res <- max((costT1T2*x - prior_H1H0*y)/(prior_H1H0+1), (prior_H1H0*y - costT1T2*x)/(prior_H1H0+1))
  } else if (error == "minimal"){
    res <- (costT1T2*x + prior_H1H0*y)/(prior_H1H0+1)
  }
}

res <- optimize(f, c(0, 1), tol = 0.00001, costT1T2 = 1, prior_H1H0 = 1, error = "minimal")




#####TESTTTTT
#####
#####

power_function <- "pwr.anova.test(f = 0.28,
                        k = 4,
                        n = 20,
                        sig.level = x)$power"

# Flexible example where we can set cost and priors for one-way ANOVA.
f <- function(x, power_function, costT1T2, prior_H1H0, error) {
  y <- 1 - eval(parse(text=paste(func)))
  print(c(x, y, x+y)) #optional: print alpha, beta, and combined error
  if(error == "balance"){
    res <- max((costT1T2*x - prior_H1H0*y)/(prior_H1H0+1), (prior_H1H0*y - costT1T2*x)/(prior_H1H0+1))
  } else if (error == "minimal"){
    res <- (costT1T2*x + prior_H1H0*y)/(prior_H1H0+1)
  }
}

res <- optimize(f, 
                c(0, 1), 
                tol = 0.00001, 
                power_function = power_function, 
                costT1T2 = 1, 
                prior_H1H0 = 1, 
                error = "minimal")




power_function <- "pwr.anova.test(f = 0.28,
                        k = 4,
                        n = 20,
                        sig.level = x)$power"

# Flexible example where we can set cost and priors for one-way ANOVA.
f <- function(x, power_function, costT1T2, prior_H1H0, error) {
  y <- 1 - eval(parse(text=paste(power_function)))
  print(c(x, y, x+y)) #optional: print alpha, beta, and combined error
  if(error == "balance"){
    res <- max((costT1T2*x - prior_H1H0*y)/(prior_H1H0+1), (prior_H1H0*y - costT1T2*x)/(prior_H1H0+1))
  } else if (error == "minimal"){
    res <- (costT1T2*x + prior_H1H0*y)/(prior_H1H0+1)
  }
}

res <- optimize(f, 
                c(0, 1), 
                tol = 0.00001, 
                power_function = "pwr.anova.test(f = 0.28, k = 4, n = 20, sig.level = x)$power", 
                costT1T2 = 1, 
                prior_H1H0 = 1, 
                error = "minimal")

#CORRECT FUNCTION
optimal_alpha <- function(power_function, costT1T2, prior_H1H0, error) {
  f = function(x, power_function, costT1T2, prior_H1H0, error) {
    y <- 1 - eval(parse(text=paste(power_function)))
    print(c(x, y, x+y)) #optional: print alpha, beta, and combined error
    if(error == "balance"){
      max((costT1T2*x - prior_H1H0*y)/(prior_H1H0+1), (prior_H1H0*y - costT1T2*x)/(prior_H1H0+1))
    } else if (error == "minimal"){
      2*(costT1T2*x + prior_H1H0*y)/(prior_H1H0+1)
    }
  }
  res <- optimize(f, 
                  c(0, 1), 
                  tol = 0.00001,
                  power_function = power_function,
                  costT1T2 = costT1T2, 
                  prior_H1H0 = prior_H1H0, 
                  error = error)
  if(error == "balance"){
    beta <- res$minimum - res$objective
  } else if (error == "minimal"){
    beta <- res$objective - res$minimum
  }
  
  invisible(list(alpha = res$minimum,
                 beta = beta,
                 tot = res$objective
  )
  )
}

x <- optimal_alpha(power_function = "pwr.anova.test(f = 0.28, k = 4, n = 20, sig.level = x)$power", 
              costT1T2 = 1, 
              prior_H1H0 = 1, 
              error = "minimal")

x <- optimal_alpha(power_function = "pwr.t.test(d=0.5, n=8, sig.level = x, type='two.sample', alternative='two.sided')$power", 
                   costT1T2 = 1, 
                   prior_H1H0 = 1, 
                   error = "balance")
x
# check
pwr.t.test(d=0.5, n=8, sig.level = 0.2963138, type='two.sample', alternative='two.sided')$power


x <- optimal_alpha(power_function = "powerTOSTtwo(alpha=x, N=108, low_eqbound_d=-0.4, high_eqbound_d=0.4)", 
                   costT1T2 = 1, 
                   prior_H1H0 = 1, 
                   error = "balance")
x$alpha + x$beta
# [1] 0.1966935 (minimal) vs [1] 0.1975292 (balance) hardly matters

