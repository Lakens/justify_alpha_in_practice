# Reproduce lowest row Table 1 Mudge et al 2012----
# First run functions_mudge.R
require("pwr")

optab(n1=100, n2=100, d=0.5, T1T2cratio=4, HaHopratio=2, type = "two.sample",
      tails = "two.tailed")

res <- optimal_alpha(power_function = "pwr.t.test(d=0.5, n=100, sig.level = x, type='two.sample', alternative='two.sided')$power", 
                     error = "minimal", 
                     priorH1H0 = 2,
                     costT1T2 = 4)

res$alpha
res$beta

1-pwr.t.test(d=0.5, n=100, sig.level = res$alpha, type='two.sample', alternative='two.sided')$power

optab(n1=100, n2=100, d=0.5, T1T2cratio=1, HaHopratio=1, type = "two.sample",
      tails = "two.tailed")

res <- optimal_alpha(power_function = "pwr.t.test(d=0.5, n=100, sig.level = x, type='two.sample', alternative='two.sided')$power", 
                     error = "minimal", 
                     priorH1H0 = 1,
                     costT1T2 = 1)

res$alpha
res$beta

optab(n1=64, n2=64, d=0.5, T1T2cratio=1, HaHopratio=1, type = "two.sample",
      tails = "two.tailed")

res <- optimal_alpha(power_function = "pwr.t.test(d=0.5, n=64, sig.level = x, type='two.sample', alternative='two.sided')$power", 
                     error = "minimal", 
                     priorH1H0 = 1,
                     costT1T2 = 1)

res$alpha
res$beta

