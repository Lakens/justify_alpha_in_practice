#Reproduce lowest row Table 1 Mudge et al 2012----

require("pwr")

optab(n1=100, n2=100, d=0.5, T1T2cratio=1, HaHopratio=1, type = "two.sample",
      tails = "two.tailed")

1-pwr.t.test(n = 5.5, d = 0.5, sig.level = 0.3228697, 
           type = "two.sample", alternative = "two.sided")$power

#Result is not equal to example, because pwr package only does equal sample sizes.

#Gpower gives power of 0.4373076 instead of 0.4639465

1-0.4373076 #which matches lowest row. 

#omega
(0.3228697 + (1-0.4373076))/2


#Random Example----

n<-100
d<-0.4
normal_alpha <- 0.05

results<-optab(n1=n, n2=n, d=d, T1T2cratio=1, HaHopratio=1, type = "two.sample",
      tails = "two.tailed")
optimal_alpha<-results$output[8]

optimal_beta<-1-pwr.t.test(n = n, d = d, sig.level = optimal_alpha, 
             type = "two.sample", alternative = "two.sided")$power

normal_beta<-1-pwr.t.test(n = n, d = d, sig.level = normal_alpha, 
                           type = "two.sample", alternative = "two.sided")$power


optimal_omega<-(optimal_beta+optimal_alpha)/2
normal_omega<-(normal_beta+normal_alpha)/2

optimal_alpha
optimal_beta
optimal_omega

normal_alpha
normal_beta
normal_omega

#error_gain
normal_omega-optimal_omega





get_alpha_beta<-function(n,d,normal_alpha){
  results<-optab(n1=n, n2=n, d=d, T1T2cratio=1, HaHopratio=1, type = "two.sample",
                 tails = "two.tailed")
  optimal_alpha<-results$output[8]
  
  optimal_beta<-1-pwr.t.test(n = n, d = d, sig.level = optimal_alpha, 
                             type = "two.sample", alternative = "two.sided")$power
  
  normal_beta<-1-pwr.t.test(n = n, d = d, sig.level = normal_alpha, 
                            type = "two.sample", alternative = "two.sided")$power
  
  
  optimal_omega<-(optimal_beta+optimal_alpha)/2
  normal_omega<-(normal_beta+normal_alpha)/2

  #error_gain
  normal_omega-optimal_omega
}


#The error gain is only a function of the power of the test, and we can plot them against each other
#Always goes to alpha/2 (since with huge power, the Type 2 error rate is tiny, and we adjust the alpha)
#But at smaller n, benefit comes from increasing alpha, which has a positive effect on the power
library(pwr)
source('ttest.R')

options(scipen=20) #disable scientific notation for numbers
n<-3000
d<-0.1
outcome<-numeric(n)
for(i in 10:n){
  outcome[i]<-get_alpha_beta(d=d, n=i, normal_alpha=0.05)
}

plot(outcome, type="l")


Power <- (function(D, n, alpha, side)
{
  ncp <- D*(n*n/(n+n))^0.5 #formula to calculate t from d from Dunlap, Cortina, Vaslow, & Burke, 1996, Appendix B
  t <- qt(1-(alpha/side),df=(n*2)-2)
  1-(pt(t,df=n*2-2,ncp=ncp)-pt(-t,df=n*2-2,ncp=ncp))
}
)
par(new=TRUE)
curve(Power(D=d, n=x, alpha=0.05, side=2), 10, n, type="l", lty=1, lwd=2, ylim=c(0,1), xlim=c(0,n), col="red", xlab='', ylab='', axes=FALSE)


plot(outcome, type="l")