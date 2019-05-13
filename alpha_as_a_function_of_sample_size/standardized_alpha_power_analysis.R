#Standardized p-value
#This is untested alpha software provided with no guarantees - use at your own risk

#Good 1982: P SQRT(N/100)
#Enter the observed p-value and the sample size

p_stan <- function(p, N){
  p * sqrt(N/100)
}

#For example:
p_stan(p = 0.05, N = 100)
p_stan(p = 0.05, N = 200)

#Calculate the standardized alpha level
a_stan <- function(p, N){
  p/sqrt(N/100)
}

#Check it yields .05 when reversing a calculation above:
a_stan(p = 0.07071068, N = 200)

#Calculate standardized alpha levels as a function of the sample size: 
a_stan(p = 0.05, N = 200)
a_stan(p = 0.05, N = 10000)


#Power analysis for an independent t-test using a standardized alpha level

#The alpha level is based on N according the a_stan formula above.
#The desired power (e.g., 90%) is based on the effect size (e.g., d = 0.3), the standarized alpha, and N. 
#The standardized alpha is also dependent on N. 
#Because of 2 unkowns, we need to iteratively tweak the alpha and and N to get the desired power.

#Load pwr package
library(pwr)

#Set the value of the unstandardized alpha level (e.g., 0.05)
alpha_unstandardized <- 0.05

#We start using the unstandardized alpha level for the standardized alpha level, and then lower it. 
alpha_standardized <- alpha_unstandardized

#Specify how close the alpha should be to the samplesize adjusted p-value
#Lower the tolerance value if you don't get a solution
tolerance <- 0.0001

#Specify the smallest effect size of interest (or, perhaps, the expected effect size, although this is rarely known)
sesoi <- 0.3
#Specify the desired power (e.g., 0.9)
desired_power <- 0.9
#Run a power analysis once just to populate the pwr.output dataframe we will use below
pwr.output <- pwr.t.test(d = sesoi, sig.level = alpha_standardized, power = desired_power, type = "two.sample", alternative="two.sided")

#While the difference between the standardized p-value (based on the sample size in the power calculation) and the alpha level is too large, do the following:
#Lower the alpha a tiny bit (by 0.00001).
#Perform a power analysis for the sesoi and the desired power, with the current alpha
#We will get a sample size. If the alpha and standardized p-value are the same, we have the alpha level that, given the sample size, yields a standardized p-value that takes into account the sample size. 
while(abs(alpha_unstandardized / sqrt(pwr.output$n/100) - alpha_standardized) > tolerance){
  alpha_standardized <- alpha_standardized - 0.00001
  pwr.output <- pwr.t.test(d = sesoi, 
                           sig.level = alpha_standardized, 
                           power = desired_power, 
                           type = "two.sample", 
                           alternative = "two.sided")
}

#What is the sample size per group we should collect?
ceiling(pwr.output$n)
#What is the standardized alpha level, adjusted for the N we will collect, that we will be using for our test?
alpha_standardized
