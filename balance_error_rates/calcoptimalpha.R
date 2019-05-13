# Calculate optimal alpha for MANYLABS2

library(pwr)
source('ttest.R')

# Example: Structure and goal pursuit (Kay et al, 2014)
# Original: N=67  (M = 4.72, SD = 1.32; t(65) = 2.00, p = 0.05, d =  0.49, 95% CI [0.001, 0.973]).
# Replication: N=6506  (M = 5.51, SD = 1.39; t(6,498.63) = -0.94, p = 0.35, d = -.02, 95% CI [-0.07, 0.03]).

# p-value of original study
p <- .005

# Orginal effect size
d <- .49

# Original sample size
n1 <- 33
n2 <- 34

# Observed power based on original observed p-value & effect size
pwr <- pwr::pwr.t2n.test(n1 = n1, n2= n2, d = d, sig.level = p, power = NULL, alternative = "two.sided")

# Test and tail type
# DL Moved this up because needed below but gave error originally.
type  = "two.sample"
tails = "two.tailed"


# Beta given p-value & observed power
b_ori  <- 1-pwr$power
b_ori2 <- beta.t.test(n1=n1,n2=n2,d=d,sig.level=p,type=type,tails=tails)

# This should be equal
b_ori==b_ori2

# Prior probability for H0 in original study (H1 = 1-q)
q <- .5

# Posterior probability of H0 and H1 in original study
postH0 <- ( 1 + ( (ifelse( p < (1/exp(1)), -exp(1)*p*log(p), 1) * q) / (1-q) )^-1 )^-1
postH1 <- 1- postH0

# Combined error probability, given posterior probabilities of original H0 and H1
omega_ori <- (p * postH0 + b_ori * postH1) / 2


# List of small, medium, large and estimate of minimum value of effect size based on replication [lower CI]
d = data.frame(small=.2,
               medium=.5,
               large=.8,
               lowerCI=.001,
               observed=.49,
               upperCI=.973)


# Cost ratio of making Typ-I versus Type-II error
T1T2cratio  = 1

# Ratio of prior probabilities for H1/H0
HaHopratio = 1

optimalAlpha <- plyr::ldply(d, function(cd){
  tmp <- optab(n1=n1,
               n2=n2,
               d=abs(cd),
               T1T2cratio=T1T2cratio,
               HaHopratio=HaHopratio,
               type = type,
               tails = tails)$output

  rn <- rownames(tmp)

  tmp <- t(tmp)

  colnames(tmp) <- rn
  return(tmp)
}
)

#
# p = 0.35
# # Prior probability for H0 in replicatio study (H1 = 1-q)
# q <- .5
#
# # Posterior probability of H0 and H1 in original study
# postH0 <- ( 1 + ( (ifelse( p < (1/exp(1)), -exp(1)*p*log(p), 1) * q) / (1-q) )^-1 )^-1
# postH1 <- 1- postH0


# Sample size of replication
n1 = 3250
n2 = 3250

# List of small, medium, large and estimate of minimum value of effect size based on replication [lower CI]
d = data.frame(small=.2,
               medium=.5,
               large=.8,
               lowerCI=-.07,
               observed=-.02,
               upperCI=.03)


# Cost ratio of making Typ-I versus Type-II error
T1T2cratio  = 1

# Ratio of prior probabilities for H1/H0, based on the original posterior probability
HaHopratio = postH1/postH0

# Test and tail type
type  = "two.sample"
tails = "two.tailed"

optimalAlpha <- plyr::ldply(d, function(cd){
  tmp <- optab(n1=n1,
               n2=n2,
               d=abs(cd),
               T1T2cratio=T1T2cratio,
               HaHopratio=HaHopratio,
               type = type,
               tails = tails)$output

  rn <- rownames(tmp)

  tmp <- t(tmp)

  colnames(tmp) <- rn
  return(tmp)
  }
)


# voorbeeld -----

# Sample size of replication
n1 = 3500
n2 = 3500

# List of small, medium, large and estimate of minimum value of effect size based on replication [lower CI]
d = data.frame(small=.2,
               medium=.5,
               large=.8)

# Cost ratio of making Typ-I versus Type-II error
T1T2cratio  = 1

# Ratio of prior probabilities for H1/H0, based on the original posterior probability
HaHopratio = postH1/postH0

# Test and tail type
type  = "two.sample"
tails = "two.tailed"


optimalAlpha <- plyr::ldply(d, function(cd){
  tmp <- optab(n1=n1,
               n2=n2,
               d=abs(cd),
               T1T2cratio=T1T2cratio,
               HaHopratio=HaHopratio,
               type = type,
               tails = tails)$output

  rn <- rownames(tmp)

  tmp <- t(tmp)

  colnames(tmp) <- rn
  return(tmp)
}
)
