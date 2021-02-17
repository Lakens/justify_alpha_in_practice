library(pwr)
library(ggplot2)

#functions to justify alpha levels. 
#first run all the code. In the beginning are functions incoroporated in later functions

#call idealalpha(n, effsize, prior, test = "t-test", type = "two.sample", alternative = "two.sided", k = 4)
#to calculate the alpha level which gives optimal learning
#the output will be alpha, power and info, which is the expected learning

#call alphaConfirmation(n, effsize, prior, postdis, test = "t-test", type = "two.sample", alternative = "two.sided", k = 4)
#to calculate the alpha level needed to achieve a given posterior after a significant result
#the output will be alpha, power, info, and confirm, which is the probability after a positive result

#alphaDisconfirmation(n, effsize, prior, postdis, test = "t-test", type = "two.sample", alternative = "two.sided", k = 4)
#to calculate the alpha level needed to achieve a given posterior after a non-significant result
#the output will be alpha power and disconfirm, which is the probability after a negative result

#probpower(effsize, prior, postdis, postconf, test = "t-test", type = "two.sample", alternative = "two.sided", k = 4)
#conduct a probability based power analysis. Here posterior after significance and posterior after non-signficance are fixed
#the function will then return the sample size and the alpha level needed
#the output will be alpha, power, disconfirm, confirm, info and n, the sample size per group



#The next functions are defined to be part of the four main functions

#calculate the probability that the hypothesis is true, after confirming
postconfirm <- function(power,alpha,prior)
                          {
        (power*prior)/(power*prior+alpha*(1-prior))
}

#calculate the probability that the hypothesis is true, after a nonsignificant finding
postdisconfirm <- function(power, alpha, prior)
{
  (1-power)*prior/((1-alpha)*(1-prior)+(1-power)*prior)
}

#define power function for t.test, anova and test for correlations based on pwr package
get_power   <- function(alpha, effsize, n, test = "t-test", type = "two.sample", alternative = "two.sided", k = 4){ #get_power for different alpha levels, two tail independent samples t-test
if(test == "t-test") {
  power <- pwr.t.test(n, effsize, alpha, type = type, alternative = alternative)
}
if(test == "anova") {
  power <- pwr.anova.test(k, n, effsize, alpha)
}
if(test == "correlation") {
  power <- pwr.r.test(n, effsize, alpha, alternative = alternative)
}
  return(power$power)
}

#calculate the expected amount of correct change in believe
information <- function(confirm, disconfirm, power, alpha, prior) #calcualte the amount of information per test as %correct change in believe about the hypothesis
            {
      rightupdate <- (prior*power*(confirm - prior) + (1-prior)*(1-alpha)*(prior -disconfirm)) 
      wrongupdate <- ((1-prior)*alpha*(confirm-prior) + (1-power)*prior*(prior - disconfirm)) #you could also update incorrectly for false postives and false negatives
      information <- rightupdate - wrongupdate #the expected correct change in belief in the hypothesis
}



#four main functions start here

#idealalpha: Calculate the alpha level that optimizes the amount of learning

idealalpha <- function(n, effsize, prior, test = "t-test", type = "two.sample", alternative = "two.sided", k = 4) #return the best alpha and power as well as the correct change in believe and plot change in believe as a function of alpha
{
  alpha   <- seq(0,1,.00001) #supply a vector of alphas
  power  <- get_power(alpha, effsize, n, test = test, type = type, alternative = alternative, k = k) # calculate power
  confirm <- postconfirm(power, alpha, prior) #posterior after confirming
  disconf <- postdisconfirm(power, alpha, prior) #posterior after disconfirming
  info <- information(confirm, disconf, power, alpha, prior) #expected correct change in believe
  d <- as.data.frame(cbind(alpha,power,info)) 
  maxinf <- max(na.omit(info))
  dmax <- as.data.frame(subset(d, info == maxinf)) #return valus with highes correct change in believe
  d05 <- as.data.frame(subset(d, alpha == 0.05)) #compare to 5% alpha
  plot(alpha,info, ylab = "Expected Correct Change in Belief", xlab = "Significance Level", type = "l", main = "Optimal Alpha", xlim = c(0,1), axes = F, ylim =c(0,round(maxinf + 0.005, digits = 2)))
  abline(v = 0.05, col = "red")
  abline(v = dmax$alpha, col = "blue")
  axis(side=2, at=c(seq(0,(maxinf+0.01),0.01)))
  axis(side=1, at=c(seq(0,1,0.05)))
  mtext(side =3, paste("Alpha = ", round(dmax$alpha, digits = 2), " Power = ", round(dmax$power, digits = 2), "Learning =", round(dmax$info*100, digits = 2), "%"))
  dmax <- cbind(dmax, d05$info)
  colnames(dmax) <- c("alpha", "power", "Learning", "Learning05")
  return(dmax)
}

#alphaConfirmation: Indicate your posterior probability after finding a positive result. 
#The function will return the alpha level needed to achieve this posterior probability

#the user specifies the probability of the hypothesis after confirming the prediction and the alpha level is adjusted accordingly
alphaConfirmation <- function(n, effsize, prior, posttrue, test = "t-test", type = "two.sample", alternative = "two.sided", k = 4)
{
  if(prior > posttrue) #stop if posterior smaller than prior
  {
    stop("Prior must be smaller than posterior")
  }
  alpha   <- seq(0,1,.00001) #define vector of alphas
  #calculate powers 
  power  <- get_power(alpha, effsize, n, test = test, type = type, alternative = alternative, k = k)
  #calculate posteriors after confirming or disconfirming and expected learning
  confirm <- postconfirm(power, alpha, prior)
  disconf <- postdisconfirm(power, alpha, prior)
  info <- information(confirm, disconf, power, alpha, prior)
  df <- as.data.frame(cbind(alpha,power,info, confirm))
  #calculate the smallest confirming posterior that is above the user specified posterior
  #the alpha level corresponding to this posterior is what we are looking for.
  maxinf <- max(na.omit(info))
  postt <- as.data.frame(subset(df, confirm > posttrue))
  postt <- min(postt$confirm)
  d <- subset(df, confirm == postt)
  #plot
  plot(alpha,info, ylab = "Expected Change in Belief", xlab = "Signifcance Level", type = "l", main = paste("Alpha for",posttrue, "Probability after Significance") , xlim = c(0,1), axes = F, ylim =c(0,round(maxinf + 0.005, digits = 2)))
  abline(v = 0.05, col = "red")
  abline(v = d$alpha, col = "blue")
  axis(side=2, at=c(seq(0,(maxinf+0.01),0.01)))
  axis(side=1, at=c(seq(0,1,0.05)))
  mtext(side =3, paste("Alpha = ", round(d$alpha, digits = 2), " Power = ", round(d$power, digits = 2), "Learning =", round(d$info*100, digits = 2), "%"))
  return(d)
}

#AlphaDisconfirmation: Indicate the posterior probability you want after a negative result
#The function will return the alpha level needed to achieve this posterior probability after testing

#now do the same and solve for postfalse
alphaDisconfirmation <- function(n, effsize, prior, postdis, test = "t-test", type = "two.sample", alternative = "two.sided", k = 4)
{
  if(postdis > prior)
  {
    stop("Prior must be higher than posterior")
  }
  alpha   <- seq(0,1,.00001)
  power  <- get_power(alpha, effsize, n, test = test, type = type, alternative = alternative, k = k) #power.prop.test(n, p1 = 0.8, p2 = 1, sig.level = alpha)$power
  confirm <- postconfirm(power, alpha, prior)
  disconf <- postdisconfirm(power, alpha, prior)
  info <- information(confirm, disconf, power, alpha, prior)
  df <- as.data.frame(cbind(alpha,power,info, disconf))
  #f <- max(na.omit(info))
  maxinf <- max(na.omit(info))
  postd <- as.data.frame(subset(df, disconf > postdis))
  postd <- min(postd$disconf)
  d <- subset(df, disconf == postd)
  #plot
  plot(alpha,info, ylab = "Expected Change in Belief", xlab = "Signifcance Level", type = "l", main = paste("Alpha for",round(postd, digits = 2), "Probability after Non-Significance") , xlim = c(0,1), axes = F, ylim =c(0,round(maxinf + 0.005, digits = 2)))
  abline(v = 0.05, col = "red")
  abline(v = d$alpha, col = "blue")
  axis(side=2, at=c(seq(0,(maxinf+0.01),0.01)))
  axis(side=1, at=c(seq(0,1,0.05)))
  mtext(side =3, paste("Alpha = ", round(d$alpha, digits = 2), " Power = ", round(d$power, digits = 2), "Learning =", round(d$info*100, digits = 2), "%"))
  return(d)
}

#Probability based power analysis: Set both the posterior probability after a positive result and the posterior probability after a negative result
#The function will calculate the required sample size and alplha level to achieve this. 

probpower <- function(effsize, prior, postdis, postconf, test = "t-test", type = "two.sample", alternative = "two.sided", k = 4)
{

  if(!postdis < prior | !prior < postconf)
  {
    stop("Posterior after disconfirming needs to be smaller than prior, posterior after confirming larger")
  }
  alpha <- seq(.0001,1,.001) 
  #for 10-20 000 particpants calculte the alpha needed to achieve the posterior after a negative result
  #stop the loop when there is an option that at the same time has the needed posterior probability after a positive result
  for (i in 10:20000) { 
  power <- get_power(alpha, effsize, i, test, type, alternative, k) 
  #calculate power for all alphas for this sample size
  disconf <- postdisconfirm(power, alpha, prior)
  confirm <- postconfirm(power, alpha, prior)
  info <- information(confirm, disconf, power, alpha, prior)
  maxinf <- max(na.omit(info))
  df <- as.data.frame(cbind(alpha,power, disconf, confirm, info))
  postd <- as.data.frame(subset(df, disconf > postdis)) 
  postd <- min(postd$disconf)
  d <- subset(df, disconf == postd)
  d <- cbind(d, i)
  postc <- d$confirm
  
  if(postd > postdis & postc > postconf)
  {
    break
  }
  }
  colnames(d) <- (c("alpha", "power", "disconf", "confirm", "info", "n"))
  plot(df$alpha,df$info, ylab = "Expected Correct Change in Belief", xlab = "Signifcance Level", type = "l", main = paste("Probability bas Power Analysis - PosteriorConfirmation:",postconf, "  Posterior Disconfirmation:", postdis) , xlim = c(0,1), axes = F, ylim =c(0,round(maxinf + 0.005, digits = 2)))
  abline(v = 0.05, col = "red")
  abline(v = d$alpha, col = "blue")
  axis(side=2, at=c(seq(0,(maxinf+0.01),0.01)))
  axis(side=1, at=c(seq(0,1,0.05)))
  mtext(side =3, paste("Sample Size", i, " Alpha = ", round(d$alpha, digits = 2), " Power = ", round(d$power, digits = 2), "Learning =", round(d$info*100, digits = 2), "%"))
  return(d)
}


