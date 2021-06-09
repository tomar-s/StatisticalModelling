setwd("D:/TCD_MScDataScience/StatisticalModelling/LabFiles")
df1 <- read.csv("simpsons_2seasons.csv")
df1$Season <- factor(df1$Season)
dim(df1)


library(ggplot2)
ggplot(df1) + geom_boxplot(aes(Season, Rating, fill = Season)) + geom_jitter(aes(Season, Rating, shape = df1$Season))
## Warning: Use of `df1$Season` is discouraged. Use `Season` instead.

tapply(df1$Rating, df1$Season, mean)
tapply(df1$Rating, df1$Season, median)
tapply(df1$Rating, df1$Season, sd)

#perform a t-test to compare the group means
t.test(Rating ~ Season, data=df1, var.equal = TRUE)

##Gibbs Sampling

compare_2_gibbs <- function(y, ind, mu0 = 5, tau0 = 1/4, del0 = 0, gamma0 = 1/4, a0 = 1, b0 = 50, maxiter = 5000)
{
  y1 <- y[ind == 2]
  y2 <- y[ind == 6]
  
  n1 <- length(y1) 
  n2 <- length(y2)
  
  ##### starting values
  mu <- (mean(y1) + mean(y2)) / 2
  del <- (mean(y1) - mean(y2)) / 2
  
  mat_store <- matrix(0, nrow = maxiter, ncol = 3)
  #####
  
  ##### Gibbs sampler
  an <- a0 + (n1 + n2)/2
  
  for(s in 1 : maxiter) 
  {
    
    ##update tau
    bn <- b0 + 0.5 * (sum((y1 - mu - del) ^ 2) + sum((y2 - mu + del) ^ 2))
    tau <- rgamma(1, an, bn)
    ##
    
    ##update mu
    taun <-  tau0 + tau * (n1 + n2)
    mun <- (tau0 * mu0 + tau * (sum(y1 - del) + sum(y2 + del))) / taun
    mu <- rnorm(1, mun, sqrt(1/taun))
    ##
    
    ##update del
    gamman <-  gamma0 + tau*(n1 + n2)
    deln <- ( del0 * gamma0 + tau * (sum(y1 - mu) - sum(y2 - mu))) / gamman
    del<-rnorm(1, deln, sqrt(1/gamman))
    ##
    
    ## store parameter values
    mat_store[s, ] <- c(mu, del, tau)
  }
  colnames(mat_store) <- c("mu", "del", "tau")
  return(mat_store)
}


library(MCMCpack)
fit <- compare_2_gibbs(df1$Rating, as.factor(df1$Season))

plot(as.mcmc(fit))

raftery.diag(as.mcmc(fit))
apply(fit, 2, mean)
 
