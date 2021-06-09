setwd("D:/TCD_MScDataScience/StatisticalModelling/LabFiles")
df1 <- read.csv("simpsons.csv")
#change the Season variable to be treated as categorical(factor)
df1$Season <- factor(df1$Season)
dim(df1)
nlevels(df1$Season)

tapply(df1$Rating, df1$Season, mean)
tapply(df1$Rating, df1$Season, median)
tapply(df1$Rating, df1$Season, sd)

#plotting to visualize the data closely
library(ggplot2)
ggplot(df1) + geom_boxplot(aes(x = reorder(Season, Rating, median), Rating, fill = reorder(Season, Rating, median)), show.legend=FALSE)

ggplot(df1, aes(x = reorder(Season, Season, length))) + stat_count()

## stat_bin() using bins = 30. Pick better value with binwidth
ggplot(df1, aes(Rating)) + stat_bin()

ggplot(data.frame(size = tapply(df1$Rating, df1$Season, length), 
                  mean_rating = tapply(df1$Rating, df1$Season, mean)), 
       aes(size, mean_rating)) + geom_point()

##gibbs sampling
#??, the overall mean across Seasons;
#??b, the precision (inverse variability) between Seasons;
#??w, the precision (inverse variability) within Seasons;
#??m, the mean Rating awarded to Season m.

compare_m_gibbs <- function(y, ind, maxiter = 5000)
{
  
  ### weakly informative priors
  a0 <- 1/2 ; b0 <- 50 ## tau_w hyperparameters
  eta0 <-1/2 ; t0 <- 50 ## tau_b hyperparameters
  mu0<-5 ; gamma0 <- 1/4
  ###
  
  ### starting values
  m <- nlevels(ind)
  ybar <- theta <- tapply(y, ind, mean)
  tau_w <- mean(1 / tapply(y, ind, var)) ##within group precision
  mu <- mean(theta)
  tau_b <-var(theta) ##between group precision
  n_m <- tapply(y, ind, length)
  an <- a0 + sum(n_m)/2
  ###
  
  ### setup MCMC
  theta_mat <- matrix(0, nrow=maxiter, ncol=m)
  mat_store <- matrix(0, nrow=maxiter, ncol=3)
  ###
  
  ### MCMC algorithm
  for(s in 1:maxiter) 
  {
    
    # sample new values of the thetas
    for(j in 1:m)
    {
      taun <- n_m[j] * tau_w + tau_b
      thetan <- (ybar[j] * n_m[j] * tau_w + mu * tau_b) / taun
      theta[j]<-rnorm(1, thetan, 1/sqrt(taun))
    }
    
    #sample new value of tau_w
    ss <- 0
    for(j in 1:m){
      ss <- ss + sum((y[ind == j] - theta[j])^2)
    }
    bn <- b0 + ss/2
    tau_w <- rgamma(1, an, bn)
    
    #sample a new value of mu
    gammam <- m * tau_b + gamma0
    mum <- (mean(theta) * m * tau_b + mu0 * gamma0) / gammam
    mu <- rnorm(1, mum, 1/ sqrt(gammam)) 
    
    # sample a new value of tau_b
    etam <- eta0 + m/2
    tm <- t0 + sum((theta - mu)^2) / 2
    tau_b <- rgamma(1, etam, tm)
    
    #store results
    theta_mat[s,] <- theta
    mat_store[s, ] <- c(mu, tau_w, tau_b)
  }
  colnames(mat_store) <- c("mu", "tau_w", "tau_b")
  return(list(params = mat_store, theta = theta_mat))
}

fit1 <- compare_m_gibbs(df1$Rating, df1$Season)
#fit$theta will have 32 sets of parameters stored (because there are 32 different seasons in the data), so we will focus first on fit$params, which contains only the main parameters ??,??b and ??w.
apply(fit1$params, 2, mean)
apply(fit1$params, 2, sd)

mean(1/sqrt(fit1$params[, 3]))
sd(1/sqrt(fit1$params[, 3]))


## reformat samples for ggplot
theta_df <- data.frame(samples = as.numeric(fit1$theta), 
                       Season = rep(1:ncol(fit1$theta), each = nrow(fit1$theta))) 

theta_med <- apply(theta_df, 2, mean) ## get basic posterior summary
sort(theta_med, decreasing = TRUE) ## which seasons did best and worst?

ggplot(theta_df) + geom_boxplot(aes(x = reorder(Season, samples, median), samples, 
                                    fill = reorder(Season, samples, median)), show.legend=FALSE)

#We can also plot the theta_hat against the relative sample size of each Season
theta_hat <- apply(fit1$theta, 2, mean)
ggplot(data.frame(size = tapply(df1$Rating, df1$Season, length), theta_hat = theta_hat), aes(size, theta_hat)) + geom_point()
