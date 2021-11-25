tempmod<-"model {
#data# temps
#data# N
#data# S
for (s in 1:S){
  means[s] ~ dunif(20, 30)
  vars[s] ~ dnorm(10, 0.5)
  for (i in 1:N){
    temps[i,s] ~ dnorm(means[s], vars[s])
  }
}
#monitor# means
#monitor# vars
#inits# means
#inits# vars
}"
library(runjags)
library(rjags)
temps <- cbind(rnorm(1000, 25, 10), rnorm(1000, 25, 10))
N <- nrow(temps)
S <- ncol(temps)
means <-list(chain1=c(25, 25), chain2=c(25, 25))
vars <- list(chain1=c(10, 10), chain2= c(10, 10))  	# a named list or a function
results <- run.jags(tempmod, burnin=5000, sample=10000, n.chains=2)
results  			# Prints summary, HPD, convergence / MC statistics
plot(results, plot.type=c('t','d'))
mcmcobj <- as.mcmc.list(results)


library(runjags)
mod<-"model {
for (i in 1:N){
#data# Data
#data# N
Data[i] ~ dnorm(mean, precision)
}
precision~ dgamma(0.001,0.001)
#data# Precision
mean~dnorm(0, Precision)
#monitor# precision
#monitor# mean
#inits# mean
#inits# precision
}"
library(runjags)
Data <- rnorm(100, mean=12, sd=2)
N <- length(Data)
Precision <- 1/(1000^2)
mean <- list(chain1=0, chain2=1)
precision <- list(chain1=0.1, chain2=0.001)
initial <- list(chain1=list(0,0.1), chain2=list(1, 0.001))
results <- run.jags(mod, burnin=5000, sample=10000, n.chains=2)
results  			# Prints summary, HPD, convergence / MC statistics
plot(results, type=c('t','d'))
mcmcobj <- as.mcmc.list(results)
