# Initial code for project 1, doing a survival oriented analysis.
# model checking: http://webpages.math.luc.edu/~ebalderama/myfiles/modelchecking101_pres.pdf
# open bugs doc: https://www.mrc-bsu.cam.ac.uk/wp-content/uploads/2021/06/OpenBUGS_Manual.pdf
# bugs issues: https://www.mrc-bsu.cam.ac.uk/wp-content/uploads/bugsbook_chapter9.pdf
# code for reparameterization of hierarchical effects: https://zenodo.org/record/3743847#.YbTydb3MKM8
library("R2OpenBUGS")
library(dplyr)
library(tidyr)
library(coda)
library(ggplot2)

setwd("C:/Users/USER/Documents/ku_leuven/courses/bayesian_II/project1")
#setwd("D:/asus_documents/ku_leuven/courses/bayesian_II/project1/")
df <- read.csv("Grubs_Nematodes.csv", header=T)
head(df)
dim(df)
df %>% select(everything()) %>% summarise_all(list(~sum(is.na(.))))
# Take midpoint between lower and upper limit, accounting for censoring
df %>%
  mutate(no_cen=ifelse(is.na(UPPERLIM),12,UPPERLIM))%>%
  mutate(mid=(no_cen+LOWERLIM)/2) -> df
head(df)
    
# Split data into interval-censored and right-censored sets
# Adjust groups for parameterization
df$GROUP <- df$GROUP - 1
df_right <- df[is.na(df$UPPERLIM), ]
df_int <- df[!is.na(df$UPPERLIM),]
n_right <- nrow(df_right)
n_int <- nrow(df_int)
# Multilevel survival model, page 396, Bayesian Biostatistics
# https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-9-9
N <- nrow(df)
J <- length(unique(df$UREPID))
set.seed(13)
# First a Weibull model with random effects but no censoring
inits <- function(){
    list(betas = rnorm(3, 0, 3), plate_unident = rnorm(J, 0, 3),
         rho=runif(1, 0, 10))
}
#inits1 <-list(betas=c(0, 0, 0), rho = 10, plate=rep(1,J))
#inits2 <-list(betas=c(1, 1, 1), rho = .1, plate=rep(0, J))
#inits3 <-list(betas=c(-1, -1, -1), rho = 1, plate=rep(-1, J))
#inits <- list(inits1, inits2, inits3)
parameters <- c("betas", "beta_star", "plate", "rho", "lambda", "ppo", "icpo")
model_file="weibull_model_no_censor_random.txt"
# Center main effect
data <- list(N=N, J=J, y=df$mid, group=df$GROUP, grubsize=df$GRUBSIZE-mean(df$GRUBSIZE), urepid=df$UREPID)
print(Sys.time())
wb_no_censor <- bugs(data, inits=inits, model.file=model_file, 
                 parameters=parameters, n.chains=3,
                    n.burnin=1000, n.iter=10000, n.thin=10, save.history=T, codaPkg=T, bugs.seed=13)

coda_wb_no_censor <- read.bugs(wb_no_censor)
#    Plot random 
matrix_coda <- as.matrix(coda_wb_no_censor)
colnames(matrix_coda)
hist(matrix_coda[,286:305])
# god I hate this package.
# Ignore betas[1] 
traceplot(coda_wb_no_censor[,c(1, 3:4)])
traceplot(coda_wb_no_censor[,286:305])
# Rho
traceplot(coda_wb_no_censor[,446])


# Check convergence
params2check <- c(1, 3:4, 286:305)
geweke.diag(coda_wb_no_censor[[1]][, params2check])
geweke.diag(coda_wb_no_censor[[2]][, params2check])
geweke.diag(coda_wb_no_censor[[3]][, params2check])
# Check with Gelman
# Looks good
gelman.diag(coda_wb_no_censor[, params2check])
gelman.plot(coda_wb_no_censor[, params2check])

crosscorr.plot(coda_wb_no_censor[, 286:305])
crosscorr.plot(coda_wb_no_censor[, c(1, 3:4)])


# Next a lognormal model no with censoring but random effects.
parameters <- c("beta_star", "betas", "plate", "ppo", "icpo", "test", "ks", "ks.rep")
model_file <- "lognormal_no_censor_random.txt"
inits <- function(){
    list(betas = rnorm(3, 0, 3), plate_unident = rnorm(J, 0, 3),
         tau=runif(1, 0, 10))
}
data <- list(N=N, J=J, logy=log(df$mid), group=df$GROUP, grubsize=df$GRUBSIZE-mean(df$GRUBSIZE), 
             urepid=df$UREPID)
print(Sys.time())
ln_no_censor <- bugs(data, inits=inits, model.file=model_file, parameters=parameters, n.chains=3,
                 n.burnin=1000, n.iter=10000, save.history=T, codaPkg=T)
coda_ln_no_censor <- read.bugs(ln_no_censor)
matrix_coda_ln_no_censor <- as.matrix(coda_ln_no_censor)
colnames(matrix_coda_ln_no_censor)
hist(matrix_coda_ln_no_censor[,148:167])
# Random eff
traceplot(coda_ln_no_censor[,148:167])
# main eff
traceplot(coda_ln_no_censor[,c(1, 3:4)])
# Formal checks
params2check <- c(1, 3:4, 148:167)
geweke.diag(coda_ln_no_censor[[1]][, params2check])
geweke.diag(coda_ln_no_censor[[2]][, params2check])
geweke.diag(coda_ln_no_censor[[3]][, params2check])
gelman.diag(coda_ln_no_censor[, params2check])
# low corr
crosscorr.plot(coda_ln_no_censor[, 148:167])
# Gelman
gelman.plot(coda_ln_no_censor[, params2check])


# Then weibull model with censoring and random effects
# Right and interval censoring are accounted for in different loops
parameters <- c("beta_star", "betas", "plate", "rho", "ppo_r", "ppo_int", "icpo_r", "icpo_int")
model_file <- "weibull_model_censor_random.txt"
data <- list(n_right=n_right, n_int=n_int, J=J, y_r=df_right$mid, y_int=df_int$mid, 
             right_cens=df_right$LOWERLIM, int_lower_cens=df_int$LOWERLIM, int_upper_cens=df_int$UPPERLIM,
             group_r=df_right$GROUP, grubsize_r=df_right$GRUBSIZE-mean(df_right$GRUBSIZE), urepid_r=df_right$UREPID,
             group_int=df_int$GROUP, grubsize_int=df_int$GRUBSIZE-mean(df_int$GRUBSIZE), urepid_int=df_int$UREPID)
# Using same inits as above
print(Sys.time())
wb_censor <- bugs(data, inits=inits, model.file=model_file, parameters=parameters, n.chains=3,
                 n.burnin=1000, n.iter=10000, n.thin=10, save.history=T, codaPkg=T, bugs.seed=13)

coda_wb_censor <- read.bugs(wb_censor)
matrix_coda_wb <- as.matrix(coda_wb_censor)
colnames(matrix_coda_wb)
hist(matrix_coda_wb[146:165])

# Random eff
traceplot(coda_wb_censor[, 146:165])
# main eff
traceplot(coda_wb_censor[,c(1, 3:4)])
# Rho
traceplot(coda_wb_censor[,306])
# Convergence checks
params2check <- c(1, 3:4, 146:165, 306)
geweke.diag(coda_wb_censor[, params2check])
gelman.diag(coda_wb_censor[, params2check])
gelman.plot(coda_wb_censor[, params2check])
# low correlation.
crosscorr.plot(coda_wb_censor[, 146:165])
# Everything has converged properly

# And finally a lognormal model with censoring and random effects.
parameters <- c("beta_star", "betas", "plate", "ppo_r", "ppo_int", "icpo_r", "icpo_int", "test", "ks", "ks.rep")
model_file <- "lognormal_censor_random.txt"
#inits1 <- list(betas = c(0, 0, 0), tau_r = 1/0.5, tau_int=1/0.5, plate=rep(1,J))
#inits2 <- list(betas = c(1, 1, 1), tau_r = 1/5, tau_int=1/5, plate=rep(0, J))
#inits3 <- list(betas = c(-1, -1, -1), tau_r = 1/0.05, tau_int=1/0.05, plate=rep(-1, J))
#inits <- list(inits1, inits2, inits3)
inits <- function(){
    list(betas = rnorm(3, 0, 3), plate_unident = rnorm(J, 0, 3),
         tau_r=runif(1, 0, 10), tau_int=runif(1, 0, 10))
}
data <- list(n_right=n_right, n_int=n_int, J=J, logy_r=log(df_right$mid), logy_int=log(df_int$mid), 
             right_cens=df_right$LOWERLIM, int_lower_cens=df_int$LOWERLIM, int_upper_cens=df_int$UPPERLIM,
             group_r=df_right$GROUP, grubsize_r=df_right$GRUBSIZE-mean(df_right$GRUBSIZE), urepid_r=df_right$UREPID,
             group_int=df_int$GROUP, grubsize_int=df_int$GRUBSIZE-mean(df_int$GRUBSIZE), urepid_int=df_int$UREPID)
print(Sys.time())
ln_censor <- bugs(data, inits=inits, model.file=model_file, parameters=parameters, n.chains=3,
                 n.burnin=1000, n.iter=10000, save.history=T, codaPkg=T)
coda_ln_censor <- read.bugs(ln_censor)
matrix_coda_ln <- as.matrix(coda_ln_censor)
colnames(matrix_coda_ln)
hist(matrix_coda_ln[,148:167])
# Random eff
traceplot(coda_ln_censor[,148:167])
# main eff
traceplot(coda_ln_censor[,c(1, 3:4)])
# Formal checks
params2check <- c(1, 3:4, 148:167)
geweke.diag(coda_ln_censor[[1]][, params2check])
geweke.diag(coda_ln_censor[[2]][, params2check])
geweke.diag(coda_ln_censor[[3]][, params2check])
gelman.diag(coda_ln_censor[, params2check])
gelman.plot(coda_ln_censor[, params2check])
# low corr
crosscorr.plot(coda_ln_censor[, 148:167])


# KS test for normality
ks <- matrix_coda_ln[,"ks"]
ks.rep <-  matrix_coda_ln[,"ks.rep"]
# P value from monte carlo test
P_ks_test <- mean(matrix_coda_ln[,"test[1]"])
# We do not reject normality
# However, test has low power because just based off of one y_i
# When there are more times when ks is smaller than ks.rep, 
# p val will be larger. When ks is larger than ks.rep most of the 
# time, we get a smaller p value.
P_ks_test

minks <- min(ks,ks.rep)
maxks <- max(ks,ks.rep)

par(pty="s")
plot(ks, ks.rep,xlim=c(minks,maxks),ylim=c(minks,maxks),
     xlab="Observed D_ks", ylab="Replicated D_ks",
     cex.lab=1.5,cex.main=1.8,cex.axis=1.3,col="dark blue",
     main="Kolmogorov-Smirnov test")
abline(0,1,lwd=2,col="red")


# PSBF
# From andy's code
icpo_ind <- 6:145
df_icpo_ln <- dplyr::bind_rows(lapply(coda_ln_censor[, icpo_ind], as.data.frame))
df_icpo_wb <- dplyr::bind_rows(lapply(coda_wb_censor[, icpo_ind], as.data.frame))
CPO_weib <- 1/colMeans(df_icpo_wb)
CPO_ln <- 1/colMeans(df_icpo_ln)
logPSBF <- sum(log(CPO_ln)) - sum(log(CPO_weib))
# Evidence for Lognormal model when considering random effects and censoring
logPSBF

