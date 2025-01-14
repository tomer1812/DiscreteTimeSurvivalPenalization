
rm(list = ls())


library(tidyverse)
library(survival)
library(FSA)
library(foreach)
library(doParallel)
library(dplyr)


outputfile <- "bias_results_comparison_different_alpha_large_d.csv"


# Number of observations to sample n:
n <- 5000


# Number of Discrete-times
d <- 50
times <- c(1:d)
maxt <- length(times)

minimal_occurrences_required <- 3
max_exact <- 7000

if (n > max_exact) {
  ours_method = "efron"
} else {
  ours_method = "exact"
}
print(paste('Using', ours_method, "in our method"))

# Definition of the true logit link function model parameters (see [1])
beta11 <- -0.5*c(log(0.8),log(3),log(3),log(2.5),log(2))
beta12 <- -0.5*c(log(1),log(3),log(4),log(3),log(2))

alpha1 <- -4 + 0.07*times
alpha2 <- -5.3 + 0.07*times
true.coef <- c(alpha1,beta11,alpha2,beta12)


# Number of features p (the length of the Z vector) 
p <- length(beta11) 


# Number of repetitions
rep <- 500
dim.time <- length(times)



sampling.discrete.competing.Cox <- function(i) {
  # Covariates sampling
  z <- runif(5,0,1) 
  
  # Beta*Z calculation
  bz1 <- as.vector(beta11%*%z)
  bz2 <- as.vector(beta12%*%z)  
  
  # Hazard function J=1
  exp1t <- exp(alpha1+bz1)/(1+exp(alpha1+bz1)) 
  
  # Hazard function J=2
  exp2t <- exp(alpha2+bz2)/(1+exp(alpha2+bz2)) 
  
  # Overall survival function
  st <- c(1,cumprod(1-exp1t-exp2t)[1:(dim.time-1)]) 
  
  # Probability of event occurrence function J=1
  prob1t <- exp1t*st 
  
  # Probability of event occurrence function J=2
  prob2t <- exp2t*st 
  
  # Marginal probability J=1
  p1 <- sum(prob1t) 
  
  # Marginal probability J=2
  p2 <- sum(prob2t)
  
  # Conditional probability for event occurrence at t given the event type J=1
  prob1t.j1 <- prob1t/p1
  
  # Conditional probability for event occurrence at t given the event type J=2
  prob2t.j2 <- prob2t/p2 
  
  # Probability for no observed event by Tmax
  p0 <- 1-p1-p2
  
  # Event type sampling
  j <- sample(c(0,1,2),1,prob=c(p0,p1,p2)) 

  # Event time sampling
  if (j==0) {
    tt <- max(times)+1
    datai <- c(tt,0,z)
  } else 
  {
    if (j==1) { 
      tt <- sample(times,1,prob = prob1t.j1)
      datai <- c(tt,j,z)
    } else {
      tt <- sample(times,1,prob = prob2t.j2)
      datai <- c(tt,j,z)
    }
  }
  
  # Right censoring time sampling - random censoring
  cen.prob <- c(rep(0.002,max(times)),1-sum(rep(0.002,max(times))))
  cen <- sample(1:(max(times)+1),1,prob=cen.prob) 

  # Update event time and type in case of censored observation
  if (cen < tt)   {
    datai <- c(cen,0,z)
  }
  return(datai)
}



data.expansion <- function(n, dat) {
  Expdata <- NULL
  for (i in 1:n) {
    for (j in 1:dat$time[i]) {
      if (dat$time[i]==j & dat$status[i]==1) {
        delta1 = 1; delta2 = 0;
      } else {
        if (dat$time[i]==j & dat$status[i]==2) {
          delta1 = 0; delta2 =1;
        } else {
          delta1 = delta2 = 0;
        }
      }
      Expdata <- rbind(Expdata,c(dat$ID[i],j,delta1,delta2,1-delta1-delta2,dat$cov1[i],dat$cov2[i],dat$cov3[i],dat$cov4[i],dat$cov5[i]))
    }
  }
  Expdata <- as.data.frame(Expdata)
  names(Expdata) <- c("ID","time","delta1","delta2","delta3","cov1","cov2","cov3","cov4","cov5")
  
  Expdata <- Expdata[order(Expdata$ID, Expdata$time), ]
  
  return(Expdata)
}


data.expansion_parallel <- function(n, dat) {
  # Register parallel backend
  num_cores <- 4  # Use one less core than available
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Perform parallel processing
  Expdata <- foreach(i = 1:nrow(dat), .combine = rbind, .packages = "base") %dopar% {
    Expdata_i <- do.call(rbind, lapply(1:dat$time[i], function(j) {
      if (dat$time[i] == j & dat$status[i] == 1) {
        delta1 <- 1; delta2 <- 0
      } else if (dat$time[i] == j & dat$status[i] == 2) {
        delta1 <- 0; delta2 <- 1
      } else {
        delta1 <- delta2 <- 0
      }
      c(dat$ID[i], j, delta1, delta2, 1 - delta1 - delta2,
        dat$cov1[i], dat$cov2[i], dat$cov3[i], dat$cov4[i], dat$cov5[i])
    }))
    Expdata_i
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  # Convert the result to a data frame and add column names
  Expdata <- as.data.frame(Expdata)
  names(Expdata) <- c("ID", "time", "delta1", "delta2", "delta3", "cov1", "cov2", "cov3", "cov4", "cov5")
  
  Expdata <- Expdata[order(Expdata$ID, Expdata$time), ]
  
  return(Expdata)
}




obj_fun <- function(a0,j,t,ndat,prob1,prob2) {
  if (j==1) {
    temp <- filter(ndat,time>=t)
    probt <- exp(a0+temp$betaZ1)/(1+exp(a0+temp$betaZ1))
    dif <- (mean(probt)-prob1[t])^2     
  } else {
    temp <- filter(ndat,time>=t)
    probt <- exp(a0+temp$betaZ2)/(1+exp(a0+temp$betaZ2))
    dif <- (mean(probt)-prob2[t])^2  
  }
  return(dif)
}



cov1.hat.simple <- matrix(0,nrow = rep,ncol = maxt+p)
cov2.hat.simple <- matrix(0,nrow = rep,ncol = maxt+p)

cov1.hat.cox.breslow <- matrix(0,nrow = rep,ncol = p)
cov1.hat.cox.breslow.bhzd <- matrix(0,nrow = rep,ncol = maxt)
cov2.hat.cox.breslow <- matrix(0,nrow = rep,ncol = p)
cov2.hat.cox.breslow.bhzd <- matrix(0,nrow = rep,ncol = maxt)

cov1.hat.cox.efron <- matrix(0,nrow = rep,ncol = p)
cov1.hat.cox.efron.bhzd <- matrix(0,nrow = rep,ncol = maxt)
cov2.hat.cox.efron <- matrix(0,nrow = rep,ncol = p)
cov2.hat.cox.efron.bhzd <- matrix(0,nrow = rep,ncol = maxt)

cov1.hat.cox.exact <- matrix(0,nrow = rep,ncol = p)
cov1.hat.cox.exact.bhzd <- matrix(0,nrow = rep,ncol = maxt)
cov2.hat.cox.exact <- matrix(0,nrow = rep,ncol = p)
cov2.hat.cox.exact.bhzd <- matrix(0,nrow = rep,ncol = maxt)

sam <- 1

# Repetitions loop
for (iii in 1:(rep+100)) {
  # Data sampling
  start_time <- Sys.time()
  print(sam)
  minimal_occurrences = 0
  while (minimal_occurrences < minimal_occurrences_required) {
    dat <- lapply(1:n, sampling.discrete.competing.Cox)
    dat <- as.data.frame(t(matrix(unlist(dat),2+p,n)))
    names(dat) <- c("time","statusJ","cov1","cov2","cov3","cov4","cov5")

    dat$ID <- 1:n
    ndat <- as_tibble(dat)
  
    table(dat$time,dat$statusJ)
    
    statusJ_range <- c(1, 2)

    # Filter the dataframe and count occurrences for each combination
    result <- dat %>%
      filter(statusJ %in% statusJ_range, time %in% times) %>%
      group_by(statusJ, time) %>%
      summarise(count = n(), .groups = "drop") %>%
      complete(statusJ = statusJ_range, time = times, fill = list(count = 0))
    
    # Find the minimal value of the count
    minimal_occurrences <- min(result$count)
    print(paste("Mid minimal occurrences: ", minimal_occurrences))

  }
  print(paste('Minimal occurrences in rep', sam, 'is', minimal_occurrences))
  
  if (ours_method == 'exact') {
    coxsurv1 <- Surv(time=dat$time, event=ifelse(dat$statusJ == 1, 1, 0))
    cox1 <- coxph(coxsurv1 ~ cov1 + cov2 + cov3 + cov4 + cov5, data=dat, ties='exact')
    base_hazard_1 <- basehaz(cox1, centered=FALSE)
    
    
    coxsurv2 <- Surv(time=dat$time, event=ifelse(dat$statusJ == 2, 1, 0))
    cox2 <- coxph(coxsurv2 ~ cov1 + cov2 + cov3 + cov4 + cov5, data=dat, ties='exact')
    base_hazard_2 <- basehaz(cox2, centered=FALSE)
    
    if (any(is.na(c(cox1$coefficients, cox2$coefficients)))) {
      print(paste("NA in exact solution of iteration", sam))
      next
    }
  
    cov1.hat.cox.exact[sam,] <- cox1$coefficients
    cov2.hat.cox.exact[sam,] <- cox2$coefficients
    cov1.hat.cox.exact.bhzd[sam,] <- qlogis(diff(c(0, base_hazard_1$hazard[1:d]))) 
    cov2.hat.cox.exact.bhzd[sam,] <- qlogis(diff(c(0, base_hazard_2$hazard[1:d]))) 
  }
  
  coxsurv1 <- Surv(time=dat$time, event=ifelse(dat$statusJ == 1, 1, 0))
  cox1 <- coxph(coxsurv1 ~ cov1 + cov2 + cov3 + cov4 + cov5, data=dat, ties='breslow')
  base_hazard_1 <- basehaz(cox1, centered=FALSE)
  
    
  coxsurv2 <- Surv(time=dat$time, event=ifelse(dat$statusJ == 2, 1, 0))
  cox2 <- coxph(coxsurv2 ~ cov1 + cov2 + cov3 + cov4 + cov5, data=dat, ties='breslow')
  base_hazard_2 <- basehaz(cox2, centered=FALSE)
  
  cov1.hat.cox.breslow[sam,] <- cox1$coefficients
  cov2.hat.cox.breslow[sam,] <- cox2$coefficients
  cov1.hat.cox.breslow.bhzd[sam,] <- qlogis(diff(c(0, base_hazard_1$hazard[1:d]))) 
  cov2.hat.cox.breslow.bhzd[sam,] <- qlogis(diff(c(0, base_hazard_2$hazard[1:d]))) 
  
  
  coxsurv1 <- Surv(time=dat$time, event=ifelse(dat$statusJ == 1, 1, 0))
  cox1 <- coxph(coxsurv1 ~ cov1 + cov2 + cov3 + cov4 + cov5, data=dat, ties='efron')
  base_hazard_1 <- basehaz(cox1, centered=FALSE)
  
  
  coxsurv2 <- Surv(time=dat$time, event=ifelse(dat$statusJ == 2, 1, 0))
  cox2 <- coxph(coxsurv2 ~ cov1 + cov2 + cov3 + cov4 + cov5, data=dat, ties='efron')
  base_hazard_2 <- basehaz(cox2, centered=FALSE)
  
  cov1.hat.cox.efron[sam,] <- cox1$coefficients
  cov2.hat.cox.efron[sam,] <- cox2$coefficients
  cov1.hat.cox.efron.bhzd[sam,] <- qlogis(diff(c(0, base_hazard_1$hazard[1:d]))) 
  cov2.hat.cox.efron.bhzd[sam,] <- qlogis(diff(c(0, base_hazard_2$hazard[1:d]))) 
  
  
  # Data Expansion
  
  print('Expansion')

  
  Expdata <- data.expansion_parallel(n, dat)
  Expdata <- as.data.frame(Expdata)
  

  # Estimation of the event specific beta_j coefficients using the approach of Meir and Gorfine (2023) [1]
  # Note that method = "approximate" refers to Breslow approximation for the likelihood function
  # Other options such as "efron" and "exact" are available
  
  print('Ours')
  fit1.short <- clogit(delta1 ~ cov1+cov2+cov3+cov4+cov5  + strata(time), data = Expdata, method = ours_method)
  fit2.short <- clogit(delta2 ~ cov1+cov2+cov3+cov4+cov5  + strata(time), data = Expdata, method = ours_method) 

  ndat <- ndat %>%
    mutate(betaZ1 =  cbind(cov1,cov2,cov3,cov4,cov5)%*%fit1.short$coefficients,
           betaZ2 =  cbind(cov1,cov2,cov3,cov4,cov5)%*%fit2.short$coefficients) %>%
    arrange(time)
  prob1 <- vector(); prob2 <- vector()
  for (t in 1:maxt) {
    n.t <- dim(filter(ndat,time>=t))[1]
    n.t.1 <- dim(filter(ndat,time==t,statusJ==1))[1]
    n.t.2 <- dim(filter(ndat,time==t,statusJ==2))[1]
    prob1[t] <- n.t.1/n.t
    prob2[t] <- n.t.2/n.t
  }

  # Estimation of the event specific alpha_jt coefficients using the approach of Meir and Gorfine (2023) [1]
  alpha.mat <- matrix(0,maxt,2)
  for (j in 1:2) {
    for (t in 1:maxt) {
      res <- nlminb(-3,obj_fun,j=j,t=t,ndat=ndat,prob1=prob1,prob2=prob2)
      alpha.mat[t,j] <- res$par
    }
  }

  cov1.hat.simple[sam,] <- c(alpha.mat[,1],fit1.short$coefficients)
  cov2.hat.simple[sam,] <- c(alpha.mat[,2],fit2.short$coefficients)

  end_time <- Sys.time()
  iteration_time <- difftime(end_time, start_time, units = "secs")[[1]]
  print(paste("Time taken for iteration", sam, ":", iteration_time))
  sam <- sam + 1
  
  if (sam > rep) {
    break
  }
}


est.coef.cox.breslow <- apply( cbind(cov1.hat.cox.breslow.bhzd, cov1.hat.cox.breslow, cov2.hat.cox.breslow.bhzd, cov2.hat.cox.breslow)  ,2,mean)
est.coef.cox.breslow.se <- apply( cbind(cov1.hat.cox.breslow.bhzd, cov1.hat.cox.breslow, cov2.hat.cox.breslow.bhzd, cov2.hat.cox.breslow)  ,2,se)
diff.coef.cox.breslow <- apply( t( t(cbind(cov1.hat.cox.breslow.bhzd, cov1.hat.cox.breslow, cov2.hat.cox.breslow.bhzd, cov2.hat.cox.breslow)) - true.coef)  ,2,mean)
diff.coef.cox.breslow.se <- apply( t( t(cbind(cov1.hat.cox.breslow.bhzd, cov1.hat.cox.breslow, cov2.hat.cox.breslow.bhzd, cov2.hat.cox.breslow)) - true.coef)  ,2,se)


est.coef.cox.efron <- apply( cbind(cov1.hat.cox.efron.bhzd, cov1.hat.cox.efron, cov2.hat.cox.efron.bhzd, cov2.hat.cox.efron)  ,2,mean)
est.coef.cox.efron.se <- apply( cbind(cov1.hat.cox.efron.bhzd, cov1.hat.cox.efron, cov2.hat.cox.efron.bhzd, cov2.hat.cox.efron)  ,2,se)
diff.coef.cox.efron <- apply( t(t(cbind(cov1.hat.cox.efron.bhzd, cov1.hat.cox.efron, cov2.hat.cox.efron.bhzd, cov2.hat.cox.efron)) - true.coef)  ,2,mean)
diff.coef.cox.efron.se <- apply( t(t(cbind(cov1.hat.cox.efron.bhzd, cov1.hat.cox.efron, cov2.hat.cox.efron.bhzd, cov2.hat.cox.efron)) - true.coef)  ,2,se)


est.coef.cox.exact <- apply( cbind(cov1.hat.cox.exact.bhzd, cov1.hat.cox.exact, cov2.hat.cox.exact.bhzd, cov2.hat.cox.exact)  ,2,mean)
est.coef.cox.exact.se <- apply( cbind(cov1.hat.cox.exact.bhzd, cov1.hat.cox.exact, cov2.hat.cox.exact.bhzd, cov2.hat.cox.exact)  ,2,se)
diff.coef.cox.exact <- apply( t(t(cbind(cov1.hat.cox.exact.bhzd, cov1.hat.cox.exact, cov2.hat.cox.exact.bhzd, cov2.hat.cox.exact))  - true.coef),2,mean)
diff.coef.cox.exact.se <- apply( t(t(cbind(cov1.hat.cox.exact.bhzd, cov1.hat.cox.exact, cov2.hat.cox.exact.bhzd, cov2.hat.cox.exact)) - true.coef)  ,2,se)


est.coef.discrete <- apply( cbind(cov1.hat.simple, cov2.hat.simple)  ,2,mean)
est.coef.discrete.se <- apply( cbind(cov1.hat.simple, cov2.hat.simple)  ,2,se)
diff.coef.discrete <- apply( t(t(cbind(cov1.hat.simple, cov2.hat.simple)) - true.coef)  ,2,mean)
diff.coef.discrete.se <- apply( t(t(cbind(cov1.hat.simple, cov2.hat.simple)) - true.coef)  ,2,se)


results.dat <- as.data.frame( cbind(true.coef, 
                                    est.coef.cox.breslow, 
                                    est.coef.cox.breslow.se, 
                                    diff.coef.cox.breslow, 
                                    diff.coef.cox.breslow.se, 
                                    est.coef.cox.efron, 
                                    est.coef.cox.efron.se,
                                    diff.coef.cox.efron, 
                                    diff.coef.cox.efron.se,
                                    est.coef.cox.exact,
                                    est.coef.cox.exact.se,
                                    diff.coef.cox.exact,
                                    diff.coef.cox.exact.se,                                    
                                    est.coef.discrete,
                                    est.coef.discrete.se,
                                    diff.coef.discrete,
                                    diff.coef.discrete.se
                                    ) )

write.csv(results.dat, file = outputfile, row.names = FALSE)

warnings()

