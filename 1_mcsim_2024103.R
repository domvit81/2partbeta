
rm(list = ls())

require(betareg)
require(simex)
require(robustbase)
require(parallel)
require(doParallel)
require(foreach)
require(lqmm) ## for positiveDefinite
options(scipen=999)

source('0_twopart_fn_20241203.R')

dir.create("scen")


nrep <- 200
n <- 10000
alpha <- c(0, 0.05, 0.4, 0.8)
beta0 <- -1.4 #corresponding to 10% ytrue and 2% of yobs when alpha = 0.8 
phi <-   c(3.5,.3,.5,-.3)
gamma <- c(.35, .7, 1.25, 3)       

beta1 <- 1
beta2 <- .5
beta3 <- -.5
theta0 <- -.7
theta1 <- .6
theta2 <- .8
theta3 <- .4

schema <- expand.grid("beta0"= beta0, "alpha"= alpha, "gamma"= gamma, "n"= n)

#-----------
# Simulation run

for(j in 1:nrow(schema)){
	
	true.par.p1 <- c("beta0"= schema$beta0[j],"beta1"= beta1, "beta2"= beta2, "beta3"= beta3, "alpha" = schema$alpha[j])
  	true.par.p2 <- c("theta0"=theta0, "theta1"= theta1, "theta2"=theta2, "theta3"=theta3, "phi0"= phi[1], "phi1"= phi[2], "phi2"= phi[3], "phi3"= phi[4], "gamma"=schema$gamma[j])

	cl <- makeCluster(detectCores()-1) # set the number of cores supported by CPU
	registerDoParallel(cl)
	invisible(clusterEvalQ(cl = cl, source("0_twopart_fn_20241203.R")))		

	result <- foreach(i = 1:nrep, .combine="rbind", .packages = c("betareg")) %dopar% {

  		dati <- data.gen(pars=c(schema$beta0[j], beta1, beta2, beta3, theta0, theta1, theta2, theta3, phi, schema$alpha[j], schema$gamma[j]), n=schema$n[j])
  		dati$Wobs[dati$Wobs>0.99] <- 0.99
  		dati$Wobs[dati$Wobs<0.01] <- 0.01
  
	 	naive.probit <- glm(yobs ~ Z1 + Z2 + Z3, data = dati, family = binomial(link="probit"))
		naive.beta <- tryCatch(betareg::betareg(Wobs ~ X1+X2+X3 | X1+X2+X3, data=dati, na.action = na.omit, link = "logit"), error=function(cond){return(NA)})
  
		p1.init <- c(naive.probit$coefficients, qnorm(0.05))
		#p1.init <- c(0,0,0,0, qnorm(0.0001))
		est.p1 <- tryCatch(optim(par = p1.init, fn=llik.part1, data=dati, 
							method="BFGS", control = list(fnscale = -1, maxit = 10000, reltol=9e-14), 
                				hessian = T) , error=function(cond){return(list("par"=rep(NA, length(p1.init)), "value"=NA))})
  		alphahat <- tryCatch(pnorm(est.p1$par[length(est.p1$par)]),  error=function(cond){return(NA)})
  		se.p1 <- c(tryCatch(sqrt(diag((-solve(est.p1$hessian)))), error=function(cond){return(rep(NA,length(p1.init)))})[1:(length(true.par.p1)-1)],
             tryCatch(sqrt(diag((-solve(est.p1$hessian)))), error=function(cond){return(rep(NA,length(p1.init)))})[length(true.par.p1)]*dnorm(pnorm(est.p1$par[length(true.par.p1)])))


		#p2.init <- c(naive.beta$start, 1)
		p2.init <- c(rep(0, 4),1,0,0,0, 3)
		est.p2 <- tryCatch(optim(par=p2.init, fn=llik.part2, data=dati, method="BFGS", control = list(fnscale = -1, maxit = 10000, reltol=9e-14), hessian = T),
       							error=function(cond){return(list("par"=rep(NA,9), "value"=NA))})
  		se.p2 <- tryCatch(sqrt(diag(-solve(est.p2$hessian))), error=function(cond){return(rep(NA, 9))})
  
  		c(true.par.p1, 
			summary(naive.probit)$coeff[,1], 
			summary(naive.probit)$coeff[,2], 
    			c(est.p1$par[-length(est.p1$par)], alphahat), 
 			se.p1,
			true.par.p2, 
			c(summary(naive.beta)$coefficients$mean[,1], summary(naive.beta)$coefficients$precision[,1]), 
			c(summary(naive.beta)$coefficients$mean[,2], summary(naive.beta)$coefficients$precision[,2]), 
			est.p2$par,
			se.p2)
		}

	colnames(result) <- c(
		names(true.par.p1), 
		paste0(names(true.par.p1)[-5], "_est_naive"), 
		paste0(names(true.par.p1)[-5], "_se_naive"), 
		paste0(names(true.par.p1), "_est_prop"), 
		paste0(names(true.par.p1), "_se_prop"),
		names(true.par.p2), 
		paste0(names(true.par.p2)[-9], "_est_naive"),
		paste0(names(true.par.p2)[-9], "_se_naive"),
		paste0(names(true.par.p2), "_est_prop"),
		paste0(names(true.par.p2), "_se_prop"))

	write.table(result, paste0("./Scen/Res_Scen_",j,".csv"), quote=FALSE, row.names=FALSE, sep=",")
	
	stopCluster(cl)

}