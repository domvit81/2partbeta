
rm(list = ls())

# --------------------------------------------
# R Code to reproduce table and figures of the manuscript entitled
# The two-part Beta regression model with mismeasured dependent variable to analyse
# quasi-formal employment in Europe
# by Arezzo MF, Guagnano G, Vitale D
# --------------------------------------------

# Likelihood function for the part one (one-side error)
llik.part1 <- function(pars.p1, data){
	beta <- pars.p1[1:(length(pars.p1)-1)]
	alpha <- pnorm(pars.p1[length(pars.p1)])
	yobs <- data$ew
	Z <- as.matrix(data.frame(data[,-(which(colnames(data)=="ew"))]))
	p2 <- (1-alpha)*pnorm(Z%*%beta)
	loglik <- (yobs==1)*log(p2) + (yobs==0)*log(1-p2)
	likf <- sum(loglik, na.rm = TRUE)
	return(likf)
	}

# Likelihood function for the part one (two-side error)
llik.part1.ext <- function(pars.p1.ext, data){
	beta <- pars.p1.ext[1:(length(pars.p1.ext)-2)]
	alpha1 <- pnorm(pars.p1.ext[length(pars.p1.ext)-1])
	alpha0 <- pnorm(pars.p1.ext[length(pars.p1.ext)])
	yobs <- data$ew
	Z <- as.matrix(data.frame(data[,-(which(colnames(data)=="ew"))]))
	p2 <- (1 - alpha1 - alpha0)*pnorm(Z%*%beta) + alpha0
	loglik <- (yobs==1)*log(p2) + (yobs==0)*log(1-p2)
	likf <- sum(loglik, na.rm = TRUE)
	return(likf)
	}


# Log-Likelihood profiling for alpha parameter
llik.part1.profile <- function(pars.p1, data){
	beta <- pars.p1[1:(length(pars.p1))]
	yobs <- data$ew
	Z <- as.matrix(data.frame(data[,-(which(colnames(data)=="ew"))]))
	p2 <- (1-pnorm(i))*pnorm(Z%*%beta)
	loglik <- (yobs==1)*log(p2) + (yobs==0)*log(1-p2)
	likf <- sum(loglik, na.rm = TRUE)
	return(likf)
	}


# Profile log-Likelihood for alpha 0 parameter
llik.part1.profile.ext <- function(pars.p1, data){
	beta <- pars.p1[1:(length(pars.p1)-1)]
	alpha0 <- pnorm(pars.p1[length(pars.p1)])
	yobs <- data$ew
	Z <- as.matrix(data.frame(data[,-(which(colnames(data)=="ew"))]))  
	p2 <- (1-pnorm(i)-alpha0)*pnorm(Z%*%beta) + alpha0 
	loglik <- (yobs==1)*log(p2) + (yobs==0)*log(1-p2)
	likf <- sum(loglik, na.rm = TRUE)
  return(likf)
}


# Profile log-Likelihood for beta parameters
llik.part1.profile.flex <- function(pars.p1, data, j){
	beta <- pars.p1[1:(length(pars.p1)-1)]
	beta[j] <- i
	alpha <- pars.p1[length(pars.p1)]
	yobs <- data$ew
	Z <- as.matrix(data.frame(data[,-(which(colnames(data)=="ew"))]))
	p2 <- (1-pnorm(alpha))*pnorm(Z%*%beta)
	loglik <- (yobs==1)*log(p2) + (yobs==0)*log(1-p2)
	likf <- sum(loglik, na.rm = TRUE)
	return(likf)
	}




#------------
# Likelihood considering UNDER-REPORTING

llik.part2.phi.cst <- function(pars.p2, data, Xmu){
	indNA <- which(!is.na(data$y))
	Wobs <- as.vector(data$y[indNA])
	X <- as.matrix(data[indNA,-1])
	k <- ncol(X)
	theta <- pars.p2[seq.int(k)]
	phi <- pars.p2[k+1]
	gamm <- pars.p2[k+2]
	mi1 <- exp(X%*%theta)/(1+exp(X%*%theta))
	mi2 <- exp(Xmu*gamm)/(1+exp(Xmu*gamm))
	mi <- mi1*mi2
	loglik <- suppressWarnings(dbeta(Wobs, shape1 = mi*phi, shape2 = (1-mi)*phi, log=TRUE)) 
	likf <- sum(loglik, na.rm = TRUE)
	return(likf)
	}

llik.part2.phi.var <- function(pars.p2, data.m, data.p, Xmu){
	indNA <- which(!is.na(data.m$y))
	Wobs <- as.vector(data.m$y[indNA])
	X <- as.matrix(data.m[indNA,-1])
	Z <- as.matrix(data.p[indNA,])
	k <- ncol(X)
	kk <- ncol(Z)
	theta <- pars.p2[seq.int(k)]
	phi <- pars.p2[(k+1):(k+kk)]
	gamm <- pars.p2[k+kk+1]
	mi1 <- exp(X%*%theta)/(1+exp(X%*%theta))
	mi2 <- exp(Xmu*gamm)/(1+exp(Xmu*gamm))
	mi <- mi1*mi2
	phi <- exp(as.vector(Z%*%phi))
	loglik <- suppressWarnings(dbeta(Wobs, shape1 = mi*phi, shape2 = (1-mi)*phi, log=TRUE))
	likf <- sum(loglik, na.rm = TRUE)
	return(likf)
	}


llik.part2.over <- function(pars.p2, data.m, data.p, Xmu){
	indNA <- which(!is.na(data.m$y))
	Wobs <- as.vector(data.m$y[indNA])
	X <- as.matrix(data.m[indNA,-1])
	Z <- as.matrix(data.p[indNA,])
	k <- ncol(X)
	kk <- ncol(Z)
	theta <- pars.p2[seq.int(k)]
	phi <- pars.p2[(k+1):(k+kk)]
	gamm <- pars.p2[k+kk+1]
	mi1 <- exp(X%*%theta)/(1+exp(X%*%theta))
	mi2 <- exp(Xmu*gamm)/(1+exp(Xmu*gamm))
	mi <- mi1/mi2
	phi <- exp(as.vector(Z%*%phi))
	loglik <- suppressWarnings(dbeta(Wobs, shape1 = mi*phi, shape2 = (1-mi)*phi, log=TRUE))
	likf <- sum(loglik, na.rm = TRUE)
	return(likf)
	}


llik.part2.profile.phi.cst <- function(pars.p2, data, Xmu){
	ind <- which(!is.na(data$y))
	Wobs <- as.vector(data$y[ind])
	X <- as.matrix(data[ind,-1])
	k <- ncol(X)
	theta <- pars.p2[seq.int(k)]
	phi <- pars.p2[k+1]
	mi1 <- exp(X%*%theta)/(1+exp(X%*%theta))
	mi2 <- exp(Xmu*i)/(1+exp(Xmu*i))
	mi <- mi1*mi2
	loglik <- suppressWarnings(dbeta(Wobs, shape1 = mi*phi, shape2 = (1-mi)*phi, log=TRUE))
	likf <- sum(loglik, na.rm = TRUE)
	return(likf)
	}


llik.part2.profile.phi.var <- function(pars.p2, data.m, data.p, Xmu){
	ind <- which(!is.na(data.m$y))
	Wobs <- as.vector(data.m$y[ind])
	X <- as.matrix(data.m[ind,-1])
	k <- ncol(X)
	Z <- as.matrix(data.p[ind,])
	kk <- ncol(Z)
	theta <- pars.p2[seq.int(k)]
	phi <- pars.p2[(k+1):(k+kk)]
	mi1 <- exp(X%*%theta)/(1+exp(X%*%theta))
	mi2 <- exp(Xmu*i)/(1+exp(Xmu*i))
	mi <- mi1*mi2
	phi1 <- exp(as.vector(Z%*%phi))
	loglik <- suppressWarnings(dbeta(Wobs, shape1 = mi*phi1, shape2 = (1-mi)*phi1, log=TRUE))
	likf <- sum(loglik, na.rm = TRUE)
	return(likf)
	}


llik.part2.profile.flex <- function(pars.p2, data.m, data.p, Xmu, j){
	ind <- which(!is.na(data.m$y))
	Wobs <- as.vector(data.m$y[ind])
	X <- as.matrix(data.m[ind,-1])
	k <- ncol(X)
	Z <- as.matrix(data.p[ind,])
	kk <- ncol(Z)
	theta <- pars.p2[seq.int(k)]
	theta[j] <- i
	gamm <- pars.p2[k+kk+1]
	phi <- pars.p2[(k+1):(k+kk)]
	mi1 <- exp(X%*%theta)/(1+exp(X%*%theta))
	mi2 <- exp(Xmu*gamm)/(1+exp(Xmu*gamm))
	mi <- mi1*mi2
	phi1 <- exp(as.vector(Z%*%phi))
	loglik <- suppressWarnings(dbeta(Wobs, shape1 = mi*phi1, shape2 = (1-mi)*phi1, log=TRUE))
	likf <- sum(loglik, na.rm = TRUE)
	return(likf)
	}

llik.part2.profile.flex.prec <- function(pars.p2, data.m, data.p, Xmu, j){
	ind <- which(!is.na(data.m$y))
	Wobs <- as.vector(data.m$y[ind])
	X <- as.matrix(data.m[ind,-1])
	k <- ncol(X)
	Z <- as.matrix(data.p[ind,])
	kk <- ncol(Z)
	theta <- pars.p2[seq.int(k)]
	gamm <- pars.p2[k+kk+1]
	phi <- pars.p2[(k+1):(k+kk)]
	phi[j] <- i
	mi1 <- exp(X%*%theta)/(1+exp(X%*%theta))
	mi2 <- exp(Xmu*gamm)/(1+exp(Xmu*gamm))
	mi <- mi1*mi2
	phi1 <- exp(as.vector(Z %*%  phi))
	loglik <- suppressWarnings(dbeta(Wobs, shape1 = mi*phi1, shape2 = (1-mi)*phi1, log=TRUE))
	likf <- sum(loglik, na.rm = TRUE)
	return(likf)
	}



#------------------------------------------------------------------------------
# DATA ANALYSIS
#------------------------------------------------------------------------------
PROFILE <- FALSE
library(ggplot2)
library(betareg)
options(scipen=999)

dati_orig <- read.csv("ebs_data.csv")

# Data set description
N <- nrow(dati_orig)
summary(dati_orig[,c("ew","y", "age", "tax_mor" )])  
table(dati_orig$east)/N
table(dati_orig$female)/N
table(dati_orig$pmi)/N
table(dati_orig$occupation)/N
table(dati_orig$form_educ)/N
table(dati_orig$det_risk)/N
table(dati_orig$sanct)/N
table(dati_orig$fin_diff)/N
table(dati_orig$resp_coop)/N

#-----------------------
#   DATA PREPARATION
#-----------------------

dati <- cbind.data.frame("uniqid" = dati_orig$uniqid, 
						 "ew" = dati_orig$ew, 
						 "y" = dati_orig$y,
                         "east" = factor(dati_orig$east), 
                         "nw" = factor(dati_orig$nw),
                         "south" = factor(dati_orig$south), 
                         "age" = dati_orig$age,
                         "female" = factor(dati_orig$female), 
                         "occupation" = factor(dati_orig$occupation),
                         "form_educ" = factor(dati_orig$form_educ), 
                         "pmi" = factor(dati_orig$pmi),
                         "det_risk" = factor(dati_orig$det_risk), 
                         "sanct" = factor(dati_orig$sanct),
                         "tax_mor" = (11 - dati_orig$tax_mor), 
                         "fin_diff" = factor(dati_orig$fin_diff),
                         "resp_coop" = 4 - dati_orig$resp_coop)
                        
Xmu <- dati$resp_coop[!is.na(dati$y)]
Xmu[which(Xmu==0)] <- 1 
table(Xmu)


#Inizialize
cov.p1 <- paste("east", "south", "age", "female", "occupation","pmi", "det_risk", "sanct", "tax_mor", "form_educ","fin_diff" , sep = "+")
cov.p2 <- paste("east", "south", "age", "female", "occupation","pmi", "det_risk", "sanct", "tax_mor", "form_educ","fin_diff", sep = "+")
cov.p2.p <- paste("east", "south",  "form_educ", "fin_diff", sep = "+")


#-----------------------
#   NAIVE MODEL
#-----------------------

naive.probit <- glm(as.formula(paste("ew ~ ", cov.p1)), data = dati, family = binomial(link="probit"))
naive.beta.const <- betareg(as.formula(paste("y ~", cov.p2)), data = dati, na.action = na.omit, link = "logit")
naive.beta <- betareg(as.formula(paste("y ~", cov.p2, "|", cov.p2.p)), data = dati, na.action = na.omit, link = "logit")

lmtest::lrtest(naive.beta.const, naive.beta)
summary(naive.beta)

#-----------------------
#   PROPOSED MODEL
#-----------------------

# Creation of dummy variables for the optim function
occ <- (model.matrix(( ~ (dati$occupation)-1)))
dr <- (model.matrix(( ~ (dati$det_risk)-1)))
san <- (model.matrix(( ~ (dati$sanct)-1)))
fd <- (model.matrix(( ~ (dati$fin_diff)-1)))
feduc <- (model.matrix(( ~ (dati$form_educ)-1)))

#----
dati.p1 <- cbind.data.frame(dati$ew, 1, as.numeric(dati$east)-1, as.numeric(dati$south)-1, dati$age, as.numeric(dati$female)-1,
                            occ[,-1], as.numeric(dati$pmi)-1, dr[,-1], san[,-1],
                            dati$tax_mor, feduc[,-1], fd[,-1])
colnames(dati.p1) <- c("ew", names(naive.probit$coefficients))

dati.p2 <- cbind.data.frame(dati$y, 1, as.numeric(dati$east)-1, as.numeric(dati$south)-1, dati$age,  as.numeric(dati$female)-1, occ[,-1], as.numeric(dati$pmi)-1, dr[,-1], san[,-1],  dati$tax_mor, feduc[,-1], fd[,-1])
colnames(dati.p2) <- c("y",names(naive.beta$coefficients$mean))

dati.p2.p <- cbind.data.frame(1, as.numeric(dati$east)-1, as.numeric(dati$south)-1, feduc[,-1], fd[,-1])
colnames(dati.p2.p) <- c(names(naive.beta$coefficients$precision))


#-------------------
# Estimation of part one
#-------------------

p1.init <- c(naive.probit$coefficients, qnorm(0.05))
est.p1 <-  optim(par = p1.init,fn=llik.part1, data=dati.p1, method="BFGS", control = list(fnscale = -1, maxit = 10000), hessian = T)
(alpha1hat <- unname(pnorm(est.p1$par[length(est.p1$par)])))
mean(pnorm(as.matrix(dati.p1[,-1])%*%est.p1$par[-length(est.p1$par)])) ## percentage of Y_true=1


#-------------------
# Estimation of part two
#--------

p2.init <- c(naive.beta$start, 3) 
est.p2 <-  optim(par=p2.init, fn=llik.part2.phi.var, data.m=dati.p2, data.p=dati.p2.p, Xmu=Xmu, method="BFGS", control = list(fnscale = -1, maxit = 5000, trace=1), hessian = T)
se.rob <-  sqrt(diag(-solve(est.p2$hess)))


if (PROFILE==TRUE){
#-------------------
# LLik profile for alpha
#-------------------
	p1.init <- c(naive.probit$coefficients)
	prof1 <- c()
	for (i in seq(-1.49,1.49,length=50)){
		prof1 <-  c(prof1, optim(par = p1.init, fn=llik.part1.profile, data=dati.p1, method="BFGS", control = list(fnscale = -1, maxit = 10000), hessian = F)$value)
	}

#-------------------
# LLik profile for gamma
#-------------------

p2.init <- c(naive.beta$start)
	prof2 <- c()
	for (i in seq(-.5,5,length=50)){
		prof2 <-  c(prof2, optim(par = p2.init, fn=llik.part2.profile.phi.var, Xmu=Xmu, data.m=dati.p2, data.p=dati.p2.p, method="BFGS", control = list(fnscale = -1, maxit = 10000, reltol=9e-14), hessian = F)$value)
	}
}


#-------------------
# Tables and Figures
#-------------------

# Part 1 - naive
out1.naive <- cbind.data.frame(c(naive.probit$coefficients, "alpha"= NA, "logLik"= logLik(naive.probit)),
									c(summary(naive.probit)$coefficients[,2], NA, NA), # St.error
									c(summary(naive.probit)$coefficients[,4], NA, NA)) # p-value
colnames(out1.naive) <- c("Estimate", "Std. Error",	"Pr(>|z|)")

# Part 1 - proposed
V1 <- (fBasics::makePositiveDefinite(-solve(est.p1$hessian)))
V1[nrow(V1),ncol(V1)] <-  V1[nrow(V1),ncol(V1)]*(dnorm(pnorm(est.p1$par[length(est.p1$par)])))^2

out1.prop <- cbind.data.frame(c(est.p1$par[-length(est.p1$par)], pnorm(est.p1$par[length(est.p1$par)]), est.p1$value),
								     c(sqrt(diag(V1)), NA), 
								     c(2*(1-pnorm(abs(est.p1$par/sqrt(diag(V1))))), NA))
colnames(out1.prop) <- c("Estimate", "Std. Error",	"Pr(>|z|)")


#Part 2 - naive
out2.naive <- cbind.data.frame(c(naive.beta$coefficients$mean, naive.beta$coefficients$precision, "gamma"= NA, naive.beta$loglik),
									c(summary(naive.beta)$coefficients$mean[,2], summary(naive.beta)$coefficients$precision[,2], "gamma"= NA, NA), # St.error
									c(summary(naive.beta)$coefficients$mean[,4], summary(naive.beta)$coefficients$precision[,4], "gamma"= NA, NA)) # p-value
colnames(out2.naive) <- c("Estimate", "Std. Error",	"Pr(>|z|)")

#Part 2 - proposed
out2.prop <- cbind.data.frame(c(est.p2$par, est.p2$value),
                              c(se.rob, NA),
                              c(2*(1-pnorm(abs(est.p2$par/se.rob))), NA))
colnames(out2.prop) <- c("Estimate", "Std. Error",	"Pr(>|z|)")

pval_gam <- pnorm((out2.prop[nrow(out2.prop)-1,1]-3)/out2.prop[nrow(out2.prop)-1,2])

out2.prop[nrow(out2.prop)-1,3] <- pval_gam
cbind(c(names(naive.beta$coeff$mean), names(naive.beta$coeff$prec), "gamma", "LLik"), out2.prop)

GAM <- out2.prop[nrow(out2.prop)-1,1] 
mean(1-exp(Xmu*GAM)/(1+exp(Xmu*GAM))) ## mean of under-reporting

UR <- 1-exp(Xmu*GAM)/(1+exp(Xmu*GAM))

p <- pnorm(as.matrix(dati.p1[,-1])%*%est.p1$par[-length(est.p1$par)])
stats::weighted.mean(p,dati_orig$wex) ## stima della proporzione di quasi formal employee di Y_true=1


## FIGURE 5
grafico <- cbind.data.frame("Resp_coop"=rev(c("Excellent", "Fair", "Bad/Average")), "Under_reporting"=rev(c(as.numeric(attributes(table(UR))$dimnames$UR))))
figUR2 <- ggplot(data=grafico, aes(x=Resp_coop, y=Under_reporting*100)) +
  		geom_bar(stat="identity",  width = 0.25) + theme_bw() + 
  		scale_x_discrete(limits=rev(c("Excellent", "Fair", "Bad/Average"))) +
  		coord_cartesian(ylim = c(0, 35)) +
  		xlab("Respondent cooperation") +
  		ylab("Under reporting (%)") +
  		theme(axis.text=element_text(size=10)) +
		theme(axis.title=element_text(size=10))

ggsave(paste0("figURbar.png"), figUR2, width=4, height=4, bg="white")


## FIGURE 6
if (PROFILE==TRUE){
png("FigIdentifiability.png", width=480*3)
par(mfrow=c(1,2), cex=2, oma=c(1,1,1,1), mar=c(3,5,1,1), cex.axis=1, cex.lab=1.25, las=1)
plot(pnorm(seq(-1.49,1.49,length=50)), prof1, type="l", lwd=5, ylab="", xlab="")
abline(h=est.p1$value - qchisq(0.95,1)/2, col=2, lty=3, lwd=3)
mtext(side=1, expression(alpha), line=2.5, cex=2.25)
mtext(side=2, "LogLik", line=4.5, las=0, cex=2)
text(pnorm(-.9), (est.p1$value- qchisq(0.95,1)/2)*.9, "95% CI", col=2)
mtext("a", side=3, line=0.2, adj=0, cex=2.5, font=2)
box(lwd=2)
plot(seq(-.5,5,length=50), prof2, ylim=c(prof2[which.max(prof2)] - (qchisq(0.95,1)), max(prof2)*1.005), type="l", ylab="", lwd=5, xlab="")
abline(h=est.p2$value - qchisq(0.95,1)/2, col=2, lty=3, lwd=3)
mtext(side=1, expression(gamma), line=2.5, cex=2.25)
mtext(side=2, "LogLik", line=3, las=0, cex=2)
#text(4, (est.p2$value - qchisq(0.95,1)/2)*1.002, "95% CI", col=2)
box(lwd=2)
mtext("b", side=3, line=0.2, adj=0, cex=2.5, font=2)
dev.off()
}

res.p1 <- round(cbind.data.frame(out1.naive," " =rep(NA, nrow(out1.naive)), out1.prop),3)
res.p2 <- round(cbind.data.frame(out2.naive," " =rep(NA, nrow(out2.naive)), out2.prop),3)

## TABLES 6 and 7 
library(xtable)
print(xtable((res.p1), type = "latex", digits=3), file = "case_study_results_p1.tex", format.args = list(big.mark = " ", decimal.mark = "."))
print(xtable((res.p2), type = "latex", digits=3), file = "case_study_results_p2.tex", format.args = list(big.mark = " ", decimal.mark = "."))




## FIGURE 7
library(car)
png("Fig7.png", width=480*4, height=480*1.25)
par(mfrow=c(1,3), cex=2, oma=c(1,1,1,1), mar=c(4,5,1,1), cex.axis=1.5, cex.lab=1.5, las=1)

# Panel a
reg1 <- as.matrix(dati.p1[,-1]) # regressori prima parte
p.true <- pnorm((reg1)%*%est.p1$par[-length(est.p1$par)])

cutoff <- seq(0.01:0.99, by=0.01)
ppc_prop <- rep(NA, length(cutoff))
ppc_naive <- rep(NA, length(cutoff))
for (i in 1:length(cutoff)) {
	ppc_naive[i] <- table(naive.probit$fitted.values[dati.p1$ew==1]>=cutoff[i])[2]/table(dati.p1$ew==1)[2]
	ppc_prop[i] <- table(p.true[dati.p1$ew==1]>=cutoff[i])[2]/table(dati.p1$ew==1)[2]
	}

plot(cutoff, ppc_naive*100, type='l', col='black', ylim = c(0,100), ylab = "", xlab="Cutoff", lwd=3) #naive
mtext(side=2, "Corrected prediction (%)", line=3.5, las=0, cex=3)
lines(cutoff, ppc_prop*100, type='l', col='black', lty=2, lwd=3) #proposed model
mtext("a", side=3, line=0.2, adj=0, cex=3, font=2)
grid(lwd=1, lty=1)
legend("topright", c("Naive", "Proposed"), lty=c(1,2), col=1, bty="n", lwd=3, cex=1.25)
box(lwd=2)

# Panel b
XX <- as.matrix(dati.p2)[!is.na(dati$y),-1]
k <- ncol(XX)
naive.hat <- as.vector(predict(naive.beta))
mi1.hat <- exp(XX%*%est.p2$par[1:k])/(1+exp(XX%*%est.p2$par[1:k])) # corrisponde alla stima dell'ytrue
mi2.hat <- exp(Xmu*est.p2$par[(length(est.p2$par))])/(1+exp(Xmu*est.p2$par[(length(est.p2$par))]))
mi.hat <- mi1.hat*mi2.hat #corrisponde alla stima dell'yobs
phi.hat <- exp(as.matrix(dati.p2.p[!is.na(dati$y),])%*%est.p2$par[ncol(dati.p2):(length(est.p2$par)-1)])
res.beta.obs <- (dati.p2$y[!is.na(dati$y)] - mi.hat)/sqrt((mi.hat*(1-mi.hat))/(1+phi.hat))
plot(res.beta.obs, pch=19, cex=.75, ylab="", xlab="i (unit)")
mtext(side=2, "Standardized residuals", line=3, las=0, cex=3)
abline(h=0, lty=2)
mtext("b", side=3, line=0.2, adj=0, cex=3, font=2)
grid(lwd=1, lty=1)
box(lwd=2)

# Panel c
qqPlot(qnorm(pbetar(na.omit(dati$y), mi.hat, phi=phi.hat)), envelope=list(level = 0.99), pch=19, cex=.5, ylab="", xlab="Normal quantiles")
mtext(side=2, "Quantile residuals", line=3, las=0, cex=3)
box(lwd=2)
mtext("c", side=3, line=0.2, adj=0, cex=3, font=2)
dev.off()

#-------------------
# LR TESTS
#-------------------

# LR test Model with precision parameter variable vs Model with precision parameter constant 
pchisq(2*(naive.beta$loglik - naive.beta.const$loglik), 6, lower.tail=FALSE)

# LR test Model with under-reporting vs Model without under-reporting (naive)
pchisq(2*(est.p2$value-naive.beta$loglik),1, lower.tail=FALSE)



## FIGURES S2, S3 and S4 of the Supplementary Material - Profile LLik for parameter identifiability

p1.init <- c(naive.probit$coefficients, qnorm(0.05))
est.p1 <-  optim(par = p1.init,fn=llik.part1, data=dati.p1, method="BFGS", control = list(fnscale = -1, maxit = 10000), hessian = T)
(alpha1hat <- unname(pnorm(est.p1$par[length(est.p1$par)])))
mean(pnorm(as.matrix(dati.p1[,-1])%*%est.p1$par[-length(est.p1$par)])) ## percentuale di Y_true=1


var.names <- c("a) Intercept", "b) East", "c) South", "d) Age", "e) Female",  "f) Occ: other white collars", "g) Occ: Manual worker",  "h) Firm size (1-99)", "i) Det risk: Fairly High", "j) Det risk: Fairly small", "k) Det risk: Very Small", "l) Sanctions: Fine", "m) Sanctions: Prison", "n) Tax morale", "o) Formal educ: 16-19 yrs", "p) Formal educ: >20 yrs",  "q) Fin diff: Sometimes", "r) Fin diff: Never")

step <- 20
png("FigIdentifiabilityBeta2.png", width=480*4, height=480*4)
par(mfrow=c(6,3), mar=c(3,5,3,5), oma=c(3,12,1,1), cex.axis=4, cex.lab=4, las=1)
for (j in 1:18){

p1.init <- c(naive.probit$coefficients, qnorm(0.05))
	prof1 <- c()
	if (j==1) for (i in seq(est.p1$par[j]-2,est.p1$par[j]+2,length=step)){
		prof1 <-  c(prof1, optim(par = p1.init, fn=llik.part1.profile.flex, data=dati.p1, j=j, method="BFGS", control = list(fnscale = -1, maxit = 10000), hessian = F)$value)
	}
	if (j==4) for (i in seq(est.p1$par[j]-.05,est.p1$par[j]+.05,length=step)){
		prof1 <-  c(prof1, optim(par = p1.init, fn=llik.part1.profile.flex, data=dati.p1, j=j, method="BFGS", control = list(fnscale = -1, maxit = 10000), hessian = F)$value)
	}
	if (j==14) for (i in seq(est.p1$par[j]-.5,est.p1$par[j]+.25,length=step)){
		prof1 <-  c(prof1, optim(par = p1.init, fn=llik.part1.profile.flex, data=dati.p1, j=j, method="BFGS", control = list(fnscale = -1, maxit = 10000), hessian = F)$value)
	}
	if (j!=1 & j!=4 & j!=14) for (i in seq(est.p1$par[j]-1,est.p1$par[j]+1,length=step)){
		prof1 <-  c(prof1, optim(par = p1.init, fn=llik.part1.profile.flex, data=dati.p1, j=j, method="BFGS", control = list(fnscale = -1, maxit = 10000), hessian = F)$value)
	}

	YUpp <- ifelse(prof1[which.max(prof1)] < est.p1$value - (qchisq(0.95,1)), prof1[which.max(prof1)]*0.8, (est.p1$value - qchisq(0.95,1))*0.99)
	YLow <- max(prof1, est.p1$value)*1.01
	if (j!=1 & j!=4 & j!=9 & j!=11 & j!=14 & j!=15) plot(seq(est.p1$par[j]-1,est.p1$par[j]+1,length=step), prof1, ylim=c(YLow, YUpp), type="l", ylab="", lwd=5, xlab="", xaxt="n")
	if (j==1) plot(seq(est.p1$par[j]-2,est.p1$par[j]+2,length=step)[-8], prof1[-8], ylim=c(YLow, YUpp), type="l", ylab="", lwd=5, xlab="", xaxt="n")
	if (j==4) plot(seq(est.p1$par[j]-.05,est.p1$par[j]+.05,length=step), prof1, ylim=c(YLow, YUpp), type="l", ylab="", lwd=5, xlab="", xaxt="n")
	if (j==9 | j==15) plot(seq(est.p1$par[j]-1,est.p1$par[j]+1,length=step)[-12], prof1[-12], ylim=c(YLow, YUpp), type="l", ylab="", lwd=5, xlab="", xaxt="n")
	if (j==11) plot(seq(est.p1$par[j]-1,est.p1$par[j]+1,length=step)[-20], prof1[-20], ylim=c(YLow, YUpp), type="l", ylab="", lwd=5, xlab="", xaxt="n")
	if (j==14) plot(seq(est.p1$par[j]-.5,est.p1$par[j]+.25,length=step), prof1, ylim=c(YLow, YUpp), type="l", ylab="", lwd=5, xlab="", xaxt="n")

	abline(h=est.p1$value - qchisq(0.95,1)/2, col=2, lty=2, lwd=3)
	mtext(side=3, var.names[j] , line=-4, cex=3.5)
	if(j==1 | j==4 | j==7 | j==10 | j==13 | j==16) mtext(side=2, "LogLik", line=11, las=0, cex=3)
	axis(1, line=1.5, col="white", cex.axis=4)
	box(lwd=3)
}
dev.off()	


step <- 30
png("FigIdentifiabilityTheta.png", width=480*4, height=480*4)
par(mfrow=c(6,3), mar=c(3,5,3,5), oma=c(3,12,1,1), cex.axis=4, cex.lab=4, las=1)
for (j in 1:18){
	p2.init <- c(naive.beta$start, 3)
	prof2 <- c()
	if (j!=1 & j!=4 & j!=14) for (i in seq(est.p2$par[j]-1, est.p2$par[j]+1, length=step)){
		prof2 <-  c(prof2, optim(par = p2.init, fn=llik.part2.profile.flex, Xmu=Xmu, j=j, data.m=dati.p2, data.p=dati.p2.p, method="BFGS", control = list(fnscale = -1, maxit = 10000, reltol=9e-14), hessian = F)$value)
	}
	if (j==1) for (i in seq(est.p2$par[j]-1.5,est.p2$par[j]+1.5,length=step)){
		prof2 <-  c(prof2, optim(par = p2.init, fn=llik.part2.profile.flex, Xmu=Xmu, j=j, data.m=dati.p2, data.p=dati.p2.p, method="BFGS", control = list(fnscale = -1, maxit = 10000, reltol=9e-14), hessian = F)$value)
	}
	if (j==4) for (i in seq(est.p2$par[j]-.05,est.p2$par[j]+.05,length=step)){
		prof2 <-  c(prof2, optim(par = p2.init, fn=llik.part2.profile.flex, Xmu=Xmu, j=j, data.m=dati.p2, data.p=dati.p2.p, method="BFGS", control = list(fnscale = -1, maxit = 10000, reltol=9e-14), hessian = F)$value)
	}
	if (j==14) for (i in seq(est.p2$par[j]-.1,est.p2$par[j]+.1,length=step)){
		prof2 <-  c(prof2, optim(par = p2.init, fn=llik.part2.profile.flex, Xmu=Xmu, j=j, data.m=dati.p2, data.p=dati.p2.p, method="BFGS", control = list(fnscale = -1, maxit = 10000, reltol=9e-14), hessian = F)$value)
	}

	YLow <- ifelse(prof2[which.max(prof2)] < est.p2$value - (qchisq(0.95,1)), prof2[which.max(prof2)]*0.8, (est.p2$value - qchisq(0.95,1))*0.99)
	YUpp <- max(prof2, est.p2$value)*1.01
	if (j!=1 & j!=4 & j!=14) plot(seq(est.p2$par[j]-1,est.p2$par[j]+1,length=step), prof2, ylim=c(YLow, YUpp), type="l", ylab="", lwd=5, xlab="", xaxt="n")
	if (j==1) plot(seq(est.p2$par[j]-1.5,est.p2$par[j]+1,length=step), prof2, ylim=c(YLow, YUpp), type="l", ylab="", lwd=5, xlab="", xaxt="n")
	if (j==4) plot(seq(est.p2$par[j]-.05,est.p2$par[j]+.05,length=step), prof2, ylim=c(YLow, YUpp), type="l", ylab="", lwd=5, xlab="", xaxt="n")
	if (j==14) plot(seq(est.p2$par[j]-.1,est.p2$par[j]+.1,length=step), prof2, ylim=c(YLow, YUpp), type="l", ylab="", lwd=5, xlab="", xaxt="n")
	abline(h=est.p2$value - qchisq(0.95,1)/2, col=2, lty=2, lwd=3)
	mtext(side=3, var.names[j] , line=-4, cex=3.5)
	if(j==1 | j==4 | j==7 | j==10 | j==13 | j==16) mtext(side=2, "LogLik", line=11, las=0, cex=3)
	axis(1, line=1, col="white", cex.axis=4)
	box(lwd=3)
}
dev.off()	


var.names <- c("a) Intercept", "b) East", "c) South", "d) Formal educ: 16-19 yrs", "e) Formal educ: >20 yrs", "f) Fin diff: Sometimes", "g) Fin diff: Never")

step <- 30
png("FigIdentifiabilityPhi.png", width=480*4, height=480*4)
par(mfrow=c(6,3), mar=c(3,5,3,5), oma=c(3,12,1,1), cex.axis=4, cex.lab=4, las=1)
for (j in 1:7){
	p2.init <- est.p2$par
	prof2 <- c()
	for (i in seq(est.p2$par[j+18]-1.25, est.p2$par[j+18]+1.25, length=step)){
		prof2 <-  c(prof2, optim(par = p2.init, fn=llik.part2.profile.flex.prec, Xmu=Xmu, j=j, data.m=dati.p2, data.p=dati.p2.p, method="BFGS", control = list(fnscale = -1, maxit = 10000, reltol=9e-14), hessian = F)$value)
	}

	YLow <- ifelse(prof2[which.max(prof2)] < est.p2$value - (qchisq(0.95,1)), prof2[which.max(prof2)]*0.8, (est.p2$value - qchisq(0.95,1))*0.99)
	YUpp <- max(prof2, est.p2$value)*1.01
	plot(seq(est.p2$par[j+18]-1,est.p2$par[j+18]+1,length=step), prof2, ylim=c(YLow, YUpp), type="l", ylab="", lwd=5, xlab="", xaxt="n")
	abline(h=est.p2$value - qchisq(0.95,1)/2, col=2, lty=2, lwd=3)
	mtext(side=3, var.names[j] , line=-4, cex=3.5)
	if(j==1 | j==4 | j==7) mtext(side=2, "LogLik", line=11, las=0, cex=3)
	axis(1, line=1, col="white", cex.axis=4)
	box(lwd=3)
}
dev.off()	
