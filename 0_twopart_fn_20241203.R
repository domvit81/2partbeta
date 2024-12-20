
# Data generating function

data.gen <- function(pars, n){
  # parameters position
  beta <- pars[1:4] 
  theta <- pars[5:8]
  phi <- pars[9:12]
  alpha <- pars[13]
  gamm <- pars[14]
  #-----------------------------------------------------------
  # Generation of part one (binary model) with misclassification 
  z1 <- rnorm(n, 0, 1)
  z2 <- rbinom(n, size=1, prob=1/3)
  z3 <- rlnorm(n, 0, 1)
  Z <- data.frame(1, z1, z2, z3)
  names(Z) <- c("cost.bin", "Z1", "Z2", "Z3")
  y <- factor(rbinom(n, 1, pnorm(beta[1] + beta[2] * z1  + beta[3]*z2 + beta[4]*z3)))
  ytrue <- y
  y <- as.factor(ytrue)
  Pi <- matrix(data = c(1, 0, alpha, 1-alpha), nrow = 2, byrow = FALSE)
  dimnames(Pi) <- list(levels(y), levels(y))
  yobs <- as.numeric(simex::misclass(data.frame(y), list(y = Pi), k = 1)[, 1])-1
  
  #---------------------------------------------------------
  # Generation of part two (beta model) with measurement error
  x1 <- z1
  x2 <- z2
  x3 <- z3
  X <- data.frame(1, x1, x2, x3)
  names(X) <- c("cost.con", "X1", "X2", "X3")
  Xmu <- 1 + rbinom(n, size=2, prob=.77) #covariata per cogliere attendibilità rispondente
  mi1 <- exp(as.matrix(X)%*%theta)/(1+exp(as.matrix(X)%*%theta))
  phi1 <- exp(as.matrix(X)%*%phi)
  p1 <- mi1*phi1
  q1 <- (1-mi1)*phi1
  SEED <-  round(sum(abs(rnorm(100)))*100000,0)
  set.seed(SEED)
  Wtrue <- rbeta(n, shape1=p1, shape2=q1)
  Wtrue[ytrue==0] <- 0
  mi2 <- exp(Xmu*gamm)/(1+exp(Xmu*gamm))
  mi <- mi1*mi2
  p <- mi*phi1
  q <- (1-mi)*phi1
  set.seed(SEED)
  Wobs <- rbeta(n, shape1=p, shape2=q)
  Wobs[yobs==0] <- NA
  db <- data.frame(yobs, Wobs, X, Z, ytrue, Wtrue, Xmu, phi1)
  return(db)
}


#*****************************************************************
# Likelihood functions

llik.part1 <- function(pars.p1, data){
	yobs <- data$yobs
	beta <- pars.p1[1:4]
	alpha <- pnorm(pars.p1[5])
	Z <- as.matrix(data[,(which(colnames(data)=="Z1")-1):which(colnames(data)=="Z3")])
	p2 <- (1-alpha)*pnorm(Z%*%beta)
	loglikpart1 <- (yobs==1)*log(p2) + (yobs==0)*log(1-p2)
	likfpart1 <- sum(loglikpart1, na.rm = TRUE)
	return(likfpart1)
	}

llik.part2 <- function(pars.p2, data){
	theta <- pars.p2[1:4]
	phi <- pars.p2[5:8]
	gamm <- pars.p2[9]
	X <- as.matrix(data[,3:6])
	Xmu <- data$Xmu
	mi1 <- exp(X[data$yobs==1,]%*%theta)/(1+exp(X[data$yobs==1,]%*%theta))
	mi2 <- exp(Xmu[data$yobs==1]*gamm)/(1+exp(Xmu[data$yobs==1]*gamm))
	phi1 <- exp(X[data$yobs==1,]%*%phi)
	mi <- mi1*mi2
	Wobs <- data$Wobs
	loglikpart2 <- dbeta(Wobs[data$yobs==1], shape1 = mi*phi1, shape2 = (1-mi)*phi1, log=TRUE)
	likfpart2 <- sum(loglikpart2, na.rm = TRUE)
	return(likfpart2)
	}

