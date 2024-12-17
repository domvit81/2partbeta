
library(ggplot2)
library(ggpubr)
library(ggbreak) 
library(xtable)
library(SimDesign)

options(scipen=999)


my_theme <-	theme(legend.text = element_text(size = 18), 
                  legend.title = element_text(size = 18), 
                  legend.key.size = unit(1.5, 'cm'))

path_in <- paste0(getwd(), "/scen")

n <- 10000
ALPHA <- c(0, 0.05, 0.4,  0.8)
BETA <- c(-1.4, 1, .5, -.5)
THETA <- c(-.7, .6, .8, .4)
PHI <-    c(3.5,.3,.5,-.3)
GAMMA <- c(.35, .7, 1.25, 3)
schema <- expand.grid("beta0"= BETA[1], "alpha"= ALPHA, "gamma"= GAMMA ,"n"= n)
nsim <- 200

lf <- list.files(path_in)

########################################################################################
##
## RELATIVE BIAS MODEL PARAMETERS
##
########################################################################################

plot_beta <- list() 
plot_theta <- list() 
plot_phi <- list() 

for (i in 1:length(lf)){
	set <- read.csv(paste0(path_in, "/Res_Scen_",i,".csv"))
	NOMICOL <- colnames(set)
	
	rb_beta <- data.frame(	
"RelBias"=c(as.vector(unlist(sweep(sweep(set[,c(6:9)], MARGIN=2, FUN="-", STATS=BETA), MARGIN=2, FUN="/", STATS=BETA)*100)),
as.vector(unlist(sweep(sweep(set[,c(14:17)], MARGIN=2, FUN="-", STATS=BETA), MARGIN=2, FUN="/", STATS=BETA)*100))),
"Coef" = rep(c(rep("B0", nsim), rep("B1", nsim), rep("B2", nsim), rep("B3", nsim)),2),
"Model"=c(rep("Naive", nsim*4), rep("Proposed", nsim*4)))

	rb_theta <- data.frame(	
"RelBias"=c(as.vector(unlist(sweep(sweep(set[,c(33:36)], MARGIN=2, FUN="-", STATS=THETA), MARGIN=2, FUN="/", STATS=THETA)*100)),
as.vector(unlist(sweep(sweep(set[,c(49:52)], MARGIN=2, FUN="-", STATS=THETA), MARGIN=2, FUN="/", STATS=THETA)*100))),
"Coef" = rep(c(rep("T0", nsim), rep("T1", nsim), rep("T2", nsim), rep("T3", nsim)),2),
"Model"=c(rep("Naive", nsim*4), rep("Proposed", nsim*4)))


	rb_phi <- data.frame(	
"RelBias"=c(as.vector(unlist(sweep(sweep(set[,c(37:40)], MARGIN=2, FUN="-", STATS=PHI), MARGIN=2, FUN="/", STATS=PHI)*100)),
as.vector(unlist(sweep(sweep(set[,c(53:56)], MARGIN=2, FUN="-", STATS=PHI), MARGIN=2, FUN="/", STATS=PHI)*100))),
"Coef" = rep(c(rep("P0", nsim), rep("P1", nsim), rep("P2", nsim), rep("P3", nsim)),2),
"Model"=c(rep("Naive", nsim*4), rep("Proposed", nsim*4)))




label <- c("a", "b", "c", "d")
plot_beta[[i]] <- ggplot(rb_beta, aes(x = Coef, y = RelBias)) +
		geom_boxplot(aes(fill = Model), position = position_dodge(0.8), outlier.size = 0.25, show.legend = TRUE, lwd=.5, outlier.shape = 1) +
		scale_fill_manual(values = c("#999999", "white")) + 
		xlab(" ") +
		ylab("Relative bias (%)") +
		ggtitle(label[i]) +
		scale_x_discrete(breaks=c("B0","B1","B2", "B3"), labels=c(expression(beta[0]), expression(beta[1]), expression(beta[2]), expression(beta[3]))) +
		theme(axis.text=element_text(size=24)) + 
		theme(axis.title=element_text(size=24)) +
		theme(plot.title = element_text(size = 30, face = "bold")) +
		geom_hline(yintercept = 0, colour="dodgerblue", lty=2) +
		my_theme


level_orderX <- c("T0","T1","T2","T3")
label <- c("d", "h", "l", "p", "c", "g", "k", "o", "b", "f", "j", "n", "a", "e", "i", "m")
plot_theta[[i]] <- ggplot(rb_theta, aes(x = Coef, y = RelBias)) +
		geom_boxplot(aes(fill = Model), position = position_dodge(0.8), outlier.size = 0.25, show.legend = TRUE, lwd=.5, outlier.shape = 1) +
		scale_fill_manual(values = c("#999999", "white")) + 
		xlab(" ") +
		ylab("Relative bias (%)") +
		ggtitle(label[i]) +
		scale_x_discrete(breaks=c("T0","T1","T2","T3"), labels=c(expression(theta[0]), expression(theta[1]), expression(theta[2]), expression(theta[3])))+
		theme(axis.text=element_text(size=24)) +
		theme(axis.title=element_text(size=24)) +
		theme(plot.title = element_text(size = 30, face = "bold")) +
        geom_hline(yintercept = 0, colour="dodgerblue", lty=2) +
		my_theme


level_orderX <- c("P0","P1","P2","P3")
label <- c("d", "h", "l", "p", "c", "g", "k", "o", "b", "f", "j", "n", "a", "e", "i", "m")
plot_phi[[i]] <- ggplot(rb_phi, aes(x = factor(Coef, level=level_orderX), y = RelBias)) + 
		geom_boxplot(aes(fill = Model), position = position_dodge(0.8), outlier.size = 0.25, show.legend = TRUE, lwd=.5, outlier.shape = 1) +
		scale_fill_manual(values = c("#999999", "white")) + 
		xlab(" ") +
		ylab("Relative bias (%)") +
		ggtitle(label[i]) +
		scale_x_discrete(breaks=c("P0","P1","P2","P3"), labels=c(expression(phi[0]), expression(phi[1]), expression(phi[2]), expression(phi[3])))+
		theme(axis.text=element_text(size=24)) +
		theme(axis.title=element_text(size=24)) +
		theme(plot.title = element_text(size = 30, face = "bold")) +
        geom_hline(yintercept = 0, colour="dodgerblue", lty=2) +
		my_theme


}


tl <- 80
fig_rb_beta <- ggarrange(
plot_beta[[1]] + coord_cartesian(ylim = c(-tl, tl)),
plot_beta[[2]] + coord_cartesian(ylim = c(-tl, tl)),
plot_beta[[3]] + coord_cartesian(ylim = c(-tl, tl)),
plot_beta[[4]] + coord_cartesian(ylim = c(-tl, tl)),
ncol = 4, nrow = 1, common.legend = TRUE, legend="bottom")
ggsave(paste0("fig_rb_beta.png"), fig_rb_beta, width=18, height=5.5, bg="white")


fig_rb_theta <- ggarrange(
plot_theta[[13]] + coord_cartesian(ylim = c(-tl, tl)),
plot_theta[[9]] + coord_cartesian(ylim = c(-tl, tl)),
plot_theta[[5]] + coord_cartesian(ylim = c(-tl, tl)),
plot_theta[[1]] + coord_cartesian(ylim = c(-tl, tl)),
plot_theta[[14]] + coord_cartesian(ylim = c(-tl, tl)),
plot_theta[[10]] + coord_cartesian(ylim = c(-tl, tl)),
plot_theta[[6]] + coord_cartesian(ylim = c(-tl, tl)),
plot_theta[[2]] + coord_cartesian(ylim = c(-tl, tl)),
plot_theta[[15]] + coord_cartesian(ylim = c(-tl, tl)),
plot_theta[[11]] + coord_cartesian(ylim = c(-tl, tl)),
plot_theta[[7]] + coord_cartesian(ylim = c(-tl, tl)),
plot_theta[[3]] + coord_cartesian(ylim = c(-tl, tl)),
plot_theta[[16]] + coord_cartesian(ylim = c(-tl, tl)),
plot_theta[[12]] + coord_cartesian(ylim = c(-tl, tl)),
plot_theta[[8]] + coord_cartesian(ylim = c(-tl, tl)),
plot_theta[[4]] + coord_cartesian(ylim = c(-tl, tl)),
ncol = 4, nrow = 4, common.legend = TRUE, legend="bottom")
ggsave(paste0("fig_rb_theta.png"), fig_rb_theta, width=18, height=18, bg="white")



fig_rb_phi <- ggarrange(
plot_phi[[13]] + coord_cartesian(ylim = c(-tl, tl)),
plot_phi[[9]] + coord_cartesian(ylim = c(-tl, tl)),
plot_phi[[5]] + coord_cartesian(ylim = c(-tl, tl)),
plot_phi[[1]] + coord_cartesian(ylim = c(-tl, tl)),
plot_phi[[14]] + coord_cartesian(ylim = c(-tl, tl)),
plot_phi[[10]] + coord_cartesian(ylim = c(-tl, tl)),
plot_phi[[6]] + coord_cartesian(ylim = c(-tl, tl)),
plot_phi[[2]] + coord_cartesian(ylim = c(-tl, tl)),
plot_phi[[15]] + coord_cartesian(ylim = c(-tl, tl)),
plot_phi[[11]] + coord_cartesian(ylim = c(-tl, tl)),
plot_phi[[7]] + coord_cartesian(ylim = c(-tl, tl)),
plot_phi[[3]] + coord_cartesian(ylim = c(-tl, tl)),
plot_phi[[16]] + coord_cartesian(ylim = c(-tl, tl)),
plot_phi[[12]] + coord_cartesian(ylim = c(-tl, tl)),
plot_phi[[8]] + coord_cartesian(ylim = c(-tl, tl)),
plot_phi[[4]] + coord_cartesian(ylim = c(-tl, tl)),
ncol = 4, nrow = 4, common.legend = TRUE, legend="bottom")
ggsave(paste0("fig_rb_phi.png"), fig_rb_phi, width=18, height=18, bg="white")



########################################################################################
##
## RELATIVE BIAS ALPHA - GAMMA
##
########################################################################################

plota <- list() 
plotg <- list() 

for (i in 1:length(lf)){
	set <- read.csv(paste0(path_in, "/Res_Scen_",i,".csv"))
	NOMICOL <- colnames(set)
	
	rba <- data.frame(RelBias=unlist(as.vector(set[,"alpha_est_prop"])), Coef=rep("alpha",nsim), Scenario=paste0("S", i)) 	  
	rbg <- data.frame(RelBias=unlist(as.vector(set[,"gamma_est_prop"])), Coef=rep("gamma",nsim), Scenario=paste0("S", i))



label <- c("a", "b", "c", "d")
plota[[i]] <- ggplot(rba, aes(x = Coef, y = RelBias)) +
		geom_boxplot(aes(fill = Scenario), outliers=FALSE, position = position_dodge(0.8), outlier.size = 0.25, show.legend = TRUE, lwd=.5, outlier.shape = 1) +
		geom_hline(yintercept=schema[i,2], col="dodgerblue", lty=2, lwd=1)+
		scale_fill_manual(values = c("white")) + 
		xlab(" ") +
		ylab(expression(alpha~~estimate)) +
		ggtitle(label[i]) +
		scale_x_discrete(breaks="alpha", labels=" ") +
		theme(axis.text=element_text(size=24)) + 
		theme(axis.title=element_text(size=24)) +
		theme(plot.title = element_text(size = 30, face = "bold")) +
		my_theme


label <- c("d", "h", "l", "p", "c", "g", "k", "o", "b", "f", "j", "n", "a", "e", "i", "m")
plotg[[i]] <- ggplot(rbg, aes(x = Coef, y = RelBias)) + 
		geom_boxplot(aes(fill = Scenario), outliers=FALSE, position = position_dodge(0.8), outlier.size = 0.25, show.legend = TRUE, lwd=.5, outlier.shape = 1) +
		geom_hline(yintercept=schema[i,3], col="dodgerblue", lty=2, lwd=1)+
		scale_fill_manual(values = c("white")) + 
		xlab(" ") +
		ylab(expression(gamma~~estimate)) +
		ggtitle(label[i]) +
		scale_x_discrete(breaks="gamma", labels=" ") +
		theme(axis.text=element_text(size=24)) +
		theme(axis.title=element_text(size=24)) +
		theme(plot.title = element_text(size = 30, face = "bold")) +
		my_theme
}

fig_alpha <- ggarrange(
plota[[1]], 
plota[[2]],
plota[[3]],
plota[[4]],
ncol = 4, nrow = 1, common.legend = TRUE, legend="none")
ggsave(paste0("fig_alpha.png"), fig_alpha, width=18, height=4, bg="white")


fig_gamma <- ggarrange(
plotg[[13]],
plotg[[9]],
plotg[[5]],
plotg[[1]],
plotg[[14]],
plotg[[10]],
plotg[[6]],
plotg[[2]],
plotg[[15]],
plotg[[11]],
plotg[[7]],
plotg[[3]],
plotg[[16]],
plotg[[12]],
plotg[[8]],
plotg[[4]],
ncol = 4, nrow = 4, common.legend = TRUE, legend="none")
ggsave(paste0("fig_gamma.png"), fig_gamma, width=18, height=16, bg="white")




########################################################################################
##
## COVERAGE
##
########################################################################################


cvg_beta <- matrix(NA, nrow=16, ncol=8)
cvg_theta <- matrix(NA, nrow=16, ncol=8)
cvg_phi <- matrix(NA, nrow=16, ncol=8)

ord_scen <- c(13,9,5,1,14,10,6,2,15,11,7,3,16,12,8,4)

for (i in 1:length(lf)){
	
	set <- read.csv(paste0(path_in, "/Res_Scen_",ord_scen[i],".csv"))
	NOMICOL <- colnames(set)
	
	cvg_beta[i,] <- c(
	sum(set[,"beta0"] >= (set[,"beta0_est_naive"] - 1.96*set[,"beta0_se_naive"]) & set[,"beta0"] <= (set[,"beta0_est_naive"] + 1.96*set[,"beta0_se_naive"]), na.rm=TRUE)/nsim*100,
	sum(set[,"beta1"] >= (set[,"beta1_est_naive"] - 1.96*set[,"beta1_se_naive"]) & set[,"beta1"] <= (set[,"beta1_est_naive"] + 1.96*set[,"beta1_se_naive"]), na.rm=TRUE)/nsim*100,
	sum(set[,"beta2"] >= (set[,"beta2_est_naive"] - 1.96*set[,"beta2_se_naive"]) & set[,"beta2"] <= (set[,"beta2_est_naive"] + 1.96*set[,"beta2_se_naive"]), na.rm=TRUE)/nsim*100,
	sum(set[,"beta3"] >= (set[,"beta3_est_naive"] - 1.96*set[,"beta3_se_naive"]) & set[,"beta3"] <= (set[,"beta3_est_naive"] + 1.96*set[,"beta3_se_naive"]), na.rm=TRUE)/nsim*100,

	sum(set[,"beta0"] >= (set[,"beta0_est_prop"] - 1.96*set[,"beta0_se_prop"]) & set[,"beta0"] <= (set[,"beta0_est_prop"] + 1.96*set[,"beta0_se_prop"]), na.rm=TRUE)/nsim*100,
	sum(set[,"beta1"] >= (set[,"beta1_est_prop"] - 1.96*set[,"beta1_se_prop"]) & set[,"beta1"] <= (set[,"beta1_est_prop"] + 1.96*set[,"beta1_se_prop"]), na.rm=TRUE)/nsim*100,
	sum(set[,"beta2"] >= (set[,"beta2_est_prop"] - 1.96*set[,"beta2_se_prop"]) & set[,"beta2"] <= (set[,"beta2_est_prop"] + 1.96*set[,"beta2_se_prop"]), na.rm=TRUE)/nsim*100,
	sum(set[,"beta3"] >= (set[,"beta3_est_prop"] - 1.96*set[,"beta3_se_prop"]) & set[,"beta3"] <= (set[,"beta3_est_prop"] + 1.96*set[,"beta3_se_prop"]), na.rm=TRUE)/nsim*100)


	cvg_theta[i,] <- c(
	sum(set[,"theta0"] >= (set[,"theta0_est_naive"] - 1.96*set[,"theta0_se_naive"]) & set[,"theta0"] <= (set[,"theta0_est_naive"] + 1.96*set[,"theta0_se_naive"]), na.rm=TRUE)/nsim*100,
	sum(set[,"theta1"] >= (set[,"theta1_est_naive"] - 1.96*set[,"theta1_se_naive"]) & set[,"theta1"] <= (set[,"theta1_est_naive"] + 1.96*set[,"theta1_se_naive"]), na.rm=TRUE)/nsim*100,
	sum(set[,"theta2"] >= (set[,"theta2_est_naive"] - 1.96*set[,"theta2_se_naive"]) & set[,"theta2"] <= (set[,"theta2_est_naive"] + 1.96*set[,"theta2_se_naive"]), na.rm=TRUE)/nsim*100,
	sum(set[,"theta3"] >= (set[,"theta3_est_naive"] - 1.96*set[,"theta3_se_naive"]) & set[,"theta3"] <= (set[,"theta3_est_naive"] + 1.96*set[,"theta3_se_naive"]), na.rm=TRUE)/nsim*100,

	sum(set[,"theta0"] >= (set[,"theta0_est_prop"] - 1.96*set[,"theta0_se_prop"]) & set[,"theta0"] <= (set[,"theta0_est_prop"] + 1.96*set[,"theta0_se_prop"]), na.rm=TRUE)/nsim*100,
	sum(set[,"theta1"] >= (set[,"theta1_est_prop"] - 1.96*set[,"theta1_se_prop"]) & set[,"theta1"] <= (set[,"theta1_est_prop"] + 1.96*set[,"theta1_se_prop"]), na.rm=TRUE)/nsim*100,
	sum(set[,"theta2"] >= (set[,"theta2_est_prop"] - 1.96*set[,"theta2_se_prop"]) & set[,"theta2"] <= (set[,"theta2_est_prop"] + 1.96*set[,"theta2_se_prop"]), na.rm=TRUE)/nsim*100,
	sum(set[,"theta3"] >= (set[,"theta3_est_prop"] - 1.96*set[,"theta3_se_prop"]) & set[,"theta3"] <= (set[,"theta3_est_prop"] + 1.96*set[,"theta3_se_prop"]), na.rm=TRUE)/nsim*100)


	cvg_phi[i,] <- c(
	sum(set[,"phi0"] >= (set[,"phi0_est_naive"] - 1.96*set[,"phi0_se_naive"]) & set[,"phi0"] <= (set[,"phi0_est_naive"] + 1.96*set[,"phi0_se_naive"]), na.rm=TRUE)/nsim*100,
	sum(set[,"phi1"] >= (set[,"phi1_est_naive"] - 1.96*set[,"phi1_se_naive"]) & set[,"phi1"] <= (set[,"phi1_est_naive"] + 1.96*set[,"phi1_se_naive"]), na.rm=TRUE)/nsim*100,
	sum(set[,"phi2"] >= (set[,"phi2_est_naive"] - 1.96*set[,"phi2_se_naive"]) & set[,"phi2"] <= (set[,"phi2_est_naive"] + 1.96*set[,"phi2_se_naive"]), na.rm=TRUE)/nsim*100,
	sum(set[,"phi3"] >= (set[,"phi3_est_naive"] - 1.96*set[,"phi3_se_naive"]) & set[,"phi3"] <= (set[,"phi3_est_naive"] + 1.96*set[,"phi3_se_naive"]), na.rm=TRUE)/nsim*100,

	sum(set[,"phi0"] >= (set[,"phi0_est_prop"] - 1.96*set[,"phi0_se_prop"]) & set[,"phi0"] <= (set[,"phi0_est_prop"] + 1.96*set[,"phi0_se_prop"]), na.rm=TRUE)/nsim*100,
	sum(set[,"phi1"] >= (set[,"phi1_est_prop"] - 1.96*set[,"phi1_se_prop"]) & set[,"phi1"] <= (set[,"phi1_est_prop"] + 1.96*set[,"phi1_se_prop"]), na.rm=TRUE)/nsim*100,
	sum(set[,"phi2"] >= (set[,"phi2_est_prop"] - 1.96*set[,"phi2_se_prop"]) & set[,"phi2"] <= (set[,"phi2_est_prop"] + 1.96*set[,"phi2_se_prop"]), na.rm=TRUE)/nsim*100,
	sum(set[,"phi3"] >= (set[,"phi3_est_prop"] - 1.96*set[,"phi3_se_prop"]) & set[,"phi3"] <= (set[,"phi3_est_prop"] + 1.96*set[,"phi3_se_prop"]), na.rm=TRUE)/nsim*100)
}    	


colnames(cvg_beta) <- c("B0_n","B1_n","B2_n","B3_n","B0_p","B1_p","B2_p","B3_p")
colnames(cvg_theta) <- c("T0_n","T1_n","T2_n", "T2_p","T0_p","T1_p","T2_p", "T3_p")
colnames(cvg_phi) <- c("P0_n","P1_n","P2_n","P3_n","P0_p","P1_p","P2_p","P3_p")
	
dcvg_beta <- data.frame(cvg_beta, row.names=paste0("S", ord_scen))
dcvg_theta <- data.frame(cvg_theta, row.names=paste0("S", ord_scen))
dcvg_phi <- data.frame(cvg_phi, row.names=paste0("S", ord_scen))

print(xtable(dcvg_beta, type = "latex", digits=1), file = "t1_cvg_beta.tex", format.args = list(big.mark = " ", decimal.mark = "."))
print(xtable(dcvg_theta, type = "latex", digits=1), file = "t2_cvg_theta.tex", format.args = list(big.mark = " ", decimal.mark = "."))
print(xtable(dcvg_phi, type = "latex", digits=1), file = "t3_cvg_phi.tex", format.args = list(big.mark = " ", decimal.mark = "."))



########################################################################################
##
## WIDTH
##
########################################################################################
width_beta <- matrix(NA, nrow=16, ncol=8)
width_theta <- matrix(NA, nrow=16, ncol=8)
width_phi <- matrix(NA, nrow=16, ncol=8)

for (i in 1:length(lf)){
	set <- read.csv(paste0(path_in, "/Res_Scen_",ord_scen[i],".csv"))
	NOMICOL <- colnames(set)
	b0w_n <- mean((set[,"beta0_est_naive"]+1.96*set[,"beta0_se_naive"]) - (set[,"beta0_est_naive"]-1.96*set[,"beta0_se_naive"]), na.rm=TRUE)
	b1w_n <- mean((set[,"beta1_est_naive"]+1.96*set[,"beta1_se_naive"]) - (set[,"beta1_est_naive"]-1.96*set[,"beta1_se_naive"]), na.rm=TRUE)
	b2w_n <- mean((set[,"beta2_est_naive"]+1.96*set[,"beta2_se_naive"]) - (set[,"beta2_est_naive"]-1.96*set[,"beta2_se_naive"]), na.rm=TRUE)
	b3w_n <- mean((set[,"beta3_est_naive"]+1.96*set[,"beta3_se_naive"]) - (set[,"beta3_est_naive"]-1.96*set[,"beta3_se_naive"]), na.rm=TRUE)

	b0w_p <- mean((set[,"beta0_est_prop"]+1.96*set[,"beta0_se_prop"]) - (set[,"beta0_est_prop"]-1.96*set[,"beta0_se_prop"]), na.rm=TRUE)
	b1w_p <- mean((set[,"beta1_est_prop"]+1.96*set[,"beta1_se_prop"]) - (set[,"beta1_est_prop"]-1.96*set[,"beta1_se_prop"]), na.rm=TRUE)
	b2w_p <- mean((set[,"beta2_est_prop"]+1.96*set[,"beta2_se_prop"]) - (set[,"beta2_est_prop"]-1.96*set[,"beta2_se_prop"]), na.rm=TRUE)
	b3w_p <- mean((set[,"beta3_est_prop"]+1.96*set[,"beta3_se_prop"]) - (set[,"beta3_est_prop"]-1.96*set[,"beta3_se_prop"]), na.rm=TRUE)

	t0w_n <- mean((set[,"theta0_est_naive"]+1.96*set[,"theta0_se_naive"]) - (set[,"theta0_est_naive"]-1.96*set[,"theta0_se_naive"]), na.rm=TRUE)
	t1w_n <- mean((set[,"theta1_est_naive"]+1.96*set[,"theta1_se_naive"]) - (set[,"theta1_est_naive"]-1.96*set[,"theta1_se_naive"]), na.rm=TRUE)
	t2w_n <- mean((set[,"theta2_est_naive"]+1.96*set[,"theta2_se_naive"]) - (set[,"theta2_est_naive"]-1.96*set[,"theta2_se_naive"]), na.rm=TRUE)
	t3w_n <- mean((set[,"theta3_est_naive"]+1.96*set[,"theta3_se_naive"]) - (set[,"theta3_est_naive"]-1.96*set[,"theta3_se_naive"]), na.rm=TRUE)

	p0w_n <- mean((set[,"phi0_est_naive"]+1.96*set[,"phi0_se_naive"]) - (set[,"phi0_est_naive"]-1.96*set[,"phi0_se_naive"]), na.rm=TRUE)
	p1w_n <- mean((set[,"phi1_est_naive"]+1.96*set[,"phi1_se_naive"]) - (set[,"phi1_est_naive"]-1.96*set[,"phi1_se_naive"]), na.rm=TRUE)
	p2w_n <- mean((set[,"phi2_est_naive"]+1.96*set[,"phi2_se_naive"]) - (set[,"phi2_est_naive"]-1.96*set[,"phi2_se_naive"]), na.rm=TRUE)
	p3w_n <- mean((set[,"phi3_est_naive"]+1.96*set[,"phi3_se_naive"]) - (set[,"phi3_est_naive"]-1.96*set[,"phi3_se_naive"]), na.rm=TRUE)

	t0w_p <- mean((set[,"theta0_est_prop"]+1.96*set[,"theta0_se_prop"]) - (set[,"theta0_est_prop"]-1.96*set[,"theta0_se_prop"]), na.rm=TRUE)
	t1w_p <- mean((set[,"theta1_est_prop"]+1.96*set[,"theta1_se_prop"]) - (set[,"theta1_est_prop"]-1.96*set[,"theta1_se_prop"]), na.rm=TRUE)
	t2w_p <- mean((set[,"theta2_est_prop"]+1.96*set[,"theta2_se_prop"]) - (set[,"theta2_est_prop"]-1.96*set[,"theta2_se_prop"]), na.rm=TRUE)
	t3w_p <- mean((set[,"theta3_est_prop"]+1.96*set[,"theta3_se_prop"]) - (set[,"theta3_est_prop"]-1.96*set[,"theta3_se_prop"]), na.rm=TRUE)

	p0w_p <- mean((set[,"phi0_est_prop"]+1.96*set[,"phi0_se_prop"]) - (set[,"phi0_est_prop"]-1.96*set[,"phi0_se_prop"]), na.rm=TRUE)
	p1w_p <- mean((set[,"phi1_est_prop"]+1.96*set[,"phi1_se_prop"]) - (set[,"phi1_est_prop"]-1.96*set[,"phi1_se_prop"]), na.rm=TRUE)
	p2w_p <- mean((set[,"phi2_est_prop"]+1.96*set[,"phi2_se_prop"]) - (set[,"phi2_est_prop"]-1.96*set[,"phi2_se_prop"]), na.rm=TRUE)
	p3w_p <- mean((set[,"phi3_est_prop"]+1.96*set[,"phi3_se_prop"]) - (set[,"phi3_est_prop"]-1.96*set[,"phi3_se_prop"]), na.rm=TRUE)

	width_beta[i,] <- c(b0w_n, b1w_n, b2w_n, b3w_n, b0w_p, b1w_p, b2w_p, b3w_p) 
	width_theta[i,] <- c(t0w_n, t1w_n, t2w_n, t3w_n, t0w_p, t1w_p, t2w_p, t3w_p) 
	width_phi[i,] <- c(p0w_n, p1w_n, p2w_n, p3w_n, p0w_p, p1w_p, p2w_p, p3w_p) 
	}

colnames(width_beta) <- c("B0_n","B1_n","B2_n","B3_n","B0_p","B1_p","B2_p","B3_p")
colnames(width_theta) <- c("T0_n","T1_n","T2_n","T3_n","T0_p","T1_p","T2_p","T3_p")
colnames(width_phi) <- c("P0_n","P1_n","P2_n", "P3_n","P0_p","P1_p","P2_p", "P3_p")

dwidth_beta <- data.frame(width_beta, row.names=paste0("S", ord_scen))
dwidth_theta <- data.frame(width_theta, row.names=paste0("S", ord_scen))
dwidth_phi <- data.frame(width_phi, row.names=paste0("S", ord_scen))

print(xtable(dwidth_beta, type = "latex", digits=2), file = "t4_width_beta.tex", format.args = list(big.mark = " ", decimal.mark = "."))
print(xtable(dwidth_theta, type = "latex", digits=2), file = "t5_width_theta.tex", format.args = list(big.mark = " ", decimal.mark = "."))
print(xtable(dwidth_phi, type = "latex", digits=2), file = "t6_width_phi.tex", format.args = list(big.mark = " ", decimal.mark = "."))


