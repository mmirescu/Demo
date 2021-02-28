# ===== Model Incorporating Fuel Price and Availability Risks and Load =====

setwd("C:/Users/Mirescu1/ucloud/Universities/TU Wien/Theses/Dissertation/Eigene/Latex")

# === Install the Necessary Packages ===

library(xts)
library(quadprog)
library(rootSolve)
library(MASS)
library(maps)       # Provides functions that let us plot the maps
library(mapdata)    # Contains the hi-resolution points that mark out the countries.
library("colorspace")
library(grDevices)
library('nloptr')
library('NlcOptim')
library(rgl)
library('sp')
library('KernSmooth')
library(kernlab)
#library('alabama')
library(latex2exp)
library(tikzDevice)
library(extrafont)
library(Cairo)

# === Define Colors to Be Used Throughout This Document ===

# Call Up on Interface Where One Can Define a Personalised Color Palette With Subsequent Color Generator in HEX-Format

#choose_palette()

# Introduce the Personalised Palettes into R

# Nikolaus Rab (nr) color palette

col_nr<-rev(c("#E5FCC2","#9DE0AD","#45ADA8","#547980","#594F4F"))

# Blue-Green (bg) color palette

col_bg<-rev(c("#00584D","#00685E","#347B73","#65948E","#A1BAB6","#A3B9BC","#689398","#387A81","#00666E","#00565F"))

# Pink-Purple (pp) color palette 

col_pp<-rev(c("#CA2497","#D456A6","#DF7CB7","#EAA4CC","#F6D5E7","#F0D6EF","#DEA7DD","#D181CF","#C65BC3","#BD2ABA"))

# Coray-Pink (cp) color palette

col_cp<-rev(c("#FF7C89","#FF97A0","#FFB0B6","#FCC8CB","#F2E0E1","#F2DFE4","#FDC7D6","#FFADC8","#FF91BA","#FF73AD"))

# Newest color palette

col_new<-c("#005D00","#0B8955","#8BB59C","#E2E2E2","#83B4B9","#008991","#00626D")

colors<-col_bg

# === Availability Factor Wind ===

source("wind_code.R")


# === Technical and Costs Parameters ===

source("tech_costs_par.R")


# === Thermal Technologies ===

source("thermal_preparation.R")


# === System Costs ===

source("system_costs.R")

slope_RC<-0


# === Load ===

N<-5

source("load_opt.R")

# Solution Resulting from the Optimal Discretisation (Using Riemann Sums) of the LDC

d
LDC.d

# Demand per Block and Resulting Overall Demand for 1 MWh

D<-rep(0,N)
for (n in 1:(N-1)) D[n]<-d[2*n]*(LDC.d[2*n-1]-LDC.d[2*n+1])
D[N]<-LDC.d[2*N-1]
D<-rev(D/sum(D))

# Load Factors

LF<-rev(d[seq(2,length(d),2)])

# === Define the Objective Function, Its Gradient As Well As the Constraints and Their Gradient ===

# == Objective Function ==

# With System Costs

fn<-function(x,beta) {
  t(1/LF)%*%matrix(x,nrow=N,ncol=J+1)%*%t(t(c(F_t,F_r)))+t(rep(1,N))%*%matrix(x,nrow=N,ncol=J+1)%*%t(t(c(VOM_t+mean_t,(slope_BC+slope_RC)*(mean(avfac))^2*t(rep(1,N))%*%matrix(x,nrow=N,ncol=J+1)[,J+1])))+beta/2*(t(rep(1,N))%*%matrix(x,nrow=N,ncol=J+1)%*%cbind(rbind(var_t,rep(0,J)),c(rep(0,J),(slope_BC+slope_RC)^2*(mean(avfac))^2*var(avfac)*(rep(1,N)%*%matrix(x,nrow=N,ncol=J+1)[,J+1])^2))%*%t(matrix(x,nrow=N,ncol=J+1))%*%t(t(rep(1,N))))
}

mu<-function(x) t(1/LF)%*%matrix(x,nrow=N,ncol=J+1)%*%t(t(c(F_t,F_r)))+t(rep(1,N))%*%matrix(x,nrow=N,ncol=J+1)%*%t(t(c(VOM_t+mean_t,(slope_BC+slope_RC)*(mean(avfac))^2*t(rep(1,N))%*%matrix(x,nrow=N,ncol=J+1)[,J+1])))

sd<-function(x) sqrt((t(rep(1,N))%*%matrix(x,nrow=N,ncol=J+1)%*%cbind(rbind(var_t,rep(0,J)),c(rep(0,J),(slope_BC+slope_RC)^2*(mean(avfac))^2*var(avfac)*(rep(1,N)%*%matrix(x,nrow=N,ncol=J+1)[,J+1])^2))%*%t(matrix(x,nrow=N,ncol=J+1))%*%t(t(rep(1,N)))))

# Without System Costs

slope_BC_wo<-0
slope_RC_wo<-0

fn_wo<-function(x,beta) {
  t(1/LF)%*%matrix(x,nrow=N,ncol=J+1)%*%t(t(c(F_t,F_r)))+t(rep(1,N))%*%matrix(x,nrow=N,ncol=J+1)%*%t(t(c(VOM_t+mean_t,(slope_BC_wo+slope_RC_wo)*(mean(avfac))^2*t(rep(1,N))%*%matrix(x,nrow=N,ncol=J+1)[,J+1])))+beta/2*(t(rep(1,N))%*%matrix(x,nrow=N,ncol=J+1)%*%cbind(rbind(var_t,rep(0,J)),c(rep(0,J),(slope_BC_wo+slope_RC_wo)^2*(mean(avfac))^2*var(avfac)*(rep(1,N)%*%matrix(x,nrow=N,ncol=J+1)[,J+1])^2))%*%t(matrix(x,nrow=N,ncol=J+1))%*%t(t(rep(1,N))))
}

mu_wo<-function(x) t(1/LF)%*%matrix(x,nrow=N,ncol=J+1)%*%t(t(c(F_t,F_r)))+t(rep(1,N))%*%matrix(x,nrow=N,ncol=J+1)%*%t(t(c(VOM_t+mean_t,(slope_BC_wo+slope_RC_wo)*(mean(avfac))^2*t(rep(1,N))%*%matrix(x,nrow=N,ncol=J+1)[,J+1])))

sd_wo<-function(x) sqrt((t(rep(1,N))%*%matrix(x,nrow=N,ncol=J+1)%*%cbind(rbind(var_t,rep(0,J)),c(rep(0,J),(slope_BC_wo+slope_RC_wo)^2*(mean(avfac))^2*var(avfac)*(rep(1,N)%*%matrix(x,nrow=N,ncol=J+1)[,J+1])^2))%*%t(matrix(x,nrow=N,ncol=J+1))%*%t(t(rep(1,N)))))

# == Gradient Objective Function ==

# With System Costs

gr<-function(x,beta) {
  as.vector(as.matrix(1/LF)%*%c(F_t,F_r))+c(rep(VOM_t+mean_t,each=N),rep(2*(slope_BC+slope_RC)*(mean(avfac))^2*rep(1,N)%*%matrix(x,nrow=N,ncol=J+1)[,J+1],N))+beta*as.vector(matrix(rep(1,N*N),nrow=N,ncol=N)%*%matrix(x,nrow=N,ncol=J+1)%*%cbind(rbind(var_t,rep(0,J)),c(rep(0,J),2*(slope_BC+slope_RC)^2*(mean(avfac))^2*var(avfac)*(rep(1,N)%*%matrix(x,nrow=N,ncol=J+1)[,J+1])^2)))
}

# Without System Costs

gr_wo<-function(x,beta) {
  as.vector(as.matrix(1/LF)%*%c(F_t,F_r))+c(rep(VOM_t+mean_t,each=N),rep(2*(slope_BC_wo+slope_RC_wo)*(mean(avfac))^2*rep(1,N)%*%matrix(x,nrow=N,ncol=J+1)[,J+1],N))+beta*as.vector(matrix(rep(1,N*N),nrow=N,ncol=N)%*%matrix(x,nrow=N,ncol=J+1)%*%cbind(rbind(var_t,rep(0,J)),c(rep(0,J),2*(slope_BC_wo+slope_RC_wo)^2*(mean(avfac))^2*var(avfac)*(rep(1,N)%*%matrix(x,nrow=N,ncol=J+1)[,J+1])^2)))
}

# == Equalities ==

# Actual

heq<-function(x) as.vector(matrix(x,nrow=N,ncol=J+1)%*%c(rep(1,J),mean(avfac))-D) 

# Jacobian

heqjac<-function(x) cbind(t(apply(diag(1,N,N),1,rep,J)),diag(mean(avfac),N,N)) 

# === Define Vector of Risk Aversion Parameters ===

beta_max<-10
beta<-seq(0,beta_max,0.01)


#############################################################################################################################################

# ===== All 4 Technologies: Gas, Coal, Nuclear and Wind =====

J<-3

# === Initialize Lists, Vectors and Matrices for the Storing of Results ===

# With System Costs

res_all<-vector(mode="list",length=length(beta))
x_all<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_all<-rep(0,length(beta))
mu_all<-rep(0,length(beta))
sd_all<-rep(0,length(beta))
it_all<-rep(0,length(beta))
status_all<-rep(0,length(beta))

# Without System Costs

res_all_wo<-vector(mode="list",length=length(beta))
x_all_wo<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_all_wo<-rep(0,length(beta))
mu_all_wo<-rep(0,length(beta))
sd_all_wo<-rep(0,length(beta))
it_all_wo<-rep(0,length(beta))
status_all_wo<-rep(0,length(beta))


# === Set Starting Value ===

x0<-rep(0.5,N*(J+1))
x0_wo<-rep(0.5,N*(J+1))

# === Compute Solution ===

for (i in 1:length(beta)) {
  f<-function(y) fn(y,beta[i])
  g<-function(y) gr(y,beta[i])
  res_all[[i]]<-nloptr(x0=x0,eval_f=f,eval_grad_f=g,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_all[i,]<-round(res_all[[i]]$solution,digits=4)
  OF_all[i]<-res_all[[i]]$objective
  it_all[i]<-res_all[[i]]$iterations
  status_all[i]<-res_all[[i]]$status
#  x0<-x_all[i,]
  ########################################################################################################################################################################################################## 
  f_wo<-function(y) fn_wo(y,beta[i])
  g_wo<-function(y) gr_wo(y,beta[i])
  res_all_wo[[i]]<-nloptr(x0=x0_wo,eval_f=f_wo,eval_grad_f=g_wo,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_all_wo[i,]<-round(res_all_wo[[i]]$solution,digits=4)
  OF_all_wo[i]<-res_all_wo[[i]]$objective
  it_all_wo[i]<-res_all_wo[[i]]$iterations
  status_all_wo[i]<-res_all_wo[[i]]$status
#  x0_wo<-x_all_wo[i,]
}

c(max(status_all),min(status_all),max(status_all_wo),min(status_all_wo))


#View(cbind(beta,x_all,x_all_wo))

# === Compute the Relative Shares of Each Technology in Each Load Block ===

x_rel_all<-round(sweep(x_all,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)
x_rel_all_wo<-round(sweep(x_all_wo,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)

#View(cbind(beta,x_rel_all,beta,x_rel_all_wo))

# === Evaluate MU and SD for Optimal Values ===

for (i in 1:length(beta)) {
  mu_all[i]<-mu(x_all[i,])
  sd_all[i]<-sd(x_all[i,])
  mu_all_wo[i]<-mu_wo(x_all_wo[i,])
  sd_all_wo[i]<-sd_wo(x_all_wo[i,])
}
#View(cbind(beta,mu_all,mu_all_wo,sd_all,sd_all_wo))
#View(cbind(beta,mu_all+beta/2*sd_all^2,OF_all))

# === Cumulative Weights per Technology (Aggregated Over All Blocks) ===

# With and Without System Costs

x_cum_all<-matrix(0,nrow=length(beta),ncol=J+1)
x_cum_all_wo<-matrix(0,nrow=length(beta),ncol=J+1)
for (i in 1:length(beta)) {
  x_cum_all[i,]<-c(sapply(1:J,function(u) sum(x_all[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_all[i,(N*J+1):(N*(J+1))]))
  x_cum_all_wo[i,]<-c(sapply(1:J,function(u) sum(x_all_wo[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_all_wo[i,(N*J+1):(N*(J+1))]))
}

#View(cbind(beta,x_cum_all,x_cum_all_wo))

# === Plots ===

# == Efficient Frontier, Expected Costs and Standard Deviation ==

# R

par(mar=c(2,2,2,2))
layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))

plot(sd_all,mu_all,type="l",col="darkblue",lwd=3,xlim=c(min(sd_all,sd_all_wo),max(sd_all,sd_all_wo)),ylim=c(min(mu_all,mu_all_wo),max(mu_all,mu_all_wo)),xlab=TeX('$\\sigma$'),ylab=TeX('$\\mu$'),main="Efficient Frontier")
lines(sd_all_wo,mu_all_wo,col="lightblue",lwd=3,type="l",lty=2)
legend("topright",c("+","-"),col=c("darkblue","lightblue"),lwd=3)

plot(beta,mu_all,type="l",col="darkblue",lwd=3,ylim=c(min(mu_all,mu_all_wo),max(mu_all,mu_all_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\mu$'),main="Expected Costs")
lines(beta,mu_all_wo,col="lightblue",lwd=3,type="l",lty=2)
legend("bottomright",c("+","-"),col=c("darkblue","lightblue"),lwd=3)

plot(beta,sd_all,type="l",col="darkblue",lwd=3,ylim=c(min(sd_all,sd_all_wo),max(sd_all,sd_all_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\sigma$'),main="Standard Deviation")
lines(beta,sd_all_wo,col="lightblue",lwd=3,type="l",lty=2)
legend("topright",c("+","-"),col=c("darkblue","lightblue"),lwd=3)

# PDF

# pdf("EF_MU_SD_all.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))
# 
# plot(sd_all,mu_all,type="l",col="darkblue",lwd=3,xlim=c(min(sd_all,sd_all_wo),max(sd_all,sd_all_wo)),ylim=c(min(mu_all,mu_all_wo),max(mu_all,mu_all_wo)),xlab=TeX('$\\sigma$'),ylab=TeX('$\\mu$'),main="Efficient Frontier")
# lines(sd_all_wo,mu_all_wo,col="lightblue",lwd=3,type="l",lty=2)
# legend("topright",c("+","-"),col=c("darkblue","lightblue"),lwd=3)
# 
# plot(beta,mu_all,type="l",col="darkblue",lwd=3,ylim=c(min(mu_all,mu_all_wo),max(mu_all,mu_all_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\mu$'),main="Expected Costs")
# lines(beta,mu_all_wo,col="lightblue",lwd=3,type="l",lty=2)
# legend("bottomright",c("+","-"),col=c("darkblue","lightblue"),lwd=3)
# 
# plot(beta,sd_all,type="l",col="darkblue",lwd=3,ylim=c(min(sd_all,sd_all_wo),max(sd_all,sd_all_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\sigma$'),main="Standard Deviation")
# lines(beta,sd_all_wo,col="lightblue",lwd=3,type="l",lty=2)
# legend("topright",c("+","-"),col=c("darkblue","lightblue"),lwd=3)
# dev.off()


# == Shares on the EF == 

# R

par(mar=c(2,2,2,2))
par(mfrow=c(1,2))
plot(beta,x_cum_all[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab="beta",ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="With System Costs")
lines(beta,x_cum_all[,1]+x_cum_all[,2],type="l")
lines(beta,x_cum_all[,1]+x_cum_all[,2]+x_cum_all[,3],type="l")
lines(beta,x_cum_all[,1]+x_cum_all[,2]+x_cum_all[,3]+x_cum_all[,4],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_all[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_all[,1],rev(x_cum_all[,1]+x_cum_all[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_all[,1]+x_cum_all[,2],rev(x_cum_all[,1]+x_cum_all[,2]+x_cum_all[,3])),col=adjustcolor(colorRampPalette(colors)(4)[3]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_all[,1]+x_cum_all[,2]+x_cum_all[,3],rev(x_cum_all[,1]+x_cum_all[,2]+x_cum_all[,3]+x_cum_all[,4])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.035,labels="gas",cex=1.2,col="black",pos=4)
text(1.2,0.036,labels="coal",cex=1.2,col="white",pos=4)
text(3.5,0.4,labels="nuclear",cex=1.2,col="white",pos=4)
text(6.5,0.83,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_all_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab="beta",ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Without System Costs")
lines(beta,x_cum_all_wo[,1]+x_cum_all_wo[,2],type="l")
lines(beta,x_cum_all_wo[,1]+x_cum_all_wo[,2]+x_cum_all_wo[,3],type="l")
lines(beta,x_cum_all_wo[,1]+x_cum_all_wo[,2]+x_cum_all_wo[,3]+x_cum_all_wo[,4],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_all_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_all_wo[,1],rev(x_cum_all_wo[,1]+x_cum_all_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_all_wo[,1]+x_cum_all_wo[,2],rev(x_cum_all_wo[,1]+x_cum_all_wo[,2]+x_cum_all_wo[,3])),col=adjustcolor(colorRampPalette(colors)(4)[3]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_all_wo[,1]+x_cum_all_wo[,2]+x_cum_all_wo[,3],rev(x_cum_all_wo[,1]+x_cum_all_wo[,2]+x_cum_all_wo[,3]+x_cum_all_wo[,4])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.035,labels="gas",cex=1.2,col="black",pos=4)
text(1.2,0.036,labels="coal",cex=1.2,col="white",pos=4)
text(2.5,0.15,labels="nuclear",cex=1.2,col="white",pos=4)
text(5.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

# # PDF
#  
# pdf("Shares_EF_all.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# par(mfrow=c(1,2))
# plot(beta,x_cum_all[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="With System Costs")
# lines(beta,x_cum_all[,1]+x_cum_all[,2],type="l")
# lines(beta,x_cum_all[,1]+x_cum_all[,2]+x_cum_all[,3],type="l")
# lines(beta,x_cum_all[,1]+x_cum_all[,2]+x_cum_all[,3]+x_cum_all[,4],type="l")
# polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_all[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_all[,1],rev(x_cum_all[,1]+x_cum_all[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_all[,1]+x_cum_all[,2],rev(x_cum_all[,1]+x_cum_all[,2]+x_cum_all[,3])),col=adjustcolor(colorRampPalette(colors)(4)[3]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_all[,1]+x_cum_all[,2]+x_cum_all[,3],rev(x_cum_all[,1]+x_cum_all[,2]+x_cum_all[,3]+x_cum_all[,4])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
# abline(h=0)
# abline(h=1)
# text(-0.05,0.035,labels="gas",cex=1.0,col="black",pos=4)
# text(1.2,0.036,labels="coal",cex=1.2,col="white",pos=4)
# text(3.5,0.4,labels="nuclear",cex=1.2,col="white",pos=4)
# text(6.5,0.83,labels="wind",cex=1.2,col="white",pos=4)
# 
# plot(beta,x_cum_all_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Without System Costs")
# lines(beta,x_cum_all_wo[,1]+x_cum_all_wo[,2],type="l")
# lines(beta,x_cum_all_wo[,1]+x_cum_all_wo[,2]+x_cum_all_wo[,3],type="l")
# lines(beta,x_cum_all_wo[,1]+x_cum_all_wo[,2]+x_cum_all_wo[,3]+x_cum_all_wo[,4],type="l")
# polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_all_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_all_wo[,1],rev(x_cum_all_wo[,1]+x_cum_all_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_all_wo[,1]+x_cum_all_wo[,2],rev(x_cum_all_wo[,1]+x_cum_all_wo[,2]+x_cum_all_wo[,3])),col=adjustcolor(colorRampPalette(colors)(4)[3]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_all_wo[,1]+x_cum_all_wo[,2]+x_cum_all_wo[,3],rev(x_cum_all_wo[,1]+x_cum_all_wo[,2]+x_cum_all_wo[,3]+x_cum_all_wo[,4])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
# abline(h=0)
# abline(h=1)
# text(-0.05,0.035,labels="gas",cex=1.0,col="black",pos=4)
# text(1.9,0.036,labels="coal",cex=1.2,col="white",pos=4)
# text(2.5,0.15,labels="nuclear",cex=1.2,col="white",pos=4)
# text(5.5,0.6,labels="wind",cex=1.2,col="white",pos=4)
# dev.off()

# == Shares in Each Load Block ==

# R

par(mar=c(2,2,2,2))
par(mfrow=c(2,2))
plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 0$ With'),cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Nuclear"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 0$ Without'),cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Nuclear"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 10$ With'),cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Nuclear","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 10$ Without'),cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Coal","Nuclear","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[3]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)


# PDF

# CairoPDF("Shares_LB_all.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# par(mfrow=c(2,2))
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 0$ With'),cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Nuclear"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3])),lwd=4)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 0$ Without'),cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Nuclear"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3])),lwd=4)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 10$ With'),cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Nuclear","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 10$ Without'),cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Coal","Nuclear","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[3]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
# dev.off()

######################################################################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################################################################


# ===== 3 Technologies: Gas, Coal and Wind =====

J<-2

IC_t<-IC_t[1:2]
FOM_t<-FOM_t[1:2]
F_t<-F_t[1:2]
VOM_t<-VOM_t[1:2]

mean_t<-mean_t[1:2]
var_t<-var_t[1:2,1:2]

# === Initialize Lists, Vectors and Matrices for the Storing of Results ===

# With System Costs

res_GCW<-vector(mode="list",length=length(beta))
x_GCW<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_GCW<-rep(0,length(beta))
mu_GCW<-rep(0,length(beta))
sd_GCW<-rep(0,length(beta))
it_GCW<-rep(0,length(beta))
status_GCW<-rep(0,length(beta))

# Without System Costs

res_GCW_wo<-vector(mode="list",length=length(beta))
x_GCW_wo<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_GCW_wo<-rep(0,length(beta))
mu_GCW_wo<-rep(0,length(beta))
sd_GCW_wo<-rep(0,length(beta))
it_GCW_wo<-rep(0,length(beta))
status_GCW_wo<-rep(0,length(beta))


# === Set Starting Value ===

x0<-rep(0.5,N*(J+1))
x0_wo<-rep(0.5,N*(J+1))

# === Compute Solution ===

for (i in 1:length(beta)) {
  f<-function(y) fn(y,beta[i])
  g<-function(y) gr(y,beta[i])
  res_GCW[[i]]<-nloptr(x0=x0,eval_f=f,eval_grad_f=g,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_GCW[i,]<-round(res_GCW[[i]]$solution,digits=4)
  OF_GCW[i]<-res_GCW[[i]]$objective
  it_GCW[i]<-res_GCW[[i]]$iterations
  status_GCW[i]<-res_GCW[[i]]$status
#  x0<-x_GCW[i,]
  ########################################################################################################################################################################################################## 
  f_wo<-function(y) fn_wo(y,beta[i])
  g_wo<-function(y) gr_wo(y,beta[i])
  res_GCW_wo[[i]]<-nloptr(x0=x0_wo,eval_f=f_wo,eval_grad_f=g_wo,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_GCW_wo[i,]<-round(res_GCW_wo[[i]]$solution,digits=4)
  OF_GCW_wo[i]<-res_GCW_wo[[i]]$objective
  it_GCW_wo[i]<-res_GCW_wo[[i]]$iterations
  status_GCW_wo[i]<-res_GCW_wo[[i]]$status
#  x0_wo<-x_GCW_wo[i,]
}

c(max(status_GCW),min(status_GCW),max(status_GCW_wo),min(status_GCW_wo))


#View(cbind(beta,x_GCW,x_GCW_wo))

# === Compute the Relative Shares of Each Technology in Each Load Block ===

x_rel_GCW<-round(sweep(x_GCW,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)
x_rel_GCW_wo<-round(sweep(x_GCW_wo,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)

#View(cbind(beta,x_rel_GCW,beta,x_rel_GCW_wo))

# === Evaluate MU and SD for Optimal Values ===

for (i in 1:length(beta)) {
  mu_GCW[i]<-mu(x_GCW[i,])
  sd_GCW[i]<-sd(x_GCW[i,])
  mu_GCW_wo[i]<-mu_wo(x_GCW_wo[i,])
  sd_GCW_wo[i]<-sd_wo(x_GCW_wo[i,])
}
#View(cbind(beta,mu_GCW,mu_GCW_wo,sd_GCW,sd_GCW_wo))
#View(cbind(beta,mu_GCW+beta/2*sd_GCW^2,OF_GCW))

# === Cumulative Weights per Technology (Aggregated Over All Blocks) ===

# With and Without System Costs

x_cum_GCW<-matrix(0,nrow=length(beta),ncol=J+1)
x_cum_GCW_wo<-matrix(0,nrow=length(beta),ncol=J+1)
for (i in 1:length(beta)) {
  x_cum_GCW[i,]<-c(sapply(1:J,function(u) sum(x_GCW[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_GCW[i,(N*J+1):(N*(J+1))]))
  x_cum_GCW_wo[i,]<-c(sapply(1:J,function(u) sum(x_GCW_wo[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_GCW_wo[i,(N*J+1):(N*(J+1))]))
}

#View(cbind(beta,x_cum_GCW,x_cum_GCW_wo))

# === Plots ===

# == Efficient Frontier, Expected Costs, Standard Deviation ==

# R

par(mar=c(2,2,2,2))
layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))
plot(sd_GCW,mu_GCW,type="l",col="darkblue",lwd=3,xlim=c(min(sd_GCW,sd_GCW_wo),max(sd_GCW,sd_GCW_wo)),ylim=c(min(mu_GCW,mu_GCW_wo),max(mu_GCW,mu_GCW_wo)),xlab=TeX('$\\sigma$'),ylab=TeX('$\\mu$'),main="Efficient Frontier")
lines(sd_GCW_wo,mu_GCW_wo,col="lightblue",lwd=3,type="l",lty=2)
legend("topright",c("+","-"),col=c("darkblue","lightblue"),lwd=3)

plot(beta,mu_GCW,type="l",col="darkblue",lwd=3,ylim=c(min(mu_GCW,mu_GCW_wo),max(mu_GCW,mu_GCW_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\mu$'),main="Expected Costs")
lines(beta,mu_GCW_wo,col="lightblue",lwd=3,type="l",lty=2)
legend("bottomright",c("+","-"),col=c("darkblue","lightblue"),lwd=3)

plot(beta,sd_GCW,type="l",col="darkblue",lwd=3,ylim=c(min(sd_GCW,sd_GCW_wo),max(sd_GCW,sd_GCW_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\sigma$'),main="Standard Deviation")
lines(beta,sd_GCW_wo,col="lightblue",lwd=3,type="l",lty=2)
legend("topright",c("+","-"),col=c("darkblue","lightblue"),lwd=3)

# # PDF
# 
# pdf("EF_MU_SD_GCW.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))
# 
# plot(sd_GCW,mu_GCW,type="l",col="darkblue",lwd=3,xlim=c(min(sd_GCW,sd_GCW_wo),max(sd_GCW,sd_GCW_wo)),ylim=c(min(mu_GCW,mu_GCW_wo),max(mu_GCW,mu_GCW_wo)),xlab=TeX('$\\sigma$'),ylab=TeX('$\\mu$'),main="Efficient Frontier")
# lines(sd_GCW_wo,mu_GCW_wo,col="lightblue",lwd=3,type="l",lty=2)
# legend("topright",c("+","-"),col=c("darkblue","lightblue"),lwd=3)
# 
# plot(beta,mu_GCW,type="l",col="darkblue",lwd=3,ylim=c(min(mu_GCW,mu_GCW_wo),max(mu_GCW,mu_GCW_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\mu$'),main="Expected Costs")
# lines(beta,mu_GCW_wo,col="lightblue",lwd=3,type="l",lty=2)
# legend("topleft",c("+","-"),col=c("darkblue","lightblue"),lwd=3)
# 
# plot(beta,sd_GCW,type="l",col="darkblue",lwd=3,ylim=c(min(sd_GCW,sd_GCW_wo),max(sd_GCW,sd_GCW_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\sigma$'),main="Standard Deviation")
# lines(beta,sd_GCW_wo,col="lightblue",lwd=3,type="l",lty=2)
# legend("topright",c("+","-"),col=c("darkblue","lightblue"),lwd=3)
# 
# dev.off()


# == Shares on the EF == 

# R

par(mar=c(2,2,2,2))
par(mfrow=c(1,2))
plot(beta,x_cum_GCW[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="With System Costs")
lines(beta,x_cum_GCW[,1]+x_cum_GCW[,2],type="l")
lines(beta,x_cum_GCW[,1]+x_cum_GCW[,2]+x_cum_GCW[,3],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_GCW[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_GCW[,1],rev(x_cum_GCW[,1]+x_cum_GCW[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_GCW[,1]+x_cum_GCW[,2],rev(x_cum_GCW[,1]+x_cum_GCW[,2]+x_cum_GCW[,3])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0.05,0.035,labels="gas",cex=1.2,col="black",pos=4)
text(2.5,0.1,labels="coal",cex=1.2,col="white",pos=4)
text(5.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_GCW_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Without System Costs")
lines(beta,x_cum_GCW_wo[,1]+x_cum_GCW_wo[,2],type="l")
lines(beta,x_cum_GCW_wo[,1]+x_cum_GCW_wo[,2]+x_cum_GCW_wo[,3],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_GCW_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_GCW_wo[,1],rev(x_cum_GCW_wo[,1]+x_cum_GCW_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_GCW_wo[,1]+x_cum_GCW_wo[,2],rev(x_cum_GCW_wo[,1]+x_cum_GCW_wo[,2]+x_cum_GCW_wo[,3])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0.05,0.035,labels="gas",cex=1.2,col="black",pos=4)
text(0.4,0.1,labels="coal",cex=1.2,col="white",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

# PDF

# pdf("Shares_EF_GCW.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# par(mfrow=c(1,2))
# plot(beta,x_cum_GCW[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="With System Costs")
# lines(beta,x_cum_GCW[,1]+x_cum_GCW[,2],type="l")
# lines(beta,x_cum_GCW[,1]+x_cum_GCW[,2]+x_cum_GCW[,3],type="l")
# polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_GCW[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_GCW[,1],rev(x_cum_GCW[,1]+x_cum_GCW[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_GCW[,1]+x_cum_GCW[,2],rev(x_cum_GCW[,1]+x_cum_GCW[,2]+x_cum_GCW[,3])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
# abline(h=0)
# abline(h=1)
# text(0.05,0.035,labels="gas",cex=1.2,col="black",pos=4)
# text(2.5,0.1,labels="coal",cex=1.2,col="white",pos=4)
# text(5.5,0.6,labels="wind",cex=1.2,col="white",pos=4)
# 
# plot(beta,x_cum_GCW_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Without System Costs")
# lines(beta,x_cum_GCW_wo[,1]+x_cum_GCW_wo[,2],type="l")
# lines(beta,x_cum_GCW_wo[,1]+x_cum_GCW_wo[,2]+x_cum_GCW_wo[,3],type="l")
# polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_GCW_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_GCW_wo[,1],rev(x_cum_GCW_wo[,1]+x_cum_GCW_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_GCW_wo[,1]+x_cum_GCW_wo[,2],rev(x_cum_GCW_wo[,1]+x_cum_GCW_wo[,2]+x_cum_GCW_wo[,3])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
# abline(h=0)
# abline(h=1)
# text(0.05,0.035,labels="gas",cex=1.2,col="black",pos=4)
# text(0.4,0.1,labels="coal",cex=1.2,col="white",pos=4)
# text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)
# dev.off()

# == Shares in Each LB == 

# R

par(mar=c(2,2,2,2))
par(mfrow=c(2,2))
plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 0$ With'),cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Coal","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 0$ Without'),cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Coal","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 10$ With'),cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-3]+x_rel_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3]+x_rel_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Coal","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 10$ Without'),cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Coal","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)


# # PDF
# 
# CairoPDF("Shares_LB_GCW.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# par(mfrow=c(2,2))
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 0$ With'),cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Coal","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 0$ Without'),cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Coal","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 10$ With'),cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-3]+x_rel_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3]+x_rel_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Coal","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 10$ Without'),cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Coal","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
# dev.off()

######################################################################################################################################################################################################################################################################################################################################################################################################
######################################################################################################################################################################################################################################################################################################################################################################################################


# ===== 2 Technologies: Gas and Wind =====

J<-1

IC_t<-IC_t[1]
FOM_t<-FOM_t[1]
F_t<-F_t[1]
VOM_t<-VOM_t[1]

mean_t<-mean_t[1]
var_t<-var_t[1,1]

# === Initialize Lists, Vectors and Matrices for the Storing of Results ===

# With System Costs

res_GW<-vector(mode="list",length=length(beta))
x_GW<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_GW<-rep(0,length(beta))
mu_GW<-rep(0,length(beta))
sd_GW<-rep(0,length(beta))
it_GW<-rep(0,length(beta))
status_GW<-rep(0,length(beta))

# Without System Costs

res_GW_wo<-vector(mode="list",length=length(beta))
x_GW_wo<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_GW_wo<-rep(0,length(beta))
mu_GW_wo<-rep(0,length(beta))
sd_GW_wo<-rep(0,length(beta))
it_GW_wo<-rep(0,length(beta))
status_GW_wo<-rep(0,length(beta))

# === Set Starting Value ===

x0<-rep(0.5,N*(J+1))
x0_wo<-rep(0.5,N*(J+1))

# Compute Results for beta in [0,beta_max]

for (i in 1:length(beta)) {
  f<-function(y) fn(y,beta[i])
  g<-function(y) gr(y,beta[i])
  res_GW[[i]]<-nloptr(x0=x0,eval_f=f,eval_grad_f=g,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4,"maxeval"=1e6))
  x_GW[i,]<-round(res_GW[[i]]$solution,digits=4)
  OF_GW[i]<-res_GW[[i]]$objective
  it_GW[i]<-res_GW[[i]]$iterations
  status_GW[i]<-res_GW[[i]]$status
#  x0<-x_GW[i,]
  ########################################################################################################################################################################################################## 
  f_wo<-function(y) fn_wo(y,beta[i])
  g_wo<-function(y) gr_wo(y,beta[i])
  res_GW_wo[[i]]<-nloptr(x0=x0_wo,eval_f=f_wo,eval_grad_f=g_wo,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4,"maxeval"=1e6))
  x_GW_wo[i,]<-round(res_GW_wo[[i]]$solution,digits=4)
  OF_GW_wo[i]<-res_GW_wo[[i]]$objective
  it_GW_wo[i]<-res_GW_wo[[i]]$iterations
  status_GW_wo[i]<-res_GW_wo[[i]]$status
#  x0_wo<-x_GW_wo[i,]
}

c(max(status_GW),min(status_GW),max(status_GW_wo),min(status_GW_wo))

# === Compute the Relative Shares of Each Technology in Each Load Block ===

x_rel_GW<-round(sweep(x_GW,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)
x_rel_GW_wo<-round(sweep(x_GW_wo,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)

# === Evaluate MU and SD for Optimal Values ===

for (i in 1:length(beta)) {
  mu_GW[i]<-mu(x_GW[i,])
  sd_GW[i]<-sd(x_GW[i,])
  mu_GW_wo[i]<-mu_wo(x_GW_wo[i,])
  sd_GW_wo[i]<-sd_wo(x_GW_wo[i,])
}

# === Cumulative Weights per Technology (Aggregated Over All Blocks) ===

# With System Costs

x_cum_GW<-matrix(0,nrow=length(beta),ncol=J+1)
x_cum_GW_wo<-matrix(0,nrow=length(beta),ncol=J+1)
for (i in 1:length(beta)) {
  x_cum_GW[i,]<-c(sapply(1:J,function(u) sum(x_GW[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_GW[i,(N*J+1):(N*(J+1))]))
  x_cum_GW_wo[i,]<-c(sapply(1:J,function(u) sum(x_GW_wo[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_GW_wo[i,(N*J+1):(N*(J+1))]))
}

# === Plots ===

# == Efficient Frontier, Expected Costs and Standard Deviation ==

# R

par(mar=c(2,2,2,2))
layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))

plot(sd_GW,mu_GW,type="l",col="darkblue",lwd=3,xlim=c(min(sd_GW,sd_GW_wo),max(sd_GW,sd_GW_wo)),ylim=c(min(mu_GW,mu_GW_wo),max(mu_GW,mu_GW_wo)),xlab=TeX('$\\sigma$'),ylab=TeX('$\\mu$'),main="Efficient Frontier")
lines(sd_GW_wo,mu_GW_wo,col="lightblue",lwd=3,type="l",lty=2)
legend("topright",c("+","-"),col=c("darkblue","lightblue"),lwd=3)

plot(beta,mu_GW,type="l",col="darkblue",lwd=3,ylim=c(min(mu_GW,mu_GW_wo),max(mu_GW,mu_GW_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\mu$'),main="Expected Costs")
lines(beta,mu_GW_wo,col="lightblue",lwd=2.5,type="l",lty=2)
legend("topleft",c("+","-"),col=c("darkblue","lightblue"),lwd=3)

plot(beta,sd_GW,type="l",col="darkblue",lwd=3,ylim=c(min(sd_GW,sd_GW_wo),max(sd_GW,sd_GW_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\sigma$'),main="Standard Deviation")
lines(beta,sd_GW_wo,col="lightblue",lwd=3,type="l",lty=2)
legend("topright",c("+","-"),col=c("darkblue","lightblue"),lwd=3)


# PDF

# pdf("EF_MU_SD_GW.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))
# plot(sd_GW,mu_GW,type="l",col="darkblue",lwd=3,xlim=c(min(sd_GW,sd_GW_wo),max(sd_GW,sd_GW_wo)),ylim=c(min(mu_GW,mu_GW_wo),max(mu_GW,mu_GW_wo)),xlab=TeX('$\\sigma$'),ylab=TeX('$\\mu$'),main="Efficient Frontier")
# lines(sd_GW_wo,mu_GW_wo,col="lightblue",lwd=3,type="l",lty=2)
# legend("topright",c("+","-"),col=c("darkblue","lightblue"),lwd=3)
# 
# plot(beta,mu_GW,type="l",col="darkblue",lwd=3,ylim=c(min(mu_GW,mu_GW_wo),max(mu_GW,mu_GW_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\mu$'),main="Expected Costs")
# lines(beta,mu_GW_wo,col="lightblue",lwd=3,type="l",lty=2)
# legend("topleft",c("+","-"),col=c("darkblue","lightblue"),lwd=3)
# 
# plot(beta,sd_GW,type="l",col="darkblue",lwd=3,ylim=c(min(sd_GW,sd_GW_wo),max(sd_GW,sd_GW_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\sigma$'),main="Standard Deviation")
# lines(beta,sd_GW_wo,col="lightblue",lwd=3,type="l",lty=2)
# legend("topright",c("+","-"),col=c("darkblue","lightblue"),lwd=3)
# dev.off()

# == Shares on the EF ==

# R

par(mar=c(2,2,2,2))
par(mfrow=c(1,2))
plot(beta,x_cum_GW[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="With System Costs")
lines(beta,x_cum_GW[,1]+x_cum_GW[,2],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_GW[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_GW[,1],rev(x_cum_GW[,1]+x_cum_GW[,2])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0.5,0.08,labels="gas",cex=1.2,col="black",pos=4)
text(4.3,0.55,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_GW_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Without System Costs")
lines(beta,x_cum_GW_wo[,1]+x_cum_GW_wo[,2],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_GW_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_GW_wo[,1],rev(x_cum_GW_wo[,1]+x_cum_GW_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0.2,0.05,labels="gas",cex=1.2,col="black",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

# PDF
 
# pdf("Shares_EF_GW.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# par(mfrow=c(1,2))
# plot(beta,x_cum_GW[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="With System Costs")
# lines(beta,x_cum_GW[,1]+x_cum_GW[,2],type="l")
# polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_GW[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_GW[,1],rev(x_cum_GW[,1]+x_cum_GW[,2])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
# abline(h=0)
# abline(h=1)
# text(0.5,0.08,labels="gas",cex=1.2,col="black",pos=4)
# text(4.3,0.55,labels="wind",cex=1.2,col="white",pos=4)
# 
# plot(beta,x_cum_GW_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Without System Costs")
# lines(beta,x_cum_GW_wo[,1]+x_cum_GW_wo[,2],type="l")
# polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_GW_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_GW_wo[,1],rev(x_cum_GW_wo[,1]+x_cum_GW_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
# abline(h=0)
# abline(h=1)
# text(0.2,0.05,labels="gas",cex=1.2,col="black",pos=4)
# text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)
# dev.off()


# == Shares in Each LB == 

# R

par(mar=c(2,2,2,2))
par(mfrow=c(2,2))
plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 0$ With'),cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 0$ Without'),cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 10$ With'),cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-3]+x_rel_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3]+x_rel_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 10$ Without'),cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_GW_wo[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_GW_wo[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_GW_wo[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_GW_wo[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)


# PDF

# CairoPDF("Shares_LB_GW.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# par(mfrow=c(2,2))
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 0$ With'),cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 0$ Without'),cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 10$ With'),cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-3]+x_rel_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3]+x_rel_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main=TeX('$\\beta = 10$ Without'),cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_GW_wo[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_GW_wo[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_GW_wo[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_GW_wo[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
# dev.off()

# ======= Overall Plots =======

# == Efficient Frontier, Expected Costs and Standard Deviation ==

# R

par(mar=c(2,2,2,2))
layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))
plot(sd_GW,mu_GW,type="l",col=col_new[6],lwd=3,xlim=c(min(sd_GW,sd_GCW,sd_all,sd_GW_wo,sd_GCW_wo,sd_all_wo),max(sd_GW,sd_GCW,sd_all,sd_GW_wo,sd_GCW_wo,sd_all_wo)),ylim=c(min(mu_GW,mu_GCW,mu_all,mu_GW_wo,mu_GCW_wo,mu_all_wo),max(mu_GW,mu_GCW,mu_all,mu_GW_wo,mu_GCW_wo,mu_all_wo)),xlab=TeX('$\\sigma$'),ylab=TeX('$\\mu$'),main="Efficient Frontier")
lines(sd_GW_wo,mu_GW_wo,col=col_new[6],lwd=3,type="l",lty=2)
lines(sd_GCW,mu_GCW,col=col_b50[40],lwd=3,type="l")
lines(sd_GCW_wo,mu_GCW_wo,col=col_b50[40],lwd=3,type="l",lty=2)
lines(sd_all,mu_all,col=col_b50[10],lwd=3,type="l")
lines(sd_all_wo,mu_all_wo,col=col_b50[10],lwd=3,type="l",lty=2)
legend("topright",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)

plot(beta,mu_GW,type="l",col=col_new[6],lwd=3,ylim=c(min(mu_GW,mu_GCW,mu_all,mu_GW_wo,mu_GCW_wo,mu_all_wo),max(mu_GW,mu_GCW,mu_all,mu_GW_wo,mu_GCW_wo,mu_all_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\mu$'),main="Expected Costs")
lines(beta,mu_GW_wo,type="l",col=col_new[6],lwd=3,lty=2)
lines(beta,mu_GCW,type="l",lwd=3,col=col_b50[40])
lines(beta,mu_GCW_wo,type="l",lwd=3,col=col_b50[40],lty=2)
lines(beta,mu_all,type="l",lwd=3,col=col_b50[10])
lines(beta,mu_all_wo,type="l",lwd=3,col=col_b50[10],lty=2)
legend("bottomright",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)

plot(beta,sd_GW,type="l",col=col_new[6],lwd=3,ylim=c(min(sd_GW,sd_GCW,sd_all,sd_GW_wo,sd_GCW_wo,sd_all_wo),max(sd_GW,sd_GCW,sd_all,sd_GW_wo,sd_GCW_wo,sd_all_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\sigma$'),main="Standard Deviation")
lines(beta,sd_GW_wo,type="l",col=col_new[6],lwd=3,lty=2)
lines(beta,sd_GCW,type="l",lwd=3,col=col_b50[40])
lines(beta,sd_GCW_wo,type="l",lwd=3,col=col_b50[40],lty=2)
lines(beta,sd_all,type="l",lwd=3,col=col_b50[10])
lines(beta,sd_all_wo,type="l",lwd=3,col=col_b50[10],lty=2)
legend("topright",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)

# PDF

pdf("Overall_EF.pdf",width=14,height=8)
par(mar=c(5,5,5,5))
par(mfrow=c(1,1))
#layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))
plot(sd_GW,mu_GW,type="l",col=col_new[6],lwd=3,xlim=c(min(sd_GW,sd_GCW,sd_all,sd_GW_wo,sd_GCW_wo,sd_all_wo),max(sd_GW,sd_GCW,sd_all,sd_GW_wo,sd_GCW_wo,sd_all_wo)),ylim=c(min(mu_GW,mu_GCW,mu_all,mu_GW_wo,mu_GCW_wo,mu_all_wo),max(mu_GW,mu_GCW,mu_all,mu_GW_wo,mu_GCW_wo,mu_all_wo)),xlab=TeX('$\\sigma$'),ylab=TeX('$\\mu$'),main="Efficient Frontier",cex.main=2.5,cex.axis=2.5,cex.lab=2.5)
grid(10,10,lwd=2)
lines(sd_GW,mu_GW,type="l",col=col_new[6],lwd=3)
lines(sd_GW_wo,mu_GW_wo,col=col_new[6],lwd=3,type="l",lty=2)
lines(sd_GCW,mu_GCW,col=col_b50[40],lwd=3,type="l")
lines(sd_GCW_wo,mu_GCW_wo,col=col_b50[40],lwd=3,type="l",lty=2)
lines(sd_all,mu_all,col=col_b50[10],lwd=3,type="l")
lines(sd_all_wo,mu_all_wo,col=col_b50[10],lwd=3,type="l",lty=2)
legend("topright",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)
dev.off()
# 
# plot(beta,mu_GW,type="l",col=col_new[6],lwd=3,ylim=c(min(mu_GW,mu_GCW,mu_all,mu_GW_wo,mu_GCW_wo,mu_all_wo),max(mu_GW,mu_GCW,mu_all,mu_GW_wo,mu_GCW_wo,mu_all_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\mu$'),main="Expected Costs")
# lines(beta,mu_GW_wo,type="l",col=col_new[6],lwd=3,lty=2)
# lines(beta,mu_GCW,type="l",lwd=3,col=col_b50[40])
# lines(beta,mu_GCW_wo,type="l",lwd=3,col=col_b50[40],lty=2)
# lines(beta,mu_all,type="l",lwd=3,col=col_b50[10])
# lines(beta,mu_all_wo,type="l",lwd=3,col=col_b50[10],lty=2)
# legend("bottomright",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)
# 
# plot(beta,sd_GW,type="l",col=col_new[6],lwd=3,ylim=c(min(sd_GW,sd_GCW,sd_all,sd_GW_wo,sd_GCW_wo,sd_all_wo),max(sd_GW,sd_GCW,sd_all,sd_GW_wo,sd_GCW_wo,sd_all_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\sigma$'),main="Standard Deviation")
# lines(beta,sd_GW_wo,type="l",col=col_new[6],lwd=3,lty=2)
# lines(beta,sd_GCW,type="l",lwd=3,col=col_b50[40])
# lines(beta,sd_GCW_wo,type="l",lwd=3,col=col_b50[40],lty=2)
# lines(beta,sd_all,type="l",lwd=3,col=col_b50[10])
# lines(beta,sd_all_wo,type="l",lwd=3,col=col_b50[10],lty=2)
# legend("topright",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)
# dev.off()

# ===== Sensitivity Analyses =====

# === 1. Average CO2-Price of 50 EUR / mt ===

source("tech_costs_par.R")

source("thermal_preparation.R")

# == All Technologies ==

J<-3

mean_t<-mean_co2_t
var_t<-var_co2_t

# === Initialize Lists, Vectors and Matrices for the Storing of Results ===

# With System Costs

res_co2_all<-vector(mode="list",length=length(beta))
x_co2_all<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_co2_all<-rep(0,length(beta))
mu_co2_all<-rep(0,length(beta))
sd_co2_all<-rep(0,length(beta))
it_co2_all<-rep(0,length(beta))
status_co2_all<-rep(0,length(beta))

# Without System Costs

res_co2_all_wo<-vector(mode="list",length=length(beta))
x_co2_all_wo<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_co2_all_wo<-rep(0,length(beta))
mu_co2_all_wo<-rep(0,length(beta))
sd_co2_all_wo<-rep(0,length(beta))
it_co2_all_wo<-rep(0,length(beta))
status_co2_all_wo<-rep(0,length(beta))


# === Set Starting Value ===

x0<-rep(0.5,N*(J+1))
x0_wo<-rep(0.5,N*(J+1))

# === Compute Solution ===

for (i in 1:length(beta)) {
  f<-function(y) fn(y,beta[i])
  g<-function(y) gr(y,beta[i])
  res_co2_all[[i]]<-nloptr(x0=x0,eval_f=f,eval_grad_f=g,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_co2_all[i,]<-round(res_co2_all[[i]]$solution,digits=4)
  OF_co2_all[i]<-res_co2_all[[i]]$objective
  it_co2_all[i]<-res_co2_all[[i]]$iterations
  status_co2_all[i]<-res_co2_all[[i]]$status
#  x0<-x_co2_all[i,]
  ########################################################################################################################################################################################################## 
  f_wo<-function(y) fn_wo(y,beta[i])
  g_wo<-function(y) gr_wo(y,beta[i])
  res_co2_all_wo[[i]]<-nloptr(x0=x0_wo,eval_f=f_wo,eval_grad_f=g_wo,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_co2_all_wo[i,]<-round(res_co2_all_wo[[i]]$solution,digits=4)
  OF_co2_all_wo[i]<-res_co2_all_wo[[i]]$objective
  it_co2_all_wo[i]<-res_co2_all_wo[[i]]$iterations
  status_co2_all_wo[i]<-res_co2_all_wo[[i]]$status
#  x0_wo<-x_co2_all_wo[i,]
}

c(max(status_co2_all),min(status_co2_all),max(status_co2_all_wo),min(status_co2_all_wo))


# === Compute the Relative Shares of Each Technology in Each Load Block ===

x_rel_co2_all<-round(sweep(x_co2_all,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)
x_rel_co2_all_wo<-round(sweep(x_co2_all_wo,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)


# === Evaluate MU and SD for Optimal Values ===

for (i in 1:length(beta)) {
  mu_co2_all[i]<-mu(x_co2_all[i,])
  sd_co2_all[i]<-sd(x_co2_all[i,])
  mu_co2_all_wo[i]<-mu_wo(x_co2_all_wo[i,])
  sd_co2_all_wo[i]<-sd_wo(x_co2_all_wo[i,])
}

# === Cumulative Weights per Technology (Aggregated Over All Blocks) ===

# With and Without System Costs

x_cum_co2_all<-matrix(0,nrow=length(beta),ncol=J+1)
x_cum_co2_all_wo<-matrix(0,nrow=length(beta),ncol=J+1)
for (i in 1:length(beta)) {
  x_cum_co2_all[i,]<-c(sapply(1:J,function(u) sum(x_co2_all[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_co2_all[i,(N*J+1):(N*(J+1))]))
  x_cum_co2_all_wo[i,]<-c(sapply(1:J,function(u) sum(x_co2_all_wo[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_co2_all_wo[i,(N*J+1):(N*(J+1))]))
}


# == 3 Technologies: Gas, Coal and Wind ==

J<-2

IC_t<-IC_t[1:2]
FOM_t<-FOM_t[1:2]
F_t<-F_t[1:2]
VOM_t<-VOM_t[1:2]
mean_t<-mean_co2_t[1:2]
var_t<-var_co2_t[1:2,1:2]

# === Initialize Lists, Vectors and Matrices for the Storing of Results ===

# With System Costs

res_co2_GCW<-vector(mode="list",length=length(beta))
x_co2_GCW<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_co2_GCW<-rep(0,length(beta))
mu_co2_GCW<-rep(0,length(beta))
sd_co2_GCW<-rep(0,length(beta))
it_co2_GCW<-rep(0,length(beta))
status_co2_GCW<-rep(0,length(beta))

# Without System Costs

res_co2_GCW_wo<-vector(mode="list",length=length(beta))
x_co2_GCW_wo<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_co2_GCW_wo<-rep(0,length(beta))
mu_co2_GCW_wo<-rep(0,length(beta))
sd_co2_GCW_wo<-rep(0,length(beta))
it_co2_GCW_wo<-rep(0,length(beta))
status_co2_GCW_wo<-rep(0,length(beta))


# === Set Starting Value ===

x0<-rep(0.5,N*(J+1))
x0_wo<-rep(0.5,N*(J+1))

# === Compute Solution ===

for (i in 1:length(beta)) {
  f<-function(y) fn(y,beta[i])
  g<-function(y) gr(y,beta[i])
  res_co2_GCW[[i]]<-nloptr(x0=x0,eval_f=f,eval_grad_f=g,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_co2_GCW[i,]<-round(res_co2_GCW[[i]]$solution,digits=4)
  OF_co2_GCW[i]<-res_co2_GCW[[i]]$objective
  it_co2_GCW[i]<-res_co2_GCW[[i]]$iterations
  status_co2_GCW[i]<-res_co2_GCW[[i]]$status
#  x0<-x_co2_GCW[i,]
  ########################################################################################################################################################################################################## 
  f_wo<-function(y) fn_wo(y,beta[i])
  g_wo<-function(y) gr_wo(y,beta[i])
  res_co2_GCW_wo[[i]]<-nloptr(x0=x0_wo,eval_f=f_wo,eval_grad_f=g_wo,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_co2_GCW_wo[i,]<-round(res_co2_GCW_wo[[i]]$solution,digits=4)
  OF_co2_GCW_wo[i]<-res_co2_GCW_wo[[i]]$objective
  it_co2_GCW_wo[i]<-res_co2_GCW_wo[[i]]$iterations
  status_co2_GCW_wo[i]<-res_co2_GCW_wo[[i]]$status
#  x0_wo<-x_co2_GCW_wo[i,]
}

c(max(status_co2_GCW),min(status_co2_GCW),max(status_co2_GCW_wo),min(status_co2_GCW_wo))


# === Compute the Relative Shares of Each Technology in Each Load Block ===

x_rel_co2_GCW<-round(sweep(x_co2_GCW,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)
x_rel_co2_GCW_wo<-round(sweep(x_co2_GCW_wo,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)


# === Evaluate MU and SD for Optimal Values ===

for (i in 1:length(beta)) {
  mu_co2_GCW[i]<-mu(x_co2_GCW[i,])
  sd_co2_GCW[i]<-sd(x_co2_GCW[i,])
  mu_co2_GCW_wo[i]<-mu_wo(x_co2_GCW_wo[i,])
  sd_co2_GCW_wo[i]<-sd_wo(x_co2_GCW_wo[i,])
}

# === Cumulative Weights per Technology (Aggregated Over All Blocks) ===

# With and Without System Costs

x_cum_co2_GCW<-matrix(0,nrow=length(beta),ncol=J+1)
x_cum_co2_GCW_wo<-matrix(0,nrow=length(beta),ncol=J+1)
for (i in 1:length(beta)) {
  x_cum_co2_GCW[i,]<-c(sapply(1:J,function(u) sum(x_co2_GCW[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_co2_GCW[i,(N*J+1):(N*(J+1))]))
  x_cum_co2_GCW_wo[i,]<-c(sapply(1:J,function(u) sum(x_co2_GCW_wo[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_co2_GCW_wo[i,(N*J+1):(N*(J+1))]))
}


# == 2 Technologies: Gas and Wind ==

J<-1

IC_t<-IC_t[1]
FOM_t<-FOM_t[1]
F_t<-F_t[1]
VOM_t<-VOM_t[1]
mean_t<-mean_co2_t[1]
var_t<-var_co2_t[1,1]

# === Initialize Lists, Vectors and Matrices for the Storing of Results ===

# With System Costs

res_co2_GW<-vector(mode="list",length=length(beta))
x_co2_GW<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_co2_GW<-rep(0,length(beta))
mu_co2_GW<-rep(0,length(beta))
sd_co2_GW<-rep(0,length(beta))
it_co2_GW<-rep(0,length(beta))
status_co2_GW<-rep(0,length(beta))

# Without System Costs

res_co2_GW_wo<-vector(mode="list",length=length(beta))
x_co2_GW_wo<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_co2_GW_wo<-rep(0,length(beta))
mu_co2_GW_wo<-rep(0,length(beta))
sd_co2_GW_wo<-rep(0,length(beta))
it_co2_GW_wo<-rep(0,length(beta))
status_co2_GW_wo<-rep(0,length(beta))


# === Set Starting Value ===

x0<-rep(0.5,N*(J+1))
x0_wo<-rep(0.5,N*(J+1))

# === Compute Solution ===

for (i in 1:length(beta)) {
  f<-function(y) fn(y,beta[i])
  g<-function(y) gr(y,beta[i])
  res_co2_GW[[i]]<-nloptr(x0=x0,eval_f=f,eval_grad_f=g,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_co2_GW[i,]<-round(res_co2_GW[[i]]$solution,digits=4)
  OF_co2_GW[i]<-res_co2_GW[[i]]$objective
  it_co2_GW[i]<-res_co2_GW[[i]]$iterations
  status_co2_GW[i]<-res_co2_GW[[i]]$status
#  x0<-x_co2_GW[i,]
  ########################################################################################################################################################################################################## 
  f_wo<-function(y) fn_wo(y,beta[i])
  g_wo<-function(y) gr_wo(y,beta[i])
  res_co2_GW_wo[[i]]<-nloptr(x0=x0_wo,eval_f=f_wo,eval_grad_f=g_wo,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_co2_GW_wo[i,]<-round(res_co2_GW_wo[[i]]$solution,digits=4)
  OF_co2_GW_wo[i]<-res_co2_GW_wo[[i]]$objective
  it_co2_GW_wo[i]<-res_co2_GW_wo[[i]]$iterations
  status_co2_GW_wo[i]<-res_co2_GW_wo[[i]]$status
#  x0_wo<-x_co2_GW_wo[i,]
}

c(max(status_co2_GW),min(status_co2_GW),max(status_co2_GW_wo),min(status_co2_GW_wo))


# === Compute the Relative Shares of Each Technology in Each Load Block ===

x_rel_co2_GW<-round(sweep(x_co2_GW,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)
x_rel_co2_GW_wo<-round(sweep(x_co2_GW_wo,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)


# === Evaluate MU and SD for Optimal Values ===

for (i in 1:length(beta)) {
  mu_co2_GW[i]<-mu(x_co2_GW[i,])
  sd_co2_GW[i]<-sd(x_co2_GW[i,])
  mu_co2_GW_wo[i]<-mu_wo(x_co2_GW_wo[i,])
  sd_co2_GW_wo[i]<-sd_wo(x_co2_GW_wo[i,])
}

# === Cumulative Weights per Technology (Aggregated Over All Blocks) ===

# With and Without System Costs

x_cum_co2_GW<-matrix(0,nrow=length(beta),ncol=J+1)
x_cum_co2_GW_wo<-matrix(0,nrow=length(beta),ncol=J+1)
for (i in 1:length(beta)) {
  x_cum_co2_GW[i,]<-c(sapply(1:J,function(u) sum(x_co2_GW[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_co2_GW[i,(N*J+1):(N*(J+1))]))
  x_cum_co2_GW_wo[i,]<-c(sapply(1:J,function(u) sum(x_co2_GW_wo[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_co2_GW_wo[i,(N*J+1):(N*(J+1))]))
}

# == Overall Plots for CO2-Price=50 EUR/mt ==

# == Efficient Frontier, Expected Costs and Standard Deviation ==

# R

par(mar=c(2,2,2,2))
layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))
plot(sd_co2_GW,mu_co2_GW,type="l",col=col_new[6],lwd=3,lty=1,xlim=c(min(sd_co2_GW,sd_co2_GCW,sd_co2_all,sd_co2_GW_wo,sd_co2_GCW_wo,sd_co2_all_wo),max(sd_co2_GW,sd_co2_GCW,sd_co2_all,sd_co2_GW_wo,sd_co2_GCW_wo,sd_co2_all_wo)),ylim=c(min(mu_co2_GW,mu_co2_GCW,mu_co2_all,mu_co2_GW_wo,mu_co2_GCW_wo,mu_co2_all_wo),max(mu_co2_GW,mu_co2_GCW,mu_co2_all,mu_co2_GW_wo,mu_co2_GCW_wo,mu_co2_all_wo)),xlab=TeX('$\\sigma$'),ylab=TeX('$\\mu$'),main="Efficient Frontier")
lines(sd_co2_GW_wo,mu_co2_GW_wo,col=col_new[6],lwd=3,type="l",lty=2)
lines(sd_co2_GCW,mu_co2_GCW,col=col_b50[40],lwd=3,lty=1,type="l")
lines(sd_co2_GCW_wo,mu_co2_GCW_wo,col=col_b50[40],lwd=3,type="l",lty=2)
lines(sd_co2_all,mu_co2_all,col=col_b50[10],lwd=3,lty=1,type="l")
lines(sd_co2_all_wo,mu_co2_all_wo,col=col_b50[10],lwd=3,type="l",lty=2)
legend("bottomright",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)

plot(beta,mu_co2_GW,type="l",col=col_new[6],lty=1,lwd=3,ylim=c(min(mu_co2_GW,mu_co2_GCW,mu_co2_all,mu_co2_GW_wo,mu_co2_GCW_wo,mu_co2_all_wo),max(mu_co2_GW,mu_co2_GCW,mu_co2_all,mu_co2_GW_wo,mu_co2_GCW_wo,mu_co2_all_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\mu$'),main="Expected Costs")
lines(beta,mu_co2_GW_wo,type="l",col=col_new[6],lwd=3,lty=2)
lines(beta,mu_co2_GCW,type="l",lty=1,lwd=3,col=col_b50[40])
lines(beta,mu_co2_GCW_wo,type="l",lwd=3,col=col_b50[40],lty=2)
lines(beta,mu_co2_all,type="l",lwd=3,lty=1,col=col_b50[10])
lines(beta,mu_co2_all_wo,type="l",lwd=3,col=col_b50[10],lty=2)
legend("bottomright",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)

plot(beta,sd_co2_GW,type="l",col=col_new[6],lty=1,lwd=3,ylim=c(min(sd_co2_GW,sd_co2_GCW,sd_co2_all,sd_co2_GW_wo,sd_co2_GCW_wo,sd_co2_all_wo),max(sd_co2_GW,sd_co2_GCW,sd_co2_all,sd_co2_GW_wo,sd_co2_GCW_wo,sd_co2_all_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\sigma$'),main="Standard Deviation")
lines(beta,sd_co2_GW_wo,type="l",col=col_new[6],lwd=3,lty=2)
lines(beta,sd_co2_GCW,type="l",lwd=3,lty=1,col=col_b50[40])
lines(beta,sd_co2_GCW_wo,type="l",lwd=3,col=col_b50[40],lty=2)
lines(beta,sd_co2_all,type="l",lwd=3,lty=1,col=col_b50[10])
lines(beta,sd_co2_all_wo,type="l",lwd=3,col=col_b50[10],lty=2)
legend("topright",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)

# PDF

# pdf("Overall_EF_MU_SD_CO2.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))
# plot(sd_co2_GW,mu_co2_GW,type="l",col=col_new[6],lwd=3,lty=1,xlim=c(min(sd_co2_GW,sd_co2_GCW,sd_co2_all,sd_co2_GW_wo,sd_co2_GCW_wo,sd_co2_all_wo),max(sd_co2_GW,sd_co2_GCW,sd_co2_all,sd_co2_GW_wo,sd_co2_GCW_wo,sd_co2_all_wo)),ylim=c(min(mu_co2_GW,mu_co2_GCW,mu_co2_all,mu_co2_GW_wo,mu_co2_GCW_wo,mu_co2_all_wo),max(mu_co2_GW,mu_co2_GCW,mu_co2_all,mu_co2_GW_wo,mu_co2_GCW_wo,mu_co2_all_wo)),xlab=TeX('$\\sigma$'),ylab=TeX('$\\mu$'),main="Efficient Frontier")
# lines(sd_co2_GW_wo,mu_co2_GW_wo,col=col_new[6],lwd=3,type="l",lty=2)
# lines(sd_co2_GCW,mu_co2_GCW,col=col_b50[40],lwd=3,lty=1,type="l")
# lines(sd_co2_GCW_wo,mu_co2_GCW_wo,col=col_b50[40],lwd=3,type="l",lty=2)
# lines(sd_co2_all,mu_co2_all,col=col_b50[10],lwd=3,lty=1,type="l")
# lines(sd_co2_all_wo,mu_co2_all_wo,col=col_b50[10],lwd=3,type="l",lty=2)
# legend("bottomright",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)
# 
# plot(beta,mu_co2_GW,type="l",col=col_new[6],lty=1,lwd=3,ylim=c(min(mu_co2_GW,mu_co2_GCW,mu_co2_all,mu_co2_GW_wo,mu_co2_GCW_wo,mu_co2_all_wo),max(mu_co2_GW,mu_co2_GCW,mu_co2_all,mu_co2_GW_wo,mu_co2_GCW_wo,mu_co2_all_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\mu$'),main="Expected Costs")
# lines(beta,mu_co2_GW_wo,type="l",col=col_new[6],lwd=3,lty=2)
# lines(beta,mu_co2_GCW,type="l",lty=1,lwd=3,col=col_b50[40])
# lines(beta,mu_co2_GCW_wo,type="l",lwd=3,col=col_b50[40],lty=2)
# lines(beta,mu_co2_all,type="l",lwd=3,lty=1,col=col_b50[10])
# lines(beta,mu_co2_all_wo,type="l",lwd=3,col=col_b50[10],lty=2)
# legend("bottomright",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)
# 
# plot(beta,sd_co2_GW,type="l",col=col_new[6],lty=1,lwd=3,ylim=c(min(sd_co2_GW,sd_co2_GCW,sd_co2_all,sd_co2_GW_wo,sd_co2_GCW_wo,sd_co2_all_wo),max(sd_co2_GW,sd_co2_GCW,sd_co2_all,sd_co2_GW_wo,sd_co2_GCW_wo,sd_co2_all_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\sigma$'),main="Standard Deviation")
# lines(beta,sd_co2_GW_wo,type="l",col=col_new[6],lwd=3,lty=2)
# lines(beta,sd_co2_GCW,type="l",lwd=3,lty=1,col=col_b50[40])
# lines(beta,sd_co2_GCW_wo,type="l",lwd=3,col=col_b50[40],lty=2)
# lines(beta,sd_co2_all,type="l",lwd=3,lty=1,col=col_b50[10])
# lines(beta,sd_co2_all_wo,type="l",lwd=3,col=col_b50[10],lty=2)
# legend("topright",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)
# dev.off()

# == Shares on the EF == 

# R

par(mar=c(2,2,2,2))
par(mfrow=c(3,2))
plot(beta,x_cum_co2_all[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal, Nuclear and Wind With")
lines(beta,x_cum_co2_all[,1]+x_cum_co2_all[,2],type="l")
lines(beta,x_cum_co2_all[,1]+x_cum_co2_all[,2]+x_cum_co2_all[,3],type="l")
lines(beta,x_cum_co2_all[,1]+x_cum_co2_all[,2]+x_cum_co2_all[,3]+x_cum_co2_all[,4],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_co2_all[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_co2_all[,1],rev(x_cum_co2_all[,1]+x_cum_co2_all[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_co2_all[,1]+x_cum_co2_all[,2],rev(x_cum_co2_all[,1]+x_cum_co2_all[,2]+x_cum_co2_all[,3])),col=adjustcolor(colorRampPalette(colors)(4)[3]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_co2_all[,1]+x_cum_co2_all[,2]+x_cum_co2_all[,3],rev(x_cum_co2_all[,1]+x_cum_co2_all[,2]+x_cum_co2_all[,3]+x_cum_co2_all[,4])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.035,labels="gas",cex=1.2,col="black",pos=4)
text(3.5,0.4,labels="nuclear",cex=1.2,col="white",pos=4)
text(6.5,0.83,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_co2_all_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal, Nuclear and Wind Without")
lines(beta,x_cum_co2_all_wo[,1]+x_cum_co2_all_wo[,2],type="l")
lines(beta,x_cum_co2_all_wo[,1]+x_cum_co2_all_wo[,2]+x_cum_co2_all_wo[,3],type="l")
lines(beta,x_cum_co2_all_wo[,1]+x_cum_co2_all_wo[,2]+x_cum_co2_all_wo[,3]+x_cum_co2_all_wo[,4],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_co2_all_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_co2_all_wo[,1],rev(x_cum_co2_all_wo[,1]+x_cum_co2_all_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_co2_all_wo[,1]+x_cum_co2_all_wo[,2],rev(x_cum_co2_all_wo[,1]+x_cum_co2_all_wo[,2]+x_cum_co2_all_wo[,3])),col=adjustcolor(colorRampPalette(colors)(4)[3]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_co2_all_wo[,1]+x_cum_co2_all_wo[,2]+x_cum_co2_all_wo[,3],rev(x_cum_co2_all_wo[,1]+x_cum_co2_all_wo[,2]+x_cum_co2_all_wo[,3]+x_cum_co2_all_wo[,4])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.035,labels="gas",cex=1.2,col="black",pos=4)
text(2.5,0.15,labels="nuclear",cex=1.2,col="white",pos=4)
text(5.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_co2_GCW[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal and Wind With")
lines(beta,x_cum_co2_GCW[,1]+x_cum_co2_GCW[,2],type="l")
lines(beta,x_cum_co2_GCW[,1]+x_cum_co2_GCW[,2]+x_cum_co2_GCW[,3],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_co2_GCW[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_co2_GCW[,1],rev(x_cum_co2_GCW[,1]+x_cum_co2_GCW[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_co2_GCW[,1]+x_cum_co2_GCW[,2],rev(x_cum_co2_GCW[,1]+x_cum_co2_GCW[,2]+x_cum_co2_GCW[,3])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0.05,0.1,labels="gas",cex=1.2,col="black",pos=4)
text(9.2,0.06,labels="coal",cex=1.2,col="white",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_co2_GCW_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal and Wind Without")
lines(beta,x_cum_co2_GCW_wo[,1]+x_cum_co2_GCW_wo[,2],type="l")
lines(beta,x_cum_co2_GCW_wo[,1]+x_cum_co2_GCW_wo[,2]+x_cum_co2_GCW_wo[,3],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_co2_GCW_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_co2_GCW_wo[,1],rev(x_cum_co2_GCW_wo[,1]+x_cum_co2_GCW_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_co2_GCW_wo[,1]+x_cum_co2_GCW_wo[,2],rev(x_cum_co2_GCW_wo[,1]+x_cum_co2_GCW_wo[,2]+x_cum_co2_GCW_wo[,3])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.05,labels="gas",cex=1.2,col="black",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_co2_GW[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas and Wind With")
lines(beta,x_cum_co2_GW[,1]+x_cum_co2_GW[,2],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_co2_GW[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_co2_GW[,1],rev(x_cum_co2_GW[,1]+x_cum_co2_GW[,2])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.08,labels="gas",cex=1.2,col="black",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_co2_GW_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas and Wind Without")
lines(beta,x_cum_co2_GW_wo[,1]+x_cum_co2_GW_wo[,2],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_co2_GW_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_co2_GW_wo[,1],rev(x_cum_co2_GW_wo[,1]+x_cum_co2_GW_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.05,labels="gas",cex=1.2,col="black",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

# PDF

# pdf("Shares_EF_CO2.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# par(mfrow=c(3,2))
# plot(beta,x_cum_co2_all[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal, Nuclear and Wind With")
# lines(beta,x_cum_co2_all[,1]+x_cum_co2_all[,2],type="l")
# lines(beta,x_cum_co2_all[,1]+x_cum_co2_all[,2]+x_cum_co2_all[,3],type="l")
# lines(beta,x_cum_co2_all[,1]+x_cum_co2_all[,2]+x_cum_co2_all[,3]+x_cum_co2_all[,4],type="l")
# polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_co2_all[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_co2_all[,1],rev(x_cum_co2_all[,1]+x_cum_co2_all[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_co2_all[,1]+x_cum_co2_all[,2],rev(x_cum_co2_all[,1]+x_cum_co2_all[,2]+x_cum_co2_all[,3])),col=adjustcolor(colorRampPalette(colors)(4)[3]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_co2_all[,1]+x_cum_co2_all[,2]+x_cum_co2_all[,3],rev(x_cum_co2_all[,1]+x_cum_co2_all[,2]+x_cum_co2_all[,3]+x_cum_co2_all[,4])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
# abline(h=0)
# abline(h=1)
# text(0,0.035,labels="gas",cex=1.2,col="black",pos=4)
# text(3.5,0.4,labels="nuclear",cex=1.2,col="white",pos=4)
# text(6.5,0.83,labels="wind",cex=1.2,col="white",pos=4)
# 
# plot(beta,x_cum_co2_all_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal, Nuclear and Wind Without")
# lines(beta,x_cum_co2_all_wo[,1]+x_cum_co2_all_wo[,2],type="l")
# lines(beta,x_cum_co2_all_wo[,1]+x_cum_co2_all_wo[,2]+x_cum_co2_all_wo[,3],type="l")
# lines(beta,x_cum_co2_all_wo[,1]+x_cum_co2_all_wo[,2]+x_cum_co2_all_wo[,3]+x_cum_co2_all_wo[,4],type="l")
# polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_co2_all_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_co2_all_wo[,1],rev(x_cum_co2_all_wo[,1]+x_cum_co2_all_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_co2_all_wo[,1]+x_cum_co2_all_wo[,2],rev(x_cum_co2_all_wo[,1]+x_cum_co2_all_wo[,2]+x_cum_co2_all_wo[,3])),col=adjustcolor(colorRampPalette(colors)(4)[3]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_co2_all_wo[,1]+x_cum_co2_all_wo[,2]+x_cum_co2_all_wo[,3],rev(x_cum_co2_all_wo[,1]+x_cum_co2_all_wo[,2]+x_cum_co2_all_wo[,3]+x_cum_co2_all_wo[,4])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
# abline(h=0)
# abline(h=1)
# text(0,0.035,labels="gas",cex=1.2,col="black",pos=4)
# text(2.5,0.15,labels="nuclear",cex=1.2,col="white",pos=4)
# text(5.5,0.6,labels="wind",cex=1.2,col="white",pos=4)
# 
# plot(beta,x_cum_co2_GCW[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal and Wind With")
# lines(beta,x_cum_co2_GCW[,1]+x_cum_co2_GCW[,2],type="l")
# lines(beta,x_cum_co2_GCW[,1]+x_cum_co2_GCW[,2]+x_cum_co2_GCW[,3],type="l")
# polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_co2_GCW[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_co2_GCW[,1],rev(x_cum_co2_GCW[,1]+x_cum_co2_GCW[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_co2_GCW[,1]+x_cum_co2_GCW[,2],rev(x_cum_co2_GCW[,1]+x_cum_co2_GCW[,2]+x_cum_co2_GCW[,3])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
# abline(h=0)
# abline(h=1)
# text(0.05,0.1,labels="gas",cex=1.2,col="black",pos=4)
# text(9.2,0.06,labels="coal",cex=1.2,col="white",pos=4)
# text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)
# 
# plot(beta,x_cum_co2_GCW_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal and Wind Without")
# lines(beta,x_cum_co2_GCW_wo[,1]+x_cum_co2_GCW_wo[,2],type="l")
# lines(beta,x_cum_co2_GCW_wo[,1]+x_cum_co2_GCW_wo[,2]+x_cum_co2_GCW_wo[,3],type="l")
# polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_co2_GCW_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_co2_GCW_wo[,1],rev(x_cum_co2_GCW_wo[,1]+x_cum_co2_GCW_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_co2_GCW_wo[,1]+x_cum_co2_GCW_wo[,2],rev(x_cum_co2_GCW_wo[,1]+x_cum_co2_GCW_wo[,2]+x_cum_co2_GCW_wo[,3])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
# abline(h=0)
# abline(h=1)
# text(0,0.05,labels="gas",cex=1.2,col="black",pos=4)
# text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)
# 
# plot(beta,x_cum_co2_GW[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas and Wind With")
# lines(beta,x_cum_co2_GW[,1]+x_cum_co2_GW[,2],type="l")
# polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_co2_GW[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_co2_GW[,1],rev(x_cum_co2_GW[,1]+x_cum_co2_GW[,2])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
# abline(h=0)
# abline(h=1)
# text(0,0.08,labels="gas",cex=1.2,col="black",pos=4)
# text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)
# 
# plot(beta,x_cum_co2_GW_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas and Wind Without")
# lines(beta,x_cum_co2_GW_wo[,1]+x_cum_co2_GW_wo[,2],type="l")
# polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_co2_GW_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
# polygon(c(beta,rev(beta)),c(x_cum_co2_GW_wo[,1],rev(x_cum_co2_GW_wo[,1]+x_cum_co2_GW_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
# abline(h=0)
# abline(h=1)
# text(0,0.05,labels="gas",cex=1.2,col="black",pos=4)
# text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)
# dev.off()

# == Shares in Each Load Block ==

# = beta=0 =

# R

par(mar=c(2,2,2,2))
par(mfrow=c(3,2))
plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Nuclear"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Nuclear"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_co2_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_co2_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_co2_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_co2_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_co2_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_co2_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_co2_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_co2_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)


# PDF

# CairoPDF("Shares_LB_CO2_0.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# par(mfrow=c(3,2))
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Nuclear"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3])),lwd=4)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Nuclear"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3])),lwd=4)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_co2_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_co2_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_co2_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_co2_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_co2_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_co2_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_co2_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_co2_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
# 
# dev.off()

# = beta=10 =

# R

par(mar=c(2,2,2,2))
par(mfrow=c(3,2))
plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Nuclear","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_co2_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_co2_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Nuclear","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-3]+x_rel_co2_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_co2_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3]+x_rel_co2_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_co2_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_co2_GCW[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_co2_GCW[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_co2_GCW[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_co2_GCW[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Coal","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_co2_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_co2_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_co2_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_co2_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-3]+x_rel_co2_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_co2_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3]+x_rel_co2_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_co2_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_co2_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_co2_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_co2_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_co2_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)


# PDF

# CairoPDF("Shares_LB_CO2_10.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# par(mfrow=c(3,2))
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Nuclear","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_co2_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_co2_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[3]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Nuclear","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-3]+x_rel_co2_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_co2_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3]+x_rel_co2_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_co2_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_co2_GCW[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_co2_GCW[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[2]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_co2_GCW[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_co2_GCW[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Coal","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_co2_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_co2_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_co2_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_co2_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-3]+x_rel_co2_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_co2_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3]+x_rel_co2_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_co2_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
# polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_co2_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_co2_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
# polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_co2_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_co2_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
# legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
# dev.off()

# === 2. Decreased Investment Costs in Wind ===

source("tech_costs_par.R")

source("thermal_preparation.R")

F_r<-F_r_ic


# == All Technologies ==

J<-3

# === Initialize Lists, Vectors and Matrices for the Storing of Results ===

# With System Costs

res_ic_all<-vector(mode="list",length=length(beta))
x_ic_all<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_ic_all<-rep(0,length(beta))
mu_ic_all<-rep(0,length(beta))
sd_ic_all<-rep(0,length(beta))
it_ic_all<-rep(0,length(beta))
status_ic_all<-rep(0,length(beta))

# Without System Costs

res_ic_all_wo<-vector(mode="list",length=length(beta))
x_ic_all_wo<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_ic_all_wo<-rep(0,length(beta))
mu_ic_all_wo<-rep(0,length(beta))
sd_ic_all_wo<-rep(0,length(beta))
it_ic_all_wo<-rep(0,length(beta))
status_ic_all_wo<-rep(0,length(beta))


# === Set Starting Value ===

x0<-rep(0.5,N*(J+1))
x0_wo<-rep(0.5,N*(J+1))

# === Compute Solution ===

for (i in 1:length(beta)) {
  f<-function(y) fn(y,beta[i])
  g<-function(y) gr(y,beta[i])
  res_ic_all[[i]]<-nloptr(x0=x0,eval_f=f,eval_grad_f=g,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_ic_all[i,]<-round(res_ic_all[[i]]$solution,digits=4)
  OF_ic_all[i]<-res_ic_all[[i]]$objective
  it_ic_all[i]<-res_ic_all[[i]]$iterations
  status_ic_all[i]<-res_ic_all[[i]]$status
#  x0<-x_ic_all[i,]
  ########################################################################################################################################################################################################## 
  f_wo<-function(y) fn_wo(y,beta[i])
  g_wo<-function(y) gr_wo(y,beta[i])
  res_ic_all_wo[[i]]<-nloptr(x0=x0_wo,eval_f=f_wo,eval_grad_f=g_wo,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_ic_all_wo[i,]<-round(res_ic_all_wo[[i]]$solution,digits=4)
  OF_ic_all_wo[i]<-res_ic_all_wo[[i]]$objective
  it_ic_all_wo[i]<-res_ic_all_wo[[i]]$iterations
  status_ic_all_wo[i]<-res_ic_all_wo[[i]]$status
#  x0_wo<-x_ic_all_wo[i,]
}

c(max(status_ic_all),min(status_ic_all),max(status_ic_all_wo),min(status_ic_all_wo))


# === Compute the Relative Shares of Each Technology in Each Load Block ===

x_rel_ic_all<-round(sweep(x_ic_all,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)
x_rel_ic_all_wo<-round(sweep(x_ic_all_wo,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)


# === Evaluate MU and SD for Optimal Values ===

for (i in 1:length(beta)) {
  mu_ic_all[i]<-mu(x_ic_all[i,])
  sd_ic_all[i]<-sd(x_ic_all[i,])
  mu_ic_all_wo[i]<-mu_wo(x_ic_all_wo[i,])
  sd_ic_all_wo[i]<-sd_wo(x_ic_all_wo[i,])
}

# === Cumulative Weights per Technology (Aggregated Over All Blocks) ===

# With and Without System Costs

x_cum_ic_all<-matrix(0,nrow=length(beta),ncol=J+1)
x_cum_ic_all_wo<-matrix(0,nrow=length(beta),ncol=J+1)
for (i in 1:length(beta)) {
  x_cum_ic_all[i,]<-c(sapply(1:J,function(u) sum(x_ic_all[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_ic_all[i,(N*J+1):(N*(J+1))]))
  x_cum_ic_all_wo[i,]<-c(sapply(1:J,function(u) sum(x_ic_all_wo[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_ic_all_wo[i,(N*J+1):(N*(J+1))]))
}

# == 3 Technologies: Gas, Coal and Wind ==

J<-2

IC_t<-IC_t[1:2]
FOM_t<-FOM_t[1:2]
F_t<-F_t[1:2]
VOM_t<-VOM_t[1:2]
mean_t<-mean_t[1:2]
var_t<-var_t[1:2,1:2]

# === Initialize Lists, Vectors and Matrices for the Storing of Results ===

# With System Costs

res_ic_GCW<-vector(mode="list",length=length(beta))
x_ic_GCW<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_ic_GCW<-rep(0,length(beta))
mu_ic_GCW<-rep(0,length(beta))
sd_ic_GCW<-rep(0,length(beta))
it_ic_GCW<-rep(0,length(beta))
status_ic_GCW<-rep(0,length(beta))

# Without System Costs

res_ic_GCW_wo<-vector(mode="list",length=length(beta))
x_ic_GCW_wo<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_ic_GCW_wo<-rep(0,length(beta))
mu_ic_GCW_wo<-rep(0,length(beta))
sd_ic_GCW_wo<-rep(0,length(beta))
it_ic_GCW_wo<-rep(0,length(beta))
status_ic_GCW_wo<-rep(0,length(beta))


# === Set Starting Value ===

x0<-rep(0.5,N*(J+1))
x0_wo<-rep(0.5,N*(J+1))

# === Compute Solution ===

for (i in 1:length(beta)) {
  f<-function(y) fn(y,beta[i])
  g<-function(y) gr(y,beta[i])
  res_ic_GCW[[i]]<-nloptr(x0=x0,eval_f=f,eval_grad_f=g,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_ic_GCW[i,]<-round(res_ic_GCW[[i]]$solution,digits=4)
  OF_ic_GCW[i]<-res_ic_GCW[[i]]$objective
  it_ic_GCW[i]<-res_ic_GCW[[i]]$iterations
  status_ic_GCW[i]<-res_ic_GCW[[i]]$status
#  x0<-x_ic_GCW[i,]
  ########################################################################################################################################################################################################## 
  f_wo<-function(y) fn_wo(y,beta[i])
  g_wo<-function(y) gr_wo(y,beta[i])
  res_ic_GCW_wo[[i]]<-nloptr(x0=x0_wo,eval_f=f_wo,eval_grad_f=g_wo,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_ic_GCW_wo[i,]<-round(res_ic_GCW_wo[[i]]$solution,digits=4)
  OF_ic_GCW_wo[i]<-res_ic_GCW_wo[[i]]$objective
  it_ic_GCW_wo[i]<-res_ic_GCW_wo[[i]]$iterations
  status_ic_GCW_wo[i]<-res_ic_GCW_wo[[i]]$status
#  x0_wo<-x_ic_GCW_wo[i,]
}

c(max(status_ic_GCW),min(status_ic_GCW),max(status_ic_GCW_wo),min(status_ic_GCW_wo))


# === Compute the Relative Shares of Each Technology in Each Load Block ===

x_rel_ic_GCW<-round(sweep(x_ic_GCW,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)
x_rel_ic_GCW_wo<-round(sweep(x_ic_GCW_wo,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)


# === Evaluate MU and SD for Optimal Values ===

for (i in 1:length(beta)) {
  mu_ic_GCW[i]<-mu(x_ic_GCW[i,])
  sd_ic_GCW[i]<-sd(x_ic_GCW[i,])
  mu_ic_GCW_wo[i]<-mu_wo(x_ic_GCW_wo[i,])
  sd_ic_GCW_wo[i]<-sd_wo(x_ic_GCW_wo[i,])
}

# === Cumulative Weights per Technology (Aggregated Over All Blocks) ===

# With and Without System Costs

x_cum_ic_GCW<-matrix(0,nrow=length(beta),ncol=J+1)
x_cum_ic_GCW_wo<-matrix(0,nrow=length(beta),ncol=J+1)
for (i in 1:length(beta)) {
  x_cum_ic_GCW[i,]<-c(sapply(1:J,function(u) sum(x_ic_GCW[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_ic_GCW[i,(N*J+1):(N*(J+1))]))
  x_cum_ic_GCW_wo[i,]<-c(sapply(1:J,function(u) sum(x_ic_GCW_wo[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_ic_GCW_wo[i,(N*J+1):(N*(J+1))]))
}


# == 2 Technologies: Gas and Wind ==

J<-1

IC_t<-IC_t[1]
FOM_t<-FOM_t[1]
F_t<-F_t[1]
VOM_t<-VOM_t[1]
mean_t<-mean_t[1]
var_t<-var_t[1,1]

# === Initialize Lists, Vectors and Matrices for the Storing of Results ===

# With System Costs

res_ic_GW<-vector(mode="list",length=length(beta))
x_ic_GW<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_ic_GW<-rep(0,length(beta))
mu_ic_GW<-rep(0,length(beta))
sd_ic_GW<-rep(0,length(beta))
it_ic_GW<-rep(0,length(beta))
status_ic_GW<-rep(0,length(beta))

# Without System Costs

res_ic_GW_wo<-vector(mode="list",length=length(beta))
x_ic_GW_wo<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_ic_GW_wo<-rep(0,length(beta))
mu_ic_GW_wo<-rep(0,length(beta))
sd_ic_GW_wo<-rep(0,length(beta))
it_ic_GW_wo<-rep(0,length(beta))
status_ic_GW_wo<-rep(0,length(beta))


# === Set Starting Value ===

x0<-rep(0.5,N*(J+1))
x0_wo<-rep(0.5,N*(J+1))

# === Compute Solution ===

for (i in 1:length(beta)) {
  f<-function(y) fn(y,beta[i])
  g<-function(y) gr(y,beta[i])
  res_ic_GW[[i]]<-nloptr(x0=x0,eval_f=f,eval_grad_f=g,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_ic_GW[i,]<-round(res_ic_GW[[i]]$solution,digits=4)
  OF_ic_GW[i]<-res_ic_GW[[i]]$objective
  it_ic_GW[i]<-res_ic_GW[[i]]$iterations
  status_ic_GW[i]<-res_ic_GW[[i]]$status
#  x0<-x_ic_GW[i,]
  ########################################################################################################################################################################################################## 
  f_wo<-function(y) fn_wo(y,beta[i])
  g_wo<-function(y) gr_wo(y,beta[i])
  res_ic_GW_wo[[i]]<-nloptr(x0=x0_wo,eval_f=f_wo,eval_grad_f=g_wo,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_ic_GW_wo[i,]<-round(res_ic_GW_wo[[i]]$solution,digits=4)
  OF_ic_GW_wo[i]<-res_ic_GW_wo[[i]]$objective
  it_ic_GW_wo[i]<-res_ic_GW_wo[[i]]$iterations
  status_ic_GW_wo[i]<-res_ic_GW_wo[[i]]$status
#  x0_wo<-x_ic_GW_wo[i,]
}

c(max(status_ic_GW),min(status_ic_GW),max(status_ic_GW_wo),min(status_ic_GW_wo))


# === Compute the Relative Shares of Each Technology in Each Load Block ===

x_rel_ic_GW<-round(sweep(x_ic_GW,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)
x_rel_ic_GW_wo<-round(sweep(x_ic_GW_wo,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)


# === Evaluate MU and SD for Optimal Values ===

for (i in 1:length(beta)) {
  mu_ic_GW[i]<-mu(x_ic_GW[i,])
  sd_ic_GW[i]<-sd(x_ic_GW[i,])
  mu_ic_GW_wo[i]<-mu_wo(x_ic_GW_wo[i,])
  sd_ic_GW_wo[i]<-sd_wo(x_ic_GW_wo[i,])
}

# === Cumulative Weights per Technology (Aggregated Over All Blocks) ===

# With and Without System Costs

x_cum_ic_GW<-matrix(0,nrow=length(beta),ncol=J+1)
x_cum_ic_GW_wo<-matrix(0,nrow=length(beta),ncol=J+1)
for (i in 1:length(beta)) {
  x_cum_ic_GW[i,]<-c(sapply(1:J,function(u) sum(x_ic_GW[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_ic_GW[i,(N*J+1):(N*(J+1))]))
  x_cum_ic_GW_wo[i,]<-c(sapply(1:J,function(u) sum(x_ic_GW_wo[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_ic_GW_wo[i,(N*J+1):(N*(J+1))]))
}

# == Overall Plots for Investment Costs of Wind = -20% ==

# == Efficient Frontier, Expected Costs and Standard Deviation ==

# R

par(mar=c(2,2,2,2))
layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))
plot(sd_ic_GW,mu_ic_GW,type="l",col=col_new[6],lwd=3,lty=1,xlim=c(min(sd_ic_GW,sd_ic_GCW,sd_ic_all,sd_ic_GW_wo,sd_ic_GCW_wo,sd_ic_all_wo),max(sd_ic_GW,sd_ic_GCW,sd_ic_all,sd_ic_GW_wo,sd_ic_GCW_wo,sd_ic_all_wo)),ylim=c(min(mu_ic_GW,mu_ic_GCW,mu_ic_all,mu_ic_GW_wo,mu_ic_GCW_wo,mu_ic_all_wo),max(mu_ic_GW,mu_ic_GCW,mu_ic_all,mu_ic_GW_wo,mu_ic_GCW_wo,mu_ic_all_wo)),xlab=TeX('$\\sigma$'),ylab=TeX('$\\mu$'),main="Efficient Frontier")
lines(sd_ic_GW_wo,mu_ic_GW_wo,col=col_new[6],lwd=3,type="l",lty=2)
lines(sd_ic_GCW,mu_ic_GCW,col=col_b50[40],lwd=3,lty=1,type="l")
lines(sd_ic_GCW_wo,mu_ic_GCW_wo,col=col_b50[40],lwd=3,type="l",lty=2)
lines(sd_ic_all,mu_ic_all,col=col_b50[10],lwd=3,lty=1,type="l")
lines(sd_ic_all_wo,mu_ic_all_wo,col=col_b50[10],lwd=3,type="l",lty=2)
legend("topleft",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)

plot(beta,mu_ic_GW,type="l",col=col_new[6],lty=1,lwd=3,ylim=c(min(mu_ic_GW,mu_ic_GCW,mu_ic_all,mu_ic_GW_wo,mu_ic_GCW_wo,mu_ic_all_wo),max(mu_ic_GW,mu_ic_GCW,mu_ic_all,mu_ic_GW_wo,mu_ic_GCW_wo,mu_ic_all_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\mu$'),main="Expected Costs")
lines(beta,mu_ic_GW_wo,type="l",col=col_new[6],lwd=3,lty=2)
lines(beta,mu_ic_GCW,type="l",lty=1,lwd=3,col=col_b50[40])
lines(beta,mu_ic_GCW_wo,type="l",lwd=3,col=col_b50[40],lty=2)
lines(beta,mu_ic_all,type="l",lwd=3,lty=1,col=col_b50[10])
lines(beta,mu_ic_all_wo,type="l",lwd=3,col=col_b50[10],lty=2)
legend("topleft",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)

plot(beta,sd_ic_GW,type="l",col=col_new[6],lty=1,lwd=3,ylim=c(min(sd_ic_GW,sd_ic_GCW,sd_ic_all,sd_ic_GW_wo,sd_ic_GCW_wo,sd_ic_all_wo),max(sd_ic_GW,sd_ic_GCW,sd_ic_all,sd_ic_GW_wo,sd_ic_GCW_wo,sd_ic_all_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\sigma$'),main="Standard Deviation")
lines(beta,sd_ic_GW_wo,type="l",col=col_new[6],lwd=3,lty=2)
lines(beta,sd_ic_GCW,type="l",lwd=3,lty=1,col=col_b50[40])
lines(beta,sd_ic_GCW_wo,type="l",lwd=3,col=col_b50[40],lty=2)
lines(beta,sd_ic_all,type="l",lwd=3,lty=1,col=col_b50[10])
lines(beta,sd_ic_all_wo,type="l",lwd=3,col=col_b50[10],lty=2)
legend("topright",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)


# PDF

pdf("Overall_EF_MU_SD_IC.pdf",width=14,height=8)
par(mar=c(5,5,5,5))
layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))
plot(sd_ic_GW,mu_ic_GW,type="l",col=col_new[6],lwd=3,lty=1,xlim=c(min(sd_ic_GW,sd_ic_GCW,sd_ic_all,sd_ic_GW_wo,sd_ic_GCW_wo,sd_ic_all_wo),max(sd_ic_GW,sd_ic_GCW,sd_ic_all,sd_ic_GW_wo,sd_ic_GCW_wo,sd_ic_all_wo)),ylim=c(min(mu_ic_GW,mu_ic_GCW,mu_ic_all,mu_ic_GW_wo,mu_ic_GCW_wo,mu_ic_all_wo),max(mu_ic_GW,mu_ic_GCW,mu_ic_all,mu_ic_GW_wo,mu_ic_GCW_wo,mu_ic_all_wo)),xlab=TeX('$\\sigma$'),ylab=TeX('$\\mu$'),main="Efficient Frontier")
lines(sd_ic_GW_wo,mu_ic_GW_wo,col=col_new[6],lwd=3,type="l",lty=2)
lines(sd_ic_GCW,mu_ic_GCW,col=col_b50[40],lwd=3,lty=1,type="l")
lines(sd_ic_GCW_wo,mu_ic_GCW_wo,col=col_b50[40],lwd=3,type="l",lty=2)
lines(sd_ic_all,mu_ic_all,col=col_b50[10],lwd=3,lty=1,type="l")
lines(sd_ic_all_wo,mu_ic_all_wo,col=col_b50[10],lwd=3,type="l",lty=2)
legend("topleft",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)

plot(beta,mu_ic_GW,type="l",col=col_new[6],lty=1,lwd=3,ylim=c(min(mu_ic_GW,mu_ic_GCW,mu_ic_all,mu_ic_GW_wo,mu_ic_GCW_wo,mu_ic_all_wo),max(mu_ic_GW,mu_ic_GCW,mu_ic_all,mu_ic_GW_wo,mu_ic_GCW_wo,mu_ic_all_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\mu$'),main="Expected Costs")
lines(beta,mu_ic_GW_wo,type="l",col=col_new[6],lwd=3,lty=2)
lines(beta,mu_ic_GCW,type="l",lty=1,lwd=3,col=col_b50[40])
lines(beta,mu_ic_GCW_wo,type="l",lwd=3,col=col_b50[40],lty=2)
lines(beta,mu_ic_all,type="l",lwd=3,lty=1,col=col_b50[10])
lines(beta,mu_ic_all_wo,type="l",lwd=3,col=col_b50[10],lty=2)
legend("topleft",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)

plot(beta,sd_ic_GW,type="l",col=col_new[6],lty=1,lwd=3,ylim=c(min(sd_ic_GW,sd_ic_GCW,sd_ic_all,sd_ic_GW_wo,sd_ic_GCW_wo,sd_ic_all_wo),max(sd_ic_GW,sd_ic_GCW,sd_ic_all,sd_ic_GW_wo,sd_ic_GCW_wo,sd_ic_all_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\sigma$'),main="Standard Deviation")
lines(beta,sd_ic_GW_wo,type="l",col=col_new[6],lwd=3,lty=2)
lines(beta,sd_ic_GCW,type="l",lwd=3,lty=1,col=col_b50[40])
lines(beta,sd_ic_GCW_wo,type="l",lwd=3,col=col_b50[40],lty=2)
lines(beta,sd_ic_all,type="l",lwd=3,lty=1,col=col_b50[10])
lines(beta,sd_ic_all_wo,type="l",lwd=3,col=col_b50[10],lty=2)
legend("topright",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)
dev.off()

# == Shares on the EF == 

# R

par(mar=c(2,2,2,2))
par(mfrow=c(3,2))
plot(beta,x_cum_ic_all[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal, Nuclear and Wind With")
lines(beta,x_cum_ic_all[,1]+x_cum_ic_all[,2],type="l")
lines(beta,x_cum_ic_all[,1]+x_cum_ic_all[,2]+x_cum_ic_all[,3],type="l")
lines(beta,x_cum_ic_all[,1]+x_cum_ic_all[,2]+x_cum_ic_all[,3]+x_cum_ic_all[,4],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_ic_all[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_all[,1],rev(x_cum_ic_all[,1]+x_cum_ic_all[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_all[,1]+x_cum_ic_all[,2],rev(x_cum_ic_all[,1]+x_cum_ic_all[,2]+x_cum_ic_all[,3])),col=adjustcolor(colorRampPalette(colors)(4)[3]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_all[,1]+x_cum_ic_all[,2]+x_cum_ic_all[,3],rev(x_cum_ic_all[,1]+x_cum_ic_all[,2]+x_cum_ic_all[,3]+x_cum_ic_all[,4])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.05,labels="gas",cex=1.2,col="black",pos=4)
text(1.5,0.03,labels="coal",cex=1.2,col="white",pos=4)
text(3.5,0.4,labels="nuclear",cex=1.2,col="white",pos=4)
text(6.5,0.83,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_ic_all_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal, Nuclear and Wind Without")
lines(beta,x_cum_ic_all_wo[,1]+x_cum_ic_all_wo[,2],type="l")
lines(beta,x_cum_ic_all_wo[,1]+x_cum_ic_all_wo[,2]+x_cum_ic_all_wo[,3],type="l")
lines(beta,x_cum_ic_all_wo[,1]+x_cum_ic_all_wo[,2]+x_cum_ic_all_wo[,3]+x_cum_ic_all_wo[,4],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_ic_all_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_all_wo[,1],rev(x_cum_ic_all_wo[,1]+x_cum_ic_all_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_all_wo[,1]+x_cum_ic_all_wo[,2],rev(x_cum_ic_all_wo[,1]+x_cum_ic_all_wo[,2]+x_cum_ic_all_wo[,3])),col=adjustcolor(colorRampPalette(colors)(4)[3]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_all_wo[,1]+x_cum_ic_all_wo[,2]+x_cum_ic_all_wo[,3],rev(x_cum_ic_all_wo[,1]+x_cum_ic_all_wo[,2]+x_cum_ic_all_wo[,3]+x_cum_ic_all_wo[,4])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.05,labels="gas",cex=1.2,col="black",pos=4)
text(0.2,0.1,labels="coal",cex=1.2,col="white",pos=4)
text(1.8,0.1,labels="nuclear",cex=1.2,col="white",pos=4)
text(5.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_ic_GCW[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal and Wind With")
lines(beta,x_cum_ic_GCW[,1]+x_cum_ic_GCW[,2],type="l")
lines(beta,x_cum_ic_GCW[,1]+x_cum_ic_GCW[,2]+x_cum_ic_GCW[,3],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_ic_GCW[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_GCW[,1],rev(x_cum_ic_GCW[,1]+x_cum_ic_GCW[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_GCW[,1]+x_cum_ic_GCW[,2],rev(x_cum_ic_GCW[,1]+x_cum_ic_GCW[,2]+x_cum_ic_GCW[,3])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.05,labels="gas",cex=1.2,col="black",pos=4)
text(0.8,0.13,labels="coal",cex=1.2,col="white",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_ic_GCW_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal and Wind Without")
lines(beta,x_cum_ic_GCW_wo[,1]+x_cum_ic_GCW_wo[,2],type="l")
lines(beta,x_cum_ic_GCW_wo[,1]+x_cum_ic_GCW_wo[,2]+x_cum_ic_GCW_wo[,3],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_ic_GCW_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_GCW_wo[,1],rev(x_cum_ic_GCW_wo[,1]+x_cum_ic_GCW_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_GCW_wo[,1]+x_cum_ic_GCW_wo[,2],rev(x_cum_ic_GCW_wo[,1]+x_cum_ic_GCW_wo[,2]+x_cum_ic_GCW_wo[,3])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.05,labels="gas",cex=1.2,col="black",pos=4)
text(0.4,0.1,labels="coal",cex=1.2,col="white",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_ic_GW[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas and Wind With")
lines(beta,x_cum_ic_GW[,1]+x_cum_ic_GW[,2],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_ic_GW[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_GW[,1],rev(x_cum_ic_GW[,1]+x_cum_ic_GW[,2])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.08,labels="gas",cex=1.2,col="black",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_ic_GW_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas and Wind Without")
lines(beta,x_cum_ic_GW_wo[,1]+x_cum_ic_GW_wo[,2],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_ic_GW_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_GW_wo[,1],rev(x_cum_ic_GW_wo[,1]+x_cum_ic_GW_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.05,labels="gas",cex=1.2,col="black",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

# PDF

pdf("Shares_EF_IC.pdf",width=14,height=8)
par(mar=c(5,5,5,5))
par(mfrow=c(3,2))
plot(beta,x_cum_ic_all[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal, Nuclear and Wind With")
lines(beta,x_cum_ic_all[,1]+x_cum_ic_all[,2],type="l")
lines(beta,x_cum_ic_all[,1]+x_cum_ic_all[,2]+x_cum_ic_all[,3],type="l")
lines(beta,x_cum_ic_all[,1]+x_cum_ic_all[,2]+x_cum_ic_all[,3]+x_cum_ic_all[,4],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_ic_all[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_all[,1],rev(x_cum_ic_all[,1]+x_cum_ic_all[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_all[,1]+x_cum_ic_all[,2],rev(x_cum_ic_all[,1]+x_cum_ic_all[,2]+x_cum_ic_all[,3])),col=adjustcolor(colorRampPalette(colors)(4)[3]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_all[,1]+x_cum_ic_all[,2]+x_cum_ic_all[,3],rev(x_cum_ic_all[,1]+x_cum_ic_all[,2]+x_cum_ic_all[,3]+x_cum_ic_all[,4])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.05,labels="gas",cex=1.2,col="black",pos=4)
text(1.5,0.03,labels="coal",cex=1.2,col="white",pos=4)
text(3.5,0.4,labels="nuclear",cex=1.2,col="white",pos=4)
text(6.5,0.83,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_ic_all_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal, Nuclear and Wind Without")
lines(beta,x_cum_ic_all_wo[,1]+x_cum_ic_all_wo[,2],type="l")
lines(beta,x_cum_ic_all_wo[,1]+x_cum_ic_all_wo[,2]+x_cum_ic_all_wo[,3],type="l")
lines(beta,x_cum_ic_all_wo[,1]+x_cum_ic_all_wo[,2]+x_cum_ic_all_wo[,3]+x_cum_ic_all_wo[,4],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_ic_all_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_all_wo[,1],rev(x_cum_ic_all_wo[,1]+x_cum_ic_all_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_all_wo[,1]+x_cum_ic_all_wo[,2],rev(x_cum_ic_all_wo[,1]+x_cum_ic_all_wo[,2]+x_cum_ic_all_wo[,3])),col=adjustcolor(colorRampPalette(colors)(4)[3]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_all_wo[,1]+x_cum_ic_all_wo[,2]+x_cum_ic_all_wo[,3],rev(x_cum_ic_all_wo[,1]+x_cum_ic_all_wo[,2]+x_cum_ic_all_wo[,3]+x_cum_ic_all_wo[,4])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.05,labels="gas",cex=1.2,col="black",pos=4)
text(0.2,0.1,labels="coal",cex=1.2,col="white",pos=4)
text(1.8,0.1,labels="nuclear",cex=1.2,col="white",pos=4)
text(5.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_ic_GCW[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal and Wind With")
lines(beta,x_cum_ic_GCW[,1]+x_cum_ic_GCW[,2],type="l")
lines(beta,x_cum_ic_GCW[,1]+x_cum_ic_GCW[,2]+x_cum_ic_GCW[,3],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_ic_GCW[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_GCW[,1],rev(x_cum_ic_GCW[,1]+x_cum_ic_GCW[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_GCW[,1]+x_cum_ic_GCW[,2],rev(x_cum_ic_GCW[,1]+x_cum_ic_GCW[,2]+x_cum_ic_GCW[,3])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.05,labels="gas",cex=1.2,col="black",pos=4)
text(0.8,0.13,labels="coal",cex=1.2,col="white",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_ic_GCW_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal and Wind Without")
lines(beta,x_cum_ic_GCW_wo[,1]+x_cum_ic_GCW_wo[,2],type="l")
lines(beta,x_cum_ic_GCW_wo[,1]+x_cum_ic_GCW_wo[,2]+x_cum_ic_GCW_wo[,3],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_ic_GCW_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_GCW_wo[,1],rev(x_cum_ic_GCW_wo[,1]+x_cum_ic_GCW_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_GCW_wo[,1]+x_cum_ic_GCW_wo[,2],rev(x_cum_ic_GCW_wo[,1]+x_cum_ic_GCW_wo[,2]+x_cum_ic_GCW_wo[,3])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.05,labels="gas",cex=1.2,col="black",pos=4)
text(0.4,0.1,labels="coal",cex=1.2,col="white",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_ic_GW[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas and Wind With")
lines(beta,x_cum_ic_GW[,1]+x_cum_ic_GW[,2],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_ic_GW[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_GW[,1],rev(x_cum_ic_GW[,1]+x_cum_ic_GW[,2])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.08,labels="gas",cex=1.2,col="black",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_ic_GW_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas and Wind Without")
lines(beta,x_cum_ic_GW_wo[,1]+x_cum_ic_GW_wo[,2],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_ic_GW_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_ic_GW_wo[,1],rev(x_cum_ic_GW_wo[,1]+x_cum_ic_GW_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.05,labels="gas",cex=1.2,col="black",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)
dev.off()

# == Shares in Each Load Block ==

# = beta=0 =

# R

par(mar=c(2,2,2,2))
par(mfrow=c(3,2))
plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Nuclear"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Nuclear"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_ic_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_ic_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_ic_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_ic_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Coal","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_ic_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_ic_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_ic_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_ic_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)


# PDF

CairoPDF("Shares_LB_IC_0.pdf",width=14,height=8)
par(mar=c(5,5,5,5))
par(mfrow=c(3,2))
plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Nuclear"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Nuclear"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_ic_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_ic_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_ic_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_ic_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Coal","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_ic_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_ic_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_ic_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_ic_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
dev.off()

# = beta=10 =

# R

par(mar=c(2,2,2,2))
par(mfrow=c(3,2))
plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_ic_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_ic_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_ic_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_ic_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Nuclear","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_ic_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_ic_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_ic_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_ic_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Coal","Nuclear","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[3]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-3]+x_rel_ic_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_ic_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3]+x_rel_ic_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_ic_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_ic_GCW[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_ic_GCW[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Coal","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_ic_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_ic_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_ic_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_ic_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-3]+x_rel_ic_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_ic_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3]+x_rel_ic_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_ic_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_ic_GW_wo[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_ic_GW_wo[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_ic_GW_wo[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_ic_GW_wo[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)


# PDF

CairoPDF("Shares_LB_IC_10.pdf",width=14,height=8)
par(mar=c(5,5,5,5))
par(mfrow=c(3,2))
plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_ic_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_ic_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_ic_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_ic_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Nuclear","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_ic_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_ic_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_ic_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_ic_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Coal","Nuclear","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[3]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-3]+x_rel_ic_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_ic_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3]+x_rel_ic_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_ic_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_ic_GCW[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_ic_GCW[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Coal","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_ic_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_ic_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_ic_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_ic_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-3]+x_rel_ic_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_ic_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3]+x_rel_ic_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_ic_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_ic_GW_wo[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_ic_GW_wo[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_ic_GW_wo[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_ic_GW_wo[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
dev.off()



# === 3. Increased Availability Factor 0.20 ===

source("tech_costs_par.R")

source("thermal_preparation.R")

avfac<-avfac_af


# == All Technologies ==

J<-3

# === Initialize Lists, Vectors and Matrices for the Storing of Results ===

# With System Costs

res_af_all<-vector(mode="list",length=length(beta))
x_af_all<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_af_all<-rep(0,length(beta))
mu_af_all<-rep(0,length(beta))
sd_af_all<-rep(0,length(beta))
it_af_all<-rep(0,length(beta))
status_af_all<-rep(0,length(beta))

# Without System Costs

res_af_all_wo<-vector(mode="list",length=length(beta))
x_af_all_wo<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_af_all_wo<-rep(0,length(beta))
mu_af_all_wo<-rep(0,length(beta))
sd_af_all_wo<-rep(0,length(beta))
it_af_all_wo<-rep(0,length(beta))
status_af_all_wo<-rep(0,length(beta))


# === Set Starting Value ===

x0<-rep(0.5,N*(J+1))
x0_wo<-rep(0.5,N*(J+1))

# === Compute Solution ===

for (i in 1:length(beta)) {
  f<-function(y) fn(y,beta[i])
  g<-function(y) gr(y,beta[i])
  res_af_all[[i]]<-nloptr(x0=x0,eval_f=f,eval_grad_f=g,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_af_all[i,]<-round(res_af_all[[i]]$solution,digits=4)
  OF_af_all[i]<-res_af_all[[i]]$objective
  it_af_all[i]<-res_af_all[[i]]$iterations
  status_af_all[i]<-res_af_all[[i]]$status
#  x0<-x_af_all[i,]
  ########################################################################################################################################################################################################## 
  f_wo<-function(y) fn_wo(y,beta[i])
  g_wo<-function(y) gr_wo(y,beta[i])
  res_af_all_wo[[i]]<-nloptr(x0=x0_wo,eval_f=f_wo,eval_grad_f=g_wo,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_af_all_wo[i,]<-round(res_af_all_wo[[i]]$solution,digits=4)
  OF_af_all_wo[i]<-res_af_all_wo[[i]]$objective
  it_af_all_wo[i]<-res_af_all_wo[[i]]$iterations
  status_af_all_wo[i]<-res_af_all_wo[[i]]$status
#  x0_wo<-x_af_all_wo[i,]
}

c(max(status_af_all),min(status_af_all),max(status_af_all_wo),min(status_af_all_wo))


# === Compute the Relative Shares of Each Technology in Each Load Block ===

x_rel_af_all<-round(sweep(x_af_all,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)
x_rel_af_all_wo<-round(sweep(x_af_all_wo,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)


# === Evaluate MU and SD for Optimal Values ===

for (i in 1:length(beta)) {
  mu_af_all[i]<-mu(x_af_all[i,])
  sd_af_all[i]<-sd(x_af_all[i,])
  mu_af_all_wo[i]<-mu_wo(x_af_all_wo[i,])
  sd_af_all_wo[i]<-sd_wo(x_af_all_wo[i,])
}

# === Cumulative Weights per Technology (Aggregated Over All Blocks) ===

# With and Without System Costs

x_cum_af_all<-matrix(0,nrow=length(beta),ncol=J+1)
x_cum_af_all_wo<-matrix(0,nrow=length(beta),ncol=J+1)
for (i in 1:length(beta)) {
  x_cum_af_all[i,]<-c(sapply(1:J,function(u) sum(x_af_all[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_af_all[i,(N*J+1):(N*(J+1))]))
  x_cum_af_all_wo[i,]<-c(sapply(1:J,function(u) sum(x_af_all_wo[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_af_all_wo[i,(N*J+1):(N*(J+1))]))
}


# == 3 Technologies: Gas, Coal and Wind ==

J<-2

IC_t<-IC_t[1:2]
FOM_t<-FOM_t[1:2]
F_t<-F_t[1:2]
VOM_t<-VOM_t[1:2]
mean_t<-mean_t[1:2]
var_t<-var_t[1:2,1:2]

# === Initialize Lists, Vectors and Matrices for the Storing of Results ===

# With System Costs

res_af_GCW<-vector(mode="list",length=length(beta))
x_af_GCW<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_af_GCW<-rep(0,length(beta))
mu_af_GCW<-rep(0,length(beta))
sd_af_GCW<-rep(0,length(beta))
it_af_GCW<-rep(0,length(beta))
status_af_GCW<-rep(0,length(beta))

# Without System Costs

res_af_GCW_wo<-vector(mode="list",length=length(beta))
x_af_GCW_wo<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_af_GCW_wo<-rep(0,length(beta))
mu_af_GCW_wo<-rep(0,length(beta))
sd_af_GCW_wo<-rep(0,length(beta))
it_af_GCW_wo<-rep(0,length(beta))
status_af_GCW_wo<-rep(0,length(beta))


# === Set Starting Value ===

x0<-rep(0.5,N*(J+1))
x0_wo<-rep(0.5,N*(J+1))

# === Compute Solution ===

for (i in 1:length(beta)) {
  f<-function(y) fn(y,beta[i])
  g<-function(y) gr(y,beta[i])
  res_af_GCW[[i]]<-nloptr(x0=x0,eval_f=f,eval_grad_f=g,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_af_GCW[i,]<-round(res_af_GCW[[i]]$solution,digits=4)
  OF_af_GCW[i]<-res_af_GCW[[i]]$objective
  it_af_GCW[i]<-res_af_GCW[[i]]$iterations
  status_af_GCW[i]<-res_af_GCW[[i]]$status
#  x0<-x_af_GCW[i,]
  ########################################################################################################################################################################################################## 
  f_wo<-function(y) fn_wo(y,beta[i])
  g_wo<-function(y) gr_wo(y,beta[i])
  res_af_GCW_wo[[i]]<-nloptr(x0=x0_wo,eval_f=f_wo,eval_grad_f=g_wo,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_af_GCW_wo[i,]<-round(res_af_GCW_wo[[i]]$solution,digits=4)
  OF_af_GCW_wo[i]<-res_af_GCW_wo[[i]]$objective
  it_af_GCW_wo[i]<-res_af_GCW_wo[[i]]$iterations
  status_af_GCW_wo[i]<-res_af_GCW_wo[[i]]$status
#  x0_wo<-x_af_GCW_wo[i,]
}

c(max(status_af_GCW),min(status_af_GCW),max(status_af_GCW_wo),min(status_af_GCW_wo))

# === Compute the Relative Shares of Each Technology in Each Load Block ===

x_rel_af_GCW<-round(sweep(x_af_GCW,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)
x_rel_af_GCW_wo<-round(sweep(x_af_GCW_wo,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)


# === Evaluate MU and SD for Optimal Values ===

for (i in 1:length(beta)) {
  mu_af_GCW[i]<-mu(x_af_GCW[i,])
  sd_af_GCW[i]<-sd(x_af_GCW[i,])
  mu_af_GCW_wo[i]<-mu_wo(x_af_GCW_wo[i,])
  sd_af_GCW_wo[i]<-sd_wo(x_af_GCW_wo[i,])
}

# === Cumulative Weights per Technology (Aggregated Over All Blocks) ===

# With and Without System Costs

x_cum_af_GCW<-matrix(0,nrow=length(beta),ncol=J+1)
x_cum_af_GCW_wo<-matrix(0,nrow=length(beta),ncol=J+1)
for (i in 1:length(beta)) {
  x_cum_af_GCW[i,]<-c(sapply(1:J,function(u) sum(x_af_GCW[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_af_GCW[i,(N*J+1):(N*(J+1))]))
  x_cum_af_GCW_wo[i,]<-c(sapply(1:J,function(u) sum(x_af_GCW_wo[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_af_GCW_wo[i,(N*J+1):(N*(J+1))]))
}


# == 2 Technologies: Gas and Wind ==

J<-1

IC_t<-IC_t[1]
FOM_t<-FOM_t[1]
F_t<-F_t[1]
VOM_t<-VOM_t[1]
mean_t<-mean_t[1]
var_t<-var_t[1,1]

# === Initialize Lists, Vectors and Matrices for the Storing of Results ===

# With System Costs

res_af_GW<-vector(mode="list",length=length(beta))
x_af_GW<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_af_GW<-rep(0,length(beta))
mu_af_GW<-rep(0,length(beta))
sd_af_GW<-rep(0,length(beta))
it_af_GW<-rep(0,length(beta))
status_af_GW<-rep(0,length(beta))

# Without System Costs

res_af_GW_wo<-vector(mode="list",length=length(beta))
x_af_GW_wo<-matrix(0,nrow=length(beta),ncol=N*(J+1))
OF_af_GW_wo<-rep(0,length(beta))
mu_af_GW_wo<-rep(0,length(beta))
sd_af_GW_wo<-rep(0,length(beta))
it_af_GW_wo<-rep(0,length(beta))
status_af_GW_wo<-rep(0,length(beta))


# === Set Starting Value ===

x0<-rep(0.5,N*(J+1))
x0_wo<-rep(0.5,N*(J+1))

# === Compute Solution ===

for (i in 1:length(beta)) {
  f<-function(y) fn(y,beta[i])
  g<-function(y) gr(y,beta[i])
  res_af_GW[[i]]<-nloptr(x0=x0,eval_f=f,eval_grad_f=g,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_af_GW[i,]<-round(res_af_GW[[i]]$solution,digits=4)
  OF_af_GW[i]<-res_af_GW[[i]]$objective
  it_af_GW[i]<-res_af_GW[[i]]$iterations
  status_af_GW[i]<-res_af_GW[[i]]$status
#  x0<-x_af_GW[i,]
  ########################################################################################################################################################################################################## 
  f_wo<-function(y) fn_wo(y,beta[i])
  g_wo<-function(y) gr_wo(y,beta[i])
  res_af_GW_wo[[i]]<-nloptr(x0=x0_wo,eval_f=f_wo,eval_grad_f=g_wo,eval_g_eq=heq,eval_jac_g_eq=heqjac,lb=rep(0,N*(J+1)),opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1e-4))
  x_af_GW_wo[i,]<-round(res_af_GW_wo[[i]]$solution,digits=4)
  OF_af_GW_wo[i]<-res_af_GW_wo[[i]]$objective
  it_af_GW_wo[i]<-res_af_GW_wo[[i]]$iterations
  status_af_GW_wo[i]<-res_af_GW_wo[[i]]$status
#  x0_wo<-x_af_GW_wo[i,]
}

c(max(status_af_GW),min(status_af_GW),max(status_af_GW_wo),min(status_af_GW_wo))


# === Compute the Relative Shares of Each Technology in Each Load Block ===

x_rel_af_GW<-round(sweep(x_af_GW,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)
x_rel_af_GW_wo<-round(sweep(x_af_GW_wo,MARGIN=2,1/rep(D,J+1)*c(rep(1,N*J),rep(mean(avfac),N)),'*'),3)


# === Evaluate MU and SD for Optimal Values ===

for (i in 1:length(beta)) {
  mu_af_GW[i]<-mu(x_af_GW[i,])
  sd_af_GW[i]<-sd(x_af_GW[i,])
  mu_af_GW_wo[i]<-mu_wo(x_af_GW_wo[i,])
  sd_af_GW_wo[i]<-sd_wo(x_af_GW_wo[i,])
}

# === Cumulative Weights per Technology (Aggregated Over All Blocks) ===

# With and Without System Costs

x_cum_af_GW<-matrix(0,nrow=length(beta),ncol=J+1)
x_cum_af_GW_wo<-matrix(0,nrow=length(beta),ncol=J+1)
for (i in 1:length(beta)) {
  x_cum_af_GW[i,]<-c(sapply(1:J,function(u) sum(x_af_GW[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_af_GW[i,(N*J+1):(N*(J+1))]))
  x_cum_af_GW_wo[i,]<-c(sapply(1:J,function(u) sum(x_af_GW_wo[i,((u-1)*N+1):(u*N)])),mean(avfac)*sum(x_af_GW_wo[i,(N*J+1):(N*(J+1))]))
}

# == Overall Plots ==

# == Efficient Frontier, Expected Costs and Standard Deviation ==

# R

par(mar=c(2,2,2,2))
layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))
plot(sd_af_GW,mu_af_GW,type="l",col=col_new[6],lwd=3,lty=1,xlim=c(min(sd_af_GW,sd_af_GCW,sd_af_all,sd_af_GW_wo,sd_af_GCW_wo,sd_af_all_wo),max(sd_af_GW,sd_af_GCW,sd_af_all,sd_af_GW_wo,sd_af_GCW_wo,sd_af_all_wo)),ylim=c(min(mu_af_GW,mu_af_GCW,mu_af_all,mu_af_GW_wo,mu_af_GCW_wo,mu_af_all_wo),max(mu_af_GW,mu_af_GCW,mu_af_all,mu_af_GW_wo,mu_af_GCW_wo,mu_af_all_wo)),xlab=TeX('$\\sigma$'),ylab=TeX('$\\mu$'),main="Efficient Frontier")
lines(sd_af_GW_wo,mu_af_GW_wo,col=col_new[6],lwd=3,type="l",lty=2)
lines(sd_af_GCW,mu_af_GCW,col=col_b50[40],lwd=3,lty=1,type="l")
lines(sd_af_GCW_wo,mu_af_GCW_wo,col=col_b50[40],lwd=3,type="l",lty=2)
lines(sd_af_all,mu_af_all,col=col_b50[10],lwd=3,lty=1,type="l")
lines(sd_af_all_wo,mu_af_all_wo,col=col_b50[10],lwd=3,type="l",lty=2)
legend("topleft",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)

plot(beta,mu_af_GW,type="l",col=col_new[6],lty=1,lwd=3,ylim=c(min(mu_af_GW,mu_af_GCW,mu_af_all,mu_af_GW_wo,mu_af_GCW_wo,mu_af_all_wo),max(mu_af_GW,mu_af_GCW,mu_af_all,mu_af_GW_wo,mu_af_GCW_wo,mu_af_all_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\mu$'),main="Expected Costs")
lines(beta,mu_af_GW_wo,type="l",col=col_new[6],lwd=3,lty=2)
lines(beta,mu_af_GCW,type="l",lty=1,lwd=3,col=col_b50[40])
lines(beta,mu_af_GCW_wo,type="l",lwd=3,col=col_b50[40],lty=2)
lines(beta,mu_af_all,type="l",lwd=3,lty=1,col=col_b50[10])
lines(beta,mu_af_all_wo,type="l",lwd=3,col=col_b50[10],lty=2)
legend("topleft",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)

plot(beta,sd_af_GW,type="l",col=col_new[6],lty=1,lwd=3,ylim=c(min(sd_af_GW,sd_af_GCW,sd_af_all,sd_af_GW_wo,sd_af_GCW_wo,sd_af_all_wo),max(sd_af_GW,sd_af_GCW,sd_af_all,sd_af_GW_wo,sd_af_GCW_wo,sd_af_all_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\sigma$'),main="Standard Deviation")
lines(beta,sd_af_GW_wo,type="l",col=col_new[6],lwd=3,lty=2)
lines(beta,sd_af_GCW,type="l",lwd=3,lty=1,col=col_b50[40])
lines(beta,sd_af_GCW_wo,type="l",lwd=3,col=col_b50[40],lty=2)
lines(beta,sd_af_all,type="l",lwd=3,lty=1,col=col_b50[10])
lines(beta,sd_af_all_wo,type="l",lwd=3,col=col_b50[10],lty=2)
legend("topright",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)


# PDF

pdf("Overall_EF_MU_SD_AF.pdf",width=14,height=8)
par(mar=c(5,5,5,5))
layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))
plot(sd_af_GW,mu_af_GW,type="l",col=col_new[6],lwd=3,lty=1,xlim=c(min(sd_af_GW,sd_af_GCW,sd_af_all,sd_af_GW_wo,sd_af_GCW_wo,sd_af_all_wo),max(sd_af_GW,sd_af_GCW,sd_af_all,sd_af_GW_wo,sd_af_GCW_wo,sd_af_all_wo)),ylim=c(min(mu_af_GW,mu_af_GCW,mu_af_all,mu_af_GW_wo,mu_af_GCW_wo,mu_af_all_wo),max(mu_af_GW,mu_af_GCW,mu_af_all,mu_af_GW_wo,mu_af_GCW_wo,mu_af_all_wo)),xlab=TeX('$\\sigma$'),ylab=TeX('$\\mu$'),main="Efficient Frontier")
lines(sd_af_GW_wo,mu_af_GW_wo,col=col_new[6],lwd=3,type="l",lty=2)
lines(sd_af_GCW,mu_af_GCW,col=col_b50[40],lwd=3,lty=1,type="l")
lines(sd_af_GCW_wo,mu_af_GCW_wo,col=col_b50[40],lwd=3,type="l",lty=2)
lines(sd_af_all,mu_af_all,col=col_b50[10],lwd=3,lty=1,type="l")
lines(sd_af_all_wo,mu_af_all_wo,col=col_b50[10],lwd=3,type="l",lty=2)
legend("topleft",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)

plot(beta,mu_af_GW,type="l",col=col_new[6],lty=1,lwd=3,ylim=c(min(mu_af_GW,mu_af_GCW,mu_af_all,mu_af_GW_wo,mu_af_GCW_wo,mu_af_all_wo),max(mu_af_GW,mu_af_GCW,mu_af_all,mu_af_GW_wo,mu_af_GCW_wo,mu_af_all_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\mu$'),main="Expected Costs")
lines(beta,mu_af_GW_wo,type="l",col=col_new[6],lwd=3,lty=2)
lines(beta,mu_af_GCW,type="l",lty=1,lwd=3,col=col_b50[40])
lines(beta,mu_af_GCW_wo,type="l",lwd=3,col=col_b50[40],lty=2)
lines(beta,mu_af_all,type="l",lwd=3,lty=1,col=col_b50[10])
lines(beta,mu_af_all_wo,type="l",lwd=3,col=col_b50[10],lty=2)
legend("topleft",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)

plot(beta,sd_af_GW,type="l",col=col_new[6],lty=1,lwd=3,ylim=c(min(sd_af_GW,sd_af_GCW,sd_af_all,sd_af_GW_wo,sd_af_GCW_wo,sd_af_all_wo),max(sd_af_GW,sd_af_GCW,sd_af_all,sd_af_GW_wo,sd_af_GCW_wo,sd_af_all_wo)),xlab=TeX('$\\beta$'),ylab=TeX('$\\sigma$'),main="Standard Deviation")
lines(beta,sd_af_GW_wo,type="l",col=col_new[6],lwd=3,lty=2)
lines(beta,sd_af_GCW,type="l",lwd=3,lty=1,col=col_b50[40])
lines(beta,sd_af_GCW_wo,type="l",lwd=3,col=col_b50[40],lty=2)
lines(beta,sd_af_all,type="l",lwd=3,lty=1,col=col_b50[10])
lines(beta,sd_af_all_wo,type="l",lwd=3,col=col_b50[10],lty=2)
legend("topright",c("GW","GCW","GCNW"),col=c(col_new[6],col_b50[40],col_b50[10]),lwd=3)
dev.off()

# == Shares on the EF == 

# R

par(mar=c(2,2,2,2))
par(mfrow=c(3,2))
plot(beta,x_cum_af_all[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal, Nuclear and Wind With")
lines(beta,x_cum_af_all[,1]+x_cum_af_all[,2],type="l")
lines(beta,x_cum_af_all[,1]+x_cum_af_all[,2]+x_cum_af_all[,3],type="l")
lines(beta,x_cum_af_all[,1]+x_cum_af_all[,2]+x_cum_af_all[,3]+x_cum_af_all[,4],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_af_all[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_all[,1],rev(x_cum_af_all[,1]+x_cum_af_all[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_all[,1]+x_cum_af_all[,2],rev(x_cum_af_all[,1]+x_cum_af_all[,2]+x_cum_af_all[,3])),col=adjustcolor(colorRampPalette(colors)(4)[3]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_all[,1]+x_cum_af_all[,2]+x_cum_af_all[,3],rev(x_cum_af_all[,1]+x_cum_af_all[,2]+x_cum_af_all[,3]+x_cum_af_all[,4])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.035,labels="gas",cex=1.2,col="black",pos=4)
text(3.5,0.4,labels="nuclear",cex=1.2,col="white",pos=4)
text(1.2,0.035,labels="coal",cex=1.2,col="white",pos=4)
text(6.5,0.83,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_af_all_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal, Nuclear and Wind Without")
lines(beta,x_cum_af_all_wo[,1]+x_cum_af_all_wo[,2],type="l")
lines(beta,x_cum_af_all_wo[,1]+x_cum_af_all_wo[,2]+x_cum_af_all_wo[,3],type="l")
lines(beta,x_cum_af_all_wo[,1]+x_cum_af_all_wo[,2]+x_cum_af_all_wo[,3]+x_cum_af_all_wo[,4],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_af_all_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_all_wo[,1],rev(x_cum_af_all_wo[,1]+x_cum_af_all_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_all_wo[,1]+x_cum_af_all_wo[,2],rev(x_cum_af_all_wo[,1]+x_cum_af_all_wo[,2]+x_cum_af_all_wo[,3])),col=adjustcolor(colorRampPalette(colors)(4)[3]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_all_wo[,1]+x_cum_af_all_wo[,2]+x_cum_af_all_wo[,3],rev(x_cum_af_all_wo[,1]+x_cum_af_all_wo[,2]+x_cum_af_all_wo[,3]+x_cum_af_all_wo[,4])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.035,labels="gas",cex=1.2,col="black",pos=4)
text(0.15,0.1,labels="coal",cex=1.2,col="white",pos=4)
text(0.8,0.1,labels="nuclear",cex=1.2,col="white",pos=4)
text(5.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_af_GCW[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal and Wind With")
lines(beta,x_cum_af_GCW[,1]+x_cum_af_GCW[,2],type="l")
lines(beta,x_cum_af_GCW[,1]+x_cum_af_GCW[,2]+x_cum_af_GCW[,3],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_af_GCW[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_GCW[,1],rev(x_cum_af_GCW[,1]+x_cum_af_GCW[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_GCW[,1]+x_cum_af_GCW[,2],rev(x_cum_af_GCW[,1]+x_cum_af_GCW[,2]+x_cum_af_GCW[,3])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.05,labels="gas",cex=1.2,col="black",pos=4)
text(0.8,0.13,labels="coal",cex=1.2,col="white",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_af_GCW_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal and Wind Without")
lines(beta,x_cum_af_GCW_wo[,1]+x_cum_af_GCW_wo[,2],type="l")
lines(beta,x_cum_af_GCW_wo[,1]+x_cum_af_GCW_wo[,2]+x_cum_af_GCW_wo[,3],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_af_GCW_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_GCW_wo[,1],rev(x_cum_af_GCW_wo[,1]+x_cum_af_GCW_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_GCW_wo[,1]+x_cum_af_GCW_wo[,2],rev(x_cum_af_GCW_wo[,1]+x_cum_af_GCW_wo[,2]+x_cum_af_GCW_wo[,3])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.05,labels="gas",cex=1.2,col="black",pos=4)
text(0.3,0.1,labels="coal",cex=1.2,col="white",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_af_GW[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas and Wind With")
lines(beta,x_cum_af_GW[,1]+x_cum_af_GW[,2],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_af_GW[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_GW[,1],rev(x_cum_af_GW[,1]+x_cum_af_GW[,2])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.08,labels="gas",cex=1.2,col="black",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_af_GW_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas and Wind Without")
lines(beta,x_cum_af_GW_wo[,1]+x_cum_af_GW_wo[,2],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_af_GW_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_GW_wo[,1],rev(x_cum_af_GW_wo[,1]+x_cum_af_GW_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.05,labels="gas",cex=1.2,col="black",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

# PDF

pdf("Shares_EF_AF.pdf",width=14,height=8)
par(mar=c(5,5,5,5))
par(mfrow=c(3,2))
plot(beta,x_cum_af_all[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal, Nuclear and Wind With")
lines(beta,x_cum_af_all[,1]+x_cum_af_all[,2],type="l")
lines(beta,x_cum_af_all[,1]+x_cum_af_all[,2]+x_cum_af_all[,3],type="l")
lines(beta,x_cum_af_all[,1]+x_cum_af_all[,2]+x_cum_af_all[,3]+x_cum_af_all[,4],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_af_all[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_all[,1],rev(x_cum_af_all[,1]+x_cum_af_all[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_all[,1]+x_cum_af_all[,2],rev(x_cum_af_all[,1]+x_cum_af_all[,2]+x_cum_af_all[,3])),col=adjustcolor(colorRampPalette(colors)(4)[3]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_all[,1]+x_cum_af_all[,2]+x_cum_af_all[,3],rev(x_cum_af_all[,1]+x_cum_af_all[,2]+x_cum_af_all[,3]+x_cum_af_all[,4])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.035,labels="gas",cex=1.2,col="black",pos=4)
text(3.5,0.4,labels="nuclear",cex=1.2,col="white",pos=4)
text(1.2,0.035,labels="coal",cex=1.2,col="white",pos=4)
text(6.5,0.83,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_af_all_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal, Nuclear and Wind Without")
lines(beta,x_cum_af_all_wo[,1]+x_cum_af_all_wo[,2],type="l")
lines(beta,x_cum_af_all_wo[,1]+x_cum_af_all_wo[,2]+x_cum_af_all_wo[,3],type="l")
lines(beta,x_cum_af_all_wo[,1]+x_cum_af_all_wo[,2]+x_cum_af_all_wo[,3]+x_cum_af_all_wo[,4],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_af_all_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_all_wo[,1],rev(x_cum_af_all_wo[,1]+x_cum_af_all_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_all_wo[,1]+x_cum_af_all_wo[,2],rev(x_cum_af_all_wo[,1]+x_cum_af_all_wo[,2]+x_cum_af_all_wo[,3])),col=adjustcolor(colorRampPalette(colors)(4)[3]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_all_wo[,1]+x_cum_af_all_wo[,2]+x_cum_af_all_wo[,3],rev(x_cum_af_all_wo[,1]+x_cum_af_all_wo[,2]+x_cum_af_all_wo[,3]+x_cum_af_all_wo[,4])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.035,labels="gas",cex=1.2,col="black",pos=4)
text(0.15,0.1,labels="coal",cex=1.2,col="white",pos=4)
text(0.8,0.1,labels="nuclear",cex=1.2,col="white",pos=4)
text(5.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_af_GCW[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal and Wind With")
lines(beta,x_cum_af_GCW[,1]+x_cum_af_GCW[,2],type="l")
lines(beta,x_cum_af_GCW[,1]+x_cum_af_GCW[,2]+x_cum_af_GCW[,3],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_af_GCW[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_GCW[,1],rev(x_cum_af_GCW[,1]+x_cum_af_GCW[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_GCW[,1]+x_cum_af_GCW[,2],rev(x_cum_af_GCW[,1]+x_cum_af_GCW[,2]+x_cum_af_GCW[,3])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.05,labels="gas",cex=1.2,col="black",pos=4)
text(0.8,0.13,labels="coal",cex=1.2,col="white",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_af_GCW_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas, Coal and Wind Without")
lines(beta,x_cum_af_GCW_wo[,1]+x_cum_af_GCW_wo[,2],type="l")
lines(beta,x_cum_af_GCW_wo[,1]+x_cum_af_GCW_wo[,2]+x_cum_af_GCW_wo[,3],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_af_GCW_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_GCW_wo[,1],rev(x_cum_af_GCW_wo[,1]+x_cum_af_GCW_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[2]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_GCW_wo[,1]+x_cum_af_GCW_wo[,2],rev(x_cum_af_GCW_wo[,1]+x_cum_af_GCW_wo[,2]+x_cum_af_GCW_wo[,3])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.05,labels="gas",cex=1.2,col="black",pos=4)
text(0.3,0.1,labels="coal",cex=1.2,col="white",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_af_GW[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas and Wind With")
lines(beta,x_cum_af_GW[,1]+x_cum_af_GW[,2],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_af_GW[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_GW[,1],rev(x_cum_af_GW[,1]+x_cum_af_GW[,2])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.08,labels="gas",cex=1.2,col="black",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)

plot(beta,x_cum_af_GW_wo[,1],type="l",cex.axis=2,cex.lab=2,cex.main=2,ylab="%",xlab=TeX('$\\beta$'),ylim=c(0,1),xlim=c(0,beta[length(beta)]),main="Gas and Wind Without")
lines(beta,x_cum_af_GW_wo[,1]+x_cum_af_GW_wo[,2],type="l")
polygon(c(beta,rev(beta)),c(rep(0,length(beta)),rev(x_cum_af_GW_wo[,1])),col=adjustcolor(colorRampPalette(colors)(4)[1]),border="white")
polygon(c(beta,rev(beta)),c(x_cum_af_GW_wo[,1],rev(x_cum_af_GW_wo[,1]+x_cum_af_GW_wo[,2])),col=adjustcolor(colorRampPalette(colors)(4)[4]),border="white")
abline(h=0)
abline(h=1)
text(0,0.05,labels="gas",cex=1.2,col="black",pos=4)
text(4.5,0.6,labels="wind",cex=1.2,col="white",pos=4)
dev.off()

# == Shares in Each Load Block ==

# = beta=0 =

# R

par(mar=c(2,2,2,2))
par(mfrow=c(3,2))
plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_af_all[1,16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_af_all[1,16]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_af_all[1,16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_af_all[1,16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Nuclear","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_af_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_af_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_af_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_af_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Coal","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_af_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_af_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_af_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_af_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)


# PDF

CairoPDF("Shares_LB_AF_0.pdf",width=14,height=8)
par(mar=c(5,5,5,5))
par(mfrow=c(3,2))
plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_af_all[1,16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_af_all[1,16]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_af_all[1,16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_af_all[1,16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Nuclear","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_af_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_af_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_af_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_af_GCW[1,11]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Coal","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_af_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_af_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_af_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_af_GW[1,6]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
dev.off()

# = beta=10 =

# R

par(mar=c(2,2,2,2))
par(mfrow=c(3,2))
plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Nuclear","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_af_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_af_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_af_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_af_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Coal","Nuclear","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[3]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-3]+x_rel_af_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_af_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3]+x_rel_af_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_af_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_af_GCW[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_af_GCW[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_af_GCW[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_af_GCW[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Coal","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_af_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_af_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_af_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_af_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Coal","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-3]+x_rel_af_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_af_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3]+x_rel_af_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_af_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_af_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_af_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_af_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_af_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)


# PDF

CairoPDF("Shares_LB_AF_10.pdf",width=14,height=8)
par(mar=c(5,5,5,5))
par(mfrow=c(3,2))
plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,1,1,0),c(LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N]+x_rel_all[length(beta),16]*(LDC.d[2*N-1]-LDC.d[2*N]),LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Nuclear","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[3]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal, Nuclear and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_af_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_af_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[3]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_af_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_af_all_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Coal","Nuclear","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[3]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-3]+x_rel_af_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_af_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3]+x_rel_af_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_af_GCW[length(beta),13]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_af_GCW[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_af_GCW[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_af_GCW[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_af_GCW[length(beta),9]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Coal","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas, Coal and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_af_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_af_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_af_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_af_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[2]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Coal","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[2]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind With",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-3]+x_rel_af_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_af_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3]+x_rel_af_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-3]+x_rel_af_GW[length(beta),8]*(LDC.d[2*N-5]-LDC.d[2*N-3]),LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)

plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2,type="l",xlab="\u2113",ylab="D(\u2113)",main="Gas and Wind Without",cex.main=2,cex.axis=2,cex.lab=2)
for (i in 1:(2*N-3))
{polygon(c(0,8760/8760*d[2*i],8760/8760*d[2*i],0),c(LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i+1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1]),LDC6.poly(d[2*i-1])/LDC6.poly(d[1])),col=col_bg[i])}
{polygon(c(0,8760/8760*d[2*N],8760/8760*d[2*N],0),c(LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1]),LDC6.poly(d[2*N-1])/LDC6.poly(d[1])),col=col_bg[i])}
polygon(c(0,1,1,0),c(LDC.d[2*N],LDC.d[2*N],LDC.d[2*N-1],LDC.d[2*N-1]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-2],d[2*N-2],0),c(LDC.d[2*N-1],LDC.d[2*N-1],LDC.d[2*N-3],LDC.d[2*N-3]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-4],d[2*N-4],0),c(LDC.d[2*N-3],LDC.d[2*N-3],LDC.d[2*N-5],LDC.d[2*N-5]),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5],LDC.d[2*N-5],LDC.d[2*N-5]+x_rel_af_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_af_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5])),col=adjustcolor(colorRampPalette(colors)(4)[4]))
polygon(c(0,d[2*N-6],d[2*N-6],0),c(LDC.d[2*N-5]+x_rel_af_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-5]+x_rel_af_GCW_wo[length(beta),14]*(LDC.d[2*N-7]-LDC.d[2*N-5]),LDC.d[2*N-7],LDC.d[2*N-7]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
polygon(c(0,d[2],d[2],0),c(LDC.d[2*N-7],LDC.d[2*N-7],LDC.d[1],LDC.d[1]),col=adjustcolor(colorRampPalette(colors)(4)[1]))
lines(seq(0,1,0.0001141553),LDC.norm,lwd=2,col="black")
lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=2)
legend("topright",c("Gas","Wind"),col=c(adjustcolor(colorRampPalette(colors)(4)[1]),adjustcolor(colorRampPalette(colors)(4)[4])),lwd=4)
dev.off()