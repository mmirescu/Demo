# ===== Optimal Discretisation of the Load Duration Curve =====

# Set working directory to source file location

setwd("C:/Users/Mirescu1/ucloud/Universities/TU Wien/Theses/Dissertation/Eigene/Latex")

# Install necessary packages

library(rootSolve)
library("colorspace")
library(grDevices)
library(GoFKernel)
library(latex2exp)
library(tikzDevice)
library(extrafont)
library(Cairo)

# === Define colors to be used throughout this document ===

# Call up on interface where one can define a personalized color palette with subsequent color generator in HEX-format

#choose_palette()

# Introduce the personalized palettes into R

# Nikolaus Rab (nr) color palette

col_nr<-c("#E5FCC2","#9DE0AD","#45ADA8","#547980","#594F4F")

# Blue-Green (bg) color palette

col_bg<-c("#00584D","#00685E","#347B73","#65948E","#A1BAB6","#A3B9BC","#689398","#387A81","#00666E","#00565F")

# Pink-Purple (pp) color palette 

col_pp<-rev(c("#CA2497","#D456A6","#DF7CB7","#EAA4CC","#F6D5E7","#F0D6EF","#DEA7DD","#D181CF","#C65BC3","#BD2ABA"))

# Coray-Pink (cp) color palette

col_cp<-rev(c("#FF7C89","#FF97A0","#FFB0B6","#FCC8CB","#F2E0E1","#F2DFE4","#FDC7D6","#FFADC8","#FF91BA","#FF73AD"))

# Blue 30 (b30) color palette

col_b30<-rev(c("#2D3184","#254289","#1C518F","#115F96","#046C9C","#0078A2","#0084A7","#068FAB","#189AAF","#28A4B3","#37ADB6","#45B6B8","#54BEBA","#62C5BC",
  "#70CCBD","#7DD2BF","#8BD8C0","#97DDC1","#A3E2C3","#AFE5C4","#BAE9C6","#C4EBC8","#CDEECA","#D6F0CD","#DEF1D0","#E4F2D3","#EAF2D6",
  "#EEF2DA","#F2F2DE","#F3F1E4"))

# Blue 50 [MAXIMUM] (b50) color palette

col_b50<-rev(c("#023FA5","#1041A5","#1A44A5","#2146A5","#2749A5","#2C4BA6","#314EA6","#3550A7","#3953A8","#3D56A9","#4158AA","#455BAB","#485DAC","#4C60AD",
               "#4F63AE","#5365AF","#5668B0","#5A6BB2","#5D6EB3","#6170B5","#6473B6","#6876B8","#6B79B9","#6E7CBB","#727FBC","#7582BE","#7985C0","#7D88C1",
               "#808BC3","#848FC5","#8892C7","#8B95C9","#8F99CB","#939CCD","#97A0CF","#9BA3D1","#9FA7D3","#A3ABD5","#A8AFD7","#ACB3DA","#B1B7DC","#B6BCDF",
               "#BBC0E1","#C0C5E4","#C5CAE7","#CBD0EA","#D2D5ED","#D9DCF1","#E1E4F5","#EEF0FC"))

# Define the colors used in the plot to be of the preferred palette

col_bg<-col_b50

# Pie Chart of German Electricity Generation by Technology

gen<-c(22.5,12.9,11.8,12.9,0.8,14.3,3,2.6,7.1,7.1,1.0,4.1)
piepercent<-round(100*gen/sum(gen),1)
lab<-c("Lignite","Coal","Nuclear","Gas","Oil","Wind Onshore","Wind Offshore","Hydro","Biomass","Solar PV","Household Wase","Other")
colors<-c(col_b50[1],col_b50[5],col_b50[9],col_b50[13],col_b50[17],col_b50[21],col_b50[25],col_b50[29],col_b50[33],col_b50[37],col_b50[41],col_b50[45])
pdf("DE_gen.pdf",width=14,height=8)
par(mar=c(5,5,5,5))
pie(gen,lab=piepercent,col=colors,cex=1.5)
legend("topright",c("Lignite","Coal","Nuclear","Gas","Oil","Wind Onshore","Wind Offshore","Hydro","Biomass","Solar PV","Household Waste","Other"),col=c(col_b50[1],col_b50[5],col_b50[9],col_b50[13],col_b50[17],col_b50[21],col_b50[25],col_b50[29],col_b50[33],col_b50[37],col_b50[41],col_b50[45]),lwd=3,cex=1.5)
dev.off()

gen<-c(49.1,11.8,17.3,7.1,9.7,5.1)
piepercent<-round(100*gen/sum(gen),1)
lab<-c("Fossil (49.1%)","Nuklear (11.8%)","Wind (17.3%)","Solar PV (7.1%)","Andere RES (9.7%)","Andere (5.1%)")
colors<-c(col_bg[1],col_bg[2],col_bg[3],col_bg[4],col_bg[5],col_bg[6])
pdf("DE_erz.pdf",width=14,height=8)
par(mar=c(5,5,5,5))
pie(gen,lab=piepercent,col=colors,cex=1.5)
legend("topright",c("Fossil (49.1%)","Nuklear (11.8%)","Wind (17.3%)","Solar PV (7.1%)","Andere RES (9.7%)","Andere (5.1%)"),col=c(col_bg[1],col_bg[2],col_bg[3],col_bg[4],col_bg[5],col_bg[6]),lwd=3,cex=1.5)
dev.off()


# === Read in the necessary data ===

# Load Germany 2015

LoadDE15<-read.csv("LoadDE15_new.csv",sep=";")
load15<-as.vector(unlist(LoadDE15))

# Load Germany 2016

LoadDE16<-read.csv("LoadDE16.csv",sep=";")
load16<-as.vector(unlist(LoadDE16))

# Load Germany 2017

LoadDE<-read.csv("LoadDE17.csv",sep=";")
load<-as.vector(unlist(LoadDE))

# === Compute the Empirical Load Duration Curve ===

# 2015

LDC15.emp<-load15[rev(order(load15))]

# 2016

LDC16.emp<-load16[rev(order(load16))]

# 2017

LDC.emp<-load[rev(order(load))]

# # Plot
# 
# # R
# 
# par(mar=c(4,4,4,4))
# par(mfrow=c(3,1))
# plot(load15*1e-3,type="l",lwd=2.5,col="darkblue",xlab="Time [h]",ylab="Load [GW]",main="2015",ylim=c(min(load15*1e-3),max(load15*1e-3)),cex.main=2.5,cex.axis=2,cex.lab=2)
# plot(load16*1e-3,type="l",lwd=2.5,col="blue",xlab="Time [h]",ylab="Load [GW]",main="2016",ylim=c(min(load16*1e-3),max(load16*1e-3)),cex.main=2.5,cex.axis=2,cex.lab=2)
# plot(load*1e-3,type="l",lwd=2.5,col="lightblue",xlab="Time [h]",ylab="Load [GW]",main="2017",ylim=c(min(load*1e-3),max(load*1e-3)),cex.main=2.5,cex.axis=2,cex.lab=2)
# 
# 
# pdf("load17.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# plot(load*1e-3,type="l",lwd=2.5,col="lightblue",xlab="Time [h]",ylab="Load [GW]",main="2017",ylim=c(min(load*1e-3),max(load*1e-3)),cex.main=2.5,cex.axis=2,cex.lab=2)
# dev.off()


# # PDF
# #
# pdf("Load_DE_151617.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# par(mfrow=c(3,1))
# plot(load15*1e-3,type="l",lwd=2.5,col="darkblue",xlab="Time [h]",ylab="Load [GW]",main="2015",ylim=c(min(load15*1e-3),max(load15*1e-3)),cex.main=3,cex.axis=2,cex.lab=2)
# plot(load16*1e-3,type="l",lwd=2.5,col="blue",xlab="Time [h]",ylab="Load [GW]",main="2016",ylim=c(min(load16*1e-3),max(load16*1e-3)),cex.main=3,cex.axis=2,cex.lab=2)
# plot(load*1e-3,type="l",lwd=2.5,col="lightblue",xlab="Time [h]",ylab="Load [GW]",main="2017",ylim=c(min(load*1e-3),max(load*1e-3)),cex.main=3,cex.axis=2,cex.lab=2)
# dev.off()

# === Estimation of a parametric polynomial Load Duration Curve ===

# Normalize the empirical Load Duration Curve

# 2015

LDC15.norm<-LDC15.emp/max(LDC15.emp)

# 2016

LDC16.norm<-LDC16.emp/max(LDC16.emp)

# 2017

LDC.norm<-LDC.emp/max(LDC.emp)

# pdf("LDC17.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# plot(LDC.emp,type="l",lwd=2.5,col="lightblue",xlab="Duration t [h]",ylab="Load [MW]",main="Empirical LDC (ELDC)",cex.main=2,cex.axis=2,cex.lab=2)
# grid(10,10,lwd=1.5)
# lines(LDC.emp,type="l",lwd=2.5,col="lightblue")
# dev.off()

# # Plot
# 
# # R
# 
# par(mfrow=c(1,2))
# plot(LDC15.emp*1e-3,type="l",lwd=2.5,col="darkblue",xlab="Duration t [h]",ylab="Load [GW]",main="(a) Empirical LDC (ELDC)",cex.main=2,cex.axis=2,cex.lab=2)
# lines(LDC16.emp*1e-3,type="l",lwd=2.5,col="blue")
# lines(LDC.emp*1e-3,type="l",lwd=2.5,col="lightblue")
# #legend("topright",c("2015","2016","2017"),col=c("darkblue","blue","lightblue"),lwd=4,cex=2)
# plot(LDC15.norm,type="l",lwd=2.5,col="darkblue",xlab="Load Factor l",ylab="Normalised Load",main="(b) Normalised ELDC",cex.main=2,cex.axis=2,cex.lab=2)
# lines(LDC16.norm,type="l",lwd=2.5,col="blue")
# lines(LDC.norm,type="l",lwd=2.5,col="lightblue")
# #legend("topright",c("2015","2016","2017"),col=c("darkblue","blue","lightblue"),lwd=4,cex=2)

# # PDF
# 
# pdf("ELDC_DE_151617.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# par(mfrow=c(1,2))
# plot(LDC15.emp*1e-3,type="l",lwd=2.5,col="darkblue",xlab="Duration t [h]",ylab="Load [GW]",main="(a) Empirical LDC (ELDC)",cex.main=2,cex.axis=2,cex.lab=2)
# lines(LDC16.emp*1e-3,type="l",lwd=2.5,col="blue")
# lines(LDC.emp*1e-3,type="l",lwd=2.5,col="lightblue")
# legend("topright",c("2015","2016","2017"),col=c("darkblue","blue","lightblue"),lwd=4,cex=2)
# plot(LDC15.norm,type="l",lwd=2.5,col="darkblue",xlab="Load Factor l",ylab="Normalised Load",main="(b) Normalised ELDC",cex.main=2,cex.axis=2,cex.lab=2)
# lines(LDC16.norm,type="l",lwd=2.5,col="blue")
# lines(LDC.norm,type="l",lwd=2.5,col="lightblue")
# legend("topright",c("2015","2016","2017"),col=c("darkblue","blue","lightblue"),lwd=4,cex=2)
# dev.off()

# Load Factors

LF<-(0:(length(LDC.emp)-1))/(length(LDC.emp)-1)

## OLS Regression assuming a linear polynomial

#LDC1.regression<-nls(LDC.norm~1+b*LF,start=list(b=1))
#LDC1.coef<-summary(LDC1.regression)$coef[,1]

## Polynomial LDC (t in [0,1])

#LDC1.poly<-function(t) max(LDC.emp)*(1+sum(sapply(1,function(i) LDC1.coef[i]*t^i)))

## Plot

#par(mfrow=c(1,1))
#plot(seq(0,1,0.01),sapply(seq(0,1,0.01),LDC1.poly)/sapply(seq(0,1,0.01),LDC1.poly)[1],col="lightblue",lwd=3,type="l",xlab="Load Factor",ylab="Normalised Load [MW]",cex.main=2,cex.axis=2,cex.lab=2)
#lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="darkblue")

## OLS Regression assuming a quadratic polynomial

#LDC2.regression<-nls(LDC.norm~1+b*LF+c*LF^2,start=list(b=1,c=1))
#LDC2.coef<-summary(LDC2.regression)$coef[,1]

## Polynomial LDC (t in [0,1])

#LDC2.poly<-function(t) max(LDC.emp)*(1+sum(sapply(1:2,function(i) LDC2.coef[i]*t^i)))

## Plot

#par(mfrow=c(1,1))
#plot(seq(0,1,0.01),sapply(seq(0,1,0.01),LDC2.poly)/sapply(seq(0,1,0.01),LDC2.poly)[1],col="lightblue",lwd=3,type="l",xlab="Load Factor",ylab="Normalised Load [MW]",cex.main=2,cex.axis=2,cex.lab=2)
#lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="darkblue")

## OLS Regression assuming a cubic polynomial

#LDC3.regression<-nls(LDC.norm~1+b*LF+c*LF^2+d*LF^3,start=list(b=1,c=1,d=1))
#LDC3.coef<-summary(LDC3.regression)$coef[,1]

## Polynomial LDC (t in [0,1])

#LDC3.poly<-function(t) max(LDC.emp)*(1+sum(sapply(1:3,function(i) LDC3.coef[i]*t^i)))

## Plot

#par(mfrow=c(1,1))
#plot(seq(0,1,0.01),sapply(seq(0,1,0.01),LDC3.poly)/sapply(seq(0,1,0.01),LDC3.poly)[1],col="lightblue",lwd=3,type="l",xlab="Load Factor",ylab="Normalised Load [MW]",cex.main=2,cex.axis=2,cex.lab=2)
#lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="darkblue")

## OLS Regression assuming a fourth-degree polynomial

#LDC4.regression<-nls(LDC.norm~1+b*LF+c*LF^2+d*LF^3+e*LF^4,start=list(b=1,c=1,d=1,e=1))
#LDC4.coef<-summary(LDC4.regression)$coef[,1]

## Polynomial LDC (t in [0,1])

#LDC4.poly<-function(t) max(LDC.emp)*(1+sum(sapply(1:4,function(i) LDC4.coef[i]*t^i)))

## Plot

#par(mfrow=c(1,1))
#plot(seq(0,1,0.01),sapply(seq(0,1,0.01),LDC4.poly)/sapply(seq(0,1,0.01),LDC4.poly)[1],col="lightblue",lwd=3,type="l",xlab="Load Factor",ylab="Normalised Load [MW]",cex.main=2,cex.axis=2,cex.lab=2)
#lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="darkblue")

## OLS Regression assuming a fifth-degree polynomial

#LDC5.regression<-nls(LDC.norm~1+b*LF+c*LF^2+d*LF^3+e*LF^4+f*LF^5,start=list(b=1,c=1,d=1,e=1,f=1))
#LDC5.coef<-summary(LDC5.regression)$coef[,1]

## Polynomial LDC (t in [0,1])

#LDC5.poly<-function(t) max(LDC.emp)*(1+sum(sapply(1:5,function(i) LDC5.coef[i]*t^i)))

## Plot

#par(mfrow=c(1,1))
#plot(seq(0,1,0.01),sapply(seq(0,1,0.01),LDC5.poly)/sapply(seq(0,1,0.01),LDC5.poly)[1],col="lightblue",lwd=3,type="l",xlab="Load Factor",ylab="Normalised Load [MW]",cex.main=2,cex.axis=2,cex.lab=2)
#lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="darkblue")


# OLS Regression assuming a sixth-degree polynomial


LDC6.regression<-nls(LDC.norm~1+b*LF+c*LF^2+d*LF^3+e*LF^4+f*LF^5+g*LF^6,start=list(b=1,c=1,d=1,e=1,f=1,g=1))
LDC6.coef<-summary(LDC6.regression)$coef[,1]

# Polynomial LDC (t in [0,1])

LDC6.poly<-function(t) max(LDC.emp)*(1+sum(sapply(1:6,function(i) LDC6.coef[i]*t^i)))


# # Plot
# 
# par(mfrow=c(1,1))
# plot(seq(0,1,0.01),sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="lightblue",lwd=3,type="l",xlab="Load Factor",ylab="Normalised Load [MW]",cex.main=2,cex.axis=2,cex.lab=2)
# lines(seq(0,1,0.0001141553),LDC15.norm,lwd=3,col="darkblue")

## OLS Regression assuming a seventh-degree polynomial

#LDC7.regression<-nls(LDC.norm~1+b*LF+c*LF^2+d*LF^3+e*LF^4+f*LF^5+g*LF^6+h*LF^7,start=list(b=1,c=1,d=1,e=1,f=1,g=1,h=1))
#LDC7.coef<-summary(LDC7.regression)$coef[,1]

## Polynomial LDC (t in [0,1])

#LDC7.poly<-function(t) max(LDC.emp)*(1+sum(sapply(1:7,function(i) LDC7.coef[i]*t^i)))

## Plot

#par(mfrow=c(1,1))
#plot(seq(0,1,0.01),sapply(seq(0,1,0.01),LDC7.poly)/sapply(seq(0,1,0.01),LDC7.poly)[1],col="lightblue",lwd=3,type="l",xlab="Load Factor",ylab="Normalised Load [MW]",cex.main=2,cex.axis=2,cex.lab=2)
#lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="darkblue")

## Plot all estimation approaches

#par(mfrow=c(1,1))
#plot(seq(0,1,0.01),sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="lightblue",lwd=3,type="l",xlab="Load Factor",ylab="Normalised Load [MW]",cex.main=2,cex.axis=2,cex.lab=2)
#lines(seq(0,1,0.01),sapply(seq(0,1,0.01),LDC5.poly)/sapply(seq(0,1,0.01),LDC5.poly)[1],col="blue",lwd=3,type="l")
#lines(seq(0,1,0.01),sapply(seq(0,1,0.01),LDC4.poly)/sapply(seq(0,1,0.01),LDC4.poly)[1],col="darkblue",lwd=3,type="l")
#lines(seq(0,1,0.01),sapply(seq(0,1,0.01),LDC3.poly)/sapply(seq(0,1,0.01),LDC3.poly)[1],col="purple",lwd=3,type="l")
#lines(seq(0,1,0.01),sapply(seq(0,1,0.01),LDC2.poly)/sapply(seq(0,1,0.01),LDC2.poly)[1],col="cyan",lwd=3,type="l")
#lines(seq(0,1,0.01),sapply(seq(0,1,0.01),LDC1.poly)/sapply(seq(0,1,0.01),LDC1.poly)[1],col="green",lwd=3,type="l")
#lines(seq(0,1,0.01),sapply(seq(0,1,0.01),LDC7.poly)/sapply(seq(0,1,0.01),LDC7.poly)[1],col="orange",lwd=3,type="l")
#lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="red")


# === Newton-Raphson Root Finding Procedure ===

# Define the parameters needed

a<-as.vector(unlist(c(1,LDC6.coef)))
D_max<-max(LDC.norm)
D_min<-min(LDC.norm)

###########################################################################################################################################################

# Determine the Desired Number of Load Blocks

#N=3  # Suppress N Here

# === Approach 1: Horizontal Equidistant Discretisation ===

# = Upper Sums =

delta_h_up=1/N

# Error

integrand<-function(x) LDC6.poly(x)/LDC6.poly(0)
err_int_up_he<-sum(sapply(1:N, function(i) delta_h_up*LDC6.poly((i-1)*delta_h_up)/LDC6.poly(0)))-integrate(Vectorize(integrand),lower=0,upper=1)$value
#print(err_int_up_he)

ell_h_up<-0
f_ell_h_up<-LDC6.poly(0)/LDC6.poly(0)
for (i in 1:N) {
  ell_h_up<-c(ell_h_up,i*delta_h_up)
  f_ell_h_up<-c(f_ell_h_up,LDC6.poly(i*delta_h_up)/LDC6.poly(0))
}

# = Alternating Sums =

delta_h_alt=1/(2*N-1)

# Error 

err_int_alt_he<-sum(sapply(1:(2*N-1), function(i) (-1)^(i+1)*delta_h_alt*LDC6.poly(2*ceiling((i-1)/2)*delta_h_alt)/LDC6.poly(0)+(-1)^i*integrate(Vectorize(integrand),lower=(i-1)*delta_h_alt,upper=i*delta_h_alt)$value))
#print(err_int_alt_he)

ell_h_alt<-0
f_ell_h_alt<-LDC6.poly(0)/LDC6.poly(0)
for (i in 1:(2*N-1)) {
  ell_h_alt<-c(ell_h_alt,i*delta_h_alt)
  f_ell_h_alt<-c(f_ell_h_alt,LDC6.poly(i*delta_h_alt)/LDC6.poly(0))
}

# === Approach 2: Vertical Equidistant Approximation of D ===

# Define the Inverse of the Polynomial

LDC6.poly.inv<-inverse(integrand,lower=0,upper=1)

# = Upper Sums =

delta_v_up=(LDC6.poly(0)-LDC6.poly(1))/(N*LDC6.poly(0))

# Error

err_int_up_ve<-sum(sapply(1:N, function(i) (LDC6.poly.inv(integrand(1)+(N-i)*delta_v_up)-LDC6.poly.inv(integrand(1)+(N-i+1)*delta_v_up))*(integrand(1)+(N-i+1)*delta_v_up)-integrate(Vectorize(integrand),lower=LDC6.poly.inv(integrand(1)+(N-i+1)*delta_v_up),upper=LDC6.poly.inv(integrand(1)+(N-i)*delta_v_up))$value))
#print(err_int_up_ve)

ell_v_up<-0
f_ell_v_up<-LDC6.poly(0)/LDC6.poly(0)
for (i in 1:N) {
  ell_v_up<-c(ell_v_up,LDC6.poly.inv(integrand(1)+(N-i)*delta_v_up))
  f_ell_v_up<-c(f_ell_v_up,1-i*delta_v_up)
}

# = Alternating Sums =

delta_v_alt=(LDC6.poly(0)-LDC6.poly(1))/((2*N-1)*LDC6.poly(0))

# Error

err_int_alt_ve<-sum(sapply(1:(2*N-1), function(i) (-1)^(i+1)*(LDC6.poly.inv(integrand(1)+(2*N-1-i)*delta_v_alt)-LDC6.poly.inv(integrand(1)+(2*N-i)*delta_v_alt))*(integrand(1)+(2*N-1-2*ceiling((i-1)/2))*delta_v_alt)+(-1)^i*integrate(Vectorize(integrand),lower=LDC6.poly.inv(integrand(1)+(2*N-i)*delta_v_alt),upper=LDC6.poly.inv(integrand(1)+(2*N-1-i)*delta_v_alt))$value))
print(err_int_alt_ve)

ell_v_alt<-0
f_ell_v_alt<-LDC6.poly(0)/LDC6.poly(0)
for (i in 1:(2*N-1)) {
  ell_v_alt<-c(ell_v_alt,LDC6.poly.inv(integrand(1)+(2*N-1-i)*delta_v_alt))
  f_ell_v_alt<-c(f_ell_v_alt,1-i*delta_v_alt)
}

# === Approach 3: Optimal Discretisation of the Load Duration Curve ===

# = Upper Sums =

# Compute the functions whose roots need to be found later via Newton-Raphson

grad_LDC_int<-function(x,j) sum(sapply(1:6,function(i) a[i+1]*(x[j-1]^i-x[j]^i+(x[j+1]-x[j])*i*x[j]^(i-1))))

grad_LDC_combined_int<-function(x) {
  f<-NULL
  for (j in 2:(N-2)) f<-c(f,grad_LDC_int(x,j))
  f}

# Compute the Needed Roots via the Newton-Raphson Procedure

if (N==2){
  model_int<-function(x) c(sum(sapply(1:6,function(i) a[i+1]*(-x[1]^i+(1-x[1])*i*x[1]^(i-1)))))
  d_int<-c(0,uniroot(f=model_int,interval=c(0.1,0.9))$root,1)
  LDC_int.d<-sapply(d_int,LDC6.poly)/LDC6.poly(0)
}

if (N==3){
  model_int<-function(x) c(sum(sapply(1:6,function(i) a[i+1]*(-x[1]^i+(x[2]-x[1])*i*x[1]^(i-1)))),
                           sum(sapply(1:6,function(i) a[i+1]*(x[N-2]^i-x[N-1]^i+(1-x[N-1])*i*x[N-1]^(i-1)))))
  d_int<-c(0,multiroot(f=model_int,start=seq(0.3,0.9,length.out=(N-1)),positive=TRUE)$root,1)
  LDC_int.d<-sapply(d_int,LDC6.poly)/LDC6.poly(0)
}

if (N>=4){
  model_int<-function(x) c(D_max-1+sum(sapply(1:6,function(i) a[i+1]*(-x[1]^i+(x[2]-x[1])*i*x[1]^(i-1)))),
                       grad_LDC_combined_int(x),
                       sum(sapply(1:6,function(i) a[i+1]*(x[N-2]^i-x[N-1]^i+(1-x[N-1])*i*x[N-1]^(i-1)))))
  d_int<-c(0,multiroot(f=model_int,start=seq(0.4,0.8,length.out=(N-1)),positive=TRUE)$root,1)
  LDC_int.d<-sapply(d_int,LDC6.poly)/LDC6.poly(0)
}

# Compute the Error Relative to the Polynomial Load Duration Curve

err_int_up_opt=sum(sapply(1:N, function(i) (d_int[i+1]-d_int[i])*LDC_int.d[i]))-integrate(Vectorize(integrand),lower=0,upper=1)$value
#print(err_int_up_opt)

ell_opt_up<-d_int
f_ell_opt_up<-LDC_int.d


# = Alternating Sums =

# Compute the functions whose roots need to be found later via Newton-Raphson

grad_LDC_int_alt<-function(x,j) c((x[j-1]-2*x[j]+x[j+1])*sum(sapply(1:6,function(i) a[i+1]*i*x[j]^(i-1))),
                                  sum(sapply(1:6,function(i) a[i+1]*(x[j]^i-2*x[j+1]^i+x[j+2]^i))))

grad_LDC_combined_int_alt<-function(x) {
  f<-NULL
  for (j in 1:(N-2)) f<-c(f,grad_LDC_int_alt(x,2*(j-1)+2))
  f}

# Compute the needed roots via the Newton-Raphson Procedure


if (N==2){
  model_int_alt<-function(x) c(sum(sapply(1:6,function(i) a[i+1]*(-2*x[1]^i+x[2]^i))),
                               (x[1]-2*x[2]+1)*sum(sapply(1:6,function(i) a[i+1]*i*x[2]^(i-1))))
  d_int_alt<-c(0,multiroot(f=model_int_alt,start=seq(0.1,0.9,length.out=2*(N-1)),positive=TRUE)$root,1)
  LDC_int_alt.d<-sapply(d_int_alt,LDC6.poly)/LDC6.poly(0)
}

if (N==3){
  model_int_alt<-function(x) c(sum(sapply(1:6,function(i) a[i+1]*(-2*x[1]^i+x[2]^i))),
                               (x[1]-2*x[2]+x[3])*sum(sapply(1:6,function(i) a[i+1]*i*x[2]^(i-1))),
                               sum(sapply(1:6,function(i) a[i+1]*(x[2]^i-2*x[3]^i+x[4]^i))),
                               (x[3]-2*x[4]+1)*sum(sapply(1:6,function(i) a[i+1]*i*x[4]^(i-1))))
  d_int_alt<-c(0,multiroot(f=model_int_alt,start=seq(0.1,0.9,length.out=2*(N-1)),positive=TRUE)$root,1)
  LDC_int_alt.d<-sapply(d_int_alt,LDC6.poly)/LDC6.poly(0)
}

if (N>=4){
  model_int_alt<-function(x) c(sum(sapply(1:6,function(i) a[i+1]*(-2*x[1]^i+x[2]^i))),
                               grad_LDC_combined_int_alt(x),
                               (x[2*(N-1)-1]-2*x[2*(N-1)]+1)*sum(sapply(1:6,function(i) a[i+1]*i*x[2*(N-1)]^(i-1))))
  d_int_alt<-c(0,multiroot(f=model_int_alt,start=seq(0.1,0.9,length.out=2*(N-1)),positive=TRUE)$root,1)
  LDC_int_alt.d<-sapply(d_int_alt,LDC6.poly)/LDC6.poly(0)
}

# Error

err_int_alt_opt<-sum(sapply(1:(2*(N-1)), function(i) (-1)^(i+1)*(d_int_alt[i+1]-d_int_alt[i])*LDC6.poly(d_int_alt[1+2*ceiling((i-1)/2)])/LDC6.poly(0)+(-1)^i*integrate(Vectorize(integrand),lower=d_int_alt[i],upper=d_int_alt[i+1])$value))
#print(err_int_alt_opt)

ell_opt_alt<-d_int_alt
f_ell_opt_alt<-LDC_int_alt.d

# = Errors of the Six Approaches =

error<-cbind(c(err_int_up_he,err_int_alt_he),c(err_int_up_ve,err_int_alt_ve),c(err_int_up_opt,err_int_alt_opt))
colnames(error)<-c("HE","VE","OPT")
rownames(error)<-c("US","AS")
#View(error)

# = Mesh Values of all Approaches =

# Upper Sums 

mesh_up<-rbind(ell_h_up,ell_v_up,ell_opt_up)
rownames(mesh_up)<-c("HE","VE","OPT")
#print(mesh_up)

# Alternating Sums

mesh_alt<-rbind(ell_h_alt,ell_v_alt,ell_opt_alt)
rownames(mesh_alt)<-c("HE","VE","OPT")
#print(mesh_alt)

d<-d_int_alt
LDC.d<-LDC_int_alt.d

area_D<-integrate(Vectorize(integrand),lower=0,upper=1)$value

# = Function of Mesh Values of all Approaches =

# Upper Sums

f_mesh_up<-rbind(f_ell_h_up,f_ell_v_up,f_ell_opt_up)
rownames(f_mesh_up)<-c("HE","VE","OPT")
#print(f_mesh_up)

# Alternating Sums

f_mesh_alt<-rbind(f_ell_h_alt,f_ell_v_alt,f_ell_opt_alt)
rownames(f_mesh_alt)<-c("HE","VE","OPT")
#print(f_mesh_alt)
#View(rbind(mesh_up,f_mesh_up))
#View(rbind(mesh_alt,f_mesh_alt))

# = Overall Plot (of all the Approaches) =

# R

# #CairoPDF("results_N_3.pdf",width=14,height=8)
# CairoPDF("N_3.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# #par(mfrow=c(3,2))
# par(mfrow=c(1,2))
# # plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3,type="l",xlab="\u2113",ylab="D(\u2113)",main="(a) Equidistant in \u2113 (Darboux)",cex.main=2,cex.axis=2,cex.lab=2)
# # #grid(20,10,lwd=1)
# # for (i in 1:N) {polygon(c(0,8760/8760*i*delta_h_up,8760/8760*i*delta_h_up,0),c(LDC6.poly(i*delta_h_up)/LDC6.poly(0),LDC6.poly(i*delta_h_up)/LDC6.poly(0),LDC6.poly((i-1)*delta_h_up)/LDC6.poly(0),LDC6.poly((i-1)*delta_h_up)/LDC6.poly(0)),col=col_b50[1])}
# # lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="black")
# # lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3)
# # #axis(1,at=seq(0,1,by=0.2),las=2)
# # #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("darkblue","lightblue"),lwd=4,cex=2)
# # 
# # plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3,type="l",xlab="\u2113",ylab="D(\u2113)",main="(b) Equidistant in \u2113 (Riemann)",cex.main=2,cex.axis=2,cex.lab=2)
# # #grid(20,10,lwd=1)
# # {polygon(c(0,8760/8760*1*delta_h_alt,8760/8760*1*delta_h_alt,0),c(LDC6.poly(2*delta_h_alt)/max(LDC.emp),LDC6.poly(2*delta_h_alt)/LDC6.poly(0),LDC6.poly(0)/LDC6.poly(0),LDC6.poly(0)/LDC6.poly(0)),col=col_b50[1])}
# # for (i in 1:(N-2))
# # {polygon(c(0,8760/8760*(2*i+1)*delta_h_alt,8760/8760*(2*i+1)*delta_h_alt,0),c(LDC6.poly(2*(i+1)*delta_h_alt)/LDC6.poly(0),LDC6.poly(2*(i+1)*delta_h_alt)/LDC6.poly(0),LDC6.poly(2*i*delta_h_alt)/LDC6.poly(0),LDC6.poly(2*i*delta_h_alt)/LDC6.poly(0)),col=col_b50[1])}
# # {polygon(c(0,8760/8760*(2*N-1)*delta_h_alt,8760/8760*(2*N-1)*delta_h_alt,0),c(LDC6.poly((2*N-1)*delta_h_alt)/LDC6.poly(0),LDC6.poly((2*N-1)*delta_h_alt)/LDC6.poly(0),LDC6.poly(2*(N-1)*delta_h_alt)/LDC6.poly(0),LDC6.poly(2*(N-1)*delta_h_alt)/LDC6.poly(0)),col=col_b50[1])}
# # lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="black")
# # lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3)
# # #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("darkblue","lightblue"),lwd=4,cex=2)
# # 
# # plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3,type="l",xlab="\u2113",ylab="D(\u2113)",main="(c) Equidistant in D (Darboux)",cex.main=2,cex.axis=2,cex.lab=2)
# # #grid(20,10,lwd=1)
# # for (i in 1:N) {polygon(c(0,8760/8760*LDC6.poly.inv(integrand(1)+(N-i)*delta_v_up),8760/8760*LDC6.poly.inv(integrand(1)+(N-i)*delta_v_up),0),c(integrand(1)+(N-(i+1)+1)*delta_v_up,integrand(1)+(N-(i+1)+1)*delta_v_up,integrand(1)+(N-i+1)*delta_v_up,integrand(1)+(N-i+1)*delta_v_up),col=col_b50[1])}
# # lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="black")
# # lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3)
# # #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("darkblue","lightblue"),lwd=4,cex=2)
# # 
# # plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3,type="l",xlab="\u2113",ylab="D(\u2113)",main="(d) Equidistant in D (Riemann)",cex.main=2,cex.axis=2,cex.lab=2)
# # #grid(20,10,lwd=1)
# # {polygon(c(0,8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-1)*delta_v_alt),8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-1)*delta_v_alt),0),c(integrand(1)+(2*N-1-2*ceiling((2-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((2-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((1-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((1-1)/2))*delta_v_alt),col=col_b50[1])}
# # for (i in 1:(N-2))
# # {polygon(c(0,8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-2*(i+1)+1)*delta_v_alt),8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-2*(i+1)+1)*delta_v_alt),0),c(integrand(1)+(2*N-1-2*ceiling((2*(i+1)-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((2*(i+1)-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((2*i-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((2*i-1)/2))*delta_v_alt),col=col_b50[1])}
# # {polygon(c(0,8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-(2*N-1))*delta_v_alt),8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-(2*N-1))*delta_v_alt),0),c(integrand(1),integrand(1),integrand(1)+(2*N-1-2*ceiling((2*(N-1)-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((2*(N-1)-1)/2))*delta_v_alt),col=col_bg[i])}
# # lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="black")
# # lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3)
# # #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("darkblue","lightblue"),lwd=4,cex=2)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3,type="l",xlab="\u2113",ylab="D(\u2113)",main="Darboux",cex.main=2,cex.axis=2,cex.lab=2)
# #grid(20,10,lwd=1)
# for (i in 1:N) {polygon(c(0,8760/8760*d_int[i+1],8760/8760*d_int[i+1],0),c(LDC6.poly(d_int[i+1])/LDC6.poly(d_int[1]),LDC6.poly(d_int[i+1])/LDC6.poly(d_int[1]),LDC6.poly(d_int[i])/LDC6.poly(d_int[1]),LDC6.poly(d_int[i])/LDC6.poly(d_int[1])),col=col_b50[1])}
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3)
# #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("darkblue","lightblue"),lwd=4,cex=2)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3,type="l",xlab="\u2113",ylab="D(\u2113)",main="Riemann",cex.main=2,cex.axis=2,cex.lab=2)
# #grid(20,10,lwd=1)
# for (i in 1:(2*N-3))
# {polygon(c(0,8760/8760*d_int_alt[2*i],8760/8760*d_int_alt[2*i],0),c(LDC6.poly(d_int_alt[2*i+1])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*i+1])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*i-1])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*i-1])/LDC6.poly(d_int_alt[1])),col=col_b50[1])}
# {polygon(c(0,8760/8760*d_int_alt[2*N],8760/8760*d_int_alt[2*N],0),c(LDC6.poly(d_int_alt[2*N])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*N])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*N-1])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*N-1])/LDC6.poly(d_int_alt[1])),col=col_b50[1])}
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="black")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3)
# #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("darkblue","lightblue"),lwd=4,cex=2)
# dev.off()

# # PDF
# 
# # N=3
# 
# pdf("results_N_3.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# par(mfrow=c(3,2))
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="lightblue",lwd=3,type="l",xlab="l",ylab="D(\u2113)",main="(a) Equidistant in l (Darboux)",cex.main=2,cex.axis=2,cex.lab=2)
# #grid(20,10,lwd=1)
# for (i in 1:N) {polygon(c(0,8760/8760*i*delta_h_up,8760/8760*i*delta_h_up,0),c(LDC6.poly(i*delta_h_up)/LDC6.poly(0),LDC6.poly(i*delta_h_up)/LDC6.poly(0),LDC6.poly((i-1)*delta_h_up)/LDC6.poly(0),LDC6.poly((i-1)*delta_h_up)/LDC6.poly(0)),col=col_bg[i])}
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="darkblue")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="lightblue",lwd=3)
# #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("darkblue","lightblue"),lwd=4,cex=2)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="lightblue",lwd=3,type="l",xlab="l",ylab="D(\u2113)",main="(b) Equidistant in l (Riemann)",cex.main=2,cex.axis=2,cex.lab=2)
# #grid(20,10,lwd=1)
# {polygon(c(0,8760/8760*1*delta_h_alt,8760/8760*1*delta_h_alt,0),c(LDC6.poly(2*delta_h_alt)/max(LDC.emp),LDC6.poly(2*delta_h_alt)/LDC6.poly(0),LDC6.poly(0)/LDC6.poly(0),LDC6.poly(0)/LDC6.poly(0)),col=col_bg[i])}
# for (i in 1:(N-2)) 
# {polygon(c(0,8760/8760*(2*i+1)*delta_h_alt,8760/8760*(2*i+1)*delta_h_alt,0),c(LDC6.poly(2*(i+1)*delta_h_alt)/LDC6.poly(0),LDC6.poly(2*(i+1)*delta_h_alt)/LDC6.poly(0),LDC6.poly(2*i*delta_h_alt)/LDC6.poly(0),LDC6.poly(2*i*delta_h_alt)/LDC6.poly(0)),col=col_bg[i])}
# {polygon(c(0,8760/8760*(2*N-1)*delta_h_alt,8760/8760*(2*N-1)*delta_h_alt,0),c(LDC6.poly((2*N-1)*delta_h_alt)/LDC6.poly(0),LDC6.poly((2*N-1)*delta_h_alt)/LDC6.poly(0),LDC6.poly(2*(N-1)*delta_h_alt)/LDC6.poly(0),LDC6.poly(2*(N-1)*delta_h_alt)/LDC6.poly(0)),col=col_bg[i])}
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="darkblue")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="lightblue",lwd=3)
# #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("darkblue","lightblue"),lwd=4,cex=2)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="lightblue",lwd=3,type="l",xlab="l",ylab="D(\u2113)",main="(c) Equidistant in D (Darboux)",cex.main=2,cex.axis=2,cex.lab=2)
# #grid(20,10,lwd=1)
# for (i in 1:N) {polygon(c(0,8760/8760*LDC6.poly.inv(integrand(1)+(N-i)*delta_v_up),8760/8760*LDC6.poly.inv(integrand(1)+(N-i)*delta_v_up),0),c(integrand(1)+(N-(i+1)+1)*delta_v_up,integrand(1)+(N-(i+1)+1)*delta_v_up,integrand(1)+(N-i+1)*delta_v_up,integrand(1)+(N-i+1)*delta_v_up),col=col_bg[i])}
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="darkblue")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="lightblue",lwd=3)
# #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("darkblue","lightblue"),lwd=4,cex=2)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="lightblue",lwd=3,type="l",xlab="l",ylab="D(\u2113)",main="(d) Equidistant in D (Riemann)",cex.main=2,cex.axis=2,cex.lab=2)
# #grid(20,10,lwd=1)
# {polygon(c(0,8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-1)*delta_v_alt),8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-1)*delta_v_alt),0),c(integrand(1)+(2*N-1-2*ceiling((2-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((2-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((1-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((1-1)/2))*delta_v_alt),col=col_bg[i])}
# for (i in 1:(N-2)) 
# {polygon(c(0,8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-2*(i+1)+1)*delta_v_alt),8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-2*(i+1)+1)*delta_v_alt),0),c(integrand(1)+(2*N-1-2*ceiling((2*(i+1)-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((2*(i+1)-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((2*i-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((2*i-1)/2))*delta_v_alt),col=col_bg[i])}
# {polygon(c(0,8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-(2*N-1))*delta_v_alt),8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-(2*N-1))*delta_v_alt),0),c(integrand(1),integrand(1),integrand(1)+(2*N-1-2*ceiling((2*(N-1)-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((2*(N-1)-1)/2))*delta_v_alt),col=col_bg[i])}
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="darkblue")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="lightblue",lwd=3)
# #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("darkblue","lightblue"),lwd=4,cex=2)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="lightblue",lwd=3,type="l",xlab="l",ylab="D(\u2113)",main="(e) Optimal (Darboux)",cex.main=2,cex.axis=2,cex.lab=2)
# #grid(20,10,lwd=1)
# for (i in 1:N) {polygon(c(0,8760/8760*d_int[i+1],8760/8760*d_int[i+1],0),c(LDC6.poly(d_int[i+1])/LDC6.poly(d_int[1]),LDC6.poly(d_int[i+1])/LDC6.poly(d_int[1]),LDC6.poly(d_int[i])/LDC6.poly(d_int[1]),LDC6.poly(d_int[i])/LDC6.poly(d_int[1])),col=col_bg[i])}
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="darkblue")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="lightblue",lwd=3)
# #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("darkblue","lightblue"),lwd=4,cex=2)
# 
# plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="lightblue",lwd=3,type="l",xlab="l",ylab="D(\u2113)",main="(f) Optimal (Riemann)",cex.main=2,cex.axis=2,cex.lab=2)
# #grid(20,10,lwd=1)
# for (i in 1:(2*N-3)) 
# {polygon(c(0,8760/8760*d_int_alt[2*i],8760/8760*d_int_alt[2*i],0),c(LDC6.poly(d_int_alt[2*i+1])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*i+1])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*i-1])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*i-1])/LDC6.poly(d_int_alt[1])),col=col_bg[i])}
# {polygon(c(0,8760/8760*d_int_alt[2*N],8760/8760*d_int_alt[2*N],0),c(LDC6.poly(d_int_alt[2*N])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*N])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*N-1])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*N-1])/LDC6.poly(d_int_alt[1])),col=col_bg[i])}
# lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="darkblue")
# lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="lightblue",lwd=3)
# #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("darkblue","lightblue"),lwd=4,cex=2)
# mtext("Title for Two Plots",outer=TRUE,cex=1.5)
# dev.off()
# 
# ################################################################################################################################################################################
# 
# # # Determine the Desired Number of Load Blocks
# # 
# N=10
# 
# # === Approach 1: Horizontal Equidistant Discretisation ===
# 
# # = Upper Sums =
# 
# delta_h_up=1/N
# 
# # Error
# 
# integrand<-function(x) LDC6.poly(x)/LDC6.poly(0)
# err_int_up_he<-sum(sapply(1:N, function(i) delta_h_up*LDC6.poly((i-1)*delta_h_up)/LDC6.poly(0)))-integrate(Vectorize(integrand),lower=0,upper=1)$value
# print(err_int_up_he)
# 
# ell_h_up<-0
# f_ell_h_up<-LDC6.poly(0)/LDC6.poly(0)
# for (i in 1:N) {
#   ell_h_up<-c(ell_h_up,i*delta_h_up)
#   f_ell_h_up<-c(f_ell_h_up,LDC6.poly(i*delta_h_up)/LDC6.poly(0))
# }
# 
# # = Alternating Sums =
# 
# delta_h_alt=1/(2*N-1)
# 
# # Error
# 
# err_int_alt_he<-sum(sapply(1:(2*N-1), function(i) (-1)^(i+1)*delta_h_alt*LDC6.poly(2*ceiling((i-1)/2)*delta_h_alt)/LDC6.poly(0)+(-1)^i*integrate(Vectorize(integrand),lower=(i-1)*delta_h_alt,upper=i*delta_h_alt)$value))
# print(err_int_alt_he)
# 
# ell_h_alt<-0
# f_ell_h_alt<-LDC6.poly(0)/LDC6.poly(0)
# for (i in 1:(2*N-1)) {
#   ell_h_alt<-c(ell_h_alt,i*delta_h_alt)
#   f_ell_h_alt<-c(f_ell_h_alt,LDC6.poly(i*delta_h_alt)/LDC6.poly(0))
# }
# 
# # === Approach 2: Vertical Equidistant Approximation of D ===
# 
# # Define the Inverse of the Polynomial
# 
# LDC6.poly.inv<-inverse(integrand,lower=0,upper=1)
# 
# # = Upper Sums =
# 
# delta_v_up=(LDC6.poly(0)-LDC6.poly(1))/(N*LDC6.poly(0))
# 
# # Error
# 
# err_int_up_ve<-sum(sapply(1:N, function(i) (LDC6.poly.inv(integrand(1)+(N-i)*delta_v_up)-LDC6.poly.inv(integrand(1)+(N-i+1)*delta_v_up))*(integrand(1)+(N-i+1)*delta_v_up)-integrate(Vectorize(integrand),lower=LDC6.poly.inv(integrand(1)+(N-i+1)*delta_v_up),upper=LDC6.poly.inv(integrand(1)+(N-i)*delta_v_up))$value))
# print(err_int_up_ve)
# 
# ell_v_up<-0
# f_ell_v_up<-LDC6.poly(0)/LDC6.poly(0)
# for (i in 1:N) {
#   ell_v_up<-c(ell_v_up,LDC6.poly.inv(integrand(1)+(N-i)*delta_v_up))
#   f_ell_v_up<-c(f_ell_v_up,1-i*delta_v_up)
# }
# 
# # = Alternating Sums =
# 
# delta_v_alt=(LDC6.poly(0)-LDC6.poly(1))/((2*N-1)*LDC6.poly(0))
# 
# # Error
# 
# err_int_alt_ve<-sum(sapply(1:(2*N-1), function(i) (-1)^(i+1)*(LDC6.poly.inv(integrand(1)+(2*N-1-i)*delta_v_alt)-LDC6.poly.inv(integrand(1)+(2*N-i)*delta_v_alt))*(integrand(1)+(2*N-1-2*ceiling((i-1)/2))*delta_v_alt)+(-1)^i*integrate(Vectorize(integrand),lower=LDC6.poly.inv(integrand(1)+(2*N-i)*delta_v_alt),upper=LDC6.poly.inv(integrand(1)+(2*N-1-i)*delta_v_alt))$value))
# print(err_int_alt_ve)
# 
# ell_v_alt<-0
# f_ell_v_alt<-LDC6.poly(0)/LDC6.poly(0)
# for (i in 1:(2*N-1)) {
#   ell_v_alt<-c(ell_v_alt,LDC6.poly.inv(integrand(1)+(2*N-1-i)*delta_v_alt))
#   f_ell_v_alt<-c(f_ell_v_alt,1-i*delta_v_alt)
# }
# 
# # === Approach 3: Optimal Discretisation of the Load Duration Curve ===
# 
# # = Upper Sums =
# 
# # Compute the functions whose roots need to be found later via Newton-Raphson
# 
# grad_LDC_int<-function(x,j) sum(sapply(1:6,function(i) a[i+1]*(x[j-1]^i-x[j]^i+(x[j+1]-x[j])*i*x[j]^(i-1))))
# 
# grad_LDC_combined_int<-function(x) {
#   f<-NULL
#   for (j in 2:(N-2)) f<-c(f,grad_LDC_int(x,j))
#   f}
# 
# # Compute the Needed Roots via the Newton-Raphson Procedure
# 
# if (N==2){
#   model_int<-function(x) c(sum(sapply(1:6,function(i) a[i+1]*(-x[1]^i+(1-x[1])*i*x[1]^(i-1)))))
#   d_int<-c(0,uniroot(f=model_int,interval=c(0.1,0.9))$root,1)
#   LDC_int.d<-sapply(d_int,LDC6.poly)/LDC6.poly(0)
# }
# 
# if (N==3){
#   model_int<-function(x) c(sum(sapply(1:6,function(i) a[i+1]*(-x[1]^i+(x[2]-x[1])*i*x[1]^(i-1)))),
#                            sum(sapply(1:6,function(i) a[i+1]*(x[N-2]^i-x[N-1]^i+(1-x[N-1])*i*x[N-1]^(i-1)))))
#   d_int<-c(0,multiroot(f=model_int,start=seq(0.3,0.9,length.out=(N-1)),positive=TRUE)$root,1)
#   LDC_int.d<-sapply(d_int,LDC6.poly)/LDC6.poly(0)
# }
# 
# if (N>=4){
#   model_int<-function(x) c(D_max-1+sum(sapply(1:6,function(i) a[i+1]*(-x[1]^i+(x[2]-x[1])*i*x[1]^(i-1)))),
#                            grad_LDC_combined_int(x),
#                            sum(sapply(1:6,function(i) a[i+1]*(x[N-2]^i-x[N-1]^i+(1-x[N-1])*i*x[N-1]^(i-1)))))
#   d_int<-c(0,multiroot(f=model_int,start=seq(0.4,0.8,length.out=(N-1)),positive=TRUE)$root,1)
#   LDC_int.d<-sapply(d_int,LDC6.poly)/LDC6.poly(0)
# }
# 
# # Compute the Error Relative to the Polynomial Load Duration Curve
# 
# err_int_up_opt=sum(sapply(1:N, function(i) (d_int[i+1]-d_int[i])*LDC_int.d[i]))-integrate(Vectorize(integrand),lower=0,upper=1)$value
# print(err_int_up_opt)
# 
# ell_opt_up<-d_int
# f_ell_opt_up<-LDC_int.d
# 
# 
# # = Alternating Sums =
# 
# # Compute the functions whose roots need to be found later via Newton-Raphson
# 
# grad_LDC_int_alt<-function(x,j) c((x[j-1]-2*x[j]+x[j+1])*sum(sapply(1:6,function(i) a[i+1]*i*x[j]^(i-1))),
#                                   sum(sapply(1:6,function(i) a[i+1]*(x[j]^i-2*x[j+1]^i+x[j+2]^i))))
# 
# grad_LDC_combined_int_alt<-function(x) {
#   f<-NULL
#   for (j in 1:(N-2)) f<-c(f,grad_LDC_int_alt(x,2*(j-1)+2))
#   f}
# 
# # Compute the needed roots via the Newton-Raphson Procedure
# 
# 
# if (N==2){
#   model_int_alt<-function(x) c(sum(sapply(1:6,function(i) a[i+1]*(-2*x[1]^i+x[2]^i))),
#                                (x[1]-2*x[2]+1)*sum(sapply(1:6,function(i) a[i+1]*i*x[2]^(i-1))))
#   d_int_alt<-c(0,multiroot(f=model_int_alt,start=seq(0.1,0.9,length.out=2*(N-1)),positive=TRUE)$root,1)
#   LDC_int_alt.d<-sapply(d_int_alt,LDC6.poly)/LDC6.poly(0)
# }
# 
# if (N==3){
#   model_int_alt<-function(x) c(sum(sapply(1:6,function(i) a[i+1]*(-2*x[1]^i+x[2]^i))),
#                                (x[1]-2*x[2]+x[3])*sum(sapply(1:6,function(i) a[i+1]*i*x[2]^(i-1))),
#                                sum(sapply(1:6,function(i) a[i+1]*(x[2]^i-2*x[3]^i+x[4]^i))),
#                                (x[3]-2*x[4]+1)*sum(sapply(1:6,function(i) a[i+1]*i*x[4]^(i-1))))
#   d_int_alt<-c(0,multiroot(f=model_int_alt,start=seq(0.1,0.9,length.out=2*(N-1)),positive=TRUE)$root,1)
#   LDC_int_alt.d<-sapply(d_int_alt,LDC6.poly)/LDC6.poly(0)
# }
# 
# if (N>=4){
#   model_int_alt<-function(x) c(sum(sapply(1:6,function(i) a[i+1]*(-2*x[1]^i+x[2]^i))),
#                                grad_LDC_combined_int_alt(x),
#                                (x[2*(N-1)-1]-2*x[2*(N-1)]+1)*sum(sapply(1:6,function(i) a[i+1]*i*x[2*(N-1)]^(i-1))))
#   d_int_alt<-c(0,multiroot(f=model_int_alt,start=seq(0.1,0.9,length.out=2*(N-1)),positive=TRUE)$root,1)
#   LDC_int_alt.d<-sapply(d_int_alt,LDC6.poly)/LDC6.poly(0)
# }
# 
# # Error
# 
# err_int_alt_opt<-sum(sapply(1:(2*(N-1)), function(i) (-1)^(i+1)*(d_int_alt[i+1]-d_int_alt[i])*LDC6.poly(d_int_alt[1+2*ceiling((i-1)/2)])/LDC6.poly(0)+(-1)^i*integrate(Vectorize(integrand),lower=d_int_alt[i],upper=d_int_alt[i+1])$value))
# print(err_int_alt_opt)
# 
# ell_opt_alt<-d_int_alt
# f_ell_opt_alt<-LDC_int_alt.d
# 
# # = Errors of the Six Approaches =
# 
# error<-cbind(c(err_int_up_he,err_int_alt_he),c(err_int_up_ve,err_int_alt_ve),c(err_int_up_opt,err_int_alt_opt))
# colnames(error)<-c("HE","VE","OPT")
# rownames(error)<-c("US","AS")
# View(error)
# 
# # = Mesh Values of all Approaches =
# 
# # Upper Sums
# 
# mesh_up<-rbind(ell_h_up,ell_v_up,ell_opt_up)
# rownames(mesh_up)<-c("HE","VE","OPT")
# #print(mesh_up)
# 
# # Alternating Sums
# 
# mesh_alt<-rbind(ell_h_alt,ell_v_alt,ell_opt_alt)
# rownames(mesh_alt)<-c("HE","VE","OPT")
# #print(mesh_alt)
# 
# d<-d_int_alt
# LDC.d<-LDC_int_alt.d
# 
# 
# # = Function of Mesh Values of all Approaches =
# 
# # Upper Sums
# 
# f_mesh_up<-rbind(f_ell_h_up,f_ell_v_up,f_ell_opt_up)
# rownames(f_mesh_up)<-c("HE","VE","OPT")
# #print(f_mesh_up)
# 
# # Alternating Sums
# 
# f_mesh_alt<-rbind(f_ell_h_alt,f_ell_v_alt,f_ell_opt_alt)
# rownames(f_mesh_alt)<-c("HE","VE","OPT")
# #print(f_mesh_alt)
# View(rbind(mesh_up,f_mesh_up))
# View(rbind(mesh_alt,f_mesh_alt))
# 
# # = Overall Plot (of all the Approaches) =
# 
# # N=10
# # 
# # CairoPDF("results_N_10.pdf",width=14,height=8)
# # par(mar=c(5,5,5,5))
# # par(mfrow=c(3,2))
# # plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3,type="l",xlab="\u2113",ylab="D(\u2113)",main="(a) Equidistant in \u2113 (Darboux)",cex.main=2,cex.axis=2,cex.lab=2)
# # #grid(20,10,lwd=1)
# # for (i in 1:N) {polygon(c(0,8760/8760*i*delta_h_up,8760/8760*i*delta_h_up,0),c(LDC6.poly(i*delta_h_up)/LDC6.poly(0),LDC6.poly(i*delta_h_up)/LDC6.poly(0),LDC6.poly((i-1)*delta_h_up)/LDC6.poly(0),LDC6.poly((i-1)*delta_h_up)/LDC6.poly(0)),col=col_b50[1])}
# # lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="black")
# # lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3)
# # #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("black","yellow"),lwd=4,cex=2)
# # 
# # plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3,type="l",xlab="\u2113",ylab="D(\u2113)",main="(b) Equidistant in \u2113 (Riemann)",cex.main=2,cex.axis=2,cex.lab=2)
# # #grid(20,10,lwd=1)
# # {polygon(c(0,8760/8760*1*delta_h_alt,8760/8760*1*delta_h_alt,0),c(LDC6.poly(2*delta_h_alt)/max(LDC.emp),LDC6.poly(2*delta_h_alt)/LDC6.poly(0),LDC6.poly(0)/LDC6.poly(0),LDC6.poly(0)/LDC6.poly(0)),col=col_b50[1])}
# # for (i in 1:(N-2))
# # {polygon(c(0,8760/8760*(2*i+1)*delta_h_alt,8760/8760*(2*i+1)*delta_h_alt,0),c(LDC6.poly(2*(i+1)*delta_h_alt)/LDC6.poly(0),LDC6.poly(2*(i+1)*delta_h_alt)/LDC6.poly(0),LDC6.poly(2*i*delta_h_alt)/LDC6.poly(0),LDC6.poly(2*i*delta_h_alt)/LDC6.poly(0)),col=col_b50[1])}
# # {polygon(c(0,8760/8760*(2*N-1)*delta_h_alt,8760/8760*(2*N-1)*delta_h_alt,0),c(LDC6.poly((2*N-1)*delta_h_alt)/LDC6.poly(0),LDC6.poly((2*N-1)*delta_h_alt)/LDC6.poly(0),LDC6.poly(2*(N-1)*delta_h_alt)/LDC6.poly(0),LDC6.poly(2*(N-1)*delta_h_alt)/LDC6.poly(0)),col=col_b50[1])}
# # lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="black")
# # lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3)
# # #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("black","yellow"),lwd=4,cex=2)
# # 
# # plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3,type="l",xlab="\u2113",ylab="D(\u2113)",main="(c) Equidistant in D (Darboux)",cex.main=2,cex.axis=2,cex.lab=2)
# # #grid(20,10,lwd=1)
# # for (i in 1:N) {polygon(c(0,8760/8760*LDC6.poly.inv(integrand(1)+(N-i)*delta_v_up),8760/8760*LDC6.poly.inv(integrand(1)+(N-i)*delta_v_up),0),c(integrand(1)+(N-(i+1)+1)*delta_v_up,integrand(1)+(N-(i+1)+1)*delta_v_up,integrand(1)+(N-i+1)*delta_v_up,integrand(1)+(N-i+1)*delta_v_up),col=col_b50[1])}
# # lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="black")
# # lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3)
# # #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("black","yellow"),lwd=4,cex=2)
# # 
# # plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3,type="l",xlab="\u2113",ylab="D(\u2113)",main="(d) Equidistant in D (Riemann)",cex.main=2,cex.axis=2,cex.lab=2)
# # #grid(20,10,lwd=1)
# # {polygon(c(0,8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-1)*delta_v_alt),8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-1)*delta_v_alt),0),c(integrand(1)+(2*N-1-2*ceiling((2-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((2-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((1-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((1-1)/2))*delta_v_alt),col=col_b50[1])}
# # for (i in 1:(N-2))
# # {polygon(c(0,8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-2*(i+1)+1)*delta_v_alt),8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-2*(i+1)+1)*delta_v_alt),0),c(integrand(1)+(2*N-1-2*ceiling((2*(i+1)-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((2*(i+1)-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((2*i-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((2*i-1)/2))*delta_v_alt),col=col_b50[1])}
# # {polygon(c(0,8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-(2*N-1))*delta_v_alt),8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-(2*N-1))*delta_v_alt),0),c(integrand(1),integrand(1),integrand(1)+(2*N-1-2*ceiling((2*(N-1)-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((2*(N-1)-1)/2))*delta_v_alt),col=col_b50[1])}
# # lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="black")
# # lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3)
# # #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("black","yellow"),lwd=4,cex=2)
# # 
# # plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3,type="l",xlab="\u2113",ylab="D(\u2113)",main="(e) Optimal (Darboux)",cex.main=2,cex.axis=2,cex.lab=2)
# # #grid(20,10,lwd=1)
# # for (i in 1:N) {polygon(c(0,8760/8760*d_int[i+1],8760/8760*d_int[i+1],0),c(LDC6.poly(d_int[i+1])/LDC6.poly(d_int[1]),LDC6.poly(d_int[i+1])/LDC6.poly(d_int[1]),LDC6.poly(d_int[i])/LDC6.poly(d_int[1]),LDC6.poly(d_int[i])/LDC6.poly(d_int[1])),col=col_b50[1])}
# # lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="black")
# # lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3)
# # #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("black","yellow"),lwd=4,cex=2)
# # 
# # plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3,type="l",xlab="\u2113",ylab="D(\u2113)",main="(f) Optimal (Riemann)",cex.main=2,cex.axis=2,cex.lab=2)
# # #grid(20,10,lwd=1)
# # for (i in 1:(2*N-3))
# # {polygon(c(0,8760/8760*d_int_alt[2*i],8760/8760*d_int_alt[2*i],0),c(LDC6.poly(d_int_alt[2*i+1])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*i+1])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*i-1])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*i-1])/LDC6.poly(d_int_alt[1])),col=col_b50[1])}
# # {polygon(c(0,8760/8760*d_int_alt[2*N],8760/8760*d_int_alt[2*N],0),c(LDC6.poly(d_int_alt[2*N])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*N])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*N-1])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*N-1])/LDC6.poly(d_int_alt[1])),col=col_b50[1])}
# # lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="black")
# # lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3)
# # #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("black","yellow"),lwd=4,cex=2)
# # dev.off()
# # 
# # ################################################################################################################################################################
# # 
# # Determine the Desired Number of Load Blocks
# 
# N=50
# 
# # === Approach 1: Horizontal Equidistant Discretisation ===
# 
# # = Upper Sums =
# 
# delta_h_up=1/N
# 
# # Error
# 
# integrand<-function(x) LDC6.poly(x)/LDC6.poly(0)
# err_int_up_he<-sum(sapply(1:N, function(i) delta_h_up*LDC6.poly((i-1)*delta_h_up)/LDC6.poly(0)))-integrate(Vectorize(integrand),lower=0,upper=1)$value
# print(err_int_up_he)
# 
# ell_h_up<-0
# f_ell_h_up<-LDC6.poly(0)/LDC6.poly(0)
# for (i in 1:N) {
#   ell_h_up<-c(ell_h_up,i*delta_h_up)
#   f_ell_h_up<-c(f_ell_h_up,LDC6.poly(i*delta_h_up)/LDC6.poly(0))
# }
# 
# # = Alternating Sums =
# 
# delta_h_alt=1/(2*N-1)
# 
# # Error
# 
# err_int_alt_he<-sum(sapply(1:(2*N-1), function(i) (-1)^(i+1)*delta_h_alt*LDC6.poly(2*ceiling((i-1)/2)*delta_h_alt)/LDC6.poly(0)+(-1)^i*integrate(Vectorize(integrand),lower=(i-1)*delta_h_alt,upper=i*delta_h_alt)$value))
# print(err_int_alt_he)
# 
# ell_h_alt<-0
# f_ell_h_alt<-LDC6.poly(0)/LDC6.poly(0)
# for (i in 1:(2*N-1)) {
#   ell_h_alt<-c(ell_h_alt,i*delta_h_alt)
#   f_ell_h_alt<-c(f_ell_h_alt,LDC6.poly(i*delta_h_alt)/LDC6.poly(0))
# }
# 
# # === Approach 2: Vertical Equidistant Approximation of D ===
# 
# # Define the Inverse of the Polynomial
# 
# LDC6.poly.inv<-inverse(integrand,lower=0,upper=1)
# 
# # = Upper Sums =
# 
# delta_v_up=(LDC6.poly(0)-LDC6.poly(1))/(N*LDC6.poly(0))
# 
# # Error
# 
# err_int_up_ve<-sum(sapply(1:N, function(i) (LDC6.poly.inv(integrand(1)+(N-i)*delta_v_up)-LDC6.poly.inv(integrand(1)+(N-i+1)*delta_v_up))*(integrand(1)+(N-i+1)*delta_v_up)-integrate(Vectorize(integrand),lower=LDC6.poly.inv(integrand(1)+(N-i+1)*delta_v_up),upper=LDC6.poly.inv(integrand(1)+(N-i)*delta_v_up))$value))
# print(err_int_up_ve)
# 
# ell_v_up<-0
# f_ell_v_up<-LDC6.poly(0)/LDC6.poly(0)
# for (i in 1:N) {
#   ell_v_up<-c(ell_v_up,LDC6.poly.inv(integrand(1)+(N-i)*delta_v_up))
#   f_ell_v_up<-c(f_ell_v_up,1-i*delta_v_up)
# }
# 
# # = Alternating Sums =
# 
# delta_v_alt=(LDC6.poly(0)-LDC6.poly(1))/((2*N-1)*LDC6.poly(0))
# 
# # Error
# 
# err_int_alt_ve<-sum(sapply(1:(2*N-1), function(i) (-1)^(i+1)*(LDC6.poly.inv(integrand(1)+(2*N-1-i)*delta_v_alt)-LDC6.poly.inv(integrand(1)+(2*N-i)*delta_v_alt))*(integrand(1)+(2*N-1-2*ceiling((i-1)/2))*delta_v_alt)+(-1)^i*integrate(Vectorize(integrand),lower=LDC6.poly.inv(integrand(1)+(2*N-i)*delta_v_alt),upper=LDC6.poly.inv(integrand(1)+(2*N-1-i)*delta_v_alt))$value))
# print(err_int_alt_ve)
# 
# ell_v_alt<-0
# f_ell_v_alt<-LDC6.poly(0)/LDC6.poly(0)
# for (i in 1:(2*N-1)) {
#   ell_v_alt<-c(ell_v_alt,LDC6.poly.inv(integrand(1)+(2*N-1-i)*delta_v_alt))
#   f_ell_v_alt<-c(f_ell_v_alt,1-i*delta_v_alt)
# }
# 
# # === Approach 3: Optimal Discretisation of the Load Duration Curve ===
# 
# # = Upper Sums =
# 
# # Compute the functions whose roots need to be found later via Newton-Raphson
# 
# grad_LDC_int<-function(x,j) sum(sapply(1:6,function(i) a[i+1]*(x[j-1]^i-x[j]^i+(x[j+1]-x[j])*i*x[j]^(i-1))))
# 
# grad_LDC_combined_int<-function(x) {
#   f<-NULL
#   for (j in 2:(N-2)) f<-c(f,grad_LDC_int(x,j))
#   f}
# 
# # Compute the Needed Roots via the Newton-Raphson Procedure
# 
# if (N==2){
#   model_int<-function(x) c(sum(sapply(1:6,function(i) a[i+1]*(-x[1]^i+(1-x[1])*i*x[1]^(i-1)))))
#   d_int<-c(0,uniroot(f=model_int,interval=c(0.1,0.9))$root,1)
#   LDC_int.d<-sapply(d_int,LDC6.poly)/LDC6.poly(0)
# }
# 
# if (N==3){
#   model_int<-function(x) c(sum(sapply(1:6,function(i) a[i+1]*(-x[1]^i+(x[2]-x[1])*i*x[1]^(i-1)))),
#                            sum(sapply(1:6,function(i) a[i+1]*(x[N-2]^i-x[N-1]^i+(1-x[N-1])*i*x[N-1]^(i-1)))))
#   d_int<-c(0,multiroot(f=model_int,start=seq(0.3,0.9,length.out=(N-1)),positive=TRUE)$root,1)
#   LDC_int.d<-sapply(d_int,LDC6.poly)/LDC6.poly(0)
# }
# 
# if (N>=4){
#   model_int<-function(x) c(D_max-1+sum(sapply(1:6,function(i) a[i+1]*(-x[1]^i+(x[2]-x[1])*i*x[1]^(i-1)))),
#                            grad_LDC_combined_int(x),
#                            sum(sapply(1:6,function(i) a[i+1]*(x[N-2]^i-x[N-1]^i+(1-x[N-1])*i*x[N-1]^(i-1)))))
#   d_int<-c(0,multiroot(f=model_int,start=seq(0.4,0.8,length.out=(N-1)),positive=TRUE)$root,1)
#   LDC_int.d<-sapply(d_int,LDC6.poly)/LDC6.poly(0)
# }
# 
# # Compute the Error Relative to the Polynomial Load Duration Curve
# 
# err_int_up_opt=sum(sapply(1:N, function(i) (d_int[i+1]-d_int[i])*LDC_int.d[i]))-integrate(Vectorize(integrand),lower=0,upper=1)$value
# print(err_int_up_opt)
# 
# ell_opt_up<-d_int
# f_ell_opt_up<-LDC_int.d
# 
# 
# # = Alternating Sums =
# 
# # Compute the functions whose roots need to be found later via Newton-Raphson
# 
# grad_LDC_int_alt<-function(x,j) c((x[j-1]-2*x[j]+x[j+1])*sum(sapply(1:6,function(i) a[i+1]*i*x[j]^(i-1))),
#                                   sum(sapply(1:6,function(i) a[i+1]*(x[j]^i-2*x[j+1]^i+x[j+2]^i))))
# 
# grad_LDC_combined_int_alt<-function(x) {
#   f<-NULL
#   for (j in 1:(N-2)) f<-c(f,grad_LDC_int_alt(x,2*(j-1)+2))
#   f}
# 
# # Compute the needed roots via the Newton-Raphson Procedure
# 
# 
# if (N==2){
#   model_int_alt<-function(x) c(sum(sapply(1:6,function(i) a[i+1]*(-2*x[1]^i+x[2]^i))),
#                                (x[1]-2*x[2]+1)*sum(sapply(1:6,function(i) a[i+1]*i*x[2]^(i-1))))
#   d_int_alt<-c(0,multiroot(f=model_int_alt,start=seq(0.1,0.9,length.out=2*(N-1)),positive=TRUE)$root,1)
#   LDC_int_alt.d<-sapply(d_int_alt,LDC6.poly)/LDC6.poly(0)
# }
# 
# if (N==3){
#   model_int_alt<-function(x) c(sum(sapply(1:6,function(i) a[i+1]*(-2*x[1]^i+x[2]^i))),
#                                (x[1]-2*x[2]+x[3])*sum(sapply(1:6,function(i) a[i+1]*i*x[2]^(i-1))),
#                                sum(sapply(1:6,function(i) a[i+1]*(x[2]^i-2*x[3]^i+x[4]^i))),
#                                (x[3]-2*x[4]+1)*sum(sapply(1:6,function(i) a[i+1]*i*x[4]^(i-1))))
#   d_int_alt<-c(0,multiroot(f=model_int_alt,start=seq(0.1,0.9,length.out=2*(N-1)),positive=TRUE)$root,1)
#   LDC_int_alt.d<-sapply(d_int_alt,LDC6.poly)/LDC6.poly(0)
# }
# 
# if (N>=4){
#   model_int_alt<-function(x) c(sum(sapply(1:6,function(i) a[i+1]*(-2*x[1]^i+x[2]^i))),
#                                grad_LDC_combined_int_alt(x),
#                                (x[2*(N-1)-1]-2*x[2*(N-1)]+1)*sum(sapply(1:6,function(i) a[i+1]*i*x[2*(N-1)]^(i-1))))
#   d_int_alt<-c(0,multiroot(f=model_int_alt,start=seq(0.1,0.9,length.out=2*(N-1)),positive=TRUE)$root,1)
#   LDC_int_alt.d<-sapply(d_int_alt,LDC6.poly)/LDC6.poly(0)
# }
# 
# # Error
# 
# err_int_alt_opt<-sum(sapply(1:(2*(N-1)), function(i) (-1)^(i+1)*(d_int_alt[i+1]-d_int_alt[i])*LDC6.poly(d_int_alt[1+2*ceiling((i-1)/2)])/LDC6.poly(0)+(-1)^i*integrate(Vectorize(integrand),lower=d_int_alt[i],upper=d_int_alt[i+1])$value))
# print(err_int_alt_opt)
# 
# ell_opt_alt<-d_int_alt
# f_ell_opt_alt<-LDC_int_alt.d
# 
# # = Errors of the Six Approaches =
# 
# error<-cbind(c(err_int_up_he,err_int_alt_he),c(err_int_up_ve,err_int_alt_ve),c(err_int_up_opt,err_int_alt_opt))
# colnames(error)<-c("HE","VE","OPT")
# rownames(error)<-c("US","AS")
# View(error)
# 
# # = Mesh Values of all Approaches =
# 
# # Upper Sums
# 
# mesh_up<-rbind(ell_h_up,ell_v_up,ell_opt_up)
# rownames(mesh_up)<-c("HE","VE","OPT")
# #print(mesh_up)
# 
# # Alternating Sums
# 
# mesh_alt<-rbind(ell_h_alt,ell_v_alt,ell_opt_alt)
# rownames(mesh_alt)<-c("HE","VE","OPT")
# #print(mesh_alt)
# 
# d<-d_int_alt
# LDC.d<-LDC_int_alt.d
# 
# 
# # = Function of Mesh Values of all Approaches =
# 
# # Upper Sums
# 
# f_mesh_up<-rbind(f_ell_h_up,f_ell_v_up,f_ell_opt_up)
# rownames(f_mesh_up)<-c("HE","VE","OPT")
# #print(f_mesh_up)
# 
# # Alternating Sums
# 
# f_mesh_alt<-rbind(f_ell_h_alt,f_ell_v_alt,f_ell_opt_alt)
# rownames(f_mesh_alt)<-c("HE","VE","OPT")
# #print(f_mesh_alt)
# View(rbind(mesh_up,f_mesh_up))
# View(rbind(mesh_alt,f_mesh_alt))
# 
# # # = Overall Plot (of all the Approaches) =
# # 
# # N=50
# 
# col_b50[1]<-col_bg[25]
# 
# # CairoPDF("results_N_50.pdf",width=14,height=8)
# # par(mar=c(5,5,5,5))
# # par(mfrow=c(3,2))
# # plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3,type="l",xlab="\u2113",ylab="D(\u2113)",main="(a) Equidistant in \u2113 (Darboux)",cex.main=2,cex.axis=2,cex.lab=2)
# # #grid(20,10,lwd=1)
# # for (i in 1:N) {polygon(c(0,8760/8760*i*delta_h_up,8760/8760*i*delta_h_up,0),c(LDC6.poly(i*delta_h_up)/LDC6.poly(0),LDC6.poly(i*delta_h_up)/LDC6.poly(0),LDC6.poly((i-1)*delta_h_up)/LDC6.poly(0),LDC6.poly((i-1)*delta_h_up)/LDC6.poly(0)),col=col_b50[1])}
# # lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="black")
# # lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3)
# # #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("black","yellow"),lwd=4,cex=2)
# # 
# # plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3,type="l",xlab="\u2113",ylab="D(\u2113)",main="(b) Equidistant in \u2113 (Riemann)",cex.main=2,cex.axis=2,cex.lab=2)
# # #grid(20,10,lwd=1)
# # {polygon(c(0,8760/8760*1*delta_h_alt,8760/8760*1*delta_h_alt,0),c(LDC6.poly(2*delta_h_alt)/max(LDC.emp),LDC6.poly(2*delta_h_alt)/LDC6.poly(0),LDC6.poly(0)/LDC6.poly(0),LDC6.poly(0)/LDC6.poly(0)),col=col_b50[1])}
# # for (i in 1:(N-2))
# # {polygon(c(0,8760/8760*(2*i+1)*delta_h_alt,8760/8760*(2*i+1)*delta_h_alt,0),c(LDC6.poly(2*(i+1)*delta_h_alt)/LDC6.poly(0),LDC6.poly(2*(i+1)*delta_h_alt)/LDC6.poly(0),LDC6.poly(2*i*delta_h_alt)/LDC6.poly(0),LDC6.poly(2*i*delta_h_alt)/LDC6.poly(0)),col=col_b50[1])}
# # {polygon(c(0,8760/8760*(2*N-1)*delta_h_alt,8760/8760*(2*N-1)*delta_h_alt,0),c(LDC6.poly((2*N-1)*delta_h_alt)/LDC6.poly(0),LDC6.poly((2*N-1)*delta_h_alt)/LDC6.poly(0),LDC6.poly(2*(N-1)*delta_h_alt)/LDC6.poly(0),LDC6.poly(2*(N-1)*delta_h_alt)/LDC6.poly(0)),col=col_b50[1])}
# # lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="black")
# # lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3)
# # #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("black","yellow"),lwd=4,cex=2)
# # 
# # plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3,type="l",xlab="\u2113",ylab="D(\u2113)",main="(c) Equidistant in D (Darboux)",cex.main=2,cex.axis=2,cex.lab=2)
# # #grid(20,10,lwd=1)
# # for (i in 1:N) {polygon(c(0,8760/8760*LDC6.poly.inv(integrand(1)+(N-i)*delta_v_up),8760/8760*LDC6.poly.inv(integrand(1)+(N-i)*delta_v_up),0),c(integrand(1)+(N-(i+1)+1)*delta_v_up,integrand(1)+(N-(i+1)+1)*delta_v_up,integrand(1)+(N-i+1)*delta_v_up,integrand(1)+(N-i+1)*delta_v_up),col=col_b50[1])}
# # lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="black")
# # lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3)
# # #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("black","yellow"),lwd=4,cex=2)
# # 
# # plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3,type="l",xlab="\u2113",ylab="D(\u2113)",main="(d) Equidistant in D (Riemann)",cex.main=2,cex.axis=2,cex.lab=2)
# # #grid(20,10,lwd=1)
# # {polygon(c(0,8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-1)*delta_v_alt),8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-1)*delta_v_alt),0),c(integrand(1)+(2*N-1-2*ceiling((2-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((2-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((1-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((1-1)/2))*delta_v_alt),col=col_b50[1])}
# # for (i in 1:(N-2))
# # {polygon(c(0,8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-2*(i+1)+1)*delta_v_alt),8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-2*(i+1)+1)*delta_v_alt),0),c(integrand(1)+(2*N-1-2*ceiling((2*(i+1)-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((2*(i+1)-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((2*i-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((2*i-1)/2))*delta_v_alt),col=col_b50[1])}
# # {polygon(c(0,8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-(2*N-1))*delta_v_alt),8760/8760*LDC6.poly.inv(integrand(1)+(2*N-1-(2*N-1))*delta_v_alt),0),c(integrand(1),integrand(1),integrand(1)+(2*N-1-2*ceiling((2*(N-1)-1)/2))*delta_v_alt,integrand(1)+(2*N-1-2*ceiling((2*(N-1)-1)/2))*delta_v_alt),col=col_b50[1])}
# # lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="black")
# # lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3)
# # #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("black","yellow"),lwd=4,cex=2)
# # 
# # plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3,type="l",xlab="\u2113",ylab="D(\u2113)",main="(e) Optimal (Darboux)",cex.main=2,cex.axis=2,cex.lab=2)
# # #grid(20,10,lwd=1)
# # for (i in 1:N) {polygon(c(0,8760/8760*d_int[i+1],8760/8760*d_int[i+1],0),c(LDC6.poly(d_int[i+1])/LDC6.poly(d_int[1]),LDC6.poly(d_int[i+1])/LDC6.poly(d_int[1]),LDC6.poly(d_int[i])/LDC6.poly(d_int[1]),LDC6.poly(d_int[i])/LDC6.poly(d_int[1])),col=col_b50[1])}
# # lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="black")
# # lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3)
# # #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("black","yellow"),lwd=4,cex=2)
# # 
# # plot(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3,type="l",xlab="\u2113",ylab="D(\u2113)",main="(f) Optimal (Riemann)",cex.main=2,cex.axis=2,cex.lab=2)
# # #grid(20,10,lwd=1)
# # for (i in 1:(2*N-3))
# # {polygon(c(0,8760/8760*d_int_alt[2*i],8760/8760*d_int_alt[2*i],0),c(LDC6.poly(d_int_alt[2*i+1])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*i+1])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*i-1])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*i-1])/LDC6.poly(d_int_alt[1])),col=col_b50[1])}
# # {polygon(c(0,8760/8760*d_int_alt[2*N],8760/8760*d_int_alt[2*N],0),c(LDC6.poly(d_int_alt[2*N])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*N])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*N-1])/LDC6.poly(d_int_alt[1]),LDC6.poly(d_int_alt[2*N-1])/LDC6.poly(d_int_alt[1])),col=col_b50[1])}
# # lines(seq(0,1,0.0001141553),LDC.norm,lwd=3,col="black")
# # lines(seq(0,1,0.01)*8760/8760,sapply(seq(0,1,0.01),LDC6.poly)/sapply(seq(0,1,0.01),LDC6.poly)[1],col="yellow",lwd=3)
# # #legend("topright",c("Empirical LDC","Polynomial LDC"),col=c("black","yellow"),lwd=4,cex=2)
# # dev.off()
# # 
# # ############################################################################################################################################
# # 
# N_max=50
# E<-matrix(0,N_max-1,6)
# colnames(E)<-c("HE-US","HE-AS","VE-US","VE-AS","OPT-US","OPT-AS")
# rownames(E)<-seq(2,N_max,1)
# for (j in 2:N_max) {
#   delta_h_up<-1/j
#   delta_h_alt<-1/(2*j-1)
#   delta_v_up<-(LDC6.poly(0)-LDC6.poly(1))/(j*LDC6.poly(0))
#   delta_v_alt<-(LDC6.poly(0)-LDC6.poly(1))/((2*j-1)*LDC6.poly(0))
#   if (j==2){
#     model_int<-function(x) c(D_max-1+sum(sapply(1:6,function(i) a[i+1]*(-x[1]^i+(1-x[1])*i*x[1]^(i-1)))))
#     d_int<-c(0,uniroot(f=model_int,interval=c(0.1,0.9))$root,1)
#     LDC_int.d<-sapply(d_int,LDC6.poly)/LDC6.poly(0)
#     model_int_alt<-function(x) c(D_max-1+sum(sapply(1:6,function(i) a[i+1]*(-2*x[1]^i+x[2]^i))),
#                                  (x[1]-2*x[2]+1)*sum(sapply(1:6,function(i) a[i+1]*i*x[2]^(i-1))))
#     d_int_alt<-c(0,multiroot(f=model_int_alt,start=seq(0.1,0.9,length.out=2*(j-1)),positive=TRUE)$root,1)
#     LDC_int_alt.d<-sapply(d_int_alt,LDC6.poly)
#   }
#   if (j==3){
#     model_int<-function(x) c(D_max-1+sum(sapply(1:6,function(i) a[i+1]*(-x[1]^i+(x[2]-x[1])*i*x[1]^(i-1)))),
#                              sum(sapply(1:6,function(i) a[i+1]*(x[j-2]^i-x[j-1]^i+(1-x[j-1])*i*x[j-1]^(i-1)))))
#     d_int<-c(0,multiroot(f=model_int,start=seq(0.3,0.9,length.out=(j-1)),positive=TRUE)$root,1)
#     LDC_int.d<-sapply(d_int,LDC6.poly)/LDC6.poly(0)
#     model_int_alt<-function(x) c(D_max-1+sum(sapply(1:6,function(i) a[i+1]*(-2*x[1]^i+x[2]^i))),
#                                  (x[1]-2*x[2]+x[3])*sum(sapply(1:6,function(i) a[i+1]*i*x[2]^(i-1))),
#                                  sum(sapply(1:6,function(i) a[i+1]*(x[2]^i-2*x[3]^i+x[4]^i))),
#                                  (x[3]-2*x[4]+1)*sum(sapply(1:6,function(i) a[i+1]*i*x[4]^(i-1))))
#     d_int_alt<-c(0,multiroot(f=model_int_alt,start=seq(0.3,0.9,length.out=2*(j-1)),positive=TRUE)$root,1)
#     LDC_int_alt.d<-sapply(d_int_alt,LDC6.poly)
#   }
#   if (j>=4){
#     grad_LDC_combined_int<-function(x) {
#       f<-NULL
#       for (k in 2:(j-2)) f<-c(f,grad_LDC_int(x,k))
#       f}
#     model_int<-function(x) c(D_max-1+sum(sapply(1:6,function(i) a[i+1]*(-x[1]^i+(x[2]-x[1])*i*x[1]^(i-1)))),
#                              grad_LDC_combined_int(x),
#                              sum(sapply(1:6,function(i) a[i+1]*(x[j-2]^i-x[j-1]^i+(1-x[j-1])*i*x[j-1]^(i-1)))))
#     d_int<-c(0,multiroot(f=model_int,start=seq(0.4,0.8,length.out=(j-1)),positive=TRUE)$root,1)
#     LDC_int.d<-sapply(d_int,LDC6.poly)/LDC6.poly(0)
#     grad_LDC_combined_int_alt<-function(x) {
#       f<-NULL
#       for (k in 1:(j-2)) f<-c(f,grad_LDC_int_alt(x,2*(k-1)+2))
#       f}
#     model_int_alt<-function(x) c(D_max-1+sum(sapply(1:6,function(i) a[i+1]*(-2*x[1]^i+x[2]^i))),
#                                  grad_LDC_combined_int_alt(x),
#                                  (x[2*(j-1)-1]-2*x[2*(j-1)]+1)*sum(sapply(1:6,function(i) a[i+1]*i*x[2*(j-1)]^(i-1))))
#     d_int_alt<-c(0,multiroot(f=model_int_alt,start=seq(0.2,0.9,length.out=2*(j-1)),positive=TRUE)$root,1)
#     LDC_int_alt.d<-sapply(d_int_alt,LDC6.poly)
#   }
#   E[j-1,]<-c(sum(sapply(1:j, function(i) delta_h_up*LDC6.poly((i-1)*delta_h_up)/LDC6.poly(0)))-integrate(Vectorize(integrand),lower=0,upper=1)$value,sum(sapply(1:(2*j-1), function(i) (-1)^(i+1)*delta_h_alt*LDC6.poly(2*ceiling((i-1)/2)*delta_h_alt)/LDC6.poly(0)+(-1)^i*integrate(Vectorize(integrand),lower=(i-1)*delta_h_alt,upper=i*delta_h_alt)$value)),sum(sapply(1:j, function(i) (LDC6.poly.inv(integrand(1)+(j-i)*delta_v_up)-LDC6.poly.inv(integrand(1)+(j-i+1)*delta_v_up))*(integrand(1)+(j-i+1)*delta_v_up)-integrate(Vectorize(integrand),lower=LDC6.poly.inv(integrand(1)+(j-i+1)*delta_v_up),upper=LDC6.poly.inv(integrand(1)+(j-i)*delta_v_up))$value)),sum(sapply(1:(2*j-1), function(i) (-1)^(i+1)*(LDC6.poly.inv(integrand(1)+(2*j-1-i)*delta_v_alt)-LDC6.poly.inv(integrand(1)+(2*j-i)*delta_v_alt))*(integrand(1)+(2*j-1-2*ceiling((i-1)/2))*delta_v_alt)+(-1)^i*integrate(Vectorize(integrand),lower=LDC6.poly.inv(integrand(1)+(2*j-i)*delta_v_alt),upper=LDC6.poly.inv(integrand(1)+(2*j-1-i)*delta_v_alt))$value)),sum(sapply(1:j, function(i) (d_int[i+1]-d_int[i])*LDC_int.d[i]))-integrate(Vectorize(integrand),lower=0,upper=1)$value,sum(sapply(1:(2*(j-1)), function(i) (-1)^(i+1)*(d_int_alt[i+1]-d_int_alt[i])*LDC6.poly(d_int_alt[1+2*ceiling((i-1)/2)])/LDC6.poly(0)+(-1)^i*integrate(Vectorize(integrand),lower=d_int_alt[i],upper=d_int_alt[i+1])$value)))
# }
# View(E)
# 
# # # R Plot
# # 
# # par(mfrow=c(1,2))
# # plot(E[,1],type="l",col="darkblue",lwd=2.5,xlab="N",ylab="Error",main="Darboux Sums")
# # lines(E[,3],type="l",col="blue",lwd=2.5)
# # lines(E[,5],type="l",col="lightblue",lwd=2.5)
# # grid(20,10,lwd=1)
# # legend("topright",c("Equidistant in l","Equidistant in D","Optimal"),col=c("darkblue","blue","lightblue"),lwd=2)
# # plot(E[,2],type="l",col="darkblue",lwd=2.5,xlab="N",ylab="Error",main="Riemann Sums")
# # lines(E[,4],type="l",col="blue",lwd=2.5)
# # lines(E[,6],type="l",col="lightblue",lwd=2.5)
# # grid(20,10,lwd=1)
# # legend("topright",c("Equidistant in l","Equidistant in D","Optimal"),col=c("darkblue","blue","lightblue"),lwd=2)
# # 
# # PDF Plot
# 
# # CairoPDF("results_errors_N.pdf",width=14,height=8)
# # par(mar=c(5,5,5,5))
# # par(mfrow=c(1,2))
# # plot(E[,1],type="l",col="darkblue",lwd=2.5,xlab="Number of Load Blocks",ylab="Error",main="Darboux Sums",cex.main=2,cex.axis=2,cex.lab=2)
# # grid(20,10,lwd=1.5)
# # lines(E[,3],type="l",col="blue",lwd=2.5)
# # lines(E[,5],type="l",col="lightblue",lwd=2.5)
# # legend("topright",c("Equidistant in \u2113","Equidistant in D","Optimal"),col=c("darkblue","blue","lightblue"),lwd=2)
# # 
# # plot(E[,2],type="l",col="darkblue",lwd=2.5,xlab="Number of Load Blocks",ylab="Error",main="Riemann Sums",cex.main=2,cex.axis=2,cex.lab=2)
# # grid(20,10,lwd=1.5)
# # lines(E[,4],type="l",col="blue",lwd=2.5)
# # lines(E[,6],type="l",col="lightblue",lwd=2.5)
# # legend("topright",c("Equidistant in \u2113","Equidistant in D","Optimal"),col=c("darkblue","blue","lightblue"),lwd=2)
# # dev.off()
# # 
# # 
# # par(mfrow=c(1,1))
# # plot(E[,1],type="l",col="darkblue",lwd=2.5,xlab="Number of Load Blocks",ylab="Error",main="Error Assuming Equidistant l")
# # lines(E[,2],type="l",col="blue",lwd=2.5)
# # grid(20,10,lwd=1)
# # legend("topright",c("Upper Sums","Alternating Sums"),col=c("darkblue","lightblue"),lwd=2)
# # 
# # par(mfrow=c(1,1))
# # plot(E[,3],type="l",col="darkgreen",lwd=2.5,xlab="Number of Load Blocks",ylab="Error",main="Error Assuming Equidistant D")
# # lines(E[,4],type="l",col="green",lwd=2.5)
# # grid(20,10,lwd=1)
# # legend("topright",c("Upper Sums","Alternating Sums"),col=c("darkgreen","green"),lwd=2)
# # 
# # par(mfrow=c(1,1))
# # plot(E[,5],type="l",col="darkred",lwd=2.5,xlab="Number of Load Blocks",ylab="Error",main="Error with Optimal l")
# # lines(E[,6],type="l",col="red",lwd=2.5)
# # grid(20,10,lwd=1)
# # legend("topright",c("Upper Sums","Alternating Sums"),col=c("darkred","red"),lwd=2)
# # 
# # # R
# # 
# # par(mfrow=c(1,1))
# # plot(E[1:9,1],type="l",col=col_b50[10],lwd=2.5,xlab="N",ylab="Error",main="All Six Approaches",ylim=c(0.01,0.14))
# # lines(E[1:9,2],type="l",col=col_b50[18],lwd=2.5)
# # lines(E[1:9,3],type="l",col=col_b50[25],lwd=2.5)
# # lines(E[1:9,4],type="l",col=col_b50[30],lwd=2.5)
# # lines(E[1:9,5],type="l",col=col_b50[40],lwd=2.5)
# # lines(E[1:9,6],type="l",col=col_b50[50],lwd=2.5)
# # grid(20,10,lwd=1)
# # legend("topright",c("Equidistant in l (D)","Equidistant in l (R)","Equidistant in D (D)","Equidistant in D (R)","Optimal (D)","Optimal (R)"),col=c(col_b50[10],col_b50[18],col_b50[25],col_b50[30],col_b50[40],col_b50[50]),lwd=2)
# # 
# # # PDF
# # 
# # CairoPDF("comparison_approaches.pdf",width=14,height=8)
# # par(mar=c(5,5,5,5))
# # par(mfrow=c(1,1))
# # grid(20,10,lwd=1.5)
# # plot(E[1:9,1],type="l",col=col_b50[10],lwd=2.5,xlab="Number of Load Blocks",ylab="Error",main="Comparison of All Approaches",cex.main=2,cex.axis=2,cex.lab=2,ylim=c(0.01,0.14))
# # lines(E[1:9,2],type="l",col=col_b50[18],lwd=2.5)
# # lines(E[1:9,3],type="l",col=col_b50[25],lwd=2.5)
# # lines(E[1:9,4],type="l",col=col_b50[30],lwd=2.5)
# # lines(E[1:9,5],type="l",col=col_b50[40],lwd=2.5)
# # lines(E[1:9,6],type="l",col=col_b50[50],lwd=2.5)
# # legend("topright",c("Equidistant in \u2113 (D)","Equidistant in \u2113 (R)","Equidistant in D (D)","Equidistant in D (R)","Optimal (D)","Optimal (R)"),col=c(col_b50[10],col_b50[18],col_b50[25],col_b50[30],col_b50[40],col_b50[50]),lwd=4,cex=2)
# # dev.off()
# # # 
# # Compare the Errors Using Upper and Alternating Sums for the Same Number of Variables
# 
# E_U<-matrix(0,2*N_max-1-1,3)
# colnames(E_U)<-c("HE","VE","OPT")
# rownames(E_U)<-seq(2,2*N_max-1,1)
# # View(E_U)
# 
# 
# for (j in 2:(2*N_max-1)) {
#   delta_h_up<-1/j
#   delta_v_up<-(LDC6.poly(0)-LDC6.poly(1))/(j*LDC6.poly(0))
#   if (j==2){
#     model_int<-function(x) c(D_max-1+sum(sapply(1:6,function(i) a[i+1]*(-x[1]^i+(1-x[1])*i*x[1]^(i-1)))))
#     d_int<-c(0,uniroot(f=model_int,interval=c(0.1,0.9))$root,1)
#     LDC_int.d<-sapply(d_int,LDC6.poly)/LDC6.poly(0)
#   }
#   if (j==3){
#     model_int<-function(x) c(D_max-1+sum(sapply(1:6,function(i) a[i+1]*(-x[1]^i+(x[2]-x[1])*i*x[1]^(i-1)))),
#                              sum(sapply(1:6,function(i) a[i+1]*(x[j-2]^i-x[j-1]^i+(1-x[j-1])*i*x[j-1]^(i-1)))))
#     d_int<-c(0,multiroot(f=model_int,start=seq(0.3,0.9,length.out=(j-1)),positive=TRUE)$root,1)
#     LDC_int.d<-sapply(d_int,LDC6.poly)/LDC6.poly(0)
#   }
#   if (j>=4){
#     grad_LDC_combined_int<-function(x) {
#       f<-NULL
#       for (k in 2:(j-2)) f<-c(f,grad_LDC_int(x,k))
#       f}
#     model_int<-function(x) c(D_max-1+sum(sapply(1:6,function(i) a[i+1]*(-x[1]^i+(x[2]-x[1])*i*x[1]^(i-1)))),
#                              grad_LDC_combined_int(x),
#                              sum(sapply(1:6,function(i) a[i+1]*(x[j-2]^i-x[j-1]^i+(1-x[j-1])*i*x[j-1]^(i-1)))))
#     d_int<-c(0,multiroot(f=model_int,start=seq(0.4,0.8,length.out=(j-1)),positive=TRUE)$root,1)
#     LDC_int.d<-sapply(d_int,LDC6.poly)/LDC6.poly(0)
#   }
#   E_U[j-1,]<-c(sum(sapply(1:j, function(i) delta_h_up*LDC6.poly((i-1)*delta_h_up)/LDC6.poly(0)))-integrate(Vectorize(integrand),lower=0,upper=1)$value,sum(sapply(1:j, function(i) (LDC6.poly.inv(integrand(1)+(j-i)*delta_v_up)-LDC6.poly.inv(integrand(1)+(j-i+1)*delta_v_up))*(integrand(1)+(j-i+1)*delta_v_up)-integrate(Vectorize(integrand),lower=LDC6.poly.inv(integrand(1)+(j-i+1)*delta_v_up),upper=LDC6.poly.inv(integrand(1)+(j-i)*delta_v_up))$value)),sum(sapply(1:j, function(i) (d_int[i+1]-d_int[i])*LDC_int.d[i]))-integrate(Vectorize(integrand),lower=0,upper=1)$value)
# }
# View(E_U)
# 
# # # Plots of the Errors Using Upper and Alternating Sums for the Same Number of Variables
# # 
# # # Equidistant in l
# # 
# # par(mfrow=c(1,1))
# # plot(E[,2],type="l",col="darkblue",lwd=2.5,xlab="Number of Variables",ylab="Error",main="Error for the Same Number of Variables (Equidistant l)")
# # lines(E_U[seq(2,2*N_max-1,2),1],type="l",col="blue",lwd=2.5)
# # legend("topright",c("Alternating Sums","Upper Sums"),col=c("darkblue","blue"),lwd=2)
# # 
# # # Equidistant in D
# # 
# # plot(E[,4],type="l",col="darkblue",lwd=2.5,xlab="Number of Variables",ylab="Error",main="Error for the Same Number of Variables (Equidistant D)")
# # lines(E_U[seq(2,2*N_max-1,2),2],type="l",col="blue",lwd=2.5)
# # legend("topright",c("Alternating Sums","Upper Sums"),col=c("darkblue","blue"),lwd=2)
# # 
# # # Optimal
# # 
# # plot(E[,6],type="l",col="darkblue",lwd=2.5,xlab="Number of Variables",ylab="Error",main="Error for the Same Number of Variables (Optimal)")
# # lines(E_U[seq(2,2*N_max-1,2),3],type="l",col="blue",lwd=2.5)
# # legend("topright",c("Alternating Sums","Upper Sums"),col=c("darkblue","blue"),lwd=2)
# # grid(20,10,lwd=1)
# # 
# # # Plots of the Difference Between the Errors Using Upper and Alternating Sums for the Same Number of Variables
# # 
# # # R
# # 
# # par(mfrow=c(3,1))
# # plot(E_U[seq(2,2*N_max-1,2),1]-E[,2],type="l",col="darkblue",lwd=2.5,xlab="Number of Variables",ylab="Error Difference",main="Equidistant in l")
# # grid(20,10,lwd=1.5)
# # plot(E_U[seq(2,2*N_max-1,2),2]-E[,4],type="l",col="darkblue",lwd=2.5,xlab="Number of Variables",ylab="Error Difference",main="Equidistant in D")
# # grid(20,10,lwd=1.5)
# # plot(E_U[seq(2,2*N_max-1,2),3]-E[,6],type="l",col="darkblue",lwd=2.5,xlab="Number of Variables",ylab="Error Difference",main="Optimal")
# # grid(20,10,lwd=1.5)
# # # 
# # # # PDF
# # # 
# # CairoPDF("results_errors_variables.pdf",width=14,height=8)
# # par(mar=c(5,5,5,5))
# # par(mfrow=c(3,1))
# # plot(E_U[seq(2,2*N_max-1,2),1]-E[,2],type="l",col="darkblue",lwd=2.5,xlab="Number of Variables",ylab="Error Difference",main="Equidistant in \u2113",cex.main=2.5,cex.axis=2,cex.lab=2)
# # grid(20,10,lwd=2)
# # plot(E_U[seq(2,2*N_max-1,2),2]-E[,4],type="l",col="darkblue",lwd=2.5,xlab="Number of Variables",ylab="Error Difference",main="Equidistant in D",cex.main=2.5,cex.axis=2,cex.lab=2)
# # grid(20,10,lwd=2)
# # plot(E_U[seq(2,2*N_max-1,2),3]-E[,6],type="l",col="darkblue",lwd=2.5,xlab="Number of Variables",ylab="Error Difference",main="Optimal",cex.main=2.5,cex.axis=2,cex.lab=2)
# # grid(20,10,lwd=2)
# # dev.off()
