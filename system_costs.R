# ===== System Costs =====

setwd("C:/Users/Mirescu1/ucloud/Universities/TU Wien/Theses/Dissertation/Eigene/Latex")

# === Install the necessary packages ===

library(xts)
library(quadprog)
library(rootSolve)
library(MASS)
library("colorspace")
library(grDevices)

# === Define colors to be used throughout this document ===

# Call up on interface where one can define a personalized color palette with subsequent color generator in HEX-format

#choose_palette()

# Introduce the personalized palettes into R

# Nikolaus Rab (nr) color palette

col_nr<-rev(c("#E5FCC2","#9DE0AD","#45ADA8","#547980","#594F4F"))

# Blue-Green (bg) color palette

col_bg<-rev(c("#00584D","#00685E","#347B73","#65948E","#A1BAB6","#A3B9BC","#689398","#387A81","#00666E","#00565F"))

# Pink-Purple (pp) color palette 

col_pp<-rev(c("#CA2497","#D456A6","#DF7CB7","#EAA4CC","#F6D5E7","#F0D6EF","#DEA7DD","#D181CF","#C65BC3","#BD2ABA"))

# Coray-Pink (cp) color palette

col_cp<-rev(c("#FF7C89","#FF97A0","#FFB0B6","#FCC8CB","#F2E0E1","#F2DFE4","#FDC7D6","#FFADC8","#FF91BA","#FF73AD"))

# Blue 50 [MAXIMUM] (b50) color palette

col_b50<-rev(c("#023FA5","#1041A5","#1A44A5","#2146A5","#2749A5","#2C4BA6","#314EA6","#3550A7","#3953A8","#3D56A9","#4158AA","#455BAB","#485DAC","#4C60AD",
               "#4F63AE","#5365AF","#5668B0","#5A6BB2","#5D6EB3","#6170B5","#6473B6","#6876B8","#6B79B9","#6E7CBB","#727FBC","#7582BE","#7985C0","#7D88C1",
               "#808BC3","#848FC5","#8892C7","#8B95C9","#8F99CB","#939CCD","#97A0CF","#9BA3D1","#9FA7D3","#A3ABD5","#A8AFD7","#ACB3DA","#B1B7DC","#B6BCDF",
               "#BBC0E1","#C0C5E4","#C5CAE7","#CBD0EA","#D2D5ED","#D9DCF1","#E1E4F5","#EEF0FC"))

# Define the colors used in the plot to be of the preferred palette

colors<-col_b50

# === Read in the Necessary Data ===

balancing<-read.csv("balancing.csv",header=T,sep=",",dec=".")
#capacity<-read.csv("capacity.csv",header=T,sep=",",dec=".")

# === Linear and Non-Linear Estimation ===

# = Balancing Costs =

# Linear Model Without Intercept

slope_BC<-as.vector(lm(balancing[,2]~balancing[,1]-1)$coef)*100

# Linear Model With Intercept

lmi_BC<-as.vector(lm(balancing[,2]~balancing[,1])$coef)

# Quadratic Polynomial Without Intercept

BC_1<-balancing[,1]
BC_2<-balancing[,2]

bal_sq<-nls(BC_2~0+b*BC_1+c*BC_1^2,start=list(b=1,c=1))
sq_BC<-as.vector(summary(bal_sq)$coef[,1])

# Quadratic Polynomial With Intercept

bali_sq<-nls(BC_2~a+b*BC_1+c*BC_1^2,start=list(a=1,b=1,c=1))
sqi_BC<-summary(bali_sq)$coef[,1]

# = Capacity / Reliability Costs =

# Linear Model Without Intercept

#slope_RC<-mean(capacity[,2])*1e-2*100

# Quadratic Polynomial Without Intercept

#sq_RC<-lm(capacity[,2]~capacity[,1]-1)$coef/2

# Plot

# Plot the 4 Approximation Possibilities for Wind BC for Increasing Capacity

# # R Plot
# 
# par(mar=c(5,5,5,5))
# plot(balancing[1:4,1],balancing[1:4,2],col=col_nr[1],xlim=c(0,32),ylim=c(0,6),cex.axis=1.5,cex.lab=2,cex.main=2,xlab="Share of Wind Production [%]",ylab="System Costs [EUR/MWh_e]",main="Additional System Costs via Wind Integration",lwd=2,pch=3,cex.lab=1.2,cex.main=1.2)
# points(balancing[5:7,1],balancing[5:7,2],col=col_bg[1],lwd=2,pch=3)
# points(balancing[8,1],balancing[8,2],col=col_nr[3],lwd=2,pch=3)
# points(balancing[9:10,1],balancing[9:10,2],col=col_bg[5],lwd=2,pch=3)
# points(balancing[11:13,1],balancing[11:13,2],col=col_nr[4],lwd=2,pch=3)
# points(balancing[14:15,1],balancing[14:15,2],col=col_b50[30],lwd=2,pch=3)
# points(balancing[16,1],balancing[16,2],col=col_nr[5],lwd=2,pch=3)
# points(balancing[17:18,1],balancing[17:18,2],col=col_b50[50],lwd=2,pch=3)
# points(balancing[19:20,1],balancing[19:20,2],col="blue",lwd=2,pch=3)
# lines(sapply(seq(0,30,1),function(x) slope_BC*x/100),type="l",lwd=2,col=colors[50])
# #lines(sapply(seq(0,30,1),function(x) lmi_BC[1]+lmi_BC[2]*x),type="l",lwd=2,col=colors[37])
# lines(sapply(seq(0,30,1),function(x) sq_BC[1]*x+sq_BC[2]*x^2),type="l",lwd=2,col=colors[24])
# #lines(sapply(seq(0,30,1),function(x) sqi_BC[1]+sqi_BC[2]*x+sqi_BC[3]*x^2),type="l",lwd=2,col=colors[11])
# legend("bottomright",c("linear","quadratic"),col=c(col_b50[50],col_b50[24]),lwd=3,cex=1.1)
# 
# # PDF Plot
# 
# pdf("system_costs.pdf",width=8,height=5)
# par(mar=c(5,5,5,5))
# plot(balancing[1:4,1],balancing[1:4,2],col=col_nr[1],xlim=c(0,32),ylim=c(0,6),cex.axis=1.5,cex.lab=2,cex.main=2,xlab="Share of Wind Production [%]",ylab="System Costs [EUR/MWh_e]",main="Additional System Costs via Wind Integration",lwd=2,pch=3,cex.lab=1.2,cex.main=1.2)
# points(balancing[5:7,1],balancing[5:7,2],col=col_bg[1],lwd=2,pch=3)
# points(balancing[8,1],balancing[8,2],col=col_nr[3],lwd=2,pch=3)
# points(balancing[9:10,1],balancing[9:10,2],col=col_bg[5],lwd=2,pch=3)
# points(balancing[11:13,1],balancing[11:13,2],col=col_nr[4],lwd=2,pch=3)
# points(balancing[14:15,1],balancing[14:15,2],col=col_b50[30],lwd=2,pch=3)
# points(balancing[16,1],balancing[16,2],col=col_nr[5],lwd=2,pch=3)
# points(balancing[17:18,1],balancing[17:18,2],col=col_b50[50],lwd=2,pch=3)
# points(balancing[19:20,1],balancing[19:20,2],col="blue",lwd=2,pch=3)
# lines(sapply(seq(0,30,1),function(x) slope_BC*x/100),type="l",lwd=2,col=colors[50])
# #lines(sapply(seq(0,30,1),function(x) lmi_BC[1]+lmi_BC[2]*x),type="l",lwd=2,col=colors[37])
# lines(sapply(seq(0,25,1),function(x) sq_BC[1]*x+sq_BC[2]*x^2),type="l",lwd=2,col=colors[24])
# #lines(sapply(seq(0,30,1),function(x) sqi_BC[1]+sqi_BC[2]*x+sqi_BC[3]*x^2),type="l",lwd=2,col=colors[11])
# legend("bottomright",c("linear","linear-quadratic"),col=c(col_b50[50],col_b50[24]),lwd=3,cex=1.1)
# dev.off()
