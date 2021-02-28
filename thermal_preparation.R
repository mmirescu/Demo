# ===== Preparating and Analysing the Fuel Price Data for the Thermal Technologies =====

setwd("C:/Users/Mirescu1/ucloud/Universities/TU Wien/Theses/Dissertation/Eigene/Latex")

# === Install the necessary packages ===

library(xts)
library(quadprog)
library(rootSolve)
library(MASS)
library("colorspace")
library(grDevices)

# === Define Colors to Be Used Throughout This Document ===

# Call Up an Interface Where One Can Define a Personalised Color Palette With Subsequent Color Generator in HEX-Format

#choose_palette()

# Introduce the personalised Palettes into R

# Nikolaus Rab (nr) Color Palette

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

col_b50<-c("#023FA5","#1041A5","#1A44A5","#2146A5","#2749A5","#2C4BA6","#314EA6","#3550A7","#3953A8","#3D56A9","#4158AA","#455BAB","#485DAC","#4C60AD",
           "#4F63AE","#5365AF","#5668B0","#5A6BB2","#5D6EB3","#6170B5","#6473B6","#6876B8","#6B79B9","#6E7CBB","#727FBC","#7582BE","#7985C0","#7D88C1",
           "#808BC3","#848FC5","#8892C7","#8B95C9","#8F99CB","#939CCD","#97A0CF","#9BA3D1","#9FA7D3","#A3ABD5","#A8AFD7","#ACB3DA","#B1B7DC","#B6BCDF",
           "#BBC0E1","#C0C5E4","#C5CAE7","#CBD0EA","#D2D5ED","#D9DCF1","#E1E4F5","#EEF0FC")

# Define the Colors Used in the Plot to Be of the Preferred Palette

col_bg<-col_b50


# === Read in the Necessary Data ===

data<-read.csv("data.csv",header=T,sep=",",dec=".")
parameters<-read.csv("data_parameters.csv",header=T,sep=",",dec=".")
mean_co2<-mean(data[80:240,7])
co2_plus<-50
mean_co2_plus<-mean_co2+co2_plus

# Plot the CO2 Prices

# # R Plot
# 
# plot(data[80:240,7]*data[80:240,5],type="l",lwd=4,col=col_b30[30],axes=F,main="Evolution of the Monthly CO2 Futures Price (August 2005 - December 2018)",xlab="Time [Years]",ylab="CO2 Price [EUR/mt]")
# axis(1,at=seq(1,161,12),labels=c("05","06","07","08","09","10","11","12","13","14","15","16","17","18"))
# axis(2)
# abline(h=mean(data[80:240,7]*data[80:240,5]),lwd=3,lty=2,col=col_b30[30])
# abline(h=sd(data[80:240,7]*data[80:240,5]),lwd=3,lty=3,col=col_b30[30])
# box()
# 
# # PDF Plot
# 
# pdf("prices_co2.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# plot(data[80:240,7]*data[80:240,5],type="l",lwd=4,col=col_b30[30],axes=F,main="Evolution of the Monthly CO2 Futures Price (August 2005 - December 2018)",xlab="Time [Years]",ylab="CO2 Price [EUR/mt]",cex.main=2.5,cex.axis=2,cex.lab=2)
# axis(1,at=seq(1,161,12),labels=c("05","06","07","08","09","10","11","12","13","14","15","16","17","18"))
# axis(2)
# abline(h=mean(data[80:240,7]*data[80:240,5]),lwd=3,lty=2,col=col_b30[30])
# abline(h=sd(data[80:240,7]*data[80:240,5]),lwd=3,lty=3,col=col_b30[30])
# box()
# dev.off()

# Define the Efficiency / Energy Conversion Factors 

eta<-parameters[1:(dim(parameters)[1]-1),dim(parameters)[2]-4]

# === Data Preparation ===

# Transform Everything in EUR/MWh_e

# = Gas =

# 1 MWh: 0.0036 TJ ; H_s in H_i: 0.903 ; H_i in H_s: 1/0.903 ; HICP: data[,5] ; efficiency factor: eta[1]
# incl. CO2: data[,7]; HICP: data[,5]; Emissions in kgCO2/kWh: 0.2; efficiency factor: eta[1]

gas_no<-data[,2]*0.0036*1/0.901*1/data[,5]*1/eta[1]
gas<-gas_no+data[,7]*data[,5]*0.2*1/eta[1]
gas_co2<-gas_no+(data[,7]*data[,5]+co2_plus)*0.2*1/eta[1]

# Plot gas prices in EUR/MWh:
# 1. H_i with HICP
# 2. H_i without HICP
# 3. H_s with HICP
# 4. H_s without HICP

# # R Plot
# 
# par(mfrow=c(1,1))
# plot(gas,type="l",lwd=4,col=col_b30[30],axes=F,main="Evolution of the Monthly Gas Cross-Border Price in Germany (January 1999 - December 2018)",xlab="Time [Years]",ylim=c(10,max(gas)),ylab="Gas Price [EUR/MWh]")
# axis(1,at=seq(6,dim(data)[1],12),labels=c("99","00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18"))
# axis(2)
# lines(data[,2]*0.0036*1/0.901*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1],lwd=4,col=col_b30[27])
# lines(data[,2]*0.0036*1/data[,5]*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1],lwd=4,col=col_b30[24])
# lines(data[,2]*0.0036*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1],lwd=4,col=col_b30[21])
# abline(h=mean(gas),lwd=3,lty=1,col=col_b30[30])
# abline(h=sd(gas),lwd=3,lty=3,col=col_b30[30])
# box()
# legend("topleft",c("LHV+HICP","LHV-HICP","HHV+HICP","HHV-HICP"),col=c(col_b30[30],col_b30[27],col_b30[24],col_b30[21]),lwd=4)
# 
# # PDF Plot
# 
# pdf("prices_gas.pdf",width=14,height=8)
# plot(gas,type="l",lwd=4,col=col_b30[30],axes=F,main="Evolution of the Monthly Gas Cross-Border Price in Germany (January 1999 - December 2018)",xlab="Time [Years]",ylim=c(10,max(gas)),ylab="Gas Price [EUR/MWh]")
# axis(1,at=seq(6,dim(data)[1],12),labels=c("99","00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18"))
# axis(2)
# lines(data[,2]*0.0036*1/0.901*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1],lwd=4,col=col_b30[27])
# lines(data[,2]*0.0036*1/data[,5]*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1],lwd=4,col=col_b30[24])
# lines(data[,2]*0.0036*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1],lwd=4,col=col_b30[21])
# abline(h=mean(gas),lwd=3,lty=1,col=col_b30[30])
# abline(h=sd(gas),lwd=3,lty=3,col=col_b30[30])
# box()
# legend("topleft",c("LHV+HICP","LHV-HICP","HHV+HICP","HHV-HICP"),col=c(col_b30[30],col_b30[27],col_b30[24],col_b30[21]),lwd=4)
# dev.off()


# # Boxplot(s)
# 
# # R Plot
# 
# boxplot(gas,data[,2]*0.0036*1/0.901*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1],data[,2]*0.0036*1/data[,6]*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1],data[,2]*0.0036*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1],axes=F,xlab="Different Possibilities",ylab="Gas Price [EUR/MWh]",main="Monthly Gas Cross-Border Price in Germany (January 1999 - December 2018)",col=c(col_b30[30],col_b30[27],col_b30[24],col_b30[21]))
# means<-c(mean(gas),mean(data[,2]*0.0036*1/0.901*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1]),mean(data[,2]*0.0036*1/data[,6]*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1]),mean(data[,2]*0.0036*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1]))
# points(means,col="blue",pch=18)
# sds<-c(sd(gas),sd(data[,2]*0.0036*1/0.901*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1]),sd(data[,2]*0.0036*1/data[,6]*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1]),sd(data[,2]*0.0036*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1]))
# points(sds,col="lightblue",pch=18)
# axis(1,at=seq(4),labels=c("LHV+HICP","LHV-HICP","HHV+HICP","HHV-HICP"))
# axis(2)
# box()
# legend("topright",c("Mean","Standard Deviation"),col=c("blue","lightblue"),lwd=4)
# 
# # PDF Plot
# 
# pdf("box_gas.pdf",width=14,height=8)
# boxplot(gas,data[,2]*0.0036*1/0.901*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1],data[,2]*0.0036*1/data[,6]*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1],data[,2]*0.0036*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1],axes=F,xlab="Different Possibilities",ylab="Gas Price [EUR/MWh]",main="Monthly Gas Cross-Border Price in Germany (January 1999 - December 2018)",col=c(col_b30[30],col_b30[27],col_b30[24],col_b30[21]))
# means<-c(mean(gas),mean(data[,2]*0.0036*1/0.901*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1]),mean(data[,2]*0.0036*1/data[,6]*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1]),mean(data[,2]*0.0036*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1]))
# points(means,col="blue",pch=18)
# sds<-c(sd(gas),sd(data[,2]*0.0036*1/0.901*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1]),sd(data[,2]*0.0036*1/data[,6]*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1]),sd(data[,2]*0.0036*1/eta[1]+data[,7]*data[,5]*0.2*1/eta[1]))
# points(sds,col="lightblue",pch=18)
# axis(1,at=seq(4),labels=c("LHV+HICP","LHV-HICP","HHV+HICP","HHV-HICP"))
# axis(2)
# box()
# legend("topright",c("Mean","Standard Deviation"),col=c("blue","lightblue"),lwd=4)
# dev.off()

# = Coal =

# South African Coal Spot Price ; 1 MWh: 0.122835033 mt ; H_s in H_i: 0.98 ; H_i in H_s: 1/0.98 ; HICP: data[,6] ; USD/EUR: data[,7] ; efficiency factor: eta[2]
# incl. CO2: data[,7]; HICP: data[,5]; Emissions in kgCO2/kWh: 0.34; efficiency factor: eta[2]
# Konstantin: 1 mt=> MWh_th: 1/7=0.1428571

coal_no<-1/eta[2]*data[,3]*0.122835033*1/data[,5]*data[,6]
coal<-coal_no+data[,7]*data[,5]*0.34*1/eta[2]
coal_co2<-coal_no+(data[,7]*data[,5]+co2_plus)*0.34*1/eta[2]

# Plot coal prices in EUR/MWh:
# 1. H_i with HICP
# 2. H_i without HICP
# 3. H_s with HICP
# 4. H_s without HICP

# # R Plot
# 
# plot(coal,type="l",lwd=4,col=col_b30[30],axes=F,main="Evolution of the Monthly Coal Price in Germany (January 1999 - December 2018)",xlab="Time [Years]",ylim=c(min(sd(coal),coal,1/eta[2]*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7]),1/eta[2]*0.98*c(data[1:193,3]*0.122835033*1/data[1:193,6],data[194:221,4]*0.122835033*1/data[194:221,6]*data[194:221,7]),1/eta[2]*0.98*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7])),max(coal,1/eta[2]*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7]),1/eta[2]*0.98*c(data[1:193,3]*0.122835033*1/data[1:193,6],data[194:221,4]*0.122835033*1/data[194:221,6]*data[194:221,7]),1/eta[2]*0.98*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7]))),ylab="Coal Price [EUR/MWh]")
# axis(1,at=seq(6,dim(data)[1],12),labels=c("99","00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18"))
# axis(2)
# #lines(1/eta[2]*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7]),lwd=4,col=col_b30[27])
# #lines(1/eta[2]*0.98*c(data[1:193,3]*0.122835033*1/data[1:193,6],data[194:221,4]*0.122835033*1/data[194:221,6]*data[194:221,7]),lwd=4,col=col_b30[24])
# #lines(1/eta[2]*0.98*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7]),lwd=4,col=col_b30[21])
# abline(h=mean(coal),lwd=3,lty=1,col=col_b30[30])
# abline(h=sd(coal),lwd=3,lty=3,col=col_b30[30])
# box()
# legend("topleft",c("LHV+HICP","LHV-HICP","HHV+HICP","HHV-HICP"),col=c(col_b30[30],col_b30[27],col_b30[24],col_b30[21]),lwd=4)
# 
# # PDF Plot
# 
# pdf("prices_coal.pdf",width=14,height=8)
# plot(coal,type="l",lwd=4,col=col_b30[30],axes=F,main="Evolution of the Monthly Coal Price in Germany (January 1999 - December 2018)",xlab="Time [Years]",ylim=c(min(sd(coal),coal,1/eta[2]*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7]),1/eta[2]*0.98*c(data[1:193,3]*0.122835033*1/data[1:193,6],data[194:221,4]*0.122835033*1/data[194:221,6]*data[194:221,7]),1/eta[2]*0.98*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7])),max(coal,1/eta[2]*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7]),1/eta[2]*0.98*c(data[1:193,3]*0.122835033*1/data[1:193,6],data[194:221,4]*0.122835033*1/data[194:221,6]*data[194:221,7]),1/eta[2]*0.98*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7]))),ylab="Coal Price [EUR/MWh]")
# axis(1,at=seq(6,dim(data)[1],12),labels=c("99","00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18"))
# axis(2)
# #lines(1/eta[2]*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7]),lwd=4,col=col_b30[27])
# #lines(1/eta[2]*0.98*c(data[1:193,3]*0.122835033*1/data[1:193,6],data[194:221,4]*0.122835033*1/data[194:221,6]*data[194:221,7]),lwd=4,col=col_b30[24])
# #lines(1/eta[2]*0.98*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7]),lwd=4,col=col_b30[21])
# abline(h=mean(coal),lwd=3,lty=1,col=col_b30[30])
# abline(h=sd(coal),lwd=3,lty=3,col=col_b30[30])
# box()
# legend("topleft",c("LHV+HICP","LHV-HICP","HHV+HICP","HHV-HICP"),col=c(col_b30[30],col_b30[27],col_b30[24],col_b30[21]),lwd=4)
# dev.off()

## Boxplot(s)

## R Plot

#boxplot(coal,1/eta[2]*c(data[1:193,3]*0.122835033*1/0.98,data[194:221,4]*0.122835033*1/0.98*data[194:221,7]),1/eta[2]*c(data[1:193,3]*0.122835033*1/data[1:193,6],data[194:221,4]*0.122835033*1/data[194:221,6]*data[194:221,7]),1/eta[2]*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7]),axes=F,xlab="Different Possibilities",ylab="Coal Price [EUR/MWh]",main="Coal Price in Germany (January 1999 - May 2017)",col=c(col_b30[30],col_b30[27],col_b30[24],col_b30[21]),ylim=c(3,27))
#means_c<-c(mean(coal),mean(1/eta[2]*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7])),mean(1/eta[2]*0.98*c(data[1:193,3]*0.122835033*1/data[1:193,6],data[194:221,4]*0.122835033*1/data[194:221,6]*data[194:221,7])),mean(1/eta[2]*0.98*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7])))
#points(means_c,col="blue",pch=18)
#sds_c<-c(sd(coal),sd(1/eta[2]*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7])),sd(1/eta[2]*0.98*c(data[1:193,3]*0.122835033*1/data[1:193,6],data[194:221,4]*0.122835033*1/data[194:221,6]*data[194:221,7])),sd(1/eta[2]*0.98*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7])))
#points(sds_c,col="lightblue",pch=18)
#axis(1,at=seq(4),labels=c("LHV+HICP","LHV-HICP","HHV+HICP","HHV-HICP"))
#axis(2)
#box()
#legend("topright",c("Mean","Standard Deviation"),col=c("blue","lightblue"),lwd=4)

## PDF Plot

#pdf("box_coal.pdf",width=14,height=8)
#boxplot(coal,1/eta[2]*c(data[1:193,3]*0.122835033*1/0.98,data[194:221,4]*0.122835033*1/0.98*data[194:221,7]),1/eta[2]*c(data[1:193,3]*0.122835033*1/data[1:193,6],data[194:221,4]*0.122835033*1/data[194:221,6]*data[194:221,7]),1/eta[2]*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7]),axes=F,xlab="Different Possibilities",ylab="Coal Price [EUR/MWh]",main="Coal Price in Germany (January 1999 - May 2017)",col=c(col_b30[30],col_b30[27],col_b30[24],col_b30[21]),ylim=c(3,27))
#means_c<-c(mean(coal),mean(1/eta[2]*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7])),mean(1/eta[2]*0.98*c(data[1:193,3]*0.122835033*1/data[1:193,6],data[194:221,4]*0.122835033*1/data[194:221,6]*data[194:221,7])),mean(1/eta[2]*0.98*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7])))
#points(means_c,col="blue",pch=18)
#sds_c<-c(sd(coal),sd(1/eta[2]*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7])),sd(1/eta[2]*0.98*c(data[1:193,3]*0.122835033*1/data[1:193,6],data[194:221,4]*0.122835033*1/data[194:221,6]*data[194:221,7])),sd(1/eta[2]*0.98*c(data[1:193,3]*0.122835033,data[194:221,4]*0.122835033*data[194:221,7])))
#points(sds_c,col="lightblue",pch=18)
#axis(1,at=seq(4),labels=c("LHV+HICP","LHV-HICP","HHV+HICP","HHV-HICP"))
#axis(2)
#box()
#legend("topright",c("Mean","Standard Deviation"),col=c("blue","lightblue"),lwd=4)
#dev.off()

# = Uranium =

# 1 kg : 2.20462262 pounds ; 1 MBtu: 1/180 lbs ; 1 MWh: 3.41214163 MBtu ; HICP: data[,5] ; USD/EUR: data[,6]
# Konstantin transformation: 1 lbs->MWh_th: 55.6

#old
uranium<-(data[,4])*3.41214163*2.20462262*1/180*1/data[,5]*data[,6]*1/eta[3]
#new
#uranium<-(data[,4]+49.1)*1/55.6*1/data[,5]*data[,6]*1/eta[3]

# Plot uranium prices in EUR/MWh:
# 1. Price with HICP
# 2. Price without HICP

# # R plot
# 
# plot(uranium,type="l",lwd=4,col=col_b30[30],axes=F,main="Evolution of the Uranium Price in Germany (January 1999 - December 2018)",xlab="Time [Years]",ylim=c(min(uranium,data[,4]*3.41214163*2.20462262*1/180*data[,6]*1/eta[3]),max(uranium,data[,4]*3.41214163*2.20462262*1/180*data[,6]*1/eta[3])),ylab="Uranium Price [EUR/MWh]")
# axis(1,at=seq(1,dim(data)[1],12),labels=c("99","00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18"))
# axis(2)
# lines(data[,4]*3.41214163*2.20462262*1/180*data[,6]*1/eta[3],lwd=4,col=col_b30[27])
# abline(h=mean(uranium),lwd=3,lty=1,col=col_b30[30])
# abline(h=sd(uranium),lwd=3,lty=3,col=col_b30[30])
# box()
# legend("topleft",c("+HICP","-HICP"),col=c(col_b30[30],col_b30[27]),lwd=4)
# 
# # PDF Plot
# 
# pdf("prices_uranium.pdf",width=14,height=8)
# plot(uranium,type="l",lwd=4,col=col_b30[30],axes=F,main="Evolution of the Uranium Price in Germany (January 1999 - December 2018)",xlab="Time [Years]",ylim=c(min(uranium,data[,4]*3.41214163*2.20462262*1/180*data[,6]*1/eta[3]),max(uranium,data[,4]*3.41214163*2.20462262*1/180*data[,6]*1/eta[3])),ylab="Uranium Price [EUR/MWh]")
# axis(1,at=seq(1,dim(data)[1],12),labels=c("99","00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18"))
# axis(2)
# lines(data[,4]*3.41214163*2.20462262*1/180*data[,6]*1/eta[3],lwd=4,col=col_b30[27])
# abline(h=mean(uranium),lwd=3,lty=1,col=col_b30[30])
# abline(h=sd(uranium),lwd=3,lty=3,col=col_b30[30])
# box()
# legend("topleft",c("+HICP","-HICP"),col=c(col_b30[30],col_b30[27]),lwd=4)
# dev.off()
# 
# # Boxplot(s)
# 
# # R Plot
# 
# boxplot(uranium,data[,4]*3.41214163*2.20462262*1/180*data[,6]*1/eta[3],axes=F,xlab="Different Possibilities",ylab="Uranium Price [EUR/MWh]",main="Uranium Price in Germany (January 1999 - May 2017)",col=c(col_b30[30],col_b30[27]))
# means_u<-c(mean(uranium),mean(data[,4]*3.41214163*2.20462262*1/180*data[,6]*1/eta[3]))
# points(means_u,col="blue",pch=18)
# sds_u<-c(sd(uranium),sd(data[,4]*3.41214163*2.20462262*1/180*data[,6]*1/eta[3]))
# points(sds_u,col="lightblue",pch=18)
# axis(1,at=seq(2),labels=c("+HICP","-HICP"))
# axis(2)
# box()
# legend("topright",c("Mean","Standard Deviation"),col=c("blue","lightblue"),lwd=4)
# 
# # PDF Plot
# 
# pdf("box_uranium.pdf",width=14,height=8)
# boxplot(uranium,data[,4]*3.41214163*2.20462262*1/180*data[,6]*1/eta[3],axes=F,xlab="Different Possibilities",ylab="Uranium Price [EUR/MWh]",main="Uranium Price in Germany (January 1999 - May 2017)",col=c(col_b30[30],col_b30[27]))
# means_u<-c(mean(uranium),mean(data[,4]*3.41214163*2.20462262*1/180*data[,6]*1/eta[3]))
# points(means_u,col="blue",pch=18)
# sds_u<-c(sd(uranium),sd(data[,4]*3.41214163*2.20462262*1/180*data[,6]*1/eta[3]))
# points(sds_u,col="lightblue",pch=18)
# axis(1,at=seq(2),labels=c("+HICP","-HICP"))
# axis(2)
# box()
# legend("topright",c("Mean","Standard Deviation"),col=c("blue","lightblue"),lwd=4)
# dev.off()

# === Create the Matrix of the Thermal Technology as Time Series in EUR/MWh ===

therm_no<-cbind(gas_no,coal_no,uranium)
colnames(therm_no)<-c("Gas","Coal","Uranium")
rownames(therm_no)<-data[,1]

therm<-cbind(gas,coal,uranium)
colnames(therm)<-c("Gas","Coal","Uranium")
rownames(therm)<-data[,1]

therm_co2<-cbind(gas_co2,coal_co2,uranium)
colnames(therm_co2)<-c("Gas","Coal","Uranium")
rownames(therm_co2)<-data[,1]

# Plot the Prices of All Three Thermal Technologies 

# # No CO2 Prices
# 
# # R Plot
# 
# par(mar=c(5,5,5,5))
# plot(therm_no[,1],type="l",lwd=4,col=col_bg[1],cex.axis=2,cex.lab=2,cex.main=2,axes=F,main="Evolution of the Input Prices in Germany (January 1999 - December 2018)",xlab="Time [Years]",ylim=c(min(therm_no[,1],therm_no[,2],therm_no[,3]),max(therm_no[,1],therm_no[,2],therm_no[,3])),ylab="Input Prices [EUR/MWh]",cex.main=2,cex.axis=2,cex.lab=2)
# axis(1,at=seq(6,dim(data)[1],12),labels=c("99","00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18"))
# axis(2)
# lines(therm_no[,2],lwd=4,col=col_bg[21])
# lines(therm_no[,3],lwd=4,col=col_bg[41])
# abline(h=mean(therm_no[,1]),lwd=3,lty=2,col=col_bg[1])
# abline(h=mean(therm_no[,2]),lwd=3,lty=2,col=col_bg[21])
# abline(h=mean(therm_no[,3]),lwd=3,lty=2,col=col_bg[41])
# abline(h=sd(therm[,1]),lwd=3,lty=3,col=col_bg[1])
# abline(h=sd(therm[,2]),lwd=3,lty=3,col=col_bg[21])
# abline(h=sd(therm[,3]),lwd=3,lty=3,col=col_bg[41])
# box()
# legend("topleft",c("Gas","Coal","Uranium"),col=c(col_bg[1],col_bg[21],col_bg[41]),lwd=4,cex=2)
# 
# # # PDF Plot
# # 
# pdf("prices_thermal_no.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# plot(therm_no[,1],type="l",lwd=4,col=col_bg[1],cex.axis=2,cex.lab=2,cex.main=2,axes=F,main="Evolution of the Input Prices in Germany (01.1999 - 12.2018)",xlab="Time [Years]",ylim=c(min(therm_no[,1],therm_no[,2],therm_no[,3]),max(therm_no[,1],therm_no[,2],therm_no[,3])),ylab="Input Prices [EUR/MWh_e]")
# #axis(1,at=seq(6,dim(data)[1],12),labels=c("99","00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18"))
# #axis(2)
# lines(therm_no[,2],lwd=4,col=col_bg[21])
# lines(therm_no[,3],lwd=4,col=col_bg[41])
# #abline(h=mean(therm_no[,1]),lwd=3,lty=2,col=col_bg[1])
# #abline(h=mean(therm_no[,2]),lwd=3,lty=2,col=col_bg[21])
# #abline(h=mean(therm_no[,3]),lwd=3,lty=2,col=col_bg[41])
# #abline(h=sd(therm[,1]),lwd=3,lty=3,col=col_bg[1])
# #abline(h=sd(therm[,2]),lwd=3,lty=3,col=col_bg[21])
# #abline(h=sd(therm[,3]),lwd=3,lty=3,col=col_bg[41])
# #box()
# legend("topleft",c("Gas","Coal","Uranium"),col=c(col_bg[1],col_bg[21],col_bg[41]),lwd=4,cex=2)
# dev.off()
# 
# # With CO2 Prices
# 
# # R Plot
# 
# par(mar=c(5,5,5,5))
# plot(therm[,1],type="l",lwd=4,col=col_bg[1],cex.axis=2,cex.lab=2,cex.main=2,axes=F,main="Evolution of the Input Prices in Germany (01.1999 - 12.2018)",xlab="Time [Years]",ylim=c(min(therm[,1],therm[,2],therm[,3]),max(therm[,1],therm[,2],therm[,3])),ylab="Input Prices [EUR/MWh]",cex.main=2,cex.axis=2,cex.lab=2)
# #axis(1,at=seq(6,dim(data)[1],12),labels=c("99","00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18"))
# #axis(2)
# lines(therm[,2],lwd=4,col=col_bg[21])
# lines(therm[,3],lwd=4,col=col_bg[41])
# abline(h=mean(therm[,1]),lwd=3,lty=2,col=col_bg[1])
# abline(h=mean(therm[,2]),lwd=3,lty=2,col=col_bg[21])
# abline(h=mean(therm[,3]),lwd=3,lty=2,col=col_bg[41])
# abline(h=sd(therm[,1]),lwd=3,lty=3,col=col_bg[1])
# abline(h=sd(therm[,2]),lwd=3,lty=3,col=col_bg[21])
# abline(h=sd(therm[,3]),lwd=3,lty=3,col=col_bg[41])
# box()
#legend("topleft",c("Gas","Coal","Uranium"),col=c(col_bg[1],col_bg[21],col_bg[41]),lwd=4,cex=2)

# PDF Plot
# 
# pdf("prices_thermal.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# plot(therm[,1],type="l",lwd=4,col=col_bg[1],cex.axis=2,cex.lab=2,cex.main=2,axes=F,main="Evolution of the Monthly Input Prices in Germany (January 1999 - December 2018)",xlab="Time [Years]",ylim=c(min(therm[,1],therm[,2],therm[,3]),max(therm[,1],therm[,2],therm[,3])),ylab="Input Prices [EUR/MWh_e]",cex.main=2,cex.axis=2,cex.lab=2)
# axis(1,at=seq(6,dim(data)[1],12),labels=c("99","00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18"))
# axis(2)
# lines(therm[,2],lwd=4,col=col_bg[21])
# lines(therm[,3],lwd=4,col=col_bg[41])
# abline(h=mean(therm[,1]),lwd=3,lty=2,col=col_bg[1])
# abline(h=mean(therm[,2]),lwd=3,lty=2,col=col_bg[21])
# abline(h=mean(therm[,3]),lwd=3,lty=2,col=col_bg[41])
# abline(h=sd(therm[,1]),lwd=3,lty=3,col=col_bg[1])
# abline(h=sd(therm[,2]),lwd=3,lty=3,col=col_bg[21])
# abline(h=sd(therm[,3]),lwd=3,lty=3,col=col_bg[41])
# box()
# legend("topleft",c("Gas","Coal","Uranium"),col=c(col_bg[1],col_bg[21],col_bg[41]),lwd=4,cex=2)
# dev.off()

# # R Plot
# 
# par(mar=c(5,5,5,5))
# plot(therm[,1],type="l",lwd=4,col=col_bg[1],cex.axis=2,cex.lab=2,cex.main=2,axes=F,main="Evolution of the Input Prices in Germany (January 1999 - December 2018)",xlab="Time [Years]",ylim=c(min(therm[,1],therm[,2],therm[,3]),max(therm[,1],therm[,2],therm[,3])),ylab="Input Prices [EUR/MWh]",cex.main=2,cex.axis=2,cex.lab=2)
# axis(1,at=seq(6,dim(data)[1],12),labels=c("99","00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18"))
# axis(2)
# lines(therm[,2],lwd=4,col=col_bg[21])
# lines(therm[,3],lwd=4,col=col_bg[41])
# abline(h=mean(therm[,1]),lwd=3,lty=3,col=col_bg[1])
# abline(h=mean(therm[,2]),lwd=3,lty=3,col=col_bg[21])
# abline(h=mean(therm[,3]),lwd=3,lty=3,col=col_bg[41])
# #abline(h=sd(therm[,1]),lwd=3,lty=3,col=col_b30[30])
# #abline(h=sd(therm[,2]),lwd=3,lty=3,col=col_b30[27])
# #abline(h=sd(therm[,3]),lwd=3,lty=3,col=col_b30[24])
# box()
# #legend("topleft",c("Gas","Coal","Uranium"),col=c(col_bg[1],col_bg[21],col_bg[41]),lwd=4,cex=2)
# 
# # PDF Plot
# 
# pdf("prices_th.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# par(mfrow=c(2,1))
# plot(therm_no[,1],type="l",lwd=4,col=col_bg[1],cex.axis=2,cex.lab=2,cex.main=2,axes=F,main="(a) No EUAs Included",xlab="Time [Years]",ylim=c(min(therm_no[,1],therm_no[,2],therm_no[,3]),max(therm_no[,1],therm_no[,2],therm_no[,3])),ylab="Input Prices [EUR/MWh]",cex.main=2,cex.axis=2,cex.lab=2)
# axis(1,at=seq(6,dim(data)[1],12),labels=c("99","00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18"))
# axis(2)
# lines(therm_no[,2],lwd=4,col=col_bg[21])
# lines(therm_no[,3],lwd=4,col=col_bg[41])
# abline(h=mean(therm_no[,1]),lwd=3,lty=3,col=col_bg[1])
# abline(h=mean(therm_no[,2]),lwd=3,lty=3,col=col_bg[21])
# abline(h=mean(therm_no[,3]),lwd=3,lty=3,col=col_bg[41])
# #abline(h=sd(therm[,1]),lwd=3,lty=3,col=col_b30[30])
# #abline(h=sd(therm[,2]),lwd=3,lty=3,col=col_b30[27])
# #abline(h=sd(therm[,3]),lwd=3,lty=3,col=col_b30[24])
# box()
# #legend("topleft",c("Gas","Coal","Uranium"),col=c(col_bg[1],col_bg[21],col_bg[41]),lwd=4,cex=2)
# plot(therm[,1],type="l",lwd=4,col=col_bg[1],cex.axis=2,cex.lab=2,cex.main=2,axes=F,main="(b) EUAs Included",xlab="Time [Years]",ylim=c(min(therm[,1],therm[,2],therm[,3]),max(therm[,1],therm[,2],therm[,3])),ylab="Input Prices [EUR/MWh]",cex.main=2,cex.axis=2,cex.lab=2)
# axis(1,at=seq(6,dim(data)[1],12),labels=c("99","00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18"))
# axis(2)
# lines(therm[,2],lwd=4,col=col_bg[21])
# lines(therm[,3],lwd=4,col=col_bg[41])
# abline(h=mean(therm[,1]),lwd=3,lty=3,col=col_bg[1])
# abline(h=mean(therm[,2]),lwd=3,lty=3,col=col_bg[21])
# abline(h=mean(therm[,3]),lwd=3,lty=3,col=col_bg[41])
# #abline(h=sd(therm[,1]),lwd=3,lty=3,col=col_b30[30])
# #abline(h=sd(therm[,2]),lwd=3,lty=3,col=col_b30[27])
# #abline(h=sd(therm[,3]),lwd=3,lty=3,col=col_b30[24])
# box()
# #legend("topleft",c("Gas","Coal","Uranium"),col=c(col_bg[1],col_bg[21],col_bg[41]),lwd=4,cex=2)
# dev.off()

# # Plot the Normed Prices
# 
# # With CO2 Prices 
# 
# # R Plot
# 
# plot(therm_no[,1]/therm_no[1,1],type="l",lwd=4,col=col_b30[30],axes=F,main="Evolution of the Thermal Input Prices in Germany (January 1999 - December 2018)",xlab="Time [Years]",ylim=c(min(therm_no[,1]/therm_no[1,1],therm_no[,2]/therm_no[1,2],therm_no[,3]/therm_no[1,3],sd(therm_no[,1]/therm_no[1,1]),sd(therm_no[,2]/therm_no[1,2]),sd(therm_no[,3]/therm_no[1,3])),max(therm_no[,1]/therm_no[1,1],therm_no[,2]/therm_no[1,2],therm_no[,3]/therm_no[1,3])),ylab="Input Prices [EUR/MWh]")
# axis(1,at=seq(6,dim(data)[1],12),labels=c("99","00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18"))
# axis(2)
# lines(therm_no[,2]/therm_no[1,2],lwd=4,col=col_b30[27])
# lines(therm_no[,3]/therm_no[1,3],lwd=4,col=col_b30[24])
# #bline(h=mean(therm[,1]/therm[1,1]),lwd=3,lty=1,col=col_b30[30])
# #abline(h=mean(therm[,2]/therm[1,2]),lwd=3,lty=1,col=col_b30[27])
# #abline(h=mean(therm[,3]/therm[1,3]),lwd=3,lty=1,col=col_b30[24])
# abline(h=sd(therm_no[,1]/therm_no[1,1]),lwd=3,lty=3,col=col_b30[30])
# abline(h=sd(therm_no[,2]/therm_no[1,2]),lwd=3,lty=3,col=col_b30[27])
# abline(h=sd(therm_no[,3]/therm_no[1,3]),lwd=3,lty=3,col=col_b30[24])
# box()
# legend("topleft",c("Gas","Coal","Uranium"),col=c(col_b30[30],col_b30[27],col_b30[24]),lwd=4)
# 
# # PDF Plot
# 
# pdf("prices_all_normed_no.pdf",width=14,height=8)
# plot(therm_no[,1]/therm_no[1,1],type="l",lwd=4,col=col_b30[30],axes=F,main="Evolution of the Thermal Input Prices in Germany (January 1999 - December 2018)",xlab="Time [Years]",ylim=c(min(therm_no[,1]/therm_no[1,1],therm_no[,2]/therm_no[1,2],therm_no[,3]/therm_no[1,3],sd(therm_no[,1]/therm_no[1,1]),sd(therm_no[,2]/therm_no[1,2]),sd(therm_no[,3]/therm_no[1,3])),max(therm_no[,1]/therm_no[1,1],therm_no[,2]/therm_no[1,2],therm_no[,3]/therm_no[1,3])),ylab="Input Prices [EUR/MWh]")
# axis(1,at=seq(6,dim(data)[1],12),labels=c("99","00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18"))
# axis(2)
# lines(therm_no[,2]/therm_no[1,2],lwd=4,col=col_b30[27])
# lines(therm_no[,3]/therm_no[1,3],lwd=4,col=col_b30[24])
# #bline(h=mean(therm[,1]/therm[1,1]),lwd=3,lty=1,col=col_b30[30])
# #abline(h=mean(therm[,2]/therm[1,2]),lwd=3,lty=1,col=col_b30[27])
# #abline(h=mean(therm[,3]/therm[1,3]),lwd=3,lty=1,col=col_b30[24])
# abline(h=sd(therm_no[,1]/therm_no[1,1]),lwd=3,lty=3,col=col_b30[30])
# abline(h=sd(therm_no[,2]/therm_no[1,2]),lwd=3,lty=3,col=col_b30[27])
# abline(h=sd(therm_no[,3]/therm_no[1,3]),lwd=3,lty=3,col=col_b30[24])
# box()
# legend("topleft",c("Gas","Coal","Uranium"),col=c(col_b30[30],col_b30[27],col_b30[24]),lwd=4)
# dev.off()
# 
# # With CO2 Prices 
# 
# # R Plot
# 
# plot(therm[,1]/therm[1,1],type="l",lwd=4,col=col_b30[30],axes=F,main="Evolution of the Thermal Input Prices in Germany (January 1999 - December 2018)",xlab="Time [Years]",ylim=c(min(therm[,1]/therm[1,1],therm[,2]/therm[1,2],therm[,3]/therm[1,3],sd(therm[,1]/therm[1,1]),sd(therm[,2]/therm[1,2]),sd(therm[,3]/therm[1,3])),max(therm[,1]/therm[1,1],therm[,2]/therm[1,2],therm[,3]/therm[1,3])),ylab="Input Prices [EUR/MWh]")
# axis(1,at=seq(1,dim(data)[1],12),labels=c("99","00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18"))
# axis(2)
# lines(therm[,2]/therm[1,2],lwd=4,col=col_b30[27])
# lines(therm[,3]/therm[1,3],lwd=4,col=col_b30[24])
# #bline(h=mean(therm[,1]/therm[1,1]),lwd=3,lty=1,col=col_b30[30])
# #abline(h=mean(therm[,2]/therm[1,2]),lwd=3,lty=1,col=col_b30[27])
# #abline(h=mean(therm[,3]/therm[1,3]),lwd=3,lty=1,col=col_b30[24])
# abline(h=sd(therm[,1]/therm[1,1]),lwd=3,lty=3,col=col_b30[30])
# abline(h=sd(therm[,2]/therm[1,2]),lwd=3,lty=3,col=col_b30[27])
# abline(h=sd(therm[,3]/therm[1,3]),lwd=3,lty=3,col=col_b30[24])
# box()
# legend("topleft",c("Gas","Coal","Uranium"),col=c(col_b30[30],col_b30[27],col_b30[24]),lwd=4)
# 
# # PDF Plot
# 
# pdf("prices_all_normed.pdf",width=14,height=8)
# plot(therm[,1]/therm[1,1],type="l",lwd=4,col=col_b30[30],axes=F,main="Evolution of the Monthly Input Prices in Germany (January 1999 - December 2018)",xlab="Time [Years]",ylim=c(min(therm[,1]/therm[1,1],therm[,2]/therm[1,2],therm[,3]/therm[1,3],sd(therm[,1]/therm[1,1]),sd(therm[,2]/therm[1,2]),sd(therm[,3]/therm[1,3])),max(therm[,1]/therm[1,1],therm[,2]/therm[1,2],therm[,3]/therm[1,3])),ylab="Input Prices [EUR/MWh]",cex.main=2.5,cex.axis=2,cex.lab=2)
# axis(1,at=seq(6,dim(data)[1],12),labels=c("99","00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18"))
# axis(2)
# lines(therm[,2]/therm[1,2],lwd=4,col=col_b30[27])
# lines(therm[,3]/therm[1,3],lwd=4,col=col_b30[24])
# #bline(h=mean(therm[,1]/therm[1,1]),lwd=3,lty=1,col=col_b30[30])
# #abline(h=mean(therm[,2]/therm[1,2]),lwd=3,lty=1,col=col_b30[27])
# #abline(h=mean(therm[,3]/therm[1,3]),lwd=3,lty=1,col=col_b30[24])
# abline(h=sd(therm[,1]/therm[1,1]),lwd=3,lty=3,col=col_b30[30])
# abline(h=sd(therm[,2]/therm[1,2]),lwd=3,lty=3,col=col_b30[27])
# abline(h=sd(therm[,3]/therm[1,3]),lwd=3,lty=3,col=col_b30[24])
# box()
# legend("topleft",c("Gas","Coal","Uranium"),col=c(col_b30[30],col_b30[27],col_b30[24]),lwd=4)
# dev.off()
# 
# # Boxplot(s)
# 
# # No CO2
# 
# # R Plot
# 
# boxplot(gas_no,coal_no,uranium,axes=F,xlab="Different Possibilities",ylab="Price [EUR/MWh]",main="Monthly Fuel Input Prices in Germany (January 1999 - December 2018)",col=c(col_b30[30],col_b30[27],col_b30[24]))
# means_a_no<-c(mean(gas_no),mean(coal_no),mean(uranium))
# points(means_a_no,col="blue",pch=18)
# sds_a_no<-c(sd(gas_no),sd(coal_no),sd(uranium))
# points(sds_a_no,col="lightblue",pch=18)
# axis(1,at=seq(3),labels=c("Gas","Coal","Uranium"))
# axis(2)
# box()
# legend("topright",c("Mean","Standard Deviation"),col=c("blue","lightblue"),lwd=4)
# 
# # PDF Plot
# 
# pdf("box_all_no.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# boxplot(gas_no,coal_no,uranium,axes=F,ylab="Price [EUR/MWh_e]",main="Monthly Fuel Input Prices in Germany (January 1999 - December 2018)",col=c(col_bg[1],col_bg[21],col_bg[41]),cex.main=2.5,cex.axis=2,cex.lab=2)
# means_a_no<-c(mean(gas_no),mean(coal_no),mean(uranium))
# points(means_a_no,col="blue",pch=18)
# sds_a_no<-c(sd(gas_no),sd(coal_no),sd(uranium))
# points(sds_a_no,col="lightblue",pch=18)
# axis(1,at=seq(3),labels=c("Gas","Coal","Uranium"))
# axis(2)
# box()
# legend("topright",c("Mean","Standard Deviation"),col=c("blue","lightblue"),lwd=4)
# dev.off()
# 
# # With CO2
# 
# # R Plot
# 
# boxplot(gas,coal,uranium,axes=F,xlab="Different Possibilities",ylab="Price [EUR/MWh]",main="Monthly Fuel Input Prices in Germany (January 1999 - December 2018)",col=c(col_b30[30],col_b30[27],col_b30[24]))
# means_a<-c(mean(gas),mean(coal),mean(uranium))
# points(means_a,col="blue",pch=18)
# sds_a<-c(sd(gas),sd(coal),sd(uranium))
# points(sds_a,col="lightblue",pch=18)
# axis(1,at=seq(3),labels=c("Gas","Coal","Uranium"))
# axis(2)
# box()
# legend("topright",c("Mean","Standard Deviation"),col=c("blue","lightblue"),lwd=4)
# 
# # PDF Plot
# 
# pdf("box_all.pdf",width=14,height=8)
# par(mar=c(5,5,5,5))
# boxplot(gas,coal,uranium,axes=F,ylab="Price [EUR/MWh_e]",main="Monthly Fuel Input Prices in Germany (January 1999 - December 2018)",col=c(col_bg[1],col_bg[21],col_bg[41]),cex.main=2.5,cex.axis=2,cex.lab=2)
# means_a<-c(mean(gas),mean(coal),mean(uranium))
# points(means_a,col="blue",pch=18)
# sds_a<-c(sd(gas),sd(coal),sd(uranium))
# points(sds_a,col="lightblue",pch=18)
# axis(1,at=seq(3),labels=c("Gas","Coal","Uranium"))
# axis(2)
# box()
# legend("topright",c("Mean","Standard Deviation"),col=c("blue","lightblue"),lwd=4)
# dev.off()

# = Compute the Mean, Standard Deviation, Variance and Correlation of the Thermal Technology Matrix =

mean_monthly_t<-apply(therm,2,mean)
sd_monthly_t<-apply(therm,2,sd)
var_monthly_t<-cov(therm)
corr_monthly_t<-cor(therm)

mean_monthly_co2_t<-apply(therm_co2,2,mean)
sd_monthly_co2_t<-apply(therm_co2,2,sd)
var_monthly_co2_t<-cov(therm_co2)
corr_monthly_co2_t<-cor(therm_co2)

# === Compute the Yearly Thermal Matrix ===

thermal<-matrix(0,20,3)
thermal_co2<-matrix(0,20,3)
for (j in 1:3){
  for (i in 1:20) {
    thermal[i,j]<-mean(therm[12*i-11:11,j])
    thermal_co2[i,j]<-mean(therm_co2[12*i-11:11,j])
  }
}
rownames(thermal)<-c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018")
colnames(thermal)<-c("Gas","Coal","Uranium")
#View(thermal)

rownames(thermal_co2)<-c("1999","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017","2018")
colnames(thermal_co2)<-c("Gas","Coal","Uranium")
#View(thermal)

# plot(thermal[,1],type="l",col="blue",lwd=2.5,ylim=c(1,68))
# lines(thermal[,2],type="l",col="darkblue",lwd=2.5)
# lines(thermal[,3],type="l",col="lightblue",lwd=2.5)

# = Compute the Mean, Standard Deviation, Variance and Correlation of the Thermal Technology Matrix =

mean_t<-apply(thermal,2,mean)
#sd_t<-apply(thermal,2,sd)
var_t<-cov(thermal)
#corr_t<-cor(thermal)

mean_co2_t<-apply(thermal_co2,2,mean)
#sd_co2_t<-apply(thermal_co2,2,sd)
var_co2_t<-cov(thermal_co2)
#corr_co2_t<-cor(thermal_co2)


