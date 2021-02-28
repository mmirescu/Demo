# ===== Technical and Costs Parameters =====

setwd("C:/Users/Mirescu1/ucloud/Universities/TU Wien/Theses/Dissertation/Eigene/Latex")

# === Install the necessary packages ===

library(xts)
library(quadprog)
library(rootSolve)
library(MASS)

# === Read in the Necessary Data ===

parameters<-read.csv("data_parameters_old.csv",header=T,sep=",",dec=".")
#parameters<-read.csv("data_parameters.csv",header=T,sep=",",dec=".")


# === Capital Recovery Factor ===

r<-0.0642
L_t<-parameters[1:(dim(parameters)[1]-1),dim(parameters)[2]-3]
L_r<-parameters[dim(parameters)[1],dim(parameters)[2]-3]

CRF_t<-sapply(L_t, function(L) r/(1-(1+r)^(-L)))
CRF_r<-sapply(L_r, function(L) r/(1-(1+r)^(-L)))

# === Fixed Costs ===

# Investment Costs [Euro per MW]

IC_t<-parameters[1:(dim(parameters)[1]-1),dim(parameters)[2]-2]
IC_r<-parameters[dim(parameters)[1],dim(parameters)[2]-2]
IC_r_ic<-0.8*IC_r

# Fixed Operations and Maintenance Costs [Euro per MW per Year ]

FOM_t<-parameters[1:(dim(parameters)[1]-1),dim(parameters)[2]-1]
FOM_r<-parameters[dim(parameters)[1],dim(parameters)[2]-1]

# Overall Fixed Costs 

F_t<-(IC_t*CRF_t+FOM_t)/8760
#F_t<-F_t[1]
F_r<-(IC_r*CRF_r+FOM_r)/8760
F_r_ic<-(IC_r_ic*CRF_r+FOM_r)/8760

# === Variable Costs ===

VOM_t<-parameters[1:(dim(parameters)[1]-1),dim(parameters)[2]]
#VOM_t<-VOM_t[1]
VOM_r<-parameters[(dim(parameters)[1]),dim(parameters)[2]]
