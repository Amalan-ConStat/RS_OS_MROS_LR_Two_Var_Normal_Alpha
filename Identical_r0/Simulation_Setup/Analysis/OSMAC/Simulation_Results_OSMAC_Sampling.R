# Run for OSMAC Sampling----
library(compiler)
library(here)
library(LaplacesDemon)

enableJIT(1)

# Load the No Correlated Data----
load(here("Identical_r0","Simulation_Setup","Analysis","OSMAC","Init.RData"))
#load(here("Init.RData"))

# Load the OSMAC Sample----
load(here("Identical_r0","Simulation_Setup","Analysis","OSMAC","Run_OSMAC.RData"))
#load(here("Run_OSMAC.RData"))

# Model 1 and Model 2 ----
# Generate for Random sample of 1000 different times ---
## Pi=MMSE and Pi=MVC Proportional probabilities ---

Final_Parameter<-Run_OSMAC(Replicates = Replicates, r1 = r0, r2 = c(100*(2:15)),
                           Y = as.matrix(Simulated_Data$All_Data[,1]), 
                           X = Simulated_Data$All_Data[,-1],
                           Real_Data = Simulated_Data$All_Data[,colnames(Simulated_Data$All_Data) %in% c("Y_Data",Simulated_Data$All_Models$Real_Model)],
                           N = Simulated_Data$N,
                           alpha = list("Equal"=Equal_Alpha,"UnEqual"=Alpha),
                           combs = Simulated_Data$All_Models[-1],
                           All_Covariates = colnames(Simulated_Data$All_Data)[-1],
                           Theta=Simulated_Data$Theta)
