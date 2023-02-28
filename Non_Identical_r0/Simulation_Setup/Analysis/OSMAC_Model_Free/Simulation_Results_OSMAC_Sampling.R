# Run for OSMAC Sampling----
library(compiler)
library(here)
library(LaplacesDemon)

enableJIT(1)

# Load the No Correlated Data----
load(here("Non_Identical_r0","Simulation_Setup","Analysis","OSMAC_Model_Free","Init.RData"))

# Load the OSMAC Sample----
load(here("Non_Identical_r0","Simulation_Setup","Analysis","OSMAC_Model_Free","Run_OSMAC.RData"))
 
# Generate for Random sample of 1000 different times ----
## Pi=MMSE and Pi=MVC Proportional probabilities ---
Final_Parameter<-Run_OSMAC_MF(Replicates=Replicates,r1=r0,r2=c(100*(2:15)),
                              Y=as.matrix(Simulated_Data$All_Data[,1]),
                              X=as.matrix(Simulated_Data$All_Data[,colnames(Simulated_Data$All_Data) %in%
                                                                    paste0("X",0:4)]),
                              alpha = list("Equal"=Equal_Alpha,"UnEqual"=Alpha), 
                              N=Simulated_Data$N,
                              All_Covariates=paste0("X",0:4),
                              combs=Simulated_Data$All_Models,
                              Theta=Simulated_Data$Theta)
