# Run for OSMAC Sampling----
library(compiler)
library(here)
library(LaplacesDemon)

enableJIT(1)

# Load the No Correlated Data----
load(here("Non_Identical_r0","Simulation_Setup","Analysis","OSMAC","Init.RData"))

# Load the OSMAC Sample----
load(here("Non_Identical_r0","Simulation_Setup","Analysis","OSMAC","Run_OSMAC.RData"))

print("Real Model")
# Real Model ----
# Generate for Random sample of 1000 different times ---
## Pi=MMSE and Pi=MVC Proportional probabilities ---
Final_Parameter<-Run_OSMAC(Replicates=Replicates,r1=r0,r2=c(100*(2:15)),
                           Y=as.matrix(Simulated_Data$All_Data[,1]),
                           X=as.matrix(Simulated_Data$All_Data[,colnames(Simulated_Data$All_Data) %in%
                                                                 Simulated_Data$All_Models$Real_Model]),
                           Real_Data=as.matrix(Simulated_Data$All_Data[,colnames(Simulated_Data$All_Data) %in%
                                                                         c("Y_Data",Simulated_Data$All_Models$Real_Model)]),
                           N=Simulated_Data$N,Model="Real_Model/OSMAC_output.RData",
                           Theta=Simulated_Data$Theta)

# Assumed Model 1 to 15 ----

for (j in 2:length(Simulated_Data$All_Models)) 
{
  print(paste0("Assumed Model ",j-1))
  # Generate for Random sample of 1000 different times ---
  ## Pi=MMSE and Pi=MVC Proportional probabilities ---
  Final_Parameter<-Run_OSMAC(Replicates=Replicates,r1=r0,r2=c(100*(2:15)),
                             Y=as.matrix(Simulated_Data$All_Data[,1]),
                             X=as.matrix(Simulated_Data$All_Data[,colnames(Simulated_Data$All_Data) %in% 
                                                                   Simulated_Data$All_Models[[j]]]),
                             Real_Data = as.matrix(Simulated_Data$All_Data[,colnames(Simulated_Data$All_Data) %in%
                                                                             c("Y_Data",Simulated_Data$All_Models$Real_Model)]),
                             N=Simulated_Data$N, Model=paste0("Assumed_Model/OSMAC_output_",j-1,".RData"),
                             Theta=Simulated_Data$Theta)
}
