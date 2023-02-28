# Run for Rare Event Random Sampling----
library(compiler)
library(here)
library(LaplacesDemon)

enableJIT(1)

# Load the No Correlated Data----
load(here("Non_Identical_r0","Simulation_Setup","Analysis","RE_Random_Sampling","Init.RData"))

# Load the Rare Event Random Sample----
load(here("Non_Identical_r0","Simulation_Setup","Analysis","RE_Random_Sampling","RE_RandomSample.RData"))

# Real Model ----
# Generate for Rare Evnet Random sample of 1000 different times
Final_Parameter<-Run_RE_RandomSample(Replicates = Replicates, 
                                     FullData = Simulated_Data$All_Data[,colnames(Simulated_Data$All_Data) %in%
                                                                          c("Y_Data",Simulated_Data$All_Models$Real_Model)],
                                     N=Simulated_Data$N,
                                     Subsample_Size=Subsample_Size,Choices=c(100*(2:15)),
                                     No_of_models = length(Simulated_Data$All_Models),
                                     Theta=Simulated_Data$Theta)

Est_Param_RE_RandomSample<-Final_Parameter$EstPar[[1]]
Utility_RE_RandomSample<-Final_Parameter$Utility[[1]]
Bias_RE_RandomSample<-Final_Parameter$Bias[[1]]
#SelectedData_RE_RandomSample<-Final_Parameter$Subsampled_Data[[1]]

# Save the Results ---
save(Est_Param_RE_RandomSample,Utility_RE_RandomSample,Bias_RE_RandomSample,#SelectedData_RE_RandomSample,
     file = here("Non_Identical_r0","Simulation_Setup","Analysis","RE_Random_Sampling",
                 "Results","Real_Model","RE_Random_Sample_output.RData"))

save(Est_Param_RE_RandomSample,Utility_RE_RandomSample,Bias_RE_RandomSample,#SelectedData_RE_RandomSample,
     file = here("Non_Identical_r0","Outputs",Model_Path,"RE_Random_Sampling",
                 "Real_Model","RE_Random_Sample_output.RData"))

# Assumed Model 1----
# Generate for Rare Evnet Random sample of 1000 different times

for(j in 2:length(Simulated_Data$All_Models))
{
  Est_Param_RE_RandomSample<-Final_Parameter$EstPar[[j]]
  Utility_RE_RandomSample<-Final_Parameter$Utility[[j]]
  Bias_RE_RandomSample<-Final_Parameter$Bias[[j]]
  #SelectedData_RE_RandomSample<-Final_Parameter$Subsampled_Data
  
  # Save the Results ---
  save(Est_Param_RE_RandomSample,Utility_RE_RandomSample,Bias_RE_RandomSample,#SelectedData_RE_RandomSample,
       file = here("Non_Identical_r0","Simulation_Setup","Analysis","RE_Random_Sampling",
                   "Results","Assumed_Model",paste0("RE_Random_Sample_output_",j-1,".RData")))
  
  save(Est_Param_RE_RandomSample,Utility_RE_RandomSample,Bias_RE_RandomSample,#SelectedData_RE_RandomSample,
       file = here("Non_Identical_r0","Outputs",Model_Path,"RE_Random_Sampling",
                   "Assumed_Model",paste0("RE_Random_Sample_output_",j-1,".RData")))
}

