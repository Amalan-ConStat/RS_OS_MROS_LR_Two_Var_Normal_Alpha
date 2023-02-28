# Run for Random Sampling----
library(compiler)
library(here)
library(LaplacesDemon)

enableJIT(1)

# Load the No Correlated Data----
load(here("Non_Identical_r0","Simulation_Setup","Analysis","Random_Sampling","Init.RData"))

# Load the Random Sample----
load(here("Non_Identical_r0","Simulation_Setup","Analysis","Random_Sampling","Run_RandomSample.RData"))

# Real Model ----
# Generate for Random sample of 1000 different times
Final_Parameter<-Run_RandomSample(Replicates = Replicates, 
                                  FullData = Simulated_Data$All_Data[,colnames(Simulated_Data$All_Data) %in%
                                                                       c("Y_Data",Simulated_Data$All_Models$Real_Model)],
                                  N=Simulated_Data$N,
                                  Subsample_Size=Subsample_Size,Choices=c(100*(2:15)),
                                  No_of_Models = length(Simulated_Data$All_Models),
                                  Theta=Simulated_Data$Theta)

Est_Param_RandomSample<-Final_Parameter$EstPar[[1]]
Utility_RandomSample<-Final_Parameter$Utility[[1]]
#SelectedData_RandomSample<-Final_Parameter$Subsampled_Data[[1]]
Bias_RandomSample<-Final_Parameter$Bias[[1]]

# Save the Results ----
save(Est_Param_RandomSample,Utility_RandomSample,Bias_RandomSample,#SelectedData_RandomSample,
     file = here("Non_Identical_r0","Simulation_Setup","Analysis","Random_Sampling",
                 "Results","Real_Model","Random_Sample_output.RData"))

save(Est_Param_RandomSample,Utility_RandomSample,Bias_RandomSample,#SelectedData_RandomSample,
     file = here("Non_Identical_r0","Outputs",Model_Path,"Random_Sampling",
                 "Real_Model","Random_Sample_output.RData"))

# Assumed Models 1 to 15 ----
# Generate for Random sample of 1000 different times

for (j in 2:length(Simulated_Data$All_Models)) 
{
  Est_Param_RandomSample<-Final_Parameter$EstPar[[j]]
  Utility_RandomSample<-Final_Parameter$Utility[[j]]
  #SelectedData_RandomSample<-Final_Parameter$Subsampled_Data[[j]]
  Bias_RandomSample<-Final_Parameter$Bias[[j]]
  
  # Save the Results ---
  save(Est_Param_RandomSample,Utility_RandomSample,Bias_RandomSample,#SelectedData_RandomSample,
       file = here("Non_Identical_r0","Simulation_Setup","Analysis","Random_Sampling",
                   "Results","Assumed_Model",paste0("Random_Sample_output_",j-1,".RData")))
  
  save(Est_Param_RandomSample,Utility_RandomSample,Bias_RandomSample,#SelectedData_RandomSample,
       file = here("Non_Identical_r0","Outputs",Model_Path,"Random_Sampling",
                   "Assumed_Model",paste0("Random_Sample_output_",j-1,".RData")))
}

