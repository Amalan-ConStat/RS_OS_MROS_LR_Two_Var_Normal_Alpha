library(here)
library(compiler)
library(Zelig)
library(LaplacesDemon)
library(psych)

# Using Rare event Random Sampling Sub-sample from Big Data----
Run_RE_RandomSample <- function(Replicates,FullData,Subsample_Size,N,Choices,No_of_models,Theta)
{
  #RandomSample<-list()
  #Sample_RE_RS<-list()
  
  All_Parameter<-list()
  Est_Parameter <- matrix(NA,nrow = length(Choices),ncol = ncol(FullData))

  All_Bias<-list()
  Bias <- matrix(NA,nrow = length(Choices),ncol = ncol(FullData))

  All_optimality<-list()
  optimality <- matrix(NA,nrow = length(Choices),ncol = 3)

  Final_Parameter<-Final_optimality<-Final_Bias<-list()
  
  for (k in 1:No_of_models)
  {
    cat("Model :",k,"\n")
    for(i in 1:Replicates)
    {
      Sampled<-sample(1:N,size = Subsample_Size)
      Temp_Data<-as.data.frame(FullData[Sampled,])
      #RandomSample[[i]]<-Temp_Data
      
      for (j in 1:length(Choices))
      {
        Temp_Temp_Data<-Temp_Data[1:Choices[j],]
        Results <- Zelig::relogit(Y_Data~.-1,data=Temp_Temp_Data,family = binomial(),
                                  tau=mean(FullData[,1]),case.control = "weighting",
                                  bias.correct = TRUE)
        
        Temp_results<-Zelig::relogit(Y_Data~.-1,data=Temp_Temp_Data,family = binomial(),
                                     tau=mean(FullData[,1]),case.control = "weighting",
                                     bias.correct = FALSE)
        
        Est_Parameter[j,]<-c(Choices[j],Results$coefficients)
        Bias[j,]<-c(Choices[j],Temp_results$coefficients-Results$coefficients)
        
        pi_1<-invlogit(as.matrix(Temp_Temp_Data[,-1])%*% Theta)
        W_1<-diag(as.vector(pi_1*(1-pi_1)))
        Mx_1<-t(as.matrix(Temp_Temp_Data[,-1])) %*% W_1 %*% as.matrix(Temp_Temp_Data[,-1])
        
        optimality[j,]<-c(Choices[j],(nrow(Temp_Temp_Data)/(nrow(Temp_Temp_Data)+length(Results$coefficients)))^2*tr(vcov(Results)),
                          det(Mx_1) )
      }
      
      All_Parameter[[i]]<-Est_Parameter
      All_optimality[[i]]<-optimality
      All_Bias[[i]]<-Bias
      #Sample_RE_RS[[i]]<-cbind(i,RandomSample[[i]])
      
      if(i%% 100 ==0){cat(i,":")}
    }
    
    Final_Parameter[[k]]<-do.call(rbind,All_Parameter)
    Final_optimality[[k]]<-do.call(rbind,All_optimality)
    Final_Bias[[k]]<-do.call(rbind,All_Bias)
    #Final_Sample_RE_RS<-do.call(rbind,Sample_RE_RS)
  }
  
  return(list("EstPar"=Final_Parameter,"Utility"=Final_optimality,
              "Bias"=Final_Bias#,"Subsampled_Data"=Final_Sample_RE_RS
              ))
}

Run_RE_RandomSample<-cmpfun(Run_RE_RandomSample)

# Save the Rare Event Random Sample function for No Correlated Covariate Data----
save(Run_RE_RandomSample,
     file=here("Non_Identical_r0","Simulation_Setup","Analysis",
               "RE_Random_Sampling","RE_RandomSample.RData"))

rm(list = ls())

# Run the Rare Event Random sampling method ----
source(here("Non_Identical_r0","Simulation_Setup","Analysis","RE_Random_Sampling",
            "Simulation_Results_RE_Random_Sampling.R"))
