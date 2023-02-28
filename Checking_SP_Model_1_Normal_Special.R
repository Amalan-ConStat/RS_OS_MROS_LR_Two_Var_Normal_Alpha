Generate_Data<-function(Mean_X,Stdev_X,Theta,N,All_Models)
{
  X<-rmvn(n = N, mu= Mean_X,Sigma=Stdev_X)
  Complete_Data<-cbind(1,X)
  colnames(Complete_Data)<-c(paste0("X",0:length(Mean_X)))
  
  Pi_Data <- invlogit(Complete_Data[,colnames(Complete_Data) %in% All_Models$Real_Model]%*%Theta)
  Y_Data <- rbern(N,Pi_Data)
  
  Outputs<-list()
  for (i in 1:length(All_Models)) 
  {
    All_Data<-cbind(Y=Y_Data,Complete_Data[,colnames(Complete_Data) %in% All_Models[[i]] ])
    Results<-glm(Y~.-1,data=data.frame(All_Data),family="binomial")
    Outputs[[i]]<-coef(Results)
  }
  
  Simulated_Data<-list("N"=N,
                       "Theta"=Theta,
                       "Mean_X"=Mean_X,
                       "Stdev_X"=Stdev_X,
                       "Pi"=Pi_Data,
                       "Outputs"=Outputs,
                       "All_Data"=cbind(Y_Data,Complete_Data),
                       "All_Models"=All_Models)
  
  return(Simulated_Data)
}

N <- 10000
Mean_X<-c(0,0,0,0) 
Stdev_X<-matrix(c(1.5,0,0,0,
                  0,1.5,0,0,
                  0,0,1.5,0,
                  0,0,0,1.5),nrow = 4)
Subsample_Size<-1500; r0<-Nc_size<-100; Replicates <- 1000

# Model 1 ----
# Real Model :$$ log\bigg(\frac{p}{1-p}\bigg) = \beta_0 + \beta_1 X_1 + \beta_2 X_2 $$
No_of_Covariates<-4
All_Models<-list(Real_Model=c("X0","X1","X2"),
                 Assumed_Model_1= c("X0","X1","X2","X3","X4"),
                 Assumed_Model_2= c("X0","X1","X2","X3"),
                 Assumed_Model_3= c("X0","X1","X2","X4"),
                 Assumed_Model_4= c("X0","X1","X3","X4"),
                 Assumed_Model_5= c("X0","X2","X3","X4"),
                 Assumed_Model_6= c("X0","X1","X3"),
                 Assumed_Model_7= c("X0","X1","X4"),
                 Assumed_Model_8= c("X0","X2","X3"),
                 Assumed_Model_9= c("X0","X2","X4"),
                 Assumed_Model_10=c("X0","X3","X4"),
                 Assumed_Model_11=c("X0","X1"),
                 Assumed_Model_12=c("X0","X2"),
                 Assumed_Model_13=c("X0","X3"),
                 Assumed_Model_14=c("X0","X4"),
                 Assumed_Model_15=c("X0"))

# Set parameters
Theta<-c(-1,0.8,0.5)
Simulated_Data<-Generate_Data(Mean_X = Mean_X,Stdev_X = Stdev_X,
                              Theta = Theta,N=N,All_Models =All_Models)
Model_Path<-"Model_1"; 

Equal_Alpha<-rep(1/length(All_Models),length(All_Models))
Alpha<-NULL

for (i in 1:length(All_Models)) 
{
  Alpha[i]<-(5*choose(4,length(Simulated_Data$All_Models[[i]])-1))^(-1)
}

save(list = c("Simulated_Data","Subsample_Size","Replicates","Model_Path","Nc_size","r0","Equal_Alpha","Alpha"),
     file=here("Identical_r0","Generate_Big_Data",Model_Path,"No_Correlated_Covariate.RData"))

save(list = c("Simulated_Data","Subsample_Size","Replicates","Model_Path","Nc_size","r0","Equal_Alpha","Alpha"),
     file=here("Non_Identical_r0","Generate_Big_Data",Model_Path,"No_Correlated_Covariate.RData"))

# Model 2 ----
# Real Model :$$ log\bigg(\frac{p}{1-p}\bigg) = \beta_0 + \beta_1 X_1 + \beta_2 X_2 + \beta_3 X_3$$

No_of_Covariates<-4

All_Models<-list(Real_Model=c("X0","X1","X2","X3"),
                 Assumed_Model_1= c("X0","X1","X2","X3","X4"),
                 Assumed_Model_2= c("X0","X1","X2","X4"),
                 Assumed_Model_3= c("X0","X1","X3","X4"),
                 Assumed_Model_4= c("X0","X2","X3","X4"),
                 Assumed_Model_5= c("X0","X1","X2"),
                 Assumed_Model_6= c("X0","X1","X3"),
                 Assumed_Model_7= c("X0","X1","X4"),
                 Assumed_Model_8= c("X0","X2","X3"),
                 Assumed_Model_9= c("X0","X2","X4"),
                 Assumed_Model_10=c("X0","X3","X4"),
                 Assumed_Model_11=c("X0","X1"),
                 Assumed_Model_12=c("X0","X2"),
                 Assumed_Model_13=c("X0","X3"),
                 Assumed_Model_14=c("X0","X4"),
                 Assumed_Model_15=c("X0"))

# Set parameters
Theta<-c(-1,0.8,0.7,0.5)
Simulated_Data<-Generate_Data(Mean_X = Mean_X,Stdev_X = Stdev_X,
                              Theta = Theta,N=N,All_Models =All_Models)

Equal_Alpha<-rep(1/length(All_Models),length(All_Models))
Alpha<-NULL

for (i in 1:length(All_Models)) 
{
  Alpha[i]<-(5*choose(4,length(Simulated_Data$All_Models[[i]])-1))^(-1)
}

Model_Path<-"Model_2";

save(list = c("Simulated_Data","Subsample_Size","Replicates","Model_Path","Nc_size","r0","Equal_Alpha","Alpha"),
     file=here("Identical_r0","Generate_Big_Data",Model_Path,"No_Correlated_Covariate.RData"))

save(list = c("Simulated_Data","Subsample_Size","Replicates","Model_Path","Nc_size","r0","Equal_Alpha","Alpha"),
     file=here("Non_Identical_r0","Generate_Big_Data",Model_Path,"No_Correlated_Covariate.RData"))

rm(list = ls())
