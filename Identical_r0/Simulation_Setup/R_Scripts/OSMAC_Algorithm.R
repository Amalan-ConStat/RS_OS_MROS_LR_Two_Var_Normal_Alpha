# Cordeiro Bias Estimation ----
Cordeiro<-function(XData,With_bias)
{
  p <- as.vector(invlogit(XData%*%as.vector(With_bias)))
  W <- diag(p*(1-p))
  inverse_term <- solve(t(XData)%*%W%*%XData)
  
  Term1 <- inverse_term%*%t(XData)%*%W
  Term2 <- diag(diag(XData%*%(inverse_term)%*%t(XData))) %*% as.matrix(p-0.5)
  
  bias <- as.vector(Term1%*%Term2)
  return(bias)
} 

# Newton Rhapson Algorithm for MLE ----
getMLE <- function(x, y, w) {
  beta <- rep(0, ncol(x))
  loop  <- 1
  Loop  <- 200
  msg <- "NA"
  while (loop <= Loop) {
    pr <- c(1 - 1 / (1 + exp(x %*% beta)))
    H <- t(x) %*% (pr * (1 - pr) * w * x)
    S <- colSums((y - pr) * w * x)
    tryCatch(
      {shs <- NA
      shs <- solve(H, S) }, 
      error=function(e){
        cat("\n ERROR :", loop, conditionMessage(e), "\n")})
    if (is.na(shs[1])) {
      msg <- "Not converge"
      beta <- loop <- NA
      break
    }
    beta.new <- beta + shs
    tlr  <- sum((beta.new - beta)^2)
    beta  <- beta.new
    if(tlr < 0.000001) {
      msg <- "Successful convergence"
      break
    }
    if (loop == Loop)
      warning("Maximum iteration reached")
    loop  <- loop + 1
  }
  list(par=beta, message=msg, iter=loop)
}

# Two step OSMAC ----
AlgTwoStp <- function(r1=r1, r2=r2,Y,X,n,Real_Data,alpha,combs,All_Covariates,Theta){
    Y_Real<-Real_Data[,1] #  Real Data
    X_Real<-Real_Data[,-1] # Real Data
    
    n1 <- sum(Y)
    n0 <- n - n1
    PI.prop <- rep(1/(2*n0), n)
    PI.prop[Y==1] <- 1/(2*n1)
    idx.prop <- sample(1:n, r1, T, PI.prop)
    
    x.prop<-lapply(1:length(combs),function(a){
      as.matrix(X[idx.prop,All_Covariates %in% combs[[a]]]) # Assumed Data 
    })
    y.prop <- Y[idx.prop,]  # Assumed Data 
    
    x_Real.prop<-X_Real[idx.prop,] # Real Data
    y_Real.prop<-Y_Real[idx.prop] # Real Data
    
    pinv.prop <- n
    pinv.prop <- 1/PI.prop[idx.prop]
    
    fit.prop <- lapply(1:length(combs), function(a){
      getMLE(x=x.prop[[a]], y=y.prop, w=pinv.prop) # Assumed Data
    })
    fit_Real.prop <- getMLE(x=x_Real.prop, y=y_Real.prop, w=pinv.prop) # Real Data
    
    beta.prop<-list()
    for (j in 1:length(combs)) 
      {
      beta.prop[[j]] <- fit.prop[[j]]$par # Assumed Data
      if(anyNA(beta.prop[[j]]))
        {
        return(list(opt=NA, msg="first stage not converge"))
        }
      }
    beta_Real.prop <- fit_Real.prop$par # Real Data
    
    if(anyNA(beta_Real.prop[1]))
      return(list(opt=NA, msg="first stage not converge"))
    
    P.prop  <- lapply(1:length(combs),function(a){
      1 - 1 / (1 + exp(as.matrix(X[,All_Covariates %in% combs[[a]] ]) %*% beta.prop[[a]])) # Assumed Data
    })
    P_Real.prop  <- 1 - 1 / (1 + exp(X_Real %*% beta_Real.prop)) # Real Data
    
    # For the Real Model 
    beta.mVc_Real<-matrix(nrow = length(r2),ncol = ncol(X_Real) )
    Utility_mVc_Real<-matrix(nrow = length(r2),ncol = 3 )
    Bias_mVc_Real<-matrix(nrow = length(r2),ncol = ncol(X_Real) )
    
    beta.mMSE_Real<-matrix(nrow = length(r2),ncol = ncol(X_Real) )
    Utility_mMSE_Real<-matrix(nrow = length(r2),ncol = 3 )
    Bias_mMSE_Real<-matrix(nrow = length(r2),ncol = ncol(X_Real) )
    
    #Sample.mMSE_Real<-list()
    #Sample.mVc_Real<-list()
    
    beta.mVc_Old<-Utility_mVc_Old<-Bias_mVc_Old<-list()
    beta.mMSE_Old<-Utility_mMSE_Old<-Bias_mMSE_Old<-list()
    #Sample.mMSE_Assumed_Old<-list()
    #Sample.mVc_Assumed_Old<-list()
    
    # For the Assumed Model with Already Available Sub-sampling probabilities
    for (a in 1:length(combs)) 
    {
      beta.mVc_Old[[a]]<-matrix(nrow = length(r2),ncol = ncol(X_Real) )
      Utility_mVc_Old[[a]]<-matrix(nrow = length(r2),ncol = 3 )
      Bias_mVc_Old[[a]]<-matrix(nrow = length(r2),ncol = ncol(X_Real) )

      beta.mMSE_Old[[a]]<-matrix(nrow = length(r2),ncol = ncol(X_Real) )
      Utility_mMSE_Old[[a]]<-matrix(nrow = length(r2),ncol = 3 )
      Bias_mMSE_Old[[a]]<-matrix(nrow = length(r2),ncol = ncol(X_Real) )
      
      #Sample.mMSE_Assumed_Old[[a]]<-list()
      #Sample.mVc_Assumed_Old[[a]]<-list()
    }
    
    # For the Real Model with joined Sub-sampling probabilities Equal
    beta.mVc_join<-matrix(nrow = length(r2),ncol = ncol(X_Real) )
    Utility_mVc_join<-matrix(nrow = length(r2),ncol = 3 )
    Bias_mVc_join<-matrix(nrow = length(r2),ncol = ncol(X_Real) )
    
    beta.mMSE_join<-matrix(nrow = length(r2),ncol = ncol(X_Real) )
    Utility_mMSE_join<-matrix(nrow = length(r2),ncol = 3 )
    Bias_mMSE_join<-matrix(nrow = length(r2),ncol = ncol(X_Real) )
    
    #Sample.mMSE_join<-list()
    #Sample.mVc_join<-list()
    
    # For the Real Model with joined Sub-sampling probabilities Un-Equal
    beta.mVc_join_e<-matrix(nrow = length(r2),ncol = ncol(X_Real) )
    Utility_mVc_join_e<-matrix(nrow = length(r2),ncol = 3 )
    Bias_mVc_join_e<-matrix(nrow = length(r2),ncol = ncol(X_Real) )
    
    beta.mMSE_join_e<-matrix(nrow = length(r2),ncol = ncol(X_Real) )
    Utility_mMSE_join_e<-matrix(nrow = length(r2),ncol = 3 )
    Bias_mMSE_join_e<-matrix(nrow = length(r2),ncol = ncol(X_Real) )
    
    #Sample.mMSE_join_e<-list()
    #Sample.mVc_join_e<-list()
    
    ## mVc
    PI_Real.mVc <- sqrt((Y - P_Real.prop)^2 * rowSums(X_Real^2)) # Real Data
    PI_Real.mVc <- PI_Real.mVc / sum(PI_Real.mVc) # Real Data
    
    PI_Assumed_Old.mVc <- lapply(1:length(combs), function(a){
      sqrt((Y - P.prop[[a]])^2 * rowSums(as.matrix(X[,All_Covariates %in% combs[[a]] ])^2)) # Assumed Data Old probabilities
    } )
    
    PI_Assumed_Old.mVc <- lapply(1:length(combs), function(a){
      PI_Assumed_Old.mVc[[a]] / sum(PI_Assumed_Old.mVc[[a]]) # Assumed Data Old probabilities
    })
    
    PI_join.mVc<-rowSums2(cbind(alpha$Equal[1]*PI_Real.mVc,
                                do.call(cbind,PI_Assumed_Old.mVc)%*%diag(alpha$Equal[-1])) )
    
    PI_join.mVc_e<-rowSums2(cbind(alpha$UnEqual[1]*PI_Real.mVc,
                                  do.call(cbind,PI_Assumed_Old.mVc)%*%diag(alpha$UnEqual[-1])) )
    
    ## mMSE
    p_Real.prop <- P_Real.prop[idx.prop] # Real data
    w_Real.prop <- p_Real.prop * (1 - p_Real.prop) # Real data
    W_Real.prop <- solve(t(x_Real.prop) %*% (x_Real.prop * w_Real.prop * pinv.prop)) # Real data
    
    p_Assumed.prop <- lapply(1:length(combs),function(a){
      P.prop[[a]][idx.prop] # Assumed data
    })
    w_Assumed.prop <- lapply(1:length(combs),function(a){
      p_Assumed.prop[[a]] * (1 - p_Assumed.prop[[a]]) # Assumed data
    })
    W_Assumed.prop <- lapply(1:length(combs),function(a){
      solve(t(x.prop[[a]]) %*% (x.prop[[a]] * w_Assumed.prop[[a]] * pinv.prop)) # Assumed data
    })
    
    PI_Real.mMSE <- sqrt((Y_Real - P_Real.prop)^2 * rowSums((X_Real%*%W_Real.prop)^2)) # Real data
    PI_Real.mMSE <- PI_Real.mMSE / sum(PI_Real.mMSE) # Real data
    
    PI_Assumed_Old.mMSE <- lapply(1:length(combs),function(a){
      sqrt((Y - P.prop[[a]])^2 * rowSums((as.matrix(X[,All_Covariates %in% combs[[a]]])%*%W_Assumed.prop[[a]])^2)) # Assumed data
    })
    
    PI_Assumed_Old.mMSE <- lapply(1:length(combs),function(a){
      PI_Assumed_Old.mMSE[[a]] / sum(PI_Assumed_Old.mMSE[[a]]) # Assumed data
    })
    
    PI_join.mMSE<-rowSums2(cbind(alpha$Equal[1]*PI_Real.mMSE,(do.call(cbind,PI_Assumed_Old.mMSE)%*%diag(alpha$Equal[-1]))))
    PI_join.mMSE_e<-rowSums2(cbind(alpha$UnEqual[1]*PI_Real.mMSE,(do.call(cbind,PI_Assumed_Old.mMSE)%*%diag(alpha$UnEqual[-1]))))
    
    for (i in 1:length(r2)) 
    {
      # mVc
      idx_Real.mVc <- sample(1:n, r2[i]-r1, T, PI_Real.mVc) # Real Data 
      idx_Assumed.mVc <- lapply(1:length(combs), function(a){
        sample(1:n, r2[i]-r1, T, PI_Assumed_Old.mVc[[a]]) # Assumed Data Old probabilities
      })
      idx_join.mVc <- sample(1:n, r2[i]-r1, T, PI_join.mVc) # Joined Real Data Equal
      idx_join.mVc_e <- sample(1:n, r2[i]-r1, T, PI_join.mVc_e) # Joined Real Data Un Equal
      
      x_Real.mVc <- X_Real[c(idx_Real.mVc, idx.prop),] # Real Data
      y_Real.mVc <- Y_Real[c(idx_Real.mVc, idx.prop)] # Real Data
      
      x_Assumed.mVc <-lapply(1:length(combs),function(a){
        as.matrix(X_Real[c(idx_Assumed.mVc[[a]], idx.prop),]) # Assumed Data
      }) 
        
      y_Assumed.mVc <- lapply(1:length(combs),function(a){
        Y_Real[c(idx_Assumed.mVc[[a]], idx.prop)] # Assumed Data
      })
        
      x_join.mVc <- X_Real[c(idx_join.mVc, idx.prop),] # Joined Data Equal 
      y_join.mVc <- Y_Real[c(idx_join.mVc, idx.prop)] # Joined Data Equal
      
      x_join.mVc_e <- X_Real[c(idx_join.mVc_e, idx.prop),] # Joined Data Un-Equal
      y_join.mVc_e <- Y_Real[c(idx_join.mVc_e, idx.prop)] # Joined Data Un-Equal
      
      fit_Real.mVc <- getMLE(x=x_Real.mVc, y=y_Real.mVc, 
                             w=c(1 / PI_Real.mVc[idx_Real.mVc], pinv.prop)) # Real Data
      
      fit_Assumed_Old.mVc <-lapply(1:length(combs), function(a){
        getMLE(x=x_Assumed.mVc[[a]], y=y_Assumed.mVc[[a]], 
               w=c(1 / PI_Assumed_Old.mVc[[a]][idx_Assumed.mVc[[a]]], pinv.prop)) # Assumed Data Old probabilities
      }) 
      
      fit_join.mVc <- getMLE(x=x_join.mVc, y=y_join.mVc, 
                             w=c(1 / PI_join.mVc[idx_join.mVc], pinv.prop)) # Joined Data Equal
      
      fit_join.mVc_e <- getMLE(x=x_join.mVc_e, y=y_join.mVc_e, 
                               w=c(1 / PI_join.mVc_e[idx_join.mVc_e], pinv.prop)) # Joined Data Un-Equal
      
      # Sample.mVc_Real[[i]]<-cbind(r2[i],y_Real.mVc,x_Real.mVc,
      #                             c(PI_Real.mVc[idx_Real.mVc], 1 / pinv.prop)) # Real Data
      
      # for (j in 1:length(combs)) 
      # {
      #   Sample.mVc_Assumed_Old[[j]][[i]]<-cbind(r2[i],y_Assumed.mVc[[j]],x_Assumed.mVc[[j]],
      #                                           c(PI_Assumed_Old.mVc[[j]][idx_Assumed.mVc[[j]] ], 1 / pinv.prop)) 
      # Assumed Data Old probabilities
      # 
      # }
      # 
      # Sample.mVc_join[[i]]<-cbind(r2[i],y_join.mVc,x_join.mVc,
      #                             c(PI_join.mVc[idx_join.mVc], 1 / pinv.prop)) # Joined Data
      
      if(anyNA(fit_Real.mVc$par) || anyNA(fit_Assumed_Old.mVc$par) || anyNA(fit_join.mVc$par) || anyNA(fit_join.mVc_e$par))
      {
        stop("There are NA or NaN values")
      }
        
      beta.mVc_Real[i,] <- fit_Real.mVc$par # Real Data
      
      for (j in 1:length(combs)) 
      {
        beta.mVc_Old[[j]][i,] <- fit_Assumed_Old.mVc[[j]]$par # Assumed Data Old probabilities 
      }
      
      beta.mVc_join[i,] <- fit_join.mVc$par # Joined Data Equal
      beta.mVc_join_e[i,] <- fit_join.mVc_e$par # Joined Data Un Equal
      
      # Real Data
      pi<- invlogit(x_Real.mVc %*% beta.mVc_Real[i,])
      W<-diag(as.vector(pi*(1-pi)*c(1 / PI_Real.mVc[idx_Real.mVc], pinv.prop)))
      Mx<-(t(x_Real.mVc) %*% W %*% x_Real.mVc)
      Mx<-solve(Mx)
      Middle<-diag(((as.vector(y_Real.mVc)-as.vector(pi))*
                      as.vector(c(1 / PI_Real.mVc[idx_Real.mVc], pinv.prop)))^2)
      V_Temp<-(t(x_Real.mVc)%*%Middle%*%x_Real.mVc)
      V_Final<-Mx %*% V_Temp %*% Mx
      
      pi_1<- invlogit(x_Real.mVc %*% Theta)
      W_1<-diag(as.vector(pi_1*(1-pi_1)))
      Mx_1<-(t(x_Real.mVc) %*% W_1 %*% x_Real.mVc)
      
      Utility_mVc_Real[i,]<-cbind(r2[i],tr(V_Final),det(Mx_1))
      Bias_mVc_Real[i,]<-Cordeiro(XData=x_Real.mVc,With_bias = beta.mVc_Real[i,])
      
      # Assumed Data Old probabilities
      pi<- lapply(1:length(combs), function(a){
        invlogit(x_Assumed.mVc[[a]] %*% beta.mVc_Old[[a]][i,])
      })
        
      W<-lapply(1:length(combs),function(a){
        diag(as.vector(pi[[a]]*(1-pi[[a]])*c(1 / PI_Assumed_Old.mVc[[a]][idx_Assumed.mVc[[a]]], pinv.prop)))
      })
        
      Mx<-lapply(1:length(combs),function(a){
        solve((t(x_Assumed.mVc[[a]]) %*% W[[a]] %*% x_Assumed.mVc[[a]]))
      }) 
      Middle<-lapply(1:length(combs),function(a){
        diag(((as.vector(y_Assumed.mVc[[a]])-as.vector(pi[[a]]))*
                as.vector(c(1 / PI_Assumed_Old.mVc[[a]][idx_Assumed.mVc[[a]]], pinv.prop)))^2)
      }) 
      V_Temp<-lapply(1:length(combs), function(a){
        (t(x_Assumed.mVc[[a]])%*%Middle[[a]]%*%x_Assumed.mVc[[a]])
      })
        
      V_Final<-lapply(1:length(combs),function(a){
        Mx[[a]] %*% V_Temp[[a]] %*% Mx[[a]]
      }) 
      
      pi_1<- lapply(1:length(combs), function(a){
        invlogit(x_Assumed.mVc[[a]] %*% Theta)
      })
      
      W_1<-lapply(1:length(combs),function(a){
        diag(as.vector(pi_1[[a]]*(1-pi_1[[a]])))
      })
      
      Mx_1<-lapply(1:length(combs),function(a){
        (t(x_Assumed.mVc[[a]]) %*% W_1[[a]] %*% x_Assumed.mVc[[a]])
      })
      
      for (j in 1:length(combs)) 
      {
        Utility_mVc_Old[[j]][i,]<-cbind(r2[i],tr(V_Final[[j]]),det(Mx_1[[j]]))
        Bias_mVc_Old[[j]][i,]<-Cordeiro(XData=x_Assumed.mVc[[j]],With_bias = beta.mVc_Old[[j]][i,])    
      }
      
      # Assumed Data Joined Probabilities
      pi<- invlogit(x_join.mVc %*% beta.mVc_join[i,])
      W<-diag(as.vector(pi*(1-pi)*c(1 / PI_join.mVc[idx_join.mVc], pinv.prop)))
      Mx<-(t(x_join.mVc) %*% W %*% x_join.mVc)
      Mx<-solve(Mx)
      Middle<-diag(((as.vector(y_join.mVc)-as.vector(pi))*
                      as.vector(c(1 / PI_join.mVc[idx_join.mVc], pinv.prop)))^2)
      V_Temp<-(t(x_join.mVc)%*%Middle%*%x_join.mVc)
      V_Final<-Mx %*% V_Temp %*% Mx
      
      pi_1<- invlogit(x_join.mVc %*% Theta)
      W_1<-diag(as.vector(pi_1*(1-pi_1)))
      Mx_1<-(t(x_join.mVc) %*% W_1 %*% x_join.mVc)
      
      Utility_mVc_join[i,]<-cbind(r2[i],tr(V_Final),det(Mx_1))
      Bias_mVc_join[i,]<-Cordeiro(XData=x_join.mVc,With_bias = beta.mVc_join[i,])
      
      # Assumed Data Joined Probabilities Un Equal 
      pi<- invlogit(x_join.mVc_e %*% beta.mVc_join_e[i,])
      W<-diag(as.vector(pi*(1-pi)*c(1 / PI_join.mVc_e[idx_join.mVc_e], pinv.prop)))
      Mx<-(t(x_join.mVc_e) %*% W %*% x_join.mVc_e)
      Mx<-solve(Mx)
      Middle<-diag(((as.vector(y_join.mVc_e)-as.vector(pi))*
                      as.vector(c(1 / PI_join.mVc_e[idx_join.mVc_e], pinv.prop)))^2)
      V_Temp<-(t(x_join.mVc_e)%*%Middle%*%x_join.mVc_e)
      V_Final<-Mx %*% V_Temp %*% Mx
      
      pi_1<- invlogit(x_join.mVc_e %*% Theta)
      W_1<-diag(as.vector(pi_1*(1-pi_1)))
      Mx_1<-(t(x_join.mVc_e) %*% W_1 %*% x_join.mVc_e)
      
      Utility_mVc_join_e[i,]<-cbind(r2[i],tr(V_Final),det(Mx_1))
      Bias_mVc_join_e[i,]<-Cordeiro(XData=x_join.mVc_e,With_bias = beta.mVc_join_e[i,])    
      
      ## mMSE
      idx_Real.mMSE <- sample(1:n, r2[i]-r1, T, PI_Real.mMSE) # Real data
      idx_Assumed.mMSE <- lapply(1:length(combs),function(a){
        sample(1:n, r2[i]-r1, T, PI_Assumed_Old.mMSE[[a]]) # Assumed data
      }) 
      idx_join.mMSE <- sample(1:n, r2[i]-r1, T, PI_join.mMSE) # Joined Data Equal 
      idx_join.mMSE_e <- sample(1:n, r2[i]-r1, T, PI_join.mMSE_e) # Joined Data Un-Equal
      
      x_Real.mMSE <- X_Real[c(idx_Real.mMSE, idx.prop),] # Real Data
      y_Real.mMSE <- Y_Real[c(idx_Real.mMSE, idx.prop)] # Real Data
      
      x_Assumed.mMSE <- lapply(1:length(combs),function(a){
        as.matrix(X_Real[c(idx_Assumed.mMSE[[a]], idx.prop),]) # Assumed Data
      })
      y_Assumed.mMSE <- lapply(1:length(combs),function(a){
        Y_Real[c(idx_Assumed.mMSE[[a]], idx.prop)] # Assumed Data
      })
      
      x_join.mMSE <- X_Real[c(idx_join.mMSE, idx.prop),] # Joined Data Equal
      y_join.mMSE <- Y_Real[c(idx_join.mMSE, idx.prop)] # Joined Data Equal
      
      x_join.mMSE_e <- X_Real[c(idx_join.mMSE_e, idx.prop),] # Joined Data Un-Equal
      y_join.mMSE_e <- Y_Real[c(idx_join.mMSE_e, idx.prop)] # Joined Data Un-Equal
      
      fit_Real.mMSE <- getMLE(x=x_Real.mMSE, y=y_Real.mMSE, 
                              w=c(1 / PI_Real.mMSE[idx_Real.mMSE], pinv.prop)) # Real Data
      fit_Assumed_Old.mMSE <- lapply(1:length(combs),function(a){
        getMLE(x=x_Assumed.mMSE[[a]], y=y_Assumed.mMSE[[a]], 
               w=c(1 / PI_Assumed_Old.mMSE[[a]][idx_Assumed.mMSE[[a]]], pinv.prop)) # Assumed Data Old probabilities
      })
        
      fit_join.mMSE <- getMLE(x=x_join.mMSE, y=y_join.mMSE, 
                             w=c(1 / PI_join.mMSE[idx_join.mMSE], pinv.prop)) # Joined Data Equal
      fit_join.mMSE_e <- getMLE(x=x_join.mMSE_e, y=y_join.mMSE_e, 
                                w=c(1 / PI_join.mMSE_e[idx_join.mMSE_e], pinv.prop)) # Joined Data Un-Equal
      
      # Sample.mMSE_Real[[i]]<-cbind(r2[i],y_Real.mMSE,x_Real.mMSE,
      #                             c(PI_Real.mMSE[idx_Real.mMSE], 1 / pinv.prop)) # Real Data
      
      # for (j in 1:length(combs))
      # {
      #   Sample.mMSE_Assumed_Old[[j]][[i]]<-cbind(r2[i],y_Assumed.mMSE[[j]],x_Assumed.mMSE[[j]],
      #                                            c(PI_Assumed_Old.mMSE[[j]][idx_Assumed.mMSE[[j]]], 1 / pinv.prop)) 
      # Assumed Data Old probabilities
      #   
      # }
      # 
      # Sample.mMSE_join[[i]]<-cbind(r2[i],y_join.mMSE,x_join.mMSE,
      #                             c(PI_join.mMSE[idx_join.mMSE], 1 / pinv.prop)) # Joined Data
      
      if(anyNA(fit_Real.mMSE$par) || anyNA(fit_Assumed_Old.mMSE$par) || anyNA(fit_join.mMSE$par) || anyNA(fit_join.mMSE_e$par))
      {
        stop("There are NA or NaN values")
      }
      
      beta.mMSE_Real[i,] <- fit_Real.mMSE$par # Real Data
      for (j in 1:length(combs)) 
      {
        beta.mMSE_Old[[j]][i,] <- fit_Assumed_Old.mMSE[[j]]$par # Assumed Data Old probabilities 
      }
      beta.mMSE_join[i,] <- fit_join.mMSE$par # Joined Data Equal
      beta.mMSE_join_e[i,] <- fit_join.mMSE_e$par # Joined Data Un-Equal
      
      # Real Data
      pi<- invlogit(x_Real.mMSE %*% beta.mMSE_Real[i,])
      W<-diag(as.vector(pi*(1-pi)*c(1 / PI_Real.mMSE[idx_Real.mMSE], pinv.prop)))
      Mx<-(t(x_Real.mMSE) %*% W %*% x_Real.mMSE)
      Mx<-solve(Mx)
      Middle<-diag(((as.vector(y_Real.mMSE)-as.vector(pi))*
                      as.vector(c(1 / PI_Real.mMSE[idx_Real.mMSE], pinv.prop)))^2)
      V_Temp<-(t(x_Real.mMSE)%*%Middle%*%x_Real.mMSE)
      V_Final<-Mx %*% V_Temp %*% Mx
      
      pi_1<- invlogit(x_Real.mMSE %*% Theta)
      W_1<-diag(as.vector(pi_1*(1-pi_1)))
      Mx_1<-(t(x_Real.mMSE) %*% W_1 %*% x_Real.mMSE)
      
      Utility_mMSE_Real[i,]<-cbind(r2[i],tr(V_Final),det(Mx_1))
      Bias_mMSE_Real[i,]<-Cordeiro(XData=x_Real.mMSE,With_bias = beta.mMSE_Real[i,])
      
      # Assumed Data Old probabilities
      pi<-lapply(1:length(combs),function(a){
        invlogit(x_Assumed.mMSE[[a]] %*% beta.mMSE_Old[[a]][i,])
      }) 
      W<-lapply(1:length(combs),function(a){
        diag(as.vector(pi[[a]]*(1-pi[[a]])*c(1 / PI_Assumed_Old.mMSE[[a]][idx_Assumed.mMSE[[a]]], pinv.prop)))
      }) 
      Mx<-lapply(1:length(combs),function(a){
        solve((t(x_Assumed.mMSE[[a]]) %*% W[[a]] %*% x_Assumed.mMSE[[a]]))
      }) 
      Middle<-lapply(1:length(combs),function(a){
        diag(((as.vector(y_Assumed.mMSE[[a]])-as.vector(pi[[a]]))*
                as.vector(c(1 / PI_Assumed_Old.mMSE[[a]][idx_Assumed.mMSE[[a]]], pinv.prop)))^2)
      })
      V_Temp<-lapply(1:length(combs),function(a){
        (t(x_Assumed.mMSE[[a]])%*%Middle[[a]]%*%x_Assumed.mMSE[[a]])
        }) 
      V_Final<-lapply(1:length(combs),function(a){
        Mx[[a]] %*% V_Temp[[a]] %*% Mx[[a]]
      }) 
      
      pi_1<-lapply(1:length(combs),function(a){
        invlogit(x_Assumed.mMSE[[a]] %*% Theta)
      }) 
      W_1<-lapply(1:length(combs),function(a){
        diag(as.vector(pi_1[[a]]*(1-pi_1[[a]])))
      }) 
      Mx_1<-lapply(1:length(combs),function(a){
        (t(x_Assumed.mMSE[[a]]) %*% W_1[[a]] %*% x_Assumed.mMSE[[a]])
      })
      
      for (j in 1:length(combs)) 
      {
        Utility_mMSE_Old[[j]][i,]<-cbind(r2[i],tr(V_Final[[j]]),det(Mx_1[[j]]))
        Bias_mMSE_Old[[j]][i,]<-Cordeiro(XData=x_Assumed.mMSE[[j]],With_bias = beta.mMSE_Old[[j]][i,])    
      }
          
      # Assumed Data Joined Probabilities Equal
      pi<- invlogit(x_join.mMSE %*% beta.mMSE_join[i,])
      W<-diag(as.vector(pi*(1-pi)*c(1 / PI_join.mMSE[idx_join.mMSE], pinv.prop)))
      Mx<-(t(x_join.mMSE) %*% W %*% x_join.mMSE)
      Mx<-solve(Mx)
      Middle<-diag(((as.vector(y_join.mMSE)-as.vector(pi))*
                      as.vector(c(1 / PI_join.mMSE[idx_join.mMSE], pinv.prop)))^2)
      V_Temp<-(t(x_join.mMSE)%*%Middle%*%x_join.mMSE)
      V_Final<-Mx %*% V_Temp %*% Mx
      
      pi_1<- invlogit(x_join.mMSE %*% Theta)
      W_1<-diag(as.vector(pi_1*(1-pi_1)))
      Mx_1<-(t(x_join.mMSE) %*% W_1 %*% x_join.mMSE)
      
      Utility_mMSE_join[i,]<-cbind(r2[i],tr(V_Final),det(Mx_1))
      Bias_mMSE_join[i,]<-Cordeiro(XData=x_join.mMSE,With_bias = beta.mMSE_join[i,])    
      
      # Assumed Data Joined Probabilities Un-Equal
      pi<- invlogit(x_join.mMSE_e %*% beta.mMSE_join_e[i,])
      W<-diag(as.vector(pi*(1-pi)*c(1 / PI_join.mMSE_e[idx_join.mMSE_e], pinv.prop)))
      Mx<-(t(x_join.mMSE_e) %*% W %*% x_join.mMSE_e)
      Mx<-solve(Mx)
      Middle<-diag(((as.vector(y_join.mMSE_e)-as.vector(pi))*
                      as.vector(c(1 / PI_join.mMSE_e[idx_join.mMSE_e], pinv.prop)))^2)
      V_Temp<-(t(x_join.mMSE_e)%*%Middle%*%x_join.mMSE_e)
      V_Final<-Mx %*% V_Temp %*% Mx
      
      pi_1<- invlogit(x_join.mMSE_e %*% Theta)
      W_1<-diag(as.vector(pi_1*(1-pi_1)))
      Mx_1<-(t(x_join.mMSE_e) %*% W_1 %*% x_join.mMSE_e)
      
      Utility_mMSE_join_e[i,]<-cbind(r2[i],tr(V_Final),det(Mx_1))
      Bias_mMSE_join_e[i,]<-Cordeiro(XData=x_join.mMSE_e,With_bias = beta.mMSE_join_e[i,])    
      
    }
    
    # Full_SP_Real<-cbind(Real_Data,PI_Real.mMSE,PI_Real.mVc)
    # Full_SP_Old<-cbind(Real_Data,do.call(cbind,PI_Assumed_Old.mMSE),do.call(cbind,PI_Assumed_Old.mVc))
    # Full_SP_join<-cbind(Real_Data,PI_join.mMSE,PI_join.mVc)
    # 
    # Sample_mVc_Real<-do.call(rbind,Sample.mVc_Real)
    # for (j in 1:length(combs)) 
    # {
    #   assign(paste0("Sample_mVc_Old_",j),
    #          do.call(rbind,Sample.mVc_Assumed_Old[[j]]))
    # }
    # Sample_mVc_join<-do.call(rbind,Sample.mVc_join)
    # 
    # Sample_mMSE_Real<-do.call(rbind,Sample.mMSE_Real)
    # for (j in 1:length(combs)) 
    # {
    #   assign(paste0("Sample_mMSE_Old_",j),
    #          do.call(rbind,Sample.mMSE_Assumed_Old[[j]]))
    # }
    # Sample_mMSE_join<-do.call(rbind,Sample.mMSE_join)
    
    opt_Real<-list("Est_Theta_mMSE"=cbind(r2,beta.mMSE_Real),"Utility_mMSE"=Utility_mMSE_Real,
                   "Bias_mMSE"=cbind(r2,Bias_mMSE_Real),
                   "Est_Theta_mVc"=cbind(r2,beta.mVc_Real),"Utility_mVc"=Utility_mVc_Real,
                   "Bias_mVc"=cbind(r2,Bias_mVc_Real))
    
    for(j in 1:length(combs))
    {
      assign(paste0("opt_Old_",j),
             list("Est_Theta_mMSE"=cbind(r2,beta.mMSE_Old[[j]]),"Utility_mMSE"=Utility_mMSE_Old[[j]],
                  "Bias_mMSE"=cbind(r2,Bias_mMSE_Old[[j]]),
                  "Est_Theta_mVc"=cbind(r2,beta.mVc_Old[[j]]),"Utility_mVc"=Utility_mVc_Old[[j]],
                  "Bias_mVc"=cbind(r2,Bias_mVc_Old[[j]]))
             )
    }
    
    opt_join<-list("Est_Theta_mMSE"=cbind(r2,beta.mMSE_join),"Utility_mMSE"=Utility_mMSE_join,
                   "Bias_mMSE"=cbind(r2,Bias_mMSE_join),
                   "Est_Theta_mVc"=cbind(r2,beta.mVc_join),"Utility_mVc"=Utility_mVc_join,
                   "Bias_mVc"=cbind(r2,Bias_mVc_join))
    
    opt_join_e<-list("Est_Theta_mMSE"=cbind(r2,beta.mMSE_join_e),"Utility_mMSE"=Utility_mMSE_join_e,
                     "Bias_mMSE"=cbind(r2,Bias_mMSE_join_e),
                     "Est_Theta_mVc"=cbind(r2,beta.mVc_join_e),"Utility_mVc"=Utility_mVc_join_e,
                     "Bias_mVc"=cbind(r2,Bias_mVc_join_e))
    
    # msg <- c(fit.prop$message, fit_Real.prop$message, 
    #          fit_Real.mVc$message,fit_Assumed_Old.mVc$message,fit_Assumed_New.mVc$message,fit_join.mVc$message,
    #          fit_Real.mMSE$message,fit_Assumed_Old.mMSE$message,fit_Assumed_New.mMSE$message,fit_join.mMSE$message)
    msg <- NULL
    return(list(opt=list("opt_Real"=opt_Real,
                         "opt_Old"=mget(paste0("opt_Old_",1:length(combs))),
                         "opt_join_Equal"=opt_join,"opt_join_UnEqual"=opt_join_e), msg=msg)) 
}
