##########################################################################
##########################################################################
###    R Code for pooling biomarkers(1-n and 2-n), no interaction term
###
###                 Yujie Wu, Abigail Sloan, Molin Wang
###
###             For any questions please contact Yujie Wu
###                  Email: yujiewu@hsph.harvard.edu
##########################################################################
##########################################################################

##########################################################################
##########################################################################
### Contents:
###   I.   Library, data sourcing, working directory
###   II.  Sourced functions
###   III. Example syntax to run code
###
##########################################################################
##########################################################################



#################################################################################
#################################################################################
###   I. Library and data sourcing

library("survival")

setwd("~/")

mydata = readRDS("yourdata.rds")

#################################################################################
#################################################################################

#################################################################################
#################################################################################
###   II. Functions
###      a. expit
###      b. fc_df:         full calibration with additional covariates
###      c. fc_df_no_cov:  internalized calibration without adidtional covariates


expit = function(x){return(exp(x)/(exp(x)+1))}

fc_df = function(data, X, S, H, W, Y, strata, nstud, nref, covar=NA){
  
  #### 1. Create matrix to store output
  output_fc = matrix(NA, ncol=3, nrow=1)
  colnames(output_fc) = c("Point estimate of betax","Estimated RR", "Estimated variance of betahatx")
  
  #### 2. Rename variables names given to function to ensure consistency
  data$X = data[[X]]
  data$S = data[[S]]
  data$H = data[[H]]
  data$W = data[[W]]
  data$Y = data[[Y]]
  data$strata = data[[strata]]
  
  #### 2b. Sort data frame, in case it wasn't already
  data = data[with(data, order(S, strata, Y)),]
  
  strata2           = c()
  ncase             = c()
  ncont             = c()
  count             = 1
  for( study in 1:nstud){
    data.temp       = subset(data, S==study)
    stratum_in_S    = max(data.temp$strata)
    for( str in 1:stratum_in_S){
      data.temp2    = subset(data.temp, strata==str)
      strata2       = append(strata2, rep(count, nrow(data.temp2)))
      ncase         = append(ncase, rep(sum(data.temp2$Y), nrow(data.temp2)))
      ncont         = append(ncont, rep((nrow(data.temp2)-sum(data.temp2$Y)), nrow(data.temp2)))
      count         = count+1
    }
  }
  data$strata2      = strata2
  data$ncase        = ncase
  data$ncont        = ncont
  
  #### 2b. Sort data frame, in case it wasn't already
  data = data[with(data, order(S, strata, Y)),]
  
  #### 3. Compute other useful quantities
  nz = length(covar)
  np = max(strata2) # total number of strata
  n1 = c(rep(NA,nstud)) # number of matched pairs in each study
  nc = nstud - nref # number of studies needing calibration
  for(k in 1:nstud){
    data.temp = subset(data, S==k)
    n1[k]     = max(data.temp$strata)
  }
  
  #### 4. Complete calibration studies and add appropriate ahat, bhat to data frame
  a_hat = c(rep(0,nref),rep(NA,nc))
  b_hat = c(rep(1,nref),rep(NA,nc))
  
  for(k in (nref+1):nstud){
    cal_data_s  = subset(data, S==k & H==1)
    fit         = lm(X~W, data=cal_data_s)
    a_hat[k]    = fit$coefficients[1]
    b_hat[k]    = fit$coefficients[2]
  }
  
  data$a_hat = a_hat[data$S] ## adding a_hat and b_hat to the dataframe
  data$b_hat = b_hat[data$S]
  
  #### 5. Create xhat_fc variable- use H==2 to indicate when using X ref lab
  data$xhat_fc  = ifelse(data$H==2, data$X, data$a_hat + data$b_hat*data$W)
  
  #### 6. Obtain point estimate from standard logistic regression
  formula         = as.formula(paste("Y~xhat_fc+strata(strata2)",paste(covar,collapse ="+"),sep="+"))
  fc_fit          = clogit(formula, data=data)
  beta_hat        = fc_fit$coefficients[1] # betahat_x value
  betazhat        = fc_fit$coefficients[-1] # extract betazhat values
  output_fc[1]   = beta_hat
  output_fc[2]   = exp(beta_hat)
  
  #### 6.b. Compute Xd for each observation, Zd, and BX (linear combo)
  
  matrix.var <- c()
  for(str2 in 1:np){
    data_s_j           = subset(data,strata2 == str2)
    if(data_s_j$ncase[1] == 1){
      ncont              = data_s_j$ncont[1]
      total_num          = choose(ncont+1,1)
      pairs              = combn(seq(1,(ncont+1),1),m=1)[,-total_num] #### all possible combinations in each strata, excluding two cases being chosen
      xd_fc              = matrix(NA, nrow = length(pairs), ncol = 1)
      Wd                 = matrix(NA, nrow = length(pairs), ncol = 1)
      Zd                 = matrix(NA, nrow = length(pairs), ncol = length(covar))
      S                  = matrix(data_s_j$S[1], nrow = length(pairs), ncol = 1)
      strata             = matrix(data_s_j$strata[1], nrow = length(pairs), ncol = 1)
      strata2            = matrix(str2, nrow = length(pairs), ncol = 1)
      for(i in 1:(total_num-1)){
        i1               = pairs[i]
        xd_fc[i]         = data_s_j[i1,]$xhat_fc - data_s_j[(ncont+1),]$xhat_fc
        Wd[i]            = data_s_j[i1,]$W - data_s_j[(ncont+1),]$W 
        col_nums         = which(names(data) %in% covar)
        for( j in 1:length(covar)){
          col_num        = col_nums[j]
          Zd[i,j]        = data_s_j[i1,col_num] - data_s_j[(ncont+1),col_num] 
        }
      }
      BXd                = beta_hat*(xd_fc) + (Zd) %*%betazhat
      matrix.temp        = cbind(xd_fc, Zd, Wd, BXd, S, strata, strata2)
      matrix.var         = rbind(matrix.var, matrix.temp)
    } 
    if(data_s_j$ncase[1] == 2){
      ncont              = data_s_j$ncont[1]
      total_num          = choose(ncont+2,2)
      pairs              = combn(seq(1,(ncont+2),1),m=2)[,-total_num] #### all possible combinations in each strata, excluding two cases being chosen
      xd_fc              = matrix(NA, nrow = ncol(pairs), ncol = 1)
      Wd                 = matrix(NA, nrow = ncol(pairs), ncol = 1)
      Zd                 = matrix(NA, nrow = ncol(pairs), ncol = length(covar))
      S                  = matrix(data_s_j$S[1], nrow = ncol(pairs), ncol = 1)
      strata             = matrix(data_s_j$strata[1], nrow = ncol(pairs), ncol = 1)
      strata2            = matrix(str2, nrow = ncol(pairs), ncol = 1)
      for(i in 1:(total_num-1)){
        i1               = pairs[,i][1]
        i2               = pairs[,i][2]
        xd_fc[i]        = (data_s_j[i1,]$xhat_fc + data_s_j[i2,]$xhat_fc) - 
                           (data_s_j[(ncont+1),]$xhat_fc + data_s_j[(ncont+2),]$xhat_fc)
        Wd[i]            = (data_s_j[i1,]$W + data_s_j[i2,]$W) - (data_s_j[(ncont+1),]$W + data_s_j[(ncont+2),]$W)
        col_nums         = which(names(data) %in% covar)
        for( j in 1:length(covar)){
          col_num        = col_nums[j]
          Zd[i,j]        = (data_s_j[i1,col_num] + data_s_j[i2,col_num]) - (data_s_j[(ncont+1),col_num] + data_s_j[(ncont+2),col_num])
        }
      }
      BXd                = beta_hat*(xd_fc) + (Zd) %*%betazhat
      matrix.temp        = cbind(xd_fc, Zd, Wd, BXd, S, strata, strata2)
      matrix.var         = rbind(matrix.var, matrix.temp)
    } 
  }
  covard               = paste(covar,"d",sep = "")
  colnames(matrix.var) = c("xd_fc",covard,"Wd","BXd","S","strata","strata2")
  data.var             = data.frame(matrix.var)
  
  #### 7. Compute variance; prepare matrices 
  dim_sand  = 2*nstud+1+nz
  A         = matrix(0, ncol = dim_sand, nrow = dim_sand)
  B         = matrix(0, ncol = dim_sand, nrow = dim_sand)
  
  ###################### A MATRIX #######################################
  
  #### 7.a.i: The upper block diagonals of A
  for(k in (nref+1):nstud){
    cal_data_s               = subset(data, (H==1 & S==k)) # Specific cal data with only controls
    A[2*k-1, 2*k-1]          = (1/np)*sum((cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
    A[2*k-1,(2*k-1)+1]       = (1/np)*sum(cal_data_s$W*(cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
    A[(2*k-1)+1,(2*k-1)]     = A[2*k-1, (2*k-1)+1]
    A[(2*k-1)+1,(2*k-1)+1]   = (1/np)*sum((cal_data_s$W^2)*(cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
  }
  
  #### 7.a.ii. Upper right entries of A (A12)
  A12    = matrix(NA, ncol=(nz+1), nrow=(2*nstud))
  
  for(k in 1:nstud){
    
    ### Betax and a,b ###
    ### we have to use for loop to calculate psi1,psi2,psiR for each strata ###
    data_s_h     = subset(data, S==k) # the kth study
    var_s_h      = subset(data.var, S==k)
    psi1         = rep(NA, n1[k])
    psi2         = rep(NA, n1[k])
    psiR         = rep(NA, n1[k])
    for(str in 1:n1[k]){
      data_s_h.temp     = subset(data_s_h, strata==str)
      data_s_h.temp2    = subset(data_s_h.temp, H==1)
      psi1[str]         = sum(data_s_h.temp2$X-a_hat[k]-b_hat[k]*data_s_h.temp2$W) # estimating equ 1
      psi2[str]         = sum((data_s_h.temp2$X-a_hat[k]-b_hat[k]*data_s_h.temp2$W)*data_s_h.temp2$W) # estimating equ 2
      
      var_s_h.temp      = subset(var_s_h, strata==str)
      Xd_s_h            = var_s_h.temp$xd_fc
      BXd_s_h           = var_s_h.temp$BXd
      psiR[str]         = -sum(Xd_s_h*exp(BXd_s_h))/(sum(exp(BXd_s_h))+1) # last estimating equation
    }
    
    A12[(2*k-1),1]      = (1/np)*sum(psi1*psiR)
    A12[(2*k),1]        = (1/np)*sum(psi2*psiR)
    
    ### Betaz and a,b ###
    ## Grab the Z covariate vector of interest; do computations with that particular Z
    for(i in 1:nz){ 
      col_num           = which(names(var_s_h) == covard[i])
      psiZ              = rep(NA, n1[k])
      for(str in 1:n1[k]){
        var_s_h.temp    = subset(var_s_h, strata==str)
        ZZ              = var_s_h.temp[,col_num]
        BXd_s_h         = var_s_h.temp$BXd
        psiZ[str]       = -sum(ZZ*exp(BXd_s_h))/(1+sum(exp(BXd_s_h)))
      }
      A12[(2*k-1),(1+i)]  = (1/np)*sum(psi1*psiZ)
      A12[(2*k),(1+i)]    = (1/np)*sum(psi2*psiZ) 
    }
    
  }
  
  
  A[1:(2*nstud), (2*nstud+1):dim_sand] = A12 # Fill appropriate piece of A matrix
  A21                                  = t(A12) # Reflect entries for A21
  A[(2*nstud+1):dim_sand,1:(2*nstud)]  = A21
  
  
  ## Now compute A22 entries
  A22    = matrix(NA, ncol=(nz+1), nrow=(nz+1))
  
  ## Betax^2 entry
  psiX            = rep(NA, np)
  for(str2 in 1:np){
    var_s_h       = subset(data.var, strata2==str2)
    Xd            = var_s_h$xd_fc
    BXd           = var_s_h$BXd
    psiX[str2]    = -(sum(Xd*exp(BXd)))/(sum(exp(BXd))+1) #psiX for each strata in each study
  }
  
  A22[1,1]        = (1/np)*sum(psiX^2)
  
  ### Computations in A22 involving betax and betaz, betaz squared
  
  for(i in 1:nz){ 
    col_num         = which(names(data.var) == covard[i])
    psiZ            = rep(NA, np)
    for(str2 in 1:np){
      var_s_h      = subset(data.var, strata2==str2)
      ZZ            = var_s_h[,col_num]
      BXd           = var_s_h$BXd
      psiZ[str2]    = -(sum(ZZ*exp(BXd)))/(sum(exp(BXd))+1) #psiX for each strata in each study
    }
    A22[(1),(1+i)]      = (1/np)*sum(psiX*psiZ)
    A22[(1+i),(1)]      = (1/np)*sum(psiX*psiZ)
    A22[(1+i),(1+i)]    = (1/np)*sum(psiZ*psiZ) 
  }
  
  ### Computations involving cross terms of Z
  
  ## Interact the covariate psi terms (Z1 and Z2...) for all pairs
  
  pairs = combn(seq(1,nz,1),m=2)
  for(j in 1:choose(nz,2)){ #for each pairwise combination
    pair = pairs[,j]
    
    # Obtain appropriate matrix column vector based on pairwise selection
    col_num1 = which(names(data.var) == covard[pair[1]])
    col_num2 = which(names(data.var) == covard[pair[2]])
    
    psi_betaz1            = rep(NA, np)
    psi_betaz2            = rep(NA, np)
    for(str2 in 1:np){
      var_s_h             = subset(data.var, strata2==str2)
      ZZ1                 = var_s_h[,col_num1]
      ZZ2                 = var_s_h[,col_num2]
      BXd                 = var_s_h$BXd
      psi_betaz1[str2]    = -(sum(ZZ1*exp(BXd)))/(sum(exp(BXd))+1) #psi_betaz for each strata in each study
      psi_betaz2[str2]    = -(sum(ZZ2*exp(BXd)))/(sum(exp(BXd))+1) 
    }
    
    #Fill matrix elements
    A22[(1+pair[2]),(1+pair[1])] = (1/np)*sum(psi_betaz1*psi_betaz2)
    A22[(1+pair[1]),(1+pair[2])] = (1/np)*sum(psi_betaz1*psi_betaz2)
  }
  
  
  ## Now place A22 in the main A matrix
  A[(2*nstud+1):dim_sand, (2*nstud+1):dim_sand] = A22
  
  ## Eliminate extra entries from A not associated with calibration studies
  A = A[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  #################### B Matrix ############################  
  
  ### 8.a.i. B11 entries: derivatives associated with calibration studies
  for(k in (nref+1):nstud){
    cal_data_s               = subset(data, S==k & H==1) # Specific cal data with only controls
    B[2*k-1,2*k-1]           = -1
    B[2*k-1,(2*k-1)+1]       = (-1/np)*sum(cal_data_s$W)
    B[(2*k-1)+1,(2*k-1)]     = B[2*k-1,(2*k-1)+1]
    B[(2*k-1)+1,(2*k-1)+1]   = (-1/np)*sum(cal_data_s$W^2)
  }
  
  ### 8.b. B12 entries: 
  
  B12 = matrix(0, ncol=(1+nz), nrow=(2*nstud))
  
  for(k in (nref+1):nstud){
    var_s_j               = subset(data.var, S==k)
    # Odd entries (a and betax), internalized only
    B12[(2*k-1), 1]       = 0
    
    # a and betaz, which exists for internalized only
    
    for(i in 1:nz){ 
      B12[(2*k-1),(1+i)]  = 0
    }
    
    # Even entries (b and betax)
    
    frac1                     = rep(NA, n1[k])
    for( str in 1:n1[k]){
      #### for observations outside calibration study
      var_s_x                 = subset(var_s_j, strata == str)
      BXd                     = var_s_x$BXd #get all BX_d
      Xd                      = var_s_x$xd_fc
      
      Wd                      = var_s_x$Wd  #W_d 
      numerator1              = (sum(Wd*exp(BXd))+
                                   sum(Xd*beta_hat*Wd*exp(BXd)))*(1+sum(exp(BXd)))
      numerator2              = sum(Xd*exp(BXd))*sum(beta_hat*Wd*exp(BXd))
      denominator             = (1+sum(exp(BXd)))^2
      frac1[str]              = -(numerator1-numerator2)/denominator
    }
    ## Sum all components
    B12[(2*k), 1]   = (1/np)*(sum(frac1))
    
    
    ### Betaz and b : full calibration and internalized pairs
    
    for(i in 1:nz){ 
      col_num                   = which(names(data.var) == covard[i])
      frac1                     = rep(NA, n1[k])
      
      for( str in 1:n1[k]){
        #### for observations outside calibration study
        var_s_x                 = subset(var_s_j, strata == str)
        BXd                     = var_s_x$BXd #get all BX_d
        Zd                      = var_s_x[,col_num]
        
        Wd                      = var_s_x$Wd  #W_d 
        numerator1              = sum(Zd*beta_hat*Wd*exp(BXd))*(1+sum(exp(BXd)))
        numerator2              = sum(Zd*exp(BXd))*sum(beta_hat*Wd*exp(BXd))
        denominator             = (1+sum(exp(BXd)))^2
        frac1[str]              = -(numerator1-numerator2)/denominator
      }
      B12[(2*k),(1+i)]  = (1/np)*(sum(frac1))
    }
  }
  
  #### Replace B12 and its transpose in the B matrix
  B[1:(2*nstud),(2*nstud+1):dim_sand] = B12
  B[(2*nstud+1):dim_sand,1:(2*nstud)] = t(B12)
  B[1:(2*nstud),(2*nstud+1):dim_sand] = 0
  ### B22 entry : Completed component-wise #######
  
  B22    = matrix(NA, ncol=(nz+1), nrow=(nz+1))
  
  ## Betax^2 entry
  frac <- rep(NA, np)
  for(str2  in 1:np){
    var_s_j      = subset(data.var, strata2 == str2)
    Xd_s_j       = var_s_j$xd_fc
    BXd_s_j      = var_s_j$BXd
    numerator    = sum((Xd_s_j)^2*exp(BXd_s_j))*(1+sum(exp(BXd_s_j)))-(sum(Xd_s_j*exp(BXd_s_j)))^2
    denominator  = (1+sum(exp(BXd_s_j)))^2
    frac[str2]   = -numerator/denominator
  }
  B22[1,1]  = (1/np)*sum(frac)
  
  ### Computations in B22 involving betax and betaz, betaz squared
  
  for(i in 1:nz){ 
    col_num  = which(names(data.var) == covard[i])
    
    fracZZ    = rep(NA,np)
    fracXZ    = rep(NA,np)
    for(str2 in 1:np){
      var_s_j       = subset(data.var, strata2 == str2)
      Xd_s_j        = var_s_j$xd_fc
      Zd_s_j        = var_s_j[,col_num]
      BXd_s_j       = var_s_j$BXd
      numerator1    = sum(Xd_s_j*Zd_s_j*exp(BXd_s_j))*(1+sum(exp(BXd_s_j)))-
        (sum(Xd_s_j*exp(BXd_s_j)))*(sum(Zd_s_j*exp(BXd_s_j)))
      denominator1  = (1+sum(exp(BXd_s_j)))^2
      fracXZ[str2]  = -numerator1/denominator1
      
      numerator2    = sum(Zd_s_j^2*exp(BXd_s_j))*(1+sum(exp(BXd_s_j)))-(sum(Zd_s_j*exp(BXd_s_j)))^2
      denominator2  = denominator1
      fracZZ[str2]  = -numerator2/denominator2
    }
    
    B22[(1),(1+i)]      = (1/np)*sum(fracXZ)
    B22[(1+i),(1)]      = (1/np)*sum(fracXZ)
    B22[(1+i),(1+i)]    = (1/np)*sum(fracZZ) 
  }
  
  ### Computations involving cross terms of Z
  
  ## Interact the covariate psi terms (Z1 and Z2...) for all pairs
  
  pairs = combn(seq(1,nz,1),m=2)
  for(j in 1:choose(nz,2)){ #for each pairwise combination
    pair = pairs[,j]
    
    # Obtain appropriate matrix column vector based on pairwise selection
    col_num1 = which(names(data.var) == covard[pair[1]])
    
    col_num2 = which(names(data.var) == covard[pair[2]])
    
    fracZZ         = rep(NA, np)
    for(str2 in 1:np){
      var_s_j      = subset(data.var, strata2 == str2)
      ZZ1          = var_s_j[,col_num1]
      ZZ2          = var_s_j[,col_num2]
      BXd_s_j      = var_s_j$BXd
      
      numerator    = sum(ZZ1*ZZ2*exp(BXd_s_j))*(1+sum(exp(BXd_s_j)))-
        sum(ZZ1*exp(BXd_s_j))*sum(ZZ2*exp(BXd_s_j))
      denominator  = (1+sum(exp(BXd_s_j)))^2
      
      fracZZ[str2] = -numerator/denominator
    }
    
    #Fill matrix elements
    B22[(1+pair[2]),(1+pair[1])] = (1/np)*sum(fracZZ)
    B22[(1+pair[1]),(1+pair[2])] = (1/np)*sum(fracZZ)
  }
  
  
  ## Now place B22 in the main B matrix
  B[(2*nstud+1):dim_sand, (2*nstud+1):dim_sand] = B22
  
  ## Eliminate extra entries from B associated with noncalibration studies
  B = B[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  ### Compute sandwich variance estimator and place variance estimate in output
  V                = solve(B)%*%A%*%t(solve(B))
  output_fc[3]    = (1/np)*V[(2*nc+1),(2*nc+1)]
  
  #### 99. Return appropriate output
  return(output_fc)
  
}


fc_df_no_cov = function(data, X, S, H, W, Y, strata, nstud, nref){
  
  #### 1. Create matrix to store output
  output_fc = matrix(NA, ncol=3, nrow=1)
  colnames(output_fc) = c("Point estimate of betax","Estimated RR", "Estimated variance of betahatx")
  
  #### 2. Rename variables names given to function to ensure consistency
  data$X = data[[X]]
  data$S = data[[S]]
  data$H = data[[H]]
  data$W = data[[W]]
  data$Y = data[[Y]]
  data$strata = data[[strata]]
  
  #### 2b. sort the data in case it is not ready
  data = data[with(data, order(S, strata, Y)),]
  
  strata2           = c()
  ncase             = c()
  ncont             = c()
  count             = 1
  for( study in 1:nstud){
    data.temp       = subset(data, S==study)
    stratum_in_S    = max(data.temp$strata)
    for( str in 1:stratum_in_S){
      data.temp2    = subset(data.temp, strata==str)
      strata2       = append(strata2, rep(count, nrow(data.temp2)))
      ncase         = append(ncase, rep(sum(data.temp2$Y), nrow(data.temp2)))
      ncont         = append(ncont, rep((nrow(data.temp2)-sum(data.temp2$Y)), nrow(data.temp2)))
      count         = count+1
    }
  }
  data$strata2      = strata2
  data$ncase        = ncase
  data$ncont        = ncont
  
  #### 2b. Sort data frame, in case it wasn't already
  data = data[with(data, order(S, strata, Y)),]
  
  #### 3. Compute other useful quantities
  np = max(strata2) # total number of strata
  n1 = c(rep(NA,nstud)) # number of matched pairs in each study
  nc = nstud - nref # number of studies needing calibration
  for(k in 1:nstud){
    data.temp = subset(data, S==k)
    n1[k]     = max(data.temp$strata)
  }
  
  #### 4. Complete calibration studies and add appropriate ahat, bhat to data frame
  a_hat = c(rep(0,nref),rep(NA,nc))
  b_hat = c(rep(1,nref),rep(NA,nc))
  
  for(k in (nref+1):nstud){
    cal_data_s  = subset(data, S==k & H==1)
    fit         = lm(X~W, data=cal_data_s)
    a_hat[k]    = fit$coefficients[1]
    b_hat[k]    = fit$coefficients[2]
  }
  
  data$a_hat = a_hat[data$S] ## adding a_hat and b_hat to the dataframe
  data$b_hat = b_hat[data$S]
  
  #### 5. Create xhat_fc variable- use H==2 to indicate when using X ref lab
  data$xhat_fc  = ifelse(data$H==2, data$X, data$a_hat + data$b_hat*data$W)

  
  #### 6. Obtain point estimate from standard logistic regression
  formula         = as.formula(paste("Y~xhat_fc+strata(strata2)",sep="+"))
  fc_fit          = clogit(formula, data=data)
  beta_hat        = fc_fit$coefficients[1] # betahat_x value
  betazhat        = fc_fit$coefficients[-1] # extract betazhat values
  output_fc[1]   = beta_hat
  output_fc[2]   = exp(beta_hat)
  
  #### 6.b. Compute Xd for each observation, Zd, and BX (linear combo)
  
  matrix.var <- c()
  for(str2 in 1:np){
    data_s_j           = subset(data,strata2 == str2)
    if(data_s_j$ncase[1] == 1){
      ncont              = data_s_j$ncont[1]
      total_num          = choose(ncont+1,1)
      pairs              = combn(seq(1,(ncont+1),1),m=1)[,-total_num] #### all possible combinations in each strata, excluding two cases being chosen
      xd_fc              = matrix(NA, nrow = length(pairs), ncol = 1)
      Wd                 = matrix(NA, nrow = length(pairs), ncol = 1)
      S                  = matrix(data_s_j$S[1], nrow = length(pairs), ncol = 1)
      strata             = matrix(data_s_j$strata[1], nrow = length(pairs), ncol = 1)
      strata2            = matrix(str2, nrow = length(pairs), ncol = 1)
      for(i in 1:(total_num-1)){
        i1               = pairs[i]
        xd_fc[i]         = data_s_j[i1,]$xhat_fc - data_s_j[(ncont+1),]$xhat_fc
        Wd[i]            = data_s_j[i1,]$W - data_s_j[(ncont+1),]$W 
      }
      BXd                = beta_hat*(xd_fc) 
      matrix.temp        = cbind(xd_fc, Wd, BXd, S, strata, strata2)
      matrix.var         = rbind(matrix.var, matrix.temp)
    } 
    if(data_s_j$ncase[1] == 2){
      ncont              = data_s_j$ncont[1]
      total_num          = choose(ncont+2,2)
      pairs              = combn(seq(1,(ncont+2),1),m=2)[,-total_num] #### all possible combinations in each strata, excluding two cases being chosen
      xd_fc              = matrix(NA, nrow = ncol(pairs), ncol = 1)
      Wd                 = matrix(NA, nrow = ncol(pairs), ncol = 1)
      S                  = matrix(data_s_j$S[1], nrow = ncol(pairs), ncol = 1)
      strata             = matrix(data_s_j$strata[1], nrow = ncol(pairs), ncol = 1)
      strata2            = matrix(str2, nrow = ncol(pairs), ncol = 1)
      for(i in 1:(total_num-1)){
        i1               = pairs[,i][1]
        i2               = pairs[,i][2]
        xd_fc[i]        = (data_s_j[i1,]$xhat_fc + data_s_j[i2,]$xhat_fc) - 
          (data_s_j[(ncont+1),]$xhat_fc + data_s_j[(ncont+2),]$xhat_fc)
        Wd[i]            = (data_s_j[i1,]$W + data_s_j[i2,]$W) - (data_s_j[(ncont+1),]$W + data_s_j[(ncont+2),]$W)
      }
      BXd                = beta_hat*(xd_fc) 
      matrix.temp        = cbind(xd_fc, Wd, BXd, S, strata, strata2)
      matrix.var         = rbind(matrix.var, matrix.temp)
    } 
  }
  colnames(matrix.var) = c("xd_fc","Wd","BXd","S","strata","strata2")
  data.var             = data.frame(matrix.var)
  
  #### 7. Compute variance; prepare matrices 
  dim_sand  = 2*nstud+1
  A         = matrix(0, ncol = dim_sand, nrow = dim_sand)
  B         = matrix(0, ncol = dim_sand, nrow = dim_sand)
  
  ###################### A MATRIX #######################################
  
  #### 7.a.i: The upper block diagonals of A
  for(k in (nref+1):nstud){
    cal_data_s               = subset(data, (H==1 & S==k)) # Specific cal data with only controls
    A[2*k-1, 2*k-1]          = (1/np)*sum((cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
    A[2*k-1,(2*k-1)+1]       = (1/np)*sum(cal_data_s$W*(cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
    A[(2*k-1)+1,(2*k-1)]     = A[2*k-1, (2*k-1)+1]
    A[(2*k-1)+1,(2*k-1)+1]   = (1/np)*sum((cal_data_s$W^2)*(cal_data_s$X-a_hat[k]-b_hat[k]*cal_data_s$W)^2)
  }
  
  #### 7.a.ii. Upper right entries of A (A12)
  A12    = matrix(NA, ncol=1, nrow=(2*nstud))
  
  for(k in 1:nstud){
    
    ### Betax and a,b ###
    ### we have to use for loop to calculate psi1,psi2,psiR for each strata ###
    data_s_h     = subset(data, S==k) # the kth study
    var_s_h      = subset(data.var, S==k)
    psi1         = rep(NA, n1[k])
    psi2         = rep(NA, n1[k])
    psiR         = rep(NA, n1[k])
    for(str in 1:n1[k]){
      data_s_h.temp     = subset(data_s_h, strata==str)
      data_s_h.temp2    = subset(data_s_h.temp, H==1)
      psi1[str]         = sum(data_s_h.temp2$X-a_hat[k]-b_hat[k]*data_s_h.temp2$W) # estimating equ 1
      psi2[str]         = sum((data_s_h.temp2$X-a_hat[k]-b_hat[k]*data_s_h.temp2$W)*data_s_h.temp2$W) # estimating equ 2
      
      var_s_h.temp      = subset(var_s_h, strata==str)
      Xd_s_h            = var_s_h.temp$xd_fc
      BXd_s_h           = var_s_h.temp$BXd
      psiR[str]         = -sum(Xd_s_h*exp(BXd_s_h))/(sum(exp(BXd_s_h))+1) # last estimating equation
    }
    A12[(2*k-1),1]      = (1/np)*sum(psi1*psiR)
    A12[(2*k),1]        = (1/np)*sum(psi2*psiR)
  }
  
  
  A[1:(dim_sand-1), dim_sand] = A12 # Fill appropriate piece of A matrix
  A21                         = t(A12) # Reflect entries for A21
  A[dim_sand,1:(dim_sand-1)]  = A21
  
  
  ## Betax^2 entry
  psiX            = rep(NA, np)
  for(str2 in 1:np){
    var_s_h       = subset(data.var, strata2==str2)
    Xd            = var_s_h$xd_fc
    BXd           = var_s_h$BXd
    psiX[str2]    = -(sum(Xd*exp(BXd)))/(sum(exp(BXd))+1) #psiX for each strata in each study
  }
  
  A[dim_sand,dim_sand] = (1/np)*sum( psiX^2 )
  
  
  ## Eliminate extra entries from A associated with calibration studies
  A = A[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  #################### B Matrix ############################  
  
  ### 8.a.i. B11 entries: derivatives associated with calibration studies
  for(k in (nref+1):nstud){
    cal_data_s               = subset(data, S==k & H==1) # Specific cal data with only controls
    B[2*k-1,2*k-1]           = -1
    B[2*k-1,(2*k-1)+1]       = (-1/np)*sum(cal_data_s$W)
    B[(2*k-1)+1,(2*k-1)]     = B[2*k-1,(2*k-1)+1]
    B[(2*k-1)+1,(2*k-1)+1]   = (-1/np)*sum(cal_data_s$W^2)
  }
  
  ### 8.b. B21 and B12 entries: place directly in B matrix
  
  B12 = matrix(0, ncol=1, nrow=(2*nstud))
  
  for(k in (nref+1):nstud){
    var_s_j               = subset(data.var, S==k)
    # Odd entries (a and betax), internalized only
    B12[(2*k-1), 1]       = 0
    
    # Even entries (b and betax)
    
    frac1                     = rep(NA, n1[k])
    for( str in 1:n1[k]){
      #### for observations outside calibration study
      var_s_x                 = subset(var_s_j, strata == str)
      BXd                     = var_s_x$BXd #get all BX_d
      Xd                      = var_s_x$xd_fc
      
      Wd                      = var_s_x$Wd  #W_d 
      numerator1              = (sum(Wd*exp(BXd))+
                                   sum(Xd*beta_hat*Wd*exp(BXd)))*(1+sum(exp(BXd)))
      numerator2              = sum(Xd*exp(BXd))*sum(beta_hat*Wd*exp(BXd))
      denominator             = (1+sum(exp(BXd)))^2
      frac1[str]              = -(numerator1-numerator2)/denominator
    }
    ## Sum all components
    B12[(2*k), 1]   = (1/np)*(sum(frac1))
  }
  
  
  B[1:(2*nstud), dim_sand] = 0
  B[dim_sand, 1:(2*nstud)] = t(B12)
  ### Single corner entry for betax
  frac <- rep(NA, np)
  for(str2  in 1:np){
    var_s_j      = subset(data.var, strata2 == str2)
    Xd_s_j       = var_s_j$xd_fc
    BXd_s_j      = var_s_j$BXd
    numerator    = sum((Xd_s_j)^2*exp(BXd_s_j))*(1+sum(exp(BXd_s_j)))-(sum(Xd_s_j*exp(BXd_s_j)))^2
    denominator  = (1+sum(exp(BXd_s_j)))^2
    frac[str2]   = -numerator/denominator
  }
  B[dim_sand,dim_sand] = (1/np)*sum( frac)
  
  ## Eliminate extra entries from B associated with calibration studies
  B = B[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  ### Compute sandwich variance estimator and place variance estimate in output
  V                = solve(B)%*%A%*%t(solve(B))
  output_fc[3]    = (1/np)*V[(2*nc+1),(2*nc+1)]
  
  #### 99. Return appropriate output
  return(output_fc)
  
}

#######################################################
#######################################################
#######################################################
#######################################################
###   III. Running functions on actual data

## Example syntax: fc_df(data=mydata, X="ref",S="S",H="H",W="local",Y="Y", strata="strata", nref=2, nstud=5,
##                       covar=c("age","bmi"))

########### The function arguments of main_x are as follows:
## mydata: a dataframe of data that you have formatted in advance.
## nstud:  is the total number of studies contributing to the analysis
## nref:   is the total number of studies that used the reference lab initially (if none, then nref=0)
## H:      is an indicator variable saying whether that observation was part of the calibration subset (0=no, 1=yes, 2=NA because
##the associated study used reference lab for all measurements and thus does not calibrate)
## S:      is the study number (numbering begins at 1 and increases incrementally by 1. Studies using the reference lab initially
## for all observations must be numbered before those that did not use the reference lab)
## X:      is the reference lab measurement (NA if not available)
## W:      is the local laboratory measurement (NA if all measurements taken at reference lab initially)
## Y:      Binary outcome data (1/0 for yes/no)
## covar:  is a list of covariate names.
## strata: is the stratum number in each study. In each study, it should starts from 1 and increases
##         incrementally by 1. 


#######################################################
#######################################################



