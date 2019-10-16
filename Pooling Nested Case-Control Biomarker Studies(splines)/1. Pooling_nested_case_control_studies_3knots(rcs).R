#################################################################################
#################################################################################
###   R Code for pooling biomarkers using restricted cubic spline functions
###   (three knots)
###
###                         Yujie Wu and Molin Wang
###
###             For any questions please contact Yujie Wu
###                  Email: yujiewu@hsph.harvard.edu
#################################################################################
#################################################################################

#################################################################################
#################################################################################
###   I. Library and data sourcing
library("splines")
library("MASS")
library("Hmisc")
library("survival")

setwd("~/")

mydata = readRDS("yourdata.rds")


#################################################################################
#################################################################################
###   II. Functions
###      a. sp.basis:       returns restricted cubic spline functions
###      b. deriv.spline:   returns derivatives of cubic spline functions
###      c. int_df:         Internalized calibration with additional covariates
###      d. int_df_no_cov:  Internalized calibration without additional covariates
###      e. fc_df:          Full calibration with additional covariates
###      f. fc_df_no_cov:   Full calibration without additional covariates


#### returns restricted cubic spline basis matrix
sp.basis         = function(knots, x){
  f1 <- x
  
  f2 <- 0 * (x <= knots[1]) + ( (x - knots[1])^3 )*( x > knots[1] & x <= knots[2]) +
    ( (x - knots[1])^3-((x - knots[2])^3*(knots[3]-knots[1]))/(knots[3]-knots[2]))*(x > knots[2] & x <= knots[3]) +
    ( (x - knots[1])^3-((x - knots[2])^3*(knots[3]-knots[1]))/(knots[3]-knots[2]) + 
        ((x-knots[3])^3*(knots[2]-knots[1]))/(knots[3]-knots[2]))*(x > knots[3])
  
  spline.basis <- data.frame(f1 = f1, f2 = f2)
  return(spline.basis)
}

deriv.spline <- function(knots,x){
  f1.prime <- 1
  
  f2.prime <- 0*(x<=knots[1]) + (3*(x-knots[1])^2)*(x>knots[1] & x<=knots[2]) +
    (3*(x-knots[1])^2 - (3*(x-knots[2])^2*(knots[3]-knots[1]))/(knots[3]-knots[2]))*(x>knots[2] & x<=knots[3])+
    (3*(x-knots[1])^2 - (3*(x-knots[2])^2*(knots[3]-knots[1]))/(knots[3]-knots[2])+
       (3*(x-knots[3])^2*(knots[2]-knots[1]))/(knots[3]-knots[2]))*(x>knots[3])
  
  spline.der <- data.frame(f1.prime = f1.prime, f2.prime = f2.prime)
  return(spline.der)
}

int_df = function(data, X, S, H, W, Y, strata, nstud, nref, knots, covar=NA){
  
  #### 1. Create matrix to store output
  output_int = matrix(NA, ncol=6, nrow=1) #### to be generalized
  colnames(output_int) = c("point estimate of beta_fx1","point estimate of beta_fx2",
                           "Estimated RR_fx1", "Estimated RR_fx2", 
                           "Estimated variance of beta_fx1","Estimated variance of beta_fx2")
  
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
  
  #### 5. Create xhat_int 
  data$xhat_fc          = ifelse(data$H==2, data$X, data$a_hat + data$b_hat*data$W)
  data$xhat_int         = ifelse(data$H==0, data$xhat_fc, data$X)
  
  knots2                = quantile(data$xhat_int,probs=knots) ####to be generalized
  fxhat_int             = sp.basis(knots2, data$xhat_int)
  ns                    = ncol(fxhat_int) # number of spline functions
  fxhat_int_name        = paste0("fxhat_int_",1:ns,sep="")  
  colnames(fxhat_int)   = fxhat_int_name
  data                  = cbind(data, fxhat_int)
  
  fxhat_int.p           = deriv.spline(knots2, data$xhat_int)
  fxhat_int_name.p      = paste0("fxhat_int_",1:ns,".p",sep="")
  colnames(fxhat_int.p) = fxhat_int_name.p
  data                  = cbind(data, fxhat_int.p)
  
  #### 6. Obtain point estimate from standard logistic regression
  formula          = as.formula(paste("Y~strata(strata2)",paste(fxhat_int_name,collapse = "+"), 
                                      paste(covar,collapse ="+"),sep="+"))
  int_fit          = clogit(formula, data=data)
  beta_hat         = int_fit$coefficients[1:ns] # betahat_x value 
  betazhat         = int_fit$coefficients[-(1:ns)] # extract betazhat values 
  output_int[1:ns] = beta_hat 
  output_int[(ns+1):(2*ns)]  = exp(beta_hat) 
  
  #### 6.b. Compute Xd for each observation, Zd, and BX (linear combo)
  
  matrix.var <- c()
  for(str2 in 1:np){
    data_s_j           = subset(data,strata2 == str2)
    if(data_s_j$ncase[1] == 1){
      ncont              = data_s_j$ncont[1]
      total_num          = choose(ncont+1,1)
      pairs              = combn(seq(1,(ncont+1),1),m=1)[,-total_num] #### all possible combinations in each strata, excluding two cases being chosen
      fxd_int            = matrix(NA, nrow = length(pairs), ncol = ns)
      fp.d               = matrix(NA, nrow = length(pairs), ncol = ns)
      fpW.d              = matrix(NA, nrow = length(pairs), ncol = ns)
      Zd                 = matrix(NA, nrow = length(pairs), ncol = nz)
      S                  = matrix(data_s_j$S[1], nrow = length(pairs), ncol = 1)
      strata             = matrix(data_s_j$strata[1], nrow = length(pairs), ncol = 1)
      strata2            = matrix(str2, nrow = length(pairs), ncol = 1)
      for(i in 1:(total_num-1)){
        i1               = pairs[i]
        col_fx           = which(names(data_s_j) %in% fxhat_int_name)
        for(j in 1:length(col_fx)){
          col_num        = col_fx[j]
          fxd_int[i,j]   = data_s_j[i1, col_num] - data_s_j[(ncont+1), col_num]
        }
        
        col_fpx          = which(names(data_s_j) %in% fxhat_int_name.p)
        for(j in 1:length(col_fpx)){
          col_num        = col_fpx[j]
          fp.d[i,j]      = data_s_j[i1, col_num]*(data_s_j[i1,]$H==0) - data_s_j[(ncont+1), col_num]
        }
        
        for(j in 1:length(col_fpx)){
          col_num        = col_fpx[j]
          fpW.d[i,j]     = data_s_j[i1, col_num]*data[i1,]$W*(data_s_j[i1,]$H==0) - data_s_j[(ncont+1), col_num]*data[(ncont+1),]$W
        }
        
        col_nums         = which(names(data_s_j) %in% covar)
        for( j in 1:length(covar)){
          col_num        = col_nums[j]
          Zd[i,j]        = data_s_j[i1,col_num] - data_s_j[(ncont+1),col_num] 
        }
      }
      BXd                = fxd_int%*%beta_hat + (Zd) %*%betazhat
      matrix.temp        = cbind(fxd_int, Zd, fp.d, fpW.d, BXd, S, strata, strata2)
      matrix.var         = rbind(matrix.var, matrix.temp)
    } 
    if(data_s_j$ncase[1] == 2){
      ncont              = data_s_j$ncont[1]
      total_num          = choose(ncont+2,2)
      pairs              = combn(seq(1,(ncont+2),1),m=2)[,-total_num] #### all possible combinations in each strata, excluding two cases being chosen
      fxd_int            = matrix(NA, nrow = ncol(pairs), ncol = ns)
      fp.d               = matrix(NA, nrow = ncol(pairs), ncol = ns)
      fpW.d              = matrix(NA, nrow = ncol(pairs), ncol = ns)
      Zd                 = matrix(NA, nrow = ncol(pairs), ncol = nz)
      S                  = matrix(data_s_j$S[1], nrow = ncol(pairs), ncol = 1)
      strata             = matrix(data_s_j$strata[1], nrow = ncol(pairs), ncol = 1)
      strata2            = matrix(str2, nrow = ncol(pairs), ncol = 1)
      for(i in 1:(total_num-1)){
        i1               = pairs[,i][1]
        i2               = pairs[,i][2]
        
        col_fx           = which(names(data_s_j) %in% fxhat_int_name)
        for(j in 1:length(col_fx)){
          col_num        = col_fx[j]
          fxd_int[i,j]   = (data_s_j[i1, col_num]+data_s_j[i2, col_num]) -
            (data_s_j[(ncont+1), col_num] + data_s_j[(ncont+2), col_num])
        }
        
        col_fpx          = which(names(data_s_j) %in% fxhat_int_name.p)
        for(j in 1:length(col_fpx)){
          col_num        = col_fpx[j]
          fp.d[i,j]      = (data_s_j[i1, col_num]*(data_s_j[i1,]$H==0)+data_s_j[i2, col_num]*(data_s_j[i2,]$H==0)) - 
            (data_s_j[(ncont+1), col_num] + data_s_j[(ncont+2), col_num])
        }      
        for(j in 1:length(col_fpx)){
          col_num        = col_fpx[j]
          fpW.d[i,j]     = (data_s_j[i1, col_num]*data[i1,]$W*(data_s_j[i1,]$H==0) +data_s_j[i2, col_num]*data[i2,]$W*(data_s_j[i2,]$H==0)) - 
            (data_s_j[(ncont+1), col_num]*data[(ncont+1),]$W + data_s_j[(ncont+2), col_num]*data[(ncont+2),]$W)
        }
        
        
        col_nums         = which(names(data) %in% covar)
        for( j in 1:length(covar)){
          col_num        = col_nums[j]
          Zd[i,j]        = (data_s_j[i1,col_num] + data_s_j[i2,col_num]) - (data_s_j[(ncont+1),col_num] + data_s_j[(ncont+2),col_num])
        }
      }
      BXd                = fxd_int%*%beta_hat + (Zd) %*%betazhat
      matrix.temp        = cbind(fxd_int, Zd, fp.d, fpW.d, BXd, S, strata, strata2)
      matrix.var         = rbind(matrix.var, matrix.temp)
    } 
  }
  fxd_int.name         = paste(fxhat_int_name,"d",sep="")
  covard               = paste(covar,"d",sep = "")
  fpd.name             = paste(fxhat_int_name,"p.d",sep="")
  fpWd.name            = paste(fxhat_int_name,"pW.d",sep="")
  colnames(matrix.var) = c(fxd_int.name,covard,fpd.name,fpWd.name,"BXd","S","strata","strata2")
  data.var             = data.frame(matrix.var)
  
  #### 7. Compute variance; prepare matrices 
  dim_sand  = 2*nstud+nz+ns  #### to be generalized
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
  A12    = matrix(NA, ncol=(nz+ns), nrow=(2*nstud))
  
  for(k in 1:nstud){
    
    ### fx and a,b ###
    ### we have to use for loop to calculate psi1,psi2,psiR for each strata ###
    data_s_h     = subset(data, S==k) # the kth study
    var_s_h      = subset(data.var, S==k)
    
    psi1         = rep(NA, n1[k])
    psi2         = rep(NA, n1[k])
    col_fxd      = which(names(var_s_h) %in% fxd_int.name)
    for(i in 1:length(fxd_int.name)){
      col_num      = col_fxd[i]
      psiR         = rep(NA, n1[k])
      for(str in 1:n1[k]){
        data_s_h.temp     = subset(data_s_h, strata==str)
        data_s_h.temp2    = subset(data_s_h.temp, H==1)
        psi1[str]         = sum(data_s_h.temp2$X-a_hat[k]-b_hat[k]*data_s_h.temp2$W) # estimating equ 1
        psi2[str]         = sum((data_s_h.temp2$X-a_hat[k]-b_hat[k]*data_s_h.temp2$W)*data_s_h.temp2$W) # estimating equ 2
        
        var_s_h.temp      = subset(var_s_h, strata==str)
        fXd_s_h           = var_s_h.temp[, col_num]
        BXd_s_h           = var_s_h.temp$BXd
        psiR[str]         = -sum(fXd_s_h*exp(BXd_s_h))/(sum(exp(BXd_s_h))+1) # last estimating equation
      }
      
      A12[(2*k-1),i]      = (1/np)*sum(psi1*psiR)
      A12[(2*k),i]        = (1/np)*sum(psi2*psiR)
    }
    
    
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
      A12[(2*k-1),(ns+i)]  = (1/np)*sum(psi1*psiZ)
      A12[(2*k),(ns+i)]    = (1/np)*sum(psi2*psiZ) 
    }
  }
  
  
  A[1:(2*nstud), (2*nstud+1):dim_sand] = A12 # Fill appropriate piece of A matrix
  A21                                  = t(A12) # Reflect entries for A21
  A[(2*nstud+1):dim_sand,1:(2*nstud)]  = A21
  
  
  ## Now compute A22 entries
  A22    = matrix(NA, ncol=(nz+ns), nrow=(nz+ns))
  
  ## fx_1^2 entry
  psiX            = rep(NA, np)
  col_num         = which(names(data.var) == fxd_int.name[1])
  for(str2 in 1:np){
    var_s_h       = subset(data.var, strata2==str2)
    fXd           = var_s_h[,col_num]
    BXd           = var_s_h$BXd
    psiX[str2]    = -(sum(fXd*exp(BXd)))/(sum(exp(BXd))+1) #psiX for each strata in each study
  }
  
  A22[1,1]        = (1/np)*sum(psiX^2)
  
  ### Computations in A22 involving interaction of fx_1 between fx_2, fx_3 and betaz, 
  ### and their squared
  covard.combo      = c(fxd_int.name[-1], covard)
  for(i in 1:length(covard.combo)){ 
    col_num         = which(names(data.var) == covard.combo[i])
    psiZ            = rep(NA, np)
    for(str2 in 1:np){
      var_s_h       = subset(data.var, strata2==str2)
      ZZ            = var_s_h[,col_num]
      BXd           = var_s_h$BXd
      psiZ[str2]    = -(sum(ZZ*exp(BXd)))/(sum(exp(BXd))+1) #psiX for each strata in each study
    }
    A22[(1),(1+i)]      = (1/np)*sum(psiX*psiZ)
    A22[(1+i),(1)]      = (1/np)*sum(psiX*psiZ)
    A22[(1+i),(1+i)]    = (1/np)*sum(psiZ*psiZ) 
  }
  
  ### Computations involving cross terms of Z
  
  ## Interact the covariate psi terms (fx_2, fx_3 , Z1 and Z2...) for all pairs
  
  pairs = combn(seq(1,length(covard.combo),1),m=2)
  for(j in 1:choose(length(covard.combo),2)){ #for each pairwise combination
    pair = pairs[,j]
    
    # Obtain appropriate matrix column vector based on pairwise selection
    col_num1 = which(names(data.var) == covard.combo[pair[1]])
    col_num2 = which(names(data.var) == covard.combo[pair[2]])
    
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
  
  B12 = matrix(0, ncol=(nz+ns), nrow=(2*nstud))
  
  for(k in (nref+1):nstud){
    data_s_j               = subset(data.var, S==k)
    
    # (a and fx1,fx2,fx3
    for(i in 1:length(fxd_int.name)){ 
      col_num  = which(names(data_s_j) == fxd_int.name[i])
      col_num2 = which(names(data_s_j) == fpd.name[i])
      
      frac                  = rep(0, n1[k])
      for( str in 1:n1[k]){
        data_s_j.temp       = subset(data_s_j, strata == str)
        BXd                 = data_s_j.temp$BXd
        fkxd                = data_s_j.temp[,col_num]
        fpkxd               = data_s_j.temp[,col_num2]
        fpxd                = as.matrix(data_s_j.temp[,which(names(data_s_j)%in%fpd.name)])
        
        numerator1          = (sum(fpkxd*exp(BXd))+sum(fkxd*(fpxd%*%beta_hat)*exp(BXd)))*(1+sum(exp(BXd)))
        numerator2          = sum(fkxd*exp(BXd))*sum((fpxd%*%beta_hat)*exp(BXd))
        denominator         = (1+sum(exp(BXd)))^2
        frac[str]           = -(numerator1-numerator2)/denominator
      }
      B12[(2*k-1),i]  = (1/np)*sum(frac)
    }
    
    
    # a and betaz
    for(i in 1:nz){ 
      col_num               = which(names(data_s_j) == covard[i])
      frac                  = rep(0, n1[k])
      for( str in 1:n1[k]){
        data_s_j.temp       = subset(data_s_j, strata ==  str)
        BXd                 = data_s_j.temp$BXd
        Zkd                 = data_s_j.temp[,col_num]
        fpxd                = as.matrix(data_s_j.temp[,which(names(data_s_j)%in%fpd.name)])
        
        numerator1          = sum(Zkd*(fpxd%*%beta_hat)*exp(BXd))*(1+sum(exp(BXd)))
        numerator2          = sum(Zkd*exp(BXd))*sum((fpxd%*%beta_hat)*exp(BXd))
        denominator         = (1+sum(exp(BXd)))^2
        frac[str]           = -(numerator1-numerator2)/denominator
      }
      B12[(2*k-1),(ns+i)]  = (1/np)*sum(frac)
    }
    
    
    # (b and fx1,fx2,fx3
    for(i in 1:ns){ 
      col_num  = which(names(data_s_j) == fxd_int.name[i])
      col_num2 = which(names(data_s_j) == fpWd.name[i])
      frac                  = rep(0, n1[k])
      for( str in 1:n1[k]){
        data_s_j.temp       = subset(data_s_j, strata ==  str)
        BXd                 = data_s_j.temp$BXd
        fkxd                = data_s_j.temp[,col_num]
        fkpWd               = data_s_j.temp[,col_num2]
        fpWd                = as.matrix(data_s_j.temp[,which(names(data_s_j)%in%fpWd.name)])
        
        numerator1          = (sum(fkpWd*exp(BXd))+sum(fkxd*(fpWd%*%beta_hat)*exp(BXd)))*(1+sum(exp(BXd)))
        numerator2          = sum(fkxd*exp(BXd))*sum((fpWd%*%beta_hat)*exp(BXd))
        denominator         = (1+sum(exp(BXd)))^2
        frac[str]           = -(numerator1-numerator2)/denominator
      }
      B12[(2*k),i]  = (1/np)*sum(frac)
    }
    
    # b and betaz
    for(i in 1:nz){ 
      col_num               = which(names(data_s_j) == covard[i])
      frac                  = rep(0, n1[k])
      for( str in 1:n1[k]){
        data_s_j.temp       = subset(data_s_j, strata ==  str)
        BXd                 = data_s_j.temp$BXd
        Zkd                 = data_s_j.temp[,col_num]
        fpWd                = as.matrix(data_s_j.temp[,which(names(data_s_j)%in%fpWd.name)])
        
        numerator1          = sum(Zkd*(fpWd%*%beta_hat)*exp(BXd))*(1+sum(exp(BXd)))
        numerator2          = sum(Zkd*exp(BXd))*sum((fpWd%*%beta_hat)*exp(BXd))
        denominator         = (1+sum(exp(BXd)))^2
        frac[str]           = -(numerator1-numerator2)/denominator
      }
      B12[(2*k),(ns+i)]  = (1/np)*sum(frac)
    }
  }
  
  #### Replace B12 and its transpose in the B matrix
  B[1:(2*nstud),(2*nstud+1):dim_sand] = B12
  B[(2*nstud+1):dim_sand,1:(2*nstud)] = t(B12)
  B[1:(2*nstud),(2*nstud+1):dim_sand] = 0
  ### B22 entry : Completed component-wise #######
  
  B22    = matrix(NA, ncol=(nz+ns), nrow=(nz+ns))
  
  ## fx_1^2 entry
  frac <- rep(NA, np)
  col_num = which(names(data.var)==fxd_int.name[1])
  for(str2  in 1:np){
    var_s_j      = subset(data.var, strata2 == str2)
    fXd_s_j      = var_s_j[,col_num]
    BXd_s_j      = var_s_j$BXd
    numerator    = sum((fXd_s_j)^2*exp(BXd_s_j))*(1+sum(exp(BXd_s_j)))-(sum(fXd_s_j*exp(BXd_s_j)))^2
    denominator  = (1+sum(exp(BXd_s_j)))^2
    frac[str2]   = -numerator/denominator
  }
  B22[1,1]  = (1/np)*sum(frac)
  
  ### Computations in B22 involving fx_1 with fx_2, fx_3 and betaz, and their squared 
  
  for(i in 1:length(covard.combo)){ 
    col_num  = which(names(data.var) == covard.combo[i])
    col_num2 = which(names(data.var) == fxd_int.name[1])
    
    fracZZ    = rep(NA,np)
    fracXZ    = rep(NA,np)
    for(str2 in 1:np){
      var_s_j       = subset(data.var, strata2 == str2)
      fXd_s_j       = var_s_j[,col_num2]
      Zd_s_j        = var_s_j[,col_num]
      BXd_s_j       = var_s_j$BXd
      numerator1    = sum(fXd_s_j*Zd_s_j*exp(BXd_s_j))*(1+sum(exp(BXd_s_j)))-
        (sum(fXd_s_j*exp(BXd_s_j)))*(sum(Zd_s_j*exp(BXd_s_j)))
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
  
  pairs = combn(seq(1,length(covard.combo),1),m=2)
  for(j in 1:choose(length(covard.combo),2)){ #for each pairwise combination
    pair = pairs[,j]
    
    # Obtain appropriate matrix column vector based on pairwise selection
    col_num1 = which(names(data.var) == covard.combo[pair[1]])
    
    col_num2 = which(names(data.var) == covard.combo[pair[2]])
    
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
  V                             = solve(B)%*%A%*%t(solve(B))
  output_int[(2*ns+1):(3*ns)]   = (1/np)*diag(V)[(2*nc+1):(2*nc+ns)]
  
  #### 99. Return appropriate output
  return(output_int)
  
}

int_df_no_cov = function(data, X, S, H, W, Y, strata, nstud, nref, knots){
  
  #### 1. Create matrix to store output
  output_int = matrix(NA, ncol=6, nrow=1) #### to be generalized
  colnames(output_int) = c("point estimate of beta_fx1","point estimate of beta_fx2",
                           "Estimated RR_fx1", "Estimated RR_fx2", 
                           "Estimated variance of beta_fx1","Estimated variance of beta_fx2")
  
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
  
  #### 5. Create xhat_int 
  data$xhat_fc          = ifelse(data$H==2, data$X, data$a_hat + data$b_hat*data$W)
  data$xhat_int         = ifelse(data$H==0, data$xhat_fc, data$X)
  
  knots2                = quantile(data$xhat_int,probs=knots) ####to be generalized
  fxhat_int             = sp.basis(knots2, data$xhat_int)
  ns                    = ncol(fxhat_int)
  fxhat_int_name        = paste0("fxhat_int_",1:ns,sep="")  
  colnames(fxhat_int)   = fxhat_int_name
  data                  = cbind(data, fxhat_int)
  
  fxhat_int.p           = deriv.spline(knots2, data$xhat_int)
  fxhat_int_name.p      = paste0("fxhat_int_",1:ns,".p",sep="")
  colnames(fxhat_int.p) = fxhat_int_name.p
  data                  = cbind(data, fxhat_int.p)
  
  #### 6. Obtain point estimate from standard logistic regression
  formula                    = as.formula(paste("Y~strata(strata2)",paste(fxhat_int_name,collapse = "+"), sep="+"))
  int_fit                    = clogit(formula, data=data)
  beta_hat                   = int_fit$coefficients[1:ns] # betahat_x value 
  output_int[1:ns]           = beta_hat #### to be generalized
  output_int[(ns+1):(2*ns)]  = exp(beta_hat) #### to be generalized
  
  #### 6.b. Compute Xd for each observation, Zd, and BX (linear combo)
  
  matrix.var <- c()
  for(str2 in 1:np){
    data_s_j             = subset(data,strata2 == str2)
    if(data_s_j$ncase[1] == 1){
      ncont              = data_s_j$ncont[1]
      total_num          = choose(ncont+1,1)
      pairs              = combn(seq(1,(ncont+1),1),m=1)[,-total_num] #### all possible combinations in each strata, excluding two cases being chosen
      fxd_int            = matrix(NA, nrow = length(pairs), ncol = ns)
      fp.d               = matrix(NA, nrow = length(pairs), ncol = ns)
      fpW.d              = matrix(NA, nrow = length(pairs), ncol = ns)
      S                  = matrix(data_s_j$S[1], nrow = length(pairs), ncol = 1)
      strata             = matrix(data_s_j$strata[1], nrow = length(pairs), ncol = 1)
      strata2            = matrix(str2, nrow = length(pairs), ncol = 1)
      for(i in 1:(total_num-1)){
        i1               = pairs[i]
        col_fx           = which(names(data_s_j) %in% fxhat_int_name)
        for(j in 1:length(col_fx)){
          col_num        = col_fx[j]
          fxd_int[i,j]   = data_s_j[i1, col_num] - data_s_j[(ncont+1), col_num]
        }
        
        col_fpx          = which(names(data_s_j) %in% fxhat_int_name.p)
        for(j in 1:length(col_fpx)){
          col_num        = col_fpx[j]
          fp.d[i,j]      = data_s_j[i1, col_num]*(data_s_j[i1,]$H==0) - data_s_j[(ncont+1), col_num]
        }
        
        for(j in 1:length(col_fpx)){
          col_num        = col_fpx[j]
          fpW.d[i,j]     = data_s_j[i1, col_num]*data[i1,]$W*(data_s_j[i1,]$H==0) - data_s_j[(ncont+1), col_num]*data[(ncont+1),]$W
        }
      }
      BXd                = fxd_int%*%beta_hat 
      matrix.temp        = cbind(fxd_int, fp.d, fpW.d, BXd, S, strata, strata2)
      matrix.var         = rbind(matrix.var, matrix.temp)
    } 
    if(data_s_j$ncase[1] == 2){
      ncont              = data_s_j$ncont[1]
      total_num          = choose(ncont+2,2)
      pairs              = combn(seq(1,(ncont+2),1),m=2)[,-total_num] #### all possible combinations in each strata, excluding two cases being chosen
      fxd_int            = matrix(NA, nrow = ncol(pairs), ncol = ns)
      fp.d               = matrix(NA, nrow = ncol(pairs), ncol = ns)
      fpW.d              = matrix(NA, nrow = ncol(pairs), ncol = ns)
      S                  = matrix(data_s_j$S[1], nrow = ncol(pairs), ncol = 1)
      strata             = matrix(data_s_j$strata[1], nrow = ncol(pairs), ncol = 1)
      strata2            = matrix(str2, nrow = ncol(pairs), ncol = 1)
      for(i in 1:(total_num-1)){
        i1               = pairs[,i][1]
        i2               = pairs[,i][2]
        
        col_fx           = which(names(data_s_j) %in% fxhat_int_name)
        for(j in 1:length(col_fx)){
          col_num        = col_fx[j]
          fxd_int[i,j]   = (data_s_j[i1, col_num]+data_s_j[i2, col_num]) -
            (data_s_j[(ncont+1), col_num] + data_s_j[(ncont+2), col_num])
        }
        
        col_fpx          = which(names(data_s_j) %in% fxhat_int_name.p)
        for(j in 1:length(col_fpx)){
          col_num        = col_fpx[j]
          fp.d[i,j]      = (data_s_j[i1, col_num]*(data_s_j[i1,]$H==0)+data_s_j[i2, col_num]*(data_s_j[i2,]$H==0)) - 
            (data_s_j[(ncont+1), col_num] + data_s_j[(ncont+2), col_num])
        }      
        for(j in 1:length(col_fpx)){
          col_num        = col_fpx[j]
          fpW.d[i,j]     = (data_s_j[i1, col_num]*data[i1,]$W*(data_s_j[i1,]$H==0) +data_s_j[i2, col_num]*data[i2,]$W*(data_s_j[i2,]$H==0)) - 
            (data_s_j[(ncont+1), col_num]*data[(ncont+1),]$W + data_s_j[(ncont+2), col_num]*data[(ncont+2),]$W)
        }
      }
      BXd                = fxd_int%*%beta_hat 
      matrix.temp        = cbind(fxd_int, fp.d, fpW.d, BXd, S, strata, strata2)
      matrix.var         = rbind(matrix.var, matrix.temp)
    } 
  }
  fxd_int.name         = paste(fxhat_int_name,"d",sep="")
  fpd.name             = paste(fxhat_int_name,"p.d",sep="")
  fpWd.name            = paste(fxhat_int_name,"pW.d",sep="")
  colnames(matrix.var) = c(fxd_int.name,fpd.name,fpWd.name,"BXd","S","strata","strata2")
  data.var             = data.frame(matrix.var)
  
  #### 7. Compute variance; prepare matrices 
  dim_sand  = 2*nstud+ns  #### to be generalized
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
  A12    = matrix(NA, ncol=ns, nrow=(2*nstud))
  
  for(k in 1:nstud){
    
    ### fx and a,b ###
    ### we have to use for loop to calculate psi1,psi2,psiR for each strata ###
    data_s_h     = subset(data, S==k) # the kth study
    var_s_h      = subset(data.var, S==k)
    
    psi1         = rep(NA, n1[k])
    psi2         = rep(NA, n1[k])
    col_fxd      = which(names(var_s_h) %in% fxd_int.name)
    for(i in 1:length(fxd_int.name)){
      col_num      = col_fxd[i]
      psiR         = rep(NA, n1[k])
      for(str in 1:n1[k]){
        data_s_h.temp     = subset(data_s_h, strata==str)
        data_s_h.temp2    = subset(data_s_h.temp, H==1)
        psi1[str]         = sum(data_s_h.temp2$X-a_hat[k]-b_hat[k]*data_s_h.temp2$W) # estimating equ 1
        psi2[str]         = sum((data_s_h.temp2$X-a_hat[k]-b_hat[k]*data_s_h.temp2$W)*data_s_h.temp2$W) # estimating equ 2
        
        var_s_h.temp      = subset(var_s_h, strata==str)
        fXd_s_h           = var_s_h.temp[, col_num]
        BXd_s_h           = var_s_h.temp$BXd
        psiR[str]         = -sum(fXd_s_h*exp(BXd_s_h))/(sum(exp(BXd_s_h))+1) # last estimating equation
      }
      
      A12[(2*k-1),i]      = (1/np)*sum(psi1*psiR)
      A12[(2*k),i]        = (1/np)*sum(psi2*psiR)
    }
  }
  
  
  A[1:(2*nstud), (2*nstud+1):dim_sand] = A12 # Fill appropriate piece of A matrix
  A21                                  = t(A12) # Reflect entries for A21
  A[(2*nstud+1):dim_sand,1:(2*nstud)]  = A21
  
  
  ## Now compute A22 entries
  A22    = matrix(NA, ncol=(ns), nrow=(ns))
  
  ## fx_1^2 entry
  psiX            = rep(NA, np)
  col_num         = which(names(data.var) == fxd_int.name[1])
  for(str2 in 1:np){
    var_s_h       = subset(data.var, strata2==str2)
    fXd           = var_s_h[,col_num]
    BXd           = var_s_h$BXd
    psiX[str2]    = -(sum(fXd*exp(BXd)))/(sum(exp(BXd))+1) #psiX for each strata in each study
  }
  
  A22[1,1]        = (1/np)*sum(psiX^2)
  
  ### Computations in A22 involving fx_2, fx_3 and  their squared
  covard.combo      = c(fxd_int.name[-1])
  for(i in 1:length(covard.combo)){ 
    col_num         = which(names(data.var) == covard.combo[i])
    psiZ            = rep(NA, np)
    for(str2 in 1:np){
      var_s_h       = subset(data.var, strata2==str2)
      ZZ            = var_s_h[,col_num]
      BXd           = var_s_h$BXd
      psiZ[str2]    = -(sum(ZZ*exp(BXd)))/(sum(exp(BXd))+1) #psiX for each strata in each study
    }
    A22[(1),(1+i)]      = (1/np)*sum(psiX*psiZ)
    A22[(1+i),(1)]      = (1/np)*sum(psiX*psiZ)
    A22[(1+i),(1+i)]    = (1/np)*sum(psiZ*psiZ) 
  }
  
  ### Computations involving cross terms of fx2 fx3
  
  ## Interact the covariate psi terms (Z1 and Z2...) for all pairs
  if(length(covard.combo) >= 2){
    pairs = combn(seq(1,length(covard.combo),1),m=2)
    for(j in 1:choose(length(covard.combo),2)){ #for each pairwise combination
      pair = pairs[,j]
      
      # Obtain appropriate matrix column vector based on pairwise selection
      col_num1 = which(names(data.var) == covard.combo[pair[1]])
      col_num2 = which(names(data.var) == covard.combo[pair[2]])
      
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
  
  B12 = matrix(0, ncol=(ns), nrow=(2*nstud))
  
  for(k in (nref+1):nstud){
    data_s_j               = subset(data.var, S==k)
    
    # (a and fx1,fx2,fx3
    for(i in 1:ns){ 
      col_num  = which(names(data_s_j) == fxd_int.name[i])
      col_num2 = which(names(data_s_j) == fpd.name[i])
      
      frac                  = rep(0, n1[k])
      for( str in 1:n1[k]){
        data_s_j.temp       = subset(data_s_j, strata ==  str)
        BXd                 = data_s_j.temp$BXd
        fkxd                = data_s_j.temp[,col_num]
        fpkxd               = data_s_j.temp[,col_num2]
        fpxd                = as.matrix(data_s_j.temp[,which(names(data_s_j)%in%fpd.name)])
        
        numerator1          = (sum(fpkxd*exp(BXd))+sum(fkxd*(fpxd%*%beta_hat)*exp(BXd)))*(1+sum(exp(BXd)))
        numerator2          = sum(fkxd*exp(BXd))*sum((fpxd%*%beta_hat)*exp(BXd))
        denominator         = (1+sum(exp(BXd)))^2
        frac[str]           = -(numerator1-numerator2)/denominator
      }
      B12[(2*k-1),i]  = (1/np)*sum(frac)
    }
    
    # (b and fx1,fx2,fx3
    for(i in 1:ns){ 
      col_num  = which(names(data_s_j) == fxd_int.name[i])
      col_num2 = which(names(data_s_j) == fpWd.name[i])
      frac                  = rep(0, n1[k])
      for( str in 1:n1[k]){
        data_s_j.temp       = subset(data_s_j, strata ==  str)
        BXd                 = data_s_j.temp$BXd
        fkxd                = data_s_j.temp[,col_num]
        fkpWd               = data_s_j.temp[,col_num2]
        fpWd                = as.matrix(data_s_j.temp[,which(names(data_s_j)%in%fpWd.name)])
        
        numerator1          = (sum(fkpWd*exp(BXd))+sum(fkxd*(fpWd%*%beta_hat)*exp(BXd)))*(1+sum(exp(BXd)))
        numerator2          = sum(fkxd*exp(BXd))*sum((fpWd%*%beta_hat)*exp(BXd))
        denominator         = (1+sum(exp(BXd)))^2
        frac[str]           = -(numerator1-numerator2)/denominator
      }
      B12[(2*k),i]  = (1/np)*sum(frac)
    }
  }
  
  #### Replace B12 and its transpose in the B matrix
  B[1:(2*nstud),(2*nstud+1):dim_sand] = B12
  B[(2*nstud+1):dim_sand,1:(2*nstud)] = t(B12)
  B[1:(2*nstud),(2*nstud+1):dim_sand] = 0
  ### B22 entry : Completed component-wise #######
  
  B22    = matrix(NA, ncol=(ns), nrow=(ns))
  
  ## fx_1^2 entry
  frac <- rep(NA, np)
  col_num = which(names(data.var)==fxd_int.name[1])
  for(str2  in 1:np){
    var_s_j      = subset(data.var, strata2 == str2)
    fXd_s_j      = var_s_j[,col_num]
    BXd_s_j      = var_s_j$BXd
    numerator    = sum((fXd_s_j)^2*exp(BXd_s_j))*(1+sum(exp(BXd_s_j)))-(sum(fXd_s_j*exp(BXd_s_j)))^2
    denominator  = (1+sum(exp(BXd_s_j)))^2
    frac[str2]   = -numerator/denominator
  }
  B22[1,1]  = (1/np)*sum(frac)
  
  ### Computations in B22 involving fx_1 with fx_2, fx_3 and betaz, and their squared 
  
  for(i in 1:length(covard.combo)){ 
    col_num  = which(names(data.var) == covard.combo[i])
    col_num2 = which(names(data.var) == fxd_int.name[1])
    
    fracZZ    = rep(NA,np)
    fracXZ    = rep(NA,np)
    for(str2 in 1:np){
      var_s_j       = subset(data.var, strata2 == str2)
      fXd_s_j       = var_s_j[,col_num2]
      Zd_s_j        = var_s_j[,col_num]
      BXd_s_j       = var_s_j$BXd
      numerator1    = sum(fXd_s_j*Zd_s_j*exp(BXd_s_j))*(1+sum(exp(BXd_s_j)))-
        (sum(fXd_s_j*exp(BXd_s_j)))*(sum(Zd_s_j*exp(BXd_s_j)))
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
  if(length(covard.combo) >= 2){
    pairs = combn(seq(1,length(covard.combo),1),m=2)
    for(j in 1:choose(length(covard.combo),2)){ #for each pairwise combination
      pair = pairs[,j]
      
      # Obtain appropriate matrix column vector based on pairwise selection
      col_num1 = which(names(data.var) == covard.combo[pair[1]])
      
      col_num2 = which(names(data.var) == covard.combo[pair[2]])
      
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
  }

  ## Now place B22 in the main B matrix
  B[(2*nstud+1):dim_sand, (2*nstud+1):dim_sand] = B22
  
  ## Eliminate extra entries from B associated with noncalibration studies
  B = B[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  ### Compute sandwich variance estimator and place variance estimate in output
  V                = solve(B)%*%A%*%t(solve(B))
  output_int[(2*ns+1):(3*ns)]   = (1/np)*diag(V)[(2*nc+1):(2*nc+ns)]
  
  #### 99. Return appropriate output
  return(output_int)
}

fc_df = function(data, X, S, H, W, Y, strata, nstud, nref, knots, covar=NA){
  
  #### 1. Create matrix to store output
  output_fc = matrix(NA, ncol=6, nrow=1) #### to be generalized
  colnames(output_fc) = c("point estimate of beta_fx1","point estimate of beta_fx2",
                          "Estimated RR_fx1", "Estimated RR_fx2", 
                          "Estimated variance of beta_fx1","Estimated variance of beta_fx2")
  
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
  data$xhat_fc        = ifelse(data$H==2, data$X, data$a_hat + data$b_hat*data$W)
  
  knots2              = quantile(data$xhat_fc,probs=knots) ####to be generalized
  fxhat_fc            = sp.basis(knots2, data$xhat_fc)
  ns                  = ncol(fxhat_fc)
  fxhat_fc_name       = paste0("fxhat_fc_",1:ns,sep="")  
  colnames(fxhat_fc)  = fxhat_fc_name
  data                = cbind(data, fxhat_fc)
  
  fxhat_fc.p          = deriv.spline(knots2, data$xhat_fc)
  fxhat_fc_name.p     = paste0("fxhat_fc_",1:ns,".p",sep="")
  colnames(fxhat_fc.p)= fxhat_fc_name.p
  data                = cbind(data, fxhat_fc.p)
  
  #### 6. Obtain point estimate from standard logistic regression
  formula          = as.formula(paste("Y~strata(strata2)",paste(fxhat_fc_name,collapse = "+"), 
                                      paste(covar,collapse ="+"),sep="+"))
  fc_fit           = clogit(formula, data=data)
  beta_hat         = fc_fit$coefficients[1:ns] # betahat_x value #### to be generalized
  betazhat         = fc_fit$coefficients[-(1:ns)] # extract betazhat values #### to be generalized
  output_fc[1:ns]   = beta_hat #### to be generalized
  output_fc[(ns+1):(2*ns)]   = exp(beta_hat) #### to be generalized
  
  #### 6.b. Compute Xd for each observation, Zd, and BX (linear combo)
  
  matrix.var <- c()
  for(str2 in 1:np){
    data_s_j           = subset(data,strata2 == str2)
    if(data_s_j$ncase[1] == 1){
      ncont              = data_s_j$ncont[1]
      total_num          = choose(ncont+1,1)
      pairs              = combn(seq(1,(ncont+1),1),m=1)[,-total_num] #### all possible combinations in each strata, excluding two cases being chosen
      fxd_fc             = matrix(NA, nrow = length(pairs), ncol = ns)
      fp.d               = matrix(NA, nrow = length(pairs), ncol = ns)
      fpW.d              = matrix(NA, nrow = length(pairs), ncol = ns)
      Zd                 = matrix(NA, nrow = length(pairs), ncol = nz)
      S                  = matrix(data_s_j$S[1], nrow = length(pairs), ncol = 1)
      strata             = matrix(data_s_j$strata[1], nrow = length(pairs), ncol = 1)
      strata2            = matrix(str2, nrow = length(pairs), ncol = 1)
      for(i in 1:(total_num-1)){
        i1               = pairs[i]
        
        col_fx           = which(names(data_s_j) %in% fxhat_fc_name)
        for(j in 1:length(col_fx)){
          col_num        = col_fx[j]
          fxd_fc[i,j]    = data_s_j[i1, col_num] - data_s_j[(ncont+1), col_num]
        }
        
        col_fpx          = which(names(data_s_j) %in% fxhat_fc_name.p)
        for(j in 1:length(col_fpx)){
          col_num        = col_fpx[j]
          fp.d[i,j]      = data_s_j[i1, col_num] - data_s_j[(ncont+1), col_num]
        }
        
        for(j in 1:length(col_fpx)){
          col_num        = col_fpx[j]
          fpW.d[i,j]     = data_s_j[i1, col_num]*data[i1,]$W - data_s_j[(ncont+1), col_num]*data[(ncont+1),]$W
        }
        
        col_nums         = which(names(data_s_j) %in% covar)
        for( j in 1:length(covar)){
          col_num        = col_nums[j]
          Zd[i,j]        = data_s_j[i1,col_num] - data_s_j[(ncont+1),col_num] 
        }
      }
      BXd                = fxd_fc%*%beta_hat + (Zd) %*%betazhat
      matrix.temp        = cbind(fxd_fc, Zd, fp.d, fpW.d, BXd, S, strata, strata2)
      matrix.var         = rbind(matrix.var, matrix.temp)
    } 
    if(data_s_j$ncase[1] == 2){
      ncont              = data_s_j$ncont[1]
      total_num          = choose(ncont+2,2)
      pairs              = combn(seq(1,(ncont+2),1),m=2)[,-total_num] #### all possible combinations in each strata, excluding two cases being chosen
      fxd_fc             = matrix(NA, nrow = ncol(pairs), ncol = ns)
      fp.d               = matrix(NA, nrow = ncol(pairs), ncol = ns)
      fpW.d              = matrix(NA, nrow = ncol(pairs), ncol = ns)
      Zd                 = matrix(NA, nrow = ncol(pairs), ncol = nz)
      S                  = matrix(data_s_j$S[1], nrow = ncol(pairs), ncol = 1)
      strata             = matrix(data_s_j$strata[1], nrow = ncol(pairs), ncol = 1)
      strata2            = matrix(str2, nrow = ncol(pairs), ncol = 1)
      for(i in 1:(total_num-1)){
        i1               = pairs[,i][1]
        i2               = pairs[,i][2]
        
        col_fx           = which(names(data_s_j) %in% fxhat_fc_name)
        for(j in 1:length(col_fx)){
          col_num        = col_fx[j]
          fxd_fc[i,j]    = (data_s_j[i1, col_num]+data_s_j[i2, col_num]) -
            (data_s_j[(ncont+1), col_num] + data_s_j[(ncont+2), col_num])
        }
        
        col_fpx          = which(names(data_s_j) %in% fxhat_fc_name.p)
        for(j in 1:length(col_fpx)){
          col_num        = col_fpx[j]
          fp.d[i,j]      = (data_s_j[i1, col_num]+data_s_j[i2, col_num]) - 
            (data_s_j[(ncont+1), col_num] + data_s_j[(ncont+2), col_num])
        }      
        for(j in 1:length(col_fpx)){
          col_num        = col_fpx[j]
          fpW.d[i,j]     = (data_s_j[i1, col_num]*data[i1,]$W +data_s_j[i2, col_num]*data[i2,]$W) - 
            (data_s_j[(ncont+1), col_num]*data[(ncont+1),]$W + data_s_j[(ncont+2), col_num]*data[(ncont+2),]$W)
        }
        
        
        col_nums         = which(names(data) %in% covar)
        for( j in 1:length(covar)){
          col_num        = col_nums[j]
          Zd[i,j]        = (data_s_j[i1,col_num] + data_s_j[i2,col_num]) - (data_s_j[(ncont+1),col_num] + data_s_j[(ncont+2),col_num])
        }
      }
      BXd                = fxd_fc%*%beta_hat + (Zd) %*%betazhat
      matrix.temp        = cbind(fxd_fc, Zd, fp.d, fpW.d, BXd, S, strata, strata2)
      matrix.var         = rbind(matrix.var, matrix.temp)
    } 
  }
  fxd_fc.name          = paste(fxhat_fc_name,"d",sep="")
  covard               = paste(covar,"d",sep = "")
  fpd.name             = paste(fxhat_fc_name,"p.d",sep="")
  fpWd.name            = paste(fxhat_fc_name,"pW.d",sep="")
  colnames(matrix.var) = c(fxd_fc.name,covard,fpd.name,fpWd.name,"BXd","S","strata","strata2")
  data.var             = data.frame(matrix.var)
  
  #### 7. Compute variance; prepare matrices 
  dim_sand  = 2*nstud+nz+ns  #### to be generalized
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
  A12    = matrix(NA, ncol=(nz+ns), nrow=(2*nstud))
  
  for(k in 1:nstud){
    
    ### fx and a,b ###
    ### we have to use for loop to calculate psi1,psi2,psiR for each strata ###
    data_s_h     = subset(data, S==k) # the kth study
    var_s_h      = subset(data.var, S==k)
    
    psi1         = rep(NA, n1[k])
    psi2         = rep(NA, n1[k])
    col_fxd      = which(names(var_s_h) %in% fxd_fc.name)
    for(i in 1:ns){
      col_num      = col_fxd[i]
      psiR         = rep(NA, n1[k])
      for(str in 1:n1[k]){
        data_s_h.temp     = subset(data_s_h, strata==str)
        data_s_h.temp2    = subset(data_s_h.temp, H==1)
        psi1[str]         = sum(data_s_h.temp2$X-a_hat[k]-b_hat[k]*data_s_h.temp2$W) # estimating equ 1
        psi2[str]         = sum((data_s_h.temp2$X-a_hat[k]-b_hat[k]*data_s_h.temp2$W)*data_s_h.temp2$W) # estimating equ 2
        
        var_s_h.temp      = subset(var_s_h, strata==str)
        fXd_s_h           = var_s_h.temp[, col_num]
        BXd_s_h           = var_s_h.temp$BXd
        psiR[str]         = -sum(fXd_s_h*exp(BXd_s_h))/(sum(exp(BXd_s_h))+1) # last estimating equation
      }
      
      A12[(2*k-1),i]      = (1/np)*sum(psi1*psiR)
      A12[(2*k),i]        = (1/np)*sum(psi2*psiR)
    }
    
    
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
      A12[(2*k-1),(ns+i)]  = (1/np)*sum(psi1*psiZ)
      A12[(2*k),(ns+i)]    = (1/np)*sum(psi2*psiZ) 
    }
  }
  
  
  A[1:(2*nstud), (2*nstud+1):dim_sand] = A12 # Fill appropriate piece of A matrix
  A21                                  = t(A12) # Reflect entries for A21
  A[(2*nstud+1):dim_sand,1:(2*nstud)]  = A21
  
  
  ## Now compute A22 entries
  A22    = matrix(NA, ncol=(nz+ns), nrow=(nz+ns))
  
  ## fx_1^2 entry
  psiX            = rep(NA, np)
  col_num         = which(names(data.var) == fxd_fc.name[1])
  for(str2 in 1:np){
    var_s_h       = subset(data.var, strata2==str2)
    fXd           = var_s_h[,col_num]
    BXd           = var_s_h$BXd
    psiX[str2]    = -(sum(fXd*exp(BXd)))/(sum(exp(BXd))+1) #psiX for each strata in each study
  }
  
  A22[1,1]        = (1/np)*sum(psiX^2)
  
  ### Computations in A22 involving fx_2, fx_3 and betaz, and their squared
  covard.combo      = c(fxd_fc.name[-1], covard)
  for(i in 1:length(covard.combo)){ 
    col_num         = which(names(data.var) == covard.combo[i])
    psiZ            = rep(NA, np)
    for(str2 in 1:np){
      var_s_h       = subset(data.var, strata2==str2)
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
  
  pairs = combn(seq(1,length(covard.combo),1),m=2)
  for(j in 1:choose(length(covard.combo),2)){ #for each pairwise combination
    pair = pairs[,j]
    
    # Obtain appropriate matrix column vector based on pairwise selection
    col_num1 = which(names(data.var) == covard.combo[pair[1]])
    col_num2 = which(names(data.var) == covard.combo[pair[2]])
    
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
  
  B12 = matrix(0, ncol=(nz+ns), nrow=(2*nstud))
  
  for(k in (nref+1):nstud){
    data_s_j               = subset(data.var, S==k)
    
    # (a and fx1,fx2,fx3
    for(i in 1:ns){ 
      col_num  = which(names(data_s_j) == fxd_fc.name[i])
      col_num2 = which(names(data_s_j) == fpd.name[i])
      
      frac                  = rep(0, n1[k])
      for( str in 1:n1[k]){
        data_s_j.temp       = subset(data_s_j, strata ==  str)
        BXd                 = data_s_j.temp$BXd
        fkxd                = data_s_j.temp[,col_num]
        fpkxd               = data_s_j.temp[,col_num2]
        fpxd                = as.matrix(data_s_j.temp[,which(names(data_s_j)%in%fpd.name)])
        
        numerator1          = (sum(fpkxd*exp(BXd))+sum(fkxd*(fpxd%*%beta_hat)*exp(BXd)))*(1+sum(exp(BXd)))
        numerator2          = sum(fkxd*exp(BXd))*sum((fpxd%*%beta_hat)*exp(BXd))
        denominator         = (1+sum(exp(BXd)))^2
        frac[str]           = -(numerator1-numerator2)/denominator
      }
      B12[(2*k-1),i]  = (1/np)*sum(frac)
    }
    
    
    # a and betaz
    for(i in 1:nz){ 
      col_num               = which(names(data_s_j) == covard[i])
      frac                  = rep(0, n1[k])
      for( str in 1:n1[k]){
        data_s_j.temp       = subset(data_s_j, strata ==  str)
        BXd                 = data_s_j.temp$BXd
        Zkd                 = data_s_j.temp[,col_num]
        fpxd                = as.matrix(data_s_j.temp[,which(names(data_s_j)%in%fpd.name)])
        
        numerator1          = sum(Zkd*(fpxd%*%beta_hat)*exp(BXd))*(1+sum(exp(BXd)))
        numerator2          = sum(Zkd*exp(BXd))*sum((fpxd%*%beta_hat)*exp(BXd))
        denominator         = (1+sum(exp(BXd)))^2
        frac[str]           = -(numerator1-numerator2)/denominator
      }
      B12[(2*k-1),(ns+i)]  = (1/np)*sum(frac)
    }
    
    
    # (b and fx1,fx2,fx3
    for(i in 1:length(fxd_fc.name)){ 
      col_num  = which(names(data_s_j) == fxd_fc.name[i])
      col_num2 = which(names(data_s_j) == fpWd.name[i])
      frac                  = rep(0, n1[k])
      for( str in 1:n1[k]){
        data_s_j.temp       = subset(data_s_j, strata ==  str)
        BXd                 = data_s_j.temp$BXd
        fkxd                = data_s_j.temp[,col_num]
        fkpWd               = data_s_j.temp[,col_num2]
        fpWd                = as.matrix(data_s_j.temp[,which(names(data_s_j)%in%fpWd.name)])
        
        numerator1          = (sum(fkpWd*exp(BXd))+sum(fkxd*(fpWd%*%beta_hat)*exp(BXd)))*(1+sum(exp(BXd)))
        numerator2          = sum(fkxd*exp(BXd))*sum((fpWd%*%beta_hat)*exp(BXd))
        denominator         = (1+sum(exp(BXd)))^2
        frac[str]           = -(numerator1-numerator2)/denominator
      }
      B12[(2*k),i]  = (1/np)*sum(frac)
    }
    
    # b and betaz
    for(i in 1:nz){ 
      col_num               = which(names(data_s_j) == covard[i])
      frac                  = rep(0, n1[k])
      for( str in 1:n1[k]){
        data_s_j.temp       = subset(data_s_j, strata ==  str)
        BXd                 = data_s_j.temp$BXd
        Zkd                 = data_s_j.temp[,col_num]
        fpWd                = as.matrix(data_s_j.temp[,which(names(data_s_j)%in%fpWd.name)])
        
        numerator1          = sum(Zkd*(fpWd%*%beta_hat)*exp(BXd))*(1+sum(exp(BXd)))
        numerator2          = sum(Zkd*exp(BXd))*sum((fpWd%*%beta_hat)*exp(BXd))
        denominator         = (1+sum(exp(BXd)))^2
        frac[str]           = -(numerator1-numerator2)/denominator
      }
      B12[(2*k),(ns+i)]  = (1/np)*sum(frac)
    }
  }
  
  #### Replace B12 and its transpose in the B matrix
  B[1:(2*nstud),(2*nstud+1):dim_sand] = B12
  B[(2*nstud+1):dim_sand,1:(2*nstud)] = t(B12)
  B[1:(2*nstud),(2*nstud+1):dim_sand] = 0
  ### B22 entry : Completed component-wise #######
  
  B22    = matrix(NA, ncol=(nz+ns), nrow=(nz+ns))
  
  ## fx_1^2 entry
  frac <- rep(NA, np)
  col_num = which(names(data.var)==fxd_fc.name[1])
  for(str2  in 1:np){
    var_s_j      = subset(data.var, strata2 == str2)
    fXd_s_j      = var_s_j[,col_num]
    BXd_s_j      = var_s_j$BXd
    numerator    = sum((fXd_s_j)^2*exp(BXd_s_j))*(1+sum(exp(BXd_s_j)))-(sum(fXd_s_j*exp(BXd_s_j)))^2
    denominator  = (1+sum(exp(BXd_s_j)))^2
    frac[str2]   = -numerator/denominator
  }
  B22[1,1]  = (1/np)*sum(frac)
  
  ### Computations in B22 involving fx_1 with fx_2, fx_3 and betaz, and their squared 
  
  for(i in 1:length(covard.combo)){ 
    col_num  = which(names(data.var) == covard.combo[i])
    col_num2 = which(names(data.var) == fxd_fc.name[1])
    
    fracZZ    = rep(NA,np)
    fracXZ    = rep(NA,np)
    for(str2 in 1:np){
      var_s_j       = subset(data.var, strata2 == str2)
      fXd_s_j       = var_s_j[,col_num2]
      Zd_s_j        = var_s_j[,col_num]
      BXd_s_j       = var_s_j$BXd
      numerator1    = sum(fXd_s_j*Zd_s_j*exp(BXd_s_j))*(1+sum(exp(BXd_s_j)))-
        (sum(fXd_s_j*exp(BXd_s_j)))*(sum(Zd_s_j*exp(BXd_s_j)))
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
  
  pairs = combn(seq(1,length(covard.combo),1),m=2)
  for(j in 1:choose(length(covard.combo),2)){ #for each pairwise combination
    pair = pairs[,j]
    
    # Obtain appropriate matrix column vector based on pairwise selection
    col_num1 = which(names(data.var) == covard.combo[pair[1]])
    
    col_num2 = which(names(data.var) == covard.combo[pair[2]])
    
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
  output_fc[(2*ns+1):(3*ns)]   = (1/np)*diag(V)[(2*nc+1):(2*nc+ns)]
  
  #### 99. Return appropriate output
  return(output_fc)
  
}

fc_df_no_cov = function(data, X, S, H, W, Y, strata, nstud, nref, knots){
  
  #### 1. Create matrix to store output
  output_fc = matrix(NA, ncol=6, nrow=1) #### to be generalized
  colnames(output_fc) = c("point estimate of beta_fx1","point estimate of beta_fx2",
                          "Estimated RR_fx1", "Estimated RR_fx2", 
                          "Estimated variance of beta_fx1","Estimated variance of beta_fx2")
  
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
  data$xhat_fc        = ifelse(data$H==2, data$X, data$a_hat + data$b_hat*data$W)
  
  knots2              = quantile(data$xhat_fc,probs=knots) ####to be generalized
  fxhat_fc            = sp.basis(knots2, data$xhat_fc)
  ns                  = ncol(fxhat_fc)
  fxhat_fc_name       = paste0("fxhat_fc_",1:ns,sep="")  ####to be generalized
  colnames(fxhat_fc)  = fxhat_fc_name
  data                = cbind(data, fxhat_fc)
  
  fxhat_fc.p          = deriv.spline(knots2, data$xhat_fc)
  fxhat_fc_name.p     = paste0("fxhat_fc_",1:ns,".p",sep="")
  colnames(fxhat_fc.p)= fxhat_fc_name.p
  data                = cbind(data, fxhat_fc.p)
  
  #### 6. Obtain point estimate from standard logistic regression
  formula                    = as.formula(paste("Y~strata(strata2)",paste(fxhat_fc_name,collapse = "+"), sep="+"))
  fc_fit                     = clogit(formula, data=data)
  beta_hat                   = fc_fit$coefficients[1:ns] # betahat_x value #### to be generalized
  output_fc[1:ns]            = beta_hat #### to be generalized
  output_fc[(ns+1):(2*ns)]   = exp(beta_hat) #### to be generalized
  
  #### 6.b. Compute Xd for each observation, Zd, and BX (linear combo)
  
  matrix.var <- c()
  for(str2 in 1:np){
    data_s_j           = subset(data,strata2 == str2)
    if(data_s_j$ncase[1] == 1){
      ncont              = data_s_j$ncont[1]
      total_num          = choose(ncont+1,1)
      pairs              = combn(seq(1,(ncont+1),1),m=1)[,-total_num] #### all possible combinations in each strata, excluding two cases being chosen
      fxd_fc             = matrix(NA, nrow = length(pairs), ncol = ns)
      fp.d               = matrix(NA, nrow = length(pairs), ncol = ns)
      fpW.d              = matrix(NA, nrow = length(pairs), ncol = ns)
      S                  = matrix(data_s_j$S[1], nrow = length(pairs), ncol = 1)
      strata             = matrix(data_s_j$strata[1], nrow = length(pairs), ncol = 1)
      strata2            = matrix(str2, nrow = length(pairs), ncol = 1)
      for(i in 1:(total_num-1)){
        i1               = pairs[i]
        
        col_fx           = which(names(data_s_j) %in% fxhat_fc_name)
        for(j in 1:length(col_fx)){
          col_num        = col_fx[j]
          fxd_fc[i,j]    = data_s_j[i1, col_num] - data_s_j[(ncont+1), col_num]
        }
        
        col_fpx          = which(names(data_s_j) %in% fxhat_fc_name.p)
        for(j in 1:length(col_fpx)){
          col_num        = col_fpx[j]
          fp.d[i,j]      = data_s_j[i1, col_num] - data_s_j[(ncont+1), col_num]
        }
        
        for(j in 1:length(col_fpx)){
          col_num        = col_fpx[j]
          fpW.d[i,j]     = data_s_j[i1, col_num]*data[i1,]$W - data_s_j[(ncont+1), col_num]*data[(ncont+1),]$W
        }
      }
      BXd                = fxd_fc%*%beta_hat 
      matrix.temp        = cbind(fxd_fc, fp.d, fpW.d, BXd, S, strata, strata2)
      matrix.var         = rbind(matrix.var, matrix.temp)
    } 
    if(data_s_j$ncase[1] == 2){
      ncont              = data_s_j$ncont[1]
      total_num          = choose(ncont+2,2)
      pairs              = combn(seq(1,(ncont+2),1),m=2)[,-total_num] #### all possible combinations in each strata, excluding two cases being chosen
      fxd_fc             = matrix(NA, nrow = ncol(pairs), ncol = ns)
      fp.d               = matrix(NA, nrow = ncol(pairs), ncol = ns)
      fpW.d              = matrix(NA, nrow = ncol(pairs), ncol = ns)
      S                  = matrix(data_s_j$S[1], nrow = ncol(pairs), ncol = 1)
      strata             = matrix(data_s_j$strata[1], nrow = ncol(pairs), ncol = 1)
      strata2            = matrix(str2, nrow = ncol(pairs), ncol = 1)
      for(i in 1:(total_num-1)){
        i1               = pairs[,i][1]
        i2               = pairs[,i][2]
        
        col_fx           = which(names(data_s_j) %in% fxhat_fc_name)
        for(j in 1:length(col_fx)){
          col_num        = col_fx[j]
          fxd_fc[i,j]    = (data_s_j[i1, col_num]+data_s_j[i2, col_num]) -
            (data_s_j[(ncont+1), col_num] + data_s_j[(ncont+2), col_num])
        }
        
        col_fpx          = which(names(data_s_j) %in% fxhat_fc_name.p)
        for(j in 1:length(col_fpx)){
          col_num        = col_fpx[j]
          fp.d[i,j]      = (data_s_j[i1, col_num]+data_s_j[i2, col_num]) - 
            (data_s_j[(ncont+1), col_num] + data_s_j[(ncont+2), col_num])
        }      
        for(j in 1:length(col_fpx)){
          col_num        = col_fpx[j]
          fpW.d[i,j]     = (data_s_j[i1, col_num]*data[i1,]$W +data_s_j[i2, col_num]*data[i2,]$W) - 
            (data_s_j[(ncont+1), col_num]*data[(ncont+1),]$W + data_s_j[(ncont+2), col_num]*data[(ncont+2),]$W)
        }
      }
      BXd                = fxd_fc%*%beta_hat 
      matrix.temp        = cbind(fxd_fc, fp.d, fpW.d, BXd, S, strata, strata2)
      matrix.var         = rbind(matrix.var, matrix.temp)
    } 
  }
  fxd_fc.name          = paste(fxhat_fc_name,"d",sep="")
  fpd.name             = paste(fxhat_fc_name,"p.d",sep="")
  fpWd.name            = paste(fxhat_fc_name,"pW.d",sep="")
  colnames(matrix.var) = c(fxd_fc.name,fpd.name,fpWd.name,"BXd","S","strata","strata2")
  data.var             = data.frame(matrix.var)
  
  #### 7. Compute variance; prepare matrices 
  dim_sand  = 2*nstud+ns  #### to be generalized
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
  A12    = matrix(NA, ncol=(ns), nrow=(2*nstud))
  
  for(k in 1:nstud){
    
    ### fx and a,b ###
    ### we have to use for loop to calculate psi1,psi2,psiR for each strata ###
    data_s_h     = subset(data, S==k) # the kth study
    var_s_h      = subset(data.var, S==k)
    
    psi1         = rep(NA, n1[k])
    psi2         = rep(NA, n1[k])
    col_fxd      = which(names(var_s_h) %in% fxd_fc.name)
    for(i in 1:ns){
      col_num      = col_fxd[i]
      psiR         = rep(NA, n1[k])
      for(str in 1:n1[k]){
        data_s_h.temp     = subset(data_s_h, strata==str)
        data_s_h.temp2    = subset(data_s_h.temp, H==1)
        psi1[str]         = sum(data_s_h.temp2$X-a_hat[k]-b_hat[k]*data_s_h.temp2$W) # estimating equ 1
        psi2[str]         = sum((data_s_h.temp2$X-a_hat[k]-b_hat[k]*data_s_h.temp2$W)*data_s_h.temp2$W) # estimating equ 2
        
        var_s_h.temp      = subset(var_s_h, strata==str)
        fXd_s_h           = var_s_h.temp[, col_num]
        BXd_s_h           = var_s_h.temp$BXd
        psiR[str]         = -sum(fXd_s_h*exp(BXd_s_h))/(sum(exp(BXd_s_h))+1) # last estimating equation
      }
      
      A12[(2*k-1),i]      = (1/np)*sum(psi1*psiR)
      A12[(2*k),i]        = (1/np)*sum(psi2*psiR)
    }
  }
  
  A[1:(2*nstud), (2*nstud+1):dim_sand] = A12 # Fill appropriate piece of A matrix
  A21                                  = t(A12) # Reflect entries for A21
  A[(2*nstud+1):dim_sand,1:(2*nstud)]  = A21
  
  
  ## Now compute A22 entries
  A22    = matrix(NA, ncol=(ns), nrow=(ns))
  
  ## fx_1^2 entry
  psiX            = rep(NA, np)
  col_num         = which(names(data.var) == fxd_fc.name[1])
  for(str2 in 1:np){
    var_s_h       = subset(data.var, strata2==str2)
    fXd           = var_s_h[,col_num]
    BXd           = var_s_h$BXd
    psiX[str2]    = -(sum(fXd*exp(BXd)))/(sum(exp(BXd))+1) #psiX for each strata in each study
  }
  
  A22[1,1]        = (1/np)*sum(psiX^2)
  
  ### Computations in A22 involving fx_2, fx_3  and their squared
  covard.combo      = c(fxd_fc.name[-1])
  for(i in 1:length(covard.combo)){ 
    col_num         = which(names(data.var) == covard.combo[i])
    psiZ            = rep(NA, np)
    for(str2 in 1:np){
      var_s_h       = subset(data.var, strata2==str2)
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
  if(length(covard.combo) >= 2){
    pairs = combn(seq(1,length(covard.combo),1),m=2)
    for(j in 1:choose(length(covard.combo),2)){ #for each pairwise combination
      pair = pairs[,j]
      
      # Obtain appropriate matrix column vector based on pairwise selection
      col_num1 = which(names(data.var) == covard.combo[pair[1]])
      col_num2 = which(names(data.var) == covard.combo[pair[2]])
      
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
  
  B12 = matrix(0, ncol=(ns), nrow=(2*nstud))
  
  for(k in (nref+1):nstud){
    data_s_j               = subset(data.var, S==k)
    
    # (a and fx1,fx2,fx3
    for(i in 1:length(fxd_fc.name)){ 
      col_num  = which(names(data_s_j) == fxd_fc.name[i])
      col_num2 = which(names(data_s_j) == fpd.name[i])
      
      frac                  = rep(0, n1[k])
      for( str in 1:n1[k]){
        data_s_j.temp       = subset(data_s_j, strata ==  str)
        BXd                 = data_s_j.temp$BXd
        fkxd                = data_s_j.temp[,col_num]
        fpkxd               = data_s_j.temp[,col_num2]
        fpxd                = as.matrix(data_s_j.temp[,which(names(data_s_j)%in%fpd.name)])
        
        numerator1          = (sum(fpkxd*exp(BXd))+sum(fkxd*(fpxd%*%beta_hat)*exp(BXd)))*(1+sum(exp(BXd)))
        numerator2          = sum(fkxd*exp(BXd))*sum((fpxd%*%beta_hat)*exp(BXd))
        denominator         = (1+sum(exp(BXd)))^2
        frac[str]           = -(numerator1-numerator2)/denominator
      }
      B12[(2*k-1),i]  = (1/np)*sum(frac)
    }
    
    # (b and fx1,fx2,fx3
    for(i in 1:length(fxd_fc.name)){ 
      col_num  = which(names(data_s_j) == fxd_fc.name[i])
      col_num2 = which(names(data_s_j) == fpWd.name[i])
      frac                  = rep(0, n1[k])
      for( str in 1:n1[k]){
        data_s_j.temp       = subset(data_s_j, strata ==  str)
        BXd                 = data_s_j.temp$BXd
        fkxd                = data_s_j.temp[,col_num]
        fkpWd               = data_s_j.temp[,col_num2]
        fpWd                = as.matrix(data_s_j.temp[,which(names(data_s_j)%in%fpWd.name)])
        
        numerator1          = (sum(fkpWd*exp(BXd))+sum(fkxd*(fpWd%*%beta_hat)*exp(BXd)))*(1+sum(exp(BXd)))
        numerator2          = sum(fkxd*exp(BXd))*sum((fpWd%*%beta_hat)*exp(BXd))
        denominator         = (1+sum(exp(BXd)))^2
        frac[str]           = -(numerator1-numerator2)/denominator
      }
      B12[(2*k),i]  = (1/np)*sum(frac)
    }
  }
  
  #### Replace B12 and its transpose in the B matrix
  B[1:(2*nstud),(2*nstud+1):dim_sand] = B12
  B[(2*nstud+1):dim_sand,1:(2*nstud)] = t(B12)
  B[1:(2*nstud),(2*nstud+1):dim_sand] = 0
  ### B22 entry : Completed component-wise #######
  
  B22    = matrix(NA, ncol=(ns), nrow=(ns))
  
  ## fx_1^2 entry
  frac <- rep(NA, np)
  col_num = which(names(data.var)==fxd_fc.name[1])
  for(str2  in 1:np){
    var_s_j      = subset(data.var, strata2 == str2)
    fXd_s_j      = var_s_j[,col_num]
    BXd_s_j      = var_s_j$BXd
    numerator    = sum((fXd_s_j)^2*exp(BXd_s_j))*(1+sum(exp(BXd_s_j)))-(sum(fXd_s_j*exp(BXd_s_j)))^2
    denominator  = (1+sum(exp(BXd_s_j)))^2
    frac[str2]   = -numerator/denominator
  }
  B22[1,1]  = (1/np)*sum(frac)
  
  ### Computations in B22 involving fx_1 with fx_2, fx_3, and their squared 
  
  for(i in 1:length(covard.combo)){ 
    col_num  = which(names(data.var) == covard.combo[i])
    col_num2 = which(names(data.var) == fxd_fc.name[1])
    
    fracZZ    = rep(NA,np)
    fracXZ    = rep(NA,np)
    for(str2 in 1:np){
      var_s_j       = subset(data.var, strata2 == str2)
      fXd_s_j       = var_s_j[,col_num2]
      Zd_s_j        = var_s_j[,col_num]
      BXd_s_j       = var_s_j$BXd
      numerator1    = sum(fXd_s_j*Zd_s_j*exp(BXd_s_j))*(1+sum(exp(BXd_s_j)))-
        (sum(fXd_s_j*exp(BXd_s_j)))*(sum(Zd_s_j*exp(BXd_s_j)))
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
  if(length(covard.combo) >= 2){
    pairs = combn(seq(1,length(covard.combo),1),m=2)
    for(j in 1:choose(length(covard.combo),2)){ #for each pairwise combination
      pair = pairs[,j]
      
      # Obtain appropriate matrix column vector based on pairwise selection
      col_num1 = which(names(data.var) == covard.combo[pair[1]])
      
      col_num2 = which(names(data.var) == covard.combo[pair[2]])
      
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
  }

  
  ## Now place B22 in the main B matrix
  B[(2*nstud+1):dim_sand, (2*nstud+1):dim_sand] = B22
  
  ## Eliminate extra entries from B associated with noncalibration studies
  B = B[(2*nref+1):dim_sand,(2*nref+1):dim_sand]
  
  ### Compute sandwich variance estimator and place variance estimate in output
  V                = solve(B)%*%A%*%t(solve(B))
  output_fc[(2*ns+1):(3*ns)]   = (1/np)*diag(V)[(2*nc+1):(2*nc+ns)]
  
  #### 99. Return appropriate output
  return(output_fc)
}

#######################################################
#######################################################
#######################################################
#######################################################
###   III. Running functions on actual data

## Example syntax: fc_df(data=mydata, X="ref",S="S",H="H",W="local",Y="Y", nref=2, nstud=5, knots=c(0.25,0.5,0.75),
##                       covar=c("age","bmi"))
##                 fc_df_no_cov(data=mydata, X="ref",S="S",H="H",W="local",Y="Y", nref=2, nstud=5, knots=c(0.25,0.5,0.75))

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
## knots:  The place where the knots are placed. K knots yield K-1 spline terms. For example
##knots=c(0.25,0.5,0.75)
## covar:  is a list of covariate names.


#######################################################
#######################################################









