# NEW PALMS

multi_signal_lasso <- function(Y,
                               phi,
                               lam,
                               etaini=NULL,
                               admmAbsTol=1e-2,
                               admmRelTol=1e-2){
  
  f <- function(phi, Y, x, lam){ # objective function
    #Y=n*1
    #X=p*1
    #phi=n*p
    1/2*norm(Y -  phi %*% matrix(x,ncol=1)  , "2")^2 + lam*sum(
      apply(data.frame(abs(x),abs(x-1)),1,min)
    ) # 修改为annie Qu老师的惩罚
  }
  
  S <- function(z,theta1,theta2){ # the soft threshold function
    if(z<=0){return(min(z,0))}
    if(0<z & z<=1+theta2){return(max(0,z-theta2))}
    if(z>1+theta2){max(1,z-theta1)}
  }
  
  require(quadprog) # solving the Quadratic Programming Problem
  t1 <- proc.time()#timer
  p <- ncol(phi)
  n <- length(Y)
  #if(nrow(phi)!=n | ncol(phi)!=p){print("输入矩阵维数错误！")
  # return(NULL)}
  # basic setting &  transformation
  if(is.null(etaini)){
    x  <- rep(0.5,p)
  }else{
    x  <- etaini
  }
  eta <- x   
  nu <- eta-eta ##  dual variable
  indt <- 1 # number of repeats totally
  eps <- sqrt(length(x))+2
  
  while(eps > 0.1 & indt < 50){ #ind : ADMM max repeats
    #print(indt)
    old.obj <- f(phi, Y, x, lam)
    #update x ; coordinate decent 
    old.x = x
    for(j in 1:p){
      #print(paste("p=",j))
      epsil = 1 # initial residuals
      r <- Y  - as.matrix(phi[,-j],ncol=ncol(phi)-1) %*% matrix(x[-j],ncol=1)
      tem <-  phi[,j] %*% r / sum(phi[,j]^2 ) # 判别，即xhat
      # print(tem)
      if(abs(tem)<=abs(tem-1)){
        x[j] <- sign(tem)*max(0,abs(tem) - n * lam/sum(phi[,j]^2))
      }
      if(abs(tem>abs(tem-1))){
        x[j] <- 1 + sign(tem-1)*max(0, abs(tem-1) - n * lam/sum(phi[,j]^2)) 
      }
    }
    eps <- sum(abs(x-old.x))
    #print(eps)
    indt <- indt + 1
  }
  
  # Return and timer
  t2 <- proc.time();t=t2-t1#timer
  #print(paste("【计算完成】，已返回求解结果。用时：",t[3][[1]],'秒'))
  return(as.numeric(x))
}

 

PALMS <- function(trainY,
                  trainPhi,
                  lam,  
                 etaini=NULL,rho=10,
                 admmAbsTol=0.1,admmRelTol=0.1,
                 Gamma_tuning=1){
  
  t1 <- proc.time() #timer
  
  if(nrow(trainPhi)>ncol(trainPhi)){
    tem <- lm(trainY~trainPhi)$coef[-1]
    #hist(tem)
    w <- tem
  }
  
  if(nrow(trainPhi)<=ncol(trainPhi)){
    w <- as.numeric(multi_signal_lasso(trainY,trainPhi,
                                                lam, 
                                                etaini=etaini,
                                                admmAbsTol=1,admmRelTol=1))
  }
  
  t2 <- proc.time()
  
  wei <- as.numeric(apply(cbind(abs(w),abs(w-1)),1,min)) # 距离0和1最小的距离
  
  # if one step over, refit and return
  iiind <- 0
  if(sum(wei)==0){#print("MDSL已经全部压缩为01");
    #print("ALL 01")
    iiind <- 1
    return(list(w,
                w,
                t2-t1,
                t2-t1,
                iiind))} # 全部完全压缩为0和1，则直接返回结果
  if(sum(wei)!=0){
    #print("【NOT】———————————— 01")
    ind <- which(wei==0)  #  否则就标记出来，将非0，1的系数进行惩罚回归
    indi <- which(wei!=0) #   ind for Fixed coef; indi for non-fixed coef     
    phi2 <- t( t(trainPhi[,indi]) * wei[indi]^Gamma_tuning )  # Fix the coefficients has already been estimated as 0 or 1
    if(length(indi)==1){
      # phi2 <- matrix( phi2, ncol = 1 )
      # r <- Y- as.matrix(phi[,ind]) %*% matrix(x[ind],ncol=1) # residuals
      return(list(w,
                  w,
                  t2-t1,
                  t2-t1,
                  iiind))
    } #防止单列情况报错，进行格式调整
    
    if(length(ind)==1){
      Y2 <- trainY- matrix(trainPhi[,ind],ncol=1) %*% matrix(w[ind],ncol=1) 
    }
    
    if(length(ind)!=1){
      Y2 <- trainY- as.matrix(trainPhi[,ind]) %*% matrix(w[ind],ncol=1)  
    }
    
    t3 <- proc.time()
    x2 <- constrained_MDSL(as.numeric(Y2),
                           phi2,
                           lam,
                           rho = rho, #用极大的权重限制估计区间，因为加权后的惩罚比较大
                           CMALweights = 1 / (wei[indi]^Gamma_tuning ) ,       ###!!!!!重要
                           etaini =  etaini[indi] , # 
                           admmAbsTol=1e-2,
                           admmRelTol=1e-2)
    t4 <- proc.time()
    final <- c()
    final[indi] <- x2   *   wei[indi]^Gamma_tuning   ## [这里也做相应的调整]
    final[ind] <- w[ind]
    # final[which(final!=1 & final !=0)] <- one.step.x[which(final!=1 & final !=0)]
    
    final
    
    return(list(as.numeric(w),
                as.numeric(final),
                t2-t1,
                t4-t1,
                iiind))}
}









PALMS_cv <- function(Y,
                    phi,
                    lam.list=seq(0.01,10,length = 50),Gamma_tuning=1,
                    etaini,rho,x,ct="cv",c2=1 # ct for crateria 
){
  require(sBIC)
  n = 5
  
  require(glmnet)
  group <- c(1:n) # 为了防止有一组没有观测
  
  
  mse <- c()
  for(lam in 1:length(lam.list)){
    mse.tem <- c()
    for(g in 1:n){
      trainY <- Y[-which(group==g)]
      testY  <- Y[which(group==g)]
      trainPhi <- phi[-which(group==g),]
      testPhi <- phi[which(group==g),] 
      etaini=rep(0.5,ncol(trainPhi))
      xhat <- PALMS(trainY,trainPhi,
                   lam=lam.list[lam],
                   rho=rho,
                   etaini=etaini,Gamma_tuning=Gamma_tuning)[[2]]
      #mse.tem[g] <- mean((x-xhat)^2)
      if(ct=="cv"){mse.tem[g] <- mean((testY-testPhi%*%xhat)^2)}
      if(ct=="bic"){
        # modified BIC cerition
        sbic = log(sum((testY-testPhi%*%xhat)^2)) + c2* log(length(testY))*(length(xhat)-sum(xhat==0|xhat ==1))
        mse.tem[g] <- sbic}
      #print(mse.tem)
    }
    mse[lam] <- mean(na.omit(mse.tem))
  }
  (lambda_final <- min(which(mse==min(mse))))
  #print(lam.list[lambda_final])
  xhat_final <- PALMS(Y,phi,rho = rho,
                     lam=lam.list[lambda_final],
                     etaini=etaini)
  return(xhat_final)
}



