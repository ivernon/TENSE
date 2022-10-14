############################################################################################################
#########               TENSE Torn Embedding Non-Stationary Emulation Functions                    #########
############################################################################################################



### Standard simple Bayes Linear/Gaussian Process Stationary Emulator ###
simple_BL_emulator <- function(x1,y1,xp,theta=1,sig,nugget=0.1,nu=5/2,just_var=1,em_mcmc_notf=1,gaus_or_Mater=1,mu=0){
  D <- y1
  x1 <- t(t(x1)/theta)				# rescale inputs and prediction inputs by theta
  xp <- t(t(xp)/theta)
  
  if(gaus_or_Mater) Cmat <- function(d) sig^2 * exp(-d^2) * (1-nugget) else								# exponential structure
  if(!gaus_or_Mater) Cmat <- function(d) return( Matern(d, range = 1, nu= nu, phi=sig^2) * (1-nugget))		# Matern structure
  
  covBD <- Cmat( as.matrix(pdist(xp,x1)) )
  varD  <- Cmat( as.matrix(dist(x1)) ) + diag(sig^2*nugget,nrow=nrow(x1))
  varB  <- Cmat( as.matrix(dist(xp)) ) + em_mcmc_notf*diag(sig^2*nugget,nrow=nrow(xp)) # em_mcmc_notf = 1 => variance includes noise
  EB    <- mu																			 # em_mcmc_notf = 0 => emul f(x) itself not
  ED 	  <- mu
  varD_inv <- chol2inv(chol(varD))
  ED_B   <- EB + covBD %*% varD_inv %*% (D - ED)
  VarD_B <- varB - covBD %*% varD_inv %*% t(covBD)
  if(just_var) return(list(ED_B,diag(VarD_B))) else
  if(!just_var) return(list(ED_B=ED_B,VarD_B=VarD_B,varD=varD,varB=varB,covBD=covBD,corD=varD/sig^2,corB=varB/sig^2,corBD=covBD/sig^2))
}



### Bayes Linear Non-stationary Emulator ###
simple_Non_Stationary_emulator <- function(xp,x1,y1,mu=0,sig=1,nugget=10^(-6),just_var=1,just_varD=0,matrix_out=0){

	D <- y1 
	np <- nrow(xp)
	n1 <- nrow(x1)
	
	# construct cov matrix element wise, done some efficiency improvements.
	varD <- diag(sig^2/2,n1)				# half sig^2 as gets doubled when transpose added
	covBD <- matrix(0,nrow=np,ncol=n1)

	for(i in 1:(n1-1)) for(j in (i+1):n1)  varD[i,j]  <-  (1-nugget)*k_NS(x=x1[i,],y=x1[j,],sig=sig)		# only calc upper diagonal
	for(i in 1:np) 	   for(j in 1:n1)	    covBD[i,j]  <-  (1-nugget)*k_NS(x=xp[i,],y=x1[j,],sig=sig)		# only calc upper diagonal
	varD <- varD + t(varD)					# add transpose to form full matrix, now with correct sig^2 diagonal
	if(just_varD==1) return(varD)

	if(just_var==1){
		diag_varB <- rep(sig^2,np)		# if just_var=1 we only need just the variances of varB
	}
	if(just_var==0){
		varB <- diag(sig^2/2,np)			# half sig^2 as gets doubled when transpose added
		for(i in 1:(np-1)) for(j in (i+1):np)  varB[i,j]  <-  (1-nugget)*k_NS(x=xp[i,],y=xp[j,],sig=sig)		# only calc upper diagonal
		varB <- varB + t(varB)				# add transpose to form full matrix, now with correct sig^2 diagonal
	}
	
	### evaluated Bayes Linear Update ###
	EB    <- mu												
	ED 	  <- mu
	covBD_varDinv <- covBD %*% chol2inv(chol(varD))			# need this for full VarD_B calc so may as well use it for ED_B
	ED_B  <- EB + covBD_varDinv %*% (D - ED) 				
	if(just_var==1) diag_VarD_B	<- diag_varB -  apply(covBD_varDinv * covBD,1,sum)		# fast way of just getting the diagonal elements
	if(just_var==0) VarD_B <- varB - covBD_varDinv %*% t(covBD)				                # if full covariance matrix is required	
	
	### Return results ###
	if(just_var & !matrix_out) return(list(ED_B=c(as.matrix(ED_B)),diag_VarD_B=diag_VarD_B)) else 
	if(!just_var & !matrix_out) return(list(ED_B=ED_B,VarD_B=VarD_B,varB=varB,varD=varD,covBD=covBD)) else
	if(just_var & matrix_out) return(cbind("ED_B"=c(ED_B),"diag_VarD_B"=diag_VarD_B))
	print("Non-Stationary emulator Complete")
}



### Batch wrapper for Non-stationary BL emulator. Splits evaluation points into batches for efficiency gains ###
batch_Non_Stationary_emulator <- function(x1,y1,xp,pts_per_batch=100,mu=0,sig=0.7,nugget=0.00001){
	
	np <- nrow(xp)
	em_mean <- 1:np								# set up emulator mean and variance vectors
	em_var <- 1:np
	
	batches <- ceiling(np / pts_per_batch)		# number of batches to perform
	
	for(m in 1:batches){
		ind <- ( 1 + (m-1)*pts_per_batch ):( min(m*pts_per_batch,np) )			#	 set up indices including last tricky shorter set
		GPem <- simple_Non_Stationary_emulator(x1=x1,y1=y1,xp=xp[ind,,drop=FALSE],mu=mu,sig=sig,nugget=nugget,just_var=1)  
		em_mean[ind] <- c(GPem[[1]])							# emulator mean mu(x) for input points xp
		em_var[ind]  <-   GPem[[2]]								# emulator var sigma^2(x) for input points xp	
		cat("Completed batch ",m," out of ",batches,"\n")
	}
	return(list(em_mean=em_mean,em_var=em_var))
}



### Bayes Linear Non-stationary Emulator in batches using futures for parallel computation ###
### batch up the emulator evaluation points into list and use FUTURE lapply to evaluate the emulator ### 
batch_future_lapply_Non_Stationary_emulator <- function(x1,y1,xp,pts_per_batch=100,mu=0,sig=0.7,nugget=0.00001){
  
  np <- nrow(xp)
  batches <- ceiling(np / pts_per_batch)		# number of batches to perform
  xp_batch_list <- vector("list",batches)   # split evaluation points xp into smaller batches in a list 
  cat(1)
  for(m in 1:batches){
    ind <- ( 1 + (m-1)*pts_per_batch ):( min(m*pts_per_batch,np) )			#	 set up indices including last tricky shorter set
    xp_batch_list[[m]] <- xp[ind,,drop=FALSE]               # create smaller batches of evaluation points in list
  }
  cat(2)
  
  ### parallel futures applied across list, (may not fit into array as eval no. not divisible by core no.) ###
  GPem <- future_lapply(xp_batch_list,simple_Non_Stationary_emulator,x1=x1,y1=y1,mu=mu,sig=sig,nugget=nugget,just_var=1)    
  cat(3)
  
  ### convert batch list back to preallocated vectors keep memory under control ###
  em_mean <- 1:np					  # set up emulator mean and variance vectors
  em_var <- 1:np  
  for(m in 1:batches){
    ind <- ( 1 + (m-1)*pts_per_batch ):( min(m*pts_per_batch,np) )			#	 set up indices including last tricky shorter set
    em_mean[ind] <- GPem[[m]]$ED_B								      # batch emulator mean mu(x) for input points xp
    em_var[ind]  <- GPem[[m]]$diag_VarD_B								# batch emulator var sigma^2(x) for input points xp	
  }
  cat(4)
  
  return(list(em_mean=em_mean,em_var=em_var))
}




































