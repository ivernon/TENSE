##################################################################################################################
######                  Emulators with Torn Embeddings for Discontinuities                 ######
##################################################################################################################

### Created using R version 4.2.1 (2022-06-23) -- "Funny-Looking Kid" ###

### Ensure current working directory is set to that of "TENSE Emulating functions with multiple discontinuities.R" ###

### Ensure plotting directory exists called "plots/" ###

setwd("/Users/ianvernon/Work/Oil/Roxar/Roxar Analysis 2018 - 2021/Well Placement/Project with Jonathan Carter/TNO challenge 2/Example Code for toy models in Paper TNO Challenge II")

library("pdist")
library("fields")
library("viridisLite")

### If you are happy to use "future" package (for parallel evaluation using futures) set 'use_futures' to TRUE
### Only used for final example (4 non-linear discontinuities, figure 3) 
use_futures <- FALSE      # set to TRUE to use futures and parallelisation 
if(use_futures){
  library("future")               
  library("future.apply")           
}

### source main emulation functions ###
source("TENSE Torn Embedding Non-Stationary Emulation Functions.R")


##################################################################################################################
### Simple 2d function with single discontinuity along straight line ###

#################################################
### Plot Figure 1: simple 2d toy function

axis_labels <- c(expression(x),expression(y))

pdf("plots/example_1discontinuity_figure1.pdf",height=8,width=9)

	### define function f(x) with discontinuity ###
	sim2d_fun <- function(x,y) 0.4*sin(5*x) + 0.4*cos(5*y) + 0.8*(x>0.75)*(x-0.75)^2*sign(y-1)

	### plot true function f(x) in Figure 1(a) ###
	xseq <- seq(0,2,len=80)
	xp <- as.matrix(expand.grid(xseq,xseq))    # define grid of points for plotting
	y_true <- sim2d_fun(xp[,1],xp[,2])
	y_true_mat <- matrix(y_true,nrow=length(xseq),ncol=length(xseq))
	filled.contour(x=xseq,y=xseq,z=y_true_mat,color.palette=terrain.colors,levels=pretty(range(y_true_mat),20),xlab=axis_labels[1],ylab=axis_labels[2],
						plot.axes={axis(1);axis(2)
									lines(c(0.75,3),c(1,1),lwd=4,col=1)
									})
	
	### embed in higher dimension function ###
	raise_to_higher_dimension <- function(x,y) {
		xh <- cbind(x,y,0)										            # raise 2d x into 3d xh space			
		xh[,3] <- -0.4*(x>0.75)*(x-0.75)^2*sign(y-1)			# embed in 3D using v(x,y) 
		return(xh)
	}
	
	### design simple 4x4 grid of points (simple design to highlight effects of emulation with embedding)
	x1seq <- seq(0.25,1.75,len=4)
	x1 <- as.matrix(expand.grid(x1seq,x1seq))
	y1 <- sim2d_fun(x1[,1],x1[,2])              # evaluate f(x) at 16 design points
	
	### embed in 3D ###
	xh1 <- raise_to_higher_dimension(x1[,1],x1[,2])   # embed 2D design points in 3D
	xhp <- raise_to_higher_dimension(xp[,1],xp[,2])   # embed 2D evaluation points in 3D
	
	### evaluate simple emulator in 3D ###
	GPem <- simple_BL_emulator(x1=xh1,y1=y1,xp=xhp,theta=0.5,sig=0.7,nugget=0.00001,gaus_or_Mater=1,mu=0)
	yp <- c(GPem[[1]])						# emulator mean = E_D[f(x)] for input points xp
	ys <- sqrt(GPem[[2]])					# emulator sd = sqrt(Var_D[f(x)]) for input points xp
	
	### sort into matrices for plotting ###
	yp_mat <- matrix(yp,nrow=length(xseq),ncol=length(xseq))
	ys_mat <- matrix(ys,nrow=length(xseq),ncol=length(xseq))
	yembed_mat <- matrix(xhp[,3],nrow=length(xseq),ncol=length(xseq))
	
	### plot embedding surface v(x,y) ###
	filled.contour(x=xseq,y=xseq,z=yembed_mat,color.palette=magma,levels=pretty(range(yembed_mat),20),xlab=axis_labels[1],ylab=axis_labels[2],
	               plot.axes={axis(1);axis(2)
	                 lines(c(0.75,3),c(1,1),lwd=4,col=1)
	               })
	### plot emulator mean E_D[f(x)] ###
	filled.contour(x=xseq,y=xseq,z=yp_mat,color.palette=terrain.colors,levels=pretty(range(y_true_mat),20),xlab=axis_labels[1],ylab=axis_labels[2],
						plot.axes={axis(1);axis(2)
									lines(c(0.75,3),c(1,1),lwd=4,col=1)
									points(x1,pch=16)
									})
	### plot emulator sd sqrt(Var_D[f(x)]) ###
	filled.contour(x=xseq,y=xseq,z=ys_mat,color.palette=cm.colors,xlab=axis_labels[1],ylab=axis_labels[2],
						plot.axes={axis(1);axis(2)
									lines(c(0.75,3),c(1,1),lwd=4,col=1)
									points(x1,pch=16)
									})

dev.off()


###################################################################
### Plot figure 8, appendix C: realisations from the emulator 

library(MASS)

### define lower resolution grid ###
# xseq <- seq(0,2,len=60)     # slow detailed version
xseq <- seq(0,2,len=30)       # fast, coarse version
xp <- as.matrix(expand.grid(xseq,xseq))
xhp <- raise_to_higher_dimension(xp[,1],xp[,2])

### evaluate simple emulator in 3D, but request full covariance matrices ###
mu <- 0
GPem <- simple_BL_emulator(x1=xh1,y1=y1,xp=xhp,theta=0.5,sig=0.7,nugget=0.00001,gaus_or_Mater=1,mu=mu,just_var=0)
varB <- GPem$varB          # prior covariance matrix Var[f(xp)]
EB <- rep(mu,nrow(xhp))    # prior expectation vector E[f(xp)]
VarD_B <- GPem$VarD_B      # adjusted covariance matrix Var_D[f(xp)]
ED_B   <- GPem$ED_B        # adjusted expectation vector E_D[f(xp)]
set.seed(3)
ndraw <- 4                 # just do 4 realisations
y_draw_prior <- mvrnorm(n=ndraw,mu=EB,Sigma=varB)       # evaluate prior realisations using normal assumption
y_draw_post  <- mvrnorm(n=ndraw,mu=ED_B,Sigma=VarD_B)   # evaluate adjusted realisations using normal assumption

### Plot Figure 8 in appendix C ### 
pdf("plots/example_1discontinuity_realisations_figure8.pdf",height=8,width=9)
for(j in 1:2){    # loop over prior then adjusted
	for(i in 1:ndraw) {
		if(j==1) y_draw_mat <- matrix(y_draw_prior[i,],nrow=length(xseq),ncol=length(xseq))
		if(j==2) y_draw_mat <- matrix(y_draw_post[i,],nrow=length(xseq),ncol=length(xseq))
		filled.contour(x=xseq,y=xseq,z=y_draw_mat,color.palette=inferno,main=paste("Realisations from",c("prior","posterior")[j]),xlab=axis_labels[1],ylab=axis_labels[2],
							plot.axes={axis(1);axis(2)
										lines(c(0.75,3),c(1,1),lwd=4,col=1)
										if(j==2) points(x1,pch=16)
										})
	}
}
dev.off()

### End Simple 2d function with single discontinuity along straight line ###
##################################################################################################################



########################################################################################################################################
### Two discontinuities of differing length, using full TENSE structure

### define function with two discontinuities ###
yrift <- c(0.75,1.25)						# y-coordinate of rifts (i.e. discontinuities)
xrift <- c(0.6,1)								# x-coordinate of left edge of rifts (i.e. discontinuities)
sim2d_fun <- function(x,y) 0.4*sin(5*x) + 0.4*cos(5*y) + 
  1.2*(x>xrift[2])*(x-xrift[2])^2*(y>yrift[2]) - 0.6*(x>xrift[1])*(x-xrift[1])^2*(y<yrift[1])

### place poins close to rift to help clarity of plotting ###
rdev <- 0.005					# distance of some points from the rift 
xseq <- seq(0,2,len=40)
yseq <- sort(c(xseq,yrift[1] + rdev,yrift[1] - rdev,yrift[2] + rdev,yrift[2] - rdev)) 

### look at true function ###
xp <- as.matrix(expand.grid(xseq,yseq))
y_true <- sim2d_fun(xp[,1],xp[,2])
y_true_mat <- matrix(y_true,nrow=length(xseq),ncol=length(yseq))
filled.contour(x=xseq,y=yseq,z=y_true_mat,color.palette=terrain.colors,levels=pretty(range(y_true_mat),20),xlab=axis_labels[1],ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 lines(c(xrift[2],3),c(yrift[2],yrift[2]),lwd=3,col=1)		# two rifts
                 lines(c(xrift[1],3),c(yrift[1],yrift[1]),lwd=3,col=1)		# two rifts
               })

### embed in higher dimension using v(x,y) ###
raise_to_higher_dimension <- function(x,y) {
  xh <- cbind(x,y,0)										# raise 2d x into 3d xh space			
  gradient_between_rift_endpoints <- (yrift[2]-yrift[1]) / (xrift[2]-xrift[1])
  xrift_bet <- xrift[1] + (y-yrift[1]) / gradient_between_rift_endpoints 			# x point on line between rift left endpoints.
  xh[,3] <- 0.6*(x-xrift_bet)^2*(x>xrift_bet)*(y<yrift[2])*(y>yrift[1]) 	-		# between two rifts
    0.6*(x-xrift[1])^2*(x>xrift[1])*(y<yrift[1])
  return(xh)
}

### function to calculate partial derivatives of v(x,y) surface ###
deriv_of_v_surface <- function(x,y,partial_derv="x1") {			# partial derivatives of embedded surface v(x,y) for Sigma function
  gradient_between_rift_endpoints <- (yrift[2]-yrift[1]) / (xrift[2]-xrift[1])
  xrift_bet <- xrift[1] + (y-yrift[1]) / gradient_between_rift_endpoints 			# x point on line between rift left endpoints.
  if(partial_derv=="x1") return( 0.6*2*(x-xrift_bet)*(x>xrift_bet)*(y<yrift[2])*(y>yrift[1]) - 
                                   0.6*2*(x-xrift[1])*(x>xrift[1])*(y<yrift[1]) )
  if(partial_derv=="x2") return( 0.6*2*(x-xrift_bet)*(-1/gradient_between_rift_endpoints)*(x>xrift_bet)*(y<yrift[2])*(y>yrift[1]) ) 	
}

### Design and perform a set of 12 runs ###
x1seq <- seq(0.25,1.75,len=4)
y1seq <- seq(0.375,1.625,len=3)
x1 <- as.matrix(expand.grid(x1seq,y1seq))
y1 <- sim2d_fun(x1[,1],x1[,2])

### Raise runs and evaluation grid into 3rd dimension using v(x,y) ###
xh1 <- raise_to_higher_dimension(x1[,1],x1[,2])
xhp <- raise_to_higher_dimension(xp[,1],xp[,2])
raised_mat <- matrix(xhp[,3],nrow=length(xseq),ncol=length(yseq))

### TENSE Non-Stationary Covariance Setup ###
Q 	  <- function(x,y) (x-y) %*% solve( (Sigma(x) + Sigma(y))/2 ) %*% (x-y)   # Quadratic Form, section 3.2
k_S   <- function(d,sig=1) sig^2*exp(-d^2)     # squared exponential covariance function
k_NS  <- function(x,y,sig=1) {                 # non-stationary covariance function, section 3.2
  p <- length(x)
  2^(p/2) * det(Sigma(x))^(1/4) * det(Sigma(y))^(1/4)  /  det(Sigma(x) + Sigma(y))^(1/2)  *  k_S( d=sqrt(Q(x,y)),sig=sig )
}

### Now define two choices for Sigma_3D(x) ###
Sigma_diag <- function(x) diag(theta(x)^2)						  # Simple isotropic diagonal matrix for Sigma_3D (i.e. no correction for warping)
Sigma_real <- function(x,theta=0.5,alph3=0.25) {				# The full Sigma_3D using TENSE framework equation (3.28) 
  v_x <- deriv_of_v_surface(x=x[1],y=x[2],partial_derv="x1")
  v_y <- deriv_of_v_surface(x=x[1],y=x[2],partial_derv="x2")
  w_3 <- c(-v_x,-v_y,1)
  r2 <- v_x^2 + v_y^2
  return(										
    theta^2 * matrix( c(    1,   0,   v_x,
                            0,   1,   v_y,
                            v_x, v_y, r2   ),  nrow=3, ncol=3, byrow=TRUE)   +   alph3^2 / (r2+1) * w_3 %*% t(w_3)
  )
}

pdf("plots/example_2unequal_discontinuities_figure2.pdf",height=8,width=9)

### plot real function f(x,y) ###
filled.contour(x=xseq,y=yseq,z=y_true_mat,color.palette=terrain.colors,levels=pretty(range(y_true_mat),20),xlab=axis_labels[1],ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 lines(c(xrift[2],3),c(yrift[2],yrift[2]),lwd=3,col=1)		# two rifts
                 lines(c(xrift[1],3),c(yrift[1],yrift[1]),lwd=3,col=1)		# two rifts
               })

### plot embedding into 3rd dimension ###
filled.contour(x=xseq,y=yseq,z=raised_mat,color.palette=magma,levels=pretty(range(raised_mat),50),xlab=axis_labels[1],ylab=axis_labels[2],
               plot.axes={axis(1);axis(2)
                 lines(c(xrift[2],3),c(yrift[2],yrift[2]),lwd=3,col=1)		# two rifts
                 lines(c(xrift[1],3),c(yrift[1],yrift[1]),lwd=3,col=1)		# two rifts
               })

### loop over k=1: simple diagonal Sigma_3D, k=2: full TENSE Sigma_3D with warping correction ###
for(k in 1:2){				
  if(k==1) {
    theta <- function(x) 0.5*c(1,1,0.7)
    Sigma <- Sigma_diag							# Simple diagonal matrix for Sigma_3D (i.e. no correction for warping)
  }
  if(k==2) Sigma <- Sigma_real			# The full Sigma_3D(x) using TENSE framework (i.e. with warping correction)
  
  ### perform deep GP emulation, using batches ###	
  GPem   <- batch_Non_Stationary_emulator(x1=xh1,y1=y1,xp=xhp[,],pts_per_batch=min(1000,nrow(xhp)),mu=0,sig=0.7,nugget=0.00001)
  yp_mat <- matrix(GPem[[1]],nrow=length(xseq),ncol=length(yseq))			    # matrix of emulator mean mu(x) for input points xp
  ys_mat <- matrix(sqrt(GPem[[2]]),nrow=length(xseq),ncol=length(yseq))   # matrix of emulator sd sigma(x) for input points xp
  
  ### plot emulator mean ###
  filled.contour(x=xseq,y=yseq,z=yp_mat,color.palette=terrain.colors,levels=pretty(range(yp_mat),20),
                  plot.axes={axis(1);axis(2);points(x1,pch=16)
                   lines(c(xrift[2],3),c(yrift[2],yrift[2]),lwd=3,col=1)		# two rifts
                   lines(c(xrift[1],3),c(yrift[1],yrift[1]),lwd=3,col=1)		# two rifts
                 } )
  ### plot emulator sd ###
  filled.contour(x=xseq,y=yseq,z=ys_mat,color.palette=cm.colors,levels=pretty(c(0,0.7),14),xlab=axis_labels[1],ylab=axis_labels[2],
                 plot.axes={axis(1);axis(2);points(x1,pch=16)
                   lines(c(xrift[2],3),c(yrift[2],yrift[2]),lwd=3,col=1)		# two rifts
                   lines(c(xrift[1],3),c(yrift[1],yrift[1]),lwd=3,col=1)		# two rifts
                 } )
}

dev.off()

### End Two discontinuities of differing length, using full TENSE structure ###
########################################################################################################################################



########################################################################################################################################
### More general example with 4 non-linear discontinuities


### define function to classify region ###
determine_region <- function(x,y,circ_rad=0.4){    # the region identifier function m(x) in paper
  r <- sqrt(x^2 + y^2)
  if(r < circ_rad) reg1 <- 0 else { 
    if(x<=0 & y<0 ) if( (x^2 + (y+1)^2) < 1) reg1 <- 4 else reg1 <- 3
    if(x<0 & y>=0 ) if( ((x+1)^2 + y^2) < 1) reg1 <- 3 else reg1 <- 2
    if(x>=0 & y>0 ) if( (x^2 + (y-1)^2) < 1) reg1 <- 2 else reg1 <- 1
    if(x>0 & y<=0 ) if( ((x-1)^2 + y^2) < 1) reg1 <- 1 else reg1 <- 4
  }
  return(reg1)   # return which region point is in
}

### define function f(x,y) with curved discontinuities (sim_not_embed==1), and embedding v(x,y) (sim_not_embed==0) ###
sim2d_fun_swirl <- function(xvec,sim_not_embed=0,em_power=2,em_mult=-0.5) {
  x <- xvec[1]
  y <- xvec[2]
  r <- sqrt(x^2 + y^2)
  circ_rad <- 0.4

  reg1 <- determine_region(x,y,circ_rad=circ_rad)
  
  base_func <- 0.5*(sin(3*x) + cos(3.5*y))
  if(sim_not_embed==1) if(reg1==0) full_fun <- base_func else full_fun <- base_func + 2*(reg1 %% 2 - 0.5)*(r-circ_rad)^2    # returns simulator f(x,y)
  if(sim_not_embed==0) if(reg1==0) full_fun <- 0 else full_fun <- em_mult*c(1,0,-1,0)[reg1]*(r-circ_rad)^em_power           # returns embedding surface v(x,y)
  full_fun
}

### setup grid of evaluation points ###
xseq <- seq(-1,1,len=100)
yseq <- xseq

### define points along curved rifts ###
circ_rad <- 0.4
rift_pts1 <- cbind(xseq, sqrt(1-(xseq+1)^2))
rift_pts1 <- rift_pts1[!is.nan(rift_pts1[,2]),]
rift_pts1 <- rift_pts1[apply(rift_pts1^2,1,sum)>(circ_rad^2),]

### rotate first set of rift points to get 3 other sets of rift points ###
rot_mat <- matrix(c(0,1, -1,0),byrow=TRUE,ncol=2)
rift_pts2 <- rift_pts1%*%rot_mat
rift_pts3 <- rift_pts2%*%rot_mat
rift_pts4 <- rift_pts3%*%rot_mat

### function for plotting rift points ###
plot_rift_points <- function(){
  lines(rift_pts1[,1],rift_pts1[,2],lwd=5,col=1)		# two rifts
  lines(rift_pts2[,1],rift_pts2[,2],lwd=5,col=1)		# two rifts
  lines(rift_pts3[,1],rift_pts3[,2],lwd=5,col=1)		# two rifts
  lines(rift_pts4[,1],rift_pts4[,2],lwd=5,col=1)		# two rifts
}

### look at true function ###
xp <- as.matrix(expand.grid(xseq,yseq))
y_true <- apply(xp,1,sim2d_fun_swirl,sim_not_embed=1)     # evaluate f(x,y) on xp points
y_true_mat <- matrix(y_true,nrow=length(xseq),ncol=length(yseq))
filled.contour(x=xseq,y=yseq,z=y_true_mat,color.palette=terrain.colors,levels=pretty(range(y_true_mat),20),xlab=axis_labels[1],ylab=axis_labels[2],
				plot.axes={axis(1); axis(2); plot_rift_points()}
)

### look at embedding ###
y_embed <- apply(xp,1,sim2d_fun_swirl,sim_not_embed=0)    # evaluate v(x,y) on xp points
y_embed_mat <- matrix(y_embed,nrow=length(xseq),ncol=length(yseq))     
filled.contour(x=xseq,y=yseq,z=y_embed_mat,color.palette=magma,levels=pretty(range(y_embed_mat),20),xlab=axis_labels[1],ylab=axis_labels[2],
               plot.axes={axis(1); axis(2); plot_rift_points()}
)

### embed in higher dimension ###
raise_to_higher_dimension <- function(x,y) {
  xh <- cbind(x,y,0)										# raise 2d x into 3d xh space			
  xh[,3] <- apply(xh[,1:2],1,sim2d_fun_swirl,sim_not_embed=0)
  return(xh)
}

### function to calculate partial derivatives of v(x,y) surface ###
deriv_of_v_surface <- function(x,y,partial_derv="x1",circ_rad=0.4,em_power=2,em_mult=-0.5) {			# partial derivatives of embedded surface v(x,y) for Sigma function
     r <- sqrt(x^2 + y^2)
     reg1 <- determine_region(x,y,circ_rad=circ_rad)
  if(partial_derv=="x1") if(reg1==0) return(0) else return(em_mult*c(1,0,-1,0)[reg1]*em_power*(r-circ_rad)^(em_power-1) * x / r )
  if(partial_derv=="x2") if(reg1==0) return(0) else return(em_mult*c(1,0,-1,0)[reg1]*em_power*(r-circ_rad)^(em_power-1) * y / r ) 	
}

### wrapper for embedding ###
deriv_of_v_surface_wrapper <- function(x1,partial_derv="x1",circ_rad=0.4,em_power=2,em_mult=-0.5){
  x <- x1[1]
  y <- x1[2]
  return(deriv_of_v_surface(x,y,partial_derv=partial_derv,circ_rad=circ_rad,em_power=em_power,em_mult=em_mult))
  }

### Design and perform a set of runs ###
x1seq <- seq(-0.95,0.95,len=4)
y1seq <- seq(-0.95,0.95,len=4)
x1 <- as.matrix(expand.grid(x1seq,y1seq))
y1 <- apply(x1,1,sim2d_fun_swirl,sim_not_embed=1)   # evaluate f(x,y) on x1 points

### Raise runs and grid into 3rd dimension using v(x,y) ###
xh1 <- raise_to_higher_dimension(x1[,1],x1[,2])
xhp <- raise_to_higher_dimension(xp[,1],xp[,2])
raised_mat <- matrix(xhp[,3],nrow=length(xseq),ncol=length(yseq))

### TENSE Non-Stationary Covariance Setup ###
Q 	  <- function(x,y) (x-y) %*% solve( (Sigma(x) + Sigma(y))/2 ) %*% (x-y)      # Quadratic Form, section 3.2
k_S   <- function(d,sig=1) sig^2*exp(-d^2)            # squared exponential covariance function
k_NS  <- function(x,y,sig=1) {                        # non-stationary covariance function, section 3.2
	p <- length(x)
	2^(p/2) * det(Sigma(x))^(1/4) * det(Sigma(y))^(1/4)  /  det(Sigma(x) + Sigma(y))^(1/2)  *  k_S( d=sqrt(Q(x,y)),sig=sig )
}

### Now define the full Sigma_3D(x) using TENSE framework  ###
Sigma_real <- function(x,theta=0.5,alph3=0.25) {			# The full Sigma_3D(x) using TENSE framework equation (3.28) 
	v_x <- deriv_of_v_surface(x=x[1],y=x[2],partial_derv="x1")
	v_y <- deriv_of_v_surface(x=x[1],y=x[2],partial_derv="x2")
	w_3 <- c(-v_x,-v_y,1)
	r2 <- v_x^2 + v_y^2
	return(										
	theta^2 * matrix( c(    1,   0,   v_x,
						            	0,   1,   v_y,
							           v_x, v_y,  r2  ),  nrow=3, ncol=3, byrow=TRUE)  +  alph3^2 / (r2+1) * w_3 %*% t(w_3)
		  )
}

pdf("plots/example_4curved_discontinuites_figure3.pdf",height=8,width=9)

### plot embedding into 3rd dimension ###
filled.contour(x=xseq,y=yseq,z=raised_mat,color.palette=magma,levels=pretty(range(raised_mat),50),xlab=axis_labels[1],ylab=axis_labels[2],
               plot.axes={axis(1);axis(2);plot_rift_points()}
)
### plot real function ###
filled.contour(x=xseq,y=yseq,z=y_true_mat,color.palette=terrain.colors,levels=pretty(range(y_true_mat),20),xlab=axis_labels[1],ylab=axis_labels[2],
               plot.axes={axis(1);axis(2);plot_rift_points()}
)

Sigma <- Sigma_real								# The full Sigma_3D(x) using TENSE framework (i.e. with warping correction)

### perform non-stationary emulation, using batches, not using future parallelisation ###	
if(!use_futures) GPem <- batch_Non_Stationary_emulator(x1=xh1,y1=y1,xp=xhp[,],pts_per_batch=min(1000,nrow(xhp)),mu=0,sig=0.7,nugget=0.00001)

### perform non-stationary emulation, using futures for multicore (multisession) parallel computation ###
if(use_futures){
  ncores <- parallel::detectCores() - 2        # use total number of cores minus two
  plan(multisession, workers = ncores)         # set up parallel session using ncores
  GPem <- batch_future_lapply_Non_Stationary_emulator(x1=xh1,y1=y1,xp=xhp,pts_per_batch=min(ceiling(nrow(xhp)/ncores),nrow(xhp)),mu=0,sig=0.7,nugget=0.00001)
  plan(sequential)          # cancel multicore workers: return to single core usage
}

yp_mat <- matrix(GPem[[1]],nrow=length(xseq),ncol=length(yseq))			# matrix of emulator mean mu(x) for input points xp
ys_mat <- matrix(sqrt(GPem[[2]]),nrow=length(xseq),ncol=length(yseq))   # matrix of emulator sd sigma(x) for input points xp


### plot emulator mean ###
filled.contour(x=xseq,y=yseq,z=yp_mat,color.palette=terrain.colors,levels=pretty(range(yp_mat),20),xlab=axis_labels[1],ylab=axis_labels[2],
               plot.axes={axis(1);axis(2);plot_rift_points()
                 points(x1,pch=16)
               }
)

### plot emulator sd ###
filled.contour(x=xseq,y=yseq,z=ys_mat,color.palette=cm.colors,levels=pretty(c(0,0.7),14),xlab=axis_labels[1],ylab=axis_labels[2],
               plot.axes={axis(1);axis(2);plot_rift_points()
                 points(x1,pch=16)
               }
)

dev.off()


### End attempt proper fix for the projection compression issues for 2d fold in x1 and x2 simultaneously, two unequal rifts ###
########################################################################################################################################
















































