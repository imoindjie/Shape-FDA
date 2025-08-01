# Set working directory to file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

###################################################################################
### Part 1: Import and visualization

n_sh="butterfly" ## planar curves names ## allowed names in list.files('data')

source('functions/packages.R')

X_obs=readRDS(paste0('data', '/', n_sh, '.RDS') ) ## load the n_sh collections (from '/data')

source('functions/format=shape.R') ## Set of functions (as.shape(s), plot, etc) for the format "shape" and "shapes"

plot(X_obs, xlab=expression(X^{(1)}), ylab=expression(X^{(2)})) 

source('functions/fda for curves.R') ## Load fda related functions: approximation of fourier functions, plot of the dimension of X 

## M Effect of the number of basis functions for smoothing curves
if(T){ 
  for(M in 2*seq(1, 20) ){ 
    plot(get_alpha(get_ind(X_obs, rep(1, 2) ), M)$sm, main=paste0('M=', M, ' Fourier functions'), xlim=c(-2, 2), ylim=c(-2, 2))
    Sys.sleep(1)
  }
}


###################################################################################
### Part 2: Alignment


source('functions/translation and scaling.R') # Set of functions for estimating translation and scaling

std=stand(X_obs) # get the unstranslated and unscaled version of X
plot(std$t_vec, pch=16, cex=std$s_vec) ## The translation and the scaling factors
X_star=std$x_star 
plot(X_star) ## The un-scaled and un-translated observations 


## Coordinate functions 
plot(to_fda(X_obs, 1:2), xlab='t', ylab=c(expression(X^{(1)}), expression(X^{(2)})) )
plot(to_fda(X_star, 1:2), xlab='t', ylab=c(expression(X^{(1)}), expression(X^{(2)})) )


## Transformations and optimization 
source('functions/Transformations.R') ## For operation of shape variables 
source('functions/optimization_deltas and thetas.R') ## Functions for estimating delta and thetas

out=get_alpha(X_star, M=51) ## Approximation into Fourier functions 
u=rbind(out$alpha_1[1, -1], 
            out$alpha_2[1, -1]) # The coefficient of the template 

hat_d=hat_t=NULL ## Estimation of delta (hat_d) and theta (hat_t)

S_0=seq(0, 1, by=0.01) ## Considered grid for estimating delta 
for(i in 1:length(X_star)){
  message(paste0('Image ',i, ' out of ', length(X_star) ) )
  
  alpha=rbind(out$alpha_1[i, -1], 
              out$alpha_2[i, -1]) ## The score of the i-th individual 
  
  Oout=Optf(alpha, u, S_0) ## Optimization 
  
  hat_d=c(hat_d, Oout$delta)
  hat_t=c(hat_t, Oout$theta)
  
}

## The alignment  
X_tilde=X_star

for(i in 1:length(X_star)){
  X_i=as.shape(X_star[[i]], attr(X_star, 'time')) ## extract the i-th individual from X_star
  X_tilde[[i]]=rot(param(X_i, hat_d[i] ), -hat_t[i]) 
}
## From list of shape to shapes 

X_tilde=as.shapes(X_tilde)

plot(X_tilde, xlab= expression(X^{(1)}),ylab=expression(X^{(2)})) ## Plot the planar curves 

plot(to_fda(X_tilde, 1:2), xlab='t', ylab=list(expression(X^{(1)}),expression(X^{(2)})) ) ## Plot the coordinate functions 


####################################################
#### Part 3: Principal component analysis 

Z_2=cbind(std$t_vec, map(cbind(hat_d, -hat_t))) ## The deformation matrix: Translation, MAP(delta, theta). MAP is a bijective maping from [0, 1] to [0, 2\pi] to \mathbf{R}^2  

#### The Functional part Z_1=rho \tilde{X}

Z_1=X_tilde
for(i in 1:length(X_star)){ 
  Z_1[[i]]=std$s_vec[i]*X_tilde[[i]]
}

source('functions/pcas.R') ## Functions for shape based PCAs

## PCA of forms: Perform MFPCA on Z_1 

np=length(Z_1)-1 ## Maximum number of components 

G=pca.shape( get_alpha(Z_1, M = 51),ncomp=np) ## The pca object

np=min(which(G$var_exp>=.9)) ## Keep the components explaining 90\% of totat variance

plot(G$var_exp, type='b', pch=16, xlab='components', ylab='explained variance')
abline(v=np, col='red')

G=pca.shape( get_alpha(Z_1, M = 51),ncomp=np) ## Estimate Final PCA with np explains 90 %

## Figs of the first 5 PCs and associated variations mean-0.5*phi_k, where phi_k=get_ind(G$f_propres$shape, k) ) is the k-th eigen function

for(k in 1:np){ 
  plot(G$mean, xlab=expression(X^{(1)}), ylab=expression(X^{(2)}), xlim=c(-0.7, 0.7), ylim=c(-0.7, 0.7), arr=F, main=paste0('Comp: ' , k) )
  points(add(G$mean, 0.05*get_ind(G$f_propres$shape, k) ), xlab=expression(X^{(1)}), ylab=expression(X^{(2)}), pch=16, col='blue', lty=2)
  points(add(G$mean, -0.05*get_ind(G$f_propres$shape, k) ), xlab=expression(X^{(1)}), ylab=expression(X^{(2)}), pch=16, col='red',lty=2)
  Sys.sleep(2)
  dev.off()
}


plot(G$Scores) ## Present the scores 

## Joint PCA : Z_1 and Z_2 (shapes and deformations)

Xi_all=cbind(G$Scores, Z_2 ) ## Scores of Shapes PCAs and Z_2 (the deformation variables)

Sigma=var(Xi_all) ## Estimate the covariance matrix of scores and Z_2

## Generate new observations X_gen 

n_gen=10 ## Number of generated curves

Xi_gen=mvrnorm(n_gen,apply(Xi_all, mean, MARGIN=2), Sigma = Sigma) ## Generation of samples using the covariance matrice (Sigma) of the scores

Z_gen_1=generate(G, Xi_gen[, 1:ncol(G$Scores)]) ## Generated the curves 

## Generated the deformations 

Z_gen_2=Xi_gen[, -(1:ncol(G$Scores)) ]
Def=cbind(Z_gen_2[, 1:2], mapinv(Z_gen_2[, 3:4]))

##  Generate the contours 
X_gen=NULL
for(i in 1:nrow(Z_gen_2)){ 
  x_i=as.shape(Z_gen_1[[i]], attr(Z_gen_1, 'time'))
  x_i=tran(rot(param(x_i, Def[i,3]),Def[i,4]), Def[i, 1:2])
  X_gen=append(X_gen, list(x_i))
  }

X_gen=as.shapes(X_gen) 

## and plot them

for(k in 1:n_gen){ 
  plot(get_ind(X_gen, k), xlab=expression(X^{(1)}), ylab=expression(X^{(2)}), xlim=c(-2, 2), ylim=c(-2, 2), main=paste0('Ind: ' , k) ) ## With deformations 
  Sys.sleep(2)
  dev.off()
  }
