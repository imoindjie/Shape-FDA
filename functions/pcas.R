## The Principal component analysis of shapes (inputs: the Fourier scores (get_alpha function), and ncomp (number of components))

pca.shape=function(alpha, ncomp=3){ 
  
  ti=attr(alpha$sm, "time")
  
  ## Principal component analysis by dimensions: based pca.fd of fda package 
  
  fpca1=pca.fd(alpha$fdobj[[1]]$fd, nharm=ncol(alpha$alpha_1) ) 
  fpca2=pca.fd(alpha$fdobj[[2]]$fd, nharm=ncol(alpha$alpha_2) ) 
  A=cbind(fpca1$scores, fpca2$scores)
  
  ## Principal component analysis of scores 
  
  U=eigen(var(A))
  
  S=U$vectors[, 1:ncomp]
  
  ## Eigen functions 
  
  Phi_1=fpca1$harmonics
  Phi_2=fpca2$harmonics
  Phi_1$coefs=fpca1$harmonics$coefs%*%S[1:ncol(fpca1$scores), ]
  Phi_2$coefs=fpca2$harmonics$coefs%*%S[-(1:ncol(fpca1$scores)), ]
  
  ## Shapes 
  m_Phi=list(eval.fd(ti, Phi_1 ), eval.fd(ti, Phi_2)) ## Functions to matrix 
  S_phi=NULL
  for(i in 1:ncomp){
    s_i=t(cbind(m_Phi[[1]][, i], m_Phi[[2]][, i]))
    attr(s_i, 'time')=attr(alpha$sm, 'time')
    class(s_i)='shape'
    S_phi=append(S_phi, list(s_i))
  }  
  
  S_phi=as.shapes(S_phi)
  
  mu_S=t(cbind(eval.fd(ti, fpca1$meanfd), eval.fd(ti, fpca2$meanfd)))
  attr(mu_S, 'time')=attr(alpha$sm, 'time')
  class(mu_S)='shape'
  
  prop=cumsum(U$values)/sum(U$values)
  
  Scores=A%*%U$vectors[, 1:ncomp]
  
  res=list(var_exp=prop[1:ncomp], f_propres=list(shape=S_phi, fd=list(Phi_1, Phi_2) ), mean=mu_S, Scores=Scores, pcaobj=list(fpca1, fpca2, U))
  class(res)='f_shape'
  res
}

## predict the PCA scores of new shapes newx

predict.f_shape=function(pca, newx){ 
  ti=attr(pca$f_propres$shape, 'time')
  Us=list(pca$pcaobj[[1]]$harmonics$coefs, pca$pcaobj[[2]]$harmonics$coefs)
  Ms=list(pca$pcaobj[[1]]$meanfd$coefs, pca$pcaobj[[2]]$meanfd$coefs)
  Bs=list(pca$pcaobj[[1]]$harmonics$basis, pca$pcaobj[[2]]$harmonics$basis)
  
  coefs=lapply(1:2, function(x)smooth.basis(ti, t(get_dim(newx, x)), Bs[[x]])$fd$coefs-Ms[[x]][, rep(1, length(new))])  
  Scores=cbind(t(coefs[[1]])%*%Us[[1]],t(coefs[[2]])%*%Us[[2]])%*%pca$pcaobj[[3]]$vectors
  Scores[, 1:length(pca$f_propres)]
}

## From scores, generate new curves using pca objects 

generate=function(pca, scores){ 
  
  if(length(pca$f_propres[[1]])==ncol(scores)){ 
    ti=attr(pca$f_propres[[1]], 'time')
    shapes=NULL
    for(i in 1:nrow(scores)){ 
      shape_i=pca$mean
      for(k in 1:ncol(scores)){
        shape_i=shape_i+scores[i, k]*pca$f_propres[[1]][[k]]
      }
      
      attr(shape_i, 'time')=ti
      class(shape_i)='shape'
      shapes=append(shapes, list(shape_i))
      
    }
    if(length(shapes)==1){ 
      shape_i
    }else{
      as.shapes(shapes)
      
    }
    
  }
}

