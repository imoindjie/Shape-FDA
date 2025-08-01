## The Beta(delta) function with M fourier functions 
Beta=function(delta=0, M){ 
  I_M=1:M
  I_M=I_M[(I_M+1)%%2==0 ]
  B=NULL
  for(k in 1:length(I_M)){ 
    B_k=Ort(-pi*(I_M[k]+1)*delta )
    B=append(B, list(B_k))
  }
  as.matrix(Matrix::bdiag(B))
}

## The norm  \norm{X\circ \gamma_\delta- \mu}_H^2, using their scores X=alpha \phi and Mu= mu \phi
F_Opt=function(alpha, mu, delta, theta=0){ 
  norm(alpha%*%Beta(delta, ncol(alpha))-Ort(theta)%*%mu, type='F')^2
}


## The solution of the Procrustres: tan(theta)= R_1 
R_1=function(alpha, mu, delta){
  M=ncol(alpha)
  N_1= tcrossprod(alpha[2, ]%*%Beta(delta, M), mu[1, ]) -tcrossprod(alpha[1, ]%*%Beta(delta, M), mu[2, ])
  D_1= tcrossprod(alpha[1, ]%*%Beta(delta, M), mu[1, ]) +tcrossprod(alpha[2, ]%*%Beta(delta, M), mu[2, ])
  
  N_1/D_1
  
}

## Get the solution of delta and theta 
Optf=function(alpha, mu,s_deltas, info=F ){
  
  values=c()
  thetas=c()
  for(j in s_deltas){
    r_1=as.numeric(R_1(alpha, mu, j))
    t_j=atan(r_1)+c(0, pi)
    thetas=rbind(thetas,t_j )
    values=rbind(values, c(F_Opt(alpha, mu,j, t_j[1]), F_Opt(alpha, mu,j, t_j[2])) )
  }
  rownames(values)=s_deltas
  colnames(values)=c('theta_0', 'theta_1')
  
  t_star=thetas[which(values==min(values), arr.ind = T)][1]
  
  d_star=round(s_deltas[which(values==min(values), arr.ind = T)[1]], 2)
  
  if(info){ 
    list(delta=d_star, theta=t_star,final_loss=values, angles=thetas) 
  }else{ 
    list(delta=d_star, theta=t_star) 
  }
  
}
