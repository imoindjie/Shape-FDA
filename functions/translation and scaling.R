## Convert shapes to functional data (MFPCA package)
to_fda=function(X_t, col){ 
  ti=attr(X_t, 'time')
  n=length(X_t)
  if(length(col)==1){ 
    funData(ti, t(simplify2array(lapply(1:n, function(x)X_t[[x]][col, ]))))
  }else{
    out=NULL
    for( cc in col){ 
      out=append(out, list(funData(ti,
                                   t(simplify2array(lapply(1:n, function(x)X_t[[x]][cc, ]))))) )
    }
    multiFunData(out)
  }
  
}

## From shapes get the translation vectors 
get_trans=function(shapes){ 
  n=length(shapes)
  T_1=unlist(lapply(1:n, function(x)integrate(to_fda(shapes, 1))[x]))
  T_2=unlist(lapply(1:n, function(x)integrate(to_fda(shapes, 2))[x]))
  
  cbind(T_1, T_2)
  
}

## From shapes get the scaling factors 
get_rho=function(shapes){ 
  X_1d=to_fda(shapes, 1:2)
  
  Rho=NULL
  for(i in 1:length(shapes)){ 
    Rho=c(Rho, 
          sqrt(scalarProduct(X_1d[i], X_1d[i])
          ) )
  }
  Rho
}

# Gives the X_star the untranslated and unscaled version of X: the preshape 
stand=function(X_obs, details=T){ 
  Tran=get_trans(X_obs)
  
  X_0=X_obs
  for(i in 1:length(X_obs)){ 
    X_0[[i]]=X_0[[i]]-matrix(Tran[i, ], ncol=1)[, rep(1, length(attr(X_obs, 'time')) )]
  }
  
  R=get_rho(X_0)
  for(i in 1:length(X_obs)) X_0[[i]]=1/R[i]*X_0[[i]]
  if(details){ 
    list(x_star=X_0, t_vec=Tran, s_vec=R)
  }else{ 
    X_0
      }

}
