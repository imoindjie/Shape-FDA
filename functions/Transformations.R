## For angle in [0, 2\pi], gives the Orthogonal matrix 
Ort=function(angle){ 
  if(class(angle)=='numeric')class(angle)='rotation'
  
  if(class(angle)=='rotation'){ 
    matrix(c(cos(angle),sin(angle), -sin(angle), cos(angle)), nrow=2)
  }else{ 
    matrix(c(cos(angle),sin(angle), sin(angle), -cos(angle)), nrow=2)
    }
  
}

## rotate shape with angle
rot=function(shape, angle){ 
  res=Ort(angle)%*%shape
  class(res)='shape'
  attr(res, 'time')=attr(shape, 'time')
  res
}
## translate shape 
tran=function(shape, vec){ 
  res=shape+matrix(rep(vec, length(attr(shape, 'time'))), nrow=2)
  class(res)='shape'
  attr(res, 'time')=attr(shape, 'time')
  res
}
## substrate shape_1 and shape_2: shape_1-shape_2
diff.shape=function(shape_1, shape_2){ 
  res=shape_1
  if(class(shape_1)==class(shape_2)){ 
    
    for(i in 1:length(shape_1)){ 
      
      res[[i]]=shape_1[[i]]-shape_2[[i]]
      attr(res[[i]], 'time')=NULL
      attr(res[[i]], 'class')=NULL
    }
    class(res)=class(shape_1)
    
  }else{ 
    if(class(shape_2)=='shape' ){ 
      for(i in 1:length(shape_1)){ 
        res[[i]]=shape_1[[i]]-shape_2
        attr(res[[i]], 'time')=NULL
      }
      res=as.shapes(res)
    }else{ 
      message('Erreur')  
      }  
    
    }
  res
}
## summation of shapes: shape_1+shape_2  
add=function(shape_1, shape_2){ 
  diff(shape_1, - shape_2)
}

## multiply the shape with rho: rho*shape 
scal=function(shape, rho){ 
  res=rho*shape
  class(res)='shape'
  attr(res, 'time')=attr(shape, 'time')
  res
}
## scalar product of shapes: shape_1 and shape_2, inner_product(shape_1, shape_2)
prod_scal=function(shape_1, shape_2){ 
  if(length(shape_1)==1){ 
    t(get_dim(shape_1, 1))%*%get_dim(shape_2, 1)+t(get_dim(shape_1, 2))%*%get_dim(shape_2, 1)
  }else{ 
    get_dim(shape_1, 1)%*%get_dim(shape_2, 1)+get_dim(shape_1, 2)%*%get_dim(shape_2, 2)
    }
  
  }
## Linear geodesics of two vectors score_1 and score_2: x*score_1+(1-x)*score_2
line=function(score_1, score_2, len=100){ 
  
  lambda=seq(0, 1, by=1/len)
  rlist::list.rbind(lapply(lambda, function(x)x*score_1+(1-x)*score_2))
  
  
}

## Parametrize a shape: shape \circ \gamma_\delta
param=function(shape, delta){ 
  if(delta<0) delta=h(delta)
  ti=attr(shape, 'time')
  if(is.null(ti)) ti=seq(0, 1, length.out=ncol(shape))
  xi=which.min(abs(delta-ti))
  ind_i=c(xi:(length(ti)-1),1:(xi) ) 
  
  res=shape
  if(!(xi %in%  c(1,length(ti)))  ){ 
    for(i in 1:nrow(res)){ 
      res[i, ]= res[i, ind_i]
    }
  }
  res
}

## Transformation of the map:[0, 2\pi] \times [0, 1] \to \mathbf{R}^2
map=function(x)cbind(tan(pi/2*(x[, 1]-0.5)) ,tan(1/4*(0.5*x[, 2]-pi) ))  


## Inverse transformation of map, map^{-1}: \mathbf{R}^2 \to [0, 2\pi] \times [0, 1]
mapinv=function(x){ 
  res=cbind(atan(x[, 1])*2/pi+0.5, (4*atan(x[, 2])+pi)*2)
  for(i in 1:nrow(res)){ 
    if( res[i, 1]<0 |res[i, 1]>1){ 
      res[i, 1]=res[i, 1]%%1
    }
    if(res[i, 2]<0){ 
      res[i, 2]=res[i, 2]+2*pi
    }
    
  }
  res
  
}

## Transform angle in R to [0, 2\pi]
good_angle=function(x){ 
  ifelse(x<0, x+2*pi, ifelse(x>2*pi,x-2*pi,  x) )
}

## The modulo 1 function: h(ti)=mod(ti, 1)
h=function(ti)ti%%1

## For the shapes, get the translation vectors

get_trans=function(shapes){ 
  n=length(shapes)
  T_1=unlist(lapply(1:n, function(x)integrate(to_fda(shapes, 1))[x]))
  T_2=unlist(lapply(1:n, function(x)integrate(to_fda(shapes, 2))[x]))
  
  cbind(T_1, T_2)
  
}

## For the shapes, get the the norm 

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

## For the shapes, get  the norm

get_scale=function(X_fin, m_def){ 
  n=length(X_fin)
  X_mean=get_ind(X_fin, 1)
  for( i in 2:n){ 
    X_mean=X_mean+X_fin[[i]]
  }
  X_mean=1/n*X_mean
  X_s=X_fin
  for( i in 1:n){ 
    X_s[[i]]=X_s[[i]]-X_mean
  }
  num=sum(get_rho(X_s)^2)
  
  dd=scale(m_def, center=T, F)
  num/norm(t(dd)%*%dd, type='F')^2
}


## Transform X curves to preshapes

untr=function(X_obs, scale=F){ 
  Tran=get_trans(X_obs)
  
  X_0=X_obs
  for(i in 1:length(X_obs)){ 
    X_0[[i]]=X_0[[i]]-matrix(Tran[i, ], ncol=1)[, rep(1, length(attr(X_obs, 'time')) )]
  }
  if(scale){ 
    R=get_rho(X_0)
    for(i in 1:length(X_obs)) X_0[[i]]=1/R[i]*X_0[[i]]
    }
 
  X_0
}

