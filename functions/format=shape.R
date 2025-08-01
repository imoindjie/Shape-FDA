## This function takes a matrix of 2\times k dimensions: XX and associated ti of k-elements representing t \in [0,1] 
## The output is "shape" object, which has associated functions: plot, etc. 
as.shape=function(XX, ti){ 
  class(XX)='shape'
  attr(XX, 'time')=ti
  XX
}

### Takes a list of shape objects, obtained with the last function: shapes= list(shape_1, shape_2, ...)
## The output is an object "shapes", which has also some associated functions 
as.shapes=function(shapes){ 
  ti=attr(shapes[[1]], 'time')
  for(i in 1:length(shapes)){ 
    attr(shapes[[i]], 'time')=NULL
    attr(shapes[[i]], 'class')=NULL
  }
  class(shapes)='shapes'
  attr(shapes, 'time')=ti
  shapes
}

## for (a) shape(s) X, gives X^{(dim)}
get_dim=function(shape, dim){ 
  if(all(class(shape)=='shape') ){ 
    as.matrix(shape[dim, ])
  }else{ 
    matrix(unlist(lapply(1:length(shape), function(x)shape[[x]][dim, ])), nrow=length(shape), byrow = T )  
  }
}
## for shapes X, get the observations ind: X_i, i \in ind
get_ind=function(shapes, ind){ 
  if(length(ind)>1){ 
    res=as.shapes(lapply(ind, function(x)shapes[[x]]))
    attr(res, 'time')=attr(shapes, 'time')
    class(res)='shapes'
    res
  }else{ 
    out=shapes[[ind]]
    ti=attr(shapes, 'time')
    as.shape(out, ti)
  }
} 

## Associated plot functions for shapes and shape
plot.shapes<-function(ss, cols=NULL, types=NULL,pch=16, xlim=NULL, ylim=NULL, ...){ 
  xs=get_dim(ss, 1)
  ys=get_dim(ss, 2)
  if(is.null(types)) types='l'
  if(length(cols)==1 ){ 
    c_0=cols[[1]]
    for(i in 1:nrow(xs))cols=append(cols, list(c_0))
  }
  if(is.null(cols)) for(i in 1:length(ss))cols=append(cols, list(i))
  if(is.null(xlim))xlim=c(min(xs), max(xs))
  if(is.null(ylim))ylim=c(min(ys), max(ys))
  plot(NA, xlim=xlim, ylim=ylim, ...)
  for(i in 1:length(ss)){
    points(ss[[i]][1, ], ss[[i]][2, ], col=cols[[i]], type=types, pch=pch)
  }
}


plot.shape<-function(ss, cols=NULL, types=NULL,pch=16, xlim=NULL, ylim=NULL, arr=T, ...){ 
  xs=ss[1, ]
  ys=ss[2, ]
  if(is.null(types)) types='l'
  if(is.null(cols))cols=1
  if(is.null(xlim) & is.null(xlim)){ 
    xlim=c(min(xs), max(xs)); ylim=c(min(ys), max(ys))
    }  
  plot(xs, ys, xlim=xlim, ylim=ylim, col=cols, type=types, ...)
  if(arr) arrows(x0 = xs[1],y0 = ys[1], x1 = xs[2],y1 = ys[2], col='red', lwd=2) ## The parametrization starting point 
}

points.shape=function(ss, cols=NULL, types=NULL,pch=16, ...){ 
  xs=ss[1, ]
  ys=ss[2, ]
  if(is.null(types)) types='l'
  if(is.null(cols))cols=1
  
  points(xs, ys, xlim=c(min(xs), max(xs)), ylim=c(min(ys), max(ys)), col=cols, type=types, ...)
  #arrows(x0 = xs[1],y0 = ys[1], x1 = xs[2],y1 = ys[2], col='red', lwd=2) ## The parametrization starting point 
}
