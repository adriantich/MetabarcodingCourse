
std <- function(x) sd(x)/sqrt(length(x))

abund_rel <- function(x){
  if (sum(x)>0) {
    x <- x/sum(x)
  } else{
    x<-x
  }
}

abund_rel_hellinger <- function(x){
  if (sum(x)>0) {
    x <- sqrt(x/sum(x))
    x <- x/sum(x)
  } else{
    x<-x
  }
  return(x)
}
hellinger_transf <- function(x){
  if (sum(x)>0) {
    x <- sqrt(x/sum(x))
  } else{
    x<-x
  }
  return(x)
}

idem <- function(x){
  x <- x 
}

arrel_quarta <- function(x){
  x <- sqrt(sqrt(x))
}

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

sq <- function(x,a) {
  if (x==0) x<-0 else
    if (x<=a[1]) x<-1 else
      if (x<=a[2]) x<-2 else
        if (x<=a[3]) x<-3 else
          x<-4
}
semi_quant <- function(x){
  a <- unname(quantile(x[x>0], probs = c(.5,.75,.9)))
  x<- sapply(x, sq, a = a)
}

order_factor_levels <- function(x){
  if (!is.numeric(x)&!is.logical(x)) {
    x <- factor(x,levels = unique(x))
  } else {
    x
  }
}

order_factors <- function(x) { 
  as.data.frame(lapply(x, order_factor_levels))
}

in2mm <- function(x) { 
  x*25.4
}

mm2in <- function(x) { 
  x/25.4
}

points2size <- function(x) { 
  x/2.845276
}
