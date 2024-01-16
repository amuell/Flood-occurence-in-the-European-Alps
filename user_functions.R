#' Bandwidth smoothing
#' 
#' @param data Matrix of observations, can have missing values.
#' @param h Bandwidth for the smoothing.
#' 
#' @return Returns smoothed curves (adapted moving average estimates).
#' 
bin.smooth <- function(data,h=40,na.replace = F, min.h = 20){
  n <- nrow(data)
  d <- ncol(data)
  s <- matrix(NA,nrow = n, ncol = d)
  
  for(j in 1:d){
    obs <- !is.na(data[,j]) # observed points of column j
    vec <- data[obs,j] # observations for lake j
    n.j <- length(vec) # number of observations for lake j
    s.j <- numeric(n.j) # vector of size n.j
    
    for(i in 1:n.j){
      y <- 0
      
      if(i< min.h){
        y <- vec[max(i-min.h+1,1):(i+min.h-1)]
        s.j[i] <- sum(y)/length(y)
        
      }else if(i>(n.j-min.h)){
        y <- vec[min(i-min.h+1,n.j-min.h):n.j]
        s.j[i] <- mean(y)
      }else{
        
        y <- vec[max(i-h+1,1):min(i+h-1,n.j)]
        s.j[i] <- mean(y)
        
        #s.j[i] <- sum(vec[(i-h):n.j])/(2*h)
      } 
    }
    s[obs,j] <- s.j
  }
  
  if(na.replace){
    for(i in 1:ncol(s)){
      if(any(is.na(s[,i]))){ s[is.na(s[,i]),i] <- mean(s[,i],na.rm=T) }
    }
  }
  
  s
}

#' Produces an array with all evaluations of f1+f2 (combined spatial predictors)
#' 
#' @param coef Estimated predictors. 
#' 
#' @return Returns an array of size nn^3.
#'
spatial.array <- function(coef){
  n <- nrow(coef$smterms$`Intercept(yindex)`$coef)
  
  
  B <- array(NA, dim = c(n,n,n))
  
  for(t in 1:n){
    a <- coef$smterms$`s(long,basis)`$coef %>% 
      filter(yindex.vec == unique(yindex.vec)[t])
    
    b <- coef$smterms$`s(lat,basis)`$coef %>% 
      filter(yindex.vec == unique(yindex.vec)[t])
    
    c <- coef$smterms$`Intercept(yindex)`$coef$value[t]
    
    A <- outer(X = a$value, Y = b$value, FUN = "+") 
    B[,,t] <- A
  }
  B
}