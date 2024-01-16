#' This script procudes subfigures in Figure 9
library(refund)
library(reshape2)
library(tidyverse)
source("./user_functions.R") #some useful functions
load("./sediment_pffr.RData")
nn <- 300 # precision of bivariate plots (decrease this number if you don't have enough memory)
coef <- refund:::coef.pffr(pffr, n1 = nn, n2 = nn)

total.spatial.time.effect <- function(predictors = coef, time = 1, nlev = 100, contour=F, X = X){
  
  a <- predictors$smterms$`s(lat,basis)`$coef
  b <- predictors$smterms$`s(long,basis)`$coef
  
  lat <- unique(a$lat)  
  long <- unique(b$long)  
  year <- unique(a$yindex.vec)
  t <- which(year >= time)[1]
  if(is.na(t)){stop("choose a different time")}
  
  A <- spatial.array(predictors)
  
  At <- A[,,t]; colnames(At) <- long; rownames(At) <- lat
  
  df <- melt(At) %>%
    rename(lat = "Var1", lon = "Var2")
  
  df <- df[order(df$lon, df$lat, decreasing = c(F,F)),]
  
  year <- round(year[t])

  lev <- unique(c(seq(min(A), 0, length.out = nlev/2),
                  seq(0, max(A), length.out = nlev/2)))
  
  lev.sd <- (lev-min(lev))/(diff(range(lev)))

  color <-  c("purple","blue","lightblue","white","yellow","orange", "red")
  
  filled.contour(x = long,
                 y = lat,
                 z = A[ , , t],
                 levels = lev,
                 color.palette = colorRampPalette(color),
                 xlab = "long",
                 ylab = "lat",
                 main = paste("t = ",year, sep=""),
                 plot.axes={points(X$long,X$lat, pch = 16, cex = 0.75); axis(1); axis(2)})
}

total.spatial.time.effect(coef, time = 1100, nlev = 100, X = X)
p.1450 <- total.spatial.time.effect(coef, time = 1450, nlev = 100, X = X)
p.1700 <- total.spatial.time.effect(coef, time = 1700, nlev = 100, X = X)
p.1950 <- total.spatial.time.effect(coef, time = 1950, nlev = 100, X = X)

