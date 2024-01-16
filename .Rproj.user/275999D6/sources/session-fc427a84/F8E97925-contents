#' This scripts performs the leave one out validation and produces Figures 13 and 14


library(mgcv)
library(refund)
library(readr)
library(latex2exp)
library(RColorBrewer)
library(reshape2)
library(mgcViz)
library(tidyverse)
rm(list=ls())

load("./sediment_pffr.RData")

#' Compute the overall out of sample prediction error
#'
#' @param data Functional sample 
#' @param data.binned Two sided moving average estimator
#' @param Xcov Data frame of covariates
#' @param yind Grid of time arguments
#' @param kint Number of basis functions for regression terms that vary in the same index as the functional intercept.
#' @param kbiv Number of functions for marginal basis that do not vary in time 
#'
#' @returns A matrix of out of sample predictions
#'
out_sample_error <- function(data = data.s,
                             data.binned = data.bin,
                             Xcov = X,
                             yind = year.s){
  # Out of sample
  err.mat <- matrix(NA,ncol = ncol(data.binned),nrow = nrow(data.binned))
  
  for(i in 1:ncol(data)){
    #data without lake i
    data.sv <- data[,-i]
    X1v <- Xcov[-i,] 

    rownames(X1v) <- 1:ncol(data.sv)
    
    ydata.i <- melt(as.matrix(data.sv))%>%
      mutate(.obs = as.numeric(Var2))%>%
      select(.obs, .index = Var1, .value = value)%>%
      filter(!is.na(.value))
    
    # check if all id's corresponds
    if(!all(ydata.i$.obs %in% rownames(X1v))){stop("issue with identifiability of covariates")}
    
    #Regression formula
    f <- Y ~  c(s(ca, bs = "cr")) +
      s(long, bs = "cr") +
      s(lat, bs = "cr") +
      c(s(id, bs = "re"))
    
    yindex.i <- unique(sort(ydata.i$.index)) # time arguments for the regression
    
    # fix knots positions to quantiles of data (if there are large places of the domain where we have 0 data points)
    knot.list <- list(ca = quantile(X1v$ca,probs = seq(0, 1,length.out = 10)))
    
    pffr <- pffr(formula = f,
                 yind = yindex.i,
                 data = X1v,
                 ydata = ydata.i,
                 method = "REML",
                 knots = knot.list,
                 bs.int = list(bs = "cr", k = 20, m = c(2, 1)),
                 bs.yindex = list(bs = "cr", k = 5, m = c(2, 1))
                 );gc()
    
    # Predict the missing lake
    new.data <- Xcov[i,] 
        
    pffr1 <- pffr
    pffr1$pffr$nyindex <- length(yindex.i)
    pffr1$pffr$yind <- as.numeric(year.s)
    
    pred <- 1/(1+exp(-t(predict(pffr1, type="response", exclude = "c(s(id))",
                                newdata = new.data, newdata.guaranteed = T))))
    
    err.mat[,i] <- 1/(1+exp(-data[,i])) - pred
    
    print(paste("Out of sample:",i))
  }
  err.mat
}

mat <- out_sample_error(data = data.s,
                        data.binned = data.bin,
                        Xcov = X,
                        yind = year.s)

colnames(mat) <- names(data1)

mean(colSums(mat^2, na.rm = T))

# plot prediction errors
library(ggplot2)
library(reshape2) 
library(tidyverse)
#source("./theme_publication.R")
exception <- c("OES", "LED", "LDB", "AMM","MON", "ALO", "BAR")

df <- melt(mat) %>%
  select(t = Var1, lake = Var2, error = value) %>%
  filter(!is.na(error)) %>% mutate(t = t)

df.exception <- df %>% filter(lake %in% exception)
theme_set(theme_bw()+ 
  theme(text = element_text(size=11)) +
  theme(legend.title.align = 0.5) +
  theme(axis.text=element_text(size=rel(1.2))))

error <- ggplot()+ 
  geom_point(data = df, aes(t, error), shape = 20, size = 1) +
  geom_point(data = df.exception, aes(t, error, color = lake), shape = 20, size = 1.5) +
  #geom_label(data = filter(df, t == 1500),aes(label = lake), nudge_x = 0.35, size = 4, show.legend = F) +
  guides(color=guide_legend("Lake:"))+
  ylab(TeX("$\\widehat{\\gamma}_{i}(t) - \\widehat{\\gamma}_{-i}(t, X_{r,-i})"))+
  xlim(1000,2000) +
  xlab(TeX("t"))

error

#' plot errors over time
#' Going towards red means that the GFAMM is over-estimating the probability
#' The colour green is the reference for a difference close to 0
#' Going towards purple means that the GFAMM is under-estimating the probability

#' Produces a map with the out of sample predictions errors
#' 
#' @param df A data.frame that provides informations on latitude, longitude and predictions error, for a specific time year t.
#' 
error.map <- function(df,scale = 100, size = 30, nlev = 100, bins = 11, min.val, max.val){
  library(leaflet)
  lonr <- range(df$long)
  latr <- range(df$lat)
  
  color.scale <- rev(c( "darkred","red","orange", "yellow", "green","lightblue", "aquamarine","blue", "darkblue"))

  scale.lev <- unique(c(seq(min.val,0, length.out = nlev/2),
                        seq(0, max.val, length.out = nlev/2)))
  
  pal <- colorNumeric(color.scale, scale.lev)
  
  color <- pal(df$val)

  m <- leaflet(data = df,
               padding = 1) %>%
    fitBounds(lng1 = lonr[1],
              lng2 = lonr[2],
              lat1 = latr[1],
              lat2 = latr[2]
    ) %>%
    addCircleMarkers(~long,
                     ~lat,
                     radius = size,
                     label = ~round(val,4),
                     stroke = F,
                     fillOpacity = 0.5,
                     labelOptions = labelOptions(permanent = F, direction = "bottom") ,
                     color=color) %>%
    addLabelOnlyMarkers(~long,
                        ~lat,
                        label = ~id,
                        labelOptions = labelOptions(noHide = T,
                                                    opacity = 1,
                                                    textsize = "20px",
                                                    textOnly = T,
                                                    direction = "center",
                                                    style = list(
                                                      "color" = "black")
                        )) %>%
    addProviderTiles(providers$Esri.WorldTerrain) %>%
    addLegend(position = "bottomright",
              pal = pal,
              values = scale.lev,
              bins = 10,
              title = "Prediction error",
              opacity = 1)
  
  return(m)
}

# plot errors for a specific year t
error.map.time <- function(data = mat, t = 1, W = X, scale = 30, size = 50, nlev = 101){
  df <- W %>% mutate(val = data[t,]) %>% filter(!is.na(val))
  
  min.val <- -1
  max.val <- 1
  
  error.map(df, scale, size, nlev, min.val = min.val, max.val = max.val)
}

library(mapview)
size = 30
m1 <- error.map.time(data = mat, t = 100, W = X, size)
m2 <- error.map.time(data = mat, t = 300, W = X, size)
m3 <- error.map.time(data = mat, t = 500, W = X, size)
m4 <- error.map.time(data = mat, t = 700, W = X, size)
m5 <- error.map.time(data = mat, t = 900, W = X, size)
m6 <- error.map.time(data = mat, t = 1000, W = X, size)

m1
m2
m3
m4
m5
m6

# # Uncomment these lines if you want to save plots
# width <- 1504
# height <- 1155
# 
# mapshot(x = m1, file = "./plot/error_h15_t1100.png", vwidth = width, vheight = height)
# mapshot(x = m2, file = "./plot/error_h15_t1300.png", vwidth = width, vheight = height)
# mapshot(x = m3, file = "./plot/error_h15_t1500.png", vwidth = width, vheight = height)
# mapshot(x = m4, file = "./plot/error_h15_t1700.png", vwidth = width, vheight = height)
# mapshot(x = m5, file = "./plot/error_h15_t1900.png", vwidth = width, vheight = height)
# mapshot(x = m6, file = "./plot/error_h15_t2000.png", vwidth = width, vheight = height)

