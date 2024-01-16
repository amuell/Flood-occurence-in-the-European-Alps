#' This script produces an animation of combined spatial predictors f1+f2 across time
#' 
#' 
library(refund)
library(reshape2)
library(plotly)
library(tidyverse)

source("./user_functions.R") #some useful functions
load("./sediment_pffr.RData")
nn <- 300 # precision of bivariate plots (decrease this number if you don't have enough memory)
coef <- refund:::coef.pffr(pffr, n1 = nn, n2 = nn)

time <- coef$smterms$`Intercept(yindex)`$x
lon <- coef$smterms$`s(long,basis)`$x
lat <- coef$smterms$`s(lat,basis)`$x

A <- spatial.array(coef) 
colnames(A) <- lat; rownames(A) <- lon

B <- tibble(melt(A)) %>%
  mutate(time = time[Var3]) %>%
  rename(lon = Var1, lat = Var2) %>%
  select(-Var3) %>%
  filter(time < 2000)

nlev <- 100
lev <- unique(c( seq(min(B$value), 0, length.out = nlev/2),
                 seq(0, max(B$value), length.out = nlev/2)))
lev.sd <- (lev-min(lev))/(diff(range(lev)))

colour <-  c("purple","blue","lightblue","white","yellow","orange", "red")

fig <- B %>%
  plot_ly( x = ~lon,
           y = ~lat,
           z = ~value,
           frame = ~time,
           type = "contour",
           colors = colour,
           contours = list(contours = lev,
                           start = min(B$value) ,
                           end = max(B$value),
                           size = diff(range(B$value))/nlev)) %>%
  layout(title = TeX("$ \\widehat{f}_1(t,\\text{lon})+\\widehat{f}_2(t,\\text{lat})"),
         xaxis = list(title = "lon"),
         yaxis = list(title = "lat")) %>%
  animation_opts(frame = 300, easing = "quad") %>%
  animation_slider( currentvalue = list(prefix = "t = ", font = list(color = "black"))) %>%
  config(.Last.value, mathjax = 'cdn')

fig