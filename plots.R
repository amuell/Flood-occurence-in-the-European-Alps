#' This script Produces all the remaining plots needed in the paper
library(refund)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(latex2exp)
library(reshape2)
library(tidyverse)

load("sediment_pffr.RData") #load all necessary data

nn <- 200 # precision of plots
coef <- refund:::coef.pffr(pffr, n1 = nn, n2 = nn) # extract estimated predictors 

user.theme <- theme_bw()+ 
  theme(text = element_text(size=11),
        legend.title.align = 0.5,
        axis.text = element_text(size=rel(1.2)),
        legend.text = element_text(size = 10))

#' Diagnostic plot
#' 
#' @param text Should Lake names be printed as well?
#' 
#' @return Creates Figure 12 
fpc.plot <- function(fpca.res,
                     fpca.fit,
                     gg.theme = user.theme,
                     lab = lakes,
                     text = F){
  
  rownames(fpca.res$scores) <- rownames(fpca.fit$scores) <- lab
  
  A <- melt(fpca.res$scores, value.name = "residual") %>%
    mutate(Var1 = as.factor(Var1)) ; names(A)[1:2] <- c("lake", "resFpc")
  
  B <- melt(fpca.fit$score, value.name = "fitted")%>%
    mutate(Var1 = as.factor(Var1)); names(B)[1:2] <- c("lake", "fitFpc")
  
  df <- A %>%
    full_join(y = B, by  = join_by(lake), relationship = "many-to-many")
  
  theme_set(gg.theme)
  
  g <- ggplot(df, aes( x= fitted, y = residual)) +
    geom_point(aes(color = lake), size = 2, shape = 20) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "black") +
    facet_grid(resFpc ~ fitFpc, scales = "free", shrink = F) +
    ylab(TeX("$ \\widehat{\\xi}_{ik}^R"))+
    xlab(TeX("$ \\widehat{\\xi}_{ik}^f"))+
    guides(color=guide_legend("Lake:"))
  
  if(text){
    g <- g + geom_text(aes(label = lake), hjust = 0, nudge_x = 0.075, size = 1.5)
  }
  
  g
  
}

#' Produce all estimated predictors plots in the paper
#' 
#' @param nlev Number of levels used for contour plots.
#' @param qn Quantile used for confidence interval plots.
#' @param gg.theme User defined theme for ggplot plots.
#' @param width.biv Width of bivariate figures (for saving plots).
#' @param height.biv Height of bivariate figures (for saving plots).
#' @param save Boolean. Should figures be saved?
#' @param path.plot Path to save plots.
#' @param width.pdf Width of other figures (for saving plots).
#' @param height.pdf Height of other figures (for saving plots).
#' @param unit Dimension unit for saving plots.
#' @param res Resolution of saved figures.
#' 
#' @return Prints all figures of estimated predictors present in the paper. They can also be saved in a local directory.
#' 
pffr.paper.plot <- function(coef,
                       nlev = 50,
                       qn = qnorm(0.975),
                       gg.theme = user.theme,
                       width.biv = 10,
                       height.biv = 4,
                       save = F,
                       path.plot = "./plot",
                       width.pdf = 7,
                       height.pdf = 4,
                       unit = "in",
                       res = 300){

  name <- names(coef$smterms)
  
  # plot the functional intercept (Figure 4)
  df.int <- coef$smterms$`Intercept(yindex)`$coef %>%
    mutate(value = value +coef$pterms[1], year = yindex.vec) %>%
    mutate(up = value + qn*se,
           low = value - qn*se) 
  
  theme_set(gg.theme)
  
  int.plot <- ggplot( df.int, aes(x = year)) + 
    geom_line(aes(y = value), color = "black") +
    geom_line(aes(y = up), color = "red", linetype = "dashed") +
    geom_line(aes(y = low), color = "red", linetype = "dashed") +
    ylab(TeX("$ \\widehat{\\beta}_0 + \\widehat{\\beta}(t)"))+
    xlim(1000,2000) +
    xlab(TeX("$ t"))  
  
  print(int.plot)

  # plot the effect of the catchment area (Figure 6)
  df.ca <- coef$smterms$`c(s(ca,basis))`$coef %>%
    mutate(up = value + qn*se,
           low = value - qn*se)

  theme_set(gg.theme)
  ca.plot <- ggplot( df.ca, aes(x = ca)) +
    geom_line(aes(y = value), color = "black") +
    geom_line(aes(y = up), color = "red", linetype = "dashed") +
    geom_line(aes(y = low), color = "red", linetype = "dashed") +
    ylab(TeX("$ \\widehat{f}_3(\\log(ca_i))"))+
    xlab(TeX("$ \\log(ca_i)"))

  print(ca.plot)
  
  plot.fpc <- fpc.plot(fpca.res, fpca.fit, text = T)
  print(plot.fpc)
  
  # make the bivariate plots together
    
  df1 <- coef$smterms$`s(long,basis)`$coef %>% mutate( up = value + qn*se,
                                                       low = value - qn*se,
                                                       sig = ifelse(up<0 | low>0, 1, NA)) 
  
  df2 <- coef$smterms$`s(lat,basis)`$coef %>% mutate( up = value + qn*se,
                                                      low = value - qn*se,
                                                      sig = ifelse(up<0 | low>0, 1, NA))
  
  min_value <- min(df1$value, df2$value)
  max_value <- max(df1$value, df2$value)
      
  # set contour levels such that white is the colour used for null values
  lev <- unique(c( seq(min_value, 0, length.out = nlev/2),
                   seq(0, max_value, length.out = nlev/2)))
  
  lev.sd <- (lev-min(lev))/(diff(range(lev)))
  
  color <-  c("purple","blue","lightblue","white","yellow","orange", "red")
      
  theme_set(gg.theme+ theme(legend.key.height= unit(2.5, units = "lines")))
  
  # Create your first plot (plot1)
  plot1 <- ggplot(df1, aes(x = yindex.vec, y = long, fill = value)) + 
    ylab(TeX("$ \\widehat{f}_1(t, long_i)"))+
    xlab(TeX("$ t")) +  
    xlim(1000,2000) +
    geom_tile() +
    scale_fill_gradientn(limits = c(min_value, max_value), colors = color, values = lev.sd) 
  
    # Create your second plot (plot2)
  plot2 <- ggplot(df2, aes(x = yindex.vec, y = lat, fill = value)) + 
    ylab(TeX("$ \\widehat{f}_2(t, lat_i)"))+
    xlab(TeX("$ t")) +
    geom_tile() +
    xlim(1000,2000) +
    scale_fill_gradientn(limits = c(min_value, max_value), colors = color, values = lev.sd)

  # plot both spatial predictors (Figure 7)
  p3 <- ggarrange(plot1, plot2,align = "hv", nrow = 1, ncol = 2, common.legend = T, legend = "right",
                  font.label = list(size = 11, color = "black", face = "bold", family = NULL),
                  labels = list("(A)","(B)"))
  
  print(p3)

  # Plot confidence interval of predictors f1 and f2
  lev.long <- unique(c( seq(min(df1$low), 0, length.out = nlev/2),
                        seq(0, max(df1$up), length.out = nlev/2)))
    
  lev.long.sd <- (lev.long-min(lev.long))/(diff(range(lev.long)))
    
  p.long <- ggplot(filter(df1, !is.na(sig)), aes(x = yindex.vec, y = long, fill = value)) + 
    ylab(TeX("$ \\widehat{f}_1(t, long_i)"))+
    xlab("t") +    xlim(1000,2000) +
    geom_tile() +
    scale_fill_gradientn(limits = c(min(df1$low), max(df1$up)), colors = color, values = lev.long.sd) 
    
  p.long.up <- ggplot(df1, aes(x = yindex.vec, y = long, fill = up)) + 
    ylab("upper interval")+
    xlab("t") +  
    xlim(1000,2000) +
    geom_tile() +
    scale_fill_gradientn(limits = c(min(df1$low), max(df1$up)), colors = color, values = lev.long.sd) 
  
  p.long.low <- ggplot(df1, aes(x = yindex.vec, y = long, fill = low)) + 
    ylab("lower interval")+
    xlab("t") +  
    xlim(1000,2000) +
    geom_tile() +
    scale_fill_gradientn(limits = c(min(df1$low), max(df1$up)), colors = color, values = lev.long.sd) 
  
  # Figure 8a  
  p.long.ci <- ggarrange(p.long.low, p.long,  p.long.up, align = "hv", nrow = 1, ncol = 3, common.legend = T, legend = "right",
                         font.label = list(size = 11, color = "black", face = "bold", family = NULL),
                         legend.grob = get_legend(p.long, position = "right"))
  
  lev.lat <- unique(c( seq(min(df2$low), 0, length.out = nlev/2),
                       seq(0, max(df2$up), length.out = nlev/2)))
  
  lev.lat.sd <- (lev.lat-min(lev.lat))/(diff(range(lev.lat)))
  
  p.lat <- ggplot(filter(df2, !is.na(sig)), aes(x = yindex.vec, y = lat, fill = value)) + 
    ylab(TeX("$ \\widehat{f}_2(t, lat_i)"))+
    xlab("t") +  
    xlim(1000,2000) +
    geom_tile() +
    scale_fill_gradientn(limits = c(min(df2$low), max(df2$up)), colors = color, values = lev.lat.sd) 
    
  p.lat.up <- ggplot(df2, aes(x = yindex.vec, y = lat, fill = up)) + 
    ylab("upper interval")+
    xlab("t") +  
    xlim(1000,2000) +
    geom_tile() +
    scale_fill_gradientn(limits = c(min(df2$low), max(df2$up)), colors = color, values = lev.lat.sd) 
    
  p.lat.low <- ggplot(df2, aes(x = yindex.vec, y = lat, fill = low)) + 
    ylab("lower interval")+
    xlab("t") +  
    xlim(1000,2000) +
    geom_tile() +
    scale_fill_gradientn(limits = c(min(df2$low), max(df2$up)), colors = color, values = lev.lat.sd) 
  
  # Figure 8b  
  p.lat.ci <- ggarrange(p.lat.low, p.lat,  p.lat.up, align = "hv", nrow = 1, ncol = 3, common.legend = T, legend = "right",
                        font.label = list(size = 11, color = "black", face = "bold", family = NULL),
                        legend.grob = get_legend(p.lat, position = "right"))
  print(p.long.ci)
  print(p.lat.ci)
  
  # plot random effect (Figure 5)

   df.re <- unique(coef$smterms$`c(s(id))`$coef) %>%
    mutate(up = value + qn*se,
           low = value - qn*se) 
  
  theme_set(gg.theme + theme(axis.text.x = element_text(angle = 90)))
  re.plot <- ggplot(df.re, aes(x = id, y = value, color = id))+
    geom_errorbar(aes(ymin = low, ymax = up), width = 0.75)+
    geom_point(size = 2) +
    xlab(TeX("$ i")) +
    ylab(TeX("$ \\widehat{u}_i")) +
    theme(legend.position = "none")
  
  print(re.plot)

  # Save plots
  if(save){
    ggsave(filename = "functional_mean.pdf", path = path.plot, plot = int.plot, width = width.pdf, height = height.pdf, units = unit, dpi = res)
    ggsave(filename = "f_ca_biv.pdf",path = path.plot, plot = ca.plot, device = "pdf", width = width.pdf, height = height.pdf, units = unit, dpi = res)
    ggsave(filename = "f12.pdf", path = path.plot, plot = p3, device = "pdf", width = width.biv, height = height.biv, units = unit, dpi = res)
    ggsave(filename = "f1_CI.pdf", path = path.plot, plot = p.long.ci, device = "pdf", width = width.biv, height = height.biv, units = unit, dpi = res)
    ggsave(filename = "f2_CI.pdf", path = path.plot, plot = p.lat.ci, device = "pdf", width = width.biv, height = height.biv, units = unit, dpi = res)
    ggsave(filename = "f_re.pdf",path = path.plot, plot = re.plot, device = "pdf", width = width.pdf, height = height.pdf, units = unit, dpi = res)
    ggsave(filename = "FPC_RE.pdf", plot = plot.fpc, path = path.plot, width = 11.27, height = 8.27, unit = "in", dpi = res)
  }
}

pffr.paper.plot(coef, save = F, res = 72)