#' In the following code you find the all the GFAMM routine, as well as necessary steps to perform the analysis presented in the paper
#' You can find a runned version of the script by loading sediment_pffr.RData (to save time)
#'  
# Preamble ----
library(latex2exp)
library(reshape2)
library(tidyverse)
rm(list=ls())
#load data ----
library(readr)

# Load the data and transform 0 1 into a probability. This is done by dynamic smoothing.
source("./user_functions.R") #some useful functions
data <- read.table("./data/paleoflood_count_data.txt", header=TRUE, quote="\"")
data$Year <- c(sort(0:2000,T),-sort(1:7050,F))
data <- data[order(data[,1],decreasing = F),] #order by increasing time and remove years

lake.name <- c("Blanc Belledonne",
               "Inferiore di Laures",
               "Bourget",
               "Anterne",
               "Blanc Aig. Rouges",
               "Foréant",
               "Allos",
               "Savine",
               "Oeschinen",
               "Mondsee",
               "Ammersee",
               "Iffig",
               "Glattalp",
               "Alzasca",
               "Cadagno",
               "Thun",
               "Garlate",
               "Grimsel",
               "Baldegg",
               "Faelen",
               "Hinterburg",
               "Hinterer Schw.",
               "Lauerz",
               "Seelsberg",
               "Trueb",
               "Ledro",
               "di Braies")
lake.name.abbr <- paste(lake.name," (", names(data[,-1]),")", sep ="")
lakes <- names(data[,-1])

# visualise missing values
library(naniar)
df.miss <- tibble( lake = factor(lake.name.abbr, levels = lake.name.abbr),
                   miss = round(colSums(is.na(data[,-1]))/nrow(data),3)*100)

g.miss <- ggplot(df.miss, aes(x = miss, y = lake)) + geom_col() + theme_bw()+ 
  theme(text = element_text(size=11),
        legend.title.align = 0.5,
        axis.text = element_text(size=rel(1)),
        #legend.justification = "left",
        legend.text = element_text(size = 10)
        #legend.key.height= unit(2.5, units = "lines"),
  ) + xlab("% of missing observations") + ylab("Lake:")
  
colnames(data) <-c("Year", lake.name.abbr)
rownames(data) <- data$Year

miss <- as.data.frame(as.matrix(data[,-1]))
names(miss) <- lake.name.abbr

g.miss <- 
  vis_miss(miss, show_perc_col = F) +
  theme_bw()+ 
  theme(text = element_text(size=12),
        legend.position = "none",
        legend.title.align = 0.5,
        axis.text = element_text(size=rel(1)),
        legend.text = element_text(size = 10))+
  scale_y_continuous(breaks = c(0,2500,5000,7550, nrow(data)), labels = c("7050 BC", "4550 BC", "2050 BC", "1500 AD", "2000 AD")) + 
  coord_flip() +
  scale_x_discrete(position = "bottom") +
  ylab("") 

g.miss

#ggsave("missing_values.pdf", path = "./plot", device = "pdf", dpi = 72, plot = g.miss, units = "in", width = 11, height = 10)

# Load the covariate information
X.base <- read_csv("data/flood_covariates.csv", col_types = c("f","d","d","d","d","d","f","f", "f")) %>%
  mutate(full_name = as.factor(lake.name))

X <- X.base %>% 
  mutate(ca = log(ca),
         caalt = log(caalt),
         alt = log(alt))
X

# temporal selection of data, how many years do we consider in the analysis? 
to <- 1000 # select last "to" years of observations

data1 <- data %>%
  filter(Year>=(2000-to+1)) %>%
  select(-Year)

year.s <- data$Year[data$Year>=(2000-to+1)] #vector of years of observations

# Mirror last years of data
mirror <- 16
data1 <- rbind(data1, tail(data1, mirror)[1:(mirror-1),])
year.s <- c(year.s, max(year.s) + 1:(mirror-1))

rownames(data1) <- year.s

# PRELIMINARY: Two sided moving average non-parametric smoothing ----
library(mgcv)
h <- 15 # half-bandwith size for smoothing
data.bin <- bin.smooth(data1, h = h, na.replace = F, min.h = floor(h)); colnames(data.bin) <- lakes

# smooth the discrete data, keep the logit scale ----
l.s <- list()
data.s <- data1
for(i in 1:ncol(data.bin)){
  
  k.gam <- 60 # number of basis functions for the GAM
  obs <- !is.na(data.bin[,i])
  id <- lake.name[i]
  
  if(id %in% c("GRI","BRA")){k.gam <- 5} # we use a lower number of basis functions since these lakes have only a few years of observations
  
  l.s[[i]]<-gam(data.bin[obs,i]~ s(year.s[obs], bs = "cr", m = c(3,2), k = k.gam),
                family = betar(link = "logit",eps = 0.01),
                control = gam.control(epsilon = 1e-4,maxit = 10000))
  
  data.s[obs,i] <- predict(l.s[[i]]) # prediction of additive score
  
  print(paste("Lake",i))
};gc()

rownames(data.s) <- year.s
colnames(data.s) <- colnames(data.bin) <- lakes

# graphical comparison of GAM prediction and non-parametric estimators
df.s <- melt(as.matrix(data.s)) %>%
  rename(year = Var1,
         id = Var2 ) %>%
  filter(!is.na(value)) %>%
  mutate(fit = 1/(1+exp(-value))) %>%
  select(-value)

df.bin <- melt(as.matrix(data.bin)) %>%
  rename(year = Var1,
         id = Var2,
         bin = value) %>%
  filter(!is.na(bin))

df.smooth <- bind_cols(df.s, df.bin %>% select(-year, -id)) 

theme_set(theme_bw())
smooth.plot <- ggplot(df.smooth, aes(x = year, color = id)) +
  geom_line(aes(y = fit), linewidth = 0.75) +
  geom_point(aes(y = bin), size = 0.25) +
  theme(legend.position = "right" ) +
  guides(color=guide_legend("Lake:")) +
  theme(text = element_text(size = 11),
        axis.text = element_text(size = rel(1)),
        legend.text = element_text(size = 10)) +
  xlab(TeX("t")) + 
  xlim(2000-to,2000) +
  ylab(TeX("$ \\widehat{\\gamma}_{i}(t) \\; and \\; W_{ij}"))

smooth.plot
#ggsave("smooth_sample.pdf", path = "./plot", plot = smooth.plot, units = "in", width = 11, height = 5)

# Functional Regression Model ----
# covariates for each lake
library(refund)
data.cov <- data2 <- list( data = X %>%
                             mutate(g = as.factor(g),
                                    id = as.factor(id),
                                    zone = as.factor(zone), 
                                    season = as.factor(season)))
rownames(data.cov$data) <- rownames(data2$data) <- 1:ncol(data.s)

ydata <- melt(as.matrix(data.s))%>%
  mutate(.obs = as.numeric(Var2))%>%
  select(.obs, .index = Var1, .value = value)%>%
  filter(!is.na(.value))

all(ydata$.obs %in% rownames(data2$data)) # check if all lakes id's corresponds

# Model specifications ----
basis = "cr" # type of basis function

yindex <- unique(sort(ydata$.index))

knot.list <- list(ca = quantile(data2$data$ca, probs = seq(0, 1, length.out = 10)))

#' regression formula (Model (21))
f <- Y ~ c(s(ca, bs = basis)) +
  s(long, bs = basis) +
  s(lat, bs = basis) +
  c(s(id, bs ="re"))

# functional regression using gam smoothed curves
pffr <- pffr(formula = f,
             yind = yindex,
             data = data2$data,
             ydata = ydata,
             method = "REML",
             knots = knot.list,
             bs.int = list(bs = "cr", k = 20, m = c(2, 1)),
             bs.yindex = list(bs = "cr", k = 5, m = c(2, 1))
);gc()

# assess the quality of the regression model
library(mgcViz)
summary(pffr)

pffr1 <- pffr
class(pffr1)[1] <- "gam"
b <- mgcViz::getViz(pffr1,nsim = 500)
check(b, pch = 16, cex = .8, type = "response")

# In sample predictions
pffr$pffr$nyindex <- length(yindex)
pffr$pffr$yind <- as.numeric(year.s)

fit <- predict(pffr,type="response")%>%
  mutate(.value1 = 1/(1+exp(-.value)),
         .bin = na.omit(c(data.bin)),
         .smooth = na.omit(c(as.matrix(data.s))),
         .smooth1 = 1/(1+exp(-.smooth)),
         .res = .value1-.bin,
         .res1 = .value1 - 1/(1+exp(-.smooth)),
         id = lakes[.obs])

# FPCA of residual curves and fitted predictions
fpca.res <- fpca.sc( ydata = data.frame(.id = fit$.obs,
                                        .index = fit$.index,
                                        .value = fit$.res1),
                     nbasis = 10,
                     pve = 0.95);gc()

fpca.fit <- fpca.sc( ydata = data.frame(.id = fit$.obs,
                                        .index = fit$.index,
                                        .value = fit$.value1),
                     nbasis = 10,
                     pve = 0.95);gc()
#' you should run the following statement the first time you run main.R
#' this way you don't need to run the functional regression each time you want to use other scripts
#save.image("sediment_pffr.RData")
