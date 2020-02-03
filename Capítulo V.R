library(stringr)
library(spdep)
library(rgdal)
library(magrittr)
library(ggplot2)
library(sf)



#======================================================

# Load data
map <- read.table("C:/Users/franc/Dropbox/Franca/Doctorado/tmedia_vs_Kd/Muestreo.txt", sep = "\t", header = T)
colnames(map)
map <- st_as_sf(map, coords = c("Xt", "Yt"), crs = 32720)

# map <- readOGR(dsn="R:/Dropbox/Dout/Data Dout", layer="test_map")
# head(map@data)

# Variables to use in the correlation: white and black population in each census track
x <- map$Kda
y <- map$tmedia




#======================================================
# Programming some functions

# Bivariate Moran's I
moran_I <- function(x, y = NULL, W){
  if(is.null(y)) y = x
  
  xp <- scale(x)[, 1]
  yp <- scale(y)[, 1]
  W[which(is.na(W))] <- 0
  n <- nrow(W)
  
  global <- (xp%*%W%*%yp)/(n - 1)
  local  <- (xp*W%*%yp)
  
  list(global = global, local  = as.numeric(local))
}


# Permutations for the Bivariate Moran's I
simula_moran <- function(x, y = NULL, W, nsims = 2000){
  
  if(is.null(y)) y = x
  
  n   = nrow(W)
  IDs = 1:n
  
  xp <- scale(x)[, 1]
  W[which(is.na(W))] <- 0
  
  global_sims = NULL
  local_sims  = matrix(NA, nrow = n, ncol=nsims)
  
  ID_sample = sample(IDs, size = n*nsims, replace = T)
  
  y_s = y[ID_sample]
  y_s = matrix(y_s, nrow = n, ncol = nsims)
  y_s <- (y_s - apply(y_s, 1, mean))/apply(y_s, 1, sd)
  
  global_sims  <- as.numeric( (xp%*%W%*%y_s)/(n - 1) )
  local_sims  <- (xp*W%*%y_s)
  
  list(global_sims = global_sims,
       local_sims  = local_sims)
}

#library("UScensus2000tract")
#data("oregon.tract")

#======================================================
# Adjacency Matrix (Queen)
misPoligonos <- st_make_grid(map)
plot(misPoligonos)
#as(spgrdWithin, "SpatialPolygons")


#nb <- nb2listw(dnearneigh(map, 0, all.linked), style = "B", zero.policy = TRUE)  
#Wij <- as.matrix( as(nb, "symmetricMatrix") ) 	






nb1 <- nbdists(knn2nb(knearneigh(st_coordinates(map), k = 25)),
               st_coordinates(map))
all.linked <-100000 #max(unlist(nbdists(knn2nb(knearneigh(st_coordinates(map))), 
# st_coordinates(map))))  
nb1 <- dnearneigh(st_coordinates(map), 0, all.linked)
lw <- nb2listw(nb1, style = "B", zero.policy = T)
W  <- as(lw, "symmetricMatrix")
W  <- as.matrix(W/rowSums(W))
W[which(is.na(W))] <- 0


#======================================================
# Calculating the index and its simulated distribution
# for global and local values

m <- moran_I(x, y, W)

# Global Moral
global_moran <- m[[1]][1]
#> 0.2218409

# Local values
m_i <- m[[2]] 

# local simulations
sims <- simula_moran(x, y, W, nsims = 100)

local_sims <- sims$local_sims


# global pseudo p-value  
# get all simulated global moran
global_sims <- sims$global_sims

# Proportion of simulated global values taht are higher (in absolute terms) than the actual index 
moran_pvalue <- sum(abs(global_sims) > abs( global_moran )) / length(global_sims)
#> 0


# Identifying the significant values 
alpha <- .05  # for a 95% confidence interval
probs <- c(alpha/2, 1-alpha/2)
intervals <- t( apply(local_sims, 1, function(x) quantile(x, probs=probs)))
sig       <- ( m_i < intervals[,1] )  | ( m_i > intervals[,2] )
sign <- as.data.frame(sig)


#======================================================
# Preparing for plotting


# Convert shape file into sf object
map_sf     <- st_as_sf(map)
map_sf$sig <- sig


# Identifying the LISA clusters
xp <- scale(x)[,1]
yp <- scale(y)[,1]


patterns <- as.character( interaction(xp > 0, W%*%yp > 0) )
patterns <- patterns %>% 
  str_replace_all("TRUE","High") %>% 
  str_replace_all("FALSE","Low")

patterns[map_sf$sig==0] <- "No significant"
map_sf$patterns <- patterns


# Rename LISA clusters
map_sf$patterns2 <- factor(map_sf$patterns, levels=c("High.High", "High.Low", "Low.High", "Low.Low", "Not significant"),
                           labels=c("High income - High access gain", "High income - Low access gain", "Low income - High access gain","Low income - Low access gain", "Not significant"))





### PLOT

ggplot() +
  geom_sf(data=map_sf, aes(fill=patterns, color=patterns)) +
  scale_color_manual(values = c("red", "pink", "light blue", "dark blue", "grey80")) +
  scale_fill_manual(values = c("red", "pink", "light blue", "dark blue", "grey80")) +
  guides(fill = guide_legend(title="LISA clusters"),
         color = guide_legend(title="LISA clusters"))# +
# theme_minimal()
