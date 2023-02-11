######## Example code for point pattern with spde #########
######## A case study of Nigeria ##########################
library(viridis)
library(rgdal)
library(readxl)
library(tidyverse)
library(rgeos)
library(INLA)
set.seed(12345)
# Simulate violent event location
longitude  <- sample(seq( 4.710416,14.65429,length=20000),2000,replace = T)
latitude   <- sample(seq( 4.286446,14.37073,length=20000),2000,replace = T)
dat   <- data.frame(longitude,latitude)
n     <- nrow(dat) # get total number of event
xy    <- cbind(dat$longitude,dat$latitude) # extract long & lat

# Roughly determine nigeria boundary

loc.2D <- cbind(c(min(dat$longitude),
                 max(dat$longitude),
                 max(dat$longitude),
                 min(dat$longitude),
                 min(dat$longitude)),
               ##############
                c(min( dat$latitude),
                 min(dat$latitude),
                 max(dat$latitude),
                 max(dat$latitude),
                 min(dat$latitude))
)
##### Determine structure of mesh 
mesh <- inla.mesh.2d(loc.domain = loc.2D, offset = c(1, 1), 
                     max.edge = c(1, 1), cutoff =0.7)
nfield <- mesh$n

spde <- inla.spde2.pcmatern(mesh = mesh,
                            # PC-prior for spatial range and vriance
                            prior.range = c(0.05, 0.01),
                            prior.sigma = c(1, 0.01)) 

####  Dual mesh function, sourced from R-INLA (Lindgren, F., Rue, H., and Lindstrom, J. (2011))
####  included here for self containment

book.mesh.dual <- function(mesh) {
  if (mesh$manifold=='R2') {
    ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
      colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
    library(parallel)
    pls <- mclapply(1:mesh$n, function(i) {
      p <- unique(Reduce('rbind', lapply(1:3, function(k) {
        j <- which(mesh$graph$tv[,k]==i)
        if (length(j)>0) 
          return(rbind(ce[j, , drop=FALSE],
                       cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1], 
                             mesh$loc[mesh$graph$tv[j, k], 2] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
        else return(ce[j, , drop=FALSE])
      })))
      j1 <- which(mesh$segm$bnd$idx[,1]==i)
      j2 <- which(mesh$segm$bnd$idx[,2]==i)
      if ((length(j1)>0) | (length(j2)>0)) {
        p <- unique(rbind(mesh$loc[i, 1:2], p,
                          mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2, 
                          mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
        yy <- p[,2]-mean(p[,2])/2-mesh$loc[i, 2]/2
        xx <- p[,1]-mean(p[,1])/2-mesh$loc[i, 1]/2
      }
      else {
        yy <- p[,2]-mesh$loc[i, 2]
        xx <- p[,1]-mesh$loc[i, 1]
      }
      Polygon(p[order(atan2(yy,xx)), ])
    })
    return(SpatialPolygons(lapply(1:mesh$n, function(i)
      Polygons(list(pls[[i]]), i))))
  }
  else stop("It only works for R2!")
}


dmesh        <- book.mesh.dual(mesh=mesh)
window.polys <- Polygons(list(Polygon(loc.2D)), '0')
windowSP     <- SpatialPolygons(list(window.polys))

weight <- sapply(1:length(dmesh), function(i) {
  if (gIntersects(dmesh[i, ], windowSP))
    return(gArea(gIntersection(dmesh[i, ], windowSP)))
  else return(0)
})

response <- rep(0:1, c(nfield, n))
offset.response <- c(weight, rep(0, n)) 
unobserved.A <- Diagonal(nfield, rep(1, nfield))
observed.A <- inla.spde.make.A(mesh, xy)
A <- rbind(unobserved.A, observed.A)
datastack <- inla.stack(
  data = list(y = response, e = offset.response), 
  A = list(1, A),
  effects = list(list(intercept = rep(1, nfield + n)), list(i = 1:nfield)),
  tag = 'pp')
# Estimation
Model <- inla(y ~ 0 + intercept + f(i, model = spde), 
                        family = 'poisson', data = inla.stack.data(datastack), 
                        control.predictor = list(A = inla.stack.A(datastack)), 
                        E = inla.stack.data(datastack)$e,
              control.compute = list(config = TRUE,dic = TRUE,waic=TRUE,cpo=TRUE))

# Read in boundary file
boundary <- read.csv("boundary.csv") # Available in this repository.

# Projection into unsample regions
stepsize <- 0.1 # the smaller the value the higher the smoothness
nxy <- round(
  c(diff(range( boundary$X1)), 
    diff(range( boundary$X2))) / stepsize)

projgrid <- inla.mesh.projector(
  mesh, xlim = range(boundary$X1), 
  ylim = range(boundary$X2), dims = nxy)

xmean <- inla.mesh.project(
  projgrid, Model$summary.random$i$mean)

df <-  expand.grid(x = projgrid$x, y = projgrid$y)
df$mean <- as.vector(xmean)
# delete point outside Nigeria map
ind <- point.in.polygon(
  df$x, df$y,
  boundary$X1, boundary$X2
)
dff <- df[which(ind == 1), ]

# Plot
ggplot(dff, aes(x=x, y=y)) +
  geom_tile(aes(fill = mean)) +
  scale_fill_viridis_c(option = "H", direction = 1,name="") +
  labs(title ="",
       y = "",x="") +
  theme_light()

