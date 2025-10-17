library(sp)
library(sf)
library(rnaturalearth)
library(fmesher)
library(INLA)
library(INLAspacetime)
library(inlabru)
library(ggplot2)
library(fields)

setwd(here::here())

### CRS PROJ4 string to project the Saudi Arabia map
### https://epsg.io/20440
crsKM <- st_crs("+proj=utm +zone=40 +ellps=intl +towgs84=-143,-236,7,0,0,0,0 +units=km +no_defs +type=crs")

## get the world map
worldMapll <- ne_countries()

## project the world map
worldMap <- st_transform(
    worldMapll, crsKM)

## read the data for the analysis
tavg_day <- read.csv("tavg_day.csv")

## convert the dataset into spatial dataset and project it
dTavg <- st_transform(
    st_as_sf(
        x = tavg_day, 
        coords = c(3, 2),
        crs = st_crs(worldMapll)),
    crsKM
)

## data colums refer to temperature for weeks 21 to 34 (14 weeks)
head(dTavg, 2)
jjd <- 5:ncol(tavg_day)

Dates <- as.Date(gsub("X", "", colnames(tavg_day)[5:ncol(tavg_day)]), "%Y%m%d")

plot(Dates, t(tavg_day[, jjd])[, 1], type = "n",
     ylim = range(tavg_day[, jjd], na.rm = TRUE))
for(k in 1:nrow(tavg_day))
    lines(Dates, t(tavg_day[, jjd])[, k], col=k)

par(mfrow = c(1, 1), mar = c(0,0,0,0))
plot(st_geometry(dTavg))
stlines(t(tavg_day[,jjd]), as(dTavg, "Spatial"))
plot(worldMap, add = TRUE, col = 'transparent')

(nt <- length(jjd))
(ns <- nrow(tavg_day))

## prepare the data in the long format
colnames(tavg_day) ## data is in columns 5 to 18
longdf <- data.frame(
    xloc = rep(st_coordinates(dTavg)[, 1], nt),
    yloc = rep(st_coordinates(dTavg)[, 2], nt),
    time = rep(1:nt, each = ns),
    elevation = rep(tavg_day$elevation, nt),
    temp = unlist(tavg_day[, jjd])
)

## define the data model (with inlabru)
pprec <- list(prec = list(prior = "pc.prec", param = c(1, 0.05)))
data_model <- bru_obs(
  formula = temp ~ ., 
  family = "gaussian",
  control.family = list(hyper = pprec), 
  data = longdf)

## unite polygons to create a boundary
bnd <- worldMap[worldMap$name == "Saudi Arabia", ]
bb <- st_bbox(bnd)

## define the mesh around the Ethiopia boundary
mesh <- fm_mesh_2d(
    boundary = bnd,
    max.edge = c(100, 300),
    offset = c(100, 400), 
    cutoff = 50
)
mesh$crs <- st_crs(bnd)

mesh$n

ggplot() + theme_minimal() +
    geom_sf(data = worldMap, fill = 'transparent') +
    xlim(bb[c(1,3)] + c(-1, 1)*500) + 
    ylim(bb[c(2,4)] + c(-1, 1)*500) +
    gg(mesh)

plot(mesh)
plot(st_geometry(worldMap), add = TRUE)
stlines(t(tavg_day[,jjd]), as(dTavg, "Spatial"))


## to fit a spacetime model, we need also a mesh over time
nt
tmesh <- fm_mesh_1d(loc = seq(1, nt, 3))
tmesh$loc
tmesh$n

## define the spacetime SPDE model
stModel <- stModel.define(
    smesh = mesh,
    tmesh = tmesh,
    model = "102",
    control.priors = list(
        prs = c(30, 0.05), ## P(spatial range < 30) = 0.05
        prt = c(1, 0.05),  ## P(temporal range < 1) = 0.05
        psigma = c(1, 0.05)   ## P(sigma > 1) = 0.05
    )
)

## model components
mcomps <- ~ 1 + elevation +
    spacetime(list(space = cbind(xloc, yloc), 
                   time = time),
              model = stModel)


## fit the model using inlabru
stfit <- bru(mcomps, data_model,
             options = list(verbose = TRUE))

stfit$cpu.used

## summary of the fixed effects: intercept and elevation/1000
stfit$summary.fixed

## summary of the error noise and the SPDE model parameters
stfit$summary.hyperpar

exp(stfit$summary.hyperpar$mean[2:4])

## plot posterior mean of the fitted values against observed
plot(stfit$summary.fitted.values$mean[1:(nt*ns)], longdf$temp)
cor(stfit$summary.fitted.values$mean[1:(nt*ns)], longdf$temp, use = 'pair')^2

## create a grid/projector/mapper to visualization of the
## spatial effect, and later the fitted values, in a map
rr <- apply(matrix(bb, 2), 1, diff)
rr

projGrid <- inla.mesh.projector(
    mesh = mesh,
    xlim = bb[c(1,3)], 
    ylim = bb[c(2,4)],
    dims = round(100 * rr/rr[1])
)

str(projGrid)

## conver the grid locations into an sf object
sfGrid <- st_as_sf(
    as.data.frame(projGrid$lattice$loc),
    coords = 1:2,
    crs = crsKM
)

## find which grid points are outside the Ethiopia boundary
igrid.out <- which(sapply(st_within(sfGrid, bnd), length)==0)

## organize the spacetime effect at the mesh nodes in a matrix
st.mean <- matrix(stfit$summary.random$spacetime$mean, ncol = tmesh$n)

## project and organize the spacetime posterior mean into an array
## which represent the grid at each time (an 3D array)
st.array <- array(
    NA, 
    c(length(projGrid$x), length(projGrid$y), tmesh$n)
)
for(k in 1:tmesh$n) {
    st.array[, , k] <- inla.mesh.project(projGrid, st.mean[, k])
    st.array[,,k][igrid.out] <- NA ## values outside boundary st to NA
}

## visualize the spacetime effect
par(mfrow = c(4,5), mar = c(0,0,0,0))
for(k in 1:tmesh$n) {
    image.plot(
        x = projGrid$x,
        y = projGrid$y,
        st.array[, , k],
        asp = 1,
        axes = FALSE
    )
    plot(st_geometry(worldMap), add = TRUE)
}


## download the TOPO dataset: elevation (on land) or depth (on ocean)
if(!file.exists("ETOPO2.RData")) {
    download.file(
        url = paste0(
            "http://leesj.sites.oasis.unc.edu/",
            "FETCH/GRAB/RPACKAGES/",
            "ETOPO2.RData"
        ),
        destfile = "ETOPO2.RData"
    )
}

### load the TOPO dataset
load("ETOPO2.RData")

## extract the longitude (and fix it), and latitude
elon <- attr(ETOPO2, "lon")
elon[elon >= 180] <- 180 - rev(elon[elon >= 180])
elat <- attr(ETOPO2, "lat")

### fix the order of the columns (swap)
ETOPO2 <- ETOPO2[, ncol(ETOPO2):1]

## transform the grid coordinates into latlong
sfGrid.ll <- st_transform(sfGrid, st_crs(worldMapll))
## extract it as a matrix
glocs <- st_coordinates(sfGrid.ll)

## find which pixels in the long/lat of TOPO data belong each grid location
ij <- list(
  i = findInterval(glocs[, 1], c(-180, elon + 1 / 60)),
  j = findInterval(glocs[, 2], elat)
)

## extract the TOPO data values at the grid locations
etopo.grid <- sapply(1:nrow(glocs), function(i) ETOPO2[ij$i[i], ij$j[i]])
summary(etopo.grid)

### make NA the grid locations outside the boundary
etopo.grid[igrid.out] <- NA
summary(etopo.grid)

## visualize
ggplot() +
    geom_sf(aes(color = etopo.grid), data = sfGrid, cex = 2) +
    scale_color_distiller(palette = "RdBu", direction = 1, na.value = 'transparent')

## in order to compute the uncertainty for the predictions
## make a prediction scenario and re-evaluate the model

## make a prediction dataset
head(longdf)
preddf <- data.frame(
    xloc = rep(projGrid$lattice$loc[, 1], tmesh$n),
    yloc = rep(projGrid$lattice$loc[, 2], tmesh$n),
    time = rep(tmesh$loc, nrow(projGrid$lattice$loc)),
    elevation = rep(etopo.grid, tmesh$n),
    temp = NA ## set temp to NA so to predict it
)

## re-define the data model
data_model_pred <- bru_obs(
  formula = temp ~ ., 
  family = "gaussian",
  control.family = list(hyper = pprec), 
  data = rbind(longdf, preddf)
)


## consider the fitted mode
stfit$mode$theta

## re-evaluate the model at the fitted mode
stfit.pred <- bru(
    mcomps,
    data_model_pred,
    options = list(
        verbose = TRUE,
        control.mode = list(
            theta = stfit$mode$theta,
            fixed = TRUE
        )
    )
)

stfit.pred$cpu.used

## create an index for the prediction
i.pred <- nrow(longdf) + 1:nrow(preddf)
length(i.pred)
length(projGrid$x) * length(projGrid$y) * tmesh$n

## organize the predicted values (posterior mean) into an array
## which represent the grid at each time (an 3D array)
pred.array <- array(
  stfit.pred$summary.fitted.values$mean[i.pred], 
  c(length(projGrid$x), length(projGrid$y), tmesh$n)
)
for(k in 1:tmesh$n)
    pred.array[,,k][igrid.out] <- NA ## values outside boundary st to NA

## organize the standard error of the predictions to visualize
sd.pred.array <- array(
  stfit.pred$summary.fitted.values$sd[i.pred], 
  c(length(projGrid$x), length(projGrid$y), tmesh$n)
)
for(k in 1:tmesh$n)
    sd.pred.array[,,k][igrid.out] <- NA ## values outside boundary st to NA

### visualize the predicted (posterior mean) temperature
par(mfrow = c(4,5), mar = c(0,0,0,0))
for(k in 1:tmesh$n) {
    image.plot(
        x = projGrid$x,
        y = projGrid$y,
        pred.array[, , k],
        asp = 1,
        axes = FALSE
    )
    plot(st_geometry(bnd), add = TRUE)
    plot(st_geometry(worldMap), add = TRUE)
}

## visualize the standard error of the prediction (for one time, other times are similar)
image.plot(
    x = projGrid$x,
    y = projGrid$y,
    sd.pred.array[, , 1],
    asp = 1,
    axes = FALSE
)
plot(st_geometry(bnd), add = TRUE)
plot(st_geometry(worldMap), add = TRUE)
points(st_geometry(dTavg), pch = 8)

## NOTICE that the error (std dev of the prediction)
## near the stations is smaller 
