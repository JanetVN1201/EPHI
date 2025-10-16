library(sf)
library(rnaturalearth)
library(fmesher)
library(INLA)
library(ggplot2)

### https://epsg.io/20138
crsEthiopia <- st_crs(
    "+proj=utm +zone=38 +a=6378249.145 +rf=293.465 +towgs84=-165,-11,206,0,0,0,0 +units=km +no_defs +type=crs")

worldMap <- st_transform(
    ne_countries(scale = "medium", returnclass = "sf"),
    crsEthiopia)

map2 <- st_read("map2EthiopiaData.shp")
map2km <- st_transform(map2, crsEthiopia)

tavg <- read.csv("tavg.csv")

head(tavg)

summary(lm(X29 ~ elevation, data = tavg))
summary(lm(X29 ~ I(elevation/1000), data = tavg))

inla(X29 ~ I(elevation/1000), data = tavg)$summary.fixed

Tavg <- st_transform(
    st_as_sf(
        x = tavg, 
        coords = c(3, 2),
        crs = st_crs(map2)),
    crsEthiopia
)

bb <- st_bbox(map2km)
bb

ggplot() + theme_minimal() +
    geom_sf(data = worldMap,
            fill = rgb(165/256,142/256,42/256))+
    geom_sf(data = map2) +
    geom_sf(aes(color = X29,
                size = elevation), Tavg) +
    xlim(c(-1, 1)*1200) +
    ylim(c(0, 2000)) + 
    scale_color_distiller(
        palette = "RdBu"
    )

bnd <- st_union(map2km)

mesh <- fm_mesh_2d(
    ##    loc.domain = cbind(mean(bb[c(1,3)]), mean(bb[c(2,4)])),
    boundary = bnd,
    max.edge = c(70, 300),
    offset = c(700, 2500), 
    cutoff = 30
)
mesh$crs <- st_crs(bnd)

mesh$n

plot(mesh)

ggplot() + theme_minimal() +
    geom_sf(data = worldMap,
            fill = rgb(165/256,142/256,42/256))+
    geom_sf(data = map2km) +
    inlabru::gg(mesh) +
    geom_sf(aes(color = X29,
                size = elevation), Tavg) +
    scale_color_distiller(
        palette = "RdBu"
    ) +
    xlim(c(-3500, 3500)) + 
    ylim(c(-2200, 4300)) 
    

spdeModel <- inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = c(30, 0.05), ## P(range < 30) = 0.05
    prior.sigma = c(1, 0.05)   ## P(sigma > 1) = 0.05
)

Amap <- inla.spde.make.A(
    mesh = mesh,
    loc = st_coordinates(Tavg)
)

smodel <- X29 ~ I(elevation/1000) +
    f(spatial, model = spdeModel, A.local = Amap)

dataf <- data.frame(
    tavg[c("X29", "elevation")], spatial = NA
)

sfit <- inla(
    formula = smodel,
    data = dataf
)

sfit$summary.fixed

sfit$summary.hyperpar

plot(sfit$summary.fitted.values$mean, dataf$X29)
cor(sfit$summary.fitted.values$mean, dataf$X29, use = 'pair')^2

bb

projGrid <- inla.mesh.projector(
    mesh = mesh,
    xlim = c(-850, 850),
    ylim = c(350, 1700),
    dims = c(170, 205)
)

str(projGrid)

sfGrid <- st_as_sf(
    as.data.frame(projGrid$lattice$loc),
    coords = 1:2,
    crs = crsEthiopia
)
igrid.out <- which(sapply(st_within(sfGrid, bnd), length)==0)

gridSmean <- inla.mesh.project(
    projGrid,
    sfit$summary.random$spatial$mean)
gridSmean[igrid.out] <- NA

library(fields)

image.plot(
    x = projGrid$x,
    y = projGrid$y,
    gridSmean,
    asp = 1
)
plot(st_geometry(map2km), add = TRUE)
plot(st_geometry(worldMap), add = TRUE)


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

### extract the longitude (and fix) and latitude
load("ETOPO2.RData")

elon <- attr(ETOPO2, "lon")
elon[elon >= 180] <- 180 - rev(elon[elon >= 180])
elat <- attr(ETOPO2, "lat")


### fix the order of the lines
ETOPO2 <- ETOPO2[, ncol(ETOPO2):1]

sfGrid.ll <- st_transform(sfGrid, st_crs(map2))
glocs <- st_coordinates(sfGrid.ll)

ij <- list(
  i = findInterval(glocs[, 1], c(-180, elon + 1 / 60)),
  j = findInterval(glocs[, 2], elat)
)

etopoll <- sapply(1:nrow(glocs), function(i) ETOPO2[ij$i[i], ij$j[i]])
etopoll[igrid.out] <- NA

summary(etopoll)

ggplot() +
    geom_sf(aes(color = etopoll), data = sfGrid, cex = 2) +
    scale_color_distiller(palette = "RdBu", direction = 1)

gridfitt <- sfit$summary.fixed$mean[1] +
    (etopoll/1000) * sfit$summary.fixed$mean[2] +
    gridSmean

str(gridfitt)

par(mfrow = c(1,1), mar = c(0,0,0,0))
image.plot(
    x = projGrid$x,
    y = projGrid$y,
    gridfitt,
    asp = 1
)
plot(st_geometry(map2km), add = TRUE)
plot(st_geometry(worldMap), add = TRUE)
