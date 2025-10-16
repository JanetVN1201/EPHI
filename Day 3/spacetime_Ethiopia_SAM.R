
library(sf)
library(ggplot2)
library(ggpubr)
library(INLA)

getwd() ## the current working directory
setwd("put the folder where map2EthiopiaData is") ## change the working directory

## read the shapefile (map + data)
map2 <- st_read("map2EthiopiaData.shp")

## easy map plot
plot(map2[,"SAM"])

## ggplot setup to map
gg0 <- ggplot() + theme_minimal()
gg2 <- gg0 + scale_fill_distiller(
  palette = "RdBu", 
  transf = "log", 
  breaks = sort(c(3 * 10^(-4:1), 1 * 10^(-4:2))),
  labels = function(x) format(x, digits = 3, scientific = FALSE)
  )

## visualize SAM/Pop0to4 map (sum from January to October)
gg2 +  geom_sf(aes(fill = SAM / Pop0to4), map2)

## visualize SAM map at each month
ggarrange(
    gg2 +  geom_sf(aes(fill = SAM01 / Pop0to4), map2),
    gg2 +  geom_sf(aes(fill = SAM02 / Pop0to4), map2),
    gg2 +  geom_sf(aes(fill = SAM03 / Pop0to4), map2),
    gg2 +  geom_sf(aes(fill = SAM04 / Pop0to4), map2),
    gg2 +  geom_sf(aes(fill = SAM05 / Pop0to4), map2),
    gg2 +  geom_sf(aes(fill = SAM06 / Pop0to4), map2),
    gg2 +  geom_sf(aes(fill = SAM07 / Pop0to4), map2),
    gg2 +  geom_sf(aes(fill = SAM08 / Pop0to4), map2),
    gg2 +  geom_sf(aes(fill = SAM09 / Pop0to4), map2),
    gg2 +  geom_sf(aes(fill = SAM10 / Pop0to4), map2),
    gg2 +  geom_sf(aes(fill = SAM / Pop0to4), map2)
)

summary(ppHealthFac <- map2$Popul / map2$nHlthFc)

gg2 + geom_sf(aes(fill = ppHealthFac), data = map2)

rate <- sum(map2$SAM) / sum(map2$Pop0to4)
rate

map2$ESAM <- rate * map2$Pop0to4
sum(map2$SAM) == sum(map2$ESAM)

plot(map2$ESAM, map2$SAM, pch = 8, bty = 'n')

plot(map2$ESAM, map2$SAM + (map2$SAM==0)*0.5, pch = 8, bty = 'n', log = 'xy')

plot(ppHealthFac,
     map2$SAM / map2$ESAM,
     log = 'x', pch = 8, bty = 'n',
     xlab = 'Population / Health Facilities',
     ylab = 'Observed standardized rate')

plot(ppHealthFac,
     ifelse(map2$SAM == 0, 0.5, map2$SAM) / map2$ESAM,
     log = 'xy', pch = 8, bty = 'n',
     xlab = 'Population / Health Facilities',
     ylab = 'Observed standardized rate')

n <- nrow(map2)
nt <- 10

head(st_drop_geometry(map2)[, 18:27])

longdf <- data.frame(
    i = rep(1:n, nt),
    t = rep(1:nt, each = n),
    hf = rep(log(ppHealthFac), nt),
    E = rep(map2$ESAM/nt, nt),
    y = unlist(st_drop_geometry(map2)[, 18:27])
)

head(longdf)

ppc <- list(prior = "pc.prec", param = c(0.5, 0.01))

m0 <- y ~ hf +
    f(t, model = 'rw1', scale.model = TRUE, hyper = list(prec = ppc)) +
    f(i, model = 'iid', hyper = list(prec = ppc))


fit0 <- inla(
    formula = m0,
    data = longdf,
    family = 'poisson', E = E
    )

fit0$summary.fixed

exp(fit0$summary.fixed[2, c(1,3,5)])

par(mar = c(4,4,1,1))
plot(inla.smarginal(fit0$marginals.fixed$hf), type = "l", bty = 'n')
abline(v = 0)

fit0$summary.hyperpar

nblist <- spdep::poly2nb(map2)
spdep::nb2INLA(file = 'nb2', nb = nblist)

m1 <- y ~ hf +
    f(t, model = 'rw1', scale.model = TRUE, hyper = list(prec = ppc)) +
    f(i, model = 'besagproper', graph = 'nb2', hyper = list(prec = ppc))

fit1 <- inla(
    formula = m1,
    data = longdf,
    family = 'poisson', E = E
)

fit1$summary.fixed

plot(inla.smarginal(fit1$marginals.fixed$hf), type = "l", bty = 'n')
lines(inla.smarginal(fit0$marginals.fixed$hf), lwd = 2, lty = 2)
abline(v = 0)

fit1$summary.hyperpar

plot.ts(fit1$summary.random$t$mean)

## spacetime
pcar1 <- list(theta = list(prior = 'pc.cor1', param = c(0.5, 0.7)))
m.st <- y ~ hf +
    f(i, model = 'besagproper', graph = 'nb2', hyper = list(prec = ppc), 
      group = t, control.group = list(model = 'ar1', hyper = pcar1))

fit.st <- inla(
    formula = m.st,
    data = longdf,
    family = 'poisson', E = E
)

fit.st$summary.fixed

par(mar = c(4,4,1,1))
plot(inla.smarginal(fit.st$marginals.fixed$hf), type = "l", bty = 'n')
lines(inla.smarginal(fit1$marginals.fixed$hf), lwd = 2, lty = 2)
lines(inla.smarginal(fit0$marginals.fixed$hf), lwd = 2, lty = 3)
abline(v = 0)

fit.st$summary.hyperpar

##
m.st2 <- y ~ hf +
##  f(t, model = 'rw1', scale.model = TRUE, hyper = list(prec = ppc)) +
  f(i, model = 'besagproper', graph = 'nb2', hyper = list(prec = ppc)) +
  f(i2, model = 'besagproper', graph = 'nb2', hyper = list(prec = ppc), 
    group = t2, control.group = list(model = 'ar1', hyper = pcar1))

longdf$i2 <- longdf$i
longdf$t2 <- longdf$t

fit.st2 <- inla(
  formula = m.st2,
  data = longdf,
  family = 'poisson', E = E
)

fit.st2$cpu.used

fit.st2$summary.fixed

plot(inla.smarginal(fit.st2$marginals.fixed$hf), type = "l", bty = 'n')
lines(inla.smarginal(fit.st$marginals.fixed$hf), lwd = 2, lty = 2, col = 'red')
lines(inla.smarginal(fit1$marginals.fixed$hf), lwd = 2, lty = 2)
lines(inla.smarginal(fit0$marginals.fixed$hf), lwd = 2, lty = 3)
abline(v = 0)

fit.st2$summary.hyperpar

map2$spatial_common <- fit.st2$summary.random$i$mean

plot(map2[,"spatial_common"])

#fit.st2$summary.random$t
#plot(fit.st2$summary.random$t$mean)

strfit <- matrix(fit.st2$summary.random$i2$mean, n)

library(sp)
par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
plot(as(map2, 'Spatial'))
INLAspacetime::stlines(t(strfit), 
        as(map2, "Spatial"), lwd = 3)

estfit <- exp(strfit)

ggarrange(
    gg2 +  geom_sf(aes(fill = estfit[, 1]), map2),
    gg2 +  geom_sf(aes(fill = estfit[, 2]), map2),
    gg2 +  geom_sf(aes(fill = estfit[, 3]), map2),
    gg2 +  geom_sf(aes(fill = estfit[, 4]), map2),
    gg2 +  geom_sf(aes(fill = estfit[, 5]), map2),
    gg2 +  geom_sf(aes(fill = estfit[, 6]), map2),
    gg2 +  geom_sf(aes(fill = estfit[, 7]), map2),
    gg2 +  geom_sf(aes(fill = estfit[, 8]), map2),
    gg2 +  geom_sf(aes(fill = estfit[, 9]), map2),
    gg2 +  geom_sf(aes(fill = estfit[, 10]), map2),
    gg2 +  geom_sf(aes(fill = SAM / Pop0to4), map2)
)

##
cv0 <- inla.group.cv(fit0, num.level.sets = 5)
cv1 <- inla.group.cv(fit1, num.level.sets = 5)
cv.st <- inla.group.cv(fit.st, num.level.sets = 5)
cv.st2 <- inla.group.cv(fit.st2, num.level.sets = 5)

str(cv0,1)
cv0$groups[1:2]

longdf[cv0$groups[[1]]$idx,]

c(fit0$mlik[[1]], fit1$mlik[[1]], fit.st$mlik[[1]], fit.st2$mlik[[1]])
c(sum(cv0$cv), sum(cv1$cv), sum(cv.st$cv), sum(cv.st2$cv))
