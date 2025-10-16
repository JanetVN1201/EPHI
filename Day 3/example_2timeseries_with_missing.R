nmat <- structure(
    c(29843L, 2134L, 29881L, 2254L, 29049L, 2140L, 29620L, 
      2151L, 29031L, 2067L, 27108L, 2008L, 26371L, 2028L, 24856L, 1939L, 
      25209L, 1889L, 24547L, 1909L, 24735L, 1883L, 24379L, 1932L, 25329L, 
      1918L, 24868L, 2010L, 25345L, 2089L, 25249L, 2249L, 25079L, 2285L, 
      24932L, 2336L, 24815L, 2465L, 24506L, 2460L, 23217L, 2255L, 22745L, 
      2358L, 22112L, 2207L, 21393L, 2262L, 19728L, 2188L, 18575L, 2034L, 
      18396L, 1951L, 17978L, 2009L),
    dim = c(2L, 28L),
    dimnames = list(
        c("Curitiba", "Araucária"),
        c("1996", "1997", "1998", "1999", "2000", "2001", "2002",
          "2003", "2004", "2005", "2006", "2007", "2008", "2009",
          "2010", "2011", "2012", "2013", "2014", "2015", "2016",
          "2017", "2018", "2019", "2020", "2021", "2022", "2023"
          )))

dmat <- structure(
    c(538L, 60L, 476L, 44L, 483L, 47L, 436L, 49L, 436L, 
      31L, 370L, 27L, 311L, 29L, 308L, 38L, 283L, 35L, 293L, 31L, 259L, 
      19L, 256L, 18L, 250L, 28L, 223L, 23L, 231L, 22L, 222L, 33L, 238L, 
      31L, 220L, 21L, 192L, 31L, 220L, 21L, 201L, 16L, 189L, 26L, 183L, 
      26L, 139L, 22L, 140L, 14L, 135L, 22L, 158L, 16L, 142L, 19L),
    dim = c(2L, 28L),
    dimnames = list(
        c("Curitiba", "Araucária"),
        c("1996", "1997", "1998", "1999", "2000", "2001", "2002", 
          "2003", "2004", "2005", "2006", "2007", "2008", "2009", 
          "2010", "2011", "2012", "2013", "2014", "2015", "2016",
          "2017", "2018", "2019", "2020", "2021", "2022", "2023")))

nmat
dmat

yy <- 1996:2023

sum(dmat)/sum(nmat)

par(mfrow=c(1,1), mar=c(3,4,0,0), mgp=c(3,0.5,0))
plot(yy, dmat[1, ] / nmat[1, ], axes=FALSE, type='o',
     xlab='', ylab='Deaths <1 year old / Born alive',
     pch=19, las=1, ylim=c(0, 0.03))
points(yy, dmat[2, ] / nmat[2, ], col=2, pch=8, type='o')
axis(1, yy, las=2)
axis(2, las=1)
legend('topright', c('Curitiba', 'Araucária'),
       col=1:2, pch=c(19, 8), lty=1, title='1996 - 2023')


na <- 2
nt <- ncol(nmat)

dat <- data.frame(
    iarea = rep(1:na, each = nt),
    itime = rep(1:nt, na), 
    nborn = c(nmat[1, ], nmat[2, ]),
    death = c(dmat[1, ], dmat[2, ])
)

dat

library(INLA)

res <- inla(
    formula = death ~ 0 + f(itime, model='rw2', replicate = iarea, constr = FALSE),
    family='poisson', E = nborn,
    data = dat, 
    control.predictor = list(compute = TRUE),
    control.compute = list(return.marginals.predictor = TRUE))

fitt <- res$summary.fitted.val


par(mfrow=c(1,1), mar=c(3,4,0,0), mgp=c(3,0.5,0))
plot(yy, dmat[1, ] / nmat[1, ], axes=FALSE, type='o',
     xlab='', ylab='Deaths <1 year old / Born alive',
     pch=19, las=1, ylim=c(0, 0.03))
polygon(c(yy, rev(yy), yy[1]),
        c(fitt[1:nt, 3], fitt[nt:1,5], fitt[1,3]),
        col=gray(0.7, 0.5), border=gray(0.7, 0.5))
polygon(c(yy, rev(yy), yy[1]),
        c(fitt[nt+1:nt, 3], fitt[nt+nt:1,5], fitt[nt+1,3]),
        col=rgb(1,0.7,0.5,0.5), border=rgb(1,0.7,0.5,0.5))
points(yy, dmat[2, ] / nmat[2, ], col=2, pch=8, type='o')
axis(1, yy, las=2)
axis(2, las=1)
legend('topright', c('Curitiba', 'Araucária'),
       col=1:2, pch=c(19, 8), lty=1, title='1996 - 2023')

res$summary.fitted.values[1:2, ]

res$summary.fitted.values[c(1, nt+1), ]

res$summary.fitted.values[c(nt, 2*nt), ]

par(mfrow=c(1,1), mar=c(3,4,0,0), mgp=c(3,0.5,0))
plot(inla.smarginal(res$marginals.fitted.val[[nt]], factor=5),
     xlim=c(0.0055, 0.011),
     xlab='Taxa', ylab='', type='l', axes=FALSE)
lines(inla.smarginal(res$marginals.fitted.val[[nt+nt]], factor=5), col=2)
axis(1)
abline(v=dat$y[nt]/dat$n[nt])
abline(v=dat$y[nt*2]/dat$n[nt*2], col=2)
legend('topright', c('Curitiba', 'Araucária'),
       col=1:2, lty=1, title='Infant Mortality rate in 2023')


## suppose we have missing data
dat2 <- dat
nt
dat2$death[5:24] <- NA
dat2$death[nt + 11:22] <- NA
dat2

res2 <- inla(
    formula = death ~ 0 + f(itime, model='rw2', replicate = iarea, constr = FALSE),
    family='poisson', E = nborn,
    data = dat2, 
    control.predictor = list(compute = TRUE,
                             link = 1), ## add to predict at the response scale
    control.compute = list(return.marginals.predictor = TRUE))

fitt2 <- res2$summary.fitted.val


par(mfrow=c(1,1), mar=c(3,4,0,0), mgp=c(3,0.5,0))
plot(yy, dmat[1, ] / nmat[1, ], axes=FALSE, type='o',
     xlab='', ylab='Deaths <1 year old / Born alive',
     pch = ifelse(is.na(dat2$death[1:nt]), 8, 19), las=1, ylim=c(0, 0.03))
polygon(c(yy, rev(yy), yy[1]),
        c(fitt2[1:nt, 3], fitt2[nt:1,5], fitt2[1,3]),
        col=gray(0.7, 0.5), border=gray(0.7, 0.5))
polygon(c(yy, rev(yy), yy[1]),
        c(fitt2[nt+1:nt, 3], fitt2[nt+nt:1,5], fitt2[nt+1,3]),
        col=rgb(1,0.7,0.5,0.5), border=rgb(1,0.7,0.5,0.5))
points(yy, dmat[2, ] / nmat[2, ], col=2, ifelse(is.na(dat2$death[nt+1:nt]), 8, 19), type='o')
axis(1, yy, las=2)
axis(2, las=1)
legend('topright', c('Curitiba', 'Araucária'),
       col=1:2, pch=c(19, 8), lty=1, title='1996 - 2023')
