## GOAL: 
## 1. implement different disease mapping models using
## - a linear covariate
## - a spatial random effect 

### load map
library(rgdal)
map <- readOGR('./', 'Scotland')

class(map)
names(map)
slotNames(map)
names(map@data)
head(map@data)

plot(map)

n <- nrow(map)
n

## computing the SMR
## (ratio between the observed and expected number of cases)
map$SMR <- map$O/map$E

library(gridExtra)
grid.arrange(spplot(map, 'X'), 
             spplot(map, 'SMR'), nrow=1)

plot(map$X, map$SMR)

## library spdep for the neighborhood list
library(spdep)
nb <- poly2nb(map)
nb
nb[[1]]

par(mar=c(0,0,0,0))
plot(map)
plot(nb, coordinates(map), add=TRUE)

args(nb2INLA)
nb2INLA('graph.Scotland', nb)

## model 0
f0 <- O ~ 1

## model 1
## log(smr_i) ~ beta_0 + beta_1 X_i
f1 <- O ~ X

## model 2
## log(smr_i) ~ beta_0 + s_i
map@data$area <- 1:nrow(map)
f2 <- O ~ f(area, model='bym2', graph='graph.Scotland',
      hyper=list(theta1=list(
                     prior='pc.prec',
                     param=c(0.5, 0.1)),
                 theta2=list(
                     prior='pc',
                     param=c(0.5, 0.5))))

## model 3
## log(smr_i) ~ beta_0 + beta_1 X_i + s_i 
f3 <- O ~ X +
    f(area, model='bym2', graph='graph.Scotland',
      hyper=list(theta1=list(
                     prior='pc.prec',
                     param=c(0.5, 0.1)),
                 theta2=list(
                     prior='pc',
                     param=c(0.5, 0.5))))

### actually fit the models
library(INLA)

### some INLA settings
inla.setOption(
    smtp='pardiso',
    pardiso.license='~/.pardiso.lic')

r0 <- inla(f0, family='poisson', data=map@data, E=E,
           control.compute=list(cpo=TRUE))

r1 <- inla(f1, family='poisson', data=map@data, E=E,
           control.compute=list(cpo=TRUE))

r2 <- inla(f2, family='poisson', data=map@data, E=E,
           control.compute=list(cpo=TRUE))

r3 <- inla(f3, family='poisson', data=map@data, E=E,
           control.compute=list(cpo=TRUE))

### CPO from each model 
c(r0=-sum(log(r0$cpo$cpo)),
  r1=-sum(log(r1$cpo$cpo)),
  r2=-sum(log(r2$cpo$cpo)),
  r3=-sum(log(r3$cpo$cpo)))

### fixed effects summary from some models
round(r1$summary.fixed, 2)
round(r3$summary.fixed, 2)

### hyperparameter summary from some models
round(r2$summary.hy, 2)
round(r3$summary.hy, 2)

### plot the covariate effect
par(mar=c(3,3,1,1), mgp=c(2,1,0))
plot(r1$marginals.fixed[[2]],
     xlim=c(0.01, 0.08), type='l')
lines(r3$marginals.fixed[[2]], col='red')

### visualize the spatial log risk
map$spatial <-
    r3$summary.random$area$mean[1:n]
spplot(map, 'spatial')

### the fitted relative risk
map$rr.mean <-
    r3$summary.fitted.values$mean
spplot(map, 'rr.mean')

###
with(map@data, plot(SMR, rr.mean, pch=19, log='xy'))
abline(0:1)

### the PIT values map 
map$r3.pit <- r3$cpo$pit 
spplot(map, 'r3.pit')

map$rr.smr <- map$rr.mean/map$SMR

grid.arrange(
    spplot(map, 'X', main='X'),
    spplot(map, 'SMR', main='SMR'),
    spplot(map, 'spatial', main='spatial'),
    spplot(map, 'rr.mean', main='rr.mean'),
    spplot(map, 'r3.pit', main='r3.pit'),
    spplot(map, 'rr.smr', main='rr.smr'), ncol=3)

