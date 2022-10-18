## https://github.com/eliaskrainski/ektutorials
### load script to build the mesh
system.time(source('us_mesh.R'))

### load script to get the data
system.time(source('us_get_tmax1day.R'))

### summary of the data
summary(tmax1day)
sd(tmax1day$tmax, na.rm=TRUE)

### Construct latent model components
matern <- inla.spde2.pcmatern(
    mesh=mesh, alpha=2, 
    prior.sigma = c(5, 0.01), 
    prior.range = c(25, 0.01))

### load inlabru (easier to code the model) 
library(inlabru)

### define the model
model <- tmax ~ Intercept(1) +
    field(main=coordinates, model = matern)

### prior for the likelihood parameter 
lik.prec <- list(prec=list(prior='pc.prec', param=c(5, 0.01)))

### set some INLA parameters
inla.setOption(
    inla.mode='experimental',
    num.threads='4:-1',
    smtp='pardiso', 
    pardiso.license='~/.pardiso.lic')

### fit the model using inlabru
fit <- bru(
    model, tmax1day, family='gaussian',
    options=list(
        verbose=TRUE,
        control.family=list(hyper=lik.prec)))

### some summary 
fit$summary.fix
fit$summary.hyperpar

### consider the posterior mean of the random field 
s.mean <- fit$summary.ran$field$mean

### project it into a grid (for plotting) 
y.m <- inla.mesh.project(grid.proj, field=s.mean)
y.m[id.grid.out] <- NA

library(fields)

### visualize the random field + b0
par(mfrow=c(1,1), mar=c(0,0,0,0))
image.plot(
    x=grid.proj$x,
    y=grid.proj$y,
    z=y.m+fit$summary.fix$mean[1], asp=1)
points(tmax1day, cex=0.05, pch=8)
plot(map.moll, add=TRUE, border=gray(0.3,0.5))

### group cross validation
system.time(gcpo <- inla.group.cv(fit, 20))

### selected locations to visualize
isel <- c(867, 1349, 1914, 2114, 2618, 3055, 3658, 4608, 4666, 5060, 5348)

### number of neighbors (with m=10) at the selected data locations
nnb <- sapply(gcpo$groups[isel], function(x) length(x$idx)-1)
nnb

### plot the neighbors for some data points
locs <- coordinates(tmax1day)
par(mfrow=c(1,1), mar=c(0,0,0,0))
image.plot(
    x=grid.proj$x,
    y=grid.proj$y,
    z=y.m+fit$summary.fix$mean[1], asp=1)
plot(map.moll, add=TRUE, border=gray(0.3,0.5))
points(tmax1day, cex=0.5, pch=8)
for(i in isel) {
    jj <- gcpo$groups[[i]]$idx[-1]
    segments(locs[i, 1], locs[i, 2], locs[jj, 1], locs[jj, 2])
    points(locs[jj, ], pch=19, cex=1, col='white')
}
points(locs[isel, ], pch=19, cex=3, col='white')
text(locs[isel, 1], locs[isel, 2], paste(nnb), col='blue3', cex=.8)

if(FALSE) {

    ll <- locator()
    isel <- sapply(1:length(ll[[1]]), function(i)
        which.min(sqrt((locs[,1]-ll$x[i])^2 +
                   (locs[,2]-ll$y[i])^2)))
    isel <- sort(isel)
    isel
}
