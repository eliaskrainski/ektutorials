## https://eliaskrainski.github.io/tutorials/ICoMCoS2022.zip
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
