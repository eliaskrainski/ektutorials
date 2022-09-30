if(!any(ls()=='plotting'))
    plotting <- TRUE

### load package to fast read large dataset
library(data.table)

### stations data
### https://www.ncei.noaa.gov/pub/data/ghcn/daily/
url0 <- 'https://www.ncei.noaa.gov/pub/data/ghcn/daily/'
stfl <- 'ghcnd-stations.txt'

if(FALSE) 
    browseURL(url0)

if(!file.exists(stfl)) {
    options(timeout=10*60) ### for bad Internet connection
    download.file(paste0(url0, stfl), stfl)
}

### width of the colums in the file:          
(ws <- diff(c(0,11,20,30,37,71,75,79,85)))

### read station information: longitute, latitude & altitude information
stations <- read.fwf(stfl, ws[1:4])
colnames(stations) <- c('station', 'latitude', 'longitude', 'elevation')

###
cat('Readed', ncol(stations), 'variables from',
    nrow(stations), 'stations worldwide\n')

### deal with sp and projection
library(sp)
coordinates(stations) <- ~ longitude + latitude
stations@proj4string <- CRS('+proj=longlat +datum=WGS84')

### index of stations with 'US' code
ii0us <- which(substr(stations$station, 1, 2)=='US')
table(substr(stations$station, 1, 3)[ii0us])

### define a box around US main territory
rect.ll <- SpatialPolygons(list(Polygons(list(Polygon(
    cbind(c(-130, -60, -60, -130, -130),
          c(50, 50, 23, 23, 50)))), '0')),
    proj4string=stations@proj4string)

### select stations with 'US' code and inside the rectangle 
ii1us <- ii0us[which(!is.na(over(stations[ii0us,], rect.ll)))]

### projection with units in km
stations.mkm <- spTransform(
    stations[ii1us, ], CRS('+proj=moll +units=km'))

##
cat('Found', nrow(stations.mkm), 'stations in US mainland\n')

### Visualize the selected stations
if(plotting) {
    par(mar=c(0,0,0,0))
    plot(stations.mkm, pch=8, cex=0.1)
}

### download data
dfl <- '2022.csv.gz'
if(!file.exists(dfl)) {
    options(timeout=30*60) ### for bad Internet connection
    download.file(paste0(url0, 'by_year/', dfl), dfl)
}

### read the data
system.time(alld <- as.data.frame(fread(dfl)))
dim(alld)
head(alld,2)

###
cat('Readed', nrow(alld), 'weather data at stations worldwide\n')

### id of data on TMAX
id0 <- which((alld$V6 == '') &
             (alld$V3 == 'TMAX'))

### id to stations in the rectangle
id1 <- id0[which(alld$V1[id0]%in%stations.mkm$station)]

###
cat('Found', length(id1), 'TMAX data in US mainland\n')

### reshape the data to select stations looking over time series
system.time(dataw <- reshape(
                data=alld[id1, c('V1','V2','V4')],
                idvar='V1', timevar='V2', direction='wide'))

dim(dataw)
dataw[1:3, 1:5]

colnames(dataw)[1] <- 'station'
colnames(dataw) <- gsub('V4.', '', colnames(dataw), fixed=TRUE)
dataw[1:3, 1:5]

### some stations have no good quality data
### usually these stations also have several missing data
### select stations with at least 200 days without missing data
ndays <- ncol(dataw)-1
nds <- rowSums(!is.na(dataw[,-1]))

ndt <- colSums(!is.na(dataw[,-1]))

ndt[201:length(ndt)]

dates <- as.Date('2021-12-31') + 1:ndays

summary(nds)

tail(table(nds), 15)

id.wsel <- which(nds>240)
ns <- length(id.wsel)
c(length(nds), ns)

### a time series plot
tmax.t <- colMeans(dataw[id.wsel, -1], na.rm=TRUE)/10

if(plotting) {
    par(mfrow=c(1,1), mar=c(3,3,0.5,0.5), mgp=c(2,1,0))
    plot(dates, tmax.t,     
         type='l', xlab='', ylab='Daily average of TMAX')
}

dates[tsel <- which.max(tmax.t)]
tmax.t[-5:5 + tsel]

summary(dataw[id.wsel, 1+tsel]/10)

### organize one day data
##which(stations.mkm$station %in% dataw[id.wsel,1])
id.ss <- pmatch(dataw[id.wsel,1], stations.mkm$station)
summary(id.ss)

tmax1day <- stations.mkm[id.ss, ]
tmax1day$tmax <- dataw[id.wsel, 1+tsel+1]/10

###
cat('Selected TMAX from', nrow(tmax1day), 'stations\n')


if(plotting) {
    spplot(tmax1day, 'elevation')
    spplot(tmax1day, 'tmax')
}

obj.not.rm <- c('map.moll', 'mesh', 'grid.proj', 'id.grid.out',
                'dataw', 'tmax1day')
if(!FALSE)
    rm(list=setdiff(ls(), obj.not.rm))

gc(reset=TRUE)
