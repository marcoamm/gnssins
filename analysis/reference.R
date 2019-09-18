rm(list=ls()) 
options(digits=12)
#setwd('/home/marco/Documents/code/gnssins/analysis')


xyz2lla <- function(x, y, z) {
    dtr <-  pi / 180.0

    EARTH <- wgs84()
    esq <- EARTH$Esq

    rp <- sqrt(x ^ 2 + y ^ 2 + z ^ 2)

    flatgc <- asin(z / rp) / dtr

    testval <- abs(x) + abs(y)
    if (testval < 1.0e-10)
        flon <- 0.0
    else
        flon <- atan2(y, x) / dtr

    if (flon < 0.0)
        flon <- flon + 360.0

    p <- sqrt(x ^ 2 + y ^ 2)

    # on pole special case
    if (p < 1.0e-10) {
        flat <- 90.0

        if (z < 0.0)
            flat <- -90.0

        altkm <- rp - rearth(flat)

        return(cbind(flat, flon, altkm*1000))  # i dont even know how to test this
    }

    # first iteration, use flatgc to get altitude
    # and alt needed to convert gc to gd lat.

    rnow  <- rearth(flatgc)
    altkm <- rp - rnow
    flat  <- gc2gd(flatgc, altkm)

    rrnrm <- radcur(flat)
    rn    <- rrnrm[2]

    #for (i in 1:5) {  # why only five times?
    while (TRUE) {
        slat  <- sin(dtr * flat)
        tangd <- (z + rn * esq * slat) / p
        flatn <- atan(tangd) / dtr

        dlat <- flatn - flat
        flat <- flatn
        clat <- cos(dtr * flat)

        rrnrm <- radcur(flat)
        rn    <- rrnrm[2]

        altkm <- (p / clat) - rn

        if (abs(dlat) < 1.0e-12)
            break
    }

    return(cbind(flat, flon, altkm*1000))
}

llh2xyz <- function(lat,long, h) {

    a = 6378137;
    e = 0.0818191910435;

    lat = lat/180*pi;   # converting to radians
    long = long/180*pi; # converting to radians

    e2 = e^2;

    chi = sqrt(1-e2*(sin(lat))^2); 
    X = (a/chi +h)*cos(lat)*cos(long); 
    Y = (a/chi +h)*cos(lat)*sin(long); 
    Z = (a*(1-e2)/chi + h)*sin(lat);

    return(cbind(X,Y,Z))

}

xyz2enu <- function(phi,lambda,x,y,z) {

    phi = phi*pi/180;
    lambda = lambda*pi/180;
    cl = cos(lambda);  
    sl = sin(lambda);
    cb = cos(phi);
    sb = sin(phi);

    F = rbind(c(-sl, -sb*cl, cb*cl), c(cl, -sb*sl, cb*sl), c(0, cb, sb))

    local_vect = t(F)%*%rbind(x, y, z);
    e = local_vect[1,];
    n = local_vect[2,];
    u = local_vect[3,];

    return(cbind(e,n,u))

}

ref <- read.table('../data/19032019/reference.pos');
pos <- read.table('../out/out_PVA.txt');

#Clean PVA for integration epochs only
pos <- pos[pos[,11]==1,];

ref.time <- round(ref[,2]);
pos.time <- round(pos[,1]);

#Match time series by time stamps
index <- match(round(pos[,1]),round(ref[,2]))
ref <- ref[index,];
ref.time <- round(ref[,2]);



#Use each point as llh reference 
#(accounting for earth's curvature in long trajectories)
refxyz <- llh2xyz(ref[,3],ref[,4], ref[,5])
posxyz <- llh2xyz(pos[,2],pos[,3], pos[,4])

#ENU differences from solution to reference

enu <- xyz2enu(mean(ref[,3]),mean(ref[,4]),posxyz[,1]-refxyz[,1],posxyz[,2]-refxyz[,2],posxyz[,3]-refxyz[,3])
#phi=mean(ref[,3]); lambda=mean(ref[,4]); x=posxyz[,1]-refxyz[,1]; y=posxyz[,2]-refxyz[,2]; z=posxyz[,3]-refxyz[,3];

hz=sqrt(enu[,1]^2+enu[,2]^2)
up=enu[,3]

jpeg("Horizontal.jpg", width = 1400, height = 700)
plot(ref.time,hz,type="l",lwd=3,col="blue",ylim=c(mean(hz)-3,mean(hz)+3),
     main="Horizontal differences",ylab="m", xlab="GPS time (SoW)",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
grid()
dev.off();

jpeg("Hist_Horizontal.jpg", width = 700, height = 700)
hist(hz[hz<=30],col="blue",main="Horizontal differences",xlab="m", xlim=c(mean(hz)-7,mean(hz)+7),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,breaks=30)
grid()
dev.off();

jpeg("Vertical.jpg", width = 1400, height = 700)
plot(ref.time,up,type="l",lwd=3,col="red",ylim=c(mean(up)-3,mean(up)+3),
     main="Vertical differences",ylab="m", xlab="GPS time (SoW)",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
grid()
dev.off();

jpeg("Hist_Vertical.jpg", width = 700, height = 700)
hist(up,col="red",main="Vertical differences",xlab="m", xlim=c(mean(up)-7,mean(up)+7),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,breaks=30)
grid()
dev.off();


write("Analysis summary:\n", "summary.txt")
write("Horizontal average bias: ", "summary.txt",append=TRUE)
write(paste(mean(hz)," m\n"), "summary.txt",append=TRUE)
write("Horizontal RMS: ", "summary.txt",append=TRUE)
write(paste(sd(hz)," m"), "summary.txt",append=TRUE)

