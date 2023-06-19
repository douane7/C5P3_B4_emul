#10AUG14

"cartosurf" <- function(lon,lat,champ,limits=NA,zlim=NA,option="p",coldir=1,cex=NA,name="",screen=TRUE)
# lon, lat two vector of longitudes and latitudes for plotting on the map
# champ = vector of z-values to plot with color code on the map
# plot a map with surface contours (option="s") or coloured points (option="p") surface contours + points (option="sp")
# if coldir=1: the palette goes from blue to red                               
# if coldir=-1: the palette goes from red to blue                              
# if coldir=2: the palette goes from blue to red through white                               
# if coldir=-2: the palette goes from red to blue through white                           
# limits[4] = vector of limits of the map (minx, maxx, miny, maxy)             
# 	if this vector is NULL, they are calculated on the data (lon, lat)           
# zlim = two-values of min and max for champ
# cex = relative size of the coloured points: 1 = standard size, any value between 0.1 and 10
# name = character chain (between quotes) for the title of the graphics and for the graphic name if screen=FALSE
# screen=TRUE graphics on the screen, otherwise on a pdf file
{
require(fields)
# elimination des données manquantes
lon=lon[!is.na(champ)]
lat=lat[!is.na(champ)]
champ=champ[!is.na(champ)]

# initialisation
if (is.na(limits[1])) {
    minx<-min(lon)
    miny<-min(lat)
    maxx<-max(lon)
    maxy<-max(lat)
} else {
    minx<-limits[1]
    maxx<-limits[2]  
    miny<-limits[3]   
    maxy<-limits[4]
}
if (is.na(zlim[1])) {
    minz<-min(champ)
    maxz<-max(champ)
    zlim=c(minz,maxz)
} else {
    minz<-zlim[1]
    maxz<-zlim[2]  
}

# couleurs
ncolor<-50
if (!screen) pdf(file=paste("map_",name,".pdf",sep=""))
if (abs(coldir)==1) {coli<-rainbow(ncolor,start=0,end=0.8)} 
if (abs(coldir)==2) {coli<-colorRampPalette(c("dark red","red","orange","white","white","cyan","blue","dark blue"))(ncolor)}
indx<-seq(ncolor,1,by=-1)
cold<-coli[indx]

if (coldir>=0) {cols=cold}
else {cols=coli}

# draw surface
ccex=min(sqrt((max(lon)-min(lon))*(max(lat)-min(lat))/length(champ)),1.8)
if (is.na(cex))cex=ccex
xmoy=champ
xmoy[xmoy>maxz]=maxz
xmoy[xmoy<minz]=minz
qq=seq(zlim[1],zlim[2],length=length(cols))
l.col<-c()
for(i in 1:length(champ)) {
	l.col=c(l.col,which(xmoy[i]>=qq[1:(length(qq)-1)] & xmoy[i]<=qq[2:length(qq)]))
}
if (option=="p" ) {
    plot(seq(minx,maxx,length=10),seq(miny,maxy,length=10), type="n",xlab="longitude",ylab="latitude",main=name)
    points(lon,lat,cex=cex,col=cols[l.col],bg=cols[l.col],pch=19)
    world(ylim=c(miny,maxy),xlim=c(minx,maxx),add=TRUE)
  	image.plot(xmoy,col=cols,legend.only=TRUE, graphics.reset=T,zlim=range(qq))
} else {
    fit<-Tps(cbind(lon,lat),champ)
    out.p=predictSurface(fit)
    out.p$z[out.p$z>maxz]=maxz
    out.p$z[out.p$z<minz]=minz
    surface(out.p,type="I",xlab="longitude",ylab="latitude", main=name,ylim=c(miny,maxy),xlim=c(minx,maxx), zlim=c(minz,maxz), graphics.reset=T,col=cols)
    world(ylim=c(miny,maxy),xlim=c(minx,maxx),add=TRUE)
    if (option=="sp") {points(lon,lat,cex=cex,col='black',bg=cols[l.col],pch=21)}
}
message="well done"
if (!screen) dev.off()
return(message)
}
