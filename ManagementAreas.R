###############################################################
#
#  Author:        Matt Grinnel
#  Contact:       grinnellmatt@gmail.com
#  Institution:      
#  Project:       
#  Date started:  May 11, 2011
#  Date modified: 
#  Goal:          
#
###############################################################


########################
##### Housekeeping #####
########################

rm( list=ls() )
graphics.off( ) 
sTime <- Sys.time( )

tolerance = 0.5   # in km

library( PBSmapping ) ;  library( maptools ) ;  library( sp )
#install.packages( c("PBSmapping", "maptools", "sp") )
data(nepacLL)

# get the contour line data
#write.table( x=bathy, append=F, sep="\t", col.names=T, row.names=F,
#  file="bcBath.txt" )
#bathy = read.table( file="bcBath.txt", head=TRUE, sep="\t" )
#bathyTop = makeTopography( bathy )
#save( bathyTop, file="bathyTop.RData" )
load( file="Figures/bathyTop.RData" )
bathContours = c( 50, 100, 200, 500, 1000 ) 
bathyCon = contourLines( bathyTop, levels=bathContours )
bathyCP  = convCP( bathyCon )    
bcBathy  = bathyCP$PolySet
attr( bcBathy, "projection" ) = "LL"

# Get canada data
canadaLL <- importShapefile( fn="Figures/global_worldborders_shp/world_borders", 
    readDBF=FALSE, projection="LL" )
#canadaUTM <- convUL( xydata=canadaLL, km=TRUE )

# entire NE pacific - this is embedded
nePac = function( ) {
  
  nAmCoast = thinPolys( nepacLL, tol=tolerance*5, filter=3 )
  nAmCoast = as.PolySet( nAmCoast, projection="UTM", zone=9 )
  
  plotPolys( canadaLL, col="grey", bg="white", ann=FALSE, xlab="", ylab="", 
      xaxt="n", yaxt="n", xlim=c(-180, max(nepacLL$X)), 
      ylim=c(40, max(nepacLL$Y)), border="black" )
  
  text( x=-150, y=65.5, "Alaska,\nU.S.A.", cex=0.8 )
  text( x=-118, y=54.2, "British\nColumbia,\nCanada", cex=0.8 )
  text( x=-144, y=55, "Gulf of\nAlaska", cex=0.8 )
  text( x=-162, y=46, "North-East\nPacific Ocean", cex=1 )
  text( x=-133, y=45.5, "Area\nenlarged", cex=0.7 )
  text( x=-117, y=43.5, "U.S.A.", cex=0.8 )
  lines( x=c(lef,-125), y=c(top,top), lwd=2 )
  lines( x=c(lef,rig), y=c(bot,bot), lwd=2 )
  lines( x=c(lef,lef), y=c(bot,top), lwd=2 )
  lines( x=c(rig,rig), y=c(bot,50.50), lwd=2 )
  box( lwd=2, col="black" )
  #addLines( canadaUTM, col="black" )
}

# get the trawl area data; either from the shp file, or the saved data object
#bcAreas = importShapefile( fn="TrawlAreasSmall/Trawl_Areas.shp", 
#  projection="UTM", zone=9 )
#bcAreas = convUL( bcAreas, km=TRUE )
#save( bcAreas, file="bcAreas.RData" )
load( file="Figures/bcAreas.RData" )

# remove data that is right along the coast
bcCoast = thinPolys( nepacLL, tol=tolerance, filter=3 )
bcCoast = as.PolySet( bcCoast, projection="LL", zone=9 )

# set area
top=54.75;  bot=47.75;  lef=-135;  rig=-119

# Distance scale (100 km)
distLL <- data.frame( PID=c(1, 1), POS=c(1, 2), X=c(-120.6, -119.21), Y=c(49.6, 49.6) )
distLL <- as.PolySet( x=distLL, projection="LL", zone=10 )

pdf( width=8, height=6, file="Figures/ManagementAreas.pdf" )
#windows( width=8, height=6 )
par( omi=c(0,0,0,0) )

#  plotMap( bcAreas, xlim=c(lef, rig), ylim=c(bot, top), lwd=2, tck=-0.015, 
#    plt=c(.07,.98,.07,.98) )
#  addLines( bcBathy, col="lightgrey" )
#  addPolys( bcCoast, col="grey" )

myPal = palette( gray(0:length(bathContours) / length(bathContours)) )
myPal = myPal[ 1:length(bathContours) ]

# Plot the map
plotMap( bcBathy, xlim=c(lef, rig), ylim=c(bot, top), lwd=2, tck=-0.015, 
    plt=c(.07,.98,.07,.98), type="n", cex.lab=1.5, las=1, xlab="", ylab="", 
    xaxt="n" )
mtext( side=1, line=2.5, expression(paste("Longitude ", "(",degree,"W)", 
            sep="")) )
mtext( side=2, line=2.7, expression(paste("Latitude ", "(",degree,"N)", 
            sep="")) )
axis( side=1, at=c(-130, -125, -120), tick=TRUE, labels=c("130", "125", 
        "120"),	tck=-0.015 )
axis( side=1, at=-119:-135, labels=rep("", times=17), tck=-0.0075 ) 

addLines( bcBathy, col=myPal, lwd=1 )  
addLines( bcAreas, col="black", lwd=2 )
addPolys( bcCoast, col="grey", border="black" )
arrows( x0=distLL$X[1], y0=distLL$Y[1], x1=distLL$X[2], y1=distLL$Y[2], 
    lwd=1, col="black", angle=90, code=3, length=0.05 )
#addLines( distLL, col="black", lwd=2 )
text( x=-119.925, y=49.7, "100 kilometres", cex=0.52 )

# Two legends (directly on top) to fit the title, which has 2 lines
legend( "bottomright", legend=rep(max(bathContours), times=length(myPal)+1), 
    title="bathymetry", bg="white", col="white", lty="solid", cex=0.75, 
    text.col="white", title.col="white" )
legend( "bottomright", legend=bathContours, title="bathymetry\n(metres)", 
    col=myPal, lty="solid", cex=0.75, bty="n" )
box( lwd=2, lty="solid", col="black" )
#plotMap( bcCoast, xlim=c(lef, rig), ylim=c(bot, top), lwd=1, tck=-0.015, 
#  plt=c(.07,.98,.07,.98), col="grey" )
#points( x=-125.6, y=49, pch="x", col="red", cex=3 )
#text( x=-128, y=48.7, "Pirate Booty", cex=2, font=2, col="red" )

# Management areas
text( x=-128.5, y=48.4, "3C", cex=1.5, font=2, col="black" )
text( x=-130.5, y=49.7, "3D", cex=1.5, font=2, col="black" )
text( x=-122.9, y=50.2, "4B", cex=1.5, font=2, col="black" )
text( x=-131, y=50.9, "5A", cex=1.5, font=2, col="black" )
text( x=-131.5, y=51.6, "5B", cex=1.5, font=2, col="black" )
text( x=-134, y=53, "5E", cex=1.5, font=2, col="black" )

# Rectangles under 5C and 5D
rect( xleft=-130.85, ybottom=52.55, xright=-130.15, ytop=52.85, 
    col=rgb(1,1,1,0.9), border=NA )
text( x=-130.5, y=52.7, "5C", cex=1.5, font=2, col="black" )
rect( xleft=-131.85, ybottom=54.25, xright=-131.15, ytop=54.55, 
    col=rgb(1,1,1,0.9), border=NA )
text( x=-131.5, y=54.4, "5D", cex=1.5, font=2, col="black" )



# Vancouver Island
text( x=-125.45, y=49.54, "Vancouver    Isl.", srt=-43, cex=1.5, col="black" )

# North arrow
xAr = -119.9 ;  yAr = 50.65
points( x=xAr, y=yAr, cex=3, pch=24, col="black", bg="black" )
lines( x=c(xAr, xAr), y=c(yAr-0.5, yAr), lwd=3, col="black" )
text( x=xAr, y=yAr+0.5, col="black", "N", cex=1.5, font=2 )

# Point to 4B area
arrows( x0=-123.6, x1=-123.7, y0=49.35, y1=49.2, lwd=4, code=2, length=0.125,
    col=rgb(1,1,1,0.9) )  
arrows( x0=-123.2, x1=-123.7, y0=50, y1=49.2, lwd=2, code=2, length=0.125,
    col="black" )

# Add inset
par( fig=c(x1=0.45, x2=0.99, y1=0.45, y2=0.98), new=T )
nePac( )
dev.off()

cat("Done\n")
