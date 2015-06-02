#---------------------------------------------------------------------
# 
# REM612 Project Script - Talk Plots
# 
# Author: Samuel Johnson
# Date: 19/1/15
# 
# The following script will be used to produce plots for the talk
# on the sablefish IBM.
# 
# ToDo: Spawning Grounds
#       Larval drift
#       Migration lines
#       Trawl polygon
#       Vector Field
#       Habitat Polygons
#       Tulip plot of habitat estimation performance
#       Whisker plot of detections.
# 
#---------------------------------------------------------------------

rm ( list = ls ( ) )

# Load required packages
require ( animation ); require ( data.table ); require ( adehabitat );
require ( parallel ); require ( dplyr );

# Load functions
source ( "functions.R" )

# Load detections control file for RData file paths
detCtl <- lisread ( "detections.ctl" )
simCtl <- lisread ( "simulator.ctl" )

# Uncomment to load data from simulator
load ( file = detCtl $ inputFish )
load ( file = detCtl $ inputBoats )
load ( file = detCtl $ inputDetections )
load ( file = "bathyTop.RData" )

# Okay, let's write some functions to estimate habitat size, and overlay plots
# on the map.

# First, get the fish numbers for each inlet (this was done poorly, but time
# is a factor)
inlets <- simCtl $ inlets
nCores <- detectCores ( ) - 1

# Run function to append inlet numbers
appendedList <- appendInlets ( fArr = F, inletLocs = inlets, 
                               detList = yearDetections )

# Pull out appended lists and arrays
F <- appendedList $ F
yearDetections <- appendedList $ detList
inletFishIdx <- appendedList $ inletFishIdx


# Okay, now let's go from the top
spawnPlot <- "bathyMapSpawning.png"
png ( filename = spawnPlot, width = 1000, height = 800 )
image ( bathyTop, col = bgPal )
polygon ( x = c ( -134, -135, -127, -126 ), y = c ( 55, 52, 47, 48 ), 
          col = "darkred", density = 25 )
dev.off ( )

# Larval Drift
inletsMap <- "inletsMap.png"
png ( filename = inletsMap, width = 1000, height = 800 )
image ( bathyTop, col = bgPal )
points ( x = inlets [ , 1], y = inlets [ , 2 ], cex = 5, col = "gold", 
         lwd = 3 )
dev.off()

# Trawl polygon with both of those
trawlPoly <- "trawlMap.png"
png ( filename = trawlPoly, width = 1000, height = 800 )
image ( bathyTop, col = bgPal )
points ( x = inlets [ , 1], y = inlets [ , 2 ], cex = 5, col = "gold", 
         lwd = 3 )
polygon ( x = c ( -134, -135, -127, -126 ), y = c ( 55, 52, 47, 48 ), 
          col = "darkred", density = 25 )
polygon ( x = c ( -126, -125, -131, -134 ), y = c ( 47, 48.5, 55, 54 ),
          col = gray70, density = 25 )
dev.off()

# Vector Field
vectorField <- "vectorField.pdf"
pdf ( file = vectorField,  width = 1000, height = 800 )
image ( bathyTop, col = bgPal )
nX <- length ( bathyTop $ x )
nY <- length ( bathyTop $ y )
xIdx <- seq ( from = 1, to = nX, by = 15 )
yIdx <- seq ( from = 1, to = nY, by = 15 )
for ( x in xIdx )
{
  for ( y in yIdx )
  {
    X <- bathyTop $ x [ x ]
    Y <- bathyTop $ y [ y ]
    dX <- -1 * bathyGrad $ dX [ x, y ]
    dY <- -1 * bathyGrad $ dY [ x, y ]

    arrows ( x0 = X, y0 = Y, x1 = X + dX, y1 = Y + dY, length = 0.05, angle = 35,
              col = "gray30", lwd = 1 )
  }
}
dev.off ( )

# Plot of detections
detections <- yearDetections [[ 1 ]]
fish <- sample ( 1:10000, size = 200 )
setkey ( detections, fish.number)
detFish <- detections [ J(fish) ]
plotDetect <- "detectionsPlot.png"
png ( filename = plotDetect, width = 1000, height = 800 )
image ( bathyTop, col = bgPal )
points ( x = detFish $ boat.lat, y = detFish $ boat.lon, cex = 0.01, 
          col = "brown" )
dev.off ()

# Now to plot average number of detections
detections <- group_by ( detections, inlet, fish.number )
numDetections <- summarise ( detections, n = n() )
par ( bg = "white" )
boxplot ( n ~ inlet, data = numDetections, main = "Detection Distributions",
          xlab = "Inlet Number", ylab = "Detections", las = 1, bg = "white" )
  abline ( h = mean ( numDetections $ n ), col = "red", lwd = 2 )
  abline ( h = median ( numDetections $ n ), col = "steelblue", lwd = 2 )


# Quantile plots of habitat estimator performance
# sampSizes <- c ( 10, 25, 50, 100, 250)
# load ( file = "sample.RData")
# quantList <- vector ( mode = "list", length = length ( sampSizes ) )
# quantArr <- array ( 0, dim = c ( length ( sampSizes ), 3, 5 ) )
# for ( i in 1:length ( sampList ) )
# {
#   ratios <- sampList [[ i ]][,3,]
#   quantList [[ i ]] <- apply (  X = ratios, MARGIN = 1, FUN = quantile, 
#                                 probs = c(0.025, 0.5, 0.975 ) )
#   quantArr [ i, , ] <- apply (  X = ratios, MARGIN = 1, FUN = quantile, 
#                                 probs = c(0.025, 0.5, 0.975 ) )
# }
# lineCols <- c ( "gray30", "gold", "darkgreen", "red", "darkblue" )
# polyCols <- c ( "gray70", "yellow", "green", "tomato", "steelblue" )
# yLim <- c ( 0, max ( quantArr ) )
# xLim <- c ( 1, 5 )
# xVtx <- c ( 1:5, 5:1 )
# par ( mfrow = c( 3, 2 ), bg = "white", xaxp = 1:5 )
# for ( i in 1:length ( sampSizes ) )
# {
#   if ( i == 1 ) title <- "Aggregate stock."
#   else title <- paste ( "Inlet ", i - 1, sep = "" )
#   plot ( x = 1:5, y = quantArr [ , 2, i ], col = lineCols [ i ], lwd = 3, 
#        xlim = xLim, ylim = yLim, type = "l", las = 1,
#        ylab = "Percentage Overlap", xlab = "Sample Size", main = title, 
#        xaxt = "n" )
#   axis ( side = 1, at = 1:5, labels = sampSizes )
#   polygon ( x = xVtx, y = c ( quantArr [ 1:5, 1, i], quantArr [ 5:1, 3, i ] ),
#             col = polyCols [ i ], density = 30 )  
# }
# plot ( 1, type = "n", axes = FALSE, xlab = "", ylab = "" )

# legText <- c (  "Aggregate", "Portland Inlet", "Gil Island", "Fjordland", 
#                 "King Island" )
# panLegend ( x = 0.2, y = 0.7, legTxt = legText, col = lineCols, lty = 1, 
#             lwd = 2, bty = "n" )
