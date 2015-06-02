#---------------------------------------------------------------------
# 
# REM612 Project Script - Animation script
# 
# Author: Samuel Johnson
# Date: 19/1/15
# 
# The following script takes the output of the fish and boat simulators
# and compares them at each time step. The fish and boats which come 
# close enough to each other trigger a detection. Detections are then
# collated to estimate fish habitat range and compared to the true
# simulated habitat range.
# 
# ToDo: Animate fish movement
#       Animate boat movement
#       Overlay on background
#       Animate ADEhabitat results
# 
#---------------------------------------------------------------------

# Housekeeping
rm ( list = ls ( ) )

# Load required packages
require ( animation ); require ( data.table ); require ( adehabitat )

# Load functions
source ( "functions.R" )
load ( "bathyTop.RData" )

# Load detections control file for RData file paths
detCtl <- lisread ( "detections.ctl" )
simCtl <- lisread ( "simulator.ctl" )

# Uncomment to load data for animation
load ( file = detCtl $ inputSim )
load ( file = detCtl $ inputDetections )

B <- fishingYears [[ 1 ]] $ B
nB <- dim ( B ) [ 3 ]
nT <- simCtl $ nT


# Let's try to save a short movie of the fish moving around
saveGIF ( 
{
  for ( i in 1:( nT / 50 ) ) 
  { image ( bathyTop, col = bgPal) 
    points (  x = F [ 1:simCtl $ nF, 1, ( i - 1 ) * 50 + (1:5) ],
              y = F [ 1:simCtl $ nF, 2, ( i - 1 ) * 50 + (1:5) ],
              cex = 0.2,
              col = "red" )
    points ( x = B [ 1, ( i - 1 ) * 50 + (1:5), 1:nB ],
             y = B [ 2, ( i - 1 ) * 50 + (1:5), 1:nB ],
             cex = 0.7, col = "black" )
    ani.pause ( )
  }
}, movie.name = "FishMove.gif", interval = 0.2, nmax = 500, ani.width = 1000,
ani.height = 800 )

# Load in polygon data for animations
load ( file = "polyData.RData" )
image ( bathyTop, col = bgPal )
str ( polys [[ 1 ]] )


habColours <- c ( "gray40", "gold" )
sampSizes <- c ( 10, 25, 50, 100, 250 )

saveGIF ( 
{
  for ( i in 1:length ( polys ) ) 
  { 
    titTxt <- paste ( "Habitats with ", sampSizes [ i ], 
                      " fish tagged in each inlet.", sep = "" )
    image ( bathyTop, col = bgPal, main = titTxt) 
    for ( j in 1:2 )
    {
      detPoly <- polys [[ i ]] [[ j ]] [[ 1 ]]
      simPoly <- polys [[ i ]] [[ j ]] [[ 2 ]]
      polygon ( x = detPoly $ X, y = detPoly $ Y, col = habColours [ j ],
                density = 50, border = NA )
      polygon ( x = simPoly $ X, y = simPoly $ Y,
                density = 0, border = habColours [ j ], lwd = 2 )
    }
    legText <- c (  "Aggregate", "Portland Inlet" )
    panLegend ( x = 0.2, y = 0.2, legTxt = legText, col = habColours, pch = 2, 
                bty = "0" )
    ani.pause ( )
  }
}, movie.name = "Habitats.gif", interval = 1, nmax = 500, ani.width = 1000,
ani.height = 800 )
