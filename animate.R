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
load ( file = detCtl $ inputFish )
load ( file = detCtl $ inputBoats )
load ( file = detCtl $ inputDetections )

B <- fishingYears [[ 1 ]] $ B
nB <- dim ( B ) [ 3 ]

# Let's try to save a short movie of the fish moving around
saveGIF ( 
{
  for ( i in 1:500 ) 
  { image ( bathyTop, col = bgPal) 
    points (  x = F [ 1:1000, 1, ( i - 1 ) * 20 + 1 ],
              y = F [ 1:1000, 2, ( i - 1 ) * 20 + 1 ],
              cex = 0.1,
              col = "tomato" )
    points ( x = B [ 1, ( i - 1 ) * 20 + 1, 1:nB ],
             y = B [ 2, ( i - 1 ) * 20 + 1, 1:nB ],
             cex = 0.5, col = "whitesmoke" )
    ani.pause ( )
  }
}, movie.name = "FishMove.gif", interval = 0.1, nmax = 500, ani.width = 1200,
ani.height = 960 )


