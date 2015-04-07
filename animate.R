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
require ( animation ); require ( data.table );

# Load functions
source ( "functions.R" )

# Load detections control file for RData file paths
detCtl <- lisread ( "detections.ctl" )
simCtl <- lisread ( "simulator.ctl" )

# Uncomment to load data for animation
load ( file = detCtl $ inputFish )
load ( file = detCtl $ inputBoats )
load ( file = detCtl $ inputDetections )

# Let's try to save a short movie of the fish moving around
saveGIF ( 
{
  par ( bg = NA )
  for ( i in 1:100 ) 
  { plot (  x = F[ 1:1000, 1, ( i - 1 ) * 20 + 1 ],
            y = F[ 1:1000, 2, ( i - 1 ) * 20 + 1 ],
            xlim = c ( -134, -124 ),
            ylim = c ( 48, 56 ),
            las = 1,
            cex = 0.01,
            col = "red" )
    ani.pause ( )
  }
}, movie.name = "FishOverlay.gif", interval = 0.1, nmax = 500, ani.width = 600,
ani.height = 480 )


