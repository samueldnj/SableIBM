#---------------------------------------------------------------------
# 
# REM612 Project Script - Bathymetry loading
# 
# Author: Samuel Johnson
# Date: 19/1/15
# 
# This script will load bathymetry data and turn it into a topographic
# matrix of altitiudes
# 
#---------------------------------------------------------------------

# Housekeeping
rm ( list = ls ( ) )

# require packages
require ( PBSmapping ); require ( OpenCL ); require ( RColorBrewer )

# load bathymetry data
bathyData <- read.table ( "bathymetry.txt" )
bathyData $ V1 <- bathyData $ V1 - 360

# make a topography out of it
bathyTop <- makeTopography ( bathyData )

# Colours for plotting
greenPal <- rev ( brewer.pal ( 9, "Greens" ) )
bluePal <- rev ( brewer.pal ( 9, "Blues" ) [ 4:9 ] )

# Ramp up number of colours
greenPal <- c ( colorRampPalette ( greenPal ) ( 3000 ),
                colorRampPalette ( greenPal [ 9 ] ) ( 1185 ) )
bluePal <- c ( colorRampPalette ( bluePal [ 1:5 ] ) ( 4868 ),
                colorRampPalette ( bluePal [ 6 ] ) ( 200 ) )

# Now combine into one palette
bgPal <- c ( bluePal, greenPal )

# We should pre-compute the gradients at each grid point to save computations
# later - 4 per time step per fish!!! This could be a job for openCL!
xVec <- bathyTop $ x
yVec <- bathyTop $ y
zMat <- bathyTop $ z

# Create storage for gradient components
dXMat <- matrix ( 0, nrow = length ( xVec ), ncol = length ( yVec ) )
dYMat <- matrix ( 0, nrow = length ( xVec ), ncol = length ( yVec ) )

# Loop over x and y coordinates and compute finite differences for each cell
for ( xEntry in 1 : length ( xVec ) )
{
  for ( yEntry in 1 : length( yVec ) )
  {
    if ( xEntry == 1 )
    {
      # Forward difference in the far west
      dZx <- zMat [ xEntry + 1, yEntry ] - zMat [ xEntry, yEntry ]
      dx <- ( xVec [ xEntry + 1 ] - xVec [ xEntry ] )
    }
    # Backward difference in the far east
    if ( xEntry == length ( xVec ) )
    {
      dZx <- zMat [ xEntry, yEntry ] - zMat [ xEntry - 1, yEntry ]
      dx <- ( xVec [ xEntry ] - xVec [ xEntry - 1 ] )
    }
    # and central difference elsewhere
    if ( xEntry > 1 & xEntry < length ( xVec ) )
    {
      dZx <- zMat [ xEntry + 1, yEntry ] - zMat [ xEntry - 1, yEntry ]
      dx <- ( xVec [ xEntry + 1 ] - xVec [ xEntry - 1 ] )
    }
    # Find delta x, convert to deg/deg
    dXMat [ xEntry, yEntry ] <- dZx / dx / 111000

    # Do the same for y - forward difference in the south
    if ( yEntry == 1 )
    {
      dZy <- zMat [ xEntry, yEntry + 1 ] - zMat [ xEntry, yEntry ]
      dy <- ( yVec [ yEntry + 1 ] - yVec [ yEntry ] )
    }
    # Backward difference in the far east
    if ( yEntry == length ( yVec ) )
    {
      dZy <- zMat [ xEntry, yEntry ] - zMat [ xEntry, yEntry - 1 ]
      dy <- ( yVec [ yEntry ] - yVec [ yEntry - 1 ] )
    }
    # and central difference elsewhere
    if ( yEntry > 1 & yEntry < length ( yVec ) )
    {
      dZy <- zMat [ xEntry, yEntry + 1 ] - zMat [ xEntry, yEntry - 1 ]
      dy <- ( yVec [ yEntry + 1 ] - yVec [ yEntry - 1 ] )
    }
    # Find delta y, convert to deg / deg
    dYMat [ xEntry, yEntry ] <- dZy / dy / 111000
  }
}

bathyGrad <- list ()
bathyGrad $ x <- xVec
bathyGrad $ y <- yVec
bathyGrad $ dY <- dYMat
bathyGrad $ dX <- dXMat

# Also want to compute largest euclidean radius of water (negative alt)
# in order to have a lookup table of variances for process error
zSign <- apply ( X = zMat, MARGIN = 2, FUN = sign )

# Save output to speed up computations later
save ( bathyTop, bathyGrad, bgPal, file = "bathyTop.RData" )

