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
require ( parallel )

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
# spawnPlot <- "bathyMapSpawning.png"
# png ( filename = spawnPlot, width = 1000, height = 800 )
# image ( bathyTop, col = bgPal )
# polygon ( x = c ( -134, -135, -127, -126 ), y = c ( 55, 52, 47, 48 ), 
#           col = "darkred", density = 25 )
# dev.off ( )

# Larval Drift
