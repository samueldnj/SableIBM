#---------------------------------------------------------------------
# 
# REM612 Project Script - Detection processing script
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
# ToDo: Sampling of population for bootstrapping
#       Break up analysis into inlets
# 
#---------------------------------------------------------------------

# housekeeping - if necessary
rm ( list = ls ( ) )

# load packages
require ( parallel ); require ( foreach ); require ( dplyr ); 
require ( PBSmapping ); require ( doParallel ); require ( data.table );
require ( stringr )

# load functions - mseRtools are included
source ( "functions.R" )

# Read in control parameters
detCtl <- lisread ( "detections.ctl" )
simCtl <- lisread ( "simulator.ctl" )

# Load output of simulator.R
load ( file = detCtl $ inputFish )
load ( file = detCtl $ inputBoats )

# Recover number of time steps, number of fish
nT <- simCtl $ nT
nF <- simCtl $ nF

# Create a vector list to hold the detections for each year of data
yearDetections <- vector ( mode = "list", length = length ( fishingYears ) )

# Run loop to populate list
for ( i in 1:length ( fishingYears ) )
{
  # Call detections function on each year in fishingYears
  yearDetections [[ i ]] <- detectFish (  boatDat = fishingYears [[ i ]],
                                          fArray = F, 
                                          radius = 0.005 )
}

# Save detections to output file.
# Make unique identifier from system time
splitTime <- str_split ( string = Sys.time ( ), pattern = " " )[[1]]
timeStampChar <- paste ( splitTime[1], splitTime[2], sep = "-" )
# filename <- paste ( "Output/detections-", timeStampChar, ".RData", sep = "" )
filename <- "Output/detOut.RData"
save ( yearDetections, file = filename )