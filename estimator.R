#---------------------------------------------------------------------
# 
# REM612 Project Script - estimation script
# 
# Author: Samuel Johnson
# Date: 19/1/15
# 
# The following script takes the output of the fish and boat simulators
# and estimates habitat range and movement speed. The habitat range is
# estimated based on the aggregate and each inlet's worth of fish
# for different sample sizes of the whole, and movement speed on those
# fish in each sample which have multiple detections.
# 
# ToDo: ~~ADEhabitat range~~
#       ~~sample sizes~~
#       bootstrapping
# 
#---------------------------------------------------------------------

# Housekeeping
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
load ( file = detCtl $ inputSim )
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

# areas <- compareFunc ( N = 25, detections = yearDetections [[ 1 ]],
#                           indices = inletFishIdx, states = F )

# To bootstrap, we want to pull a sample of a given size a significant number
# number of times. So, let's choose sample sizes based on expected size of
# the overlap - didn't work right now, get back to that later

# Choose sample sizes - as size increases independence drops
sampSizes <- c ( 10, 25, 50, 100, 250 )

# Choose number of samples
nSamples <- 25

# The following will perform the comparison of habitats for a given number
# of samples of a given size. It works in parallel to speed things up.
ptm <- proc.time ( )
clust <- makeCluster ( spec = min ( nCores, length ( sampSizes ) ) )
clusterExport ( cl = clust, varlist = c ( "simCtl", "detCtl" ) )

# Call wrapper functions through lapply. This will return a list with an entry
# for each sample size, each a list of nSamples results.
sampList <- parLapply ( cl = clust,  X = sampSizes, fun = sampSizeFun, 
                        nSamp = nSamples, detDT = yearDetections [[ 1 ]],
                       idx = inletFishIdx, states = F )

# Stop cluster
stopCluster ( clust )

# Comparison timing
compTime <- proc.time () - ptm
cat ( " The comparison process took: ", compTime )

# Write output to file
save ( sampList, file = "sample.RData" )


# Now compute the mcp of habitat from detections and simuations for each sample
# size
polys <- lapply ( X = sampSizes, FUN = polyFunc, 
                  detections = yearDetections [[ 1 ]],
                  indices = inletFishIdx, states = F )

save ( polys, file = "polyData.RData" )

# And get the polys and area for all fish
# Pull out detections
detections <- yearDetections [[ 1 ]]

# Make detections polygon
allPolyDet <- mcp ( xy = detections [ , c("boat.lon", "boat.lat" ), 
                                      with = FALSE ],
                    id = detections [ , agg ], percent = 100 )

# Create a data.table of the true locations
trueDT <- lapply ( X = 1:10000, FUN = trueLocs, nT = simCtl $ nT, 
                       states = F )

trueDT <- rbindlist ( trueDT )
trueDT [ , agg := 1 ]

# Make polygon for true locations
allPolyTrue <- mcp ( xy = detections [ , c("boat.lon", "boat.lat" ), 
                                      with = FALSE ],
                    id = detections [ , agg ], percent = 100 )

save ( allPolyDet, allPolyTrue, file = "allPoly.RData" )