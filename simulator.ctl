## --------------------------------------------------------------------------
## simulator.ctl -
##
## A lisread readable control file for the simulator in JohnsonProject.R
## which contains control parameters for the fish and boat simulators.
##
## Author: Samuel Johnson
## Date: April, 2015
## --------------------------------------------------------------------------

###############################
## General parameters #########
###############################
## Number of time steps
# nT
20

## Number of cores for parallelism
# nCores
7

## Random number seed
# ranSeed
612

## Recursive expressions limit for boatMoves()
# recLimit
500000

###############################
## Fish simulator parameters ##
###############################

## Number of fish to simulate
# nF
10

## Random fish jump parameters ##
## Jump distance
# distSD
0.0108108108

# distMu
0.0036036036

## Jump heading
# sigmaTheta
1.047198

## Starting locations ##
# inlets
-130.430632 54.717845
-129.281188 53.060577
-128.494292 52.403640
-127.929869 51.919249

## Potential points for movement generation ##
## Potential to overcome inlet depth
# offshorePotential
-134 52


###############################
## Boat simulator parameters ##
###############################

## Boundaries to fishing area - canadian border ##
# fBoundNorth
55

# fBoundEast
-125

