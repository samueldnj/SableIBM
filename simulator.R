#---------------------------------------------------------------------
# 
# REM612 Project Script - Simulator for fisher (like frogger)
# 
# Author: Samuel Johnson
# Date: 19/1/15
# 
# The following is a script which will run dual simulators. The first
# is a sablefish migration from nearshore inlets to offshore fishing 
# grounds using a SDE approach. The second is a DEVS of trawl boats
# fishing across their path. State variables for boats and fish are
# generated containing location and time information.
# 
# ToDo: ~~Functionalise the DEVS for the boats~~
#       ~~Parallelise DEVS and Fish~~
#       ~~Put parameters into ctl file~~
#       Add realism
#         ~~SDEs for fish movement~~ Still some bugs to iron out
#         Multiple movement modes - potential functions for home ranges
#         Fix random jump distance to logNormal dist
#         Mortality
#         Exploitation rate <- better with multiple fish
#         Fishing behaviour based on caught fish.
#         Schooling behaviour
#       Software agnostic output
# 
#---------------------------------------------------------------------

# Housekeeping
rm ( list = ls( ) )

# load packages
require ( parallel ); require ( foreach ); require ( dplyr ); 
require ( PBSmapping ); require ( doParallel ); require ( data.table );
require ( stringr )

#### Load bathymetry #### 
# if no RData file is available, source bathymetry.R - time expensive!
# source ( "bathymetry.R" )
load ( file = "bathyTop.RData" )
# If histAveByYear.RData is not available, uncomment and run histAnalysis.R
# to compute historical yearly averages for trip lengths and frequency.
# source ( "histAnalysis.R" )
load ( file = "histAveByYear.RData" )

# Read in historical fishing locations, for more realistic pressure.
trawlHist <- read.table ( "TrawlFisheries_not4B.txt", header = TRUE )
trawlHist <- as.data.table ( trawlHist )


# Source functions from separate script
source ( "functions.R" )

# Read in control list for the simulator
simCtl <- lisread ( "simulator.ctl" )

# Read number of cores from control file, or detect
# nCores <- ctlList $ nCores
nCores <- detectCores ( ) - 1

# Set seed for testing - uncomment to activate
set.seed ( simCtl $ ranSeed )

### Pseudo-code for DEVS ###
## 1. Initialise
# 1a. Set cumulative stats and clock to zero
nDetections <- 0
t <- 0

# 1b. Set number of fish and length of time from control file. Start with
# 365 days = 12410 hours, this may be tweaked in future.
nF <- simCtl $ nF 
nT <- simCtl $ nT 

# Set standard deviation and mean for fish jumps distance parameter
distSD <- simCtl $ distSD
distMu <- simCtl $ distMu

# Set standard deviation for fish jump direction parameter
sigmaTheta <- simCtl $ sigmaTheta

# A point which attracts fish out of deeper inlet waters, turns off early in
# migration
potentialPoint <- simCtl $ offshorePotential

# A matrix of tagging locations in inlets.
inlets <- simCtl $ inlets

# Fishing bounds - basically national boundaries to stop boats from straying
# outside canadian waters - will also help with plotting.
fBoundNorth <- simCtl $ fBoundNorth
fBoundEast <- simCtl $ fBoundEast


## 1c. Define state variables for boats, fish. Need to include:
#   Location as a function of time (lat, long - or grid cell reference)
#   (boats) Docked (0) or Fishing (1) TripLength (integer - days)
#   (fish) At large (0) or Caught (1) Anything else?
# Discrete events are the fishing events, so no need to include that
# info. Time spent in state will be assigned to each event (docked/towDuration). 
# Caught state is terminal for fish.

# Array to hold state variables for fish - need location and status at each
# time period for each fish.

# First, names for the dimensions of the array
Fnames <- list ( )
Fnames [[ 1 ]] <- 1:nF
Fnames [[ 2 ]] <- c ( "latitude", "longitude", "status" )
Fnames [[ 3 ]] <- 1:nT

# Then create the array
F <- array (  0, dim = c ( nF, 3, nT ), dimnames = Fnames )

# Named variables for fish and boat status
# Boats - docked or at sea
docked <- 0
fishing <- 1

# Fish statuses - no fishing mortality yet
atLarge <- 0
caught <- 1

### 2. Generate fish paths. 
# Using an SDE approach, we run fish downhill (gradient preference)
# with process error, see fishJumps() for details. We'll also time the process
# and report on
cat ( "Generating states for ", nF, " fish over ", nT, 
      " time steps (hours). \n", sep = "")
ptm <- proc.time ( )
F [ , 1:2, ] <- fishPaths ( nFish = nF, nT = nT, sigmaD = distSD, 
                            meanD = distMu, sigmaT = sigmaTheta,
                            bathy = bathyGrad )
fishTime <- proc.time () - ptm
cat ("State generation for fish took \n")
print ( fishTime )

## 3. Run boat DEVS in parallel - see functions.R for details
# GPU? Not yet.
# Increase recursive expression count
options ( expressions = simCtl $ recLimit )

# vector list to hold years of fishing
fishingYears <- vector ( mode = "list", length = ( nrow ( FishByYear ) ) )

# Create a list of variable names for export to the parallel cluster
parList <- list ( "fishing", "docked", "bathyTop", "nT" )

# Loop to populate fishingYears using function yearDEVS, which runs the DEVS
# for the historic averages. We'll also time the process.
ptm <- proc.time ( )
years <- FishByYear [ , year ]
for ( i in 1:length ( years ) )
{
  # Tell us what's going on
  cat ( "Generating fishing effort for ", years [ i ], 
        "\n", sep = "" )

  # Fill vector list with state vectors and FELs
  fishingYears [[ i ]] <- yearDEVS ( dt = FishByYear [ year == years [ i ], ], 
                                     T = nT )  
}
boatTime <- proc.time ( ) - ptm

# Check timing for the boat DEVS
cat ( "State generation for boats took \n " )
print ( boatTime )

# save output to a file for estimation procedures
# Make unique identifier for system time
splitTime <- str_split ( string = Sys.time ( ), pattern = " " )[[1]]
timeStampChar <- paste ( splitTime[1], splitTime[2], sep = "-" )
# filename <- paste ( "Output/simOut-", timeStampChar, ".RData", sep = "" )
filename <- "Output/simOut.RData"
save ( F, fishingYears, file = filename )





