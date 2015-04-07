#--------------------`-------------------------------------------------
# 
# REM612 Project Script - parallel functions
# 
# Author: Samuel Johnson
# Date: 19/1/15
# 
# This script contains the functions used inside parallel processes
# for the fishing script.
# 
#---------------------------------------------------------------------

require ( data.table )

# Function to look up bathymetric grid references for a given point.
#   Args:   lat = latitude of point
#          long = longitude of point
#         bathy = list of bathymetric data, produced by bathymetry.R
# 
#   Rets:   out = 2-ple of grid cell references for given point

bathyLookup <- function ( lat = 0, long = 0, bathy = bathyTop )
{
  # Extract vectors of grid line refs and matrix of bathymetry
  xVec <- bathy $ x
  yVec <- bathy $ y

  # Now compare Xt and Yt to the vectors of grid lines refs
  xDiff <- xVec - lat
  yDiff <- yVec - long

  # Find minimum difference
  xMin <- min ( abs ( xDiff ) )
  yMin <- min ( abs ( yDiff ) )

  # And which entry has that diff - if there are two equidistant
  # then it takes the maximum index. Unsophisticated, but it works.
  # In future, perhaps we can shift the finite difference.
  xEntry <- max ( which ( abs ( xDiff ) == xMin ) )
  yEntry <- max ( which ( abs ( yDiff ) == yMin ) )
 
  # concatenate the grid cell entries and return
  out <- c ( xEntry, yEntry )
  return ( out )
}


# This function converts a cartesian ordered pair into a polar ordered pair
#   Args: x = x-value
#         y = y-value
#   Rets:
#     (r,t) = (radius, angle)
cart2pol <- function ( x = 1, y = 1 )
{
  # Convert to a complex number
  z <- complex ( real = x, imaginary = y )

  # Extract radius and angle
  r <- Mod ( z )
  t <- Arg ( z )

  # return
  return ( c ( r, t ) )
}

# This function converts a polar ordered pair into a cartesian ordered pair
#   Args: r = radius
#         t = angle
#   Rets:
#     (x,y) = cartesian coords of point
pol2cart <- function ( r = 1, t = 0 )
{
  # convert to complex number
  z <- complex ( modulus = r, argument = t )

  # Extract real and imaginary parts
  x <- Re ( z )
  y <- Im ( z )

  # Return ordered pair
  return ( c ( x, y ) )
}

# Function to choose next location for fishing, constraints will be based
# on current boat status and whether or not it is trawling or moving. The
# function relies on recursion to make sure new locations are within
# constraints on fishing areas.
#   Args: status = fishing or docked
#            lat = latitude of current position
#           long = longitude of current position
#         towing = boolean for selecting movement speed distribution
#          eTime = amount of time movement is happening for
#          bathy = bathymetry data as produced by bathymetry.R
#   Rets: newLoc = a vector of latitude and longitude
boatMove <- function ( status = docked, lat = 0, long = 0, towing = FALSE,
                       eTime = 1, bathy = bathyTop, histData = histData,
                       recursive = FALSE )
{
  # Recover historic data from entry list.
  minDepth <- histData $ minDepth 
  maxDepth <- histData $ maxDepth 
  
  # First, if boat is docked, generate a random starting location
  if ( status == docked | ( lat == 0 & long == 0) )
  { 
    # Uniform random draw from the fishing area - very crudely drawn from
    # within a polygon around the 8 stat areas
    # Choose new latitude
    newLat <- runif ( n = 1, min = - 134, max = - 124 )
    
    # Choose new longitude within polygon
    yUB <- ( - 6 * newLat - 407 ) / 7
    yLB <- -2. * ( newLat + 125 ) / 3. + 47
    newLong <- runif ( n = 1, min = yLB, max = yUB )
  } else {
    # Draw a random direction
    tJump <- rnorm ( n = 1, mean = 0, sd = pi / 2 )

    # Draw a random motoring speed in knots
    {
      if ( towing ) speed <- rlnorm ( n = 1, mean = log ( 3 ), sd = 0.3 )
      else speed <- rlnorm ( n = 1, mean = log ( 6 ), sd = 0.5 )
    }

    # Convert speed to degrees per hour
    degSpeed <- speed * 1.852 / 111

    # Create jump distance
    rJump <- (eTime - 1) * degSpeed

    # Convert old location to a complex number for adding jump
    oldLoc <- complex ( real = lat, imaginary = long )

    # Convert jump to complex number
    jump <- complex ( modulus = rJump, argument = tJump )

    # Make new location
    newLoc <- oldLoc + jump

    # Recover new latitude and longitude
    newLat <- Re ( newLoc )
    newLong <- Im ( newLoc )
  }

  # Look up depth of new location in bathymetry table. First get gridref
  gridCell <- bathyLookup ( lat = newLat, long = newLong, bathy = bathy )
  xEntry <- gridCell [ 1 ]
  yEntry <- gridCell [ 2 ]

  # Recover depth from bathy matrix
  depth <- bathy $ z [ xEntry, yEntry ]

  # Combine new position into a vector for returning
  outVec <- c ( newLat, newLong )

  if ( recursive )
  {
    # Recurse function to correct for bad points
    while ( depth < minDepth | 
            depth > maxDepth | 
            outVec [ 2 ] > 55. |
            outVec [ 1 ] > -125. | 
            outVec [ 2 ] > ( - 6. * outVec [ 1 ]  - 407.) / 7. |
            outVec [ 2 ] < -2. * ( outVec [ 1 ] + 125 ) / 3. + 47. )
    {
      # Generate a new location
      outVec <- boatMove ( status = status, lat = lat, long = long, 
                           towing = towing, eTime = eTime, bathy = bathy,
                           histData = histData, recursive = FALSE )

      # Look up depth of new location in bathymetry table. First get gridref
      gridCell <- bathyLookup ( lat = outVec [ 1 ], long = outVec [ 2 ], 
                                bathy = bathy )
      xEntry <- gridCell [ 1 ]
      yEntry <- gridCell [ 2 ]

      # Recover depth from bathy matrix
      depth <- bathy $ z [ xEntry, yEntry ]
    }
  }
  # Return location
  return ( outVec )
}


# The following function is the event scheduler, which will take the event number
# and current state of the system and schedule the next event.
#   Arguments:   time = current time
#               event = current event number
#               state = slice of B for boat in question
#                 fel = future event list for the boat in question
# 
#   Returns:      fel = future event list with new event scheduled

eventSched <- function ( time = t, event = eCount,
                        fel = FEL, ntows = nTows, histData = histData )
{ 
  # Recover boat number, trip threshold and boat number from arguments
  boat <- fel [ event, ] $ boat.number
  thresh <- fel [ event, ] $ trip.thresh
  interFish <- histData $ interFish
  meanEventTime <- histData $ meanEventTime
  layOver <- histData $ layOver
  meanTows <- histData $ meanTows

  
  if ( time != fel [ event, ] $ time.stamp ) 
    {
      cat ( "ACK! Failure on time variable - review logic. \n
             boat = ", boat, "\n
             t = ", time, "\n
             time.stamp = ", fel [ event, ] $ time.stamp ); return ( fel )
    }

  # If boat is docked, set next fishing event
  if ( fel [ event, ] $ event.type == docked )
  {
    # Draw layover time from Poisson distribution
    tLay <- fel [ event, ] $ event.time

    # Set next time
    tNext <- min ( nT, time + tLay )

    # Create new threshold for nTows, reset nTows
    tripThresh <- rpois ( n = 1, lambda = meanTows )
    ntows <- 0

    # Create towDuration for next event
    towTime <- towDuration ( lam = meanEventTime )

    # Add event to FEL
    fel [ event + 1, ] <- c ( boat, tNext, fishing, towTime, ntows, tripThresh )
  }
  
  # If boat is fishing, check threshold variable, if below threshold then
  # set new fishing event, else set docking event.
  if ( fel [ event, ] $ event.type == fishing )
  {
    # recover tow duration for next event time
    towDur <- fel [ event, ] $ event.time

    # Set next event time by adding towDuration 
    # to current time
    tNext <- time + towDur - 1

    #   If tripLength >= thresh, dock, else fish again
    if ( ntows >= thresh )
    {
      # Generate layover time
      tLay <- layTime ( lam = layOver )

      # Add to FEL
      fel [ event + 1, ] <- c ( boat, tNext, docked, tLay, 0, 0 )
    } else {
      # Generate interfishing time
      interfishTime <- rpois ( n = 1, lambda = interFish )

      # Add to next event time
      tNext <- tNext + interfishTime

      # Generate set duration for next event
      towDurNext <- towDuration ( lam = meanEventTime )

      # Add next event to FEL and
      # Correct for going over time - dock if basically at the end of the year
      if ( tNext + towDurNext > nT )
      {
        # Set tNext to time limit
        tNext <- nT

        # update FEL to reflect end of year behaviour
        fel [ event + 1, ] <- c ( boat, tNext, docked, 0, 0, 0 )
      } else {
        fel [ event + 1, ] <- c ( boat, tNext, fishing, towDurNext, ntows + 1, 
                                  thresh )
      }
    }
  }

  # return appended future event list
  return ( fel )
}

# The following function processes discrete events according to the future
# event list (FEL) for each boat and updates the state variables to reflect
# the behaviour of the boat.
#   Arguments:    event = the counting index of the FEL entry being processed
#                   fel = the future event list
#                 state = the slice of the state variable array being updated
# 
#     Returns:    state = the updated slice of the systems state array
eventProc <- function ( event = eCount, fel = FEL, nowlat = nowLat, 
                        nowlong = nowLong, ntows = nTows, histData = histData,
                        b = B )
{ 
  # Get time stamp from FEL for logic check
  time <- fel [ event, ] $ time.stamp

  if ( time != fel [ event, ] $ time.stamp ) 
    {
      cat ( "ACK! Failure on time variable - review logic. \n
             boat = ", boat, "\n
             t = ", t, "\n" ); return ( fel )
    }

  # Recover other parameters for system state
  eTime <- fel [ event, ] $ event.time
  boat <- fel [ event, ] $ boat.number
  status <- fel [ event, ] $ event.type

  # Check event type and update state
  # First, if the boat is docked
  if ( fel [ event, ] $ event.type == docked )
  {
    # Update state: change status to docked, reset tripLength counter and set
    # location to zero
    state <- 0 # auto recycles to set state back to docked
  }
  # Next if the boat is fishing
  if ( fel [ event, ] $ event.type == fishing )
  {
    # generate event starting location
    startLoc <- boatMove ( status = status, lat = nowlat, long = nowlong,
                           towing = FALSE, eTime = interfishTime, 
                           bathy = bathyTop, histData = histData, 
                           recursive = TRUE )

    # generate event ending location
    endLoc <- boatMove (  status = status, lat = startLoc [ 1 ], 
                          long = startLoc [ 2 ],
                          towing = TRUE, eTime = eTime, bathy = bathyTop,
                          histData = histData, recursive = TRUE )

    # generate jump sizes
    if ( eTime == 1 ) jumps <- 0
    else jumps <- ( endLoc - startLoc ) / ( eTime - 1 )

    # Recover system state
    state <- b [ , time:( time + eTime - 1 ) ]

    # Fix this for small tows - not pretty
    if ( class ( state ) != "matrix" )
    {
      state <- as.matrix ( state )
    }
    # Set location - this should be made more sophisticated in future, even
    # use a Markovian approach (for tripLength > 1)
    for ( tStep in 0:( eTime - 1 ) )
    {
      state [ 1:2, tStep + 1 ] <- startLoc + tStep * jumps
    }

    # Set status to fishing
    state [ 3, ] <- fishing

    # Increment nTows and update
    state [ 4, ] <- ntows

    # Store ending location for next event
    nowLat <<- endLoc [ 1 ]
    nowLong <<- endLoc [ 2 ]
  }

  # Return state of boat during fishing event
  return ( state )
}

# Function to randomly produce a tow duration for each fishing event
#   Args: lam = parameter for Poisson distribution drawn from
#   Rets: towTime = duration of trawl event
towDuration <- function ( lam = 3 )
{
  towTime <- 0
    while ( towTime == 0 )
    {
      towTime <- rpois ( n = 1, lambda = lam)  
    }
  return ( towTime )  
}

# Function to randomly produce layover times for fishing vessels
#   Args: lam = parameter for Poisson distribution drawn from
#   Rets: lay = duration of layover
layTime <- function ( lam = layOver )
{
  lay <- rpois ( n = 1, lambda = lam )
  return ( lay )
}

# Function to count record the row numbers for non-zero
# rows in a matrix
#   Arguments:  mat = matrix to be searched
#   Returns:   ridx = vector of row indices
rowIdx <- function ( mat = boatFishingTimes )
{
  # Sum the entries in each row
  rowSums <- apply ( X = mat, MARGIN = 1, FUN = sum )

  # Get the row numbers and return
  ridx <- which ( rowSums > 0 )
  return ( ridx )
}

# Function to check distance of fish from a given centre point location,
# usually of fishing gear. Distance is computed as euclidean distance,
# as radius is usually small enough for this to be fine.
#   Arguments:     b = column number of cLoc
#               cLoc = lat/long of centre point if circle
#               fLoc = matrix of fish locations, with lat/long as rows and
#                      boats as columns
#                  t = time step for which comparison is being made
#                rad = radius of euclidean ball
# 
#   Returns: detections = a data.table with fish number, location and time
bDetections <- function ( b = 1, cLoc = boatLocs, fLoc = fishLocs, t = t,
                          rad = 0.005 )
{
  # Recover location of centrepoint from cLoc
  cLoc <- cLoc [ , b ]

  # Create matrix to hold distances from boats to fish
  diff <- matrix ( 0, nrow = nrow ( fLoc ), ncol = 2 )

  # Compare to fish locations
  diff [ , 1 ] <- fLoc [ , 1 ] - cLoc [ 1 ]
  diff [ , 2 ] <- fLoc [ , 2 ] - cLoc [ 2 ]

  # square and sum for euclidean distance
  diff2 [ , 1 ] <- diff [ , 1 ] * diff [ , 1 ]
  diff2 [ , 2 ] <- diff [ , 2 ] * diff [ , 2 ]
  diffDist <- diff2 [ ,1 ] + diff2 [ , 2 ]

  # Look for fish with diffDist < radius
  closeFish <- which ( diffDist < rad )

  # Get number of detections
  nDetections <- length ( closeFish )

  # Set up data.table to hold detections, including fish number, location
  # and time
  detections <- data.frame (  fish.number = closeFish,
                              boat.lat    = rep.int ( x = cLoc [ 1 ], 
                                                      times = nDetections ),
                              boat.lon    = rep.int ( x = cLoc [ 2 ], 
                                                      times = nDetections ),
                              time.stamp  = rep.int ( x = t, 
                                                      times = nDetections ) )

  # convert to data.table and return
  detections <- as.data.table ( detections )
  return ( detections )
}

# Function to look over a given time and record instances of boats
# and fish being within a given radius. Calls bDetect, which breaks the
# detections over boats fishing at the given time.
#   Arguments:    t = time step for detections
#            bArray = array of boat state variables
#            fArray = array of fish state variables
#            radius = radius of detection
#   Returns:
#         tDetections = a data.table of detections at time t with fish 
#                       number, location and time
tDetect <- function ( t = 1, bArray = B, fArray = fArray, 
                          radius = 0.005 )
{
  # Get boat numbers for fishing boats at each time
  boats <- which ( bArray [ 3, t, ] == 1 )

  # Make matrix of locations - rows are lat/long, cols are boats
  boatLocs <- as.matrix ( bArray [ 1:2, t, boats ] )
  
  # Get matrix of fish locations at that time
  fishLocs <- fArray [ , 1:2, t ]

  # Get number of boats fishing at that time
  nb <- length ( boats )

  # lapply bDetections to get a list of detection data.tables
  detList <- lapply ( X = 1:nb, FUN = bDetections, cLoc = boatLocs, 
                      fLoc = fishLocs, t = t, rad = radius )

  # Bind list of detection data.tables into one data.table and return
  tDetections <- rbindlist ( detList )
  return ( tDetections )
}

# Function to append inlet numbers to F array and detections data.tables
# in the output of the detections code. This will be done from the beginning
# next time, but for now we must do it like this.
#   Args: fArr = array of fish state variables
#         inletLocs = matrix of inlet tagging locations
#         detList = a list of data.tables of detections
#   Rets: outList = a list containing
#                      - list of inlet fish array indices for fArr, 
#                      - a list of appended detection data.tables
#                      - an appended version of F
appendInlets <- function ( fArr = F, inletLocs = inlets, 
                           detList = yearDetections )
{
  # Extract the indices of fish in each inlet
  # A list to hold the indices
  inletFishIdx <- vector ( mode = "list", length = 4 )
  # This loop populates the list by checking those
  # entries in F with a given starting position, and also
  # updates the array of states with the inlet of origin.
  for (i in 1:nrow ( inletLocs ) )
  {
    idx <- unique ( which ( fArr [ , 1:2, 1 ] == inletLocs [ i, ], 
                    arr.ind = T ) ) [ , 1 ]
    fArr [ idx, 3, ] <- i
    inletFishIdx [[ i ]] <- list ( idx, i )
  }

  # Now, fix the dimnames
  dimnames ( fArr ) [[ 2 ]] <- c ( "latitude", "longitude", "inlet" )

  # Now to update the detection data.tables
  # Loop over the entries in detList
  for ( i in 1:length ( detList ) )
  {
    # Recover detection dt from detList
    dt <- detList [[ i ]]

    # lapply the function to inletFishIdx
    dtJoinList <- lapply ( X = inletFishIdx, FUN = inletJoin, dt = dt )

    # bind the result together
    dt <- rbindlist ( dtJoinList )

    # Add an aggregate column
    dt [ , agg := 1 ]
    detList [[ i ]] <- dt
  }

  # Make a list for returning
  outList <- list ( )
  outList $ inletFishIdx <- inletFishIdx
  outList $ detList <- detList
  outList $ F <- fArr

  # return list of output
  return ( outList )
}

# Create a function to join the inlet numbers - path specific to the above 
# function
inletJoin <- function ( idx = inletFishIdx [[ 1 ]], dt = detList [[ 1 ]] )
{
  # Create a df of inlet numbers and fish numbers, coerce to dt
  idxDT <- data.frame ( fish.number = idx [[ 1 ]], 
                        inlet = rep.int ( x = idx [[ 2 ]], 
                                          times = length ( idx [[ 1 ]] ) ) )
  idxDT <- as.data.table ( idxDT )

  # set key of dt that you want to add inlet numbers to
  setkey ( dt, fish.number )

  # Join DTs and return the result
  dt <- dt [ idxDT, nomatch = 0 ]
  return ( dt )
}

# Function to build a data.table from fish states, for comparison
# of detected home range to actual over the course of the simulation
#   Arguments:
#          fN = fish number
#          nT = length of time to look over
#      states = array of fish states
#   Returns:
#        true = data.table of true locations of fish fN over time frame given
trueLocs <- function ( fN = 1, nT = simCtl $ nT, states = F )
{
  # Create a data.frame to hold info
  true <- data.frame ( fish.number = integer ( length = nT ),
                       fish.lat = double ( length = nT ),
                       fish.lon = double ( length = nT ),
                       time.stamp = 1:nT,
                       inlet = integer ( length = nT ) )

  # populate data frame
  # First the fish number
  true [ , 1 ] <- fN
  true [ , c( 2:3, 5 ) ] <- t ( states [ fN, 1:3, ] )

  # coerce to a data.table and return
  true <- as.data.table ( true )
  return ( true )
}


# Function to randomly sample fish and compare detections
# to true home range - for aggregate samples and inlet specific
# fish.
#   Arguments:         N = sample size
#             detections = detection data to use (based on historic fishing
#                          effort)
#                indices = indices of fish in each inlet - bad planning
#                 states = array of fish state variables
#   Returns:       areas = matrix of detected and true habitat areas, 
#                          and ratio between them.
compareFunc <- function ( N = 20, detections = yearDetections [[ 1 ]],
                          indices = inletFishIdx, states = F )
{
  # Count inlets
  nInlets <- length ( indices )

  # Take samples in each inlet
  # Make matrix to hold samples
  fishSamp <- matrix ( 0, nrow = nInlets, ncol = N )

  # Take samples
  for ( i in 1:nInlets )
  { 
    idx <- indices [[ i ]][[ 1 ]]
    fishSamp [ i, ] <- sample ( x = idx, size = N )
  }

  # Aggregate into a vector
  aggSamp <- as.vector ( fishSamp )

  # Reduce state array and detections table to sample
  sampF <- states [ aggSamp, , ]
  setkey ( detections, fish.number )
  sampDet <- detections [ J ( aggSamp ) ]

  # Need to create a data.table of locations for fish in aggSamp. Use lapply
  # on trueLocs to get one dt for each fish
  trueList <- lapply ( X = aggSamp, FUN = trueLocs, nT = simCtl $ nT, 
                       states = states )

  # Bind them together into one dt
  trueDT <- rbindlist ( trueList )
  trueDT [ , agg := 1]

  # And now we need to take the initial locations and add them to the detections
  # Recover tagging locations
  taggLocs <- trueDT [ time.stamp == 1 ]

  # Bind to sampDet
  sampDet <- rbind ( sampDet, taggLocs, use.names = FALSE )

  # Oki doke, now we can start comparing habitat areas. Polygons are for plotting
  # later, so we can keep these out of the analysis for now.
  # Find area of MCP around detections of aggregate sample
  mcpAggDetArea <- mcp.area ( xy = sampDet [ ,  c("boat.lat", "boat.lon"), 
                                                with = FALSE ], 
                     id = sampDet [ , agg ], plotit = FALSE, percent = 100 )
  # Find area of MCP of true locations of agg sample
  mcpAggTrueArea <- mcp.area (  xy = trueDT [ , c("fish.lat", "fish.lon"), 
                                                 with = FALSE ], 
                                id = trueDT [ , agg ], plotit = FALSE,
                                percent = 100 )

  areas <- matrix ( 0, nrow = nInlets + 1, ncol = 3 )

  dimnames ( areas ) [[ 1 ]] <- c (  "Aggregate", "Inlet1", "Inlet2", 
                                    "Inlet3", "Inlet4" )
  dimnames ( areas ) [[ 2 ]] <- c ( "Detections", "True", "Ratio" )

  areas [ 1, ] <- c ( as.numeric ( mcpAggDetArea ), as.numeric ( mcpAggTrueArea ),
                      as.numeric ( mcpAggDetArea ) /
                      as.numeric ( mcpAggTrueArea ) )

  
  # Function to lapply and speed up sample bootstrapping - not working right
  # now
  # inletMCP <- function ( inlet = 1, samples = sampDet, true = trueDT )
  # {
  #   # Reduce to fish from inlet i
  #   sampDetInlet <- samples [ inlet == inlet ]
  #   trueDTInlet <- true [ inlet == inlet ]

  #   # Find area of minimal convex polygon around detections
  #   mcpDetInletArea <- mcp.area ( xy = sampDetInlet [ ,  c("boat.lat", "boat.lon"),
  #                                                     with = FALSE ],
  #                                 id = sampDetInlet [ , agg ], 
  #                                 plotit = FALSE, percent = 100 )
  #   # Find area of MCP around true locations
  #   mcpTrueAreaInlet <- mcp.area (  xy = trueDTInlet [ , c( "fish.lat", 
  #                                                           "fish.lon" ),
  #                                                      with = FALSE ],
  #                                   id = trueDTInlet [ , agg ], plotit = FALSE,
  #                                   percent = 100 )


  #   areas <- c ( as.numeric ( mcpDetInletArea ), 
  #                as.numeric ( mcpTrueAreaInlet ),
  #                as.numeric ( mcpDetInletArea ) / 
  #                as.numeric ( mcpTrueAreaInlet ) )

  #   return ( areas )
  # }



  # out <-   lapply (   X = 1:nInlets, 
  #                     FUN = inletMCP,
  #                     samples = sampDet,
  #                     true = trueDT )
  for ( i in 1:nInlets )
  {
    # Reduce to fish from inlet i
    sampDetInlet <- sampDet [ inlet == i ]
    trueDTInlet <- trueDT [ inlet == i ]

    # Find area of minimal convex polygon around detections
    mcpDetInletArea <- mcp.area ( xy = sampDetInlet [ ,  c("boat.lat", "boat.lon"),
                                                      with = FALSE ],
                                  id = sampDetInlet [ , agg ], 
                                  plotit = FALSE, percent = 100 )
    # Find area of MCP around true locations
    mcpTrueAreaInlet <- mcp.area (  xy = trueDTInlet [ , c( "fish.lat", 
                                                            "fish.lon" ),
                                                       with = FALSE ],
                                    id = trueDTInlet [ , agg ], plotit = FALSE,
                                    percent = 100 )
    # Fill row of areas matrix
    areas [ i + 1, ] <- c ( as.numeric ( mcpDetInletArea ), 
                            as.numeric ( mcpTrueAreaInlet ),
                            as.numeric ( mcpDetInletArea ) / 
                            as.numeric ( mcpTrueAreaInlet ) )
  }

  return ( areas )
}

# # The following functions are wrappers for the compareFunc function,
# # in order to allow for swift parallelisation and returning
# # of lists in a decent structure for analysis.
# sampSizeFun <- function ( sampleSize = 25, nSamp = nSamples, 
#                           detDT = yearDetections [[ 1 ]],
#                           idx = inletFishIdx, states = F )
# {
#   outList <- lapply ( X = 1:nSamp, 
#                       FUN = repWrapper, 
#                       sampSize = sampleSize, det = detDT,
#                       idx = idx, stateVars = states )
#   return ( outList )
# }

repWrapper <- function ( times = 1, sampSize = 25, det = detDT, 
                         idx = idx, stateVars = states  )
{
  source ( "parFunctions.R" )
  n <- times
  areas <- compareFunc (  N = sampSize, detections = det, 
                          indices = idx, states = stateVars )

  return ( areas )
}


