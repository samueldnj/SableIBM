#---------------------------------------------------------------------
# 
# REM612 Project Script - functions
# 
# Author: Samuel Johnson
# Date: 19/1/15
# 
# This script contains the functions for the simulation model project
# 
#---------------------------------------------------------------------

# Source mseRtools for useful functions - all hail lisread
source ( "mseRtools.r" )

# The following function creates a transition matrix for a rectangular grid
# where the point will migrate towards the left edge of the grid on
# average and have non-zero probability of a stationary step.
#   Args: r = number of rows
#         c = number of columns
# 
#   Returns: westStatMat = transition matrix

westStatMat <- function ( r = rows, c = cols) 
{
  # The transition matrix is square with dimension r*c
  dim <- r*c

  # First create the matrix of zeroes
  transMat <- matrix ( data = 0, nrow = dim, ncol = dim )

  ## Then we need to populate the matrix of probabilities.
  ## First the corners of the grid, except the left corners

  # Top right
  transMat [  c, 
              c ( c - 1, 
                  2 * c - 1, 
                  c, 2 * c ) ] <- rep ( x = c(2./6., 1./6. ), 
                                        times = c( 2, 2 ) )

  # Lower right
  transMat [  dim, 
              c( dim, dim - c,
                 dim - 1, dim - 1 - c ) ] <- rep (  x = c(1./6., 2./6. ), 
                                                    times = c( 2, 2 ) )

  ## Now the edges of the grid
  # Top edge              
  for (i in 2:( c - 1) ){
    transMat [ i, c ( (i-1):(i+1), 
                      (i+c-1):(i+c+1) ) ] <- rep ( x = c ( 1./3., 
                                                           1./12., 
                                                           1./12. ),
                                                   times = 2 )
  }

  # Bottom edge
  for (i in (dim-c+2):(dim-1) ) {
    transMat [ i, 
               c( (i-1):(i+1), 
                  (i-c-1):(i-c+1) ) ] <- rep ( x = c ( 1./3., 
                                                       1./12.,
                                                       1./12. ),
                                               times = 2 )
  }

  # Left edge
  leftEdge <- seq ( from = 1, to = dim-c+1, by = c )
  for ( i in leftEdge ) {
    transMat [ i, i ] <- 1
  }

  # Right edge
  rightEdge <- seq (from = (2*c), to = (dim - c), by = c )
  for ( i in rightEdge ) {
    transMat [ i, c( i-c, 
                     i-c-1, 
                     i, i-1, i+c, i+c-1 ) ] <- rep (  x = c( 1./9., 2./9.), 
                                                                  times = 3 )
  }

  ## Now the interior cells
  for ( i in 2:(r-1) ) {
    for (j in 2:(c-1) ) {
      s <- (i - 1) * c + j 
      transMat[ s, c( (s-c-1):(s-c+1), 
                      (s-1):(s+1), 
                      (s+c-1):(s+c+1) ) ] <- rep (  x = c ( 2./9., 
                                                            1./18., 
                                                            1./18. ), 
                                                    times = 3 )
    }
  }

  return(transMat)
}


# This function creates a matrix of transition probabilities for a markov 
# chain on a rectangular grid with equal probability of moving to any nearest
# neighbour cell or staying put.
#   Args:     r = number of rows
#             c = number of columns
# 
#   Returns:  transMat = transition matrix

simpleStatMat <- function ( r = rows, c = cols ) 
{
  # The transition matrix is square with dim = r*c
  dim <- r*c

  # First create the matrix of zeroes
  transMat <- matrix ( data = 0, nrow = dim, ncol = dim )

  ## Then we need to populate the matrix of probabilities.
  ## First the corners of the grid
  # Top left
  transMat [ 1, c( 1, 2, c + 1, c + 2 ) ] <- rep ( x = 0.25, times = 4 )

  # Top right
  transMat [  c, 
              c ( c, 
                  c - 1, 
                  2 * c - 1, 2 * c - 2 ) ] <- rep ( x = 0.25, times = 4 )

  # Lower left
  # create a variable for the lower left corner
  llc <- dim - c + 1
  transMat [  llc , 
              c( llc, llc + 1, llc - c, llc - c + 1 ) ] <- rep (  x = 0.25, 
                                                                  times = 4 )

  # Lower right
  transMat [  dim, 
              c( dim, dim - 1, dim - 1 - c, dim - c ) ] <- rep (  x = 0.25, 
                                                                  times = 4 )

  ## Now the edges of the grid
  # Top edge              
  for (i in 2:( c - 1) ){
    transMat [ i, c ( (i-1):(i+1), 
                      (i+c-1):(i+c+1) ) ] <- rep ( x = 1./6., times = 6 )
  }

  # Bottom edge
  for (i in (dim-c+2):(dim-1) ) {
    transMat [ i, c( (i-1):(i+1), (i-c-1):(i-c+1) ) ] <- rep (  x = 1./6., 
                                                                times = 6 )
  }

  # Left edge
  leftEdge <- seq ( from = c+1, to = dim-2*c+1, by = c )
  for ( i in leftEdge ) {
    transMat [ i, c ( i-c, i-c+1, i, i+1, i+c, i+c+1 ) ] <- rep ( x = 1./6., 
                                                                  times = 6 )
  }

  # Right edge
  rightEdge <- seq (from = (2*c), to = (dim - c), by = c )
  for ( i in rightEdge ) {
    transMat [ i, c( i-c, i-c-1, i, i-1, i+c, i+c-1 ) ] <- rep (  x = 1./6., 
                                                                  times = 6 )
  }

  ## Now the interior cells
  for ( i in 2:(r-1) ) {
    for (j in 2:(c-1) ) {
      s <- (i - 1) * c + j 
      transMat[ s, c( (s-c-1):(s-c+1), 
                      (s-1):(s+1), 
                      (s+c-1):(s+c+1) ) ] <- rep (  x = 1./9., times = 9 )
    }
  }

  return(transMat)
}

# The following function will simulate a pair of parallel markov processes,
# a fish and a boat, moving around on a rectangular grid. The fish is 
# acoustically tagged, and the function returns the number of times the boat 
# succesfully detected the fish when they were in the same cell, with a given
# probablity of detection.
#   Args:     r = rows of the grid
#             c = columns of the grid
#             detProb = probability of detection
#             nT = number of time steps
#             print = boolean of whether to print the number of detections
# 
#   Returns:  d = number of detections
simpleSim <- function ( r = rows, c = cols, nBoats = 1,
                        detProb = p, nT = nT, print = TRUE ) 
{
  # Create transition matrix
  transitions <- simpleStatMat ( r = r, c = c )

  # make into a Markov Chain
  gridMC <- new ( "markovchain", 
                  transitionMatrix = transitions, name = "gridMC" )

  ## Simulate a boat and a fish as two parallel markov processes
  # The fishing boats
  # Create a matrix to hold fishing boats
  B <- matrix ( data = 0, nrow = nT, ncol = nBoats )

  # For loop to fill matrix
  for (boat in 1:nBoats) {
    B[ , boat ] <- as.numeric ( rmarkovchain ( n = nT, object = gridMC ) )
  }

  # The tagged fish
  F <- as.numeric ( rmarkovchain ( n = nT, object = gridMC ) )

  ## Now to use this to simulate detection using the detection efficiency
  # First, find the number of encounters between each boat and the fish
  opp <- length ( which ( B - F == 0 ) )

  # Then count the number of detections, using the probability of detection p
  d <- rbinom ( n = 1, size = opp, prob = p )

  # print the number of detections if asked
  if ( print ) {
   cat ( "The fish was detected by the boat", d, "times. \n", sep = " " )
  }

  # return the number of detections
  return ( d )

}

# This function converts body lengths per second into km/hr for a fish
# of a given length
#     Args:   blSpeed = ave. number of body lengths per second a fish swims
#             fLength = length of fish
#     Returns:  speed = ave speed in km/hr

speedConv <- function ( blSpeed = 0.6, fLength = 0.55 )
{
  # compute m/s
  msVelocity <- blSpeed * fLength

  # Convert to km/hr
  kmh <- msVelocity * 3.6

  # return
  return ( kmh )
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



# The function delta is for making a discrete time markov chain. This
# function takes a current location, grid dimensions and a transition matrix
# and returns a new location based.
#     Arguments:    loc = grid cell number of current location
#                     r = number of rows in grid
#                     c = number of columns in grid
#                 probs = transition matrix of probabilities
# 
#     Returns     delta = difference in grid refs of old and new locations
delta <- function ( loc = 1, r = r, c = c, probs = transMat ) 
{ 
  # Compute number of grid cells
  dim <- r*c
  if ( loc != 0 )
  {
    # delta ( location ) = sample ( grid cells, prob = transMat[location,] )
    newLoc <- sample ( x = 1:dim, size = 1, prob = probs [ loc, ] )
    
    # Find difference between new and old grid refs
    delta <- newLoc - loc
  } else { # loc = 0 means the fish hasn't entered grid yet.
    delta <- 0
  }

  # return difference
  return ( delta )
}

# The following function creates a matrix in which the rows are discrete
# time markov chains for fish movement through a grid.
#     Arguments:      r = rows of movement grid
#                     c = cols of movement grid
#                 nFish = number of fish to simulate
#                    nT = number of time steps
# 
#     Returns:        M = matrix with MCs of movement as rows

fishPathsOld <- function ( r = rows, c = cols, nFish = nF, nT = nT, 
                        entryCells = inlets ) 
{ 
  # Generate transition matrix
  transMat <- westStatMat ( r = r, c = c )

  # Create matrix to hold fish paths
  F <- matrix ( data = 0, nrow = nFish, ncol = nT )

  # Set random entry time for each fish
  entryTime <- sample ( x = 1:nT, replace = TRUE, size = nF )

  # Choose initial locations along the right hand edge (4 inlets) for each fish
  # Needs to be turned into an apply function - might be worth turning
  # F into a ragged array and using tapply.
  for ( f in 1:nF )
  {
    F [ f , entryTime [ f ] ] <- sample ( x = entryCells, 
                                          size = 1 )  
  }
  
  # run simulation forward for each fish, choosing next location
  # based on current location and the transition matrix - can't be
  # parallelised as not independent
  for ( t in 1:(nT - 1) )
  { 
    # The following is a little arcane. If the fish hasn't entered yet
    # F[ , t + 1] might be non-zero, but F[,t] and delta will be zero. If
    # the fish has entered, then F[,t] and delta are non-zero, but F[,t+1] 
    # is zero. Thus, while there are too many terms in the sum, it works out.
    F [ , t+1 ] <-  F [ , t + 1 ] + 
                    F [ , t ] +  
                    sapply (  X = F[ , t ], FUN = delta, r = r, c = c, 
                              probs = transMat )
  }

  # Return fish paths
  return ( F )
}

# The following function creates a matrix in which the rows are discrete
# time markov chains for fish movement through a grid.
#     Arguments:  nFish = number of fish to simulate
#                    nT = number of time steps
#             entryLocs = matrix with lat/long of initial positions
#                 sigma = sd of process error
#                   wgt = multiplicative parameter for potential function
#                 bathy = list of lat, long vectors and depth matrix
#     Returns:        P = 3D array of lat, long in time for each fish
fishPaths <- function ( nFish = nF, nT = nT, entryLocs = inlets, 
                        sigmaD = distSD, sigmaT = sigmaTheta, meanD = 0.01,
                        bathy = bathyGrad, potPoint = potentialPoint ) 
{ 
  # Create matrix to hold fish paths
  P <- array ( data = 0, dim = c ( nFish, 2, nT ) )

  # Choose initial locations along the right hand edge (4 inlets) for each fish
  # Needs to be turned into an apply function - might be worth turning
  # F into a ragged array and using tapply.
  initLoc <- sample ( x = 1:4, size = nF, replace = TRUE )
  P [ ,1:2,1] <- entryLocs [ initLoc, ]

  # run simulation forward for each fish, choosing next location
  # based on current location and the transition matrix - can't be
  # parallelised as not independent
  
  # Make a cluster
  clust <- makeCluster ( spec = nCores )


  for ( t in 1:( nT - 1 ) )
  { 
    
    # The following is a little arcane. If the fish hasn't entered yet
    # F[,,t + 1] might be non-zero, but F[,,t] and delta will be zero. If
    # the fish has entered, then F[,,t] and delta are non-zero, but F[,,t+1] 
    # is zero. Thus, while there are many terms in the sum, it works out.
    jump <- parApply (  cl = clust,  X = P[ , , t ], MARGIN = 1, 
                        FUN = fishJumps, sigmaD = sigmaD, sigmaT = sigmaT, 
                        meanD = meanD, bathy = bathy, potSource = potPoint )

    jump <- t ( jump )
    P [ , , t + 1 ] <-  P [ , , t + 1 ] + jump + P [ , , t ]
  }

  # stop Cluster
  stopCluster ( cl = clust )
    

  # Return fish paths
  return ( P )
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

# The following function will process an event from the FEL, computing the
# distance between the boat in the event and all fish at that time.
#   Arguments:   event = a row of the FEL
#                    B = the array of states for each boat
#                    F = the array of states for each fish
# 
#     Returns:    dist = a vector of length nF which contains the distances
boatFishDist <- function ( event = FEL [ 1, ], b = B, f = F )
{
  # Recover timestamp and event number
  time <- event [ 2 ]
  boat <- event [ 1 ]

  # Get location of boat at time from state array
  loc <- b [ boat, 1, time ]

  # Compare location of boat for each fish
  dist <- f [ ,1 , time ] - loc

  # Return a vector of distances
  return ( dist )
}

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


# The following function takes a current location and produces a fish move
# based on habitat (depth) preference and some process error.
#   Arguments:    Xt = latitude of current position
#                 Yt = longitude of current position
#              sigma = variance of process error
#                wgt = multiplicative weight parameter for gradient potential
#              bathy = list of bathymetry data as produced by makeTopography
#                       function in PBS mapping
#     Returns: delta = vector of jumps in lat/long
fishJumps <- function ( oldPos = c ( 0, 0 ), sigmaT = pi/3, meanD = 0.01,
                        sigmaD = 0.004, bathy = bathyGrad, potSource = potPoint )
{
  # Source parallel functions
  source ( "parFunctions.R" )

  # Extract locations from vector
  Xt <- oldPos [ 1 ]
  Yt <- oldPos [ 2 ]

  # If fish has not yet been released, just pass through
  if ( Xt == 0 | Yt == 0 )
  {
    delta <- oldPos
  } else {
    
    # Lookup grid cell reference
    gridCell <- bathyLookup ( lat = Xt, long = Yt, bathy = bathy )
    xEntry <- gridCell [ 1 ]
    yEntry <- gridCell [ 2 ]

    # Now we want to take finite differences for each coordinate, which we 
    # look up in a table in bathyGrad
    dX <- bathy $ dX [ xEntry, yEntry ]
    dY <- bathy $ dY [ xEntry, yEntry ]
    
    # Now we want to fund the compass heading that this represents, use
    # previously defined cart2pol function, which returns an ordered
    # pair (dist,heading)
    muTheta <- cart2pol ( x = -dX, y = -dY ) [ 2 ]

    # Get random jumps - fish might be able to swim in reverse, but that's
    # okay, it will simulate fish getting spooked or something - fix later
    # when you have time to think about what sigma is in a logN dist
    Theta <- rnorm ( n = 1, mean = muTheta, sd = sigmaT )
    Dist <- rnorm ( n = 1, mean = meanD, sd = sigmaD )

    # Now convert back to cartesian coordinates
    delta <- pol2cart ( r = Dist, t = Theta ) 

    # If still in inlet, add a potential to overcome depth preference
    # Border across which potential function turns off
    # y = -6/7 * x - 405/7
    if ( Yt >= -6/7 * Xt - 409 / 7 )
    {
      # Create complex numbers out of current location and the potential
      # location
      oldPointCplx <- complex ( real = Xt, imaginary = Yt )
      potPointCplx <- complex ( real = potSource [ 1 ], 
                                imaginary = potSource [ 2 ] )

      # Take the difference of the complex numbers for the direction
      potentialCplx <- potPointCplx - oldPointCplx 

      # The argument of the complex difference is the heading we want
      potTheta <- Arg ( potentialCplx )

      # Add a jump in that direction to delta
      delta <- delta + pol2cart ( r = Dist, t = potTheta )
    }

  }
  # return the new position
  return ( delta )
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

# Function for pulling out the first element of a list vector entry - 
# to be lapplied below
#     Arguments: listVec = a list vector containing lists with multiple entries
#                  entry = index of listVec from which to recover element
#                element = the element number of the list entry to be recovered
# 
#     Returns        out = the element being recovered
listElement <- function ( entry = 1, listVec, element = 1 )
{
  # Get element and return
  out <- listVec [[ entry ]] [[ element ]]
  return ( out )
}


boatDEVS <- function ( boat = 1, histData = histList, T = nT )
{
  # Source functions needed in parallel running
  source ( "parFunctions.R" )

  # Set recursion limit deeper
  options ( expressions = 500000 )

  # Recover historic data from entry list.
  layOver <- histData $ layOver
  meanEventTime <- histData $ meanEventTime 
  meanTows <- histData $ meanTows 

  # Create 2D array to hold states for each boat
  # First, names for dimensions
  Bnames <- list ( )
  Bnames [[ 1 ]] <- c ( "latitude", "longitude", "status", "tripLength" )
  Bnames [[ 2 ]] <- 1:T

  # Then make the array to hold boat state variables
  B <- array (  0, dim = c ( 4, T ), dimnames = Bnames )

  # Initialise clock
  t <- 0

  # Initialise lat, long, interFishTime and tripThresh
  nowLat <- 0
  nowLong <- 0
  interfishTime <- 0
  tripThresh <- rpois ( n = 1, lambda = meanTows )

  # Future Event List (FEL) is a data frame that contains 
  # a boat identifier and the time stamp and type of the event
  FEL <- data.frame (   boat.number = integer ( length = 0 ),
                        time.stamp = integer ( length = 0 ),
                        event.type = integer ( length = 0 ),
                        event.time = integer ( length = 0 ),
                        trip.length = integer ( length = 0 ),
                        trip.thresh = integer ( length = 0 ) )

  # 3a. Create first event, add to FEL
  # Generate a timestamp - make Poisson so it happens realistically early
  tFirst <- runif ( n = 1, min = 1, max = min ( T, 500 ) )

  # First event type
  eFirst <- fishing

  # Generate first towTime (event time variable named eTime)
  eTimeFirst <- towDuration ( lam = meanEventTime )

  # counter variable for event number
  eCount <- 1

  # set first tow counter to 1
  nTows <- 1

  # Add event to FEL - this will be the first fishing event
  FEL [ eCount, ] <- c ( boat, tFirst, eFirst, eTimeFirst, nTows, tripThresh )
  
  # 4. loop for DEVS:
  # while (t < nT) do {
  while ( t < T )
  {
    #       Forward clock to next event
    t <- FEL [ eCount , ] $ time.stamp
    eTime <- FEL [ eCount, ] $ event.time
    nTows <- FEL [ eCount, ] $ trip.length

    if ( t + eTime - 1 > T ) { break }    

    # Process next event on FEL and update boat state
    B [  
        , t:min ( T, t + eTime - 1 )
      ] <- eventProc ( 
                       event = eCount, 
                       fel = FEL,
                       nowlat = nowLat,
                       nowlong = nowLong,
                       ntows = nTows,
                       histData = histData,
                       b = B
                     )

    # Schedule a future event, based on the current event
    FEL <- eventSched ( time = t, event = eCount,
                        fel = FEL, ntows = nTows,
                        histData = histData )

    # 4d. increase event counter
    eCount <- eCount + 1
    # goto 4
  }

  # Combine boat states and FEL into a list for returning
  outList <- list ()
  outList $ state <- B
  outList $ FEL <- FEL

  # return outList
  return ( outList )
}


# Function to run a discrete event simulation given some historic
# fishing data on trip length, event duration etc. Essentially
# a wrapper for the function boatDEVS
#   Args:   dt = a line of a data.table of yeardata, as produced by
#                 histAnalysis.R
#   Rets:  out = a list containing a system state array and a FEL
yearDEVS <- function ( dt = yearAverages, T = nT )
{
  # Recover number of boats
  nB <- dt $ nVessels

  # recover historic averages to pass to the DEVS function
  histList <- list ( )
  histList $ layOver <- round ( 24 * dt [ , aveLayoverDays ] )
  histList $ interFish <- round ( dt [ , aveInterFish ] )
  histList $ meanEventTime <- round ( dt [ , aveDurHours ] )
  histList $ minDepth <- -1. *  dt [ , maxDepth ]
  histList $ maxDepth <- -1. * dt [ , minDepth ]
  histList $ meanTows <- round ( dt [ , aveTows ] )

  # Call boatDEVS function to run DEVS and return states and a FEL
  # for each boat. Parallelised
  clust <- makeCluster ( spec = nCores )

  # Export variables to clusters
  clusterExport ( cl = clust, varlist = parList )

  # Run boatDevs in parallel on each boat
  boatList <- parLapply ( cl = clust, X = 1:nB, fun = boatDEVS,
                          histData = histList, T = T )

  # Non parallel version for testing errors
  # boatList <- lapply ( X = 1:nB, FUN = boatDEVS,
  #                         histData = histList, T = nT )

  # Stop cluster to avoid ghost processes.
  stopCluster ( cl = clust )
  
  # Bind FEL and states
  # states are a little easier, as sapply has automatic simplification
  # to an array
  states <- sapply ( X = 1:nB, FUN = listElement, listVec = boatList,
                     element = 1, simplify = "array" )

  # Extract the FELs into a list.
  FELs <- lapply (  X = 1:nB, FUN = listElement, listVec = boatList, 
                    element = 2 )

  # Bind FELs into a single list.
  FEL <- rbindlist ( FELs )
  FEL <- as.data.table ( FEL )

  # Make a list for the FEL and state array, add them and return
  out <- list ()
  out $ FEL <- FEL
  out $ B <- states

  return ( out )
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
  diff2 <- diff * diff
  diffDist <- diff2 [ ,1 ] + diff2 [ , 2 ]

  # Look for fish with diffDist < radius
  closeFish <- which ( diffDist < 0.005 )

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
  # Source parallel functions
  source ( "parFunctions.R" )

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

# Function to run detections code over given fish state variables
# and a FEL and state variables for fishing effort, as produced by 
# simulator.R. The output will be used in bootstrap samples of 
# detections.
#   Argument: boatDat = a list containg a FEL and boat states
#              fArray = array of fish state variables
#   Returns:   detect = a data.table of detections including time,
#                       fish number and location in lat/long
detectFish <- function (  boatDat = fishingYears [[ 1 ]], 
                          fArray = F [ 1:1000, , ], radius = 0.005 )
{
  # Get FEL from yearData
  FEL <- boatDat [ 1 ] $ FEL

  # Recover number of boats
  nB <- max ( FEL[ , boat.number ] )

  # Get boat states from yearData - the lapply functions in simulation upset
  # the structure of the B array, so we rename the third dimension - and
  # get used to a new variable order. I'm sure this is avoidable.
  B <- boatDat [ 2 ] $ B
  dimnames ( B ) [[ 3 ]] <- 1:nB

  # Reduce to timesteps that fishing is occuring
  fishingTimes <- rowIdx ( mat = B [ 3, , ] )

  # For each timestep and each boat, compare the locations of boats and fish.
  # For any comparison that is small enough (w/in 500m = 0.005deg), record the 
  # fish number, location and time step. Add first time step of every fish to 
  # detection frame as these are tagging locations and will contribute to
  # habitat knowledge.

  # Get a list of detections, where each entry is the table of detections
  # at each time - in parallel
  clust <- makeCluster ( spec = nCores )
  detList <- parLapply (  cl = clust, X = fishingTimes, fun = tDetect, 
                          bArray = B, 
                          fArray = fArray, radius = radius )
  stopCluster ( cl = clust )

  # Bind list together
  detections <- rbindlist ( detList )

  # return detections
  return ( detections )
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

# Function to randomly sample fish and compare detections
# to true home range - for aggregate samples and inlet specific
# fish.
#   Arguments:         N = sample size
#             detections = detection data to use (based on historic fishing
#                          effort)
#                indices = indices of fish in each inlet - bad planning
#                 states = array of fish state variables
#   Returns:   areaPolys = list of estimated and simulated habitat polygons
polyFunc <- function ( N = 20, detections = yearDetections [[ 1 ]],
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
  mcpAggDet <- mcp ( xy = sampDet [ ,  c("boat.lat", "boat.lon"), 
                                                with = FALSE ], 
                     id = sampDet [ , agg ], percent = 100 )
  # Find area of MCP of true locations of agg sample
  mcpAggTrue <- mcp (  xy = trueDT [ , c("fish.lat", "fish.lon"), 
                                                 with = FALSE ], 
                                id = trueDT [ , agg ],
                                percent = 100 )

  agg <- list ( )
  agg $ detPoly <- mcpAggDet
  agg $  truePoly <- mcpAggTrue

  polyList <- vector ( mode = "list", length = nInlets + 1 )

  polyList [[ 1 ]] <- agg

  
  for ( i in 1:nInlets )
  {
    # Reduce to fish from inlet i
    sampDetInlet <- sampDet [ inlet == i ]
    trueDTInlet <- trueDT [ inlet == i ]

    # Find area of minimal convex polygon around detections
    mcpDetInlet <- mcp ( xy = sampDetInlet [ ,  c("boat.lat", "boat.lon"),
                                                      with = FALSE ],
                                  id = sampDetInlet [ , agg ], 
                                  percent = 100 )
    # Find area of MCP around true locations
    mcpTrueInlet <- mcp (  xy = trueDTInlet [ , c( "fish.lat", 
                                                            "fish.lon" ),
                                                       with = FALSE ],
                                    id = trueDTInlet [ , agg ],
                                    percent = 100 )
    # Fill row of areas matrix
    polyList [[ i + 1 ]] <- list ( mcpDetInlet, mcpTrueInlet )
  }

  return ( polyList )
}

# The following functions are wrappers for the compareFunc function,
# in order to allow for swift parallelisation and returning
# of lists in a decent structure for analysis.
sampSizeFun <- function ( sampleSize = 25, nSamp = nSamples, 
                          detDT = yearDetections [[ 1 ]],
                          idx = inletFishIdx, states = F )
{ 
  # # Count inlets
  # nInlets <- length ( idx )  

  # # Create a matrix to hold fish sample numbers
  # sampMat <- matrix ( 0, nrow = nSamples, ncol = sampleSize * nInlets )

  # for ( i in 1:nSamples )
  # {

  #   # Take samples in each inlet
  #   # Make matrix to hold samples
  #   fishSamp <- matrix ( 0, nrow = nInlets, ncol = sampleSize )

  #   # Take samples
  #   for ( i in 1:nInlets )
  #   { 
  #     ind <- indices [[ i ]][[ 1 ]]
  #     fishSamp [ i, ] <- sample ( x = ind, size = sampleSize )
  #   }

  #   # Aggregate into a vector
  #   aggSamp <- as.vector ( fishSamp )
  #   sampMat [ i, ] <- aggSamp
  # }

  # areasArray <- parApply ( cl = clust, X = sampMat, MARGIN = 1,
  #                           FUN = repWrapper, )
  source ( "parFunctions.R" )
  require ( adehabitat )
 
  areasArray <- sapply ( X = 1:nSamp, 
                            FUN = repWrapper, 
                            sampSize = sampleSize, det = detDT,
                            idx = idx, stateVars = states,
                            simplify = "array" )
  filename <- paste ( "areasArraySample", sampleSize, ".RData", sep = "" )
  save ( areasArray, file = filename )
  return ( areasArray )
}

repWrapper <- function ( times = 1, sampSize = 25, det = detDT, 
                         idx = idx, stateVars = states  )
{
  source ( "parFunctions.R" )
  require ( adehabitat )
  n <- times
  areas <- compareFunc (  N = sampSize, detections = det, 
                          indices = idx, states = stateVars )

  return ( areas )
}

