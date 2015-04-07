#---------------------------------------------------------------------
# 
# REM612 Project Script - Historical Data Analysis
# 
# Author: Samuel Johnson
# Date: 19/1/15
# 
# In order to generate realistic fishing pressure in the simulation,
# we look at the summary stats for trip length and frequency in each
# year for different boats.
# 
#---------------------------------------------------------------------

# Load required packages
require ( dplyr )

# Load the historical data
load ( file = "histData.RData" )

# Now, let's look at the number of trips per year for each vessel,
# and the average number of trip days
Fish <- group_by ( Fish, year, vesselID, tripID )

FishbyTrip <- summarise ( Fish, 
                          totalTows = max ( totalTows ), 
                          tripLength = max ( tripDays ),
                          towDuration = sum ( towDuration ),
                          nEvents = n ( ),
                          maxDepth = max ( depth ),
                          minDepth = min ( depth )
                           )

FishbyVessel <- summarise ( FishbyTrip, nTrips = n ( ),
                            totalTows = sum ( totalTows ), 
                            totalEvents = sum ( nEvents ),
                            # aveTows = mean ( totalTows ),
                            # aveLength = mean ( tripLength ),
                            seaDays = sum ( tripLength ),
                            fishTime = sum ( towDuration ),
                            maxDepth = max ( maxDepth ),
                            minDepth = min ( minDepth ) )

FishbyYear <- summarise ( FishbyVessel,
                          nVessels = n ( ),
                          meanTrips = mean ( nTrips ),
                          aveTotalDays = mean ( seaDays ),
                          aveTows = sum ( totalTows ) / sum ( nTrips ),
                          aveDuration = sum ( fishTime ) / sum ( totalEvents ),
                          maxDepth = max ( maxDepth ),
                          minDepth = min ( minDepth ) )

FishByYear <- mutate ( FishbyYear, 
                       aveTripDays = aveTotalDays / meanTrips,
                       aveLayoverDays = ( 365 - aveTotalDays ) / meanTrips,
                       aveDurHours = aveDuration / 60,
                       aveInterFish = ( aveTripDays * 24 ) / aveTows - 
                                      aveDurHours 
                     )

FishByYear <- FishByYear [ year != "2000/2001" ]

save ( FishByYear, Fish, file = "histAveByYear.RData")
