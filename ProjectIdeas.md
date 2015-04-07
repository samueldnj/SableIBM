## Project outline
21/1/15

Nearshore inlets provide a nursery area for juvenile sablefish. The yearly inlet survey tags and releases sablefish in order to track their movements as they migrate and recruit to the fisheries in offshore BC and beyond. These tags are collected when the sablefish are intercepted in the surveys and fisheries which operate on the coast and offshore and also in international waters all along the west coast. The recovery and release locations provide data about the migration habits of the sablefish, and inform movement probability and abundance estimation for the management of the population.

- Passive Tags
  - Drawbacks
    - Only two data points
      - Not a lot of resolution
    - harvest is required for detection
      - Harvest rate ~9%, though recoveries average out to 17% over the last 20 years in all fisheries, ~10% in trawl
  - Benefits:
    - Cheap

- Alternative: Acoustic Tags
  - Benefits:
    - Don't require harvest for detection
    - Can provide several data points, allowing for (hopefully) a more complex picture of migratory habits
    - Positive externality: possibility of detecting marine stage salmon, and other species with acoustic tags
  - Drawbacks
    - More costly
  - Given that the tags don't require harvest for detection, perhaps we can get away with putting out fewer tags and still retain a 17% detection rate or better.

** Question: What sample size of tags in inlet fish is needed to retain or improve a 17% detection rate (detection efficiency) over a 20 year period?**
  - Might need a simplification to a 1 year timescale, and match that efficiency.

#### Scoping
  - Start with just the trawl fishery (in Hecate Strait).
    - Coarse grid
  - Movement
    - Markov Processes
    - Fish move right to left on average
      - Enter at random times
      - Leave when they reach left edge
      - Move at a realistic speed
    - Boats move around randomly
      - Enter and leave at random times (based on realistic trip days)
  - Encounters
    - Happen when a fish and boat are in the same grid cell
  - Detections
    - When an encounter occurs, another markov process on a smaller scale is triggered
    - On a grid where cells are 500m squares (or whatever the detection radius is)
    - Boat is on a deterministic straight line set
      - Requires heading, speed
    - Fish on a markov process
    - Opportunity for detection if a fish inhabits a cell in the neighbourhood of the boat, determined by a Bernoulli trial, depends on
      - Originally:
        - engine noise
      - Ultimately:
        - habitat (coral blocking signal)
        - battery life
        - temperature and turbulence?

#### Initial simplifying assumptions:
**Encounters:**
  - Straight line tows
  - Only horizontal distance from detection gear matters (fish and nets on bottom)
  - Fish always moving at the same speed (base time steps on this)
    - Is this slower than the boat?
  - No fishing mortality or natural mortality to begin with
    - Add this later, this is necessary for publication, and it's just an exponential model
**Detections:**
- Spherical radius of detection
- Detection 100% likely, unless coral/sponge reef present.
  - Then use Bernoulli trial, with some plucked out probability
- engine noise effect only depends on depth (decays)
- no temperature and turbulence originally
- Don't allow multiple detections in DE computation


#### Future additions
Realism:
  - Fishing and natural mortality
  - Multiple detections
  - Battery failure
  - Tag pulses - detection efficiency is then detections/pulses
  - Inlets, Strait and Shelf areas
    - Different fishing techniques in each
    - Trap and longline are stationary receivers
    - Markov processes to move between areas with real movement probabilities given by passive tag data.
  - Fit Markov chains to real boats? That could be fun.

Bayesianism:
  - Is there some way to use the stationary distributions and/or a Bayesian approach?


Method: Build a model with the above scope (make modular so future realism can be added as needed), compute detection efficiency (which I'm pretty sure is just detections/tags) and run monte-carlo simulations to get a distribution.

Results and Discussion:

Future work:
- Economic analysis
  - Economic benefits? CBA? Potential Pareto Efficiency?
  - Might form part of a larger analysis for the thesis
- Using historic mark-recap data as a prior for movement probabilities, then updating with the results of a pilot study done with acoustic tags?
- Independence of tags
  - Schooling


## Project Ideas - There are no bad ideas in brainstorming

1. Individual based fishing boat model, fishing boats are individuals, changing trip movements based on habitat and non-target species encounters (rewards and penalties). Analysis of what real-time data would do, maybe collected without habitat interference or fish interception (video, hydrophone). <- STELLA

2. Simulated juvenile bycatch model, spatially auto-correlated. Might feed into 1.

3. Fishing success in individual based model feeds into a real time stock assessment, which updates an index of abundance, catch data etc and produces an updated stock assessment - perhaps for all stocks in the catch.

I fear the scope is getting away from me...

4. Power analysis for hydrophone detection of tagged fish in Coastal/Offshore fishery - markov chains for fish, boats, probability of detection etc. Make modular so you can include different gear types, more fish, more boats etc.
  - Question: Can we replace passive tags with fewer acoustic tags on juvenile sablefish to get presence absence and movement without interception, assuming that every fishing boat has detection gear attached to their fishing gear?

