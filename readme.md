<!-- A readme for the REM612 project repository -->

## REM612 simulation modelling project

1. Introduction
2. Structure
3. Future work

### 1. Introduction

Here is the repository for the REM612 simulation modelling project. The project is a power analysis for an acoustic mark-recapture program. It does so by running parallel markov processes, one for each mobile receiving station and one for the tagged animal, moving around on a rectangular grid. When the animal and a receiver are in the same grid cell, a Bernoulli process, with probability of success given by the detection efficiency of acoustic tags, determines whether there has been a succesful detection.

### 2. Structure

The project is split into two R scripts. The first, JohnsonProject.R contains the control script which sets model parameters and calls the functions. The second, functions.R, contains a function for creating a transition matrix and running the simulation.

### 3. Future Work

Future work:
  1. Add a general number of detectors as a parameter
  2. Add more tags, make independent originally, then dependent depending on herd/schooling behaoviour.
  3. Add different management areas, with a markov process determining movement within areas and another determining movement between.
  4. Tag's ability to leave the patrolled areas entirely.
  5. Investigate how the graph game "cops and robbers" can relate.