# mc-relax
Monte Carlo atomic and molecular particle simulator: Simulating the relaxation of particles, one collision at a time!
http://arxiv.org/abs/1404.5986

# version 1
Features:
* Docker image for scala-relax
* All cross sections need to be computed before hand and saved to files in the proper format
* Cross sections are read before the simulation and stored in memory for use during the simulation
* Hard sphere cross sections, universal cross sections, and quantum cross sections available from previous use

# to-do
For version 1 (v1) of mc-relax the following features need to be implemented:
* Create Dockerfile that creates an image for use with Apache Spark (latest version) and Scala
* Create functionality to check for local files containing cross sections needed for simulation in current atmosphere
  * If files do not exist, alert the user and inform them what cross sections are needed and how to compute them with scatter-calc
* Read cross sections into memory from files and use them in simulation
* Determine proper format for saving simulation results

## fortran-relax
Fortran 90 version of the Monte Carlo particle simulator. 
Uses openMPI for parallelization.


## scala-relax
Scala version of the Monte Carlo particle simulator.
Uses Apache Spark for parallelization
