# mc-relax
Monte Carlo atomic and molecular particle simulator: Simulating the relaxation of particles in astrophysical environments, one collision at a time!
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
* Evaluate use of classes in the code. Refactor as needed

  * _particle_ class
    * a particle can be an atom, ion or molecule
    * __attributes__
      * name, lab frame energy, position, velocity, mass, isotope, charge
    * __methods__
      *

  * _atmosphere_ class
    * contains information about the atmosphere being simulated: used to determine target particles, collision lengths
    * __attributes__
      * atmosphere_particles = List(_particle_)
      * atmosphere_densities = Map()
    * __methods__
      * get_atmosphere_densities(current_atmosphere_position)
      * get_cross_sections(current_atmosphere_position, current_projectile: _particle_)

  * _collision_ class
    * a class that is created when a projectile _particle_ and target _particle_ are determined
    * __attributes__
      * projectile: _particle_
      * target: _particle_
    * __methods__
      * (projectile_after: _particle_, target_after: _particle_) = collide



## scala-relax
Scala version of the Monte Carlo particle simulator.
Uses Apache Spark for parallelization. 


## fortran-relax
Fortran 90 version of the Monte Carlo particle simulator.
Uses openMPI for parallelization.

