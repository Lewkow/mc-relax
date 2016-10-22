# Fortran Monte Carlo Relax

## Input Files
* fortran-relax/Inputs/keys.in
  * What kind of simulations to do
    * Local Interstellar Medium (LISM) simulation
      * Homogeneous medium transport
    * Plantary simulation
      * Position dependent density medium (atmosphere) transport
  * What kind of secondary calculations to do (on(1) or off(0))
    * ENA Production
      * Determine altitude dependent ENA production rates through MC transport simulation using charge-exchange (CX) cross sections
      * Average Scattering Angle calculations through recording all transport parameters and performing aggregations
      * Secondary Hot (SH) MC simulation through tracking nascent SH produced through large angle collisions with MC particle
      * Calculate the escape transparency of SH
* fortran-relax/Inputs/planet_keys.in
  * Used when planet calculation is set to 1 in fortran-relax/Input/keys.in
  * Set parameters used in the simulation and output data files to write

## Usage

Clone this github repo
  ```
  $ git clone https://github.com/Lewkow/mc-relax.git
  ```

The tables used for the Monte Carlo simulation are too large to put in the github repository. They are available upon request from the author at nlewkow@gmail.com

Once the tables tarball has been downloaded, it needs to be untared and the files moved to the fortran-relax/Tables directory
  ```
  $ mv scattering_dat.tar.gz /vagrant/fortran-relax/Tables
  $ cd /vagrant/fortran-relax/Tables
  $ tar -zxf scattering_dat.tar.gz
  $ mv ./scattering_dat/* ./
  ```

If you have MPI and gfortran working already on your system (_mpif90_ and _mpirun_ specifically) you can run the default input files with 2 processors as
  ```
  $ cd fortran-relax/Execute
  $ ./BUILD 2
  ```
If you do not have these libraries properly installed there is a Vagrant box in the repository which runs Ubuntu64 which is automatically provisioned to run mc-relax.

Install Vagrant at http://www.vagrantup.com/downloads

Once installed, the vagrant box can be started
  ```
  $ cd mc-relax
  $ vagrant up
  ```

Once the Vagrant box is started you can ssh into the box
  ```
  $ cd mc-relax
  $ vagrant ssh
  ```

Now in the Vagrant box, you can run mc-relax
  ```
  $ cd /vagrant/fortran-relax/Execute
  $ ./BUILD 2
  ```

You can turn off your vagrant box with
  ```
  $ cd mc-relax
  $ vagrant halt
  ```

You can also completely destroy your vagrant box state, removing any changes you made to that box with
  ```
  $ cd mc-relax
  $ vagrant destroy
  ```

## Output data descriptions
Several data files are saved to mc-relax/fortran-relax/Data/ upon running the relaxation simulation.

### 1-Dimensional Distributions
These distributions have a single value per particle that applies to the condition (escape or thermalization).

* ```Escape_Energy_Distribution.dat```
  * [ escape energy (eV) ]
  * The energy of every particle that escapes the atmosphere
* ```Escape_NColl_Distribution.dat```
  * [ number of collisions before escape ]
  * The number of collisions before escape for every particle that escapes the atmosphere
* ```Thermalization_Height_Distribution.dat```
  * [ thermalization altitude (m) ]
  * The altitude of thermalization for every particle that thermalizes in the atmosphere
* ```Thermalization_Time_Distribution.dat```
  * [ thermalization time (s) ]
  * The amount of time, in the planet frame, that it takes to thermalize for every particle that thermalizes in the atmosphere

### Bin Variable Index Files
Indices are used for the value of the bin in 2-dimensional space.
The values for the bin indices are found in the files ```planet_X_VARIABLE.dat``` where __VARIABLE__ is the variable of interest. For example, for the file ```planet_Ncoll_vs_time.dat``` the corresponding index files are ```planet_X_Ncoll.day``` and ```planet_X_time.day```.

* All files available for bin variable indices
  * ```planet_X_Ncoll.dat```
    * Number of collisions
  * ```planet_X_energy.dat```
    * Energy (eV)
  * ```planet_X_energy_loss.dat```
    * Energy loss (eV)
  * ```planet_X_height.dat```
    * Altitude (m)
  * ```planet_X_phi.dat```
    * Phi (deg)
  * ```planet_X_theta.dat```
    * Theta (deg)
  * ```planet_X_time.dat```
    * Time (sec)
  * ```planet_X_unit_vel.dat```
    * Unit velocity
  * ```planet_X_xy.dat```
    * Position (m)

### Click Dependent 2-Dimensional Distributions
These distributions are 2-dimensional distributions based on the 2 variables assigned to the file (for example ```planet_Ncoll_vs_time.dat``` contains time dependent distributions for number of collisions vs lab frame time).

* All of the time dependent distributions are labeled in files ```planet_VARIABLE1_vs_VARIABLE2.dat``` and have the following format:
  * [ simulation click, VARIABLE1 index, VARIABLE2 index, bin value ]
  * bin value is normalized such that if you sum over the entire 2-dimensional space the result is 1, assuming at least 1 particle is within the space.
  * Particles are free to leave the 2-dimensional space and are then not counted in the data

* All files available for click dependent 2-dimensional distributions
  * ```planet_Ncoll_vs_time.dat```
  * ```planet_energy_vs_time.dat```
  * ```planet_height_vs_Ncoll.dat```
  * ```planet_height_vs_energy.dat```
  * ```planet_height_vs_energy_loss.dat```
  * ```planet_height_vs_time.dat```
  * ```planet_height_vs_vert_vel.dat```
  * ```planet_vert_vel_vs_time.dat```
