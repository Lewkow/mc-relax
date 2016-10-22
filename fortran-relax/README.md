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

* ```Escape_Energy_Distribution.dat```
  *
