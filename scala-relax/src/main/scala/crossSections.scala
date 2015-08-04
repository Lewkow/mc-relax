
package crossSections

import atmosphere.Atmosphere


class CrossSections {

  // object with current atmosphere
  var currentAtmosphere: Atmosphere = new Atmosphere

  // max energy for simulation [eV]
  var maxEnergy: Double = 0.0d

  // energy stepsize required [eV]
  var dEnergy: Double = 0.1d

  // projectile particle name
  var projectileName: String = ""


  // constructor for object
  def setParameters(projectile: String, inAtmosphere: Atmosphere, maximumEnergy: Double) {
    projectileName = projectile
    currentAtmosphere = inAtmosphere
    maxEnergy = maximumEnergy 
  }

  // check to see if cross sections needed are in database
  def haveCrossSections: Boolean = {
    true
  }

  // get scattering angle for projectile-target collision given lab energy [eV]
  // lab-frame energy [eV]
  // returns -> Random scattering angle [deg]
  def getScatteringAngle(energy: Double, target: String): Double = {
    0.5d
  }

  // get total cross section for projectile-target collision given lab energy [eV]
  // lab-frame energy [eV]
  // returns -> Total cross section [m^2]
  def getTotalCrossSection(energy: Double, target: String): Double = {
    1.0d-15
  }



}


