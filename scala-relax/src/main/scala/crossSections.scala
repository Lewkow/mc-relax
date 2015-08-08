
package crossSections

import atmosphere.Atmosphere
import scala.util.Random

class CrossSections extends Serializable {

  // random object generator
  var randy: Random = new Random()

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
    gaussianScatteringAngle(10.0d) 
  }

  def gaussianScatteringAngle(mean: Double): Double = { randy.nextGaussian() + mean }

  // get total cross section for projectile-target collision given lab energy [eV]
  // lab-frame energy [eV]
  // returns -> Total cross section [m^2]
  def getTotalCrossSection(energy: Double, target: String): Double = {
    1.0e-15
  }

}


