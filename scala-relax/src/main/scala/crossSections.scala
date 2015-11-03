
package crossSections

import atmosphere.Atmosphere
import scala.util.Random
import math._

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

  // get lambda for atom-atom and atom-molecule
  def getLambda(target: String): Double = {
    1.0
  }

  def getReducedMass(target: String): Double = {
    2.0
  }

  // universal amplitude
  def universalAmplitude(energy: Double, angle: Double, target: String): Double = {
    val tau0: Double = 50.12
    val c1: Double = -0.13
    val c2: Double = 1.0 
    val c3: Double = 2.7 
    val c4: Double = 10.0 
    val c5: Double = 2.04 
    val c6: Double = -0.03 
    val c7: Double = 32.3
    val lam: Double = getLambda(target) 
    val tau: Double = energy*angle/getReducedMass(target)
    val x: Double = math.log(tau)
    val c: Double = lam/(angle*math.sin(angle))
    var amp: Double = 0.0
    if (tau >= tau0) {
      amp = c*math.exp( c1*x + c2*x + c3 )
    }
    else {
      amp = c*c4*math.exp( c5 + c6*x ) + c7
    }
    amp
  }

  // universal scattering random angle
  def universalScatteringAngle(energy: Double, target: String): Double = {
    universalTotalCrossSection(energy, target)
  }

  def universalTotalCrossSection(energy: Double, target: String): Double = {
    val N_angles: Double = 10000
    val da: Double = math.Pi*2.0/N_angles  
    da   
  }

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


