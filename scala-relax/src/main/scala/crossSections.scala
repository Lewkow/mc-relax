
package crossSections

import atmosphere.Atmosphere
import scala.util.Random
import math._

class CrossSections extends Serializable {

  val use_Universal: Boolean = false

  // number of angles for numeric integrations
  val N_angles: Int = 10000

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

  def getAnyLambda(proj: String, targ: String): Double = {
    val projMol = currentAtmosphere.isMolecule(proj)
    val targMol = currentAtmosphere.isMolecule(targ)
    if (projMol || targMol) {
      1.4
    } else {
      1.0
    } 
  }

  // get lambda for atom-atom and atom-molecule
  def getLambda(target: String): Double = {
    getAnyLambda(projectileName, target)
  }

  def getAnyReducedMass(projectile: String, target: String): Double = {
    val projMass: Double = currentAtmosphere.getParticleMass(projectile)
    val targMass: Double = currentAtmosphere.getParticleMass(target)
    (projMass * targMass) / (projMass + targMass)
  }

  def getReducedMass(target: String): Double = {
    val projMass: Double = currentAtmosphere.getParticleMass(projectileName)
    val targMass: Double = currentAtmosphere.getParticleMass(target)
    (projMass * targMass) / (projMass + targMass)
  }

  // universal amplitude
  // Inputs:
  //   energy CM [eV]
  //   angle CM [deg]
  //   target particle
  // Outputs:
  //   differntial cross section CM [a0^2]
  def anyUniversalAmplitude(energy: Double, angle: Double, projectile: String, target: String): Double = {
    val tau0: Double = 50.12
    val c1: Double = -0.13
    val c2: Double = 1.0 
    val c3: Double = 2.7 
    val c4: Double = 10.0 
    val c5: Double = 2.04 
    val c6: Double = -0.03 
    val c7: Double = 32.3
    val lam: Double = getAnyLambda(projectile, target) 
    val tau: Double = energy*angle/getAnyReducedMass(projectile, target)
    val x: Double = math.log(tau)
    val c: Double = lam/(angle*math.sin(toRadians(angle)))
    var amp: Double = 0.0
    if (tau >= tau0) {
      amp = c*math.exp( c1*x*x + c2*x + c3 )
    }
    else {
      amp = c*c4*math.exp( c5 + c6*x ) + c7
    }
    amp
  }

  // universal amplitude
  // Inputs:
  //   energy CM [eV]
  //   angle CM [deg]
  //   target particle
  // Outputs:
  //   differntial cross section CM [a0^2]
  def universalAmplitude(energy: Double, angle: Double, target: String): Double = {
    anyUniversalAmplitude(energy, angle, projectileName, target)
  }


  def anyUniversalScatteringAngle(energy: Double, projectile: String, target: String): Double = {
    val tcs: Double = anyUniversalTotalCrossSection(energy, projectile, target)
    val theta_i: Double = 0.01
    val theta_f: Double = 170.0
    val da: Double = toRadians((theta_f-theta_i)/N_angles.toDouble)
    var rho: Double = 0.0
    var r: Double = randy.nextDouble()
    var gotIt: Boolean = false
    var i: Int = 0
    var scattering_theta: Double = 0.0
    while (rho < r) {
      scattering_theta = toDegrees(da*i + theta_i)
      i += 1
      val amp = anyUniversalAmplitude(energy, scattering_theta, projectile, target)
      rho += (da*2.0*math.Pi*math.sin(toRadians(scattering_theta))/tcs)*amp
    }
    scattering_theta
  }

  // universal scattering random angle
  // Inputs: 
  //   energy CM [eV]
  //   target particle
  // Outputs: 
  //   scattering angle CM [deg]
  def universalScatteringAngle(energy: Double, target: String): Double = {
    anyUniversalScatteringAngle(energy, projectileName, target)
  }

  def anyUniversalTotalCrossSection(energy: Double, projectile: String, target: String): Double = {
    val theta_i: Double = 0.01
    val theta_f: Double = 170.0
    val da: Double = toRadians((theta_f-theta_i)/N_angles.toDouble)
    var tcs: Double = 0.0
    var i: Int = 1
    for (i <- 0 until N_angles) {
      val theta: Double = toDegrees(da*i + theta_i)
      val amp: Double = anyUniversalAmplitude(energy, theta, projectile, target)
      tcs += da*2.0*math.Pi*math.sin(toRadians(theta))*amp
    }
    tcs
  }

  def universalTotalCrossSection(energy: Double, target: String): Double = {
    anyUniversalTotalCrossSection(energy, projectileName, target)
  }

  // check to see if cross sections needed are in database
  def haveCrossSections: Boolean = {
    true
  }

  // get scattering angle for projectile-target collision given lab energy [eV]
  // lab-frame energy [eV]
  // returns -> Random scattering angle [deg]
  def getScatteringAngle(energy: Double, target: String): Double = {
    if (use_Universal) {
      universalScatteringAngle(energy, target)
    } else {
      gaussianScatteringAngle(10.0d) 
    }
  }

  def gaussianScatteringAngle(mean: Double): Double = { randy.nextGaussian() + mean }

  // get total cross section for projectile-target collision given lab energy [eV]
  // lab-frame energy [eV]
  // returns -> Total cross section [m^2]
  def getTotalCrossSection(energy: Double, target: String): Double = {
    if (use_Universal) {
      universalTotalCrossSection(energy, target)
    } else {
      1.0e-15
    }
  }

  def toDegrees(rad: Double): Double = { rad*180.0d/math.Pi }

  def toRadians(deg: Double): Double = { deg*math.Pi/180.0d }

}


