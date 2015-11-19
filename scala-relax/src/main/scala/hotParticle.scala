package hotParticle

import scala.util.Random
import scala.math._
import scala.collection.mutable._
import atmosphere.Atmosphere
import crossSections.CrossSections

class Vector extends Serializable {
  var x: Double = 0
  var y: Double = 0
  var z: Double = 0

  def setup(x0: Double, y0: Double, z0: Double) {
    x = x0; y = y0; z = z0
  }

  override def toString = s"($x, $y, $z)"

  def update(x0: Double, y0: Double, z0: Double) { 
    // println(s"($x, $y, $z)\n($x0, $y0, $z0)")
    x = x0; y = y0; z = z0 
  }

  def toTuple: (Double, Double, Double) = {(x,y,z)}

}

class HotParticle extends Serializable {

  var projectileName: String = ""
  var projecileMass: Double = 0.0
  var currentPosition = new Vector
  var currentVelocity = new Vector
  var currentEnergy: Double = 0.0
  var currentTime: Double = 0.0
  var numberOfCollisions: Int = 0
  var numberOfClicks: Int = 0
  var collisionProbability: Double = 0.0
  var generation: Int = 0
  var randy: Random = new Random()
  var atmosphere: Atmosphere = new Atmosphere
  var crossSections: CrossSections = new CrossSections
  var trans: Transport = new Transport

  def setParameters(name: String, 
                    mass: Double, 
                    initPosition: (Double, Double, Double),
                    intitVelocity: (Double, Double, Double),
                    particleGeneration: Int,
                    inputAtmosphere: Atmosphere) {
    projectileName = name
    projecileMass = mass
    currentPosition.setup(initPosition._1,initPosition._2,initPosition._3)
    currentVelocity.setup(intitVelocity._1,intitVelocity._2,intitVelocity._3)
    currentEnergy = getEnergy
    generation = particleGeneration
    atmosphere = inputAtmosphere
    // maximum energy for cross sections
    val max_energy: Double = 10000.0d
    crossSections.setParameters(projectileName,atmosphere,max_energy)
  }

  def getEnergy: Double = {
    math.sqrt(math.pow(currentVelocity.x,2) +
              math.pow(currentVelocity.y,2) + 
              math.pow(currentVelocity.z,2))*0.5d*projecileMass
  }

  def getSpeed: Double = {
    math.sqrt(math.pow(currentVelocity.x,2) +
              math.pow(currentVelocity.y,2) + 
              math.pow(currentVelocity.z,2))
  }

  def printNumberOfCollisions {
    println(numberOfCollisions.toString + " total collisions")
  }

  def getCollisionProbability { 
    collisionProbability = 100.0d*numberOfCollisions.toDouble/numberOfClicks.toDouble
  }

  def printCollisionProbability {
    val p: Double = 100*numberOfCollisions.toDouble/numberOfClicks.toDouble
    println(p.toString + " collisions per click")
  }

  def exitConditions: Boolean = {
    if (currentEnergy < 5.0d) {
      // println("clicks     -> "+numberOfClicks.toString)
      // println("collisions -> "+numberOfCollisions.toString)
      true 
    } else {
      false
    }
    // if (numberOfCollisions > 100) true
    // else false
  }

  def collisionUpdate(theta: Double, phi: Double, dx: Double) {
    val thetaRads: Double = toRadians(theta)
    val phiRads: Double = toRadians(phi)

    // update position for straight line before collision
    currentPosition.update(currentPosition.x + currentVelocity.x*dx,
                           currentPosition.y + currentVelocity.y*dx,
                           currentPosition.z + currentVelocity.z*dx)

    if (currentVelocity.z == 1.0d) {
      currentVelocity.update(math.sin(thetaRads)*math.cos(phiRads),
                             math.sin(thetaRads)*math.sin(phiRads),
                             math.cos(thetaRads))
    } else if (currentVelocity.z == -1.0d) {
      currentVelocity.update(math.sin(thetaRads)*math.cos(phiRads),
                             -math.sin(thetaRads)*math.sin(phiRads),
                             -math.cos(thetaRads))
    } else {
      val const: Double = math.sqrt(1.0d - math.pow(currentVelocity.z,2.0))
      val x: Double = (math.sin(thetaRads)*(currentVelocity.x*currentVelocity.z*math.cos(phiRads)-currentVelocity.y*math.sin(phiRads)))/const +
                      currentVelocity.x*math.cos(thetaRads)
      val y: Double = (math.sin(thetaRads)*(currentVelocity.y*currentVelocity.z*math.cos(phiRads)+currentVelocity.x*math.sin(phiRads)))/const +
                      currentVelocity.y*math.cos(thetaRads)
      val z: Double = -const*math.sin(thetaRads)*cos(phiRads) + currentVelocity.z*math.cos(thetaRads)
      currentVelocity.update(x, y, z)
    }
  }

  def noCollisionUpdate(dx: Double) {
    // transport particle in straight line
    val velMag: Double = getSpeed

   // update velocity v1 = v0/|v|
   currentVelocity.update(currentVelocity.x/velMag,
                          currentVelocity.y/velMag,
                          currentVelocity.z/velMag)

   // update position x1 = x0 + dr*v1
   currentPosition.update(currentPosition.x + dx*currentVelocity.x,
                          currentPosition.y + dx*currentVelocity.y,
                          currentPosition.z + dx*currentVelocity.z)

  }

  // cmScatteringAngle -> center of mass scattering angle [deg]
  // targetMass -> mass of target particle
  // returns labScatteringAngle -> lab frame scattering angle [deg]
  def cmToLabAngle(cmScatteringAngle: Double, targetMass: Double): Double = {
    val ratio: Double = projecileMass/targetMass
    val tmp: Double = ( math.cos(toRadians(cmScatteringAngle)) + ratio ) / 
                      ( math.sqrt(1.0d + 2.0d*ratio*math.cos(toRadians(cmScatteringAngle)) + math.pow(ratio,2.0)) )
    toDegrees(math.acos(tmp))
  }

  def newEnergy(scatteringAngle: Double, targetMass: Double): Double = {
    val ratio: Double = projecileMass/targetMass
    var constant: Double = (1.0d + 2.0d*ratio*math.cos(toRadians(scatteringAngle)) + math.pow(ratio,2.0)) /
                           (math.pow((1.0d + ratio),2.0))
    if (constant > 1.0d) constant = 1.0d
    currentEnergy*constant
  }

  def toDegrees(rad: Double): Double = { rad*360.0d/(2.0d*math.Pi) }

  def toRadians(deg: Double): Double = { deg*2.0d*math.Pi/360.0d }

}


class Transport extends Serializable {

  def fullTransport(z: HotParticle): HotParticle = {
    var tmp: HotParticle = z
    while (!tmp.exitConditions) {
      tmp = stepTransport(tmp)
    }
    // println("This particle had " + z.numberOfCollisions.toString + " total collisions")
    tmp
  }

  /////////////////////////////////////// 
  // Transport a hot particle a single
  // step and update parameters
  /////////////////////////////////////// 
  def stepTransport(z: HotParticle): HotParticle = {

    var randy: Random = new Random()

    // get density of atmosphere at current postition [1/m^3]
    // [String, Double] -> [targName, targDensity(pos)]
    val atmosphereDensity: HashMap[String, Double] = z.atmosphere.getAtmosphereDensity(z.currentPosition.toTuple)

    // calculate total density at current postition [1/m^3]
    val totalAtmosphereDensity: Double = z.atmosphere.getTotalAtmosphereDensity(z.currentPosition.toTuple) 

    // get total cross section for projectile-atmosphere [m^2]
    // [String, Double] -> [targName, projTargTCS(pos)]
    var atmosphereTCS: HashMap[String, Double] = z.atmosphere.getTCS(z.projectileName, z.currentEnergy, z.crossSections)

    // get mean free path for projectile-atmosphere at current energy-position
    // [String, Double] -> [targName, projTargMFP(pos)]
    var atmosphereAllMFP: HashMap[String, Double] = z.atmosphere.getAllMFP(atmosphereDensity, atmosphereTCS)

    // calculate total mean free path [m]
    val atmosphereMFP: Double = z.atmosphere.getMFP(atmosphereAllMFP)

    val verbosePrinting: Boolean = false

    if (verbosePrinting) {
      println("atmosphereDensity          -> " + atmosphereDensity.toString)
      println("totalAtmosphereDensity     -> " + totalAtmosphereDensity.toString)
      println("atmosphereTCS              -> " + atmosphereTCS.toString)
      println("atmosphere all MFP         -> " + atmosphereAllMFP.toString)
      println("atmosphere MFP             -> " + atmosphereMFP.toString)
    }

    // calculate step size for transport
    // if 0.2*MFP > 1km: 1km ? 0.2*MFP
    val stepSize = if (0.2d*atmosphereMFP > 1.0e3) 1.0e3 else 0.2d*atmosphereMFP

    // draw random number for collision
    val collisionRandy: Double = randy.nextDouble()

    // calculate collision probability
    val collisionProbability: Double = math.exp(-stepSize/atmosphereMFP)

    val collisionOccurs: Boolean = {
      if (collisionRandy > collisionProbability) true
      else false
    }

    // branch for collision or not
    // collision occurs
    if (collisionOccurs) {

      // calculate collision length [m]
      val collisionLength: Double = -atmosphereMFP*math.log(collisionRandy)

      // calculate transport time [sec]
      val transportTime: Double = collisionLength/z.getSpeed

      // calculate mixing ratio for atmosphere
      val mixingRatios: HashMap[String, Double] = z.atmosphere.getMixingRatios(atmosphereDensity, totalAtmosphereDensity)

      // calculate collision target probability array
      val targetProbability: HashMap[String, Double] = z.atmosphere.getTargetProbability(atmosphereTCS, mixingRatios)

      // draw random number and get collision target
      val target: String = z.atmosphere.getTarget(targetProbability, randy.nextDouble())
      val targetMass: Double = z.atmosphere.atmosphereParticles(target)

      // calculate reduced mass for collision
      val mu: Double = targetMass*z.projecileMass/(targetMass+z.projecileMass)

      // get random angle for scattering [deg]
      val scatteringAngle: Double = z.crossSections.getScatteringAngle(z.currentEnergy, target)

      // get random phi scattring angle [deg]
      val phiScatteringAngle: Double = 360.0d*randy.nextDouble()

      // update energy after collision
      z.currentEnergy = z.newEnergy(scatteringAngle, targetMass)

      // update collision counters and nascent hot particles
      z.numberOfCollisions += 1

      // convert scattering angle 
      val scatteringAngleLabFrame: Double = z.cmToLabAngle(scatteringAngle, targetMass)

      // update energy, velocity, and position
      z.collisionUpdate(scatteringAngle, phiScatteringAngle, collisionLength)

    }
    // no collision occurs
    else {

      // calculate transport time [sec]
      val transportTime: Double = stepSize/z.getSpeed

      // update velocity and position
      z.noCollisionUpdate(stepSize)
    }
    z.numberOfClicks += 1
    z
  }
}

