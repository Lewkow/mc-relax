package hotParticle

import atmosphere.Atmosphere
import crossSections.CrossSections

class HotParticle extends Serializable {

  import scala.util.Random
  import scala.math._
  import scala.collection.mutable._

  val verbosePrinting: Boolean = false

  var projectileName: String = ""
  var projecileMass: Double = 0.0
  var currentPosition: (Double, Double, Double) = (0.0, 0.0, 0.0)
  var currentVelocity: (Double, Double, Double) = (0.0, 0.0, 0.0)
  var currentEnergy: Double = 0.0
  var currentTime: Double = 0.0
  var numberOfCollisions: Int = 0
  var numberOfClicks: Int = 0
  var generation: Int = 0
  var randy: Random = new Random()
  var atmosphere: Atmosphere = new Atmosphere
  var crossSections: CrossSections = new CrossSections


  def setParameters(name: String, 
                    mass: Double, 
                    initPosition: (Double, Double, Double),
                    intitVelocity: (Double, Double, Double),
                    particleGeneration: Int,
                    inputAtmosphere: Atmosphere) {
    projectileName = name
    projecileMass = mass
    currentPosition = initPosition
    currentVelocity = intitVelocity
    currentEnergy = getEnergy
    generation = particleGeneration
    atmosphere = inputAtmosphere
    // maximum energy for cross sections
    val max_energy: Double = 10000.0d
    crossSections.setParameters(projectileName,atmosphere,max_energy)
  }

  def getEnergy: Double = {
    math.sqrt(currentVelocity._1*currentVelocity._1 +
              currentVelocity._2*currentVelocity._2 + 
              currentVelocity._3*currentVelocity._3)*0.5d*projecileMass
  }

  def getSpeed: Double = {
    math.sqrt(currentVelocity._1*currentVelocity._1 +
              currentVelocity._2*currentVelocity._2 +
              currentVelocity._3*currentVelocity._3)
  }

  def fullTransport {
    while (!exitConditions) {
      transport
    }
  }

  def getCollisionProbability: Double = 100.0d*numberOfCollisions.toDouble/numberOfClicks.toDouble

  def printCollisionProbability {
    val p: Double = 100*numberOfCollisions.toDouble/numberOfClicks.toDouble
    println(p.toString + " collisions per click")
  }
 
  def transport {

    // get density of atmosphere at current postition [1/m^3]
    val atmosphereDensity: HashMap[String, Double] = atmosphere.getAtmosphereDensity(currentPosition)

    // calculate total density at current postition [1/m^3]
    val totalAtmosphereDensity: Double = atmosphere.getTotalAtmosphereDensity(currentPosition) 

    // get total cross section for projectile-atmosphere [m^2]
    var atmosphereTCS: HashMap[String, Double] = atmosphere.getTCS(projectileName,currentEnergy, crossSections)

    // get mean free path for projectile-atmosphere at current energy-position
    var atmosphereAllMFP: HashMap[String, Double] = atmosphere.getAllMFP(atmosphereDensity, atmosphereTCS)

    // calculate total mean free path [m]
    val atmosphereMFP: Double = atmosphere.getMFP(atmosphereAllMFP)

    if (verbosePrinting) {
      println("atmosphereDensity          -> " + atmosphereDensity.toString)
      println("totalAtmosphereDensity     -> " + totalAtmosphereDensity.toString)
      println("atmosphereTCS              -> " + atmosphereTCS.toString)
      println("atmosphere all MFP         -> " + atmosphereAllMFP.toString)
      println("atmosphere MFP             -> " + atmosphereMFP.toString)
    }

    // calculate step size for transport
    // if 0.2*MFP > 1km: 1km ? 0.2*MFP
    val stepSize = {
      if (0.2d*atmosphereMFP > 1.0e3) {
        1.0e3
      } else {
        0.2d*atmosphereMFP
      }
    }

    // draw random number for collision
    val collisionRandy: Double = randy.nextDouble()

    // calculate collision probability
    val collisionProbability: Double = exp(-stepSize/atmosphereMFP)

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
      val transportTime: Double = collisionLength/getSpeed

      // calculate mixing ratio for atmosphere

      // calculate collision probability array

      // draw random number for collision target
      val targetRandy: Double = randy.nextDouble()

      // calculate reduced mass for collision

      // get random angle for scattering

      // update energy after collision

      // update collision counters and nascent hot particles
      numberOfCollisions += 1

      // convert scattering angle 

      // update energy, velocity, and position

    }
    // no collision occurs
    else {

      // calculate transport time [sec]
      val transportTime: Double = stepSize/getSpeed

      // transport particle in straight line

      // update position and velocity 

    }
    numberOfClicks += 1

  }

  def exitConditions: Boolean = {
    if (numberOfCollisions > 100) true
    else false
  }

}