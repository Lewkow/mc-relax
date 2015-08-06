package hotParticle

import atmosphere.Atmosphere
import crossSections.CrossSections

class HotVector extends Serializable {
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

  import scala.util.Random
  import scala.math._
  import scala.collection.mutable._

  val verbosePrinting: Boolean = false

  var projectileName: String = ""
  var projecileMass: Double = 0.0
  var currentPosition = new HotVector
  var currentVelocity = new HotVector
  var currentEnergy: Double = 0.0
  var currentTime: Double = 0.0
  var numberOfCollisions: Int = 0
  var numberOfClicks: Int = 0
  var collisionProbability: Double = 0.0
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

  def fullTransport {
    while (!exitConditions) {
      transport
    }
    getCollisionProbability
  }

  def getCollisionProbability { 
    collisionProbability = 100.0d*numberOfCollisions.toDouble/numberOfClicks.toDouble
  }

  def printCollisionProbability {
    val p: Double = 100*numberOfCollisions.toDouble/numberOfClicks.toDouble
    println(p.toString + " collisions per click")
  }
 
  def transport {

    // get density of atmosphere at current postition [1/m^3]
    // [String, Double] -> [targName, targDensity(pos)]
    val atmosphereDensity: HashMap[String, Double] = atmosphere.getAtmosphereDensity(currentPosition.toTuple)

    // calculate total density at current postition [1/m^3]
    val totalAtmosphereDensity: Double = atmosphere.getTotalAtmosphereDensity(currentPosition.toTuple) 

    // get total cross section for projectile-atmosphere [m^2]
    // [String, Double] -> [targName, projTargTCS(pos)]
    var atmosphereTCS: HashMap[String, Double] = atmosphere.getTCS(projectileName,currentEnergy, crossSections)

    // get mean free path for projectile-atmosphere at current energy-position
    // [String, Double] -> [targName, projTargMFP(pos)]
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
      val mixingRatios: HashMap[String, Double] = atmosphere.getMixingRatios(atmosphereDensity, totalAtmosphereDensity)

      // calculate collision target probability array
      val targetProbability: Array[Double] = atmosphere.getTargetProbability(atmosphereTCS, mixingRatios)

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
      val velMag: Double = getSpeed

      // // update velocity v1 = v0/|v|
      currentVelocity.update(currentVelocity.x/velMag,
                             currentVelocity.y/velMag,
                             currentVelocity.z/velMag)

      // // update position x1 = x0 + dr*v1
      currentPosition.update(currentPosition.x + stepSize*currentVelocity.x,
                             currentPosition.y + stepSize*currentVelocity.y,
                             currentPosition.z + stepSize*currentVelocity.z)

    }
    numberOfClicks += 1

  }

  def exitConditions: Boolean = {
    if (numberOfCollisions > 100) true
    else false
  }

}