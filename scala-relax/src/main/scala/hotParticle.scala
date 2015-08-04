package hotParticle

import atmosphere.Atmosphere

class HotParticle extends Serializable {

  import scala.util.Random
  import scala.math._

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
    generation = particleGeneration
  }

  def getEnergy: Double = {
    math.sqrt(currentVelocity._1*currentVelocity._1 +
              currentVelocity._2*currentVelocity._2 + 
              currentVelocity._3*currentVelocity._3)*0.5d*projecileMass
  }

  def fullTransport {
    while (!exitConditions) {
      transport
    }
  }

  def printCollisionProbability {
    val p: Double = 100*numberOfCollisions.toDouble/numberOfClicks.toDouble
    println(p.toString + " collisions per click")
  }
 
  def transport {
    // get density of atmosphere at current postition [1/m^3]


    // calculate total density at current postition [1/m^3]
    
    // get total cross section for projectile-atmosphere

    // get mean free path for projectile-atmosphere at current energy-position

    // calculate total mean free path [m]
    val meanFreePath: Double = 1.0e2

    // calculate step size for transport
    // if 0.2*MFP > 1km: 1km ? 0.2*MFP
    val stepSize = {
      if (0.2d*meanFreePath > 1.0e3) {
        1.0e3
      } else {
        0.2d*meanFreePath
      }
    }

    // draw random number for collision
    val collisionRandy: Double = randy.nextDouble()

    // calculate collision probability
    val collisionProbability: Double = exp(-stepSize/meanFreePath)

    val collisionOccurs: Boolean = {
      if (collisionRandy > collisionProbability) true
      else false
    }

    // branch for collision or not
    // collision occurs
    if (collisionOccurs) {

      // calculate collision length

      // calculate transport time

      // calculate mixing ratio for atmosphere

      // calculate collision probability array

      // draw random number for collision target

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

      // calculate transport time

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