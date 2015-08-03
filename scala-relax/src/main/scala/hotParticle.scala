package hotParticle

class HotParticle extends Serializable {

  var projectileName: String = ""
  var projecileMass: Double = 0.0
  var currentPosition: (Double, Double, Double) = (0.0, 0.0, 0.0)
  var currentVelocity: (Double, Double, Double) = (0.0, 0.0, 0.0)
  var currentEnergy: Double = 0.0
  var currentTime: Double = 0.0
  var numberOfCollisions: Int = 0
  var numberOfClicks: Int = 0
  var generation: Int = 0

  def setParameters(name: String, 
                    mass: Double, 
                    initPosition: (Double, Double, Double),
                    intitVelocity: (Double, Double, Double),
                    particleGeneration: Int) {
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
  
}