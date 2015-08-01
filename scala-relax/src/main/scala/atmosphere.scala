package atmosphere

import scala.collection.mutable._

class Atmosphere {

  // atmospheric particles and densities Map[particleName, particleMass]
  val atmosphereParticles = HashMap.empty[String,Double]

  // atmospheric temperature [K]
  var temperature: Double = 0.0

  // init the Atmosphere class
  def setParameters(particles: HashMap[String,Double],
                    tmp: Double) = {
    val particleNames = particles.keysIterator
    while (particleNames.hasNext) {
      val nameNow: String = particleNames.next
      atmosphereParticles += (nameNow -> particles(nameNow))
    }
    temperature = tmp
  }

  // convert temperature from K to eV
  def temperatureEV: Double = {
    // convert temperature from K to eV
    temperature
  }



}