package atmosphere

import scala.collection.mutable._

class Atmosphere extends Serializable {

  // atmospheric particles and densities Map[particleName, particleMass]
  var atmosphereParticles = HashMap.empty[String, Double]

  // atmospheric temperature [K]
  var temperature: Double = 0.0

  // init the Atmosphere class
  def setParameters(particles: HashMap[String, Double],
                    tmp: Double) = {
    val particleNames = particles.keysIterator
    while (particleNames.hasNext) {
      val nameNow: String = particleNames.next
      atmosphereParticles += (nameNow -> particles(nameNow))
    }
    temperature = tmp
  }

  // get density of atmosphere at position (x, y, z) [m, m, m]
  def getAtmosphereDensity(position: (Double, Double, Double)): HashMap[String, Double] = {
    val atmosphereDensities = HashMap.empty[String, Double]
    atmosphereParticles foreach {case (key, value) => atmosphereDensities += (key -> 1.0e-15)}
    atmosphereParticles
  }

  // convert temperature from K to eV
  def temperatureEV: Double = {
    // convert temperature from K to eV
    temperature
  }



}