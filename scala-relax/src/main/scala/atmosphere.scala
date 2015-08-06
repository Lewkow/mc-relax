package atmosphere

import scala.collection.mutable._
import crossSections.CrossSections

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
    val dummyDensity: Double = 1.0e12
    var atmosphereDensities = HashMap.empty[String, Double]
    atmosphereParticles foreach {case (key, value) => atmosphereDensities += (key -> dummyDensity)}
    atmosphereDensities
  }

  // get total density of atmosphere at position (x, y, z) [m, m, m]
  def getTotalAtmosphereDensity(position: (Double, Double, Double)): Double = {
    val atmosphereDensities: HashMap[String, Double] = getAtmosphereDensity(position) 
    var totalDensity: Double = 0.0d
    atmosphereDensities foreach {case (key, value) => totalDensity += value}
    totalDensity
  }

  // get mixing ratios for atmosphere
  // x_density / total_density
  def getMixingRatios(densities: HashMap[String, Double], density: Double): HashMap[String, Double] = {
    val mixingRatios = HashMap.empty[String, Double]
    densities foreach {case (key, value) => mixingRatios += (key -> value/density)}
    mixingRatios
  }

  def getTCS(projectile: String, energy: Double, crossSections: CrossSections): HashMap[String, Double] = {
    var atmosphereTCS = HashMap.empty[String, Double]
    atmosphereParticles foreach {
      case (particle, mass) => atmosphereTCS += (particle -> crossSections.getTotalCrossSection(energy,particle))
    }
    atmosphereTCS
  }

  def getAllMFP(atmosphereDensities: HashMap[String, Double], 
                atmosphereTCS: HashMap[String, Double]): HashMap[String, Double] = {
    var atmosphereAllMFP = HashMap.empty[String, Double]
    atmosphereDensities foreach {
      case (particle, density) => atmosphereAllMFP += (particle -> density*atmosphereTCS(particle))
    }
    atmosphereAllMFP
  }

  def getMFP(atmosphereAllMFP: HashMap[String, Double]): Double = {
    var total: Double = 0.0
    atmosphereAllMFP foreach {
      case (particle, y) => total += y
    }
    1.0d/total
  }

  // convert temperature from K to eV
  def temperatureEV: Double = {
    // convert temperature from K to eV
    temperature
  }



}