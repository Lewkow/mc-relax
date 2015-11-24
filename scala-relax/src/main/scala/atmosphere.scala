package atmosphere

import scala.collection.mutable._
import scala.collection.immutable.ListMap
import crossSections.CrossSections

class Atmosphere extends Serializable {

  // atmospheric particles and densities Map[particleName, particleMass]
  var atmosphereParticles = HashMap.empty[String, Double]

  // atmospheric temperature [K]
  var temperature: Double = 0.0

  // init the Atmosphere class
  def setParameters(particles: HashMap[String, Double],
                    temp: Double) {
    val particleNames = particles.keysIterator
    while (particleNames.hasNext) {
      val nameNow: String = particleNames.next
      atmosphereParticles += (nameNow -> particles(nameNow))
    }
    temperature = temp
  }

  def setMarsParameters {
    val particles: HashMap[String, Double] = HashMap("H" -> 1.0e12,
                                                     "He" -> 1.0e11,
                                                     "O" -> 1.0e13,
                                                     "Ar" -> 1.0e11,
                                                     "H2" -> 1.0e12,
                                                     "N2" -> 1.0e13,
                                                     "CO" -> 5.0e12,
                                                     "CO2" -> 1.0e13)  
    val temp: Double = 100.0
    setParameters(particles,temp)
  }

  def isMolecule(particle: String): Boolean = {
    val molicules = ListMap("H"   -> false,
                            "H2"  -> true,
                            "N2"  -> true,
                            "Ar"  -> false,
                            "He"  -> false,
                            "O"   -> false,
                            "CO2" -> true,
                            "CO"  -> true)
    molicules(particle)
  }

  def getParticleMass(particle: String): Double = {
    val masses = ListMap("H"   -> 1.00782503207,
                         "H2"  -> 2.0157,
                         "N2"  -> 28.0134,
                         "Ar"  -> 39.948,
                         "He"  -> 4.00260325415,
                         "O"   -> 15.99491461956,
                         "CO2" -> 44.00964,
                         "CO"  -> 28.0147)
    masses(particle)
  }

  def getAnyReducedMass(proj: String, targ: String): Double = {
    val pMass: Double = getParticleMass(proj)
    val tMass: Double = getParticleMass(targ)
    (pMass*tMass)/(pMass+tMass)
  }

  // get density of atmosphere at position (x, y, z) [m, m, m]
  def getAtmosphereDensity(position: (Double, Double, Double)): HashMap[String, Double] = {
    val dummyDensity: Double = 1.0e12
    var atmosphereDensities = HashMap.empty[String, Double]
    var counter: Int = 1
    atmosphereParticles foreach {case (key, value) => atmosphereDensities += (key -> counter*dummyDensity); counter += 1}
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

  def getTargetProbability(tcs: HashMap[String, Double], mixingRatios: HashMap[String, Double]): HashMap[String, Double] = {
    var probs = HashMap.empty[String, Double]
    tcs foreach {case (key, value) => probs += (key -> value*mixingRatios(key))} 
    var totProbs: Double = 0.0
    probs foreach {case (key, value) => totProbs += value}
    var targetProbs = HashMap.empty[String, Double]
    for (i <- 0 until probs.size) {
      val key = probs.keys.toList(i)
      var probNow: Double = 0.0
      for (j <- 0 to i) {
        probNow += probs(key)
      }
      targetProbs += (key -> probNow/totProbs)
    }
    targetProbs
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

  def getTarget(probs: HashMap[String, Double], randy: Double): String = {
    var target: String = ""

    for (j <- 0 until probs.size) {
      if (j == 0) {
        if (randy <= probs(probs.keys.toList(j))) target = probs.keys.toList(j)
      }
      else if (j == probs.size-1) {
        if (randy >= probs(probs.keys.toList(j-1))) target = probs.keys.toList(j) 
      }
      else {
        if (randy > probs(probs.keys.toList(j-1)) && randy < probs(probs.keys.toList(j))) target = probs.keys.toList(j)
      }
    }
    target
  }

  // convert temperature from K to eV
  def temperatureEV: Double = {
    // convert temperature from K to eV
    temperature
  }



}