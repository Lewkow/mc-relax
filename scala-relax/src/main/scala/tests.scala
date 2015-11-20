package tests

import scala.collection.mutable.HashMap
import crossSections.CrossSections
import atmosphere.Atmosphere

class Tests extends Serializable {

  var atmosphere = new Atmosphere
  atmosphere.setParameters(HashMap[String, Double]("CO2" -> 44.0d, 
                                                   "O" -> 16.0d, 
                                                   "He" -> 4.0d), 
                           100.0d)
  var maxEnergy: Double = 10000.0
  var projectileName: String = "H"
  var CS: CrossSections = new CrossSections
  CS.setParameters(projectileName,atmosphere,maxEnergy)

  var targets: List[String] = List("H", "He", "O", "Ar", "H2", "N2", "CO", "CO2")

  def test_reducedMass {
    def printRM(proj: String, targ: String) {
      println("Reduced Mass: "+proj+" + "+targ+" -> "+atmosphere.getAnyReducedMass(proj,targ).toString)
    }
    for (p <- targets) {
      for (t <- targets) {
        printRM(p,t)
      }
    }
  }

}