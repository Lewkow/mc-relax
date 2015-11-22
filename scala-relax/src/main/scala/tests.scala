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
      println("Reduced Mass: "+proj+
              " + "+targ+" -> "+
              atmosphere.getAnyReducedMass(proj,targ).toString + 
              " lambda: "+CS.getAnyLambda(proj,targ))
    }
    for (p <- targets) {
      for (t <- targets) {
        printRM(p,t)
      }
    }
  }

  def assert_tests {
    assert_test_reducedMass
  }

  def assert_test_reducedMass {
    // H+H
    val x1 = ("H","H")
    val t1 = (0.503912516035, 1.0)
    val T1 = (atmosphere.getAnyReducedMass(x1._1,x1._2), CS.getAnyLambda(x1._1,x1._2))
    // He+O
    val x2 = ("He","O")
    val t2 = (3.2014621869781505, 1.0)
    val T2 = (atmosphere.getAnyReducedMass(x2._1,x2._2), CS.getAnyLambda(x2._1,x2._2))
    // Ar+CO2
    val x3 = ("Ar","CO2")
    val t3 = (20.94028725342923, 1.4)
    val T3 = (atmosphere.getAnyReducedMass(x3._1,x3._2), CS.getAnyLambda(x3._1,x3._2))
    // N2+CO
    val x4 = ("N2","CO")
    val t4 = (14.00702499245914, 1.4)
    val T4 = (atmosphere.getAnyReducedMass(x4._1,x4._2), CS.getAnyLambda(x4._1,x4._2))

    if ((t1 == T1) && (t2 == T2) && (t3 == T3) && (t4 == T4)) {
      println("Reduced Mass/Lambda test PASSED")
    } else {
      println("Reduced Mass/Lambda test FAILED")
    }
  }

}