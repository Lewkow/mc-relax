package parameters

class Parameters {

  import scala.util.Random

  var projectile: String = "H"
  var randy: Random = new Random()

  def readParameters(filename: String) {
    projectile = "He"
  }

  def getInitialPositions(N_particles: Int): Array[(Double, Double, Double)] = {
    var initPos: Array[(Double, Double, Double)] = new Array[(Double, Double, Double)](N_particles)
    for (i <- 0 until initPos.length) {
      initPos(i) = randPosition
    }
    initPos
  }

  def randPosition: (Double, Double, Double) = (randy.nextDouble*100.0d, randy.nextDouble*100.0d, randy.nextDouble*100.0d)

  def getInitialVelocities(N_particles: Int): Array[(Double, Double, Double)] = {
    var initVel: Array[(Double, Double, Double)] = new Array[(Double, Double, Double)](N_particles)
    for (i <- 0 until initVel.length) {
      initVel(i) = (100.0d - randy.nextDouble*200.0d, 
                    100.0d - randy.nextDouble*200.0d,
                    100.0d - randy.nextDouble*200.0d)
    }
    initVel 
  }

}