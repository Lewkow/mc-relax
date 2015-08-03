package parameters

class Parameters {
  var projectile: String = "H"

  def read_parameters(filename: String) {
    projectile = "He"
  }

  def getInitialPositions(N_particles: Int): Array[(Double, Double, Double)] = {
    var initPos: Array[(Double, Double, Double)] = new Array[(Double, Double, Double)](N_particles)
    for (i <- 0 until initPos.length) {
      initPos(i) = (100.0d, 0.0d, 0.0d)
    }
    initPos
  }

  def getInitialVelocities(N_particles: Int): Array[(Double, Double, Double)] = {
    var initVel: Array[(Double, Double, Double)] = new Array[(Double, Double, Double)](N_particles)
    for (i <- 0 until initVel.length) {
      initVel(i) = (-100.0d, 0.0d, 0.0d)
    }
    initVel 
  }

}