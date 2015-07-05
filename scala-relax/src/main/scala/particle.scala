package particle

class Particle(in_name: String, 
               in_mass: Double) {

  val name: String = in_name
  val mass: Double = in_mass
  var pos: (Double, Double, Double) = (0.0, 0.0, 0.0)
  var vel: (Double, Double, Double) = (0.0, 0.0, 0.0)

  
}