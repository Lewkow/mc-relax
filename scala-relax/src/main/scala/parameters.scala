package parameters

class Parameters {
  var projectile: String = "H"

  def read_parameters(filename: String): Parameters = {
    val par = new Parameters
    par.projectile = "He"
    par
  }

}