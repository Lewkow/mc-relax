
import parameters.Parameters

object mc_relax extends Serializable {

  // main relaxation function
  def main(args: Array[String]) {
    val parameters = read_inputs("test_filename.dat")
  }

  def read_inputs(filename: String): Parameters = {
    val par = new Parameters
    par.projectile = "He"
    par
  }


}