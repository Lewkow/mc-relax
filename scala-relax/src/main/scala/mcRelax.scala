
import parameters.Parameters
import hotParticle._
import atmosphere.Atmosphere
import crossSections.CrossSections
import tests.Tests

import org.apache.spark.SparkContext
import org.apache.spark.SparkConf
import org.apache.spark.rdd.RDD

import scala.collection.mutable.HashMap

import java.io._
import com.typesafe.config.{ Config, ConfigFactory }
import scala.io.Source

object mcRelax extends Serializable {

  //////////////////////////////////////////////////////////////////////
  // main relaxation function
  //////////////////////////////////////////////////////////////////////
  def main(args: Array[String]) {
    if (args contains "-test") {
      println("\nRunning Tests\n")
      val test = new Tests
    } else if (args contains "-ang") {
      average_angle(args)
    } else if (args contains "-dcs") {
      write_dcs(args)
    } else if (args contains "-tcs") {
      test_tcs(args)
    } else if (args contains "-dcs3d") {
      write_dcs3d(args)
    } else {
      println("\nRunning Monte Carlo Simulation")
      print_stupid_graphic
      mc(args)
    }
  }

  //////////////////////////////////////////////////////////////////////
  // main monte-carlo function
  //////////////////////////////////////////////////////////////////////
  def mc(args: Array[String]) {
    //////////////////////////////////////////////////////////////////////
    // read input file
    // Load our own config values from the default location, application.conf
    //////////////////////////////////////////////////////////////////////
    val conf = ConfigFactory.load()
    val N_Hots: Int = conf.getInt("relax.NumberParticles")

    //////////////////////////////////////////////////////////////////////
    // get Spark Context
    //////////////////////////////////////////////////////////////////////
    val sc = new SparkContext(new SparkConf().setAppName("mc-relax").setMaster("local[*]"))

    val test: RDD[Int] = sc.parallelize(Array(1, 2, 3, 4))
    if (test.count() != 4) println("Simple Spark test failed!")

    //////////////////////////////////////////////////////////////////////
    // read in all paremters
    // (N_part)
    // (proj, targ, start_pos_dist, start_vel_dist)
    // (system_comp, system_dist, system_escape_params)
    //////////////////////////////////////////////////////////////////////
    val parameters = new Parameters
    parameters.readParameters("test_filename.dat")
    val atmosphere: Atmosphere = new Atmosphere
    atmosphere.setParameters(HashMap[String, Double]("CO2" -> 44.0d, "O" -> 16.0d, "He" -> 4.0d), 100.0d)

    //////////////////////////////////////////////////////////////////////
    // set initial parameters
    //////////////////////////////////////////////////////////////////////
    // Rdd of initial hot particles
    val N_partitions: Int = 4
    val initPos: Array[(Double, Double, Double)] = parameters.getInitialPositions(N_Hots)
    val initVel: Array[(Double, Double, Double)] = parameters.getInitialVelocities(N_Hots)
    var gen1Hots: Array[HotParticle] = new Array[HotParticle](N_Hots)
    for (i <- 0 until gen1Hots.length) {
        gen1Hots(i) = new HotParticle
        gen1Hots(i).setParameters("He", 4.0d, initPos(i), initVel(i), 1, atmosphere)
    }
    val gen1Hots_rdd: RDD[HotParticle] = sc.parallelize(gen1Hots, N_partitions)

    println("-- " + gen1Hots_rdd.count().toString + " initial hot particles")

    //////////////////////////////////////////////////////////////////////
    // read all scattering parameters from database for simulation
    // (tcs, scattering_ang_probs)
    //////////////////////////////////////////////////////////////////////
    // if the required cross sections exist in database
    val crossSections = new CrossSections
    if (crossSections.haveCrossSections) {
        println("yay!! I have the cross sections I need!")
    }
    // if cross sections need to be calculated before mc simulation
    else {
        println("I need to do some computation")
    }

    //////////////////////////////////////////////////////////////////////
    // do transport simulation until all particles have met exit condition
    //////////////////////////////////////////////////////////////////////
    var gen1HotsTrans_rdd: RDD[HotParticle] = gen1Hots_rdd.map(x => x.trans.fullTransport(x))
    gen1HotsTrans_rdd.cache()
    var NP = gen1HotsTrans_rdd.count()
    println("Finished transporting " + NP.toString + " particles")

    // TESTS
    test_collisionCounts(gen1HotsTrans_rdd)

    // PRINT RESULTS TO SCREEN
    var collisionNumber: Int = gen1HotsTrans_rdd.map(x => x.numberOfCollisions).reduce(_+_)
    var clickNumber: Int = gen1HotsTrans_rdd.map(x => x.numberOfClicks).reduce(_+_)
    var averageScatteringAngle: Double = gen1HotsTrans_rdd.
                                           map(x => x.totScattingAngle/x.numberOfCollisions.toDouble).
                                           reduce(_+_)/NP.toDouble
    println("collisionNumber            -> " + collisionNumber.toString)
    println("clickNumber                -> " + clickNumber.toString)
    println("ave scattering angle [deg] -> " + averageScatteringAngle.toString)
    println("collision probability      -> " + (collisionNumber.toDouble/clickNumber.toDouble).toString)

    //////////////////////////////////////////////////////////////////////
    // generate statistical distributions from simulation results
    // print overall statistics to screen
    // persist all distributions to files
    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
    // stop spark
    //////////////////////////////////////////////////////////////////////
    sc.stop()
  }

  def test_collisionCounts(hotTrans: RDD[HotParticle]) {
    val allCollisionNumbers: Array[Int] = hotTrans.map(x => x.numberOfCollisions).collect()
    val collisionNumber: Int = hotTrans.map(x => x.numberOfCollisions).reduce(_+_)
    if (allCollisionNumbers.sum == collisionNumber) {
        println("Collision number test PASSSED")
    } else { println("Collision number test FAILED") }
  }

  def test_tcs(args: Array[String]) {
    val projectile: String = args(1)
    val target: String = args(2)
    val CS = new CrossSections
    val energy = 10.0 to 1000.0 by 10.0
    val file = new File("./data/universal_tcs_"+projectile+"-"+target+".txt")
    val bw = new BufferedWriter(new FileWriter(file))
    bw.write("energy,tcs\n")
    for (en <- energy) {
      val tcs = CS.anyUniversalTotalCrossSection(en, projectile, target)
      bw.write(en+","+tcs+"\n")
    }
    bw.close()
  }

  def average_angle(args: Array[String]) {
    val projectile: String = args(1)
    val target: String = args(2)
    val CS = new CrossSections
    val energy = 10.0 to 1000.0 by 50.0
    val trials: Int = 1000
    val file = new File("./data/average_scattering_angle_"+projectile+"-"+target+".txt")
    val bw = new BufferedWriter(new FileWriter(file))
    bw.write("energy,angle,dE\n")
    for (en <- energy) {
      var tot: Double = 0
      for (i <- 1 to trials) {
        tot += CS.anyUniversalScatteringAngle(en, projectile, target)
      }
      val ang: Double = tot/trials.toDouble
      val dE: Double = CS.getAverageEnergyLoss(en, projectile, target, ang)
      bw.write(en+","+ang+","+dE+"\n")
      println(en+" eV done, ave: "+ang+" dE: "+dE)
    }
    bw.close()
  }

  def write_dcs(args: Array[String]) {
    val CS = new CrossSections
    val projectile: String = args(1)
    val target: String = args(2)
    val energy: Double = args(3).toDouble
    val file = new File("./data/universal_dcs_"+projectile+"-"+target+".txt")
    val bw = new BufferedWriter(new FileWriter(file))
    bw.write("angle,dcs\n")
    val theta = 0.01 to 170.0 by 1.0
    bw.close()
  }

  def write_dcs3d(args: Array[String]) {
    val CS = new CrossSections
    val projectile: String = args(1)
    val target: String = args(2)
    val theta = 0.01 to 10.0 by 0.1
    val energy = 100.0 to 10000.0 by 100.0
    val file = new File("./data/dcs3d_"+projectile+"-"+target+".txt")
    val bw = new BufferedWriter(new FileWriter(file))
    bw.write("energy,angle,dcs\n")
    for (en <- energy) {
      for (t <- theta) {
        val amp = CS.anyUniversalAmplitude(en, t, projectile, target)
        bw.write(en+","+t+","+amp+"\n")
      }
    }
    bw.close()
  }

  def print_stupid_graphic {
    println("|    |    / |   /   .    /    /   .  / ")
    println("   . x   |   .    .  /   /      /      ")
    println("   /     .   |            |         .  ")
    println("  /  .     /        .  .    /      x   ")
    println(".      x           . x   |     .       ")
    println("~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~")
    println("~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~")
    println("                                       ")
  }

}
