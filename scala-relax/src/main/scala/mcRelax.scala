
import parameters.Parameters
import hotParticle._
import atmosphere.Atmosphere
import crossSections.CrossSections
import tests.Tests

import org.apache.spark.SparkContext
import org.apache.spark.SparkConf
import org.apache.spark.rdd.RDD

import scala.collection.mutable.HashMap

import java.io.File
import com.typesafe.config.{ Config, ConfigFactory }
import scala.io.Source

object mcRelax extends Serializable {

  //////////////////////////////////////////////////////////////////////
  // main relaxation function
  //////////////////////////////////////////////////////////////////////
  def main(args: Array[String]) {
    if (args contains "-t") {
      println("\nRunning Tests\n")
      val test = new Tests
      test.assert_tests
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

  def print_stupid_graphic {
    println("|    |    / |   /  ")
    println("   .   |   .    .  ")
    println("   /     .   |     ")
    println("  /  .     /       ")
    println(".                 .")
    println("~-~-~-~-~-~-~-~-~-~")
    println("~-~-~-~-~-~-~-~-~-~")
    println("                   ")
  }

}