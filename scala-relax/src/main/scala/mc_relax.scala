
import parameters.Parameters
import hotParticle._
import atmosphere.Atmosphere
import crossSections.CrossSections

import org.apache.spark.SparkContext
import org.apache.spark.SparkConf
import org.apache.spark.rdd.RDD

import scala.collection.mutable.HashMap

import java.io.File
import com.typesafe.config.{ Config, ConfigFactory }
import scala.io.Source

object mc_relax extends Serializable {

  //////////////////////////////////////////////////////////////////////
  // main relaxation function
  //////////////////////////////////////////////////////////////////////
  def main(args: Array[String]) {


    //////////////////////////////////////////////////////////////////////
    // read input file
    //////////////////////////////////////////////////////////////////////
    // example of how system properties override; note this
    // must be set before the config lib is used
    System.setProperty("simple-lib.whatever", "This value comes from a system property")

    // Load our own config values from the default location, application.conf
    val conf = ConfigFactory.load()
    val N_Hots: Int = conf.getInt("relax.NumberParticles")
   
    //////////////////////////////////////////////////////////////////////
    // get Spark Context
    //////////////////////////////////////////////////////////////////////
    val sc = new SparkContext(new SparkConf().setAppName("mc-relax").setMaster("local[1]"))

    val test: RDD[Int] = sc.parallelize(Array(1, 2, 3, 4))
    if (test.count() != 4) println("Simple Spark test failed!")

    //////////////////////////////////////////////////////////////////////
    // read in all paremters
    // (N_part)
    // (proj, targ, start_pos_dist, start_vel_dist)
    // (system_comp, system_dist, system_escape_params)
    //////////////////////////////////////////////////////////////////////
    val parameters = new Parameters
    parameters.read_parameters("test_filename.dat")
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
        // println("yay!! I have the cross sections I need!")
    } 
    // if cross sections need to be calculated before mc simulation
    else {
        println("I need to do some computation")
    }

    //////////////////////////////////////////////////////////////////////
    // do transport simulation until all particles have met exit condition
    //////////////////////////////////////////////////////////////////////
    var gen1HotsTrans_rdd: RDD[HotParticle] = gen1Hots_rdd.map(x => x.trans.fullTransport(x))
    var NP = gen1HotsTrans_rdd.count() 
    println("Finished transporting " + NP.toString + " particles")
    var collisionNumber: Int = gen1HotsTrans_rdd.map(x => x.numberOfCollisions).reduce(_+_)
    var clickNumber: Int = gen1HotsTrans_rdd.map(x => x.numberOfClicks).reduce(_+_)
    // println("collisionNumber -> " + collisionNumber.toString)
    // println("clickNumber     -> " + clickNumber.toString)
    println("collision probability -> "+(collisionNumber.toDouble/clickNumber.toDouble).toString)

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



}