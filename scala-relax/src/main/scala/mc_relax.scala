
import parameters.Parameters
import hotParticle._
import atmosphere.Atmosphere
import crossSections.CrossSections

import org.apache.spark.SparkContext
import org.apache.spark.SparkConf
import org.apache.spark.rdd.RDD

import scala.collection.mutable._

object mc_relax extends Serializable {

  //////////////////////////////////////////////////////////////////////
  // main relaxation function
  //////////////////////////////////////////////////////////////////////
  def main(args: Array[String]) {

    //////////////////////////////////////////////////////////////////////
    // get Spark Context
    //////////////////////////////////////////////////////////////////////
    val conf = new SparkConf().setAppName("mc-relax").setMaster("local[*]")
    val sc = new SparkContext(conf)

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
    val N_Hots: Int = 100
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
    val gen1HotsTrans_rdd = gen1Hots_rdd.map(_.fullTransport)
    println("Finished transporting " + gen1HotsTrans_rdd.count().toString + " particles")
    // val gen1HotsTrans_rdd: RDD[HotParticle] = gen1Hots_rdd.map(_.fullTransport)
    // val collisionProb_rdd: RDD[Double] = gen1HotsTrans_rdd.map(_.getCollisionProbability)
    // println("Particles collided "+collisionProb.toString+" % of the time")

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