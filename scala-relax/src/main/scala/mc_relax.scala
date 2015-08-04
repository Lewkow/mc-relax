
import parameters.Parameters
import hotParticle.HotParticle
import transport.Transport
import atmosphere.Atmosphere
import crossSections.CrossSections

import org.apache.spark.SparkContext
import org.apache.spark.SparkConf



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

    val test = sc.parallelize(Array(1, 2, 3, 4))
    if (test.count() != 4) println("Simple Spark test failed!")

    //////////////////////////////////////////////////////////////////////
    // read in all paremters
    // (N_part)
    // (proj, targ, start_pos_dist, start_vel_dist)
    // (system_comp, system_dist, system_escape_params)
    //////////////////////////////////////////////////////////////////////
    val parameters = new Parameters
    parameters.read_parameters("test_filename.dat")

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
        gen1Hots(i).setParameters("He", 4.0d, initPos(i), initVel(i), 1)
    }
    val gen1Hots_rdd = sc.parallelize(gen1Hots, N_partitions)

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
    var clickCounter: Int = 0
    var keepClicking: Boolean = true
    while (keepClicking) {
        // transport

        clickCounter += 1
        if (clickCounter > 10) {keepClicking = false}
    }

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