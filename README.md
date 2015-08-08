# mc-relax
Monte Carlo atomic and molecular particle simulator: Simulating the relaxation of particles, one collision at a time!

## Class specs

### class HotParticle
* projectileName: String
  * Name of projectile particle (He, O, CO2, H2O, etc)
* projectileMass: Double
  * Mass of the projectile [amu]
* currentAtmosphere: Atmosphere
  * Object for current atmosphere at CurrentPosition
  * Atmosphere particles, densities, temperature
* currentPosition: (Double, Double, Double)
  * Current projectile position (x, y, z) [km, km, km]
* currentVelocity: (Double, Double, Double)
  * Current projectile velocity (vx, vy, vz) [km/sec, km/sec, km/sec]
* currentEnergy: Double
  * Current projectile energy [eV]
* currentTime: Double
  * Current time since hot particle creation [sec]
* numberOfCollisions: Int
  * Collision counter
* numberOfClicks: Int
   * Click counter
* transport
  * Transport a HotParticle object one click and return (projectile, secondary) tuple which contains the updated projectile, after the transport, and the particle which was collided with if a collision occured. If no collision occured, secondary will have a currentEnergy = 0. 

### class Atmosphere
* atmosphereParticles: Map[String: Double]
  * [ParticleName -> Density[m^-3]]
* temperature
  * Atmosphere temperature [K]
* Double = temperatureEV
  * Returns temperature in eV

