  name := "mc_relax"

  version := "1.0"

  scalaVersion := "2.11.7"

  val novusRepo = "Novus Release Repository" at "http://repo.novus.com/releases/"
  val novusSnapsRepo = "Novus Snapshots Repository" at "http://repo.novus.com/snapshots/"
  val salat = "com.novus" %% "salat-core" % "1.9.9"

  libraryDependencies += "com.novus" %% "salat" % "1.9.9"
  libraryDependencies += "org.mongodb" %% "casbah" % "2.8.1"
  libraryDependencies += "org.apache.spark" %% "spark-core" % "1.4.1" // % "provided"
  libraryDependencies += "org.apache.spark" %% "spark-sql" % "1.4.1" // % "provided"
  libraryDependencies += "com.typesafe" % "config" % "1.3.0"
