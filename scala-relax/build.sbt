  name := "mc_relax"

  version := "1.0"

  scalaVersion in ThisBuild := "2.11.7"

  libraryDependencies += "org.mongodb" %% "casbah" % "2.8.1"
  libraryDependencies += "org.apache.spark" %% "spark-core" % "1.4.1" // % "provided"
  libraryDependencies += "org.apache.spark" %% "spark-sql" % "1.4.1" // % "provided"
  libraryDependencies += "com.typesafe" % "config" % "1.3.0"
