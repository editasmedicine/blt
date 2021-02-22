import java.text.SimpleDateFormat
import java.util.Date
import scala.util.Try
import com.typesafe.sbt.SbtGit.GitCommand
import sbt.Keys._
import sbt._
import sbtassembly.AssemblyKeys.assembly
import scoverage.ScoverageKeys._

import scala.sys.process.Process

////////////////////////////////////////////////////////////////////////////////////////////////
// Multi-project build file for the following projects:
// - commons: utility code and base classes used across projects
// - primer-design: tools for performing primer design
// - blt: tools for analyzing "barcoded library of targets" data
// - exco: software for calculating oligo extinction coefficients
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
// Common settings for all projects
////////////////////////////////////////////////////////////////////////////////////////////////
val htmlReportsDirectory: String = "target/test-reports"

lazy val htsjdkAndPicardExcludes = Seq(
  ExclusionRule(organization="org.apache.ant"),
  ExclusionRule(organization="gov.nih.nlm.ncbi"),
  ExclusionRule(organization="org.testng"),
  ExclusionRule(organization="com.google.cloud.genomics")
)

lazy val ToolkitVersion = {
  val date     = new SimpleDateFormat("yyyyMMdd").format(new Date())
  val hash     = Try { Process("git rev-parse --short HEAD").lineStream.head }.getOrElse("unknown")
  val modified = Try { Process("git status --porcelain").lineStream.nonEmpty }.getOrElse(false)

  date + "-" + hash + (if (modified) "-dirty" else "")
}

version in ThisBuild := ToolkitVersion

lazy val commonSettings = Seq(
  organization         := "com.editasmedicine",
  organizationName     := "Editas Medicine Inc.",
  organizationHomepage := Some(url("http://www.editasmedicine.com/")),
  homepage             := Some(url("https://github.com/editasmedicine/blt")),
  startYear            := Some(2017),
  scalaVersion         := "2.12.3",
  scalacOptions        += "-target:jvm-1.8",
  autoAPIMappings      := true,
  version              := ToolkitVersion,
  testOptions in Test  += Tests.Argument(TestFrameworks.ScalaTest, "-h", Option(System.getenv("TEST_HTML_REPORTS")).getOrElse(htmlReportsDirectory)),
  testOptions in Test  += Tests.Argument("-oDF"),
  resolvers            += Resolver.jcenterRepo,
  resolvers            += Resolver.sonatypeRepo("public"),
  resolvers            += Resolver.mavenLocal,
  shellPrompt          := { state => "%s| %s> ".format(GitCommand.prompt.apply(state), version.value) },
  updateOptions        := updateOptions.value.withCachedResolution(true),
  javaOptions in Test += "-Xmx1G",
  libraryDependencies ++= Seq(
      "org.scalatest"       %% "scalatest"     % "3.0.1" % "test->*" excludeAll ExclusionRule(organization="org.junit", name="junit"),
      "com.fulcrumgenomics" %% "sopt"          % "0.3.1",
      "com.fulcrumgenomics" %% "commons"       % "0.3.0",
      "com.fulcrumgenomics" %% "fgbio"         % "0.4.0"  excludeAll(htsjdkAndPicardExcludes:_*),
      "com.beachape"        %% "enumeratum"    % "1.5.12",
      "com.github.samtools"  % "htsjdk"        % "2.11.0" excludeAll(htsjdkAndPicardExcludes: _*)
      ),
  assemblyJarName in assembly := s"${name.value}.jar"
) ++ Defaults.coreDefaultSettings

lazy val blt = Project(id="blt", base=file("."))
  .settings(description := "Tools for processing Baroded Library of Targets data.")
  .settings(commonSettings: _*)
  .settings(libraryDependencies += "org.apache.commons" % "commons-math3" % "3.6.1")
