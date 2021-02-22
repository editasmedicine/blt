/*
 * Copyright (c) 2017-2018 Editas Medicine, Inc.
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted (subject to the limitations in the disclaimer below) provided
 * that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or
 *   other materials provided with the distribution.
 * * Neither the name of [Owner Organization] nor the names of its contributors may
 *   be used to endorse or promote products derived from this software without
 *   specific prior written permission.
 * 
 * NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY THIS
 * LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
 * TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.editasmedicine.blt

import java.util.concurrent.atomic.AtomicInteger

import com.editasmedicine.commons.clp.ClpMain
import com.editasmedicine.commons.testing.UnitSpec
import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.fastq.{FastqRecord, FastqWriter}
import com.fulcrumgenomics.util.Io
import htsjdk.samtools.SAMUtils
import org.apache.commons.math3.distribution.BinomialDistribution

import scala.collection.mutable.ListBuffer

/**
  * This test doesn't validate the outputs extensively, but is here to ensure that the whole system
  * runs end to end without problems and produces the expected set of output files.
  */
class AnalyzeExperimentTest extends UnitSpec {
  val s1 = Sample("s1", barcode="CACATACGCACTACG", guide="GGCCTCCCCAAAGCCTGGCCA", enzyme=Cas9, "GGGAGT", cut=true)
  val s2 = Sample("s2", barcode="CCTATACCCGAATCT", guide="GGCCTCCCCAAAGCCTGGCCA", enzyme=Cas9, "GGGAGT", cut=false)
  val s3 = Sample("s3", barcode="TATACAATTCGCAGC", guide="CCCAGTGTCCCCCTTCCCTAT", enzyme=Cas9, "GGGAAT", cut=true)
  val s4 = Sample("s4", barcode="CCGGAGTAGGTCCTC", guide="CCCAGTGTCCCCCTTCCCTAT", enzyme=Cas9, "GGGAAT", cut=false)

  val manifestText =
    """
      |sample	sample_barcode	guide	pam	enzyme	off_target_file	cut
      |s1	CACATACGCACTACG	GGCCTCCCCAAAGCCTGGCCA	GGGAGT	Cas9		true
      |s2	CCTATACCCGAATCT	GGCCTCCCCAAAGCCTGGCCA	GGGAGT	Cas9		false
      |s3	TATACAATTCGCAGC	CCCAGTGTCCCCCTTCCCTAT	GGGAAT	Cas9		true
      |s4	CCGGAGTAGGTCCTC	CCCAGTGTCCCCCTTCCCTAT	GGGAAT	Cas9		false
    """.stripMargin.trim

  val Rng         = new scala.util.Random(42)
  val Bases       = Array[Char]('A', 'C', 'G', 'T')
  val FullStagger = "GCGGAAGC"
  val Serial      = new AtomicInteger

  /** Generates a random stagger between 1 and 8 bases long. */
  def stagger: String = FullStagger.take(Rng.nextInt(FullStagger.length) + 1)

  /** Generates a random sequence of length n. */
  def randomBases(length: Int):String = {
    val bases = new Array[Char](length)
    forloop(from=0, until=length) { i =>
      bases(i) = Bases(Rng.nextInt(Bases.length))
    }

    new String(bases)
  }

  /** Introduces the requested number of mismatches into the given sequence. */
  def mutate(input: String, mismatches: Int): String = {
    val bases = input.toCharArray
    Rng.shuffle(bases.indices.toList).take(mismatches).foreach { i =>
      val replacement = if (bases(i) == 'A') 'C' else 'A'
      bases(i) = replacement
    }

    new String(bases)
  }

  /** Creates a fastq record with the given bases and a quality string of all `qual`. */
  private def fq(bases: String, qual: Int = 30) = {
    // Remove any layout junk from the sequence
    val bs = bases.toUpperCase.filter(b => Sequences.isValidBase(b.toByte, allowAmbiguityCodes=true))
    FastqRecord(name="q" + Serial.incrementAndGet(), bases=bs, quals=SAMUtils.phredToFastq(qual).toString * bs.length, comment=None, readNumber=None)
  }

  /** Generates a number of fastq records for a pair of cut and uncut samples. */
  def makeTestdata(cut: Sample, uncut: Sample, targets: Int = 1000): Seq[FastqRecord] = {
    require(cut.guide == uncut.guide && cut.pam == uncut.pam)
    val out = ListBuffer[FastqRecord]()
    val mmDistribution = new BinomialDistribution(cut.guide.length, 0.1)

    for (targetNumber <- 1 to targets) {
      val mismatches = mmDistribution.sample()
      val target     = mutate(uncut.guide, mismatches)

      for (umiNumber <- 1 to Rng.nextInt(5)+1) {
        val umi = randomBases(12)

        for (s <- Seq(cut, uncut); obsCount <- 1 to Rng.nextInt(3)+1) {
          val rbc = randomBases(6)
          //             stagger   a1      rbc    a2     s1 barcode              a3                a4    target     pam   gc   umi   a6    p7 adapter
          val read  = s"${stagger}-cgatct-${rbc}-tacgac-${s.barcode}-ttaccgaagatagcagcctagtggaacc-atctg-${target}-${s.pam}-GC-${umi}-tgac-agatcggaagagcac"

          for (readCount <- 1 to Rng.nextInt(3)+1) {
            out += fq(read)
          }
        }
      }
    }

    out
  }
  
  "AnalyzeExperiment" should "run end to end and generate outputs" in {
    val dir      = makeTempDir("AnalyzeExperiment.", ".testdir")
    val fastq    = dir.resolve("in.fastq")
    val manifest = dir.resolve("manifest.txt")

    // Prepare the input fastq
    val fqRecs = makeTestdata(cut=s1, uncut=s2) ++ makeTestdata(cut=s3, uncut=s4)
    val fqOut  = FastqWriter(dir.resolve("in.fastq"))
    fqRecs.foreach( fqOut.write)
    fqOut.close()

    // Write out the manifest
    Io.writeLines(manifest, manifestText.lines.toSeq)

    // Invoke the program!!
    new ClpMain("blt-test").makeItSo(Array[String]("AnalyzeExperiment", s"--input=${fastq.toAbsolutePath}", s"--sample-manifest=${manifest.toAbsolutePath}", s"--output=${dir.toAbsolutePath}"))

    // Assert that various files got created as expected:
    dir.resolve("cut_rate_by_mismatches.pdf").toFile.exists shouldBe true
    dir.resolve("demultiplexing.details.txt").toFile.exists shouldBe true
    dir.resolve("demultiplexing.summary.txt").toFile.exists shouldBe true
    dir.resolve("summary.txt").toFile.exists                shouldBe true
    dir.resolve("target_validation.txt.gz").toFile.exists   shouldBe true

    Seq(s1, s2, s3, s4).foreach { sample =>
      val sampleDir = dir.resolve(sample.name)
      sampleDir.toFile.exists      shouldBe true
      sampleDir.toFile.isDirectory shouldBe true
      sampleDir.resolve(sample.name + ".summary.txt").toFile.exists    shouldBe true
      sampleDir.resolve(sample.name + ".targets.txt.gz").toFile.exists shouldBe true
      sampleDir.resolve(sample.name + ".umis.txt.gz").toFile.exists    shouldBe true
      if (sample.cut) sampleDir.resolve(sample.name + ".pdf").toFile.exists shouldBe true
    }

    // Read back the summary demultiplexing metrics and make sure we saw the right number of reads
    val summary = com.fulcrumgenomics.util.Metric.read[DemuxMetric](dir.resolve("demultiplexing.summary.txt"))
    summary should have size 1
    summary.head.total_reads shouldBe fqRecs.size

    // Delete the temporary directory
    deletePath(dir)
  }
}
