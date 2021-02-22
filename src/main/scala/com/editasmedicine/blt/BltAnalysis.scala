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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.alignment.{Alignment, NeedlemanWunschAligner}
import com.fulcrumgenomics.commons.util.{LazyLogging, NumericCounter, SimpleCounter}
import com.fulcrumgenomics.util.{Metric => FgMetric, _}

import scala.collection.mutable
import scala.collection.mutable.{ArrayBuffer, ListBuffer}
import scala.concurrent.forkjoin.ForkJoinPool

/** Constants used in a BLT analysis. */
object BltAnalysis {
  val CutRateByMismatchRscript = "com/editasmedicine/blt/cut_rate_by_mismatches.R"
  val PerSampleRScript         = "com/editasmedicine/blt/per_sample_plots.R"
}

/**
  * Class that provides analysis methods for BLT data once the data has been extracted from the
  * fastq reads into [[BltRead]]s and demultiplexed.
  *
  * @param manifest the sample manifest that drives the analysis
  * @param minUncutReads The minimum number of uncut observations to consider a target valid
  * @param minIdenticalFraction The minimum fraction of uncut targets that are identical to consider a target valid
  * @param useCutSamplesInValidation if true, use uncut observations from cut samples along with the uncut/naive
  *                                  samples in order to do target/umi pairing and validation. If false, use only
  *                                  the uncut samples.
  * @param threads the maximum number of threads to use during the analysis
  * @param output the directory in which to create the output files
  */
class BltAnalysis(val manifest: SampleManifest,
                  val minUncutReads: Int = Defaults.MinUncutReadsPerTarget,
                  val minIdenticalFraction: Double = Defaults.MinFractionIdenticalPerTarget,
                  val useCutSamplesInValidation: Boolean = Defaults.UseCutSamplesInValidation,
                  val threads: Int,
                  val output: DirPath
                 ) extends LazyLogging {

  require(minUncutReads >= 1, "minUncutReads must be >= 1")
  require(minIdenticalFraction >= 0 && minIdenticalFraction <= 1, "minIdenticalFraction must be between 0 and 1")
  require(threads >= 1, "threads must be a positive integer.")

  private val sampleTargetMetricUmiFiles = mutable.Map[Sample, FilePath]()
  private val sampleTargetMetricFiles    = mutable.Map[Sample, FilePath]()
  private val sampleMetricFiles          = mutable.Map[Sample, FilePath]()
  private val pool = new ForkJoinPool(threads)

  /**
    * The aligner that is used to align the guide sequences to all the target sequences to
    * try and determine how many mismatches and indels there are.  The scoring parameters have
    * been picked based on the expected frequency of indels vs. mismatches in the target sequences.
    * The procedure for determining this is:
    *   1. Calculate the fraction of targets that are longer (inserted) and shorter (deleted) vs. expected length
    *   2. Use the fractions from #1 to estimate what fraction of targets at the expected length have a
    *      balancing insertion and deletion in them
    *   3. Try different values for the scoring parameters until the fraction of alignments for targets of
    *      the expected length that have indels in them is order-of-magnitude the expected fraction from #2
    */
  private val pairwiseAligner = NeedlemanWunschAligner(matchScore = 4, mismatchScore = -2, gapOpen = -5, gapExtend = -2)

  /**
    * Analyze a BLT experiment starting from demuliplexed BltReads for one or more samples.  Proceeds through:
    *   1. Collapsing raw reads into "observations" by sample & UMI
    *   2. Validating targets (finding the set of target/umi pairings that have sufficient uncut observations)
    *   3. Generating per-sample outputs of various kinds
    *   4. Generating plots
    *
    * Largely just invokes other methods in this class and then writes the results to various files!
    *
    * Note: this method will clear the input `reads` buffer once it has assembled reads into observations,
    * in order to free up memory for further processing.
    */
  def analyze(reads: mutable.Buffer[BltRead]): Unit = {
    logger.info("Building observations from raw reads.")
    val obs     = buildObservations(reads)
    logger.info(f"Created ${obs.size}%,d observations from ${reads.size}%,d reads.")
    reads.clear()  // It seems horrible to do this, but freeing the memory from the reads is really beneficial

    logger.info("Validating targets")
    val validationMetricsPath = output.resolve("target_validation.txt.gz")
    val targets               = validateTargets(obs, Some(validationMetricsPath))
    logger.info(f"Validated ${targets.size}%,d targets with ${targets.map(_.observations.length).sum}%,d observations from ${obs.size}%,d observations.")
    obs.clear()

    // Move on to per-sample stuff
    val summaries = new ListBuffer[BltMetric]
    val samplesToPlot = new ListBuffer[Sample]
    manifest.samples.foreach { sample =>
      logger.info(s"Generating sample metrics for $sample")
      val us = sampleTargetUmiMetrics(sample, targets.iterator)
      val ts = sampleTargetMetrics(us)
      val ss = sampleMetrics(sample, ts)
      summaries += bltMetric(sample, ss)

      val dir = output.resolve(sample.name)
      Io.mkdirs(dir)

      val umiOut     = dir.resolve(sample.name + ".umis.txt.gz")
      val targetOut  = dir.resolve(sample.name + ".targets.txt.gz")
      val summaryOut = dir.resolve(sample.name + ".summary.txt")
      Metric.write(umiOut,     us)
      Metric.write(targetOut,  ts)
      Metric.write(summaryOut, ss)
      this.sampleTargetMetricUmiFiles(sample) = umiOut
      this.sampleTargetMetricFiles(sample)    = targetOut
      this.sampleMetricFiles(sample)          = summaryOut

      // Only plot samples if they are cut samples with data!
      if (sample.cut) {
        if (ts.nonEmpty) samplesToPlot += sample
        else logLoudly(Seq(s"Sample ${sample.name} did not generate any usable data!"), char='-', l=logger.warning)
      }
    }

    // Write out the summary file
    Metric.write(output.resolve("summary.txt"), summaries)

    // Generate some plots
    if (samplesToPlot.nonEmpty) {
      logger.info("Generating plots via R.")
      val cutSampleMetrics = this.sampleMetricFiles.filter { case (sample, metric) => sample.cut }.values.map(_.toString).toList
      Rscript.exec(BltAnalysis.CutRateByMismatchRscript, (output.resolve("cut_rate_by_mismatches.pdf") :: cutSampleMetrics).map(_.toString): _*)

      for (sample <- samplesToPlot) {
        val out = output.resolve(sample.name).resolve(sample.name + ".pdf").toString
        val in  = this.sampleTargetMetricFiles(sample).toString
        Rscript.exec(BltAnalysis.PerSampleRScript, in, out)
      }
    }
    else {
      val lines = Seq(
        "No samples generated any usable data.  No plots will be generated.",
        "Please check sample barcodes, and PAM sequences are correct."
      )
      logLoudly(lines, l=logger.warning)
    }
  }

  /** Logs a message of one or more lines wrapped in a high visibility banner. */
  private def logLoudly(lines: Seq[String], char: Char = '#', l : String => Unit): Unit = {
    val symbol = String.valueOf(char)
    val length = lines.map(_.length).max + 4
    val banner = symbol * length
    l(banner)
    lines.foreach { line => l(s"$symbol $line".padTo(length-1, " ").mkString + symbol) }
    l(banner)
  }

  /**
    * Collapses duplicate [[BltRead]]s down into [[BltObservation]]s that share the same
    * sample, RBC and edit status.
    */
  def buildObservations(reads: Seq[BltRead]): mutable.Buffer[BltObservation] = {
    val obs    = new ArrayBuffer[BltObservation]()
    obs.sizeHint(reads.size)
    val sorted = reads.sortWith(BltRead.DuplicateOrdering.lt).iterator.bufferBetter
    while (sorted.hasNext) {
      val head  = sorted.head
      val rs    = sorted.takeWhile(BltRead.DuplicateOrdering.equiv(head,_)).toArray
      val cut   = rs.count(_.cut)

      require(cut == 0 || cut == rs.size, "Duplicate set should not have both cut and uncut sequences.")

      obs += BltObservation(
        sample          = head.sample,
        targetBarcode   = head.targetBarcode,
        targetSequences = rs.map(_.targetSequence),
        cut             = cut > 0
      )
    }

    obs
  }


  /**
    * Takes the set of all observations across all samples in an experiment and determines which targets
    * have sufficient observations of the uncut target sequence to validate a given target/umi pairing.
    *
    * Generates a set of [[TargetInfo]] objects for the valid targets, which contain all the observations
    * for those targets.
    */
  def validateTargets(observations: Seq[BltObservation], metricsPath: Option[FilePath] = None): Seq[TargetInfo] = {
    val infos   = new ArrayBuffer[TargetInfo]()
    val writer  = metricsPath.map(path => FgMetric.writer[TargetValidationMetric](path))

    observations.view
      .groupBy(o => (o.targetBarcode, o.sample.guide, o.sample.pam))
      .foreach { case ((barcode, guide, pam), obs) =>

        val uncutReads = obs.filterNot(_.cut).filter(o => this.useCutSamplesInValidation || !o.sample.cut).flatMap(_.targetSequences)

        val (target, fractionIdentical) = {
          if (uncutReads.isEmpty) {
            (None, None)
          }
          else {
            val counter = new SimpleCounter[String]()
            uncutReads.foreach(counter.count)
            val (target, count) = counter.maxBy(_._2)
            (Some(target), Some(count / counter.total.toDouble))
          }
        }

        val validated = uncutReads.size >= minUncutReads && fractionIdentical.exists(f => f >= minIdenticalFraction)

        // Only generate and write metrics if an output path was given.
        writer.foreach { w =>
          val metric = TargetValidationMetric(
            guide                        = guide,
            pam                          = pam,
            target                       = target,
            umi                          = barcode,
            valid                        = validated,
            cut_reads_in_cut_samples     = obs.filter(o =>  o.cut &&  o.sample.cut).map(_.readCount).sum,
            uncut_reads_in_cut_samples   = obs.filter(o => !o.cut &&  o.sample.cut).map(_.readCount).sum,
            cut_reads_in_naive_samples   = obs.filter(o =>  o.cut && !o.sample.cut).map(_.readCount).sum,
            uncut_reads_in_naive_samples = obs.filter(o => !o.cut && !o.sample.cut).map(_.readCount).sum,
            fraction_identical           = fractionIdentical
          )

          w.write(metric)
        }

        if (validated) {
          val targetSequence = target.getOrElse(unreachable("Target must be defined."))
          infos += TargetInfo(
            originalGuide  = guide,
            pam            = pam,
            targetSequence = targetSequence,
            targetBarcode  = barcode,
            observations   = obs.toArray,
            annotations    = annotateTarget(guide, targetSequence, obs.head.sample.enzyme)
          )
        }
    }

    writer.foreach(_.close())
    infos
  }

  /**
    * Annotates a target with information about how it matches/mismatches the guide sequence
    * @param guide the guide sequence
    * @param target the target sequence
    * @return a [[TargetAnnotation]] object with information about indels and mismatches
    */
  def annotateTarget(guide: String, target: String, enzyme: Enzyme): TargetAnnotation = {
    val indelBases        = target.length - guide.length
    val mismatchPositions = mutable.ArrayBuilder.make[Byte]()
    val alignment         = this.pairwiseAligner.align(guide, target)
    TargetAnnotation(alignment.cigar, pamIs5PrimeOfTarget=enzyme.pamIs5PrimeOfTarget)
  }


  /**
    * Generates the [[SampleTargetMetric]] objects for a sample.
    * @param sample the sample being analyzed
    * @param targets the set of targets for the sample
    */
  def sampleTargetUmiMetrics(sample: Sample, targets: Iterator[TargetInfo]): mutable.Buffer[SampleTargetMetric] = {
    val metrics = targets
      .filter(_.observations.exists(_.sample == sample))
      .map { t =>

        val obs   = t.observations.filter(_.sample == sample)
        val cut   = obs.count(_.cut)
        val uncut = obs.length - cut
        val Seq(paddedGuide, alignment, paddedTarget) = Alignment(t.originalGuide, t.targetSequence, 1, 1, t.annotations.cigar, 1).paddedString()

        SampleTargetMetric(
          sample                 = sample.name,
          guide                  = t.originalGuide,
          target                 = t.targetSequence,
          umi                    = t.targetBarcode,
          cigar                  = t.annotations.cigar.toString(),
          indel_bases            = t.annotations.indelBases,
          mismatches             = t.annotations.mismatches,
          mean_mismatch_position = t.annotations.meanMismatchPosition,
          mismatch_tuples        = t.mismatchTuples.mkString("[", ",", "]"),
          genomic_location       = sample.offTarget.getOrElse(t.targetSequence, ""),
          obs_uncut              = uncut,
          obs_cut                = cut,
          obs_total              = uncut + cut,
          cut_rate               = cut / (cut + uncut).toDouble,
          normalized_cut_rate    = 0, // to be filled in by [[normalize()]]
          padded_guide           = paddedGuide,
          alignment              = alignment,
          padded_target          = paddedTarget
        ).withAttributes(sample.extendedAttributes)
      }.toBuffer

    // Fill in the normalized cut rate, normalizing by the cut rate where there are no mismatches or indels
    normalize(metrics)
    metrics
  }

  /** Collapses observations across UMIs for the same sample/target. */
  def sampleTargetMetrics(in: Seq[SampleTargetMetric]): mutable.Buffer[SampleTargetMetric] = {
    val metrics = in.groupBy(_.target).values.map { xs => xs.reduce( (a, b) => a.mergedWith(b) ) }.toBuffer
    normalize(metrics)
    metrics
  }

  /** Re-calculates the normalized cut rate and CI for all metrics. */
  def normalize(metrics: Seq[SampleTargetMetric]): Unit = {
    val xs    = metrics.filter(m => m.indel_bases == 0 && m.mismatches == 0)
    val cut   = xs.map(_.obs_cut).sum
    val uncut = xs.map(_.obs_uncut).sum
    val rate  = cut / (cut + uncut).toDouble
    metrics.foreach(m => m.normalized(rate))
  }

  /** Constructs sample-level metric objects. */
  def sampleMetrics(sample: Sample, targetMetrics: Seq[SampleTargetMetric]): Seq[SampleMetric] = {
    type Mismatches = Int
    val targets  = new NumericCounter[Mismatches]
    val cutObs   = new NumericCounter[Mismatches]
    val uncutObs = new NumericCounter[Mismatches]

    targetMetrics.filter(_.indel_bases == 0).foreach { m =>
      targets.count(m.mismatches)
      cutObs.count(m.mismatches, m.obs_cut)
      uncutObs.count(m.mismatches, m.obs_uncut)
    }

    val maxMismatches = if (targets.nonEmpty) targets.map(_._1).max else -1
    val zeroMmCutRate = if (cutObs.countOf(0) == 0) 1 else cutObs.countOf(0) / (cutObs.countOf(0) + uncutObs.countOf(0)).toDouble

    Range.inclusive(0, maxMismatches).map { mismatches =>
      val obs     = cutObs.countOf(mismatches) + uncutObs.countOf(mismatches)
      val cutRate = if (obs == 0) 0.0 else cutObs.countOf(mismatches) / obs.toDouble

      SampleMetric(sample.name, mismatches=mismatches, targets=targets.countOf(mismatches).toInt, observations=obs.toInt,
        cut_rate=cutRate, normalized_cut_rate=cutRate / zeroMmCutRate).withAttributes(sample.extendedAttributes)
    }
  }

  /**
    * Generate the one-row-per-sample metric with the GIMP score.  The score is calculated as the area
    * under the (mismatch x cut-rate) curve for cut-rate in [1..n].
    *
    * @param sample the sample for which the metric is being calculated
    * @param sampleMetrics the set of sample metrics from which to calculate the score
    * @param n the upper bound of number of mismatches for calculating the area
    */
  def bltMetric(sample: Sample, sampleMetrics: Seq[SampleMetric], n: Int = 4): BltMetric = {
    val score = sampleMetrics
      .filter(m => m.mismatches >= 1 && m.mismatches <= n)
      .sortBy(_.mismatches)
      .sliding(2)
      .map { case Seq(a, b) => (a.normalized_cut_rate + b.normalized_cut_rate) / 2.0 }
      .sum / (n - 1)

    BltMetric(sample=sample.name, guide=sample.guide, enzyme=sample.enzyme.toString, pam=sample.pam, score=score).withAttributes(sample.extendedAttributes)
  }
}
