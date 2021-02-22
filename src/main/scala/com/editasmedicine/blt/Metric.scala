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

import java.nio.file.Path

import com.fulcrumgenomics.util.Io
import org.apache.commons.math3.stat.interval.WilsonScoreInterval

import scala.collection.immutable.SortedMap
import scala.reflect.runtime.{universe => ru}

/**
  * Base trait for all BLT metric classes to extend. Introduces a set of extended/dynamic attributes
  * that do not have to be defined at compile time.
  */
trait Metric extends com.fulcrumgenomics.util.Metric {
  var extendedAttributes: SortedMap[String, Any] = SortedMap.empty

  /** The set of all metric names, including both fixed and extended. */
  override def names: Seq[String] = super.names ++ extendedAttributes.keys.toSeq

  /** The set of all metric values, including both fixed and extended. */
  override def values: Seq[String] = super.values ++ extendedAttributes.values.map(this.formatValues).toSeq

  /** Overridden to ensure that booleans are always written as true/false. */
  override def formatValues(value: Any): String = value match {
    case b : Boolean => String.valueOf(b)
    case x           => super.formatValues(x)
  }

  /** Replaces the set of extended attributes with attributes from this map. */
  def withAttributes(attrs: Map[String,Any]): this.type = {
    attrs match {
      case sorted: SortedMap[String, Any] => this.extendedAttributes = sorted
      case map                            => this.extendedAttributes = SortedMap(map.toSeq:_*)
    }

    this
  }
}

/**
  * Used in place of fgbio's Metric object to write the metrics, in order to provide support
  * for a set of dynamic attributes in each Metric instance.
  */
object Metric {
  val Separator = "\t"

  /** Writes out a table of metrics to the Path provided. */
  def write[A <: Metric](path: Path, metrics: Seq[A])(implicit tt: ru.TypeTag[A]): Unit = {
    val out = Io.toWriter(path)

    // Write the header
    val names = if (metrics.nonEmpty) metrics.head.names else com.fulcrumgenomics.util.Metric.names[A]
    out.write(names.mkString(Separator))
    out.newLine()

    metrics.foreach { m =>
      out.write(m.values.mkString(Separator))
      out.newLine()
    }

    out.close()
  }
}

/**
  * Metrics for how the read extraction and normalization worked across all samples.
  *
  * @param total_reads the total number of reads examined (i.e. the number of reads in the run/fastqs)
  * @param reads_extracted the number of reads extracted and used in the analysis, excluding all reads accounted
  *                        for by the failure reasons given in the file
  * @param failed_to_id_landmarks the number of reads where extraction failed because landmark sequences could
  *                               not be ID'd in the sequence read with sufficient confidence
  * @param failed_to_id_sample the number of reads where extraction failed because the sample barcode did not match
  *                            any expected barcode within the given error tolerance
  * @param failed_to_extract_target the number of reads where extraction failed because the target sequences and
  *                                 UMI could not be extracted confidently
  * @param failed_quality the number of reads discarded because the extracted bases had base quality scores that
  *                       were too low
  */
case class DemuxMetric(total_reads: Long,
                       reads_extracted: Long,
                       failed_to_id_landmarks: Long,
                       failed_to_id_sample: Long,
                       failed_to_extract_target: Long,
                       failed_quality: Long
                      ) extends Metric

/**
  * Metrics about the reads attributed to each sample during extraction and de-multiplexing.
  *
  * @param sample the name of the sample
  * @param guide the guide sequence for the sample
  * @param pam the PAM sequence for the sample
  * @param total_reads the total number of reads identified as originating from the sample
  * @param reads_extracted the number of reads that were extracted and taken into analysis for the sample
  * @param failed_to_extract the number of reads where the target and/or target UMI sequence could not be extracted
  * @param failed_quality the number of reads that were discarded due to low base quality scores
  */
case class SampleDemuxMetric(sample: String,
                             guide: String,
                             pam: String,
                             total_reads: Long,
                             reads_extracted: Long,
                             failed_to_extract: Long,
                             failed_quality: Long
                            ) extends Metric

/**
  * Information about each target/UMI pair observed in the experiment, values relevant to validation of
  * the target/umi pairing and whether or not it was deemed valid.
  *
  * @param guide The guide sequence for which the target/umi observations were seen
  * @param pam The PAM sequence for which the target/umi observations were seen
  * @param target The observed target sequence. If no uncut reads were observed, or if the required fraction of
  *               agreeing observations did not agree, no sequence is stored.
  * @param umi the umi or target barcode associated with the target
  * @param cut_reads_in_cut_samples The number of cut reads observed in samples that were cut
  * @param uncut_reads_in_cut_samples The number of uncut reads observed in samples that were cut
  * @param cut_reads_in_naive_samples The number of cut reads observed in naive (i.e. uncut) samples. Should nearly always be 0.
  * @param uncut_reads_in_naive_samples The number of uncut reads observed in naive (i.e. uncut) samples.
  * @param fraction_identical The fraction of of uncut observations that were identical to the most frequently observed
  *                           uncut sequence for the target.
  */
case class TargetValidationMetric(guide:  String,
                                  pam:    String,
                                  target: Option[String],
                                  umi:    String,
                                  valid: Boolean,
                                  cut_reads_in_cut_samples: Int,
                                  uncut_reads_in_cut_samples: Int,
                                  cut_reads_in_naive_samples: Int,
                                  uncut_reads_in_naive_samples: Int,
                                  fraction_identical: Option[Double]
                                 ) extends Metric

/**
  * Metrics about how a particular target sequence performed in the assay for a given sample.
  *
  * @param sample the name of the sample
  * @param guide the guide sequence being tested
  * @param target the target sequence
  * @param umi either the UMI of the target sequence or "multiple" if observations have been combined
  * @param cigar a string describing the alignment of the guide to the target
  * @param indel_bases the number of inserted and/or deleted bases between the guide and the target
  * @param mismatches the number of mismatches between the guide and the target
  * @param mean_mismatch_position the mean position of mismatches within the target
  * @param mismatch_tuples a string of three-tuples that contain (pos,guideBase,targetBase)
  * @param genomic_location the genomic location of the target IF the target was listed in the provided off-target file
  * @param obs_uncut the number of uncut observations of the target for the sample
  * @param obs_cut the number of cut observations of the target for the sample
  * @param cut_rate the rate at which the target was cut (obs_cut / (obs_cut + obs_uncut))
  * @param normalized_cut_rate the cut_rate normalized by the cut_rate where target == guide
  */
case class SampleTargetMetric(sample:                  String,
                              guide:                   String,
                              target:                  String,
                              umi:                     String,
                              cigar:                   String,
                              indel_bases:             Int,
                              mismatches:              Int,
                              mean_mismatch_position:  Option[Double],
                              mismatch_tuples:         String,
                              genomic_location:        String,
                              obs_uncut:               Int,
                              obs_cut:                 Int,
                              obs_total:               Int,
                              cut_rate:                Double,
                              var normalized_cut_rate: Double,
                              var norm_cut_rate_ci95_low: Double  = 0,
                              var norm_cut_rate_ci95_high: Double = 0,
                              padded_guide:            String,
                              alignment:               String,
                              padded_target:           String) extends Metric {

  require(obs_total == obs_cut + obs_uncut, "obs_total doesn's match the sum of cut and uncut obs.")
  require(padded_guide.length == padded_target.length, s"padded guide $padded_guide and padded target $padded_target are different lengths.")
  require(padded_guide.length == alignment.length,     s"alignment string $alignment is different length than padded guide $padded_guide")

  /**
    * [Re-]generates the normalized cut rate and the CI given the cut rate when guide == target.
    *
    * @param rate: the cut rate for observations where the target == the guide sequence exactly
    */
  def normalized(rate: Double): this.type = {
    this.normalized_cut_rate = this.cut_rate / rate
    val ci = new WilsonScoreInterval().createInterval(obs_cut + obs_uncut, obs_cut, 0.95)
    this.norm_cut_rate_ci95_low  = ci.getLowerBound / rate
    this.norm_cut_rate_ci95_high = ci.getUpperBound / rate
    this
  }

  /** Merges two [[SampleTargetMetric]] objects that are for the same sample/target but with different UMIs. */
  def mergedWith(that: SampleTargetMetric): SampleTargetMetric = {
    val cut      = this.obs_cut   + that.obs_cut
    val uncut    = this.obs_uncut + that.obs_uncut
    val total    = cut + uncut
    val rate     = cut / total.toDouble
    copy(umi="multiple", obs_uncut=uncut, obs_cut=cut, obs_total=total, cut_rate=rate, normalized_cut_rate = -1, norm_cut_rate_ci95_low = -1, norm_cut_rate_ci95_high = -1).withAttributes(this.extendedAttributes)
  }
}

/**
  * Summary information about targets with a given number of mismatches for a single sample.
  *
  * @param sample the name of the sample
  * @param mismatches the number of mismatches between the guide and the targets
  * @param targets the number of targets observed at the mismatch distance
  * @param observations the number of observations of the targets
  * @param cut_rate the mean cut rate for targets with the given number of mismatches
  * @param normalized_cut_rate the normalized (by cut rate when target == guide) cut rate
  */
case class SampleMetric(sample: String,
                        mismatches: Int,
                        targets: Int,
                        observations: Int,
                        cut_rate: Double,
                        normalized_cut_rate: Double
                       ) extends Metric

/**
  * Sample-level summary information that gives the GIMP score per sample.
  *
  * @param sample the name of the sample
  * @param guide  the guide sequence
  * @param enzyme the enzyme used to cut
  * @param pam    the PAM sequence
  * @param score  the GIMP score
  */
case class BltMetric(sample: String, guide: String, enzyme: String, pam: String, score: Double) extends Metric
