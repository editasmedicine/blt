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
import com.fulcrumgenomics.commons.util.{Logger, SimpleCounter}
import com.fulcrumgenomics.fastq.FastqRecord
import htsjdk.samtools.SAMUtils

/**
  * Base trait for classes that extract the relevant bits of a BLT read from a fastq record.  The implementation
  * is setup with the expectation that the sequence is broken down into two halves as follows:
  *
  *                                                 <- Same for all Enzymes | Enzyme specific setup (Cas9 shown) ->
  *   stagger[1-8]-CGATCT-rbc[6]-TACGAC-sbc[15]-TTACCGAAGATAGCAGCCTAGTGGAACC-ATCTG-target-PAM-GC-umi[12]-TGAC-AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
  *
  * The parsing of the left hand side that is the same for all enzymes is done in this trait, and the parsing of the
  * right hand side is delegated to concrete implementations.
  *
  * Most of the methods and values in this class are overrideable to allow additional flexibility in implementations.
  *
  * A number of statistics are tracked which can be logged using the `printStats()` method.
  */
sealed trait ReadExtractor {
  // Counters are by by sample name
  private var failedToIdStagger   = 0L
  private var failedToDemultiplex = 0L
  private val failedToExtract     = new SimpleCounter[String]
  private val failedQuality       = new SimpleCounter[String]
  private val extracted           = new SimpleCounter[String]

  /** Case class to hold target offset and length, and PAM offset and length. */
  case class TargetAndUmiLoc(targetOffset: Int, targetLength: Int, umiOffset: Int, umiLength: Int, cut: Boolean)

  /** Case class that represents known fixed/anchor sequences and their minimum offset in the reads. */
  case class Anchor(bases: String, offset: Int) {
    val bytes: Array[Byte] = bases.getBytes
    def length: Int = bases.length
  }

  protected val anchor1: Anchor = Anchor("CGATCT", 1)
  protected val anchor2: Anchor = Anchor("TACGAC", 13)
  protected val anchor3: Anchor = Anchor("TTACCGAAGATAGCAGCCTAGTGGAACC", 34)

  /** The set of anchors that are constant and come before the target. */
  protected val anchors = Seq(anchor1, anchor2, anchor3)

  /** The maximum number of "stagger" bases at the start of the read. */
  protected val maxStagger = 8

  /** The length of the random ligation barcode used to detect unique ligation events vs. PCR dupes. */
  protected val randomBarcodeLength: Int = 6

  /** The length of the sample barcode (aka hamming code). */
  protected val sampleBarcodeLength: Int = 15

  /** The length of the target barcode (aka the UMI). */
  protected val targetBarcodeLength: Int = 12

  /** The minimum average quality of the extracted bases. */
  val minQuality: Int

  /** Must return the demultiplexer to be used during read extraction. */
  def demultiplexer: Demultiplexer

  /** Extracts the fastq into a BltRead. */
  def extract(fq: FastqRecord): Option[BltRead] = {
    val bytes = fq.bases.getBytes

    determineStagger(fq.bases, bytes) match {
      case None =>
        this.failedToIdStagger += 1
        None
      case Some(stagger) =>
        val rbcOffset       = stagger   + anchor1.bases.length
        val sbcOffset       = rbcOffset + randomBarcodeLength + anchor2.bases.length
        val targetArmOffset = sbcOffset + sampleBarcodeLength + anchor3.bases.length

        this.demultiplexer.assign(bytes, sbcOffset) match {
          case None =>
            this.failedToDemultiplex += 1
            None
          case Some(sample) =>
            locateTargetAndUmi(fq, bytes, targetArmOffset, sample) match {
              case None =>
                this.failedToExtract.count(sample.name)
                None
              case Some(locs) =>
                val q = meanQuality(fq, (rbcOffset, randomBarcodeLength), (locs.targetOffset, locs.targetLength), (locs.umiOffset, locs.umiLength))
                if (q < this.minQuality) {
                  this.failedQuality.count(sample.name)
                  None
                }
                else {
                  this.extracted.count(sample.name)
                  Some(BltRead(
                    sample         = sample,
                    stagger        = stagger.toByte,
                    randomBarcode  = fq.bases.substring(rbcOffset, rbcOffset+randomBarcodeLength),
                    targetBarcode  = fq.bases.substring(locs.umiOffset, locs.umiOffset + locs.umiLength),
                    targetSequence = fq.bases.substring(locs.targetOffset, locs.targetOffset + locs.targetLength),
                    cut            = locs.cut
                  ))
                }
            }
        }
    }
  }

  /** Determines the number of stagger bases at the start of the read. */
  def determineStagger(bases: String, bytes: Array[Byte]): Option[Int] = {
    anchors.iterator
      .map(a => bases.substring(a.offset, a.offset + maxStagger + a.bases.length).indexOf(a.bases))
      .find(_ > -1)
      .map(_ + 1) match {
      case None =>
        None
      case Some(stagger) =>
        // Validate that all three anchors are found with <=2 mismatches in the right places
        if (anchors.forall(a => Sequences.mismatches(bytes, a.offset+stagger-1, a.bytes, 0, a.length) <= 2)) {
          Some(stagger)
        }
        else {
          None
        }
    }
  }

  /**
    * Identifies the target and UMI location in the post-primer region.
    *
    * @param fq the fastq record
    * @param bytes the sequences from the fastq record as a byte array
    * @param startOffset the offset of the first base following the JT-450 primer/probe sequence
    * @param sample the sample that the read belongs to
    */
  def locateTargetAndUmi(fq: FastqRecord, bytes: Array[Byte], startOffset: Int, sample: Sample): Option[TargetAndUmiLoc]

  /** Calculates the mean quality over a set of regions of a read. */
  def meanQuality(fq: FastqRecord, offsetsAndLength: (Int, Int)*): Int = {
    val quals = SAMUtils.fastqToPhred(fq.quals)
    var total = 0
    var count = 0

    offsetsAndLength.foreach { case (offset, len) =>
      count += len
      forloop (from=offset, until=offset+len) { i=> total += quals(i) }
    }

    require(count > 0, "Invalid to try and calculate mean quality over zero bases!")
    Math.round(total / count.toDouble).toInt
  }

  /** Prints a set of statistics about the extraction and demultiplexing to the log. */
  def printStats(logger: Logger) : Unit = {
    val total = failedToIdStagger + failedToDemultiplex + failedToExtract.total + failedQuality.total + extracted.total
    logger.info(f"Total reads             : $total%,d")
    logger.info(f"Failed to ID stagger for: $failedToIdStagger%,d (${failedToIdStagger / total.toDouble}%.3f)")
    logger.info(f"Failed to demultiplex   : $failedToDemultiplex%,d (${failedToDemultiplex / total.toDouble}%.3f)")

    val padTo = this.demultiplexer.manifest.samples.map(_.name.length).max

    this.demultiplexer.manifest.samples.foreach { sample =>
      val s = sample.name
      val padded = s.padTo(padTo, ' ')
      logger.info(f"Sample $padded: extracted ${extracted.countOf(s)}%,d; failed to extract ${failedToExtract.countOf(s)}%,d; failed quality ${failedQuality.countOf(s)}%,d")
    }
  }

  /** Generates a summary metric object containing metrics about demultiplexing and extraction that can be written to disk. */
  def summaryDemuxMetrics: DemuxMetric = {
    val total = failedToIdStagger + failedToDemultiplex + failedQuality.total + failedToExtract.total + extracted.total

    DemuxMetric(
      total_reads              = total,
      reads_extracted          = this.extracted.total,
      failed_to_id_landmarks   = this.failedToIdStagger,
      failed_to_id_sample      = this.failedToDemultiplex,
      failed_to_extract_target = this.failedToExtract.total,
      failed_quality           = this.failedQuality.total
    )
  }

  /** Generates a metric object per sample with information about the reads that matched the sample's hamming code. */
  def sampleDemuxMetrics: Seq[SampleDemuxMetric] = this.demultiplexer.manifest.samples.map { sample =>
      SampleDemuxMetric(
        sample            = sample.name,
        guide             = sample.guide,
        pam               = sample.pam,
        total_reads       = extracted.countOf(sample.name) + failedToExtract.countOf(sample.name) + failedQuality.countOf(sample.name),
        reads_extracted   = extracted.countOf(sample.name),
        failed_to_extract = failedToExtract.countOf(sample.name),
        failed_quality    = failedQuality.countOf(sample.name)
      )
    }
}

/**
  * Extracts BLT reads from fastq records for the Cas9 workflow where the read is expected
  * to look like:
  *                                              JT-450 primer/probe                                            Illumina Adapter
  * uncut: stagger[1-8]-CGATCT-rbc[6]-TACGAC-sbc[15]-TTACCGAAGATAGCAGCCTAGTGGAACC-ATCTG-target-PAM-GC-umi[12]-TGAC-AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
  * cut:   stagger[1-8]-CGATCT-rbc[6]-TACGAC-sbc[15]-TTACCGAAGATAGCAGCCTAGTGGAACC----------tgt-PAM-GC-umi[12]-TGAC-AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
  */
class Cas9ReadExtractor(override val demultiplexer: Demultiplexer,
                        override val minQuality: Int = 20,
                        fixedGuideLength: Option[Int] = None) extends ReadExtractor {
  private val PreTargetAnchor = Anchor("ATCTG", 0)
  private val PostPamSequence = "GC"
  private val PostUmiAnchor   = Anchor("TGAC", 0)

  /** The maximum total indel length that should be tolerated in the target sequence. */
  protected val maximumTargetIndel: Int = 2

  /** The maximum total indel length that should be tolerated in the target barcode/umi sequence. */
  protected val maximumTargetBarcodeIndel: Int = 1
  protected val minimumTargetBarcodeLength: Int = targetBarcodeLength - maximumTargetBarcodeIndel
  protected val maximumTargetBarcodeLength: Int = targetBarcodeLength + maximumTargetBarcodeIndel

  /** The maximum allowable length for a cut target. */
  protected val maximumCutTargetLength: Int = 8

  /**
    * Implemented to locate the target and UMI by anchoring on the fixed `ATCTG` upstream of the target,
    * the PAM+GC downstream of the target and the fixed `TGAC` after the UMI.
    *
    * @param fq the fastq record
    * @param bytes the sequences from the fastq record as a byte array
    * @param startOffset the offset of the first base following the JT-450 primer/probe sequence
    * @param sample the sample that the read belongs to
    */
  override def locateTargetAndUmi(fq: FastqRecord, bytes: Array[Byte], startOffset: Int, sample: Sample): Option[TargetAndUmiLoc] = {
    val expectedTargetLength = fixedGuideLength.getOrElse(sample.guide.length)
    val minUncutTargetLength = expectedTargetLength - maximumTargetIndel
    val maxUncutTargetLength = expectedTargetLength + maximumTargetIndel
    val pamPlus = sample.pam + PostPamSequence

    fq.bases.indexOf(pamPlus, startOffset) match {
      case -1        => None
      case pamOffset =>
        val preTargetAnchorMismatches = Sequences.mismatches(bytes, startOffset, PreTargetAnchor.bytes, 0, PreTargetAnchor.length)
        val umiOffset                 = pamOffset + pamPlus.length
        val postUmiAnchorOffset       = fq.bases.indexOf(PostUmiAnchor.bases, umiOffset + minimumTargetBarcodeLength)
        val umiLength                 = if (postUmiAnchorOffset == -1) -1 else postUmiAnchorOffset - umiOffset

        // If the target is uncut, then we expect to find the PAM after the anchor and target, but allow for 1-2 bases of wiggle
        val earliestUncutPamLoc = startOffset + PreTargetAnchor.length + minUncutTargetLength
        val targetLengthIfCut   = pamOffset - startOffset

        if (postUmiAnchorOffset == -1 || umiLength < minimumTargetBarcodeLength || umiLength > maximumTargetBarcodeLength) {
          None
        }
        else if (preTargetAnchorMismatches <= 1 && pamOffset >= earliestUncutPamLoc) {
          // Uncut target
          val padding      = this.fixedGuideLength.map(_ - sample.guide.length).getOrElse(0)
          val targetOffset = startOffset + PreTargetAnchor.length + padding
          Some(TargetAndUmiLoc(targetOffset=targetOffset, targetLength=pamOffset-targetOffset, umiOffset=umiOffset, umiLength=umiLength, cut=false))
        }
        else if (targetLengthIfCut <= maximumCutTargetLength) {
          // Cut target
          val targetOffset = startOffset
          Some(TargetAndUmiLoc(targetOffset=targetOffset, targetLength=targetLengthIfCut, umiOffset=umiOffset, umiLength=umiLength, cut=true))
        }
        else {
          // Something weird that we don't want!
          None
        }
    }
  }
}
