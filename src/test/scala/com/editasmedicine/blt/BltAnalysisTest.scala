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

import com.editasmedicine.commons.testing.UnitSpec
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.alignment.Alignment
import com.fulcrumgenomics.util.{Metric => FgMetric}

class BltAnalysisTest extends UnitSpec {
  // Four samples, two pairs of cut/uncut by guide/PAM
  val s1 = Sample("g1-cut",   barcode="CACATACGCACTACG", guide="GGCCTCCCCAAAGCCTGGCCA", enzyme=Cas9, "GGGAGT", cut=true)
  val s2 = Sample("g1-naive", barcode="CCTATACCCGAATCT", guide="GGCCTCCCCAAAGCCTGGCCA", enzyme=Cas9, "GGGAGT", cut=false)
  val s3 = Sample("g2-cut",   barcode="ATTGCAAGGGCCCTT", guide="GCGGGCGGGCGGGGGAGGTG",  enzyme=Cas9, "AGGGGT", cut=true)
  val s4 = Sample("g2-naive", barcode="TCCCGTCGTCCACAA", guide="GCGGGCGGGCGGGGGAGGTG",  enzyme=Cas9, "AGGGGT", cut=false)

  val manifest = new SampleManifest(Seq(s1, s2, s3 ,s4))

  /** Method to create an analyzer with a bunch of defaults. */
  private def analyzer(out: Option[DirPath] = None, minUncut:Int=2, minIdentical:Double=0.8, useCut:Boolean=false, threads:Int=1) = {
    new BltAnalysis(
      manifest                  = manifest,
      minUncutReads             = minUncut,
      minIdenticalFraction      = minIdentical,
      useCutSamplesInValidation = useCut,
      threads                   = threads,
      output                    = out.getOrElse(makeTempDir("blt_analysis_test.", ".dir"))
    )
  }

  private val _defaultAnalyzer = analyzer()

  /** Generates a TargetInfo for use in testing. */
  private def ti(s: Sample, target: String, barcode: String, cut: Int, uncut: Int): TargetInfo = {
    val ann    = _defaultAnalyzer.annotateTarget(s.guide, target, s.enzyme)
    val cuts   = Range.inclusive(1, cut).map  (_ => BltObservation(s, barcode, Array(target.substring(0, 3)), cut=true))
    val uncuts = Range.inclusive(1, uncut).map(_ => BltObservation(s, barcode, Array(target),                 cut=false))
    TargetInfo(s.guide, s.pam, target, barcode, (cuts ++ uncuts).toArray, ann)
  }

  /** Generates a SampleTargetMetric for use in testing. */
  private def stm(s: Sample, target: String, umi: String, cut:Int, uncut: Int): SampleTargetMetric = {
    val annotation = _defaultAnalyzer.annotateTarget(s.guide, target, enzyme=Cas9)
    val alignment  = Alignment(query=s.guide, target=target, queryStart=1, targetStart=1, cigar=annotation.cigar, score=1)
    val Seq(paddedGuide, alignmentString, paddedTarget) = alignment.paddedString()


    SampleTargetMetric(
      sample                  = s.name,
      guide                   = s.guide,
      target                  = target,
      umi                     = umi,
      cigar                   = annotation.cigar.toString,
      indel_bases             = annotation.indelBases,
      mismatches              = annotation.mismatches,
      mean_mismatch_position  = annotation.meanMismatchPosition,
      mismatch_tuples         = "[]",
      genomic_location        = "",
      obs_uncut               = uncut,
      obs_cut                 = cut,
      obs_total               = cut + uncut,
      cut_rate                = cut / (cut + uncut).toDouble,
      normalized_cut_rate     = -1,
      norm_cut_rate_ci95_low  = -1,
      norm_cut_rate_ci95_high = -1,
      padded_guide            = paddedGuide,
      alignment               = alignmentString,
      padded_target           = paddedTarget
    )
  }

  //////////////////////////////////////////////////////////////////////////////
  // Basic parameter validation tests
  //////////////////////////////////////////////////////////////////////////////
  "BltAnalysis" should "throw an exception if minUncutReads < 1" in {
    val ok = analyzer(minUncut=1)
    an[Exception] shouldBe thrownBy { analyzer(minUncut = 0)}
    an[Exception] shouldBe thrownBy { analyzer(minUncut = -1)}
  }

  it should "throw an exception if minIdenticalFraction is not between 0 and 1" in {
    Seq(0.0, 0.5, 1.0).foreach { f => val ok = analyzer(minIdentical = f) }
    an[Exception] shouldBe thrownBy { analyzer(minIdentical = -0.01 )}
    an[Exception] shouldBe thrownBy { analyzer(minIdentical = 1.01 )}
  }

  it should "throw an exception if #threads is not an integer > 0" in {
    Range.inclusive(1, 16).foreach { t => val ok = analyzer(threads = t) }
    an[Exception] shouldBe thrownBy { analyzer(threads =  0)}
    an[Exception] shouldBe thrownBy { analyzer(threads = -1)}
  }


  //////////////////////////////////////////////////////////////////////////////
  // Tests for collapsing reads -> observations
  //////////////////////////////////////////////////////////////////////////////

  "BltAnalysis.buildObservations" should "collapse reads for the same sample, target barcode/umi (exact), same RBC and same stagger, cut/uncut status" in {
    val reads = Seq(
      BltRead(s1, stagger=2, randomBarcode="AAAAAA", targetBarcode="ACACACACACAC", targetSequence="ACGT", cut=false),

      BltRead(s1, stagger=2, randomBarcode="AAAAAA", targetBarcode="ACACACACACAC", targetSequence="TCTC", cut=true),
      BltRead(s1, stagger=2, randomBarcode="AAAAAA", targetBarcode="ACACACACACAC", targetSequence="TCTT", cut=true),

      BltRead(s2, stagger=2, randomBarcode="CCCCCC", targetBarcode="ACACACACACAC", targetSequence="ACGA", cut=false),
      BltRead(s2, stagger=2, randomBarcode="CCCCCC", targetBarcode="ACACACACACAC", targetSequence="ACGA", cut=false),
      BltRead(s2, stagger=2, randomBarcode="CCCCCC", targetBarcode="ACACACACACAC", targetSequence="ACGA", cut=false),

      BltRead(s3, stagger=2, randomBarcode="GGGGGG", targetBarcode="ACACACACACAC", targetSequence="AGGA", cut=false),
      BltRead(s3, stagger=2, randomBarcode="GGGGGG", targetBarcode="ACACACACACAC", targetSequence="AGGA", cut=false)
    )

    val obs = analyzer().buildObservations(reads)
    obs should have length 4
    obs.map(_.targetSequences) should contain theSameElementsAs Seq(Seq("ACGT"), Seq("TCTC", "TCTT"), Seq("ACGA", "ACGA", "ACGA"), Seq("AGGA", "AGGA"))
  }

  it should "not collapse reads if everything matches except sample" in {
    val reads = Seq(s1, s2).map(BltRead(_, stagger=2, randomBarcode="AAAAAA", targetBarcode="ACACACACACAC", targetSequence="ACGT", cut=false))
    analyzer().buildObservations(reads) should have length 2
  }

  it should "not not collapse reads if everything matches except target barcode/UMI (exact)" in {
    val reads = Seq("AAAAAAAAAAAA", "TAAAAAAAAAAA").map(umi => BltRead(s1, stagger=2, randomBarcode="AAAAAA", targetBarcode=umi, targetSequence="ACGT", cut=false))
    analyzer().buildObservations(reads) should have length 2
  }

  it should "not not collapse reads if everything matches except the random barcode (exact)" in {
    val reads = Seq("AAAAAA", "TAAAAA").map(rbc => BltRead(s1, stagger=2, randomBarcode=rbc, targetBarcode="ACACACACACAC", targetSequence="ACGT", cut=false))
    analyzer().buildObservations(reads) should have length 2
  }

  it should "not not collapse reads if everything matches except the stagger" in {
    val reads = Seq(1,2,3,4,5).map(stagger => BltRead(s1, stagger=stagger.toByte, randomBarcode="AAAAAA", targetBarcode="ACACACACACAC", targetSequence="ACGT", cut=false))
    analyzer().buildObservations(reads) should have length 5
  }

  it should "not not collapse reads if everything matches except the cut status" in {
    val reads = Seq(true, false).map(isCut => BltRead(s1, stagger=2, randomBarcode="AAAAAA", targetBarcode="ACACACACACAC", targetSequence="ACGT", cut=isCut))
    analyzer().buildObservations(reads) should have length 2
  }

  //////////////////////////////////////////////////////////////////////////////
  // Tests for collapsing reads -> observations
  //////////////////////////////////////////////////////////////////////////////

  "BltAnalysis.validateTargets" should "validate a target that has sufficient all-identical uncut observations from the same sample" in {
    val obs = Seq(
      // Two reads from a single molecule
      BltObservation(s2, targetBarcode="ACACACACACAC", targetSequences=Array("GGCCTCCCCAAAGCCTGGCCA","GGCCTCCCCAAAGCCTGGCCA"), cut=false),
      // Two reads from two different molecules
      BltObservation(s2, targetBarcode="CGCGCGCGCGCG", targetSequences=Array("GGCCTCCCCAAGCCTGGCCA"), cut=false),
      BltObservation(s2, targetBarcode="CGCGCGCGCGCG", targetSequences=Array("GGCCTCCCCAAGCCTGGCCA"), cut=false)
    )

    analyzer(minUncut=2, minIdentical=1, useCut=true).validateTargets(obs) should have size 2
  }

  it should "validate a target that has sufficient observations across samples" in {
    val obs = Seq(
      // Two reads from two different molecules across two samples, one cut sample and one uncut sample
      BltObservation(s1, targetBarcode="CGCGCGCGCGCG", targetSequences=Array("GGCCTCCCCAAAGaaTGGCCA".toUpperCase), cut=false),
      BltObservation(s2, targetBarcode="CGCGCGCGCGCG", targetSequences=Array("GGCCTCCCCAAAGaaTGGCCA".toUpperCase), cut=false)
    )

    // Should not validate if we don't allow use of cut samples (s1), but validate when using all samples
    analyzer(minUncut=2, minIdentical=1, useCut=false).validateTargets(obs) should have size 0
    analyzer(minUncut=2, minIdentical=1, useCut=true ).validateTargets(obs) should have size 1
  }

  it should "not validate a target that has sufficient observations that mismatch too much" in {
    val obs = Seq(
      // Four reads, but the last one has a mismatch as position 1
      BltObservation(s2, targetBarcode="CGCGCGCGCGCG", targetSequences=Array("GGCCTCCCCAAAGCCTGGCCA"), cut=false),
      BltObservation(s2, targetBarcode="CGCGCGCGCGCG", targetSequences=Array("GGCCTCCCCAAAGCCTGGCCA"), cut=false),
      BltObservation(s2, targetBarcode="CGCGCGCGCGCG", targetSequences=Array("GGCCTCCCCAAAGCCTGGCCA"), cut=false),
      BltObservation(s2, targetBarcode="CGCGCGCGCGCG", targetSequences=Array("AGCCTCCCCAAAGCCTGGCCA"), cut=false)
    )

    // Should not validate if we don't allow use of cut samples (s1), but validate when using all samples
    analyzer(minUncut=2, minIdentical=1,    useCut=true).validateTargets(obs) should have size 0
    analyzer(minUncut=2, minIdentical=0.8,  useCut=true).validateTargets(obs) should have size 0
    analyzer(minUncut=2, minIdentical=0.75, useCut=true).validateTargets(obs) should have size 1
  }

  // This is a bit of an odd test, because the collision is unlikely; still it's possible that two BLT
  // libraries for different guides would contain the same degenerate target, or that two BLT libraries
  // are made for the same guide but with different PAMs
  it should "not validate a target using observations across samples with different guides or PAMs" in {
    val obs = Seq(
      BltObservation(s2, targetBarcode="CGCGCGCGCGCG", targetSequences=Array("GGCCTCCCCAAAGCCTGGCCA"), cut=false),
      BltObservation(s4, targetBarcode="CGCGCGCGCGCG", targetSequences=Array("GGCCTCCCCAAAGCCTGGCCA"), cut=false)
    )

    analyzer(minUncut=2, minIdentical=1, useCut=true).validateTargets(obs) should have size 0
  }

  it should "report metrics for all observed targets, even if they did not validate" in {
    val obs = Seq(
      // Target 1: won't validate because observations are only from cut samples
      BltObservation(s1, targetBarcode="AAAAAAAAAAAA", targetSequences=Array("AAAAAAAAAAAAAAAAAAAAAA"), cut=false),
      BltObservation(s1, targetBarcode="AAAAAAAAAAAA", targetSequences=Array("AAAAAAAAAAAAAAAAAAAAAA"), cut=false),
      // Target 2: will validate because of the observations from s2
      BltObservation(s1, targetBarcode="CCCCCCCCCCCC", targetSequences=Array("CCCCCCCCCCCCCCCCCCCCCC"), cut=false),
      BltObservation(s2, targetBarcode="CCCCCCCCCCCC", targetSequences=Array("CCCCCCCCCCCCCCCCCCCCCC"), cut=false),
      BltObservation(s2, targetBarcode="CCCCCCCCCCCC", targetSequences=Array("CCCCCCCCCCCCCCCCCCCCCC"), cut=false),
      // Target 3: won't validate because the observations themselves are mostly cut
      BltObservation(s2, targetBarcode="GGGGGGGGGGGG", targetSequences=Array("GGGGGGGGGGGGGGGGGGGGGG"), cut=false),
      BltObservation(s2, targetBarcode="GGGGGGGGGGGG", targetSequences=Array("GGG"),                    cut=true),
      BltObservation(s2, targetBarcode="GGGGGGGGGGGG", targetSequences=Array("GGG"),                    cut=true),
      // Target 4: won't validate because the observations are all cut
      BltObservation(s2, targetBarcode="GAGAGAGAGAGA", targetSequences=Array("GGG"),                    cut=true),
      BltObservation(s2, targetBarcode="GAGAGAGAGAGA", targetSequences=Array("GGG"),                    cut=true),
      // Target 5: won't validate because the target sequences differ
      BltObservation(s2, targetBarcode="TTTTTTTTTTTT", targetSequences=Array("AAAAAAAAAAAACCCCCCCCCC"), cut=false),
      BltObservation(s2, targetBarcode="TTTTTTTTTTTT", targetSequences=Array("CCCCCCCCCCCCAAAAAAAAAA"), cut=false),
      // Target 6: will validate (yay!)
      BltObservation(s4, targetBarcode="CGCGCGCGCGCG", targetSequences=Array("GGCCTCACCAAAGCCTGGCCA"), cut=false),
      BltObservation(s4, targetBarcode="CGCGCGCGCGCG", targetSequences=Array("GGCCTCACCAAAGCCTGGCCA"), cut=false)
    )

    val metricsPath = makeTempFile("validation.", ".txt")
    val targets = analyzer(minUncut=2, minIdentical=0.66, useCut=false).validateTargets(obs, Some(metricsPath))

    val metrics = FgMetric.read[TargetValidationMetric](metricsPath).map(m => m.umi -> m).toMap
    metrics.size shouldBe 6

    metrics("AAAAAAAAAAAA").target                       shouldBe None // because no usable uncut reads
    metrics("AAAAAAAAAAAA").valid                        shouldBe false
    metrics("AAAAAAAAAAAA").cut_reads_in_cut_samples     shouldBe 0
    metrics("AAAAAAAAAAAA").uncut_reads_in_cut_samples   shouldBe 2
    metrics("AAAAAAAAAAAA").cut_reads_in_naive_samples   shouldBe 0
    metrics("AAAAAAAAAAAA").uncut_reads_in_naive_samples shouldBe 0
    metrics("AAAAAAAAAAAA").fraction_identical           shouldBe None // because no usable uncut reads

    metrics("CCCCCCCCCCCC").target                       shouldBe Some("CCCCCCCCCCCCCCCCCCCCCC")
    metrics("CCCCCCCCCCCC").valid                        shouldBe true
    metrics("CCCCCCCCCCCC").cut_reads_in_cut_samples     shouldBe 0
    metrics("CCCCCCCCCCCC").uncut_reads_in_cut_samples   shouldBe 1
    metrics("CCCCCCCCCCCC").cut_reads_in_naive_samples   shouldBe 0
    metrics("CCCCCCCCCCCC").uncut_reads_in_naive_samples shouldBe 2
    metrics("CCCCCCCCCCCC").fraction_identical           shouldBe Some(1.0)

    metrics("GGGGGGGGGGGG").target                       shouldBe Some("GGGGGGGGGGGGGGGGGGGGGG")
    metrics("GGGGGGGGGGGG").valid                        shouldBe false
    metrics("GGGGGGGGGGGG").cut_reads_in_cut_samples     shouldBe 0
    metrics("GGGGGGGGGGGG").uncut_reads_in_cut_samples   shouldBe 0
    metrics("GGGGGGGGGGGG").cut_reads_in_naive_samples   shouldBe 2
    metrics("GGGGGGGGGGGG").uncut_reads_in_naive_samples shouldBe 1
    metrics("GGGGGGGGGGGG").fraction_identical           shouldBe Some(1.0)

    metrics("GAGAGAGAGAGA").target                       shouldBe None
    metrics("GAGAGAGAGAGA").valid                        shouldBe false
    metrics("GAGAGAGAGAGA").cut_reads_in_cut_samples     shouldBe 0
    metrics("GAGAGAGAGAGA").uncut_reads_in_cut_samples   shouldBe 0
    metrics("GAGAGAGAGAGA").cut_reads_in_naive_samples   shouldBe 2
    metrics("GAGAGAGAGAGA").uncut_reads_in_naive_samples shouldBe 0
    metrics("GAGAGAGAGAGA").fraction_identical           shouldBe None

    metrics("TTTTTTTTTTTT").target                       shouldBe Some("AAAAAAAAAAAACCCCCCCCCC") // when 50/50 it'll just pick one
    metrics("TTTTTTTTTTTT").valid                        shouldBe false
    metrics("TTTTTTTTTTTT").cut_reads_in_cut_samples     shouldBe 0
    metrics("TTTTTTTTTTTT").uncut_reads_in_cut_samples   shouldBe 0
    metrics("TTTTTTTTTTTT").cut_reads_in_naive_samples   shouldBe 0
    metrics("TTTTTTTTTTTT").uncut_reads_in_naive_samples shouldBe 2
    metrics("TTTTTTTTTTTT").fraction_identical           shouldBe Some(0.5)

    metrics("CGCGCGCGCGCG").target                       shouldBe Some("GGCCTCACCAAAGCCTGGCCA")
    metrics("CGCGCGCGCGCG").valid                        shouldBe true
    metrics("CGCGCGCGCGCG").cut_reads_in_cut_samples     shouldBe 0
    metrics("CGCGCGCGCGCG").uncut_reads_in_cut_samples   shouldBe 0
    metrics("CGCGCGCGCGCG").cut_reads_in_naive_samples   shouldBe 0
    metrics("CGCGCGCGCGCG").uncut_reads_in_naive_samples shouldBe 2
    metrics("CGCGCGCGCGCG").fraction_identical           shouldBe Some(1.0)
  }

  //////////////////////////////////////////////////////////////////////////////
  // Tests for annotateTargets()
  //////////////////////////////////////////////////////////////////////////////

  "BltAnalysis.annotateTarget" should "correctly annotate a target that is identical to the guide" in {
    val annotation = analyzer().annotateTarget(guide="GGCCTCCCCAAAGCCTGGCCA", target="GGCCTCCCCAAAGCCTGGCCA", enzyme=Cas9)
    annotation.cigar.toString()     shouldBe "21="
    annotation.indelBases           shouldBe 0
    annotation.mismatches           shouldBe 0
    annotation.mismatchPositions    shouldBe Seq.empty
    annotation.meanMismatchPosition shouldBe None
  }

  it should "annotate a target with a few mismatches vs. the guide" in {
    val annotation = analyzer().annotateTarget(guide="GGCCTCCCCAAAGCCTGGCCA", target="GGaCTCCCCAtAGCCTGGCCg".toUpperCase, enzyme=Cas9)
    annotation.cigar.toString()     shouldBe "2=1X7=1X9=1X"
    annotation.indelBases           shouldBe 0
    annotation.mismatches           shouldBe 3
    annotation.mismatchPositions    shouldBe Seq(1, 11, 19)
    annotation.meanMismatchPosition shouldBe Some((1 + 11 + 19) / 3.0)

    annotation.copy(pamIs5PrimeOfTarget=true).mismatchPositions shouldBe Seq(3,11,21)
  }

  it should "annotate a target with adjacent mismatches vs. the guide" in {
    val annotation = analyzer().annotateTarget(guide="GGCCTCCCCAAAGCCTGGCCA", target="GGCCTCtgtAAAGCCTGGCCA".toUpperCase, enzyme=Cas9)
    annotation.cigar.toString()     shouldBe "6=3X12="
    annotation.indelBases           shouldBe 0
    annotation.mismatches           shouldBe 3
    annotation.mismatchPositions    shouldBe Seq(13,14,15)
    annotation.meanMismatchPosition shouldBe Some(14.0)

    annotation.copy(pamIs5PrimeOfTarget=true).mismatchPositions shouldBe Seq(7,8,9)
  }

  it should "annotate a target with a 1bp insertion vs. the guide" in {
    val annotation = analyzer().annotateTarget(guide="GGCCTCCCCAAAGCCTGGCCA", target="GGCCTCCCCtAAAGCCTGGCCA".toUpperCase, enzyme=Cas9)
    annotation.cigar.toString()     shouldBe "9=1D12=" // cigar is guide -> target, so shows as deletion
    annotation.indelBases           shouldBe 1
    annotation.mismatches           shouldBe 0
    annotation.mismatchPositions    shouldBe Seq.empty
    annotation.meanMismatchPosition shouldBe None
  }

  it should "annotate a target with bigger insertions vs. the guide" in {
    val annotation = analyzer().annotateTarget(guide="GGCCTCCCCAAAGCCTGGCCA", target="GGCaaCTCCCCAAAGtaCCTGGCCA".toUpperCase, enzyme=Cas9)
    annotation.cigar.toString()     shouldBe "3=2D10=2D8=" // cigar is guide -> target, so shows as deletion
    annotation.indelBases           shouldBe 4
    annotation.mismatches           shouldBe 0
    annotation.mismatchPositions    shouldBe Seq.empty
    annotation.meanMismatchPosition shouldBe None
  }

  it should "annotate a target with a deletion vs. the guide" in {
    val annotation = analyzer().annotateTarget(guide="GGCCTCCCCAAAGCCTGGCCA", target="GGCC-CCCCAAAGCCTGGCCA".replace("-", ""), enzyme=Cas9)
    annotation.cigar.toString()     shouldBe "4=1I16=" // cigar is guide -> target, so shows as insertion
    annotation.indelBases           shouldBe 1
    annotation.mismatches           shouldBe 0
    annotation.mismatchPositions    shouldBe Seq.empty
    annotation.meanMismatchPosition shouldBe None
  }

  it should "annotate a target with a deletion and compensating insertion vs. the guide" in {
    val annotation = analyzer().annotateTarget(guide="GGCCTCCCCAAAGCCTGGCCA", target="GGCC-CCCCAAAGCCTaGGCCA".replace("-", "").toUpperCase, enzyme=Cas9)
    annotation.cigar.toString()     shouldBe "4=1I11=1D5=" // cigar is guide -> target, so shows as insertion
    annotation.indelBases           shouldBe 2
    annotation.mismatches           shouldBe 0
    annotation.mismatchPositions    shouldBe Seq.empty
    annotation.meanMismatchPosition shouldBe None
  }

  it should "annotate a target with indels and mismatches vs. the guide" in {
    val annotation = analyzer().annotateTarget(guide="GGCCTCCCCAAAGCCTGGCCA", target="GGCaCTCCCCAAAGCCTGcCCA".toUpperCase, enzyme=Cas9)
    annotation.cigar.toString()     shouldBe "3=1D14=1X3="
    annotation.indelBases           shouldBe 1
    annotation.mismatches           shouldBe 1
    annotation.mismatchPositions    shouldBe Seq.empty // not given when indels present as "position" is ambiguous
    annotation.meanMismatchPosition shouldBe None
  }


  //////////////////////////////////////////////////////////////////////////////
  // Tests for sampleTargetUmiMetrics()
  //////////////////////////////////////////////////////////////////////////////

  "BltAnalysis.sampleTargetUmiMetric" should "generate metrics from simple inputs" in {
    val targets = Seq(
      ti(s1, "GGCCTCCCCAAAGCCTGGCCA".toUpperCase,  barcode="AAAAAA", cut=4, uncut=1), // target == guide
      ti(s1, "GGCCTCCgCAAAGCCTGGCCA".toUpperCase,  barcode="CCCCCC", cut=3, uncut=2), // mismatch at base 8
      ti(s1, "GGCCTCCtCCAAAGCCTGGCCA".toUpperCase, barcode="GGGGGG", cut=4, uncut=3), // 1bp insertion at base 8
      ti(s1, "GGCCTCCgCAAAGagTGGCCA".toUpperCase,  barcode="TTTTTT", cut=2, uncut=7)  // mismatches at base 8,14,15
    )

    val metrics = _defaultAnalyzer.sampleTargetUmiMetrics(s1, targets.iterator).map(m => m.umi -> m).toMap
    metrics should have size 4
    val normFactor = 4 / 5.0

    metrics("AAAAAA").cut_rate            shouldBe (4 / 5.0)
    metrics("AAAAAA").normalized_cut_rate shouldBe 1 // this entry is the only one with no edits!
    metrics("AAAAAA").mismatch_tuples     shouldBe "[]"

    metrics("CCCCCC").cut_rate            shouldBe (3 / 5.0)
    metrics("CCCCCC").normalized_cut_rate shouldBe (3 / 5.0) / normFactor
    metrics("CCCCCC").mismatch_tuples     shouldBe "[(14,C,G)]"

    metrics("GGGGGG").cut_rate            shouldBe (4 / 7.0)
    metrics("GGGGGG").normalized_cut_rate shouldBe (4 / 7.0) / normFactor
    metrics("GGGGGG").mismatch_tuples     shouldBe "[]"

    metrics("TTTTTT").cut_rate            shouldBe (2 / 9.0)
    metrics("TTTTTT").normalized_cut_rate shouldBe (2 / 9.0) / normFactor
    metrics("TTTTTT").mismatch_tuples     shouldBe "[(7,C,G),(8,C,A),(14,C,G)]"
  }

  it should "combine information across all guide==target cases for normalization" in {
    val targets = Seq(
      ti(s1, "GGCCTCCCCAAAGCCTGGCCA".toUpperCase,  barcode="AAAAAA", cut=4,  uncut=1), // target == guide
      ti(s1, "GGCCTCCCCAAAGCCTGGCCA".toUpperCase,  barcode="CCCCCC", cut=5,  uncut=1), // target == guide
      ti(s1, "GGCCTCCCCAAAGCCTGGCCA".toUpperCase,  barcode="GGGGGG", cut=12, uncut=2), // target == guide
      ti(s1, "GGCCTCCCaAAAGCCTGGCCA".toUpperCase,  barcode="TTTTTT", cut=3, uncut=2)
    )

    val metrics = _defaultAnalyzer.sampleTargetUmiMetrics(s1, targets.iterator).map(m => m.umi -> m).toMap
    metrics should have size 4
    val normFactor = (4+5+12) / (4+5+12+1+1+2).toDouble

    metrics("AAAAAA").normalized_cut_rate shouldBe ( 4 /  5.0) / normFactor
    metrics("CCCCCC").normalized_cut_rate shouldBe ( 5 /  6.0) / normFactor
    metrics("GGGGGG").normalized_cut_rate shouldBe (12 / 14.0) / normFactor
    metrics("TTTTTT").normalized_cut_rate shouldBe ( 3 /  5.0) / normFactor
  }

  it should "not combine information across samples (1/2)" in {
    val targets = Seq(
      ti(s1, "GGCCTCCCCAAAGCCTGGCCA".toUpperCase,  barcode="AAAAAA", cut=4,  uncut=1), // target == guide
      ti(s2, "GGCCTCCCCAAAGCCTGGCCA".toUpperCase,  barcode="CCCCCC", cut=5,  uncut=1), // different sample!
      ti(s1, "GGCCTCCCaAAAGCCTGGCCA".toUpperCase,  barcode="GGGGGG", cut=3, uncut=2)
    )

    val metrics = _defaultAnalyzer.sampleTargetUmiMetrics(s1, targets.iterator).map(m => m.umi -> m).toMap
    metrics should have size 2
    val normFactor = 4 / 5.0

    metrics("AAAAAA").normalized_cut_rate shouldBe (4 / 5.0) / normFactor
    metrics("GGGGGG").normalized_cut_rate shouldBe (3 / 5.0) / normFactor
  }

  it should "not combine information across samples (2/2)" in {
    val t1s1 = ti(s1, "GGCCTCCCCAAAGCCTGGCCA".toUpperCase,  barcode="AAAAAA", cut= 4,  uncut= 1) // target == guide
    val t1s2 = ti(s2, "GGCCTCCCCAAAGCCTGGCCA".toUpperCase,  barcode="AAAAAA", cut=10,  uncut=10) // target == guide
    val t2s1 = ti(s1, "GGCCTCCCaAAAGCCTGGCCA".toUpperCase,  barcode="CCCCCC", cut= 8,  uncut= 4)
    val t2s2 = ti(s2, "GGCCTCCCaAAAGCCTGGCCA".toUpperCase,  barcode="CCCCCC", cut=80,  uncut=80)

    // Targets with observations for more than one sample
    val targets = Seq(
      t1s1.copy(observations = t1s1.observations ++ t1s2.observations),
      t2s1.copy(observations = t2s1.observations ++ t2s2.observations)
    )

    val metrics = _defaultAnalyzer.sampleTargetUmiMetrics(s1, targets.iterator).map(m => m.umi -> m).toMap
    metrics should have size 2
    val normFactor = 4 / 5.0

    metrics("AAAAAA").normalized_cut_rate shouldBe (4 /  5.0) / normFactor
    metrics("CCCCCC").normalized_cut_rate shouldBe (8 / 12.0) / normFactor
  }


  //////////////////////////////////////////////////////////////////////////////
  // Tests for sampleTargetMetrics()
  //////////////////////////////////////////////////////////////////////////////

  "BltAnalysis.sampleTargetMetrics" should "roll up metrics by target sequence and re-normalize" in {
    val in = Seq(
      stm(s1, target="GGCCTCCCCAAAGCCTGGCCA",             umi="AAAAAA", cut=8, uncut=1), // target == guide
      stm(s1, target="GGCCTCCCCAAAGCCTGGCCA",             umi="CCCCCC", cut=9, uncut=1), // target == guide
      stm(s1, target="GGCCTCCCCAAAGCCTGGCCA",             umi="GGGGGG", cut=3, uncut=0), // target == guide

      stm(s1, target="GGCCTCaCCAAgaCCTGGCCA".toUpperCase, umi="ACACAC", cut=5, uncut=3),
      stm(s1, target="GGCCTCaCCAAgaCCTGGCCA".toUpperCase, umi="AGAGAG", cut=4, uncut=4),
      stm(s1, target="GGCCTCaCCAAgaCCTGGCCA".toUpperCase, umi="ATATAT", cut=2, uncut=1),

      stm(s1, target="GGCCTggggAAAGCCTGGCCA".toUpperCase, umi="CACACA", cut=2, uncut=10),
      stm(s1, target="GGCCTggggAAAGCCTGGCCA".toUpperCase, umi="CGCGCG", cut=3, uncut=4),
      stm(s1, target="GGCCTggggAAAGCCTGGCCA".toUpperCase, umi="CTCTCT", cut=1, uncut=3),

      stm(s1, target="GGCCTCCCCcccGCCTGGCCA".toUpperCase, umi="GAGAGA", cut=5, uncut=4)
    )

    val out = _defaultAnalyzer.sampleTargetMetrics(in).map(m => m.target -> m).toMap
    out should have size 4

    out("GGCCTCCCCAAAGCCTGGCCA").umi                     shouldBe "multiple"
    out("GGCCTCCCCAAAGCCTGGCCA").obs_cut                 shouldBe (8 + 9 + 3)
    out("GGCCTCCCCAAAGCCTGGCCA").obs_uncut               shouldBe (1 + 1 + 0)
    out("GGCCTCCCCAAAGCCTGGCCA").cut_rate                shouldBe (20 / 22.0)
    out("GGCCTCCCCAAAGCCTGGCCA").normalized_cut_rate     shouldBe 1.0
    out("GGCCTCCCCAAAGCCTGGCCA").norm_cut_rate_ci95_low  shouldNot be(-1.0)
    out("GGCCTCCCCAAAGCCTGGCCA").norm_cut_rate_ci95_high shouldNot be(-1.0)

    out("GGCCTCACCAAGACCTGGCCA").umi                     shouldBe "multiple"
    out("GGCCTCACCAAGACCTGGCCA").obs_cut                 shouldBe (5 + 4 + 2)
    out("GGCCTCACCAAGACCTGGCCA").obs_uncut               shouldBe (3 + 4 + 1)
    out("GGCCTCACCAAGACCTGGCCA").cut_rate                shouldBe (11 / 19.0)
    out("GGCCTCACCAAGACCTGGCCA").normalized_cut_rate     shouldBe (11/19.0) / (20/22.0)
    out("GGCCTCACCAAGACCTGGCCA").norm_cut_rate_ci95_low  shouldNot be(-1.0)
    out("GGCCTCACCAAGACCTGGCCA").norm_cut_rate_ci95_high shouldNot be(-1.0)

    out("GGCCTGGGGAAAGCCTGGCCA").umi                     shouldBe "multiple"
    out("GGCCTGGGGAAAGCCTGGCCA").obs_cut                 shouldBe ( 2 + 3 + 1)
    out("GGCCTGGGGAAAGCCTGGCCA").obs_uncut               shouldBe (10 + 4 + 3)
    out("GGCCTGGGGAAAGCCTGGCCA").cut_rate                shouldBe (6 / 23.0)
    out("GGCCTGGGGAAAGCCTGGCCA").normalized_cut_rate     shouldBe (6/23.0) / (20/22.0)
    out("GGCCTGGGGAAAGCCTGGCCA").norm_cut_rate_ci95_low  shouldNot be(-1.0)
    out("GGCCTGGGGAAAGCCTGGCCA").norm_cut_rate_ci95_high shouldNot be(-1.0)

    out("GGCCTCCCCCCCGCCTGGCCA").umi                     shouldBe "GAGAGA"
    out("GGCCTCCCCCCCGCCTGGCCA").obs_cut                 shouldBe 5
    out("GGCCTCCCCCCCGCCTGGCCA").obs_uncut               shouldBe 4
    out("GGCCTCCCCCCCGCCTGGCCA").cut_rate                shouldBe (5 / 9.0)
    out("GGCCTCCCCCCCGCCTGGCCA").normalized_cut_rate     shouldBe (5/9.0) / (20/22.0)
    out("GGCCTCCCCCCCGCCTGGCCA").norm_cut_rate_ci95_low  shouldNot be(-1.0)
    out("GGCCTCCCCCCCGCCTGGCCA").norm_cut_rate_ci95_high shouldNot be(-1.0)
  }


  //////////////////////////////////////////////////////////////////////////////
  // Tests for sampleMetrics()
  //////////////////////////////////////////////////////////////////////////////
  "BltAnalysis.sampleMetrics" should "roll up the metrics for a sample by mismatches" in {
    val in = Seq (
      stm(s1, target="GGCCTCCCCAAAGCCTGGCCA",              umi="multiple", cut=90, uncut=10), // target == guide
      stm(s1, target="GGCCaCCCCAAAGCCTGGCCA".toUpperCase,  umi="multiple", cut=80, uncut=20), // one mm
      stm(s1, target="GGCCTCCCCAAAcCCTGGCCA".toUpperCase,  umi="multiple", cut=70, uncut=30), // one mm
      stm(s1, target="GGCCTttCCAAAGCCTGGCCA".toUpperCase,  umi="multiple", cut=75, uncut=25), // two mms
      stm(s1, target="GGCCTCtCCAAAGCCaGGCCA".toUpperCase,  umi="multiple", cut=65, uncut=35), // two mms
      stm(s1, target="GGaCTCCgCAAAGgCTGGCCA".toUpperCase,  umi="multiple", cut= 6, uncut= 4), // three mms
      stm(s1, target="GGCCTtCCCAtAGCCgGGCCA".toUpperCase,  umi="multiple", cut= 5, uncut= 5), // three mms
      stm(s1, target="GGCCTgggggAAGCCTGGCCA".toUpperCase,  umi="multiple", cut=10, uncut=40), // five mms
      stm(s1, target="GGCgTCCaCAAtGCCTaaCCA".toUpperCase,  umi="multiple", cut=12, uncut=38), // five mms
      stm(s1, target="GGCgTCCCCAAGCCTGGCCA".toUpperCase, umi="multiple", cut=90, uncut=100) // one mm + one indel (should not get used)
    )

    val metrics = _defaultAnalyzer.sampleMetrics(sample=s1, targetMetrics=in).map(m => m.mismatches -> m).toMap
    metrics should have size 6 // Should have entries for zero up to max mismatches, even if there are no observations

    metrics(0).sample              shouldBe "g1-cut"
    metrics(0).mismatches          shouldBe 0
    metrics(0).targets             shouldBe 1
    metrics(0).observations        shouldBe 100
    metrics(0).cut_rate            shouldBe 0.9
    metrics(0).normalized_cut_rate shouldBe 1.0

    metrics(1).sample              shouldBe "g1-cut"
    metrics(1).mismatches          shouldBe 1
    metrics(1).targets             shouldBe 2
    metrics(1).observations        shouldBe 200
    metrics(1).cut_rate            shouldBe 0.75
    metrics(1).normalized_cut_rate shouldBe 0.75 / 0.9

    metrics(2).sample              shouldBe "g1-cut"
    metrics(2).mismatches          shouldBe 2
    metrics(2).targets             shouldBe 2
    metrics(2).observations        shouldBe 200
    metrics(2).cut_rate            shouldBe 0.7
    metrics(2).normalized_cut_rate shouldBe 0.7 / 0.9

    metrics(3).sample              shouldBe "g1-cut"
    metrics(3).mismatches          shouldBe 3
    metrics(3).targets             shouldBe 2
    metrics(3).observations        shouldBe 20
    metrics(3).cut_rate            shouldBe 0.55
    metrics(3).normalized_cut_rate shouldBe 0.55 / 0.9

    metrics(4).sample              shouldBe "g1-cut"
    metrics(4).mismatches          shouldBe 4
    metrics(4).targets             shouldBe 0
    metrics(4).observations        shouldBe 0
    metrics(4).cut_rate            shouldBe 0.0
    metrics(4).normalized_cut_rate shouldBe 0.0

    metrics(5).sample              shouldBe "g1-cut"
    metrics(5).mismatches          shouldBe 5
    metrics(5).targets             shouldBe 2
    metrics(5).observations        shouldBe 100
    metrics(5).cut_rate            shouldBe 0.22
    metrics(5).normalized_cut_rate shouldBe 0.22 / 0.9
  }

  //////////////////////////////////////////////////////////////////////////////
  // Tests for bltMetric()
  //////////////////////////////////////////////////////////////////////////////

  "BltAnalysis.bltMetric" should "calculate the GIMP score" in {
    val in = Seq(
      SampleMetric(s1.name, mismatches=0, targets=1,    observations=10,    cut_rate=0.80, normalized_cut_rate=0.80 / 0.80),
      SampleMetric(s1.name, mismatches=1, targets=10,   observations=100,   cut_rate=0.75, normalized_cut_rate=0.75 / 0.80),
      SampleMetric(s1.name, mismatches=2, targets=100,  observations=1000,  cut_rate=0.70, normalized_cut_rate=0.70 / 0.80),
      SampleMetric(s1.name, mismatches=3, targets=1000, observations=10000, cut_rate=0.55, normalized_cut_rate=0.55 / 0.80),
      SampleMetric(s1.name, mismatches=4, targets=100,  observations=100,   cut_rate=0.40, normalized_cut_rate=0.40 / 0.80),
      SampleMetric(s1.name, mismatches=5, targets=20,   observations=25,    cut_rate=0.30, normalized_cut_rate=0.30 / 0.80),
      SampleMetric(s1.name, mismatches=6, targets=10,   observations=12,    cut_rate=0.20, normalized_cut_rate=0.20 / 0.80),
      SampleMetric(s1.name, mismatches=7, targets=5,    observations=5,     cut_rate=0.25, normalized_cut_rate=0.25 / 0.80),
      SampleMetric(s1.name, mismatches=8, targets=2,    observations=2,     cut_rate=0.10, normalized_cut_rate=0.10 / 0.80)
    )

    /* Expected answers generated using the following R code:
      require(pracma)

      data = data.frame(mismatches=c(0,1,2,3,4,5,6,7,8), cut_rate=c(0.8, 0.75, 0.7, 0.55, 0.4, 0.3, 0.2, 0.25, 0.1))
      data$normalized_cut_rate = data$cut_rate / 0.8

      gimp = function(df, n) {
        ss = subset(df, df$mismatches > 0 & df$mismatches <= n)
        trapz(ss$mismatches, ss$normalized_cut_rate) / (n - 1)
      }

      for (nm in seq(3,8)) {
        print(paste("Gimp score for", nm, "is", gimp(data, nm)))
      }
    */
    _defaultAnalyzer.bltMetric(s1, in, 3).score shouldBe 0.84375 +- 0.00001
    _defaultAnalyzer.bltMetric(s1, in, 4).score shouldBe 0.76041 +- 0.00001
    _defaultAnalyzer.bltMetric(s1, in, 5).score shouldBe 0.67968 +- 0.00001
    _defaultAnalyzer.bltMetric(s1, in, 6).score shouldBe 0.60625 +- 0.00001
    _defaultAnalyzer.bltMetric(s1, in, 7).score shouldBe 0.55208 +- 0.00001
    _defaultAnalyzer.bltMetric(s1, in, 8).score shouldBe 0.50446 +- 0.00001
  }
}
