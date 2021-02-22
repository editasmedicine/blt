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
import com.fulcrumgenomics.fastq.FastqRecord
import htsjdk.samtools.SAMUtils

class ReadExtractorTest extends UnitSpec {
  val s1 = Sample("s1", barcode="CACATACGCACTACG", guide="GGCCTCCCCAAAGCCTGGCCA", enzyme=Cas9, "GGGAGT", cut=true)
  val s2 = Sample("s2", barcode="CCTATACCCGAATCT", guide="GGCCTCCCCAAAGCCTGGCCA", enzyme=Cas9, "GGGAGT", cut=true)
  val s3 = Sample("s3", barcode="TATACAATTCGCAGC", guide="CCCAGTGTCCCCCTTCCCTAT", enzyme=Cas9, "GGGAAT", cut=true)
  val s4 = Sample("s4", barcode="CCGGAGTAGGTCCTC", guide="CCCAGTGTCCCCCTTCCCTAT", enzyme=Cas9, "GGGAAT", cut=true)
  val s5 = Sample("s5", barcode="ATTGCAAGGGCCCTT", guide="GCGGGCGGGCGGGGGAGGTG",  enzyme=Cas9, "AGGGGT", cut=true)
  val s6 = Sample("s6", barcode="TCCCGTCGTCCACAA", guide="GCGGGCGGGCGGGGGAGGTG",  enzyme=Cas9, "AGGGGT", cut=true)

  // A demultiplexer that has 6 samples with three different targets/PAMs and 6 of the real sample barcodes in use
  val demuxer = Demultiplexer(SampleManifest(Seq(s1, s2, s3, s4, s5, s6)), maxMismatches = 2, minDistance   = 2)

  /** Creates a fastq record with the given bases and a quality string of all `qual`. */
  private def fq(bases: String, qual: Int = 30) = {
    // Remove any layout junk from the sequence
    val bs = bases.toUpperCase.filter(b => Sequences.isValidBase(b.toByte, allowAmbiguityCodes=true))
    FastqRecord(name="q", bases=bs, quals=SAMUtils.phredToFastq(qual).toString * bs.length, comment=None, readNumber=None)
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // The following is reproduced from the Cas9ReadExtractor doc to aid in understanding the test cases
  //
  // The structure of the reads for Cas9 BLT is expected to be as follows (where upper case sequences is FIXED):
  //   stagger[1-8]-CGATCT-rbc[6]-TACGAC-sbc[15]-TTACCGAAGATAGCAGCCTAGTGGAACC-ATCTG-target-PAM-GC-umi[12]-TGAC-AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  "Cas9ReadExtractor" should "extract a read that has no errors in it" in {
    val extractor = new Cas9ReadExtractor(demuxer, minQuality=2)
    //                              stagger   a1     rbc     a2     s1 barcode                a3                a4         target            pam   gc     umi       a6    p7 adapter
    val r1 = extractor.extract(fq("       G-cgatct-AAAAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))
    val r2 = extractor.extract(fq("      GC-cgatct-ACAAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))
    val r3 = extractor.extract(fq("     GCG-cgatct-AACAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))
    val r4 = extractor.extract(fq("    GCGG-cgatct-AAAACA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))
    val r5 = extractor.extract(fq("   GCGGA-cgatct-AAAAAC-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))
    val r6 = extractor.extract(fq("  GCGGAA-cgatct-ATAAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))
    val r7 = extractor.extract(fq(" GCGGAAG-cgatct-AATAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))
    val r8 = extractor.extract(fq("GCGGAAGC-cgatct-AAATAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))

    r1 shouldBe Some(BltRead(s1, stagger=1, randomBarcode="AAAAAA", targetBarcode="AAAAAAAAAAAA", targetSequence="GGCCTCCCCAAAGCCTGGCCA", cut=false))
    r2 shouldBe Some(BltRead(s1, stagger=2, randomBarcode="ACAAAA", targetBarcode="AAAAAAAAAAAA", targetSequence="GGCCTCCCCAAAGCCTGGCCA", cut=false))
    r3 shouldBe Some(BltRead(s1, stagger=3, randomBarcode="AACAAA", targetBarcode="AAAAAAAAAAAA", targetSequence="GGCCTCCCCAAAGCCTGGCCA", cut=false))
    r4 shouldBe Some(BltRead(s1, stagger=4, randomBarcode="AAAACA", targetBarcode="AAAAAAAAAAAA", targetSequence="GGCCTCCCCAAAGCCTGGCCA", cut=false))
    r5 shouldBe Some(BltRead(s1, stagger=5, randomBarcode="AAAAAC", targetBarcode="AAAAAAAAAAAA", targetSequence="GGCCTCCCCAAAGCCTGGCCA", cut=false))
    r6 shouldBe Some(BltRead(s1, stagger=6, randomBarcode="ATAAAA", targetBarcode="AAAAAAAAAAAA", targetSequence="GGCCTCCCCAAAGCCTGGCCA", cut=false))
    r7 shouldBe Some(BltRead(s1, stagger=7, randomBarcode="AATAAA", targetBarcode="AAAAAAAAAAAA", targetSequence="GGCCTCCCCAAAGCCTGGCCA", cut=false))
    r8 shouldBe Some(BltRead(s1, stagger=8, randomBarcode="AAATAA", targetBarcode="AAAAAAAAAAAA", targetSequence="GGCCTCCCCAAAGCCTGGCCA", cut=false))
  }

  it should "extract reads with cut and uncut targets" in {
    val extractor = new Cas9ReadExtractor(demuxer, minQuality=2)
    val r1 = extractor.extract(fq("     G-cgatct-AAAAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))
    val r2 = extractor.extract(fq("GCGGAA-cgatct-AAAAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc--------------------------CTG-GGGAGT-GC-TTTTTTTTTTTT-tgac-agatcggaagagcac"))

    r1 shouldBe Some(BltRead(s1, stagger=1, randomBarcode="AAAAAA", targetBarcode="AAAAAAAAAAAA", targetSequence="GGCCTCCCCAAAGCCTGGCCA", cut=false))
    r2 shouldBe Some(BltRead(s1, stagger=6, randomBarcode="AAAAAA", targetBarcode="TTTTTTTTTTTT", targetSequence="CTG", cut=true))
  }

  it should "extract reads with errors in anchors 1-3, so long as one a perfect match" in {
    val extractor = new Cas9ReadExtractor(demuxer, minQuality=2)
    val r1 = extractor.extract(fq("G-cgnnct-AAAAAA-tannac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))
    val r2 = extractor.extract(fq("G-cgnnct-AAAAAA-tacgac-CACATACGCACTACG-nnaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))
    val r3 = extractor.extract(fq("G-cgatct-AAAAAA-tannac-CACATACGCACTACG-nnaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))

    r1 shouldBe Some(BltRead(s1, stagger=1, randomBarcode="AAAAAA", targetBarcode="AAAAAAAAAAAA", targetSequence="GGCCTCCCCAAAGCCTGGCCA", cut=false))
    r2 shouldBe Some(BltRead(s1, stagger=1, randomBarcode="AAAAAA", targetBarcode="AAAAAAAAAAAA", targetSequence="GGCCTCCCCAAAGCCTGGCCA", cut=false))
    r3 shouldBe Some(BltRead(s1, stagger=1, randomBarcode="AAAAAA", targetBarcode="AAAAAAAAAAAA", targetSequence="GGCCTCCCCAAAGCCTGGCCA", cut=false))
  }

  it should "fail to extract when all three anchors have errors" in {
    val extractor = new Cas9ReadExtractor(demuxer, minQuality=2)
    val r1 = extractor.extract(fq("G-cgnnct-AAAAAA-tannac-CACATACGCACTACG-nnaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))
    r1 shouldBe None
  }

  it should "fail to extract reads when there are indels in anchors 1-3" in {
    val extractor = new Cas9ReadExtractor(demuxer, minQuality=2)
    val r1 = extractor.extract(fq("G-cgatctt-AAAAAA-tac-gac-CACATACGCACTACG-ttaccgaa-gatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))
    val r2 = extractor.extract(fq("G-cgatct--AAAAAA-taccgac-CACATACGCACTACG-ttaccgaa-gatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))
    val r3 = extractor.extract(fq("G-cgatct--AAAAAA-tac-gac-CACATACGCACTACG-ttaccgaaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))

    r1 shouldBe None
    r2 shouldBe None
    r3 shouldBe None
  }

  it should "fail to extract reads when there are indels in the random barcode or the sample barcode" in {
    val extractor = new Cas9ReadExtractor(demuxer, minQuality=2)
    // First read has a 7-base RBC instead of 6, second read has a 14-base sample barcode instead of 15
    val r1 = extractor.extract(fq("G-cgatct-AAAAAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))
    val r2 = extractor.extract(fq("G-cgatct-AAAAAA--tacgac-CACA-ACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))

    r1 shouldBe None
    r2 shouldBe None
  }

  it should "fail to extract when the bases are low quality" in {
    val extractor = new Cas9ReadExtractor(demuxer, minQuality=30)
    val r1 = extractor.extract(fq("G-cgatct-AAAAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac", qual=40))
    val r2 = extractor.extract(fq("G-cgatct-AAAAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac", qual=20))

    r1.isDefined shouldBe true
    r2 shouldBe None
  }

  it should "allow a 1bp insertion or deletion in the umi/target barcode but no more" in {
    val extractor = new Cas9ReadExtractor(demuxer, minQuality=2)
    val std  = extractor.extract(fq("GCG-cgatct-AACAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))
    val ins1 = extractor.extract(fq("GCG-cgatct-AACAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAAA-tgac-agatcggaagagcac"))
    val del1 = extractor.extract(fq("GCG-cgatct-AACAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAA-tgac-agatcggaagagcac"))
    val ins2 = extractor.extract(fq("GCG-cgatct-AACAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAAAA-tgac-agatcggaagagcac"))
    val del2 = extractor.extract(fq("GCG-cgatct-AACAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAA-tgac-agatcggaagagcac"))

    std  shouldBe Some(BltRead(s1, stagger=3, randomBarcode="AACAAA", targetBarcode="AAAAAAAAAAAA",  targetSequence="GGCCTCCCCAAAGCCTGGCCA", cut=false))
    ins1 shouldBe Some(BltRead(s1, stagger=3, randomBarcode="AACAAA", targetBarcode="AAAAAAAAAAAAA", targetSequence="GGCCTCCCCAAAGCCTGGCCA", cut=false))
    del1 shouldBe Some(BltRead(s1, stagger=3, randomBarcode="AACAAA", targetBarcode="AAAAAAAAAAA",   targetSequence="GGCCTCCCCAAAGCCTGGCCA", cut=false))
    ins2 shouldBe None
    del2 shouldBe None
  }

  it should "fail to extract when the exact PAM+GC sequence cannot be found" in {
    val extractor = new Cas9ReadExtractor(demuxer, minQuality=2)
    val std  = extractor.extract(fq("GCG-cgatct-AACAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))
    val err1 = extractor.extract(fq("GCG-cgatct-AACAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGgGg-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))
    val err2 = extractor.extract(fq("GCG-cgatct-AACAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-Gg-AAAAAAAAAAAA-tgac-agatcggaagagcac"))
    val ins1 = extractor.extract(fq("GCG-cgatct-AACAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GaGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))
    val del1 = extractor.extract(fq("GCG-cgatct-AACAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))

    std.isDefined shouldBe true
    err1 shouldBe None
    err2 shouldBe None
    ins1 shouldBe None
    del1 shouldBe None
  }

  it should "fail to extract when the terminal TGAC anchor cannot be found" in {
    val extractor = new Cas9ReadExtractor(demuxer, minQuality=2)
    val std  = extractor.extract(fq("GCG-cgatct-AACAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"))
    val err1 = extractor.extract(fq("GCG-cgatct-AACAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-ttac-agatcggaagagcac"))
    val err2 = extractor.extract(fq("GCG-cgatct-AACAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgat-agatcggaagagcac"))
    val ins1 = extractor.extract(fq("GCG-cgatct-AACAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tggac-agatcggaagagcac"))
    val del1 = extractor.extract(fq("GCG-cgatct-AACAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-gac-agatcggaagagcac"))

    std.isDefined shouldBe true
    err1 shouldBe None
    err2 shouldBe None
    ins1 shouldBe None
    del1 shouldBe None
  }

  it should "extract reads appropriately when the fixedGuideLength option is specified" in {
    val fixed    = new Cas9ReadExtractor(demuxer, minQuality=2, fixedGuideLength=Some(23))
    val variable = new Cas9ReadExtractor(demuxer, minQuality=2, fixedGuideLength=None)

    //            stagger   a1     rbc     a2     s6 barcode                a3                a4  pad      target          pam   gc     umi       a6    p7 adapter
    val uncut = fq("GCGGA-cgatct-AAAAAC-tacgac-TCCCGTCGTCCACAA-ttaccgaagatagcagcctagtggaacc-atctg-nnnGCGGGCGGGCGGGGGAGGTG-AGGGGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac")
    val cut   = fq("GCGGA-cgatct-AAAAAC-tacgac-TCCCGTCGTCCACAA-ttaccgaagatagcagcctagtggaacc---------------------------TAC-AGGGGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac")

    // This is good!
    fixed.extract(uncut) shouldBe Some(BltRead(sample=s6, stagger=5, randomBarcode="AAAAAC", targetBarcode="AAAAAAAAAAAA", targetSequence="GCGGGCGGGCGGGGGAGGTG", cut=false))
    fixed.extract(cut)   shouldBe Some(BltRead(sample=s6, stagger=5, randomBarcode="AAAAAC", targetBarcode="AAAAAAAAAAAA", targetSequence="TAC", cut=true))

    // This is what happens when you specified no fixed guide length, but you did pad
    variable.extract(uncut) shouldBe Some(BltRead(sample=s6, stagger=5, randomBarcode="AAAAAC", targetBarcode="AAAAAAAAAAAA", targetSequence="NNNGCGGGCGGGCGGGGGAGGTG", cut=false))
    variable.extract(cut)   shouldBe Some(BltRead(sample=s6, stagger=5, randomBarcode="AAAAAC", targetBarcode="AAAAAAAAAAAA", targetSequence="TAC", cut=true))
  }

  it should "generate metrics during extraction" in {
    val reads = Seq(
      fq(/* s1 ok     */ " GC-cgatct-ACAAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"),
      fq(/* s1 ok     */ "GCG-cgatct-AACAAA-tacgac-CACATACGCACTACG-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"),
      fq(/* s2 ok     */ "GCG-cgatct-AACAAA-tacgac-CCTATACCCGAATCT-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"),
      fq(/* s2 bad tgt*/ "GCG-cgatct-AACAAA-tacgac-CCTATACCCGAATCT-ttaccgaagatagcagcctagtggaacc-NNNNN-GGCCTCCCCAAAGCCTGGCCA-NNNNN-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"),
      fq(/* junk      */ "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"),
      fq(/* no sample */ "GCG-cgatct-AACAAA-tacgac-NNNNNNNNNNNNNNN-ttaccgaagatagcagcctagtggaacc-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac"),
      fq(/* bad stagr */ "GCG-NNNNNN-NNNNNN-NNNNNN-CACATACGCACTACG-NNNNNNNNNNNNNNNNNNNNNNNNNNNN-atctg-GGCCTCCCCAAAGCCTGGCCA-GGGAGT-GC-AAAAAAAAAAAA-tgac-agatcggaagagcac")
    )

    val extractor = new Cas9ReadExtractor(demuxer, minQuality=2, fixedGuideLength=None)
    reads.foreach(extractor.extract)

    val summaryMetrics = extractor.summaryDemuxMetrics
    summaryMetrics.total_reads              shouldBe 7
    summaryMetrics.reads_extracted          shouldBe 3
    summaryMetrics.failed_to_id_landmarks   shouldBe 2
    summaryMetrics.failed_to_id_sample      shouldBe 1
    summaryMetrics.failed_to_extract_target shouldBe 1
    summaryMetrics.failed_quality           shouldBe 0

    val sampleMetrics  = extractor.sampleDemuxMetrics.map(m => m.sample -> m).toMap
    sampleMetrics("s1").total_reads     shouldBe 2
    sampleMetrics("s1").reads_extracted shouldBe 2
    sampleMetrics("s2").total_reads     shouldBe 2
    sampleMetrics("s2").reads_extracted shouldBe 1
    Seq("s3", "s4", "s5", "s6").foreach { s =>
      sampleMetrics(s).total_reads shouldBe     0
      sampleMetrics(s).reads_extracted shouldBe 0
    }
  }
}
