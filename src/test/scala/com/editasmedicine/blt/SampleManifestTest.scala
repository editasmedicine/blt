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
import com.fulcrumgenomics.util.Io

class SampleManifestTest extends UnitSpec {
  /** A valid manifest with two samples. */
  private val GoodManifestText =
    """sample	sample_barcode	guide	enzyme	pam	off_target_file	cut	ext1	ext2
      |s1	CACATACGCACTACG	GGCCTCCCCAAAGCCTGGCCA	Cas9	GGGAGT		yes	hello	world
      |s2	ACACACACACACACA	GGCCTCCCCAAAGCCTGGCCA	Cas9	GGGAGT		true	some	thing
    """.stripMargin.trim

  "SampleManifest" should "parse a valid sample manifest" in {
    val path = makeTempFile("manifest.", ".txt")
    Io.writeLines(path, GoodManifestText.lines.toSeq)

    val manifest = SampleManifest(path)
    manifest.size shouldBe 2

    val Seq(s1, s2) = manifest.samples
    s1.name    shouldBe "s1"
    s1.barcode shouldBe "CACATACGCACTACG"
    s1.guide   shouldBe "GGCCTCCCCAAAGCCTGGCCA"
    s1.enzyme  shouldBe Cas9
    s1.pam     shouldBe "GGGAGT"
    s1.cut     shouldBe true
    s1.extendedAttributes.size shouldBe 2
    s1.extendedAttributes("ext1") shouldBe "hello"
    s1.extendedAttributes("ext2") shouldBe "world"

    s2.name    shouldBe "s2"
    s2.barcode shouldBe "ACACACACACACACA"
    s2.guide   shouldBe "GGCCTCCCCAAAGCCTGGCCA"
    s2.enzyme  shouldBe Cas9
    s2.pam     shouldBe "GGGAGT"
    s2.cut     shouldBe true
    s2.extendedAttributes.size shouldBe 2
    s2.extendedAttributes("ext1") shouldBe "some"
    s2.extendedAttributes("ext2") shouldBe "thing"
  }

  it should "fail if one or more columns are missing (or misspelled)" in {
    val lines = GoodManifestText.replace("guide", "target").lines.toSeq
    an[Exception] shouldBe thrownBy { SampleManifest(lines) }
  }

  it should "fail if there are duplicate sample names in the manifest" in {
    val lines = GoodManifestText.replace("s2", "s1").lines.toSeq
    an[Exception] shouldBe thrownBy { SampleManifest(lines) }
  }

  it should "fail if the sample barcodes are different lengths" in {
    val lines = GoodManifestText.replace("ACACACACACACACA", "ACACACACAC").lines.toSeq
    an[Exception] shouldBe thrownBy { SampleManifest(lines) }
  }

  it should "upper-case sequences as they are parsed" in {
    val lines = GoodManifestText.toLowerCase.lines.toSeq
    SampleManifest(lines).foreach { s =>
      s.barcode shouldBe s.barcode.toUpperCase
      s.guide   shouldBe s.guide.toUpperCase
      s.pam     shouldBe s.pam.toUpperCase
    }
  }

  it should "interpret true or yes or t or y case insensitive as true, all other values as false" in {
    SampleManifest(GoodManifestText.replace("true", "TRUE").replace("yes", "YES").lines.toSeq).foreach { s => s.cut shouldBe true }
    SampleManifest(GoodManifestText.replace("true", "TrUe").replace("yes", "Yes").lines.toSeq).foreach { s => s.cut shouldBe true }
    SampleManifest(GoodManifestText.replace("true", "t").replace("yes", "Y").lines.toSeq).foreach { s => s.cut shouldBe true }
    SampleManifest(GoodManifestText.replace("true", "false").replace("yes", "no").lines.toSeq).foreach { s => s.cut shouldBe false }
    SampleManifest(GoodManifestText.replace("true", "yeppers").replace("yes", "true-ish").lines.toSeq).foreach { s => s.cut shouldBe false }
  }

  it should "fail if there are invalid bases in any of the sequence fields" in {
    an[Exception] shouldBe thrownBy { SampleManifest(GoodManifestText.replace("ACACACACACACACA", "A.AC-CACACACACA").lines.toSeq) }
    an[Exception] shouldBe thrownBy { SampleManifest(GoodManifestText.replace("GGCCTCCCCAAAGCCTGGCCA", "'GGCCTCCCCAAAGCCTGGCCA'").lines.toSeq) }
    an[Exception] shouldBe thrownBy { SampleManifest(GoodManifestText.replace("GGGAGT", "NotPam").lines.toSeq) }
  }

  it should "parse a manifest that has an off-target file in it" in {
    val offTargetLines =
      """
        |GGCCTCCCCAAAGCCTGGCCANNNNNN,chr2,72934037,GGCCTCCCCAAAGCCTGGCCAGGGAGT,-,0,GGCCTCCCCAAAGCCTGGCCA,chr2:72934037
        |GGCCTCCCCAAAGCCTGGCCANNNNNN,chr3,128255199,GGCCTCCCCAAAGagTGGCCAGAGAGG,-,2,GGCCTCCCCAAAGAGTGGCCA,chr3:128255199
      """.stripMargin.trim.lines.toSeq

    val offTargetPath = makeTempFile("off_target.", ".csv")
    Io.writeLines(offTargetPath, offTargetLines)

    val manifestLines =
   s"""
      |sample	sample_barcode	guide	enzyme	pam	off_target_file	cut
      |s1	CACATACGCACTACG	GGCCTCCCCAAAGCCTGGCCA	Cas9	GGGAGT		yes
      |s2	ACACACACACACACA	GGCCTCCCCAAAGCCTGGCCA	Cas9	GGGAGT	$offTargetPath	true
    """.stripMargin.trim.lines.toSeq

    val manifest = SampleManifest(manifestLines)
    manifest.size shouldBe 2
    manifest.samples(0).offTarget shouldBe Map.empty

    val offTarget = manifest.samples(1).offTarget
    offTarget.size shouldBe 2
    offTarget("GGCCTCCCCAAAGCCTGGCCA") shouldBe "chr2:72934037"
    offTarget("GGCCTCCCCAAAGAGTGGCCA") shouldBe "chr3:128255199"
  }

  "Sample" should "throw exceptions if lower case sequence get set on it" in {
    val valid = Sample(name="hello", barcode="ACGT", guide="CCCC", enzyme=Cas9, pam="GGGG", cut=true)
    an[Exception] shouldBe thrownBy { valid.copy(barcode="acgt") }
    an[Exception] shouldBe thrownBy { valid.copy(guide  ="cccc") }
    an[Exception] shouldBe thrownBy { valid.copy(pam    ="gggg") }
  }
}
