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

class OffTargetFileTest extends UnitSpec {
  "OffTargetFile.parse" should "parse a valid file" in {
    val lines =
    """
      |GGCCTCCCCAAAGCCTGGCCANNNNNN,chr2,72934037,GGCCTCCCCAAAGCCTGGCCAGGGAGT,-,0,GGCCTCCCCAAAGCCTGGCCA,chr2:72934037
      |GGCCTCCCCAAAGCCTGGCCANNNNNN,chr3,128255199,GGCCTCCCCAAAGagTGGCCAGAGAGG,-,2,GGCCTCCCCAAAGAGTGGCCA,chr3:128255199
      |GGCCTCCCCAAAGCCTGGCCANNNNNN,chr11,44964900,GaCCTCCCCAtAGCCTGGCCAGGGAGG,-,2,GACCTCCCCATAGCCTGGCCA,chr11:44964900
      |GGCCTCCCCAAAGCCTGGCCANNNNNN,chr7,41409260,GGCCTCCtCAgAGCCTtGCCAGAGGAG,+,3,GGCCTCCTCAGAGCCTTGCCA,chr7:41409260
      |GGCCTCCCCAAAGCCTGGCCANNNNNN,chr14,24103667,attCTCCCCAAAGCCTGGCCACAGGAA,-,3,ATTCTCCCCAAAGCCTGGCCA,chr14:24103667
    """.stripMargin.trim.lines.toSeq

    val in = makeTempFile("offtarget.", ".csv")
    Io.writeLines(in, lines)
    val map = OffTargetFile.parse(in)
    map.size shouldBe 5
    map("GGCCTCCCCAAAGCCTGGCCA") shouldBe "chr2:72934037"
    map("GGCCTCCCCAAAGAGTGGCCA") shouldBe "chr3:128255199"
    map("GACCTCCCCATAGCCTGGCCA") shouldBe "chr11:44964900"
    map("GGCCTCCTCAGAGCCTTGCCA") shouldBe "chr7:41409260"
    map("ATTCTCCCCAAAGCCTGGCCA") shouldBe "chr14:24103667"
  }

  it should "parse an empty file" in {
    val in = makeTempFile("offtarget.", ".csv")
    Io.writeLines(in, Seq())
    val map = OffTargetFile.parse(in)
    map.size shouldBe 0
  }

  it should "fail if the seq column has things that don't look like sequences because column ordering is wrong" in {
    val lines =
    """
      |GGCCTCCCCAAAGCCTGGCCANNNNNN,chr2,72934037,GGCCTCCCCAAAGCCTGGCCAGGGAGT,GGCCTCCCCAAAGCCTGGCCA,-,0,chr2:72934037
      |GGCCTCCCCAAAGCCTGGCCANNNNNN,chr3,128255199,GGCCTCCCCAAAGagTGGCCAGAGAGG,GGCCTCCCCAAAGAGTGGCCA,-,2,chr3:128255199
    """.stripMargin.trim.lines.toSeq

    val in = makeTempFile("offtarget.", ".csv")
    Io.writeLines(in, lines)
    an[Exception] shouldBe thrownBy { OffTargetFile.parse(in) }
  }

  it should "fail if the location column has something that doesn't look like a location in it" in {
    val lines =
    """
      |GGCCTCCCCAAAGCCTGGCCANNNNNN,chr2,72934037,GGCCTCCCCAAAGCCTGGCCAGGGAGT,-,0,GGCCTCCCCAAAGCCTGGCCA,nope
      |GGCCTCCCCAAAGCCTGGCCANNNNNN,chr3,128255199,GGCCTCCCCAAAGagTGGCCAGAGAGG,-,2,GGCCTCCCCAAAGAGTGGCCA,nope
    """.stripMargin.trim.lines.toSeq

    val in = makeTempFile("offtarget.", ".csv")
    Io.writeLines(in, lines)
    an[Exception] shouldBe thrownBy { OffTargetFile.parse(in) }
  }
}
