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

class SequencesTest extends UnitSpec {
  "Sequences.mismatches" should "return 0 on two empty strings" in {
    Sequences.mismatches("".getBytes, "".getBytes) shouldBe 0
  }

  it should "throw an exception on two strings of different lengths" in {
    an[Exception] shouldBe thrownBy { Sequences.mismatches("AC".getBytes, "ACT".getBytes) }
    an[Exception] shouldBe thrownBy { Sequences.mismatches("AC".getBytes, "ACT".getBytes) }
  }

  it should "throw an exception if sub-ranges are invalid" in {
    an[Exception] shouldBe thrownBy { Sequences.mismatches("ACTG".getBytes,   2, "ACTGTG".getBytes, 2, 4) }
    an[Exception] shouldBe thrownBy { Sequences.mismatches("ACTGTG".getBytes, 2, "ACTG".getBytes, 2, 4) }
  }

  it should "correctly calculate mismatches" in {
    Sequences.mismatches("ACGT".getBytes, "ACGT".getBytes) shouldBe 0
    Sequences.mismatches("ACGT".getBytes, "ACCT".getBytes) shouldBe 1
    Sequences.mismatches("ACGT".getBytes, "ACCC".getBytes) shouldBe 2
    Sequences.mismatches("ACGT".getBytes, "TGCA".getBytes) shouldBe 4
    Sequences.mismatches("ACGTACGT".getBytes, 4, "AGCTAACT".getBytes, 4, 4) shouldBe 2
  }

  it should "cap the number of mismatches returned at the maximum" in {
    val (s1, s2) = ("AAAAAAAAAA".getBytes, "TTTTTTTTTT".getBytes)
    Sequences.mismatches(s1, s2) shouldBe 10
    Sequences.mismatches(s1, 0, s2, 0, s1.length) shouldBe 10
    Sequences.mismatches(s1, 0, s2, 0, s1.length, max=4) shouldBe 4
  }

  "Sequences.isValidBase" should "return true for valid bases" in {
    "ACGT".getBytes.foreach(Sequences.isValidBase(_) shouldBe true)
    "ACGTNYR".getBytes.foreach(Sequences.isValidBase(_, allowAmbiguityCodes=true) shouldBe true)
  }

  it should "return false for invalid bases" in {
    " .-+*".getBytes.foreach(Sequences.isValidBase(_) shouldBe false)
    " .-+*".getBytes.foreach(Sequences.isValidBase(_, allowAmbiguityCodes=true) shouldBe false)
  }

  "Sequences.areValidBases" should "return true on an empty sequence" in {
    Sequences.areValidBases("".getBytes) shouldBe true
  }

  it should "return true on a string of entirely valid bases" in {
    Sequences.areValidBases("TGATCGTAGCTGCGCGTCATGA".getBytes) shouldBe true
  }

  it should "return false if there are any invalid bases" in {
    Sequences.areValidBases("TGATCGTAGCTGCGCGT.ATGA".getBytes) shouldBe false
  }

  it should "return true on IUPAC codes if they are allowed and false otherwise" in {
    val fourBase = "ACGT".getBytes
    val allBases = "ACGTRYSWKMBDHVN".getBytes
    Sequences.areValidBases(fourBase, allowAmbiguityCodes=false) shouldBe true
    Sequences.areValidBases(fourBase, allowAmbiguityCodes=true)  shouldBe true
    Sequences.areValidBases(allBases, allowAmbiguityCodes=false) shouldBe false
    Sequences.areValidBases(fourBase, allowAmbiguityCodes=true)  shouldBe true
  }
}
