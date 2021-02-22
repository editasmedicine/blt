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

class DemultiplexerTest extends UnitSpec {
  /** Creates a demuxer from minimal data. */
  def buildDemuxer(mismatches:Int, distance:Int, samples: (String,String)*): Demultiplexer = {
    val ss = samples.map { case (name, barcode) => Sample(name, barcode, "AAAAAAAAAAAAAAAAAAAAAAA", Cas9, "ACGACG", cut=true) }
    new Demultiplexer(SampleManifest(ss), mismatches, distance)
  }

  "Demultiplexer.assign" should "assign obviously correct barcodes" in {
    val demuxer = buildDemuxer(0, 2, ("s1", "AAAAA"), ("s2", "CCCCC"), ("s3", "GGGGG"), ("s4", "TTTTT"))
    demuxer.assign("AAAAA".getBytes).map(_.name) shouldBe Some("s1")
    demuxer.assign("AAAAACCCCCGGGGGTTTTT".getBytes,  0).map(_.name) shouldBe Some("s1")
    demuxer.assign("AAAAACCCCCGGGGGTTTTT".getBytes,  5).map(_.name) shouldBe Some("s2")
    demuxer.assign("AAAAACCCCCGGGGGTTTTT".getBytes, 10).map(_.name) shouldBe Some("s3")
    demuxer.assign("AAAAACCCCCGGGGGTTTTT".getBytes, 15).map(_.name) shouldBe Some("s4")
  }

  it should "return None when the barcode doesn't match a valid sample" in {
    val demuxer = buildDemuxer(1, 2, ("s1", "AAAAAA"), ("s2", "CCCCCC"), ("s3", "GGGGGG"), ("s4", "TTTTTT"))
    demuxer.assign("AAAAAA".getBytes).map(_.name) shouldBe Some("s1")
    demuxer.assign("AAAAAC".getBytes).map(_.name) shouldBe Some("s1")
    demuxer.assign("AAAACC".getBytes).map(_.name) shouldBe None
    demuxer.assign("AAACCC".getBytes).map(_.name) shouldBe None
    demuxer.assign("AAGCCC".getBytes).map(_.name) shouldBe None
    demuxer.assign("AGGCCC".getBytes).map(_.name) shouldBe None
    demuxer.assign("GGGCCC".getBytes).map(_.name) shouldBe None
  }

  it should "return None when the barcode matches too closely to two samples" in {
    val demuxer = buildDemuxer(2, 2, ("s1", "ACACAC"), ("s2", "AAAAAA"), ("s3", "CCCCCC"))
    demuxer.assign("ACACAC".getBytes).map(_.name) shouldBe Some("s1")
    demuxer.assign("AAAAAA".getBytes).map(_.name) shouldBe Some("s2")
    demuxer.assign("CCCCCC".getBytes).map(_.name) shouldBe Some("s3")

    // One error to ACACAC and still 3-4 errors to other barcodes
    demuxer.assign("ACACAG".getBytes).map(_.name) shouldBe Some("s1")

    // One error to ACACAC and only two errors to AAAAAA
    demuxer.assign("ACACAA".getBytes).map(_.name) shouldBe None
  }
}
