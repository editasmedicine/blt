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
import com.fulcrumgenomics.commons.util.DelimitedDataParser
import com.fulcrumgenomics.util.Io
import htsjdk.samtools.util.SequenceUtil

/**
  * Parses off-target file.  The off-target file is expected to be comma-separated without a header row, and contains
  * the following fields (laid out with spaces and headers for easier reading):
  *
  * Guide + PAM-ish               chrom  start      off-target w/PAM              strand mm off-target              off-target-loc
  * GGCCTCCCCAAAGCCTGGCCANNNNNN   chr2   72934037   GGCCTCCCCAAAGCCTGGCCAGGGAGT   -      0  GGCCTCCCCAAAGCCTGGCCA   chr2:72934037
  * GGCCTCCCCAAAGCCTGGCCANNNNNN   chr3   128255199  GGCCTCCCCAAAGagTGGCCAGAGAGG   -      2  GGCCTCCCCAAAGAGTGGCCA   chr3:128255199
  * GGCCTCCCCAAAGCCTGGCCANNNNNN   chr11  44964900   GaCCTCCCCAtAGCCTGGCCAGGGAGG   -      2  GACCTCCCCATAGCCTGGCCA   chr11:44964900

  * Actual data will look like:
  * GGCCTCCCCAAAGCCTGGCCANNNNNN,chr2,72934037,GGCCTCCCCAAAGCCTGGCCAGGGAGT,-,0,GGCCTCCCCAAAGCCTGGCCA,chr2:72934037
  * GGCCTCCCCAAAGCCTGGCCANNNNNN,chr3,128255199,GGCCTCCCCAAAGagTGGCCAGAGAGG,-,2,GGCCTCCCCAAAGAGTGGCCA,chr3:128255199
  * GGCCTCCCCAAAGCCTGGCCANNNNNN,chr11,44964900,GaCCTCCCCAtAGCCTGGCCAGGGAGG,-,2,GACCTCCCCATAGCCTGGCCA,chr11:44964900
  */
object OffTargetFile {
  private val HeaderLine = Seq("guide_with_pam", "chrom", "pos", "off_target_with_pam", "strand", "mismatches", "off_target", "loc").mkString(",")

  /**
    * Parses an off-target file, validated that the location and sequence fields look reasonable
    * and then returns them as a Map[Sequence,Location].
    */
  def parse(path: FilePath): Map[String,String] = {
    def req(test: Boolean, msg: => String) = if (!test) throw new IllegalStateException(s"Problem with off-target file $path: $msg")
    val lines  = HeaderLine +: Io.readLines(path).toSeq
    val parser = DelimitedDataParser(lines, ',')
    parser.map { row =>
      val loc = row[String]("loc")
      val seq = row[String]("off_target").toUpperCase
      req(loc.contains(":"), s"Invalid location field: $loc")
      req(seq.forall(ch => SequenceUtil.isValidBase(ch.toByte)), s"Invalid sequence: $seq")
      seq -> loc
    }.toMap
  }
}
