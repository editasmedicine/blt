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

import java.nio.file.{Path, Paths}

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.commons.util.DelimitedDataParser
import com.fulcrumgenomics.util.Io

import scala.collection.mutable

object SampleManifest {
  // Values that are interpreted as true (case insensitive)
  val Trueish = Set("true", "yes", "t", "y")

  val ColumnNames = Seq("sample", "sample_barcode", "guide", "enzyme", "pam", "cut", "off_target_file")
  val Seq(hd_sample, hd_sample_barcode, hd_guide, hd_enzyme, hd_pam, hd_cut, hd_off_target_file) = ColumnNames

  /**
    * Parses a tab-delimited file into a SampleManifest.  The file should contain columns
    * for:
    *   sample: the name of the sample/library
    *   sample_barcode: the barcode (hamming code) for the sample
    *   guide: the original guide sequence for the sample
    *   enzyme: the enzyme in use
    *   pam: the PAM sequence for the sample
    *   cut: either true or false (or yes or no). true=sample was cut, false=Naive/uncut sample.
    *   off_target_file: [optional] the path to a file of off-target hits
    * @param path the path to the sample manifest
    */
  def apply(path: FilePath): SampleManifest = apply(Io.readLines(path).toSeq)

  /** Version of the apply method that takes in the lines as a Seq instead of reading from a file. */
  def apply(lines: Traversable[String]): SampleManifest = {
    val offTargetCache = mutable.Map[Path, Map[String,String]]()

    val parser                 = DelimitedDataParser(lines, delimiter='\t')
    val extendedAttributeNames = parser.headers.filterNot(hd => ColumnNames.contains(hd))

    val samples = parser.map { row =>
      val offTarget = row.get[String](hd_off_target_file) match {
        case None    => Map.empty[String,String]
        case Some(p) => offTargetCache.getOrElseUpdate(Paths.get(p), OffTargetFile.parse(Paths.get(p)))
      }

      val attrs : Map[String,Any] = Map(extendedAttributeNames.map(name => name -> row[String](name)):_*)

      Sample(
        name               = row[String](hd_sample),
        barcode            = row[String](hd_sample_barcode).toUpperCase,
        guide              = row[String](hd_guide).toUpperCase,
        enzyme             = Enzyme(row[String](hd_enzyme)),
        pam                = row[String](hd_pam).toUpperCase,
        cut                = Trueish.contains(row[String](hd_cut).toLowerCase),
        offTarget          = offTarget,
        extendedAttributes = attrs
      )
    }.toIndexedSeq

    // Check to make sure sample names are unique in the manifest
    val duplicateSamples = samples.groupBy(_.name).filter { case (k,v) => v.size > 1 }.keySet
    if (duplicateSamples.nonEmpty) {
      throw new IllegalStateException(s"Samples are present more than once in the manifest: ${duplicateSamples.mkString(", ")}")
    }

    // Check to make sure all barcodes are the same length!
    require(samples.map(_.barcode.length).toSet.size == 1, "Sample barcodes of different lengths in manifest.")

    SampleManifest(samples)
  }
}

/**
  * In memory representation of a sample manifest that holds the information on each
  * sample sequenced.
  */
case class SampleManifest(samples: Seq[Sample]) extends Iterable[Sample] {
  override def iterator: Iterator[Sample] = this.samples.iterator
}
