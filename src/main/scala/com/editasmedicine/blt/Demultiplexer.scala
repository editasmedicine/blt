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

/**
  * Class that provides the logic for determining which sample a sample barcode belongs to
  * and/or whether it can be assigned at all.
  *
  * Implements a simple algorithm whereby to assign a read to a sample the barcode must:
  * 1. Match the sample's barcode with <= `maxMismatches` mismatching bases
  * 2. Have at least `minDistance` more mismatches to all other sample barcodes, beyond
  *    the number of mismatches to the matched barcode
  *
  * @param manifest the manifest containing information about the samples being analyzed
  * @param maxMismatches the maximum number of mismatches to be tolerated in the barcode
  * @param minDistance the additional number of mismatches required to the next best barcode
  */
case class Demultiplexer(manifest: SampleManifest, maxMismatches: Int, minDistance: Int) {
  private val samples  = manifest.samples.toArray
  private val barcodes = samples.map(_.barcodeBytes)
  private val nSamples = samples.length

  /**
    * Attempts to assign the barcode to a sample.
    *
    * @param read a byte array with the barcode contained within it, and possibly other data
    * @param offset the starting (0-based) offset of the barcode bases
    */
  def assign(read: Array[Byte], offset: Int = 0): Option[Sample] = {
    val mismatches = this.barcodes.map(bc => Sequences.mismatches(read, offset, bc, 0, bc.length))
    val min = mismatches.min

    // If we have a sample with <= maxMismatches, and no other samples too close too it
    if (min <= maxMismatches && mismatches.count(_ < min + minDistance) == 1) {
      Some(samples(mismatches.indexOf(min)))
    }
    else {
      None
    }
  }
}
