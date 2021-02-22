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

/**
  * A sample within an experiment, as described in the manifest.
  *
  * @param name the name of the sample (should be unique within a manifest)
  * @param barcode the sample barcode (hamming code) used to sequence the sample
  * @param guide the guide sequence that is being tested
  * @param enzyme the enzyme used to cut the BLT constructs
  * @param pam the PAM that was used in the BLT constructs
  * @param cut whether or not this sample was cut with the enzyme
  * @param offTarget a map of off-target hits in the genome of the format Map[sequence, chr:pos]
  * @param extendedAttributes any attributes from the sample manifest that are not defined as fields in this class
  */
case class Sample (name: String,
                   barcode: String,
                   guide: String,
                   enzyme: Enzyme,
                   pam: String,
                   cut: Boolean,
                   offTarget: Map[String,String] = Map.empty,
                   extendedAttributes: Map[String,Any] = Map.empty) {
  require(barcode.toUpperCase == barcode, "Barcode must be upper case.")
  require(guide.toUpperCase   == guide  , "Guide must be upper case.")
  require(pam.toUpperCase     == pam    , "Pam must be upper case.")

  val barcodeBytes: Array[Byte] = barcode.getBytes
  val guideBytes:   Array[Byte] = guide.getBytes
  val pamBytes:     Array[Byte] = pam.getBytes

  require(Sequences.areValidBases(barcodeBytes), s"Barcode '$barcode' for sample $name contained invalid characters.")
  require(Sequences.areValidBases(guideBytes),   s"Guide '$guide' for sample $name contained invalid characters.")
  require(Sequences.areValidBases(pamBytes),     s"Pam '$pam' for sample $name contained invalid characters.")

  /** The equivalent of notCut, for filtering when that's the desired property. */
  def naive: Boolean = !cut

  override def toString: String = name
}
