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

import com.fulcrumgenomics.alignment.Cigar
import htsjdk.samtools.CigarOperator

import scala.collection.mutable.ListBuffer
import com.fulcrumgenomics.FgBioDef._


object BltRead {
  /**
    * Provides an ordering that will sort together reads believed to be from the same target sequence.
    * Reads are assumed to be derived from the same target sequence if they have the same target
    * barcode, the same original guide sequence (before degeneracy) and the same expected PAM.
    */
  object TargetOrdering extends Ordering[BltRead] {
    override def compare(x: BltRead, y: BltRead): Int = {
      var retval = x.targetBarcode.compareTo(y.targetBarcode)
      if (retval == 0) retval = x.sample.guide.compareTo(y.sample.guide)
      if (retval == 0) retval = x.sample.pam.compareTo(y.sample.pam)

      retval
    }
  }

  /**
    * Provides an ordering that will sort together reads believed to be duplicates of one another.
    * Reads are assumed to be duplicates if they are equal according to the [[TargetOrdering]] and
    * they have the same stagger, random barcode, are for the same sample, and are either both cut
    * or both not cut.
    */
  object DuplicateOrdering extends Ordering[BltRead] {
    override def compare(x: BltRead, y: BltRead): Int = {
      var retval = TargetOrdering.compare(x, y)
      if (retval == 0) retval = x.stagger - y.stagger
      if (retval == 0) retval = x.randomBarcode.compareTo(y.randomBarcode)
      if (retval == 0) retval = x.sample.name.compareTo(y.sample.name)
      if (retval == 0) retval = x.cut.compareTo(y.cut)

      retval
    }
  }
}


/**
  * The important parts of a BLT read extracted from the fastq read.
  *
  * @param sample the sample to which the read is attributed
  * @param stagger the length of the stagger sequence (generally 1-8bp)
  * @param randomBarcode the random-mer added during ligation that allows identification of PCR duplicates
  * @param targetBarcode the barcode or UMI that is associated with the target sequence in the BLT library
  * @param targetSequence the target sequence that was read (could be full length or cut to ~3bp)
  * @param cut true if the target appears to have been cut, false if it is full length
  * */
case class BltRead(sample: Sample, stagger: Byte, randomBarcode:  String, targetBarcode:  String, targetSequence: String, cut: Boolean)

/**
  * The results of de-duping BLT reads into unique observations.
  *
  * @param sample the Sample for which the observation is made
  * @param targetBarcode the targetBarcode sequence from the [[BltRead]]s
  * @param targetSequences the collection of target sequences (bases) observed with the targetBarcode in the sample
  * @param cut whether this represents a cut or uncut observation of the target
  * */
case class BltObservation(sample: Sample, targetBarcode: String, targetSequences: Array[String], cut: Boolean) {
  def readCount: Int = targetSequences.length
}

/* Case class to represent a mismatch or substitution in the guide vs. the target. */
case class MismatchTuple(pos: Int, guideBase:Char, targetBase: Char) {
  override def toString: String = s"($pos,$guideBase,$targetBase)"
}

/**
  * All of the information for a single target/targetBarcode across all samples.
  *
  * @param originalGuide the (common) guide sequence for the samples from the manifest
  * @param pam the (common) guide sequence for the samples from the manifest
  * @param targetSequence the validated target sequence (most common sequence meeting cutoffs)
  * @param targetBarcode the target barcode associated with the target sequence
  * @param observations the collection of observations of this target, cut and uncut, across all samples
  *                     sharing the same originalGuide and pam
  * @param annotations a set of annotations regarding how the target and guide sequences match
  * */
case class TargetInfo(originalGuide: String,
                      pam: String,
                      targetSequence: String,
                      targetBarcode: String,
                      observations: Array[BltObservation],
                      annotations: TargetAnnotation
                     ) {
  /** Returns the mismatch tuples if there are mismatches and no indels, otherwise empty. */
  def mismatchTuples: Seq[MismatchTuple] = {
    val len  = originalGuide.length
    val flip = observations.head.sample.enzyme.pamIs3PrimeOfTarget

    annotations.mismatchPositions.map { pos =>
      val posInStrings = if (flip) len - pos + 1 else pos
      val offset = posInStrings - 1
      MismatchTuple(pos, originalGuide.charAt(offset), targetSequence.charAt(offset))
    }
  }
}

/**
  * Computed annotations on how a target sequence differs from a guide.
  */
case class TargetAnnotation(cigar: Cigar, pamIs5PrimeOfTarget: Boolean) {
  /** The count of mismatches between the target and the guide. */
  def mismatches: Int = cigar.filter(_.operator == CigarOperator.X).map(_.length).sum

  /**
    * The count of inserted plus deleted bases.  If there is a single inserted base and a single
    * deleted base this will return 2, not 0.
    */
  def indelBases: Int = cigar.filter(_.operator.isIndel).map(_.length).sum

  /**
    * Returns the mismatch positions for any target that had no indels, otherwise empty. Mismatch
    * positions are expressed as 1 -based "distance from PAM".  I.e. if the guide base directly
    * adjacent to the PAM is a mismatch, the position is 1, the next base away is 2, etc.
    */
  def mismatchPositions: Seq[Int] = {
    if (cigar.exists(_.operator.isIndel)) {
      Seq.empty
    }
    else {
      var pos = 1
      val buffer = new ListBuffer[Int]
      cigar.foreach { elem =>
        if (elem.operator == CigarOperator.X) {
          forloop (from=0, until=elem.length) { i =>
            buffer += (pos + i)
          }
        }

        pos += elem.length
      }

      if (pamIs5PrimeOfTarget) {
        buffer
      }
      else {
        val len = cigar.lengthOnQuery
        buffer.map(p => len - p + 1).reverse
      }
    }
  }

  /**
    * The average mismatch position between the target and the guide. Only defined if there are no insertions
    * or deletions _and_ there is at least one mismatch.
    */
  def meanMismatchPosition: Option[Double] = mismatchPositions match {
    case Seq() => None
    case ps    => Some(ps.sum / ps.size.toDouble)
  }
}
