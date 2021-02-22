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

object Sequences {
  private val AllBases         = "ACGTRYSWKMBDHVN".getBytes
  private val UnambiguousBases = "ACGT".getBytes

  /** Turns an Array[Byte] into a String for printing. */
  private def str(bytes: Array[Byte]) = new String(bytes)

  /** Returns true if two bases are equal and false otherwise. */
  def areEqual(a: Byte, b: Byte): Boolean = a == b

  /** Computes the number of mismatches between the two arrays. */
  def mismatches(lhs: Array[Byte], rhs: Array[Byte]): Int = {
    require(lhs.length == rhs.length, s"Sequences must be the same length: ${str(lhs)} and ${str(rhs)}.")
    mismatches(lhs, 0, rhs, 0, lhs.length)
  }

  /**
    * Counts the number of mismatches in the two regions of two arrays. If a `max` is provided
    * and the regions contain `max` or more mismatches, the comparison will terminate early and
    * return `max`.
    *
    * @param lhs      the first sequence to compare
    * @param lhsStart the starting index (0-based) within the first array to compare
    * @param rhs      the second sequence to compare
    * @param rhsStart the starting index (0-based) within the second array to compare
    * @param length   the number of bases to compare
    * @param max      the maximum number of mismatches, beyond which the total is unimportant
    */
  def mismatches(lhs: Array[Byte], lhsStart: Int, rhs: Array[Byte], rhsStart: Int, length: Int, max: Int = Int.MaxValue): Int = {
    require(lhsStart + length <= lhs.length, s"Cannot compare $length bases from $lhsStart given sequence ${str(lhs)}")
    require(rhsStart + length <= rhs.length, s"Cannot compare $length bases from $rhsStart given sequence ${str(rhs)}")

    var count = 0
    forloop (from=0) (i => i < length && count < max) (i => i+1) { i =>
      if (!areEqual(lhs(lhsStart+i), rhs(rhsStart+i))) count += 1
    }

    count
  }

  /** Returns true f `b` is found within `bs`. */
  private def contains(bs: Array[Byte], b: Byte): Boolean = {
    var found = false
    forloop (from=0) (!found && _ < bs.length) (_ + 1) { i => found = bs(i) == b }
    found
  }

  /** Checks to see if a base is a valid base. */
  def isValidBase(b: Byte, allowAmbiguityCodes: Boolean = false): Boolean =
    if (allowAmbiguityCodes) contains(AllBases, b)
    else contains(UnambiguousBases, b)

  /** Checks to see if all bytes in the array are valid bases. */
  def areValidBases(bs: Array[Byte], allowAmbiguityCodes: Boolean = false): Boolean = {
    val bases = if (allowAmbiguityCodes) AllBases else UnambiguousBases
    bs.forall(b => contains(bases, b))
  }
}
