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


object Enzyme {
  /** All the Enzymes that are currently handled. */
  val values = Seq(Cas9)

  /** Returns the enzyme for the name. Throws an exception if the name is not valid. */
  def apply(s: String): Enzyme = {
    values.find(e => e.toString.toLowerCase == s.trim.toLowerCase).getOrElse(throw new IllegalArgumentException("Unrecognized enzyme: " + s))
  }
}

/** Trait representing the different kinds of cutting enzymes that are handled. */
sealed trait Enzyme extends Product {
  /** To String that ensures the result is the enzyme name. */
  override final def toString: String = productPrefix

  /** Returns the read extractor to use with this enzyme. */
  def readExtractor(demultiplexer: Demultiplexer, minQuality: Int = 20, fixedGuideLength: Option[Int] = None): ReadExtractor

  /** True if the PAM is 5' of the target, false if it is 3' of the target sequence. */
  def pamIs5PrimeOfTarget: Boolean

  /** True if the PAM is 3' of the target, false if it is 5' of the target sequence. */
  final def pamIs3PrimeOfTarget: Boolean = !pamIs5PrimeOfTarget
}

/** Enzyme object representing Cas9. */
case object Cas9 extends Enzyme {
  override val pamIs5PrimeOfTarget: Boolean = false

  /** Returns the read extractor to use with this enzyme. */
  override def readExtractor(demultiplexer: Demultiplexer, minQuality: Int, fixedGuideLength: Option[Int]): ReadExtractor = {
    new Cas9ReadExtractor(demultiplexer, minQuality, fixedGuideLength)
  }
}

