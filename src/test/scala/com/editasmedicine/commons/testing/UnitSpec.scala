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
package com.editasmedicine.commons.testing

import java.nio.file.{Files, Path}
import java.util.stream.Collectors

import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.commons.io.PathUtil
import org.scalatest.{FlatSpec, Matchers}

/** Base class for unit tests. */
class UnitSpec extends FlatSpec with Matchers {
  /** Creates a new temp file for use in testing that will be deleted when the VM exits. */
  protected def makeTempFile(prefix: String, suffix: String) : Path = {
    val path = Files.createTempFile(prefix, suffix)
    path.toFile.deleteOnExit()
    path
  }

  /** Creates a new temporary directory for use. */
  protected def makeTempDir(prefix: String, suffix: String) : Path = {
    val path = Files.createTempFile(prefix, suffix)
    Files.delete(path)
    Files.createDirectory(path)
    path.toFile.deleteOnExit()
    path
  }

  /** Deletes a path/file, recursing down if the path is directory. */
  protected def deletePath(path: Path): Unit = {
    if (Files.isDirectory(path)) {
      val childStream = Files.list(path)
      val children = childStream.collect(Collectors.toList())
      childStream.close()
      children.foreach(deletePath)
    }

    Files.deleteIfExists(path)
  }

  /** Executes a block of code with access to a temp directory that is cleaned up after the block. */
  protected def withTempDir[A](body: DirPath => A): A = {
    val dir = makeTempDir("testing.", ".dir")
    try { body(dir) }
    finally { deletePath(dir) }
  }

  /** Finds a resource on the classpath and turns it into a Path. */
  protected def testResource(name: String, pkg: Package = getClass.getPackage) : Path = {
    val packagePath = pkg.getName.replace('.', '/')
    val url = getClass.getClassLoader.getResource(packagePath + "/" + name)
    PathUtil.pathTo(url.getPath)
  }
}
