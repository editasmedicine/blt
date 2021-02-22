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

import com.editasmedicine.commons.clp.{ClpGroups, EditasTool}
import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.fastq.FastqSource
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, ProgressLogger}

@clp(group=ClpGroups.Blt, description=
  """
    |Analyzes a BLT experiment starting from one or more fastq files.  Takes in a set of fastq files and a
    |sample manifest then proceeds to:
    |  1. Demultiplex the reads and extract relevant sequences from them
    |  2. Validate the set of targets and target barcodes (UMIs) with sufficient uncut observations
    |  3. Generate information on the target observations for each sample
    |  4. Emits various text and PDF files with data about the experiment
    |
    |Within the output folder a number of top level summary files are created:
    |  1. demultiplexing.summary.txt - summary information about the number of reads extracted vs. discarded
    |  2. demultiplexing.details.txt - information on extracted/discarded reads per sample
    |  3. summary.txt - high level summary for each sample including the score for each guide
    |  4. cut_rate_by_mismatches.pdf - plots of cut rate by #mismatches for all samples
    |
    |In addition a directory per sample is created with the following files:
    |  1. {sample}.umis.txt.gz - detailed infomration about each target barcode/umi observed
    |  2. {sample}.targets.txt.gz - similar to {sample}.umis.txt.gz but rolled up by target sequence
    |  3. {sample}.summary.txt - summary information by #mismatches
    |  4. {sample}.pdf - (cut samples only) a plot of the distribution of cut rates by #mismatches
  """)
class AnalyzeExperiment
( @arg(flag='i', doc="Input fastq file(s), optionally gzipped.") val input: Seq[PathToFastq],
  @arg(flag='s', doc="Sample manifest file.") val sampleManifest: FilePath,
  @arg(flag='o', doc="Output directory") val output: DirPath,
  @arg(flag='m', doc="Maximum mismatches allowed between expected and observed sample barcodes.") val maxMismatches: Int = 2,
  @arg(flag='d', doc="Minimum distance (in mismatches) between sample barcodes.")    val minDistance: Int = 2,
  @arg(flag='q', doc="Minimum average quality of extracted bases from a read.")      val minQuality: Int = 20,
  @arg(flag='u', doc="Minimum number of un-cut observations to keep a UMI/target.")  val minUncutReads: Int = Defaults.MinUncutReadsPerTarget,
  @arg(flag='f', doc="Minimum fraction of reads for a UMI that must be identical.")  val minIdenticalFraction: Double = Defaults.MinFractionIdenticalPerTarget,
  @arg(flag='c', doc="Include cut samples when validating targets/umis.")            val useCutSamplesInValidation: Boolean = Defaults.UseCutSamplesInValidation,
  @arg(flag='l', doc="Length if guides/targets were padded to a fixed length.")      val fixedGuideLength: Option[Int] = None,
  @arg(flag='t', doc="Number of threads to use.") val threads: Int = 4
) extends EditasTool {
  Io.assertReadable(input)
  Io.mkdirs(output)
  Io.assertWritableDirectory(output)
  validate(maxMismatches >= 0, "--max-mismatches must be a positive number")
  validate(minDistance   >= 0, "--min-distance must be >= 1")
  validate(minUncutReads >= 1, "--min-uncut-reads must be >= 1")

  override def execute(): Unit = {
    val fqIterator = input.iterator.flatMap(FastqSource(_))
    val manifest   = SampleManifest(sampleManifest)
    val demuxer    = Demultiplexer(manifest, maxMismatches=maxMismatches, minDistance=minDistance)
    val extractor  = new Cas9ReadExtractor(demultiplexer=demuxer, minQuality=minQuality, fixedGuideLength=fixedGuideLength)

    // Validate that the fixed guide length is >= all the actual guide lengths if it was supplied
    for (len <- fixedGuideLength; sample <- manifest) {
      require(len >= sample.guide.length, s"Sample ${sample.name} had a guide (${sample.guide}) longer than the specified fixed guide length.")
    }

    // Extract the BltReads
    logger.info("Extracting information from fastq reads.")
    val progress = ProgressLogger(logger, noun="reads", verb="Processed", unit=2.5e6.toInt)
    val reads = fqIterator.flatMap { fq => yieldAndThen(extractor.extract(fq)) { progress.record() } }.toBuffer
    extractor.printStats(logger)

    // Output demultiplexing metrics
    Metric.write(output.resolve("demultiplexing.summary.txt"), Seq(extractor.summaryDemuxMetrics))
    Metric.write(output.resolve("demultiplexing.details.txt"), extractor.sampleDemuxMetrics)

    // Run the main analysis
    val analyzer = new BltAnalysis(
      manifest                  = manifest,
      minUncutReads             = minUncutReads,
      minIdenticalFraction      = minIdenticalFraction,
      useCutSamplesInValidation = useCutSamplesInValidation,
      threads                   = threads,
      output                    = output)

    analyzer.analyze(reads)
  }
}
