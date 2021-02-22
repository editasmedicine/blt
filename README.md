# BLT - Barcoded Library Of Targets

This repository contains software for analyzing data generated using the Barcoded Library of Targets (BLT) protocol.

## License


## Pre-requisites

In order to build and run the BLT software you will need:

1. A Java Development Kit v8 or higher available from [the Oracle Java site ](http://www.oracle.com/technetwork/java/javase/downloads/index.html).
2. A working installation of the Scala Built Tool aka _sbt_. Instructions for installing sbt can be found on the [SBT website](https://www.scala-sbt.org/download.html).
3. A working internet connection which SBT will use to automatically download dependencies
1. R version 3.2 or higher available from [the R website](https://cloud.r-project.org/)
2. R package `ggplot2`. ggplot2 can be installed from within R, or by running `echo 'install.packages("ggplot2", repos="http://cran.us.r-project.org", dependencies=TRUE)' | R --no-save`

## Building

The software has been built and tested on Linux and Mac OS X.  The software should build and run on any operating system that supports Java >= 8 and SBT, but other platforms have not been tested.

To compile and test the software run the following in the same directory as this README file:

```
sbt clean test
```

This process may take several minutes the first time it is run as sbt will download needed software libraries from the internet.  Repeated executions should take 15-30s.  A successful run of the tests will result in two lines similar to the following being printed before sbt exits:

```
[info] All tests passed.
[info] Passed: Total 69, Failed 0, Errors 0, Passed 69
[success] Total time: 20 s, completed Feb 8, 2018 1:44:01 PM
```

The BLT software is compiled into a standalone [JAR](https://en.wikipedia.org/wiki/JAR_(file_format)) file for use.  To build the JAR file run:

```
sbt assembly
```

As with the tests, the last line printed should contain `[success]`, e.g.: 

```
[success] Total time: 17 s, completed Feb 8, 2018 1:53:53 PM
```

This will produce a JAR file located at `./target/scala-2.12/blt.jar`.

## Running

The software is run by invoking the following command with appropriate options:

```
java -Xmx8g -jar target/scala-2.12/blt.jar AnalyzeExperiment
```

For example, in a directory containing multiple fastqs for analysis and a file called `samples.txt` describing the samples within the experiment the following command could be used:

```
java -Xmx8g -jar target/scala-2.12/blt.jar AnalyzeExperiment -s samples.txt -i *.fastq.gz
```

Runtime is highly dependent on a) the amount of input data, b) the quality of input data, c) the diversity of the BLT library.  For high quality data runtime is on the order of 10-30 seconds per million reads in the fastq file. 

### AnalyzeExperiment Usage

Analyzes a BLT experiment starting from one or more fastq files.  Takes in a set of fastq files and a
sample manifest then proceeds to:

  1. Demultiplex the reads and extract relevant sequences from them
  2. Validate the set of targets and target barcodes (UMIs) with sufficient uncut observations
  3. Generate information on the target observations for each sample
  4. Emits various text and PDF files with data about the experiment

Within the output folder a number of top level summary files are created:

  1. demultiplexing.summary.txt - summary information about the number of reads extracted vs. discarded
  2. demultiplexing.details.txt - information on extracted/discarded reads per sample
  3. summary.txt - high level summary for each sample including the score for each guide
  4. cut\_rate\_by\_mismatches.pdf - plots of cut rate by #mismatches for all samples

In addition a directory per sample is created with the following files:

  1. {sample}.umis.txt.gz - detailed infomration about each target barcode/umi observed
  2. {sample}.targets.txt.gz - similar to {sample}.umis.txt.gz but rolled up by target sequence
  3. {sample}.summary.txt - summary information by #mismatches
  4. {sample}.pdf - (cut samples only) a plot of the distribution of cut rates by #mismatches

### AnalyzeExperiment Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToFastq|Input fastq file(s), optionally gzipped.|Required|Unlimited||
|sample-manifest|s|FilePath|Sample manifest file.|Required|1||
|output|o|DirPath|Output directory|Required|1||
|max-mismatches|m|Int|Maximum mismatches allowed between expected and observed sample barcodes.|Optional|1|2|
|min-distance|d|Int|Minimum distance (in mismatches) between sample barcodes.|Optional|1|2|
|min-quality|q|Int|Minimum average quality of extracted bases from a read.|Optional|1|20|
|min-uncut-reads|u|Int|Minimum number of un-cut observations to keep a UMI/target.|Optional|1|3|
|min-identical-fraction|f|Double|Minimum fraction of reads for a UMI that must be identical.|Optional|1|0.9|
|use-cut-samples-in-validation|c|Boolean|Include cut samples when validating targets/umis.|Optional|1|false|
|fixed-guide-length|l|Int|Length if guides/targets were padded to a fixed length.|Optional|1||
|threads|t|Int|Number of threads to use.|Optional|1|4|

### Sample Manifest Description

The sample manifest is a tab-delimited text file describing the samples to be analyzed.  An [example](./example_data/inputs/sample-manifest.txt) is provided in the [example_data](./example_data) directory.

|column name|description|example value(s)|
|-----------|-----------|-------------|
|sample|The name of the sample. Should be unique within the manifest.|"s1", "control"|
|sample_barcode|The molecular barcode associated with the sample.|CACATACGCACTACG|
|guide|The guide sequence associated with the sample.|ATGATGTATGCGCTATACGTCA|
|enzyme|The name of the enzyme used to cut.|Cas9|
|pam|The PAM sequence used in the assay.|GGGAGT|
|cut|Whether the sample is cut or an uncut control|yes, no, true, false|
|off\_target\_file|Optional path to a file of predicted off-target sequences.|s1-off-target.txt|

### Off-Target File Description

The off-target file describes the predicted off-target sites for a guide. Off-target information is used solely to highlight predicted off-target sequences when plotting.  The file should be in CSV format and should contain the following columns:

|column name|description|example value(s)|
|-----------|-----------|-------------|
|loc|The location of the off-target hit in the format `chr:pos`.|chr1:1234567|
|off_target|The genomic sequence at the off-target site. Must be only [ACGTagct].| ATGATcTATGCGCgAaACGTCA |

## Example Data

A small example data set is contained within the [example_data](example_data) directory.   The example dataset can be processed using the following command:

```
java -jar target/scala-2.12/blt.jar AnalyzeExperiment \
  -i example_data/inputs/reads.fastq.gz \
  -s example_data/inputs/sample-manifest.txt \
  --fixed-guide-length=23 \
  -o example_analysis_outputs \
  --min-uncut-reads=1 \
  --use-cut-samples-in-validation=true 
```

The `sample-manifest.txt` contains relative paths to the off-target file(s); the command given above must be executed from the main project directory, or the paths in the manifest must be adjusted.

A set of outputs from running the above command on the example inputs is available in the [example_data/outputs](example_data/outputs) directory.
