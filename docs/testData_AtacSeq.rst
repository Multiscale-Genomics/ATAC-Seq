.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Test Data for ATAC-seq pipeline
================================

The following document is for the preparation of data set required for testing
the ATAC-seq pipeline. The document has been written with macOS Sierra in mind,
although many of the commands are cross platform (\*nix) compliant.

You would need to have the tools listed in "Prerequisites" installed on your system.
For more details on installing the tools for this pipeline please refer to

http://multiscale-genomics.readthedocs.io/projects/mg-process-fastq/en/latest/full_installation.html

If you already have certain packages installed feel free to skip over certain
steps. Likewise the bin, lib and code directories are relative to the home dir;
if this is not the case for your system then make the required changes when
running these commands.

Prerequisites
-------------

   - Cutadapt
   - Trim Galore
   - Bowtie2
   - Samtools
   - macs2

Data set for genome file
------------------------

Filtering for required coverage
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Download the genome file from

.. code-block:: none

   wget "http://www.ebi.ac.uk/ena/data/view/CM000663.2,CM000664.2,CM000665.2,CM000666.2,CM000667.2,CM000668.2,CM000669.2,CM000670.2,CM000671.2,CM000672.2,CM000673.2,CM000674.2,CM000675.2,CM000676.2,CM000677.2,CM000678.2,CM000679.2,CM000680.2,CM000681.2,CM000682.2,CM000683.2,CM000684.2,CM000685.2,CM000686.2,J01415.2&display=fasta&download=fasta&filename=entry.fasta" -O atac.Human.fasta


Checkout https://github.com/Multiscale-Genomics/mg-misc-scripts/blob/master/ATACSeq_Scripts/extract_chromosomeForATAC.py and extract chromosome 17 from the above file using the following command.

.. code-block:: none

   python extract_chromosomeForATAC.py path/to/your/input/file path/to/output/file

Download the fastq files from

.. code-block:: none

   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR165/007/ERR1659027/ERR1659027_1.fastq.gz
   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR165/007/ERR1659027/ERR1659027_2.fastq.gz

Index the fasta file

.. code-block:: none

   bowtie2-build atac.Human.fasta atac.Human.hg19

Align the fastq files

.. code-block:: none

   bowtie2 -x atac.Human.hg19 -U ERR1659027_1.fastq ERR1659027_2.fastq -S atac.Human.hg19.sam

Filter out the aligned reads from the above sam file.

.. code-block:: none

   awk '$3 != "*"' atac.Human.hg19.sam >atac.Human.hg19.filtered.sam

Sort the sam file

.. code-block:: none

   samtools sort atac.Human.hg19.filtered.sam >atac.Human.hg19.filtered.sorted.sam

Find the depths of coverage from the sorted file

.. code-block:: none

   samtools depth atac.Human.hg19.filtered.sorted.sam >atac.Human.hg19.dp

From the depth file, find regions with >= 70 depth, spanning over >=55 base pairs.
You may get the script for this from:
https://github.com/Multiscale-Genomics/mg-misc-scripts/blob/master/MNaseSeq_Scripts/traverseForCoverageRegion_MNase.py

Run it using:

.. code-block:: none

   python traverseForCoverageRegion_MNase.py path/to/atac.Human.hg19.dp

Running this script would print the spanning regions. Running this script for this data set gives multiple regions. The output is in the format : start - end - depth.  The one at the end has a maximal coverage from this data set. Since it is a continuous region, you may take the first starting base pair and the last ending base pair, as inputs for the next step. (Take out 1000 and add in 1000 to these respectively to get upstream and downstream spanning bases)

Extract the corresponding fasta sequence from the chromosome file for the positions retrieved from the above step. Checkout file from https://github.com/Multiscale-Genomics/mg-misc-scripts/blob/master/MNaseSeq_Scripts/extractChromosomalRegion.py and run using command:

.. code-block:: none

   python extractChromosomalRegion.py path/to/original/fasta/file path/to/output/file/for/region/atac.Human.hg19.fasta starting_base_position ending_base_position

Making the Fastq file
^^^^^^^^^^^^^^^^^^^^^^

Index the fasta file for the selected region

.. code-block:: none

   bowtie2-build atac.Human.hg19.fasta atac.Human.hg19

Align the fastq files

.. code-block:: none

   bowtie2 -x atac.Human.hg19 -U ERR1659027_1.fastq ERR1659027_2.fastq -S atac.Human.hg19.sam

Filter this sam file for the reads which aligned with chromosome 17 using the following command:

.. code-block:: none

   awk '$3 != "*"' atac.Human.hg19.sam >atac.Human.hg19.filtered.sam

From the filtered reads from the above output file, extract the corresponding entries in fastq file. You may do this using the file at https://github.com/Multiscale-Genomics/mg-misc-scripts/blob/master/ATACSeq_Scripts/makeFastQFiles.py

and running it via command line:

.. code-block:: none

   python makeFastQFiles.py --samfile path/to/atac.Human.hg19.filtered.sam --fastQfile /path/to/ERR1659027_1.fastq --pathToOutput /path/to/save/output/fastq/file/to/ --fastqOut ERR1659027_1_atac.fastq
      python makeFastQFiles.py --samfile path/to/atac.Human.hg19.filtered.sam --fastQfile /path/to/ERR1659027_2.fastq --pathToOutput /path/to/save/output/fastq/file/to/ --fastqOut ERR1659027_2_atac.fastq

The fastq files in the above step and fasta file atac.Human.hg19.fasta together make up the data set for MNase-seq pipeline
