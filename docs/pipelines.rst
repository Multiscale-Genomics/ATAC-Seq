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

Pipelines
=========

ATAC Seq
---------
.. automodule:: process_atac_seq

   This pipeline takes two (or optionally one) fastq files as input along with the genome file indexed with Bowtie2 
   The pipeline uses trimgalore to trim adapter sequences, bowtie2 to align reads and produce bam 
   file and biobambam to filter the bam file. It uses macs2 to map, filter and produce a bed files
   that can be used to analyze and visualize Atac Seq data sets.

   Running from the command line
   =============================

   Parameters
   ----------
   config : file
      Location of the config file for the workflow
   in_metadata : file
      Location of the input list of files required by the process
   out_metadata : file
      Location of the output results.json file for returned files

   Returns
   -------
   peaks.xls : file
      Tabular file which contains information about called peaks. You can open it in excel and sort/filter using excel functions.         Information include:
      
      .. code-block:: none
         :linenos:

         chromosome name
         start position of peak
         end position of peak
         length of peak region
         absolute peak summit position
         pileup height at peak summit, -log10(pvalue) for the peak summit (e.g. pvalue =1e-10, then this value should be 10)
         fold enrichment for this peak summit against random Poisson distribution with local lambda, -log10(qvalue) at peak summit
         
   summits.bed : file
      Contains the peak summits locations for every peaks
      
   narrowPeak : file
      BED6+4 format file which contains the peak locations together with peak summit, pvalue and qvalue

   gappedPeak : file
      BED12+3 format which contains both the broad region and narrow peaks

   broadPeak : file
      BED6+3 format which is similar to narrowPeak file, except for missing the 10th column for annotating 
      peak summits


   Example
   -------
   When using a local verion of the [COMPS virtual machine](http://www.bsc.es/computer-sciences/grid-computing/comp-superscalar/downloads-and-documentation):

   .. code-block:: none
      :linenos:

      cd /home/compss/code/ATAC-Seq/atac_seq
      runcompss --lang=python process_atac_seq.py --config /home/compss/code/ATAC-Seq/tool_config/process_atac_seq.json --in_metadata /home/compss/code/ATAC-Seq/atac_seq/tests/json/input_atac_seq.json --out_metadata /home/compss/code/ATAC-Seq/atac_seq/tests/results_atac_seq.json

   Methods
   =======
   .. autoclass:: process_atac_seq.atacSeq
      :members:
      
      
