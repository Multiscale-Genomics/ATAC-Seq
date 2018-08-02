"""
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
"""

from __future__ import print_function

import os
import shlex
import subprocess
import sys

from utils import logger

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import compss_wait_on  # pylint: disable=ungrouped-imports

from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

from mg_process_fastq.tool.aligner_utils import alignerUtils
from mg_process_fastq.tool.bowtie_indexer import bowtieIndexerTool
from mg_process_fastq.tool.bowtie_aligner import bowtie2AlignerTool
from mg_process_fastq.tool.macs2 import macs2


# ------------------------------------------------------------------------------

class atacSeqTool(Tool):  # pylint: disable=invalid-name
    """
    Tool for running pipeline over a ATAC seq data
    """

    def __init__(self, configuration=None):
        """
        Initialise the tool with its configuration.


        Parameters
        ----------
        configuration : dict
            a dictionary containing parameters that define how the operation
            should be carried out, which are specific to each Tool.
        """
        logger.info("ATAC Seq")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @task(genome_file=FILE_IN, index_loc=FILE_OUT)
    def atac_seq(self, genome_file, input_fastq1, input_fastq2, output_narrowpeak, output_summits,
                 output_broadpeak, output_gappedpeak):  # pylint: disable=unused-argument, no-self-use
        """
        Atac Seq

        Parameters
        ----------
        genome_file : str
            Location of the genome assembly FASTA file
        input_fastq1 : str
            Location of the fastq file 1
        input_fastq2 : str
            Location of the fastq file 2
        """
        
        """
        Index the genome file
        """
        
        bti_input_file = {
        "genome": genome_file
        }

        bti_output_files = {
        'index': genome_file + ".bt2.tar.gz"
        }

        metadata = {
            "genome": Metadata(
                "Assembly", "fasta", genome_file, None,
                {'assembly': 'atac'}),
        }

        bti = bowtieIndexerTool()
        bti.run(bti_input_file, metadata, bti_output_files)

        """
        Trim Adapters from fastq files. 
        """
        command_line = "cutadapt -o " +  (input_fastq1 +".trimmed ") + input_fastq1
        logger.info("Cutadapt: command_line: " + command_line)

        try:
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()
        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0}) - Cutadapt: {1}\n{2}".format(
                msg.errno, msg.strerror, command_line))
            return False

        command_line = "cutadapt -o " +  (input_fastq2 +".trimmed ") + input_fastq2
        logger.info("Cutadapt: command_line: " + command_line)

        try:
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()
        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0}) - Cutadapt: {1}\n{2}".format(
                msg.errno, msg.strerror, command_line))
            return False
        
        
        """
        Align the genome file
        """
        if ".fa" in genome_file :
            index_file = genome_file.replace("fa","idx")
        
        elif ".fasta" in genome_file :
            index_file = genome_file.replace("fasta","idx")

        command_line = "bowtie2-build " + genome_file + " " + index_file 
        logger.info("Atac-Seq, bowtie2-build: command_line: " + command_line)
        
        resource_path = os.path.join(os.path.dirname(__file__), "data/")
        genome_fa = resource_path + "bsSeeker.Mouse.GRCm38.fasta"

        fastq_file_1 = input_fastq1 #+".trimmed"
        fastq_file_2 = input_fastq2 #+".trimmed"

        input_files = {
            "genome": genome_fa,
            "index": genome_fa + ".bt2.tar.gz",
            "loc": fastq_file_1,
            "fastq2": fastq_file_2
        }

        output_files = {
            "output": fastq_file_1.replace(".fastq", "_bt2.bam")
        }

        metadata = {
            "genome": Metadata(
                "Assembly", "fasta", genome_fa, None,
                {"assembly": "test"}),
            "index": Metadata(
                "index_bwa", "", [genome_fa],
                {
                    "assembly": "test",
                    "tool": "bwa_indexer"
                }
            ),
            "loc": Metadata(
                "data_wgbs", "fastq", fastq_file_1, None,
                {"assembly": "test"}
            ),
            "fastq2": Metadata(
                "data_wgbs", "fastq", fastq_file_2, None,
                {"assembly": "test"}
            )
        }

        bowtie2_handle = bowtie2AlignerTool()
        bowtie2_handle.run(input_files, metadata, output_files)


        bam_file = input_fastq1.replace("_1.fastq",".sorted.bam")
#        command_line += " | samtools view -u - | samtools sort - >" + bam_file
#        logger.info("Atac-Seq alignment: command_line: " + command_line)

#        try:
#            args = shlex.split(command_line)
#            process = subprocess.Popen(args)
#            process.wait()
#        except (IOError, OSError) as msg:
#            logger.fatal("I/O error({0}) - Atac-Seq alignment: {1}\n{2}".format(
#                msg.errno, msg.strerror, command_line))
#            return False
        
        """
        Make bed file
        
        home = os.path.expanduser('~')
        samTobed_path = home + "/lib/ATAC_Seq/ATAC-seq/atacseq/SAMtoBED.py"
        bed_file = input_fastq1.replace("_1.fastq","bed")
        
        command_line = "samtools view -h " + bam_file + " | " + samTobed_path +" -i - -o " + bed_file 
        command_line += " -x -v"
        logger.info("Atac-Seq bed file: command_line: " + command_line)

        try:
            args = shlex.split(command_line)
            process = subprocess.Popen(args)
            process.wait()
        except (IOError, OSError) as msg:
            logger.fatal("I/O error({0}) - Atac-Seq bed file: {1}\n{2}".format(
                msg.errno, msg.strerror, command_line))
            return False
        """
        
        """
        Call peaks with macs2
        """
        
        resource_path = os.path.join(os.path.dirname(__file__), "data/")

        input_files = {
            "bam": bam_file
        }
    
        output_files = {
            "narrow_peak": output_narrowpeak,
            "summits": output_summits,
            "broad_peak": output_broadpeak,
            "gapped_peak": output_gappedpeak
        }
    
        metadata = {
            "bam": Metadata(
                "data_atacseq", "fastq", [], None,
                {'assembly' : 'test'}),
        }

        print("BAM FILE :", resource_path + bam_file)    
        print("Output_files :",output_files["narrow_peak"])    
        print("output_narrowPeak :", output_narrowpeak)    
        
        macs_handle = macs2({"macs_nomodel_param": True})
        macs_handle.run(input_files, metadata, output_files)

        return True
        

    def run(self, input_files, input_metadata, output_files):
        """
        Tool for generating bed and peak files for use with the ATAC-Seq
        data

        Parameters
        ----------
        input_files : dict
            genome : string 
                Location of the genome fasta file 
            fastq1 : string
                Location of the FASTQ file
            fastq2 : string
                Location of the paired end FASTQ file
        input_metadata : list

        Returns
        -------
        output_files = {
            "narrow_peak": resource_path + "atacseq.Human.ERR1659027_peaks.narrowPeak",
            "summits": resource_path + "atacseq.Human.ERR1659027_peaks.summits.bed",
            "broad_peak": resource_path + "atacseq.Human.ERR1659027_peaks.broadPeak",
            "gapped_peak": resource_path + "atacseq.Human.ERR1659027_peaks.gappedPeak"
        }    
        """

        # input and output share most metadata
        results = self.atac_seq(
            input_files['genome'],
            input_files['fastq1'],
            input_files['fastq2'],
            output_files['narrow_peak'],
            output_files['summits'],
            output_files['broad_peak'],
            output_files['gapped_peak']
        )
        
        print("Summits file: ", output_files["summits"])
        
        results = compss_wait_on(results)

        if results is False:
            logger.fatal("ATAC Seq: run failed")
            return {}, {}
        
        output_bed_file = output_files["summits"] 
        output_narrowPeak_file = output_files["narrow_peak"]
        
        output_metadata = {
            "narrow_peak": Metadata(
                data_type="atacseq",
                file_type="narrowpeak",
                file_path=output_files["narrow_peak"],
                #sources=[input_metadata["narrowpeak"].file_path],
                #taxon_id=input_metadata["narrowpeak"].taxon_id,
                meta_data={
                    "tool": "atac_seq"
                }
            ),
                           
            "summits": Metadata(
                data_type="atacseq",
                file_type="summits",
                file_path=output_files["summits"],
                #sources=[input_metadata["summits"].file_path],
                #taxon_id=input_metadata["summits"].taxon_id,
                meta_data={
                    "tool": "atac_seq"
                }
            )
        }
        

        return (output_files, output_metadata)

# ------------------------------------------------------------------------------
