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
from mg_process_fastq.tool.biobambam_filter import biobambam
from mg_process_fastq.tool.trimgalore import trimgalore
from mg_process_fastq.tool.macs2 import macs2


# ------------------------------------------------------------------------------

class atacSeq(Tool):  # pylint: disable=invalid-name
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

        genome_file = input_files['genome']
        input_fastq1 = input_files['fastq1']

        if "fastq2" in input_files:
            input_fastq2 = input_files['fastq2']
        output_narrowpeak = output_files['narrow_peak']
        output_summits = output_files['summits']
        output_broadpeak = output_files['broad_peak']
        output_gappedpeak = output_files['gapped_peak']

        """
            genome_file : str
                Location of the genome assembly FASTA file
            input_fastq1 : str
                Location of the fastq file 1
            input_fastq2 : str
                Location of the fastq file 2
            output_narrowpeak : str
                Location of the narrow peak output file
            output_summits : str
                Location of the summits.bed output file
            output_broadpeak : str
                Location of the broadpeak output file
            output_gappedpeak : str
                Location of the gappedpeak output file
        """
        results = {}
        resource_path = os.path.join(os.path.dirname(__file__), "data/")

        """
        Trim Adapters from fastq files.
        """

        if "fastq2" in input_files:
            files = {
                'fastq1': input_fastq1,
                'fastq2': input_fastq2
            }

            metadata = {
                "fastq1": Metadata(
                    "data_wgbs", "fastq", files['fastq1'], None,
                ),

                "fastq2": Metadata(
                    "data_wgbs", "fastq", files['fastq2'], None,
                )
            }

            files_out = {
                "fastq1_trimmed": input_fastq1 + '.trimmed',
                "fastq2_trimmed": input_fastq2 + '.trimmed',
                "fastq1_report": input_fastq1 + '.trimmed.report.txt',
                "fastq2_report": input_fastq2 + '.trimmed.report.txt'
            }

        else:
            files = {
                'fastq1': input_fastq1
            }

            metadata = {
                "fastq1": Metadata(
                    "data_wgbs", "fastq", files['fastq1'], None,
                )
            }

            files_out = {
                "fastq1_trimmed": input_fastq1 + '.trimmed',
                "fastq1_report": input_fastq1 + '.trimmed.report.txt',
            }

        self.configuration["tg_length"] = "0"

        tg_handle = trimgalore(self.configuration)
        tg_files, tg_meta = tg_handle.run(files, metadata, files_out)

        """
        Align the genome file
        """
        if ".fa" in genome_file:
            index_file = genome_file.replace("fa", "idx")

        elif ".fasta" in genome_file:
            index_file = genome_file.replace("fasta", "idx")

        fastq_file_1 = input_fastq1

        if "fastq2" in input_files:
            fastq_file_2 = input_fastq2

            input_files = {
                "genome": genome_file,
                "index": genome_file + ".bt2.tar.gz",
                "loc": fastq_file_1,
                "fastq2": fastq_file_2
            }

            metadata = {
                "genome": Metadata(
                    "Assembly", "fasta", genome_file, None,
                    {"assembly": "test"}),
                "index": Metadata(
                    "index_bwa", "", [genome_file],
                    {
                        "assembly": "test",
                        "tool": "bowtie_aligner"
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

        else:
            input_files = {
                "genome": genome_file,
                "index": genome_file + ".bt2.tar.gz",
                "loc": fastq_file_1
            }

            metadata = {
                "genome": Metadata(
                    "Assembly", "fasta", genome_file, None,
                    {"assembly": "test"}),
                "index": Metadata(
                    "index_bwa", "", [genome_file],
                    {
                        "assembly": "test",
                        "tool": "bowtie_aligner"
                    }
                ),
                "loc": Metadata(
                    "data_wgbs", "fastq", fastq_file_1, None,
                    {"assembly": "test"}
                )
            }

        output_files = {
            "output": fastq_file_1.replace(".fastq", "_bt2.bam")
        }

        bowtie2_handle = bowtie2AlignerTool()
        bowtie2_handle.run(input_files, metadata, output_files)

        bam_file = input_fastq1.replace("_1.fastq", ".bam")

        bam_filtered = bam_file + "filtered"
        input_files = {
            "input": bam_file
        }

        output_files = {
            "output": bam_filtered
        }

        metadata = {
            "input": Metadata(
                "data_chipseq", "fastq", [], None,
                {'assembly' : 'test'}),
        }

        bbb = biobambam()
        bbb.run(input_files, metadata, output_files)

        """
        Call peaks with macs2
        """

        input_files = {
            "bam": bam_filtered
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
        print("Output_files :", output_files["narrow_peak"])
        print("output_narrowPeak :", output_narrowpeak)

        self.configuration["macs_nomodel_param"] = True
        self.configuration["macs_keep-dup_param"] = "all"

        macs_handle = macs2(self.configuration)
        macs_handle.run(input_files, metadata, output_files)

        print("Summits file: ", output_files["summits"])

        results = compss_wait_on(results)

        if results is False:
            logger.fatal("ATAC Seq: run failed")
            return {}, {}

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
