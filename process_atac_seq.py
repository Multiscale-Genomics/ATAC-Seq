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
import argparse
import subprocess

from utils import logger

from basic_modules.metadata import Metadata
from basic_modules.workflow import Workflow

from mg_process_fastq.tool.bowtie_aligner import bowtie2AlignerTool
from mg_process_fastq.tool.biobambam_filter import biobambam
from mg_process_fastq.tool.trimgalore import trimgalore
from mg_process_fastq.tool.macs2 import macs2


# ------------------------------------------------------------------------------

class process_atac_seq(Workflow):  # pylint: disable=invalid-name
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

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
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
        bedpe :
            Users can specify this if the results need some changes.
            Or if the data contains singleton reads

        Returns
        -------
        output_files = {
            "narrow_peak": resource_path + "atacseq.Human.ERR1659027_peaks.narrowPeak",
            "summits": resource_path + "atacseq.Human.ERR1659027_peaks.summits.bed",
            "broad_peak": resource_path + "atacseq.Human.ERR1659027_peaks.broadPeak",
            "gapped_peak": resource_path + "atacseq.Human.ERR1659027_peaks.gappedPeak"
        }

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

        bedpe = metadata['bedpe']
        genome_file = input_files['genome']
        input_fastq1 = input_files['fastq1']

        if "fastq2" in input_files:
            input_fastq2 = input_files['fastq2']
        output_narrowpeak = output_files['narrow_peak']
        output_summits = output_files['summits']
        output_broadpeak = output_files['broad_peak']
        output_gappedpeak = output_files['gapped_peak']

        results = {}
        resource_path = os.path.join(os.path.dirname(__file__), "data/")

        # Trim Adapters from fastq files.
        
        fastq1_trimmed = input_fastq1 + '.trimmed'
        fastq2_trimmed = input_fastq2 + '.trimmed'

        if "fastq2" in input_files:
            files = {
                'fastq1': input_fastq1,
                'fastq2': input_fastq2
            }

            files_out = {
                "fastq1_trimmed": fastq1_trimmed,
                "fastq2_trimmed": fastq2_trimmed,
                "fastq1_report": input_fastq1 + '.trimmed.report.txt',
                "fastq2_report": input_fastq2 + '.trimmed.report.txt'
            }

        else:
            files = {
                'fastq1': input_fastq1
            }

            files_out = {
                "fastq1_trimmed": fastq1_trimmed,
                "fastq1_report": input_fastq1 + '.trimmed.report.txt',
            }

        self.configuration["tg_length"] = "0"

        tg_handle = trimgalore(self.configuration)
        tg_handle.run(files, metadata, files_out)

        #Align the genome file

        fastq_file_1 = input_fastq1

        if "fastq2" in input_files:
            fastq_file_2 = input_fastq2

            input_files = {
                "genome": genome_file,
                "index": genome_file + ".bt2.tar.gz",
                "loc": fastq_file_1,
                "fastq2": fastq_file_2
            }

            metadata_bowtie = {
                "genome": Metadata(
                    "Assembly", "fasta", genome_file, None,
                    {"assembly": "test"}),
                "loc": Metadata(
                    "data_atac", "fastq", fastq_file_1, None,
                    {"assembly": "test"}
                )
            }

        else:
            input_files = {
                "genome": genome_file,
                "index": genome_file + ".bt2.tar.gz",
                "loc": fastq_file_1
            }

            metadata_bowtie = {
                "genome": Metadata(
                    "Assembly", "fasta", genome_file, None,
                    {"assembly": "test"}),
                "loc": Metadata(
                    "data_atac", "fastq", fastq_file_1, None,
                    {"assembly": "test"}
                )
            }

        output_files = {
            "output": fastq_file_1.replace(".fastq", "_bt2.bam")
        }

        bowtie2_handle = bowtie2AlignerTool()
        bowtie2_handle.run(input_files, metadata_bowtie, output_files)

        bam_file = fastq_file_1.replace(".fastq", "_bt2.bam")#input_fastq1.replace("_1.fastq", ".bam")

        bam_filtered = bam_file + "filtered"
        input_files = {
            "input": bam_file
        }
        
        print ("******** Input file for BAM *******", bam_file)

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

        # Call peaks with macs2

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

        if bedpe:
            self.configuration["macs_format_param"] = "BEDPE"

        macs_handle = macs2(self.configuration)
        macs_handle.run(input_files, metadata, output_files)

        print("Summits file: ", output_files["summits"])

        #results = compss_wait_on(results)

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

def main_json(config, in_metadata, out_metadata):
    """
    Alternative main function
    -------------

    This function launches the app using configuration written in
    two json files: config.json and input_metadata.json.
    """
    # 1. Instantiate and launch the App
    print("1. Instantiate and launch the App")
    from apps.jsonapp import JSONApp
    app = JSONApp()
    result = app.launch(process_atac_seq,
                        config,
                        in_metadata,
                        out_metadata)

    # 2. The App has finished
    print("2. Execution finished; see " + out_metadata)
    print(result)

    return result


# ------------------------------------------------------------------------------

if __name__ == "__main__":

    # Set up the command line parameters
    PARSER = argparse.ArgumentParser(description="ATAC-Seq")
    PARSER.add_argument("--config", help="Configuration file")
    PARSER.add_argument("--in_metadata", help="Location of input metadata file")
    PARSER.add_argument("--out_metadata", help="Location of output metadata file")
    PARSER.add_argument("--local", action="store_const", const=True, default=False)

    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    CONFIG = ARGS.config
    IN_METADATA = ARGS.in_metadata
    OUT_METADATA = ARGS.out_metadata
    LOCAL = ARGS.local

    if LOCAL:
        import sys
        sys._run_from_cmdl = True  # pylint: disable=protected-access

    RESULTS = main_json(CONFIG, IN_METADATA, OUT_METADATA)

    print(RESULTS)
