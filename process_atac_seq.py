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

import argparse

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

        bedpe = self.configuration['macs2_bedpe']
        genome_file = input_files['genome']

        input_fastq1 = input_files['fastq1']
        fastq1_trimmed = input_fastq1 + '.trimmed'

        output_narrowpeak = output_files['narrow_peak']
        output_summits = output_files['summits']
        output_broadpeak = output_files['broad_peak']
        output_gappedpeak = output_files['gapped_peak']

        results = {}

        # Trim Adapters from fastq files.
        input_files_trim = {
            'fastq1': input_fastq1
        }

        files_out = {
            "fastq1_trimmed": fastq1_trimmed,
            "fastq1_report": input_fastq1 + '.trimmed.report.txt'
        }

        if "fastq2" in input_files:
            input_files_trim['fastq2'] = input_files['fastq2']

            files_out["fastq2_trimmed"] = input_files['fastq2'] + '.trimmed'
            files_out["fastq2_report"] = input_files['fastq2'] + '.trimmed.report.txt'

        self.configuration["tg_length"] = "0"

        tg_handle = trimgalore(self.configuration)
        tg_files, tg_meta = tg_handle.run(input_files_trim, metadata, files_out)

        # Align the genome file

        input_files_bowtie = {
            "genome": genome_file,
            "index": input_files["index"],
            "loc": tg_files["fastq1_trimmed"]
        }
        metadata_bowtie = {
            "genome": metadata["genome"],
            "index": metadata["index"],
            "loc":  Metadata(
                "data_atac", "fastq", tg_files["fastq1_trimmed"], None,
                {"assembly": "test"}
                )
        }

        if "fastq2" in input_files:
            input_files_bowtie["fastq2"] = tg_files["fastq2_trimmed"]
            metadata_bowtie["fastq2"] = tg_meta["fastq2_trimmed"]

        output_files_bowtie = {
            "output": input_fastq1.replace(".fastq", "_bt2.bam")
        }

        bowtie2_handle = bowtie2AlignerTool()
        bt_out, bt_meta = bowtie2_handle.run(input_files_bowtie, metadata_bowtie, output_files_bowtie)

        bam_file = output_files_bowtie["output"]
        bam_filtered = bam_file + "filtered"

        input_files_bbb = {
            "input": bam_file
        }

        output_files_bbb = {
            "output": bam_filtered
        }

        metadata_bbb = {
            "input": Metadata(
                "data_atacseq", "fastq", [], None,
                {'assembly': 'test'})
        }

        bbb = biobambam()
        bb_out, bb_meta = bbb.run(input_files_bbb, metadata_bbb, output_files_bbb)

        # Call peaks with macs2

        input_files["bam"] = output_files_bbb["output"]

        """output_files = {
            "narrow_peak": output_narrowpeak,
            "summits": output_summits,
            "broad_peak": output_broadpeak,
            "gapped_peak": output_gappedpeak
        }"""

        metadata_macs = {
            "bam": Metadata(
                bb_meta, [], None,
                {'assembly': 'test'}),
        }

        print("BAM FILE :", bam_file)
        print("Output_files :", output_files["narrow_peak"])
        print("output_narrowPeak :", output_narrowpeak)

        self.configuration["macs_nomodel_param"] = True
        self.configuration["macs_keep-dup_param"] = "all"

        if bedpe:
            self.configuration["macs_format_param"] = "BEDPE"

        macs_handle = macs2(self.configuration)
        macs_out, macs_meta = macs_handle.run(input_files, metadata_macs, output_files)

        output_files.update(bt_out)
        output_files.update(bb_out)
        output_files.update(macs_out)

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

        output_metadata.update(bt_meta)
        output_metadata.update(bb_meta)
        output_metadata.update(macs_meta)

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
