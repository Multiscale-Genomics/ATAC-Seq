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

import os.path
import pytest

from basic_modules.metadata import Metadata
from mg_process_fastq.tool.bowtie_aligner import bowtie2AlignerTool
from mg_process_fastq.tool.bowtie_indexer import bowtieIndexerTool
from process_atac_seq import process_atac_seq


@pytest.mark.atacseq
@pytest.mark.pipeline
def test_atac_seq():
    """
    Test case to ensure that the atac seq pipeline code works.
    """

    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    files = {
        'genome': resource_path + 'atac.Human.hg19.fasta',
        'fastq1': resource_path + 'ERR1659027_1.fastq',
        'fastq2': resource_path + 'ERR1659027_2.fastq'
    }

    metadata = {}

    files_out = {
            "narrow_peak": resource_path + "atacseq.Human.ERR1659027_peaks.narrowPeak",
            "summits": resource_path + "atacseq.Human.ERR1659027_peaks.summits.bed",
            "broad_peak": resource_path + "atacseq.Human.ERR1659027_peaks.broadPeak",
            "gapped_peak": resource_path + "atacseq.Human.ERR1659027_peaks.gappedPeak"
    }

    atac_handle = process_atac_seq()
    atac_handle.run(files, metadata, files_out)
