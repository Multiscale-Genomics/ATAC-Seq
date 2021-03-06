#!/bin/bash

# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

rc=0
pv=$(python -c 'import platform; print(platform.python_version())')

#if [[ $pv == "2.7.12" ]]; then
echo "File directory : "
pwd
echo "before tests: "
ls -l atac_seq/tests/data/
bowtie2-build atac.Human.hg19.fasta atac.Human.hg19
pytest atac_seq/tests/test_pipeline_atac_seq.py
tc=$?
rc=$(($rc + $tc))
echo "after tests: "
ls -l atac_seq/tests/data/
bash tidy_data.sh
echo "Test running"
#$fi


if [[ $rc != 0 ]]; then exit $rc; fi
