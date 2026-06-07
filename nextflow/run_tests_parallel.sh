#!/usr/bin/env bash
#
# Submit each ASAP end-to-end nf-test as its own SLURM job, running in parallel.
#
# run_tests.sh runs the whole suite serially (60+ minutes). This script instead
# loops over a list of tags — one tag per test, each matching EXACTLY ONE test
# in tests/ASAP_EtE.nf.test — and submits one `sbatch run_tests.sh --tag <tag>`
# job per test. SLURM then runs them concurrently (resources permitting), so
# the suite finishes in roughly the time of its single longest test.
#
# IMPORTANT: each tag below must resolve to exactly one test. nf-test isolates
# tests by a content-hash work directory and each test owns its own `outdir`
# under tests/test_output/<name>, but two CONCURRENT runs that both select the
# same test would collide on that hash directory and on real output files.
# If you add a new test to ASAP_EtE.nf.test, give it its own unique tag and add
# that tag here — verify uniqueness with:
#   grep -c 'tag "<tag>"' tests/ASAP_EtE.nf.test    # must print 1
#
# Usage (from the nextflow/ directory):
#   ./run_tests_parallel.sh
#
# Each line below maps a unique tag -> the test it selects:
TAGS=(
    help            # Help Documentation
    rsv             # RSV - Multi-GenBank Reference, BWA, Paired-End, Combine Output
    bowtie2         # SC2 - JSON Direct Input, Bowtie2, Paired-End, Combine Output
    sc2_se_bwa      # SC2 - Single-End Illumina, BWA
    ont             # SC2 - ONT Long Reads, Minimap2
    ivar            # SC2 - Primer Masking + iVAR Variant Calling
    primer_only     # SC2 - Primer-Only BAM Filtering
    smor            # TB - Full Feature Amplicon (primer masking, identity filter, SMOR, asaptools)
    excel_input     # TB - Excel Reference Input
    whole_genome    # TB - FASTA Reference Input, Whole Genome Mode
    tb_json_snp_aa  # TB - JSON Input, asaptools with GenBank SNP-to-AA Conversion
)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

mkdir -p logs/nf-test

echo "Submitting ${#TAGS[@]} parallel nf-test jobs (one test per job)..."
echo ""

for tag in "${TAGS[@]}"; do
    job_name="ASAP_test_${tag}"
    sbatch --job-name="${job_name}" run_tests.sh --tag "${tag}"
done

echo ""
echo "Track progress with: squeue -u \$USER --name=$(IFS=,; echo "${TAGS[*]/#/ASAP_test_}")"
echo "Logs will appear under: ${SCRIPT_DIR}/logs/nf-test/"
