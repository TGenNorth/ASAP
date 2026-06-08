#!/usr/bin/env bash
#
# Submit each ASAP end-to-end nf-test as its own SLURM job, running in parallel,
# then submit a report job that waits on all of them and summarizes PASS/FAIL.
#
# run_tests.sh runs the whole suite serially (60+ minutes). This script instead
# loops over a list of tags — one tag per test, each matching EXACTLY ONE test
# in tests/ASAP_EtE.nf.test — and submits one `sbatch run_tests.sh --tag <tag>`
# job per test. SLURM then runs them concurrently (resources permitting), so
# the suite finishes in roughly the time of its single longest test. Once every
# test job has finished (pass OR fail), a final `generate_test_report.sh` job
# runs automatically (via `--dependency=afterany:<job_ids>`) and writes a
# consolidated markdown summary under logs/nf-test/.
#
# IMPORTANT: each tag below must resolve to exactly one test. nf-test isolates
# tests by a content-hash work directory and each test owns its own `outdir`
# under tests/test_output/<name>, but two CONCURRENT runs that both select the
# same test would collide on that hash directory and on real output files.
# If you add a new test to ASAP_EtE.nf.test, give it its own unique tag and add
# both the tag AND its test name (exactly as it appears in the `test("...")`
# line — the report job matches on this) to the arrays below — verify
# uniqueness with:
#   grep -c 'tag "<tag>"' tests/ASAP_EtE.nf.test    # must print 1
#
# Usage (from the nextflow/ directory):
#   ./run_tests_parallel.sh

# Parallel arrays: TAGS[i] is the unique selector tag for NAMES[i]'s test.
TAGS=(
    help
    rsv
    bowtie2
    sc2_se_bwa
    ont
    ivar
    primer_only
    smor
    excel_input
    whole_genome
    tb_json_snp_aa
)
NAMES=(
    "Help Documentation"
    "RSV - Multi-GenBank Reference, BWA, Paired-End, Combine Output"
    "SC2 - JSON Direct Input, Bowtie2, Paired-End, Combine Output"
    "SC2 - Single-End Illumina, BWA"
    "SC2 - ONT Long Reads, Minimap2"
    "SC2 - Primer Masking + iVAR Variant Calling"
    "SC2 - Primer-Only BAM Filtering"
    "TB - Full Feature Amplicon (primer masking, identity filter, SMOR, asaptools)"
    "TB - Excel Reference Input"
    "TB - FASTA Reference Input, Whole Genome Mode"
    "TB - JSON Input, asaptools with GenBank SNP-to-AA Conversion"
)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

LOG_DIR="${SCRIPT_DIR}/logs/nf-test"
mkdir -p "${LOG_DIR}"

TIMESTAMP="$(date '+%Y-%m-%d_%H-%M-%S')"
MANIFEST_FILE="${LOG_DIR}/parallel_${TIMESTAMP}_manifest.txt"
: > "${MANIFEST_FILE}"

echo "Submitting ${#TAGS[@]} parallel nf-test jobs (one test per job)..."
echo ""

JOB_IDS=()
for i in "${!TAGS[@]}"; do
    tag="${TAGS[$i]}"
    name="${NAMES[$i]}"
    job_name="ASAP_test_${tag}"

    job_id="$(sbatch --parsable --job-name="${job_name}" run_tests.sh --tag "${tag}")"
    JOB_IDS+=("${job_id}")

    echo "  [${job_id}] ${job_name}  ->  ${name}"
    printf '%s|%s|%s\n' "${job_id}" "${tag}" "${name}" >> "${MANIFEST_FILE}"
done

# Run the report once every test job has finished, whether it passed or failed
# (afterany), so a few failures don't prevent the summary from being generated.
DEPENDENCY="afterany:$(IFS=:; echo "${JOB_IDS[*]}")"
report_job_id="$(sbatch --parsable --job-name="ASAP_test_report" --dependency="${DEPENDENCY}" generate_test_report.sh "${MANIFEST_FILE}")"

echo ""
echo "Submitted ${#JOB_IDS[@]} test jobs and report job [${report_job_id}] (depends on: ${DEPENDENCY})"
echo "Manifest: ${MANIFEST_FILE}"
echo ""
echo "Track progress with: squeue -u \$USER --name=$(IFS=,; echo "${TAGS[*]/#/ASAP_test_}"),ASAP_test_report"
echo "Logs will appear under: ${LOG_DIR}/"
