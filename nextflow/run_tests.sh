#!/usr/bin/env bash
#SBATCH --job-name=ASAP_EtE_tests
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --time=4-00:00:00
#SBATCH --output=logs/nf-test/ASAP_End_To_End_Testing_%j.out
#SBATCH --error=logs/nf-test/ASAP_End_To_End_Testing_%j.err
#
# Run ASAP end-to-end nf-tests.
#
# Direct usage (interactive / srun session):
#   ./run_tests.sh [--tag <tag>] [--keep-going] [--update-snapshot]
#
# SLURM submission — run from the nextflow/ directory:
#   cd /path/to/ASAP/nextflow
#   mkdir -p logs/nf-test   # must exist before sbatch (SLURM opens log files early)
#   sbatch run_tests.sh
#
# Tag-scoped runs:
#   sbatch --job-name=ASAP_help       run_tests.sh --tag help
#   sbatch --job-name=ASAP_rsv        run_tests.sh --tag rsv
#   sbatch --job-name=ASAP_tb         run_tests.sh --tag tb
#   sbatch --job-name=ASAP_sc2        run_tests.sh --tag sc2
#   sbatch --job-name=ASAP_ont        run_tests.sh --tag ont
#   sbatch --job-name=ASAP_ivar       run_tests.sh --tag ivar
#   sbatch --job-name=ASAP_asaptools  run_tests.sh --tag asaptools
#   sbatch --job-name=ASAP_primers    run_tests.sh --tag primer_masking
#
# NOTE: the tags above overlap (e.g. `tb` and `sc2` each match several tests,
# `ont`/`ivar`/`asaptools`/`primer_masking` are subsets of `tb`/`sc2`). Submitting
# overlapping tag-scoped jobs concurrently would run the same test twice in
# parallel and collide on its shared .nf-test work directory and outdir — only
# run one of these at a time, or run them sequentially.
#
# Run the WHOLE suite in parallel (one SLURM job per test, ~11 jobs at once):
#   ./run_tests_parallel.sh
# This loops over a curated list of tags that each match EXACTLY ONE test (see
# run_tests_parallel.sh for the tag → test mapping and the uniqueness rule to
# follow when adding new tests) and submits one `sbatch run_tests.sh --tag <tag>`
# job per test, so the suite finishes in roughly the time of its longest test.
#
# Available tags (defined in tests/ASAP_EtE.nf.test):
#   help            – CLI help output
#   rsv             – RSV amplicon tests
#   tb              – TB amplicon tests
#   sc2             – SARS-CoV-2 tests
#   bwa             – BWA aligner path
#   bowtie2         – Bowtie2 aligner path
#   minimap2        – Minimap2 / ONT path
#   multi_gb        – Multiple GenBank reference files
#   json_input      – Pre-built JSON passed directly (no conversion)
#   excel_input     – Excel → JSON conversion
#   fasta_input     – FASTA → JSON conversion
#   paired_end      – Paired-end Illumina reads
#   single_end      – Single-end Illumina reads
#   ont             – ONT long reads + fastplong
#   primer_masking  – maskPrimers.py step
#   identity_filter – identityFilter.py step
#   smor            – SMOR / SMOR-correction step
#   ivar            – iVAR trim + variant + consensus branch
#   asaptools       – R post-processing (SNP table, cov table, etc.)
#   combine_output  – OUTPUT_COMBINER + FORMAT_OUTPUT (XML → HTML)
#   whole_genome    – --whole_genome mode

set -euo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
TEST_FILE="tests/ASAP_EtE.nf.test"
LOG_DIR="${SCRIPT_DIR}/logs/nf-test"
TIMESTAMP="$(date '+%Y-%m-%d_%H-%M-%S')"
SLURM_LABEL="${SLURM_JOB_ID:+_slurm${SLURM_JOB_ID}}"
LOG_FILE="${LOG_DIR}/run_${TIMESTAMP}${SLURM_LABEL}.log"

mkdir -p "${LOG_DIR}"

# ── Parse arguments ────────────────────────────────────────────────────────────
TAG=""
PROFILE="slurm"
KEEP_GOING=false
UPDATE_SNAPSHOT=false
EXTRA_FLAGS=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tag)             TAG="$2"; shift 2 ;;
        --profile)         PROFILE="$2"; shift 2 ;;
        --keep-going)      KEEP_GOING=true; shift ;;
        --update-snapshot) UPDATE_SNAPSHOT=true; shift ;;
        *)                 EXTRA_FLAGS="${EXTRA_FLAGS} $1"; shift ;;
    esac
done

# ── Build nf-test command ──────────────────────────────────────────────────────
NF_TEST_CMD="nf-test test ${TEST_FILE} --profile ${PROFILE}"

[[ -n "${TAG}" ]]              && NF_TEST_CMD="${NF_TEST_CMD} --tag ${TAG}"
[[ "${KEEP_GOING}" == false ]] && NF_TEST_CMD="${NF_TEST_CMD} --stop-on-first-failure"
[[ "${UPDATE_SNAPSHOT}" == true ]] && NF_TEST_CMD="${NF_TEST_CMD} --update-snapshot"

NF_TEST_CMD="${NF_TEST_CMD} --verbose${EXTRA_FLAGS}"

# ── Run ───────────────────────────────────────────────────────────────────────
cd "${SCRIPT_DIR}"

echo "========================================"
echo " ASAP nf-test run"
echo " Started:  ${TIMESTAMP}"
[[ -n "${SLURM_JOB_ID:-}"   ]] && echo " SLURM ID: ${SLURM_JOB_ID}"
[[ -n "${SLURM_JOB_NAME:-}" ]] && echo " Job name: ${SLURM_JOB_NAME}"
[[ -n "${TAG}"               ]] && echo " Tag:      ${TAG}"
echo " Profile:  ${PROFILE}"
echo " Log:      ${LOG_FILE}"
echo " Command:  ${NF_TEST_CMD}"
echo "========================================"
echo ""

set +e
${NF_TEST_CMD} 2>&1 | tee "${LOG_FILE}"
EXIT_CODE=${PIPESTATUS[0]}
set -e

echo ""
echo "========================================"
if [[ ${EXIT_CODE} -eq 0 ]]; then
    echo " RESULT: ALL TESTS PASSED"
else
    echo " RESULT: FAILED (exit ${EXIT_CODE})"
    echo " Full log: ${LOG_FILE}"
fi
echo " Finished: $(date '+%Y-%m-%d_%H-%M-%S')"
echo "========================================"

exit ${EXIT_CODE}
