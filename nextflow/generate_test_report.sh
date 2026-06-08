#!/usr/bin/env bash
#SBATCH --job-name=ASAP_test_report
#SBATCH -c 1
#SBATCH --mem=512M
#SBATCH --time=00:10:00
#SBATCH --output=logs/nf-test/report_%j.out
#SBATCH --error=logs/nf-test/report_%j.err
#
# Summarize the results of a parallel ASAP nf-test run (see run_tests_parallel.sh).
#
# Submitted automatically by run_tests_parallel.sh with
# `--dependency=afterany:<job_id_1>:<job_id_2>:...`, so it only starts once
# every per-test job has finished — pass or fail. It reads the manifest written
# by run_tests_parallel.sh (job_id|tag|test_name per line), locates each job's
# `run_*_slurm<job_id>.log` under logs/nf-test/, pulls the PASSED/FAILED status
# and duration out of it, and writes a consolidated markdown report.
#
# Usage: sbatch generate_test_report.sh <manifest_file>
#   <manifest_file> — lines of `job_id|tag|test_name`, one per submitted test job

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <manifest_file>" >&2
    exit 1
fi

MANIFEST_FILE="$1"
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
LOG_DIR="${SCRIPT_DIR}/logs/nf-test"
TIMESTAMP="$(date '+%Y-%m-%d_%H-%M-%S')"
REPORT_FILE="${LOG_DIR}/report_${TIMESTAMP}.md"

if [[ ! -f "${MANIFEST_FILE}" ]]; then
    echo "Manifest file not found: ${MANIFEST_FILE}" >&2
    exit 1
fi

# Strip ANSI escape sequences (color codes, charset-select sequences, etc.) so
# the status lines and assertion messages can be matched/read as plain text.
ESC=$'\033'
strip_ansi() {
    sed -E "s/${ESC}\[[0-9;]*[A-Za-z]//g; s/${ESC}[()][A-Za-z0-9]?//g"
}

PASS_COUNT=0
FAIL_COUNT=0
MISSING_COUNT=0

SUMMARY_ROWS=""
FAILURE_DETAILS=""

while IFS='|' read -r job_id tag name; do
    [[ -z "${job_id}" ]] && continue

    log_file="$(ls -t "${LOG_DIR}"/run_*_slurm"${job_id}".log 2>/dev/null | head -n1 || true)"

    if [[ -z "${log_file}" ]]; then
        SUMMARY_ROWS+="| ${name} | \`${tag}\` | ${job_id} | ⚠️ NO LOG | — | — |"$'\n'
        MISSING_COUNT=$((MISSING_COUNT + 1))
        continue
    fi

    plain_log="$(strip_ansi < "${log_file}")"
    log_name="$(basename "${log_file}")"

    # Each job runs exactly one test, so there is exactly one status line like:
    #   PASSED (123.45s)   or   FAILED (123.45s)
    status_line="$(grep -m1 -E '^\s*(PASSED|FAILED) \([0-9.]+s\)' <<< "${plain_log}" || true)"
    status="$(grep -oE 'PASSED|FAILED' <<< "${status_line}" | head -n1 || true)"
    duration="$(grep -oE '\([0-9.]+s\)' <<< "${status_line}" | head -n1 | tr -d '()' || true)"

    if [[ "${status}" == "PASSED" ]]; then
        PASS_COUNT=$((PASS_COUNT + 1))
        SUMMARY_ROWS+="| ${name} | \`${tag}\` | ${job_id} | ✅ PASSED | ${duration:-—} | ${log_name} |"$'\n'
    elif [[ "${status}" == "FAILED" ]]; then
        FAIL_COUNT=$((FAIL_COUNT + 1))
        SUMMARY_ROWS+="| ${name} | \`${tag}\` | ${job_id} | ❌ FAILED | ${duration:-—} | ${log_name} |"$'\n'

        # Pull the assertion-failure block (between "Assertion failed:" and the
        # "Nextflow stdout:" dump that follows it in nf-test's verbose output).
        assertion="$(sed -n '/Assertion failed:/,/Nextflow stdout:/p' <<< "${plain_log}" \
            | sed '$d' | sed '/^\s*$/d')"

        FAILURE_DETAILS+="### ${name} (\`${tag}\`, job ${job_id})"$'\n\n'
        FAILURE_DETAILS+='```'$'\n'
        if [[ -n "${assertion}" ]]; then
            FAILURE_DETAILS+="${assertion}"$'\n'
        else
            FAILURE_DETAILS+="(could not extract assertion details — see ${log_name})"$'\n'
        fi
        FAILURE_DETAILS+='```'$'\n\n'
    else
        SUMMARY_ROWS+="| ${name} | \`${tag}\` | ${job_id} | ⚠️ UNKNOWN | — | ${log_name} |"$'\n'
        MISSING_COUNT=$((MISSING_COUNT + 1))
    fi
done < "${MANIFEST_FILE}"

TOTAL=$((PASS_COUNT + FAIL_COUNT + MISSING_COUNT))

{
    echo "# ASAP nf-test Parallel Run Report"
    echo ""
    echo "Generated: ${TIMESTAMP}"
    echo ""
    echo "**Result: ${PASS_COUNT}/${TOTAL} passed**, ${FAIL_COUNT} failed, ${MISSING_COUNT} missing/unknown"
    echo ""
    echo "| Test | Tag | Job ID | Status | Duration | Log |"
    echo "|---|---|---|---|---|---|"
    printf '%s' "${SUMMARY_ROWS}"

    if [[ -n "${FAILURE_DETAILS}" ]]; then
        echo ""
        echo "## Failure Details"
        echo ""
        printf '%s' "${FAILURE_DETAILS}"
    fi
} > "${REPORT_FILE}"

echo "Report written to: ${REPORT_FILE}"
cat "${REPORT_FILE}"

# Exit non-zero if anything failed or couldn't be found, so the SLURM job
# status itself reflects whether the overall suite passed.
[[ ${FAIL_COUNT} -eq 0 && ${MISSING_COUNT} -eq 0 ]]
