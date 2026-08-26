#!/bin/bash
set -euo pipefail

# ---------------------------------------------------------------------------
# Pipeline base dir = the directory where this script lives (follows symlinks).
# This way the script also works from a feature-branch checkout.
# ---------------------------------------------------------------------------
script_path=$(readlink -f "${BASH_SOURCE[0]}")
workflow_path=$(dirname "${script_path}")
pipeline_name=$(basename "${workflow_path}")

usage() {
cat <<EOF
Usage: $(basename "$0") -i <input> -o <output> -e <email> [options] [-- <nextflow args>]

Required:
  -i, --input PATH       Input dir/file
  -o, --output DIR       Output/analysis dir (basename becomes analysis_id)
  -e, --email ADDR       Mail for SLURM failure and --email

Optional (SLURM):
  -j, --job-name NAME    Job name            (default: Nextflow_${pipeline_name})
  -t, --time HH:MM:SS    Walltime            (default: 48:00:00)
  -m, --mem SIZE         Memory              (default: 10G)
      --tmpspace SIZE    tmpspace            (default: 10G)
      --account NAME     SLURM account       (default: diaggen)
  -n, --dry-run          Print the job script instead of submitting it
  -h, --help             This text

All remaining arguments are passed straight through to 'nextflow run'.
The defaults -resume, -ansi-log false, -profile singularity and
-c <pipeline>/nextflow.config are only set if you don't supply them
yourself. Use '--' if a nextflow argument clashes with one of the flags
above.

Examples:
  $(basename "$0") -i /data/run1 -o /data/analysis/run1 -e j@umcutrecht.nl
  $(basename "$0") -i in -o out -e j@x.nl -profile singularity,test --max_cpus 8
  $(basename "$0") -i in -o out -e j@x.nl -t 12:00:00 -n -- -stub-run
EOF
}

# --- defaults --------------------------------------------------------------
input=""
output=""
email=""
job_name=""
sbatch_time="48:00:00"
sbatch_mem="10G"
sbatch_tmpspace="10G"
sbatch_account="diaggen"
dry_run=false
extra_args=()

# --- argument parsing -------------------------------------------------------
while [ $# -gt 0 ]; do
    case "$1" in
        -i|--input)     input="$2";           shift 2 ;;
        -o|--output)    output="$2";          shift 2 ;;
        -e|--email)     email="$2";           shift 2 ;;
        -j|--job-name)  job_name="$2";        shift 2 ;;
        -t|--time)      sbatch_time="$2";     shift 2 ;;
        -m|--mem)       sbatch_mem="$2";      shift 2 ;;
        --tmpspace)     sbatch_tmpspace="$2"; shift 2 ;;
        --account)      sbatch_account="$2";  shift 2 ;;
        -n|--dry-run)   dry_run=true;         shift ;;
        -h|--help)      usage; exit 0 ;;
        --)             shift; extra_args+=( "$@" ); break ;;
        *)              extra_args+=( "$1" ); shift ;;
    esac
done

missing=()
[ -n "${input}" ]  || missing+=( "--input" )
[ -n "${output}" ] || missing+=( "--output" )
[ -n "${email}" ]  || missing+=( "--email" )
if [ ${#missing[@]} -gt 0 ]; then
    echo "ERROR: missing arguments: ${missing[*]}" >&2
    echo >&2
    usage >&2
    exit 1
fi

[ -n "${job_name}" ] || job_name="Nextflow_${pipeline_name}"

for f in "${workflow_path}/main.nf" "${workflow_path}/tools/nextflow"; do
    [ -e "${f}" ] || { echo "ERROR: not found: ${f}" >&2; exit 1; }
done

# --- paths -------------------------------------------------------------------
input=$(realpath "${input}")
mkdir -p "${output}"
output=$(realpath "${output}")
analysis_id=$(basename "${output}")

cd "${output}"
mkdir -p log

if [ -f 'workflow.running' ] || [ -f 'workflow.done' ] || [ -f 'workflow.failed' ]; then
    echo "Workflow job not submitted, please check ${output} for 'workflow.status' files."
    exit 0
fi

# --- build nextflow command -------------------------------------------------
# helper: is this (nextflow) flag already present in the user-supplied args?
has_flag() {
    local needle="$1" arg
    for arg in ${extra_args[@]+"${extra_args[@]}"}; do
        if [ "${arg}" = "${needle}" ] || [ "${arg#${needle}=}" != "${arg}" ]; then
            return 0
        fi
    done
    return 1
}

nf_args=( run "${workflow_path}/main.nf" )
has_flag -c        || has_flag -config || nf_args+=( -c "${workflow_path}/nextflow.config" )
nf_args+=( --input "${input}" --outdir "${output}" --analysis_id "${analysis_id}" --email "${email}" )
has_flag -resume   || nf_args+=( -resume )
has_flag -ansi-log || nf_args+=( -ansi-log false )
nf_args+=( ${extra_args[@]+"${extra_args[@]}"} )

# quote safely so it doesn't fall apart inside the heredoc
nf_cmd=$(printf '%q ' "${workflow_path}/tools/nextflow" "${nf_args[@]}")

# --- rotate nextflow_trace.txt ------------------------------------------------
trace_file="${output}/log/nextflow_trace.txt"
if [ -f "${trace_file}" ]; then
    max_suffix=0
    for f in "${output}"/log/nextflow_trace_*.txt; do
        [ -e "${f}" ] || continue
        n=$(basename "${f}" .txt); n="${n##*_}"
        case "${n}" in
            ''|*[!0-9]*) continue ;;
        esac
        [ "${n}" -gt "${max_suffix}" ] && max_suffix="${n}"
    done
    mv "${trace_file}" "${output}/log/nextflow_trace_$((max_suffix + 1)).txt"
fi

# --- job script ------------------------------------------------------------
job_script=$(cat <<EOT
#!/bin/bash
#SBATCH --time=${sbatch_time}
#SBATCH --nodes=1
#SBATCH --mem ${sbatch_mem}
#SBATCH --gres=tmpspace:${sbatch_tmpspace}
#SBATCH --job-name ${job_name}
#SBATCH -o log/slurm_${job_name}.%j.out
#SBATCH -e log/slurm_${job_name}.%j.err
#SBATCH --mail-user ${email}
#SBATCH --mail-type FAIL
#SBATCH --export=NONE
#SBATCH --account=${sbatch_account}

export NXF_JAVA_HOME='${workflow_path}/tools/java/jdk'

${nf_cmd}

if [ \$? -eq 0 ]; then
    echo "Nextflow done."

    echo "Zip work directory"
    find work -type f | egrep "\.(command|exitcode)" | zip -@ -q work.zip

    echo "Remove work directory"
    rm -r work

    echo "Creating md5sum"
    find -type f -not -iname 'md5sum.txt' -exec md5sum {} \; > md5sum.txt

    echo "${pipeline_name} workflow completed successfully."
    rm workflow.running
    touch workflow.done

    echo "Change permissions"
    chmod 775 -R ${output}

    exit 0
else
    echo "Nextflow failed"
    rm workflow.running
    touch workflow.failed

    echo "Change permissions"
    chmod 775 -R ${output}

    exit 1
fi
EOT
)

if [ "${dry_run}" = true ]; then
    echo "--- dry-run: job script (not submitted) ---"
    echo "${job_script}"
    exit 0
fi

touch workflow.running
sbatch <<< "${job_script}"
