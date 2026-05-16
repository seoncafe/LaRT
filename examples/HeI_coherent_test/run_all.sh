#!/bin/bash
# Run all 16 test cases for the HeI 10833 coherent vs. incoherent comparison.
#
# Usage:
#   bash run_all.sh           # default: NP=64 (90% of 72 cores)
#   NP=32 bash run_all.sh     # override
#
# Existing output files (*.fits.gz) are deleted before each run so the
# par%out_merge=.true. default doesn't accumulate across runs.

set -e
cd "$(dirname "$0")"

EXEC=../../LaRT.x
NP=${NP:-$(python3 -c "import os; print(max(1, int(os.cpu_count()*0.9)))")}
LOG_DIR=logs
mkdir -p "$LOG_DIR"

echo "Running with -np $NP"
echo "Executable: $EXEC"

cases=()
for src in pt un; do
  for tau in 0.1 1 10 100 1000; do
    for coh in inc coh; do
      cases+=("${src}_tau${tau}_${coh}")
    done
  done
done

start_all=$(date +%s)
for c in "${cases[@]}"; do
  in_file="${c}.in"
  log_file="$LOG_DIR/${c}.log"

  if [[ ! -f "$in_file" ]]; then
    echo "Missing $in_file -- run generate_inputs.py first"
    exit 1
  fi

  # Clean previous output so out_merge does not accumulate.
  rm -f "${c}.h5" "${c}_obs.h5" "${c}_stokes.h5" \
        "${c}.fits.gz" "${c}_obs.fits.gz" "${c}_stokes.fits.gz"

  echo "----- $(date +%H:%M:%S) running $c -----"
  t0=$(date +%s)
  mpirun -np "$NP" "$EXEC" "$in_file" > "$log_file" 2>&1
  t1=$(date +%s)
  printf "    done in %d s (log: %s)\n" $((t1 - t0)) "$log_file"
done
end_all=$(date +%s)
echo
echo "All 16 cases finished in $((end_all - start_all)) s"
