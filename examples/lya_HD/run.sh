#!/bin/bash
# Run all H+D Lyα reference cases. Uses 90% of available CPU threads.
# After completion, open inspect_HD.ipynb to inspect the spectra.

NP=$(python3 -c "import os; print(max(1, int(os.cpu_count()*0.9)))")
echo "Using NP=$NP MPI ranks."

set -e
for inp in sphere_H_only.in sphere_HD_zero.in sphere_HD_cosmic.in sphere_HD_high.in sphere_HD_dijkstra2006.in; do
  echo "===== Running $inp ====="
  mpirun -np "$NP" ../../LaRT.x "$inp"
done

echo "Done. Open inspect_HD.ipynb in Jupyter to view results."
