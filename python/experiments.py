"""Run repeated timing experiments and plot serial, OpenMP and grid scaling."""
import csv
import json
from pathlib import Path
import subprocess
import tempfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parents[1]
output = root / 'output' / 'benchmarks'
output.mkdir(parents=True, exist_ok=True)
rows = []
with tempfile.TemporaryDirectory() as tmp:
    for n in [128, 256, 512, 1024]:
        for method, threads in [('sequential', 1), ('grid', 1), ('openmp', 1), ('openmp', 2), ('openmp', 4)]:
            timings = []
            for repeat in range(3):
                subprocess.run([str(root/'gas'), '--n', str(n), '--method', method, '--threads', str(threads), '--steps', '100', '--stride', '100', '--out', tmp], check=True, stdout=subprocess.DEVNULL)
                timings.append(json.loads((Path(tmp)/'metadata.json').read_text())['compute_seconds'])
            rows.append([n, method, threads, float(np.median(timings))])
with (output/'timings.csv').open('w', newline='') as f:
    writer = csv.writer(f, lineterminator="\n"); writer.writerow(['n', 'method', 'threads', 'median_seconds']); writer.writerows(rows)
fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
for method, threads in [('sequential', 1), ('grid', 1), ('openmp', 4)]:
    data = [r for r in rows if r[1:3] == [method, threads]]
    axes[0].loglog([r[0] for r in data], [r[3] for r in data], 'o-', label=f'{method} ({threads} threads)')
for n in [128, 256, 512, 1024]:
    data = [r for r in rows if r[0] == n and r[1] == 'openmp']
    axes[1].plot([r[2] for r in data], [data[0][3]/r[3] for r in data], 'o-', label=f'N={n}')
axes[1].plot([1, 4], [1, 4], '--', color='gray', label='Ideal')
axes[0].set(xlabel='Particle count', ylabel='Seconds / 100 steps', title='Median of 3 runs; CSV I/O excluded')
axes[1].set(xlabel='Threads', ylabel='Speedup relative to OpenMP / 1 thread', title='Parallel scaling')
for ax in axes: ax.legend(); ax.grid(alpha=.2)
fig.savefig(output/'scaling.png', dpi=180)
print(f'Benchmarks saved to {output}')
