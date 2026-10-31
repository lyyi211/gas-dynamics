"""Create static diagnostics and a synchronized animation from simulator CSVs."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter

parser = argparse.ArgumentParser()
parser.add_argument('directory', nargs='?', default='output')
parser.add_argument('--format', choices=['gif', 'mp4'], default='gif')
parser.add_argument('--skip-animation', action='store_true')
args = parser.parse_args()
root = Path(args.directory)
meta = json.loads((root / 'metadata.json').read_text())
obs = np.atleast_1d(np.genfromtxt(root / 'observables.csv', delimiter=',', names=True))
raw = np.atleast_1d(np.genfromtxt(root / 'trajectory.csv', delimiter=',', names=True))
frames = raw.reshape(-1, meta['n'])
speed = np.hypot(frames['vx'], frames['vy'])
plots = root / 'plots'
plots.mkdir(exist_ok=True)
plt.style.use('dark_background')
plt.rcParams.update({'figure.facecolor': '#101827', 'axes.facecolor': '#101827', 'axes.grid': True, 'grid.alpha': .15})

fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
for key in ['kinetic', 'potential', 'total']:
    axes[0, 0].plot(obs['time'], obs[key], label=key)
axes[0, 0].set(title='Energy conservation', xlabel='Time', ylabel='Energy')
axes[0, 0].legend()
scale = max(abs(obs['total'][0]), obs['kinetic'][0], 1e-12)
axes[0, 1].plot(obs['time'], (obs['total']-obs['total'][0])/scale, color='#ffbb66')
axes[0, 1].set(title='Normalized energy error', xlabel='Time', ylabel='(E − E₀) / max(|E₀|, K₀)')
axes[1, 0].plot(obs['time'], obs['temperature'], color='#58ddcf')
axes[1, 0].set(title='Kinetic temperature', xlabel='Time', ylabel='T (reduced units)')
tail = speed[len(speed)//2:].ravel()
temperature = obs['temperature'][len(obs)//2:].mean()
values = np.linspace(0, tail.max()*1.1, 200)
axes[1, 1].hist(tail, bins=40, density=True, alpha=.7, label='Second half of trajectory')
axes[1, 1].plot(values, values/temperature*np.exp(-values**2/(2*temperature)), label='2D Maxwell reference', color='#ffbb66')
axes[1, 1].set(title='Speed distribution', xlabel='Speed', ylabel='Probability density')
axes[1, 1].legend()
fig.savefig(plots / 'diagnostics.png', dpi=180)
plt.close(fig)

if not args.skip_animation:
    fig, (ax, chart) = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw={'width_ratios': [1.4, 1]}, constrained_layout=True)
    ax.set(xlim=(0, meta['box']), ylim=(0, meta['box']), xlabel='x / σ', ylabel='y / σ', aspect='equal')
    dots = ax.scatter(frames[0]['x'], frames[0]['y'], c=speed[0], s=22, cmap='plasma', vmin=0, vmax=speed.max())
    fig.colorbar(dots, ax=ax, label='Speed')
    chart.plot(obs['time'], obs['temperature'], color='#58ddcf', alpha=.5)
    cursor = chart.axvline(0, color='#ffbb66')
    chart.set(title='Temperature', xlabel='Time', ylabel='T (reduced units)')
    title = ax.set_title('')
    def update(i):
        dots.set_offsets(np.column_stack((frames[i]['x'], frames[i]['y'])))
        dots.set_array(speed[i])
        cursor.set_xdata([obs['time'][i]]*2)
        title.set_text(f"Lennard–Jones gas · N={meta['n']} · t={obs['time'][i]:.2f}")
        return dots, cursor, title
    animation = FuncAnimation(fig, update, frames=len(frames), interval=50)
    writer = PillowWriter(fps=20) if args.format == 'gif' else FFMpegWriter(fps=20)
    animation.save(plots / f'gas_animation.{args.format}', writer=writer, dpi=90)
    plt.close(fig)
print(f'Visualizations saved to {plots}')
