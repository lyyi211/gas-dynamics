"""Structural/transport diagnostics and equal-duration time-step convergence."""
import json
import subprocess
import tempfile
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parents[1]
out = root/'output'
meta = json.loads((out/'metadata.json').read_text())
obs = np.genfromtxt(out/'observables.csv', names=True, delimiter=',')
if len(obs) < 8:
    raise ValueError('Research analysis requires at least eight saved frames.')
if not np.allclose(np.diff(obs['time']), np.diff(obs['time'])[0]):
    raise ValueError('Research analysis requires uniformly spaced saved frames (steps divisible by stride).')
f = np.genfromtxt(out/'trajectory.csv', names=True, delimiter=',').reshape(-1, meta['n'])
L, n = meta['box'], meta['n']
half = len(f)//2
pos = np.stack((f['x'], f['y']), axis=-1)
unwrapped = np.stack((f['xu'], f['yu']), axis=-1)[half:]
vel = np.stack((f['vx'], f['vy']), axis=-1)[half:]
edges = np.linspace(0, L/2, 101)
counts = np.zeros(100)
ii, jj = np.triu_indices(n, 1)
for frame in pos[half:]:
    delta = frame[ii]-frame[jj]
    delta -= L*np.round(delta/L)
    counts += np.histogram(np.linalg.norm(delta, axis=1), edges)[0]
area = np.pi*np.diff(edges**2)
g = counts / ((len(f)-half)*n*(n-1)/2*area/L**2)
r = (edges[1:]+edges[:-1])/2
np.savetxt(out/'rdf.csv', np.column_stack((r, g)), delimiter=',', header='r,g', comments='')
# Time-origin averaging uses only the second half, up to half its duration.
lags = np.arange(1, len(unwrapped)//2+1)
msd = np.array([np.mean(np.sum((unwrapped[k:]-unwrapped[:-k])**2, axis=-1)) for k in lags])
lagtime = lags*(obs['time'][1]-obs['time'][0])
vacf = np.array([np.mean(np.sum(vel[k:]*vel[:-k], axis=-1)) for k in lags])
np.savetxt(out/'transport.csv', np.column_stack((lagtime, msd, vacf)), delimiter=',', header='lag,msd,vacf', comments='')
fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
axes[0,0].plot(r,g); axes[0,0].axhline(1, ls='--', color='gray');axes[0,0].set(xlabel='r / sigma', ylabel='g(r)', title='Pair structure (second half)')
axes[0,1].plot(obs['time'],obs['pressure'],label='Total mechanical pressure')
axes[0,1].plot(obs['time'],obs['kinetic']/L**2,label='Kinetic contribution',alpha=.7)
axes[0,1].set(xlabel='Time',ylabel='Pressure',title='Mechanical pressure');axes[0,1].legend()
axes[1,0].loglog(lagtime,msd,label='Multi-origin MSD')
axes[1,0].loglog(lagtime, np.mean(np.sum(vel**2,axis=-1))*lagtime**2,'--',label='Ballistic reference')
axes[1,0].set(xlabel='Lag time',ylabel='MSD',title='Transport: no diffusion fit assumed');axes[1,0].legend()
axes[1,1].plot(lagtime,vacf);axes[1,1].axhline(0,color='gray',ls='--')
axes[1,1].set(xlabel='Lag time',ylabel='<v(t) · v(t+lag)>',title='Velocity autocorrelation')
for ax in axes.flat:ax.grid(alpha=.2)
fig.savefig(out/'plots/structure_transport.png',dpi=180);plt.close(fig)
rows=[]
with tempfile.TemporaryDirectory() as tmp:
    for dt in [.008,.004,.002,.001]:
        subprocess.run([str(root/'gas'),'--n','64','--dt',str(dt),'--steps',str(round(2/dt)),'--stride','1','--out',tmp],check=True,stdout=subprocess.DEVNULL)
        d=np.genfromtxt(Path(tmp)/'observables.csv',names=True,delimiter=',')
        error=np.max(abs(d['total']-d['total'][0]))/max(abs(d['total'][0]),d['kinetic'][0])
        rows.append([dt,error])
rows=np.array(rows)
order=float(np.polyfit(np.log(rows[:,0]),np.log(rows[:,1]),1)[0])
np.savetxt(out/'convergence.csv',rows,delimiter=',',header='dt,max_normalized_energy_error',comments='')
fig,ax=plt.subplots(figsize=(6,4),constrained_layout=True)
ax.loglog(rows[:,0],rows[:,1],'o-',label=f'Measured slope {order:.2f}')
ax.loglog(rows[:,0],rows[-1,1]*(rows[:,0]/rows[-1,0])**2,'--',label='Second-order reference')
ax.set(xlabel='Time step',ylabel='Maximum normalized energy error',title='Equal physical duration t=2; N=64')
ax.legend();ax.grid(alpha=.2);fig.savefig(out/'plots/convergence.png',dpi=180)
summary={'convergence_slope':order,'demo_max_normalized_energy_error':float(np.max(abs(obs['total']-obs['total'][0]))/max(abs(obs['total'][0]),obs['kinetic'][0])),'second_half_mean_temperature':float(obs['temperature'][half:].mean()),'second_half_mean_pressure':float(obs['pressure'][half:].mean()),'rdf_peak_r':float(r[np.argmax(g)]),'rdf_peak_height':float(max(g))}
(out/'research_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
