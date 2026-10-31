"""Physics and cross-implementation regression checks."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np

root = Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory() as tmp:
    results = {}
    for method in ['sequential', 'openmp', 'grid']:
        target = Path(tmp)/method
        subprocess.run([str(root/'gas'), '--n', '64', '--steps', '1000', '--stride', '10', '--method', method, '--out', str(target)], check=True)
        data = np.genfromtxt(target/'observables.csv', names=True, delimiter=',')
        trajectory = np.genfromtxt(target/'trajectory.csv', names=True, delimiter=',')
        results[method] = trajectory
        drift = np.max(np.abs(data['total']-data['total'][0]))/max(abs(data['total'][0]), data['kinetic'][0])
        momentum = np.max(np.hypot(data['px'], data['py']))
        assert drift < 1e-3, (method, drift)
        assert momentum < 1e-9, (method, momentum)
        print(f'{method}: maximum energy error={drift:.3g}, momentum={momentum:.3g}')
    for method in ['openmp', 'grid']:
        for field in ['x', 'y', 'vx', 'vy']:
            np.testing.assert_allclose(results[method][field], results['sequential'][field], atol=1e-8, rtol=1e-8)
print('PASS: energy, momentum, and all three implementations agree.')
# Recompute energy and pressure independently from saved coordinates.
with tempfile.TemporaryDirectory() as tmp:
    subprocess.run([str(root/'gas'), '--n', '64', '--steps', '1000', '--stride', '10', '--out', tmp], check=True)
    import json
    meta = json.loads((Path(tmp)/'metadata.json').read_text())
    frames = np.genfromtxt(Path(tmp)/'trajectory.csv', names=True, delimiter=',').reshape(-1,64)
    data = np.genfromtxt(Path(tmp)/'observables.csv', names=True, delimiter=',')
    box = meta['box']; ii,jj = np.triu_indices(64,1)
    rc=2.5; fc=24*(2/rc**13-1/rc**7); uc=4*(rc**-12-rc**-6)
    for k in [0,50,100]:
        p=np.column_stack((frames[k]['x'],frames[k]['y']))
        delta=p[ii]-p[jj]; delta-=box*np.round(delta/box)
        distance=np.linalg.norm(delta,axis=1); distance=distance[distance<rc]
        energy=np.sum(4*(distance**-12-distance**-6)-uc+(distance-rc)*fc)
        virial=np.sum((24*(2*distance**-13-distance**-7)-fc)*distance)
        np.testing.assert_allclose(energy,data['potential'][k],atol=1e-10)
        np.testing.assert_allclose((2*data['kinetic'][k]+virial)/(2*box**2),data['pressure'][k],atol=1e-10)
    for wrapped, unwrapped in [('x','xu'),('y','yu')]:
        np.testing.assert_allclose(np.mod(frames[unwrapped],box),frames[wrapped],atol=1e-10)
    assert np.any(abs(frames['xu']-frames['x'])>box/2) or np.any(abs(frames['yu']-frames['y'])>box/2)
print('PASS: independent pair energy/virial and periodic unwrapped coordinates.')
