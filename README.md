# Molecular dynamics of a two-dimensional gas

[Українська](README.uk.md) · [Physical model](docs/PHYSICS.en.md) · [Results and discussion](docs/FIGURES.en.md)

This study examines a two-dimensional particle system with Lennard–Jones interactions. It considers kinetic–potential energy exchange, the speed distribution, spatial correlations and the loss of velocity memory. Numerical accuracy is assessed through energy conservation and convergence under time-step refinement. Computational performance is compared for sequential pair summation, OpenMP and cell-based neighbor search.

## Model

Particles move in a square domain with periodic boundaries. Pair interactions are limited to $`r_c=2.5\sigma`$; a linear correction makes both potential and force vanish continuously at the cutoff. The equations of motion are integrated with velocity Verlet. External forces and a thermostat are absent, so changes in kinetic temperature arise from energy exchange between motion and interactions.

$$H=\sum_{i=1}^{N}\frac{m|\mathbf v_i|^2}{2}+\sum_{i\lt j}u_s(r_{ij}),$$

$$m\ddot{\mathbf r}_i=\sum_{j\ne i}\mathbf F_{ij}.$$

The simulation is written in C++17 with OpenMP. Trajectory analysis, figures and animation use Python. All physical quantities are expressed in Lennard–Jones reduced units.

## Experimental conditions

| Quantity | Value |
|---|---:|
| Particle count | 256 |
| Number density | 0.3 |
| Initial temperature | 1 |
| Integration step | 0.002 |
| Physical duration | 6 |
| Saved frames | 201 |
| Random seed | 42 |

Initial positions form a regular grid. Random velocities have zero total momentum and a prescribed kinetic temperature. Structural and transport observables are calculated over the second half of the trajectory.

## Particle motion

![Particle motion](output/plots/gas_animation.gif)

Color represents speed, while the vertical line on the temperature plot marks the current time. Crossing a box edge corresponds to motion into a neighboring periodic image. Interparticle forces gradually disrupt the initial spatial order.

## Energy and temperature

![Energy, temperature and speed distribution](output/plots/diagnostics.png)

Decreasing potential energy is accompanied by increasing kinetic energy. Total energy remains close to its initial value: the maximum normalized error over saved frames is $`7.8831\times10^{-5}`$. Mean kinetic temperature over times 3 to 6 is 1.47049. The speed distribution is compared with the two-dimensional Maxwell law at this temperature.

## Structure and transport

![Pair distribution, pressure, MSD and VACF](output/plots/structure_transport.png)

The pair distribution has a strongly depleted short-range region and a first maximum near $`r=1.095\sigma`$. Mean mechanical pressure over the second half is 0.546259. MSD is close to ballistic scaling at short lags, while decaying velocity autocorrelation indicates gradual loss of motion memory. The sampled duration characterizes transient dynamics but does not establish a long-time diffusion coefficient.

## Accuracy and performance

![Time-step convergence](output/plots/convergence.png)

Four integration steps are compared for a separate 64-particle system at the same final physical time. The logarithmic slope of maximum energy error against step size is 1.98528, consistent with second-order integration.

![Algorithm scaling](output/benchmarks/scaling.png)

At 1024 particles, the cell grid reduces integration time by approximately 30.19 relative to sequential summation. Four OpenMP threads give a speedup of 3.12 relative to one thread of the same implementation. Measurements are medians of three repetitions and exclude trajectory CSV writing.

## Conclusions

The system exhibits energy redistribution during relaxation, short-range positional order and gradual loss of velocity correlation. Time-step refinement is consistent with the expected integration order. Local neighbor search provides the largest performance gain among the studied implementations while retaining the same interaction model.

The conclusions concern a finite two-dimensional system and one initial state. They do not determine a phase diagram or asymptotic transport coefficients. Full derivations appear in the [physical model](docs/PHYSICS.en.md), while figure normalization and numerical evidence are presented in [results and discussion](docs/FIGURES.en.md).
