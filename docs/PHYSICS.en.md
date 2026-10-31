# Physical model and numerical analysis

[Українська](PHYSICS.uk.md) · [Overview](../README.md) · [Results](FIGURES.en.md)

## 1. Problem formulation

This study examines identical classical particles in a two-dimensional periodic domain. Its purpose is to relate microscopic interactions to energy, temperature, pressure and spatial correlations, while assessing numerical accuracy and computational cost. Each particle has mass $`m`$ and no internal degrees of freedom. Interactions are pairwise, central and conservative. External forces and heat exchange are absent.

$$A=L^2,$$

$$\rho=\frac{N}{A},$$

$$L=\sqrt{\frac{N}{\rho}}.$$

Here $`N`$ is particle count, $`L`$ is box length, $`A`$ is area and $`\rho`$ is number density. Planar geometry determines motion, statistical normalization and pressure dimensions. The Lennard–Jones potential is adopted as an interaction model; its form is not derived from two-dimensional geometry. “Gas” refers to the state under study, rather than an ideal system of noninteracting particles.

## 2. Hamiltonian and equations of motion

$$H(\lbrace \mathbf r_i,\mathbf p_i\rbrace )=\sum_{i=1}^{N}\frac{\mathbf p_i^2}{2m}+\sum_{i\lt j}u_s(r_{ij})=K+U.$$

$$\dot{\mathbf r}_i=\frac{\partial H}{\partial\mathbf p_i}=\frac{\mathbf p_i}{m},$$

$$\dot{\mathbf p}_i=-\frac{\partial H}{\partial\mathbf r_i}=\sum_{j\ne i}\mathbf F_{ij}.$$

The symbol $`u_s`$ denotes the cutoff-modified energy of one pair; $`U`$ denotes total potential energy. The distance $`r_{ij}`$ uses the nearest periodic image. Summation over $`i\lt j`$ counts each pair once. The implementation evaluates both pair directions, so its energy accumulation includes a factor of $`1/2`$.

## 3. Lennard–Jones potential and force

$$u(r)=4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12}-\left(\frac{\sigma}{r}\right)^6\right].$$

The parameter $`\epsilon`$ sets the energy scale, while $`\sigma`$ is the zero-potential distance before cutoff modification. The twelfth-power term describes steep short-range repulsion; the sixth-power term describes dispersion attraction. The exponent 12 is a convenient approximation to electronic repulsion.

$$\frac{du}{dr}=4\epsilon\left(-\frac{12\sigma^{12}}{r^{13}}+\frac{6\sigma^6}{r^7}\right),$$

$$F(r)=-\frac{du}{dr}=\frac{24\epsilon}{r}\left[2\left(\frac{\sigma}{r}\right)^{12}-\left(\frac{\sigma}{r}\right)^6\right].$$

$$\mathbf F_{ij}=F(r_{ij})\frac{\mathbf r_{ij}}{r_{ij}},$$

$$\mathbf r_{ij}=\mathbf r_i-\mathbf r_j\quad\text{(nearest image)}.$$

A positive radial force is repulsive. Setting the force to zero gives the minimum of the unmodified potential:

$$2\left(\frac{\sigma}{r_m}\right)^6=1,$$

$$r_m=2^{1/6}\sigma,$$

$$u(r_m)=-\epsilon.$$

Repulsion dominates for $`r\lt r_m`$ and attraction for $`r\gt r_m`$. The pair-energy minimum is not identical to the peak position of the pair distribution, which also depends on temperature, density and surrounding particles.

## 4. Cutoff and continuity of the interaction

The interaction radius is $`r_c=2.5\sigma`$. A linear correction to the potential ensures continuity of both energy and force:

$$u_s(r)=u(r)-u(r_c)+(r-r_c)F(r_c),\quad r\lt r_c.$$

$$u_s(r)=0,\quad r\ge r_c.$$

$$F_s(r)=-\frac{du_s}{dr}=F(r)-F(r_c),\quad r\lt r_c.$$

$$F_s(r)=0,\quad r\ge r_c.$$

$$u_s(r_c^-)=u_s(r_c^+)=0,$$

$$F_s(r_c^-)=F_s(r_c^+)=0.$$

The correction changes forces inside the cutoff, so the studied Hamiltonian differs from the full infinite-range LJ model. Energy conservation is evaluated for $`u_s`$ itself. The second derivative can be discontinuous at the cutoff. No long-range tail correction is applied.

## 5. Reduced units

The characteristic speed is set by the energy and mass scales; the characteristic time is the time to cover distance $`\sigma`$ at that speed:

$$v_0=\sqrt{\frac{\epsilon}{m}},\qquad\tau_0=\frac{\sigma}{v_0}.$$

$$\tau_0=\sigma\sqrt{\frac{m}{\epsilon}},$$

A reduced quantity is the physical quantity divided by its scale: $`r^{\ast}=1`$ means distance $`\sigma`$, and $`t^{\ast}=1`$ means time $`\tau_0`$.

$$r^{\ast}=\frac{r}{\sigma},$$

$$t^{\ast}=\frac{t}{\tau_0},$$

$$v^{\ast}=v\sqrt{\frac{m}{\epsilon}}.$$

$$T^{\ast}=\frac{k_BT}{\epsilon},$$

$$P^{\ast}=\frac{P\sigma^2}{\epsilon},$$

$$\rho^{\ast}=\rho\sigma^2,$$

$$D^{\ast}=\frac{D\tau_0}{\sigma^2}.$$

The numerical model sets $`m=\sigma=\epsilon=k_B=1`$; stars are omitted below. Pressure has units of energy per area. Conversion to kelvin or seconds would require material parameters, which are not assigned in this study. Temperatures and times in the figures are reduced quantities.

## 6. Periodic boundaries

$$\Delta r_{ij,\alpha}=r_{i,\alpha}-r_{j,\alpha}-L\mathrm{round}\left(\frac{r_{i,\alpha}-r_{j,\alpha}}{L}\right),$$

$$\alpha\in\lbrace x,y\rbrace .$$

The nearest-image displacement is computed component by component. The condition $`r_c\lt L/2`$ prevents simultaneous interaction with multiple images of the same particle. Display coordinates lie in $`[0,L)`$:

$$r_{i,\alpha}\leftarrow r_{i,\alpha}-L\left\lfloor\frac{r_{i,\alpha}}{L}\right\rfloor.$$

Continuous coordinates $`\mathbf R_i`$ are accumulated separately and are not wrapped after a boundary crossing. They determine physical displacement in MSD. Wrapped coordinates would introduce artificial jumps of one box length.

## 7. Conservation laws

$$\mathbf F_{ij}=-\mathbf F_{ji},$$

$$\frac{d}{dt}\sum_i m\mathbf v_i=\sum_{i\lt j}(\mathbf F_{ij}+\mathbf F_{ji})=\mathbf 0.$$

Differentiate kinetic energy and substitute Newton's second law:

$$\frac{dK}{dt}=\sum_i m\mathbf v_i\cdot\dot{\mathbf v}_i,\qquad m\dot{\mathbf v}_i=\mathbf F_i.$$

$$\frac{dK}{dt}=\sum_i\mathbf v_i\cdot\mathbf F_i,$$

For potential energy, use its dependence on positions:

$$\frac{dU}{dt}=\sum_i\nabla_i U\cdot\dot{\mathbf r}_i,\qquad\nabla_i U=-\mathbf F_i.$$

$$\frac{dU}{dt}=-\sum_i\mathbf v_i\cdot\mathbf F_i,$$

$$\frac{dH}{dt}=0.$$

The exact dynamics therefore conserves energy and linear momentum. Numerical energy has discretization error, while momentum has roundoff error. Ordinary total angular momentum expressed in periodic-box coordinates is not treated as an invariant: a periodic square lacks the continuous rotational symmetry of free space. Fixed $`N`$, $`A`$ and energy give the dynamics its microcanonical character; their conservation alone does not establish ergodicity.

## 8. Initial conditions and kinetic temperature

Initial coordinates lie on a square grid to avoid excessively close pairs. Gaussian random components form trial velocities $`\mathbf w_i`$. Center-of-mass motion is removed by subtracting their mean:

$$\overline{\mathbf w}=\frac1N\sum_i\mathbf w_i,$$

$$\mathbf v_i^{(0)}=\lambda(\mathbf w_i-\overline{\mathbf w}),$$

Choose the rescaling so that the initial kinetic energy corresponds to $`T_0`$:

$$K_0=\frac m2\sum_i\left\Vert\mathbf v_i^{(0)}\right\Vert^2=\frac{m\lambda^2}{2}\sum_i\left\Vert\mathbf w_i-\overline{\mathbf w}\right\Vert^2=(N-1)k_BT_0.$$

$$\lambda=\sqrt{\frac{(2N-2)k_BT_0}{m\sum_i|\mathbf w_i-\overline{\mathbf w}|^2}}.$$

The two zero-momentum constraints leave $`f=2N-2`$ thermal degrees of freedom. The kinetic temperature estimator follows from equipartition:

$$\langle K\rangle=\frac{f}{2}k_BT,$$

$$T_{\mathrm{kin}}(t)=\frac{2K(t)}{f k_B}=\frac{K(t)}{(N-1)k_B}.$$

Initial rescaling sets $`T_0`$ but does not make the configuration equilibrated. During relaxation, decreasing potential energy is accompanied by increasing kinetic energy. Temperature is not held fixed. $`T_{\mathrm{kin}}`$ is an instantaneous mechanical estimator whose interpretation as thermodynamic temperature relies on equilibrium properties.

## 9. Velocity Verlet integration

For time step $`h`$, a second-order position expansion and a trapezoidal approximation to the acceleration integral give:

$$\mathbf R_i(t+h)=\mathbf R_i(t)+h\mathbf v_i(t)+\frac{h^2}{2}\mathbf a_i(t)+O(h^3).$$

$$\mathbf v_i(t+h)=\mathbf v_i(t)+\frac h2[\mathbf a_i(t)+\mathbf a_i(t+h)]+O(h^3).$$

Dropping the remainder terms and splitting the velocity update into two half steps gives the scheme:

$$\mathbf v_{i,n+1/2}=\mathbf v_{i,n}+\frac{h}{2}\mathbf a_{i,n},$$

$$\mathbf R_{i,n+1}=\mathbf R_{i,n}+h\mathbf v_{i,n+1/2},$$

$$\mathbf r_{i,n+1}=\mathrm{wrap}(\mathbf R_{i,n+1}),$$

$$\mathbf a_{i,n+1}=\frac1m\sum_{j\ne i}\mathbf F_{ij}(\mathbf r_{n+1}),$$

$$\mathbf v_{i,n+1}=\mathbf v_{i,n+1/2}+\frac{h}{2}\mathbf a_{i,n+1}.$$

In the implementation, wrapped and unwrapped coordinates receive the same displacement; the equations above are an equivalent mathematical description. The method is a symmetric splitting into velocity and position updates. For sufficiently smooth forces, local state error is $`O(h^3)`$ and global error over a fixed time interval is $`O(h^2)`$. Symplecticity explains favorable long-time energy behavior, not exact energy conservation.

$$\delta_E(t)=\frac{E(t)-E(0)}{\max(|E(0)|,K(0))},$$

$$e(h)=\max_{0\le t\le t_f}|\delta_E(t)|.$$

$$e(h)\approx Ch^p,$$

$$p_{h,h/2}=\frac{\ln[e(h)/e(h/2)]}{\ln 2}.$$

The normalization avoids division by near-zero total energy. Runs use the same physical duration and initial state. The logarithmic regression over four step sizes gives a slope of 1.98528. This is the observed order of the energy diagnostic in the studied regime, rather than a distance to an exact trajectory.

## 10. Pressure and pair virial

The mechanical pressure tensor combines momentum transport by particle motion and the contribution of pair forces. In the periodic system, pair vectors use nearest images:

$$P_{\alpha\beta}=\frac1A\left[\sum_i m v_{i,\alpha}v_{i,\beta}+\sum_{i\lt j}r_{ij,\alpha}F_{ij,\beta}\right].$$

Add the diagonal tensor components; the sum of squared velocity components is the squared speed:

$$P_{xx}+P_{yy}=\frac1A\left[\sum_i m(v_{i,x}^2+v_{i,y}^2)+\sum_{i\lt j}\mathbf r_{ij}\cdot\mathbf F_{ij}\right]=\frac{2K+W}{A}.$$

$$P=\frac{P_{xx}+P_{yy}}2=\frac{2K+W}{2A},$$

Because the force is central, the dot product simplifies:

$$\mathbf r_{ij}\cdot\mathbf F_{ij}=\mathbf r_{ij}\cdot\left[F_s(r_{ij})\frac{\mathbf r_{ij}}{r_{ij}}\right]=r_{ij}F_s(r_{ij}).$$

$$W=\sum_{i\lt j}\mathbf r_{ij}\cdot\mathbf F_{ij}=\sum_{i\lt j}r_{ij}F_s(r_{ij}).$$

The factor of 2 in the denominator is spatial dimension. Repulsion contributes positively and attraction negatively. Without interactions, mechanical pressure is $`K/A`$. With center-of-mass momentum removed, the kinetic contribution under the adopted temperature convention is:

$$K=(N-1)k_BT_{\mathrm{kin}},\qquad\frac{N-1}{A}=\frac{N}{A}\left(1-\frac1N\right)=\rho\left(1-\frac1N\right).$$

$$P_{\mathrm{kin}}=\frac{K}{A}=\frac{N-1}{A}k_BT_{\mathrm{kin}}=\rho k_BT_{\mathrm{kin}}\left(1-\frac1N\right).$$

The usual ideal-gas form is recovered as $`N`$ becomes large. This study uses the direct mechanical estimator without replacing $`N-1`$ by $`N`$.

## 11. Two-dimensional speed distribution

In the canonical equilibrium reference without a total-momentum constraint, velocity components are independent and have Gaussian densities:

$$f_\alpha(v_\alpha)=\sqrt{\frac{m}{2\pi k_BT}}\exp\left(-\frac{mv_\alpha^2}{2k_BT}\right),\qquad\alpha\in\lbrace x,y\rbrace.$$

$$f(v_x,v_y)=f_x(v_x)f_y(v_y).$$

$$f(v_x,v_y)=\frac{m}{2\pi k_BT}\exp\left[-\frac{m(v_x^2+v_y^2)}{2k_BT}\right].$$

In polar velocity coordinates, the area element is $`v\thinspace dv\thinspace d\theta`$. Angular integration gives the speed probability density:

$$v_x=v\cos\theta,\qquad v_y=v\sin\theta,\qquad dv_x\,dv_y=v\,dv\,d\theta.$$

$$p(v)=\int_0^{2\pi}f(v\cos\theta,v\sin\theta) v d\theta=\frac{mv}{k_BT}\exp\left(-\frac{mv^2}{2k_BT}\right),$$

$$v\ge0.$$

$$\int_0^\infty p(v) dv=1,$$

Find the most probable speed from the maximum of the density:

$$\frac{dp}{dv}=\frac{m}{k_BT}\exp\left(-\frac{mv^2}{2k_BT}\right)\left(1-\frac{mv^2}{k_BT}\right)=0.$$

$$v_{\mathrm{mode}}=\sqrt{\frac{k_BT}{m}},$$

Find the mean speed as the expectation under the same density:

$$\langle v\rangle=\int_0^\infty v p(v)\,dv=\frac{m}{k_BT}\int_0^\infty v^2\exp\left(-\frac{mv^2}{2k_BT}\right)\,dv.$$

$$\langle v\rangle=\sqrt{\frac{\pi k_BT}{2m}},$$

For the root-mean-square speed, add the energies of the two components:

$$\frac m2\langle v^2\rangle=\frac m2\langle v_x^2+v_y^2\rangle=\frac12 k_BT+\frac12 k_BT=k_BT.$$

$$\sqrt{\langle v^2\rangle}=\sqrt{\frac{2k_BT}{m}}.$$

The plotted curve uses the mean kinetic temperature of the second half and serves as a thermodynamic-limit reference. A finite system with fixed momentum has corrections. Even in a canonical system with zero total momentum, one-component variance is $`(1-1/N)k_BT/m`$. In the microcanonical state, the energy constraint further affects the exact distribution. These corrections are not included in the displayed Maxwell curve.

## 12. Pair distribution and local structure

The radial function $`g(r)`$ compares pair counts at separation $`r`$ with an independent uniform arrangement. A bin $`[a_k,b_k)`$ uses annular area rather than a spherical shell:

$$\Delta A_k=\pi(b_k^2-a_k^2),$$

$$\langle H_k\rangle_{\mathrm{ideal}}=M\frac{N(N-1)}2\frac{\Delta A_k}{A},\qquad g_k=\frac{H_k}{\langle H_k\rangle_{\mathrm{ideal}}}.$$

$$g_k=\frac{H_k}{M \dfrac{N(N-1)}2 \dfrac{\Delta A_k}{A}}.$$

Here $`H_k`$ is the total unordered-pair count in a histogram bin over $`M`$ frames. This normalization gives $`g=1`$ for independent uniform points at finite $`N`$. The maximum separation $`L/2`$ keeps complete annuli within the nearest-image domain. With the same normalization, the mean number of neighbors within radius $`R`$ is:

$$n_{\mathrm{coord}}(R)=\frac{N-1}{A}\int_0^R 2\pi r g(r) dr.$$

Coordination number is given as a structural relation and was not measured separately. The figure reports $`g(r)`$: core depletion reflects repulsion, while the first peak marks preferred neighbor separations. The geometric normalization follows the two-dimensional RDF definition.

## 13. Displacement, velocity correlation and diffusion

For $`M`$ equally spaced frames in the analysis interval, lag index $`\ell`$ and sampling interval $`\Delta t_s`$, averages run over particles and valid time origins:

$$\mathrm{MSD}(\ell\Delta t_s)=\frac{1}{N(M-\ell)}\sum_{k=0}^{M-\ell-1}\sum_{i=1}^{N}\left\Vert\mathbf R_i(t_{k+\ell})-\mathbf R_i(t_k)\right\Vert^{2}.$$

$$C_v(\ell\Delta t_s)=\frac{1}{N(M-\ell)}\sum_{k=0}^{M-\ell-1}\sum_{i=1}^N\mathbf v_i(t_{k+\ell})\cdot\mathbf v_i(t_k).$$

The relation between these observables follows by integrating velocity. For a stationary process, the correlation depends only on the time difference:

$$\Delta\mathbf R(t)=\int_0^t\mathbf v(s) ds,$$

$$\mathrm{MSD}(t)=\int_0^t\int_0^t\langle\mathbf v(s)\cdot\mathbf v(s\prime)\rangle\,ds\,ds\prime,\qquad\langle\mathbf v(s)\cdot\mathbf v(s\prime)\rangle=C_v(|s-s\prime|).$$

$$\mathrm{MSD}(t)=2\int_0^t(t-s)C_v(s) ds.$$

$$\frac{d\mathrm{MSD}}{dt}=2\int_0^t C_v(s) ds,$$

$$\frac{d^2\mathrm{MSD}}{dt^2}=2C_v(t).$$

At short times velocity changes little. The mean-square displacement grows quadratically with time:

$$\mathrm{MSD}(t)\approx\langle v^2\rangle t^2.$$

The mean squared speed sets the coefficient of ballistic growth. For normal diffusion in two dimensions, the dependence becomes linear:

$$\mathrm{MSD}(t)\sim4Dt,$$

$$D=\lim_{t\to\infty}\frac{\mathrm{MSD}(t)}{4t}=\frac12\int_0^\infty C_v(s) ds.$$

The last relations require stationarity and a convergent long-time regime. The Einstein and Green–Kubo approaches describe the same transport property. This study reports finite-time MSD and VACF, but does not estimate $`D`$: a maximum lag of 1.5 does not establish an asymptotic regime. Overlapping time origins are correlated, so their count is not the number of independent measurements.

## 14. Computational methods and model validation

The sequential method evaluates all ordered pairs in $`O(N^2)`$. OpenMP distributes particles among threads: each thread writes its particles’ accelerations, while reductions accumulate energy and virial. The grid method divides the box into cells at least $`r_c`$ wide and checks the current and neighboring cells. Repeated periodic indices in small boxes are counted once. Neighbor search introduces no force approximation: all three methods implement the same Hamiltonian.

At fixed density and bounded cell occupancy, grid cost is approximately $`O(N)`$; strong clustering can restore the quadratic case. Validation covers trajectory agreement among methods, momentum conservation, energy error, independent pair-energy and pressure recomputation, and consistency of wrapped and unwrapped coordinates. Structural and transport analysis uses the second half of one run; this choice does not replace a statistical demonstration of equilibrium.
