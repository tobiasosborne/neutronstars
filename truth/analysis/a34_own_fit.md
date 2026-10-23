**Abstract.** A34 is an inference problem with substantial non-identifiability, not merely an optimizer to replace Ho et al.’s parameter table. An owned calculation must predict instrument counts and calibrated photometry from the same surface radiance used by the renderer. It must also distinguish fitting the observations from independently predicting them. The calculation below reproduces the optical excess: the Burwitz X-ray blackbody predicts a reddened flux at $2500\,\mathrm{\AA}$ smaller than Kaplan’s observationally fitted flux by **7.10055**, consistent with the published $7.1\pm0.3$. Distance rescaling cannot remove this discrepancy. However, neither a successful broadband fit nor the small pulse amplitude determines a unique surface map, redshift, or flyby appearance. Milestone 2 should establish an owned instrument-level harness and a documented single-temperature blackbody failure, while leaving absorption, extinction, unavailable measurements, and unresolved geometry explicitly conditional.

# 1. Corrections to the premise and the project plan

## 1.1 Ho et al. did not perform the proposed unrestricted joint fit

The distinction matters for both reproduction and circularity.

Ho et al. explicitly state that they do not “fit for all of the parameters in a proper sense” and instead fit a “subset of parameters while varying others manually.” For each chosen magnetic field and atmospheric column, they fit LETG data for $N_{\rm H}$, $T_\infty$, normalization expressed as $R_\infty$, and $z_g$.

Source: [refs/ho_2007_rxj1856_thin_atm.tex:580](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:580), anchors `subset of parameters while varying others manually.` and `the normalization (parameterized by \rinfty), and the redshift $\zg$`.

Their Table 1 separates **model parameters** from **fit results**. The trial fields are $3\times10^{12}$ and $4\times10^{12}\,\mathrm G$, both with $y_{\rm H}=1.2\,\mathrm{g\,cm^{-2}}$. Thus these particular $B$ and $y_{\rm H}$ values are not measurements with marginalized errors.

Source: `refs/ho_2007_rxj1856_thin_atm.tex:627`, anchors `Model Parameters` and `$\ycol$ (g cm$^{-2}$) & 1.2 & 1.2`.

The optical comparison follows the LETG fit. They fix $A_V=0.12$ after finding extinction poorly constrained, and report that their model:

> `underpredicts the optical/UV data by a factor of 15$-$20`

Source: `refs/ho_2007_rxj1856_thin_atm.tex:685`, with that anchor at line 700. Their highest-quality F606W point is above the model by $15\%$, versus a quoted photometric uncertainty of $3\%`; the broader agreement includes extinction, model, and X-ray-normalization uncertainties.

Source: `refs/ho_2007_rxj1856_thin_atm.tex:709`, anchor `is 3\% (photometric uncertainty)`.

Consequently:

- Ho’s values are legitimate **conditional reproduction targets**.
- They are not the result of a complete exploration of the parameter vector proposed for A34.
- Their reported agreement does not establish PRD R5’s literal criterion of reproducing magnitudes “within published errors” if that means the photometric errors alone.
- Reproducing their table with improved physics is not necessarily the correct acceptance test.

The last point has a concrete example: Potekhin et al. later state that Ho et al. “retained only the emission terms” in the condensed-surface boundary condition, whereas the later treatment includes reflected radiation and mode conversion.

Source: [refs/potekhin_2012_condensed_surface.tex:1104](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:1104), anchor `retained only the emission terms`.

## 1.2 R5 is narrower than the scientific inference problem

The current PRD R5 specifies HST optical/UV photometry at the adopted distance. It does **not** explicitly require an acceptable Chandra fit or pulse amplitude.

That omission permits a weak result: choose a blackbody normalization and temperature to fit the optical data alone, and potentially satisfy R5 while violating the X-ray observations. Conversely, “blackbody will fail R5” is justified for the **single-temperature blackbody constrained by X-rays**, not for every blackbody surface map.

The local literature explicitly permits a two-temperature blackbody description. Burwitz gives a hot component and a cooler component responsible for the optical flux.

Source: [refs/burwitz_2003_rxj1856_xray.tex:65](/home/tobias/Projects/neutronstars/refs/burwitz_2003_rxj1856_xray.tex:65), anchors `two-component blackbody model` and `the rest of the neutron star surface responsible for the optical flux.`

**Proposed acceptance distinction:**

1. **R5 photometric consistency:** the rendered, instrument-folded photometry agrees with the specified measurements.
2. **Joint observational consistency:** the same physical parameter set explains LETG counts, HST measurements, and timing.
3. **Predictive validation:** a declared subset of those observations was not used to fit parameters or select the model.
4. **Ownership:** every computational dependency satisfies the project’s provenance rule.

These are separate properties. Passing one does not establish the others.

# 2. Parameters and the common physical forward model

## 2.1 Use a nonredundant parameterization

For surface-model family $\mathcal M$, choose

$$
\theta_{\mathcal M}=
\left(
\rho_\infty,z_g,d,\boldsymbol\eta_T,\boldsymbol\eta_B,
y_{\rm H},i,\beta,\phi_0,
N_{\rm H},A_V,\boldsymbol\eta_{\rm ISM},
\boldsymbol\eta_{\rm cal},\boldsymbol\eta_{\rm bg}
\right),
$$

where:

- $\rho_\infty=R_\infty/d$ is the apparent angular radius of the **whole stellar surface**;
- $\boldsymbol\eta_T$ specifies a temperature map;
- $\boldsymbol\eta_B$ specifies field strength and geometry;
- $y_{\rm H}$ is fitted only for a finite hydrogen layer;
- $i$ and $\beta$ are viewing inclination and magnetic-axis obliquity;
- the final groups describe absorption, extinction, calibration, and background.

This is a **proposed modeling convention**, not a claim that the data determine every component.

Ho prints

$$
T_\infty=\frac{T_{\rm eff}}{1+z_g},
\qquad
R_\infty=R(1+z_g),
\qquad
1+z_g=\left(1-\frac{2GM}{Rc^2}\right)^{-1/2}.
$$

Source: `refs/ho_2007_rxj1856_thin_atm.tex:128`, anchors `$\tinfty=\teff/(1+\zg)$` and `$(1+\zg)=(1-2GM/Rc^2)^{-1/2}$`.

**DERIVED HERE:** defining $u=2GM/(Rc^2)$ gives

$$
u=1-(1+z_g)^{-2},
\qquad
R=\frac{d\rho_\infty}{1+z_g},
\qquad
GM=\frac{c^2Ru}{2}.
$$

Therefore $M,R,R_\infty,z_g$ must not all be independently varied.

Surface gravity must be recomputed consistently:

$$
g_s=\frac{GM(1+z_g)}{R^2}.
$$

Source: `refs/ho_2007_rxj1856_thin_atm.tex:244`, anchor `surface gravity $g$~$[=(1+\zg)GM/R^2]$`.

A critical distinction is that an X-ray blackbody **emitting radius** need not equal the apparent radius of the entire star. Equating them would turn a fitted hot-region normalization into an unsupported spacetime constraint.

## 2.2 From local intensity to Earth flux

Adopt a spherical, quasi-static surface at each rotational phase for this Earth-distance calculation. The full flyby subsequently uses the camera’s actual frequency shift and tetrad.

The radiative invariant is locally sourced:

> `Liouville's theorem, $I_\nu \propto \nu^3$`.

Source: [refs/james_2015_interstellar.tex:978](/home/tobias/Projects/neutronstars/refs/james_2015_interstellar.tex:978).

For a static surface, set $g=(1+z_g)^{-1}$. Since $E=h\nu$, the same transformation applies to intensity per unit photon energy.

**DERIVED HERE:**

$$
I_E^\infty(E)=g^3 I_E^s(E/g).
$$

Integrating rays over the apparent image gives

$$
F_E^0(E,\phi;\theta_{\mathcal M})
=
\frac{1}{d^2}
\int_{\rm image}
g^3 I^s_{E/g}
\left[
\Omega_s(b,\chi,\phi),\hat{\boldsymbol k}_s;
\theta_{\mathcal M}
\right]
b\,db\,d\chi.
$$

The image-plane measure follows directly from the local expression

> `subtends a solid angle $\dd\Omega=b\,\dd b\,\dd\phi/D^2$`.

Source: `refs/beloborodov_2002_ray_tracing.tex:573`.

Integrating in image coordinates also avoids silently dropping the absolute value of a surface-to-image Jacobian.

For a uniform isotropic blackbody, the grazing-ray relation

$$
\sin\alpha=\frac{b}{R}\sqrt{1-\frac{r_g}{R}}
$$

is printed in `refs/beloborodov_2002_ray_tracing.tex:556`, anchor `\sin\alpha=\frac{b}{R}\sqrt{1-\frac{r_g}{R}}`.

For a surface outside the photon sphere, this gives $b_{\max}=R_\infty$. Consequently,

$$
F_E^0=\pi\rho_\infty^2 B_E(T_\infty).
$$

This is **DERIVED HERE** by integrating a constant image intensity over the apparent disc. The factor $\pi$ is mandatory.

## 2.3 The three PRD R4 surface models

All three models should implement one interface:

$$
\mathcal I_{\mathcal M}:
(E,\hat{\boldsymbol k},T_s,\boldsymbol B,g_s,y_{\rm H},\ldots)
\longrightarrow
(I_E,Q_E,U_E,V_E).
$$

The fit consumes Stokes $I$; the renderer retains the full output.

### Blackbody

Use isotropic unpolarized emission,

$$
I_\lambda^s=B_\lambda(T_s),
\qquad Q_\lambda=U_\lambda=V_\lambda=0,
$$

with

$$
B_\lambda(T)=
\frac{2hc^2}{\lambda^5}
\frac{1}{\exp[hc/(\lambda k_{\rm B}T)]-1}.
$$

Source: [refs/marr_2011_planck_law.tex:42](/home/tobias/Projects/neutronstars/refs/marr_2011_planck_law.tex:42), anchor `B_\lambda =\frac{2hc^2}{\lambda^5}`.

Isotropy and zero polarization are the chosen blackbody model convention. A uniform temperature gives the minimal X-ray fit; a parameterized temperature map is a distinct, more flexible hypothesis.

### Optically thick hydrogen atmosphere

Compute the outgoing intensity from an atmosphere deep enough that the lower boundary is thermalized over the fitted energy range.

Ho’s description states that conventional models are “optically thick at all photon energies” and that their emerging spectrum is determined by the temperature profile and opacities.

Source: `refs/ho_2007_rxj1856_thin_atm.tex:313`, anchors `optically thick at all photon energies` and `determined by the temperature profile and opacities of the atmosphere.`

**DERIVED HERE:** let $y$ increase inward, $\mu>0$ denote an outward ray, $\kappa_j$ be extinction per mass, and $\eta_j$ the emission-plus-scattering source per mass. Starting from loss minus gain along a ray,

$$
\frac{dI_j}{ds}=-\rho\kappa_j I_j+\rho\eta_j,
\qquad
\frac{dy}{ds}=-\rho\mu,
$$

gives

$$
\mu\frac{\partial I_j}{\partial y}
=
\kappa_j I_j-\eta_j.
$$

The source couples directions and polarization modes. Opacities, ionization balance, hydrostatic structure, and radiative equilibrium must be computed consistently.

A fully ionized solver is defensible as a development limit. It is not an automatically adequate RX J1856 model. Ho explicitly concludes:

> `Therefore the atmospheric plasma must be partially ionized.`

Source: `refs/ho_2007_rxj1856_thin_atm.tex:268`, following the comparison of inferred temperatures with magnetic hydrogen binding energies.

### Thin hydrogen atmosphere over a condensed surface

Use the same transfer problem on a finite column $0\leq y\leq y_{\rm H}$, with the condensed material supplying thermal emission and reflection.

In compact notation, the printed boundary conditions are

$$
I^+_{E,j}
=
\frac12 J_j B_E(T_{\rm bot})
+\sum_k R_{jk}I^-_{E,k}.
$$

Source: `refs/potekhin_2012_condensed_surface.tex:1067`, with the defining text at line 1094:

> `coefficients of reflection with allowance for`

and at line 1095:

> `transformation of the incident mode $M'$ into the reflected mode`.

The atmosphere’s own downward radiation must participate in this boundary condition. Multiplying a blackbody by an emissivity curve after solving the atmosphere is not equivalent.

**DERIVED HERE — required limits:**

- As $y_{\rm H}\rightarrow0$, the outgoing solution approaches the bare condensed-surface solution.
- As optical depth becomes large over the entire fitted band, sensitivity to the bottom boundary disappears.
- For a prescribed isothermal, purely absorbing slab,

$$
I(0,\mu)
=
I_{\rm bot}e^{-\tau/\mu}
+B\left(1-e^{-\tau/\mu}\right),
$$

obtained by integrating $\mu\,dI/d\tau=I-B$.

These limits test the transfer and boundary implementation independently of a published fit.

# 3. The likelihood must live in observable space

## 3.1 LETG source and background counts

The interstellar attenuation law is

$$
F_E^{\rm obs}(E)
=
e^{-N_{\rm H}\sigma_{\rm ISM}(E)}F_E^0(E).
$$

Source: [refs/wilms_allen_mccray_2000_ism.tex:129](/home/tobias/Projects/neutronstars/refs/wilms_allen_mccray_2000_ism.tex:129), anchor `I_{\rm obs}(E) = {\rm e}^{-\sigism(E)\,N_{\rm H}}`.

For observation $r$, detector channel $i$, and diffraction order $o$, the expected source counts are **DERIVED HERE** by multiplying incident photon flux by collecting area, redistribution probability, and exposure:

$$
s_{ri}(\theta)
=
t_r\sum_o
\int
R_{rio}(E)\,A_{ro}(E)\,
e^{-N_{\rm H}\sigma_{\rm ISM}(E)}
\frac{\overline{F_E^0}(E;\theta)}{E}\,dE.
$$

Here $R$ is the RMF, $A$ the ARF, and the overbar denotes the relevant phase average. The energy in the photon-conversion denominator must use units consistent with $F_E$: dividing a per-keV energy flux by a numerical energy in keV without the keV-to-erg factor is wrong.

Ho generated responses for

> `orders $\pm 1$, $\pm 2$, and $\pm 3$`.

Source: `refs/ho_2007_rxj1856_thin_atm.tex:569`.

Thus a first-order-only response is a declared approximation requiring a convergence check, not the default reproduction of Ho.

Let $n_{ri}$ be counts in the source extraction and $k_{ri}$ counts in a background extraction. Define $\beta_{ri}$ as expected background-region counts and $\alpha_{ri}$ as the known transfer factor to the source extraction:

$$
n_{ri}\sim{\rm Poisson}(s_{ri}+\alpha_{ri}\beta_{ri}),
\qquad
k_{ri}\sim{\rm Poisson}(\beta_{ri}).
$$

Do not subtract the background and call the result Poisson.

**DERIVED HERE — Poisson likelihood and C-statistic.** Divide an exposure into $m$ independent intervals, each with event probability $\mu/m$. Then

$$
\lim_{m\rightarrow\infty}
{m\choose n}
\left(\frac{\mu}{m}\right)^n
\left(1-\frac{\mu}{m}\right)^{m-n}
=
\frac{e^{-\mu}\mu^n}{n!}.
$$

Taking minus twice its logarithm and subtracting the saturated-model value gives

$$
C(n\mid\mu)
=
2\left[\mu-n+n\ln\left(\frac n\mu\right)\right].
$$

Therefore

$$
C_X=
\sum_{ri}
\left[
C(n_{ri}\mid s_{ri}+\alpha_{ri}\beta_{ri})
+
C(k_{ri}\mid\beta_{ri})
\right].
$$

The zero-count expression is its analytic limit, not a numerical floor.

**Cash attribution: full source missing.** Cash 1979 is not locally available as a paper. Its bibliographic citation appears at `refs/vankerkwijk_kaplan_2008_rxj1856_timing.tex:518`, anchor `{Cash}, W. 1979, \apj, 228, 939`. The likelihood above is derived here; no claim about Cash’s original implementation is imported.

## 3.2 HST magnitudes

Use observed, reddened measurements. The forward model supplies

$$
F_\lambda^{\rm obs}
=
F_\lambda^0\,10^{-0.4A_\lambda}.
$$

This attenuation factor is **DERIVED HERE** from the magnitude definition: an extinction $A_\lambda$ is a positive increase in magnitude.

For a photon-counting throughput $S_b(\lambda)$, photon energy $hc/\lambda$ gives

$$
\dot N_b
=
A_{\rm tel}\int
F_\lambda^{\rm obs}S_b(\lambda)
\frac{\lambda}{hc}\,d\lambda.
$$

Equating this to the count rate from a constant-$F_\lambda$ reference spectrum yields the band-equivalent flux density

$$
\left\langle F_\lambda\right\rangle_b
=
\frac{\int \lambda S_b(\lambda)F_\lambda^{\rm obs}(\lambda)\,d\lambda}
{\int \lambda S_b(\lambda)\,d\lambda}.
$$

These are **DERIVED HERE**, conditional on $S_b$ representing photon-detection efficiency. The actual calibration files must specify whether collecting area and photon weighting are already included.

Then

$$
m_b^{\rm pred}
=
-2.5\log_{10}\left\langle F_\lambda\right\rangle_b-21.1,
$$

with flux density in $\mathrm{erg\,cm^{-2}\,s^{-1}\,\AA^{-1}}$.

Source: [refs/kaplan_2011_rxj1856_spectrum.tex:488](/home/tobias/Projects/neutronstars/refs/kaplan_2011_rxj1856_spectrum.tex:488), anchor `magnitude: ${\rm STMAG}=-2.5\log_{10} F_{\lambda}-21.1$ with`.

Adopting a Gaussian measurement model in magnitudes gives

$$
\chi_H^2
=
(\boldsymbol m^{\rm obs}-\boldsymbol m^{\rm pred})^\mathsf T
\Sigma_H^{-1}
(\boldsymbol m^{\rm obs}-\boldsymbol m^{\rm pred}).
$$

This follows directly from the Gaussian density. A fixed covariance permits omission of its normalization from optimization; a parameter-dependent covariance requires the $\ln\det\Sigma_H$ term.

Shared zeropoint and aperture uncertainties belong in correlated covariance terms or explicit nuisance parameters. They must not be counted twice or inflated after seeing a mismatch.

## 3.3 Pulse amplitude

The relevant local result is a background- and bias-corrected **sinusoidal fractional amplitude**:

$$
a_{\rm obs}=0.0096,
\qquad
\sigma_a=0.0006.
$$

Source: [refs/vankerkwijk_kaplan_2008_rxj1856_timing.tex:270](/home/tobias/Projects/neutronstars/refs/vankerkwijk_kaplan_2008_rxj1856_timing.tex:270), anchors `Correcting them for background and statistical bias` and `the best estimate of the amplitude is~$0.96\pm0.06$`.

For a model count-rate curve $C_r(\phi)$, define **DERIVED HERE**

$$
\overline C_r=\frac{1}{2\pi}\int_0^{2\pi}C_r(\phi)\,d\phi,
$$

$$
A_{c,r}=\frac1\pi\int_0^{2\pi}C_r(\phi)\cos\phi\,d\phi,
\qquad
A_{s,r}=\frac1\pi\int_0^{2\pi}C_r(\phi)\sin\phi\,d\phi,
$$

$$
a_{1,r}=
\frac{\sqrt{A_{c,r}^2+A_{s,r}^2}}{\overline C_r}.
$$

This equals the fractional semiamplitude for a sinusoid. It is not generally identical to $(C_{\max}-C_{\min})/(C_{\max}+C_{\min})$.

A summary-level constraint is

$$
\chi_a^2=
\frac{[a_{\rm model}-0.0096]^2}{0.0006^2}.
$$

However, $a_{\rm model}$ must reproduce the response-weighted estimator used for the observations. It cannot be an arbitrary bolometric pulsed fraction. The timing paper used differing event selections; for example, it lowered the PN imaging cutoff and selected events below $1\,\mathrm{keV}$.

Source: `refs/vankerkwijk_kaplan_2008_rxj1856_timing.tex:199`, anchor `default 150\,eV low-energy cutoff to 100\,eV`.

Prefer fitting the individual timing datasets with their responses and selections. If timing and spectral likelihoods reuse events, construct a joint event likelihood or a phase likelihood conditional on the spectral counts; do not assume the summaries are independent.

## 3.4 Joint objective and external information

For independent datasets, the proposed objective is

$$
Q(\theta)=
C_X+\chi_H^2+\chi_a^2
-2\ln L_{\rm astrometry}(d)
-2\ln p(\theta_{\rm nuisance}).
$$

With a nonconstant prior, minimizing this is MAP estimation, not pure maximum likelihood.

Walter’s independent astrometric distance is legitimate information. The authors explicitly adopt the individual-measurement error rather than treating several reductions of the same data as independent:

> `we adopt the mean error of the individual measurements, 0.8~mas, as`

Source: `refs/walter_2010_rxj1856_distance.tex:809`.

Their final interval is given by

> `The $\pm 1~\sigma$ confidence region encompasses 108 to 134~pc.`

Source: `refs/walter_2010_rxj1856_distance.tex:856`.

For R5, evaluate the declared reference render at the PRD’s adopted distance. For physical inference, propagate astrometric uncertainty rather than pretending that distance is exact.

There is no independent RX J1856 $N_{\rm H}$ measurement established by the supplied local sources. Burwitz’s and Ho’s columns come from fitting the same target spectra. They cannot become independent Gaussian priors on a refit of those spectra.

# 4. What the observations identify

## 4.1 X-rays constrain spectral curvature and angular normalization

For the uniform blackbody, write

$$
F_E^{\rm obs}
=
\mathcal A
\frac{2E^3}{h^3c^2}
\frac{e^{-N_{\rm H}\sigma(E)}}{e^{E/(k_{\rm B}T_\infty)}-1},
\qquad
\mathcal A=\pi\rho_\infty^2.
$$

This follows from the sourced Planck function by changing variables from frequency to energy.

Define

$$
x=\frac{E}{k_{\rm B}T_\infty},
\qquad
G(x)=\frac{x}{1-e^{-x}}.
$$

**DERIVED HERE:** logarithmic differentiation gives

$$
\delta\ln F_E
=
2\,\delta\ln\rho_\infty
+
G(x)\,\delta\ln T_\infty
-
\sigma(E)\,\delta N_{\rm H}.
$$

Thus a narrow band measures approximately one combination of normalization, temperature, and absorption. Broad energy coverage separates these parameters because their derivative shapes differ. The temperature–normalization covariance is strong, but not an exact degeneracy over a well-calibrated broad spectrum.

Using Burwitz’s $k_{\rm B}T_\infty=63.5\,\mathrm{eV}$, Python gives:

| Photon energy | $x$ | $\partial\ln F_E/\partial\ln T_\infty$ | $\delta\ln\rho_\infty/\delta\ln T_\infty$ at fixed flux and $N_{\rm H}$ |
|---:|---:|---:|---:|
| $200\,\mathrm{eV}$ | 3.14961 | 3.29067 | −1.64534 |
| $300\,\mathrm{eV}$ | 4.72441 | 4.76672 | −2.38336 |
| $500\,\mathrm{eV}$ | 7.87402 | 7.87701 | −3.93851 |

These are local spectral sensitivities, not a Fisher matrix inferred from the actual LETG counts.

To see the absorption covariance analytically, use a **declared local illustrative approximation**, $\sigma(E)\propto E^{-p}$ with $p>0$, away from edges. Then

$$
\frac{d\ln F_E}{d\ln E}
=
3-G(x)+pN_{\rm H}\sigma(E).
$$

Both increasing temperature and increasing absorption harden the observed spectrum. Consequently, along a roughly constant-hardness direction, an increased fitted column can be compensated by a decreased temperature, with normalization changing in turn. Actual edges and detector redistribution must replace this approximation in the fit.

Burwitz Table 2 gives the LETG results:

$$
N_{\rm H}=(0.95\pm0.03)\times10^{20}\,\mathrm{cm^{-2}},
$$

$$
k_{\rm B}T_\infty=63.5\pm0.2\,\mathrm{eV},
\qquad
R_{\rm bb,\infty}=(4.4\pm0.1)\left(\frac d{120\,\mathrm{pc}}\right)\mathrm{km}.
$$

Source: `refs/burwitz_2003_rxj1856_xray.tex:248`, anchor `0.95&$\pm$0.03 & 63.5&$\pm$0.2`; units appear in the preceding table header.

Python converts the printed errors into fractional scales:

- temperature: $0.3150\%$;
- emitting-radius normalization: $2.2727\%$;
- absorbing column: $3.1579\%$;
- squared angular normalization: approximately $4.5455\%$.

These are conditional fit uncertainties. They are not the uncertainty of the physical interpretation.

The same table gives EPIC-pn $N_{\rm H}=0.18\times10^{20}\,\mathrm{cm^{-2}}$, a factor **5.27778** below the LETG value. Burwitz attributes instrument-to-instrument differences to calibration systematics.

Source: `refs/burwitz_2003_rxj1856_xray.tex:249` and `:330`, anchor `are attributed to systematic uncertainties in the present`.

This historical comparison does not quantify current calibration uncertainty. It demonstrates why the small printed statistical errors cannot stand in for it.

## 4.2 Optical data measure a different temperature–area combination

Expanding the Planck denominator for $x\ll1$ gives, **DERIVED HERE**,

$$
B_\lambda
=
\frac{2ck_{\rm B}T}{\lambda^4}
\left(1-\frac{x}{2}+\cdots\right).
$$

Hence

$$
F_\lambda^{\rm obs}
\simeq
2\pi ck_{\rm B}\rho_\infty^2T_\infty
\lambda^{-4}10^{-0.4A_\lambda}.
$$

The optical amplitude primarily measures $\rho_\infty^2T_\infty$, whereas the X-ray tail is much more sensitive to temperature. Their intersection tests a one-temperature model.

For multiple blackbody regions, optical data instead constrain a weighted moment,

$$
F_\lambda^{\rm RJ}
\propto
\lambda^{-4}\sum_j \Omega_j^{\rm app}T_{\infty,j},
$$

where $\Omega_j^{\rm app}$ includes visibility and phase averaging. The Wien tail preferentially constrains the hottest regions. Optical data can reveal a cooler, larger emitting contribution without uniquely specifying its area and temperature.

For an illustrative two-component decomposition with optical excess $f$,

$$
R_{c,\infty}^2T_{c,\infty}
\simeq
(f-1)R_{h,\infty}^2T_{h,\infty}.
$$

Taking $f=7.3$, the Burwitz hot component at the adopted distance, and the illustrative $kT_{c,\infty}=26\,\mathrm{eV}$ printed by Ho, Python gives an equivalent cool-component radius of **17.6908 km**. This is a Rayleigh–Jeans projected-area calculation, not a new fit or a measurement of the stellar radius.

The cool-temperature example is sourced at `refs/ho_2007_rxj1856_thin_atm.tex:123`, anchor `$k\tinfty_{\rm opt}=26$~eV`.

Extinction remains covariant with optical slope. Kaplan explicitly says reddening is

> `highly degenerate with the power-law index`.

Source: `refs/kaplan_2011_rxj1856_spectrum.tex:596`.

For RX J1856, its table gives $d\alpha/dA_V=1.9$ and $\alpha=3.94\pm0.06$. A chosen extinction perturbation $\Delta A_V=0.025$ therefore changes the inferred slope by **0.0475**, comparable to the printed statistical error.

Source: `refs/kaplan_2011_rxj1856_spectrum.tex:1035`, anchor `$3.94\pm0.06$& 1.9`.

## 4.3 A pulse amplitude supplies one constraint, not a surface map

**DERIVED HERE:** for an axisymmetric pattern,

$$
\cos\psi(\phi)
=
\cos i\cos\beta+\sin i\sin\beta\cos\phi.
$$

This is the dot product of the line of sight and the rotating symmetry axis.

Let a visible patch contribute $K\,\mathcal G_u(\cos\psi)$ on top of a constant contribution $C_0$, where $\mathcal G_u$ includes exact bending and the emission angular law. Expanding around $C=\cos i\cos\beta$, with $S=\sin i\sin\beta$,

$$
C(\phi)\simeq
C_0+K\mathcal G_u(C)
+
KS\mathcal G_u'(C)\cos\phi,
$$

so

$$
a_1\simeq
\frac{|KS\mathcal G_u'(C)|}
{C_0+K\mathcal G_u(C)}.
$$

A small amplitude can follow from small contrast, unfavorable orientation, weak angular variation of the ray-integrated emission, or combinations of these. For this symmetric model, exchanging $i$ and $\beta$ leaves the light curve unchanged.

A uniform isotropic blackbody predicts $a_1=0$. Under the requested Gaussian summary constraint, Python gives

$$
\chi_a^2(0)=256.
$$

Thus that model fails the pulse constraint as well as the broadband spectral test. It does **not** follow that a particular temperature gradient or compactness is measured.

## 4.4 Redshift and compactness are the same degree of freedom

There is no independent “$z_g$ versus $M/R$” degeneracy: in the adopted exterior metric, $z_g$ is a one-to-one function of $M/R$.

The actual degeneracy is that a uniform blackbody measures $T_\infty$ and $R_\infty/d$, with **no separate spectral sensitivity to $z_g$**. Changing $z_g$ while adjusting local temperature, radius, and mass according to Section 2 leaves the entire distant blackbody spectrum unchanged.

Atmospheres may provide additional sensitivity through gravity, energy-dependent opacity, and spectral features. But an observed feature obeys

$$
E_{\rm obs}=\frac{E_0(B,\mathrm{composition},\ldots)}{1+z_g},
$$

so an uncertain intrinsic feature energy also couples $B$ and redshift.

Ho’s own two trial models illustrate conditional sensitivity: changing the prescribed field from $3$ to $4\times10^{12}\,\mathrm G$ changes the reported redshift from $0.27$ to $0.22$, each with formal error $0.02$. The shift is **2.5 times** either formal error.

Source: `refs/ho_2007_rxj1856_thin_atm.tex:628` and `:641`, anchor `$\zg$ & 0.27(2) & 0.22(2)`.

Using the Ho normalization only as an illustration,

$$
R_\infty=17\frac{123}{140}=14.935714\ \mathrm{km}.
$$

At this fixed apparent radius, Python gives:

| Assumed $z_g$ | $u$ | $R$ | Mass in nominal solar units |
|---:|---:|---:|---:|
| 0.10 | 0.173554 | 13.5779 km | 0.797934 |
| 0.22 | 0.328138 | 12.2424 km | 1.360260 |
| 0.40 | 0.489796 | 10.6684 km | 1.769347 |

For a uniform blackbody with fixed $T_\infty$, these choices produce the same distant spectrum. They are a degeneracy demonstration, not astrophysical mass estimates.

The nominal mass conversion uses the locally printed solar mass parameter: `refs/prsa_2016_iau_nominal.tex:135`, anchor `$1.327\,124\,4 \times 10^{20}$`.

## 4.5 Fisher structure

Taking the expected Hessian of the likelihood gives **DERIVED HERE**, for fixed nuisance parameters,

$$
\mathcal F_{ab}
=
\sum_{ri}\frac{\partial_a\mu_{ri}\,\partial_b\mu_{ri}}{\mu_{ri}}
+
(\partial_a\boldsymbol m)^\mathsf T
\Sigma_H^{-1}
(\partial_b\boldsymbol m)
+
\frac{\partial_a a_1\,\partial_b a_1}{\sigma_a^2}
+
\mathcal F^{\rm astrometry}_{ab}.
$$

The structure is more informative than a guessed numerical matrix:

- X-ray derivatives couple temperature, angular normalization, absorption, and spectral calibration.
- Optical normalization contributes a direction approximately proportional to $(2,1)$ in $(\ln\rho_\infty,\ln T_\infty)$.
- Optical colors add information on spectral departure from Rayleigh–Jeans behavior and extinction.
- A single pulse amplitude adds an outer product of **rank at most one**.
- The uniform-blackbody redshift derivative is exactly zero at fixed $T_\infty,\rho_\infty$.
- Astrometry determines the conversion from angular to physical scale.

For identifiable nuisance parameters $\eta$, profiling gives the Schur complement,

$$
\mathcal F_{\rm eff}
=
\mathcal F_{\theta\theta}
-
\mathcal F_{\theta\eta}
\mathcal F_{\eta\eta}^{-1}
\mathcal F_{\eta\theta}.
$$

Freezing calibration, absorption composition, or extinction therefore overstates information. A singular block should be reported as non-identifiability, not concealed by an unexplained numerical regularizer.

# 5. Explicit numerical optical-excess calculation

## 5.1 Constants and distance rescaling

The Python calculation used the locally printed constants:

$$
h=6.62607015\times10^{-34}\ \mathrm{J\,s},
\quad
c=299792458\ \mathrm{m\,s^{-1}},
$$

$$
k_{\rm B}=1.380649\times10^{-23}\ \mathrm{J\,K^{-1}},
\quad
1\,\mathrm{eV}=1.602176634\times10^{-19}\ \mathrm J.
$$

Source: `refs/tiesinga_2021_codata2018.txt:3403`, `:3405`, `:3420`, and `:3426`; short anchors are `299 792 458`, `6.626 070 15`, `1.380 649`, and `electron volt`.

The parsec follows the local exact definition:

> `the parsec to be exactly $648000/\pi$\,au`

with the astronomical unit specified in the same paragraph.

Source: `refs/prsa_2016_iau_nominal.tex:172`.

Computed values are

$$
1\,\mathrm{pc}=3.085677581491367\times10^{16}\ \mathrm m,
$$

$$
T_\infty=736886.9007\ \mathrm K,
\qquad
R_{\rm bb,\infty}=4.510000\ \mathrm{km},
$$

$$
\frac{R_{\rm bb,\infty}}d
=
1.18828573946\times10^{-15}.
$$

The important cancellation is

$$
\frac{4.4\,\mathrm{km}\,(d/120\,\mathrm{pc})}{d}
=
\frac{4.4\,\mathrm{km}}{120\,\mathrm{pc}}.
$$

Changing the adopted distance changes the inferred emitting radius, not the predicted Earth flux of the already fitted blackbody.

## 5.2 Comparison at $2500\,\mathrm{\AA}$

Using

$$
F_\lambda^0
=
\pi\left(\frac{R_{\rm bb,\infty}}d\right)^2B_\lambda(T_\infty)
$$

gives

$$
x_{2500}=0.0781002825,
$$

$$
F_{2500}^0
=
6.66034604\times10^{-19}
\ \mathrm{erg\,cm^{-2}\,s^{-1}\,\AA^{-1}}.
$$

The pure Rayleigh–Jeans approximation would overestimate this by **4.00869%**. That is too large to discard in a precision photometric calculation.

For this verification calculation only, use Kaplan’s tabulated $A_V=0.05$ and

$$
A_{2500}/A_V=2.32.
$$

Source: `refs/kaplan_2011_rxj1856_spectrum.tex:1035` and `:1046`, anchor `2500\,\AA\ is taken as $A_{2500\,{\rm \AA}}/A_V=2.32$.`

These are explicit extinction-model inputs, not an owned extinction calculation.

Then

$$
F_{2500}^{\rm BB,obs}
=
6.66034604\times10^{-19}\,10^{-0.4(2.32)(0.05)}
=
5.98544982\times10^{-19}.
$$

Kaplan’s observationally fitted, **absorbed** flux is

$$
F_{2500}^{\rm data\ fit}
=
(4.25\pm0.07)\times10^{-18}.
$$

Source: `refs/kaplan_2011_rxj1856_spectrum.tex:1035`, anchor `$4.25\pm0.07$`; the table note at line 1039 explicitly says `the absorbed flux density at 2500\,\AA`.

Converting both with the STMAG definition:

| Quantity | Python result |
|---|---:|
| Blackbody monochromatic STMAG | 24.457258 |
| STMAG equivalent of Kaplan’s fitted flux | $22.329028\pm0.017883$ |
| Magnitude deficit | 2.128230 |
| Data/model flux ratio | 7.100552 |
| Published excess | $7.1\pm0.3$ |
| Absolute difference from published central value | 0.000552 |
| Comparison tolerance | Published uncertainty, 0.3 |

The small central-value difference verifies the arithmetic and conventions. It is not an independent observational confirmation: the comparison intentionally uses the same published blackbody parameters and extinction prescription.

Even setting extinction to zero leaves a data/model ratio of **6.38105**. Positive extinction can only reduce the blackbody’s observed flux further. Extinction uncertainty cannot repair the single-temperature discrepancy.

## 5.3 F475W: useful numerical check, missing direct measurement

Kaplan Table 4 gives the **predicted** ACS/WFC F475W blackbody magnitude

$$
m_{\rm F475W}^{\rm BB}=27.13\pm0.05.
$$

Source: `refs/kaplan_2011_rxj1856_spectrum.tex:730`, table title anchor `Predicted \hst\ Optical Photometry from Extrapolated X-ray Blackbodies`; numerical anchor `$27.13\pm0.05$` at line 746.

Using its approximate effective wavelength $4709\,\mathrm{\AA}$ and extinction coefficient $A_b/A_V=1.24$ gives:

$$
F_{4709}^0=5.38985333\times10^{-20},
$$

$$
m_{\rm proxy}=27.133058.
$$

The difference from the published band-integrated value is **0.003058 mag**, below the published $0.05$ mag model uncertainty.

The effective wavelength and coefficient are printed at `refs/kaplan_2011_rxj1856_spectrum.tex:751`, anchors `4520 & 4709 & 5066` and `1.93 & 1.30 & 1.24`.

**This is a monochromatic proxy, not our own F475W integration.** Its close agreement does not validate a missing throughput file.

Kaplan’s $4700\,\mathrm{\AA}$ excess is $7.3\pm0.3$, from the fitted optical SED.

Source: `refs/kaplan_2011_rxj1856_spectrum.tex:1035`, anchor `$7.3\pm0.3$`.

Combining that ratio with the published F475W prediction gives an approximate implied magnitude,

$$
27.13-2.5\log_{10}(7.3)=24.971693.
$$

That number is **not a locally tabulated measured F475W STMAG**, and must not be entered as one.

The available Walter measurement is

> `The F475W magnitude of the neutron star is 25.36$\pm$0.15.`

Source: `refs/walter_2010_rxj1856_distance.tex:862`.

Its accompanying text does not establish the magnitude system needed for this comparison. It also concerns HRC photometry, whereas Kaplan’s quoted prediction is for ACS/WFC.

Therefore the requested direct measured-band comparison cannot honestly be completed from the present local files. The defensible completed comparison is the $2500\,\mathrm{\AA}$ flux-summary calculation above. The missing original photometry and throughput curves are a real Milestone 2 blocker, not permission to relabel a prediction as data.

# 6. Circularity, legitimate inputs, and the minimal honest Milestone 2

## 6.1 What is circular—and what is merely conditional

Let $D$ denote the observations and $\widehat\theta_{\rm Ho}(D)$ the published inference.

Rendering with those parameters and comparing with overlapping observations,

$$
D\longrightarrow\widehat\theta_{\rm Ho}(D)
\longrightarrow F(\widehat\theta_{\rm Ho})
\longrightarrow D,
$$

does not independently test whether the project has inferred RX J1856’s properties. It can test implementation compatibility with the published model.

However, the precise historical claim matters: Ho’s Table 1 is primarily a LETG fit with manually explored model choices, followed by an optical comparison. Calling it a fully optimized joint HST/LETG fit would be inaccurate. Its optical comparison retains some extrapolative content, while model selection and extinction choices have already used optical information.

Replacing the optimizer does not remove this logical issue:

$$
D\longrightarrow\widehat\theta_{\rm ours}(D)
\longrightarrow F(\widehat\theta_{\rm ours})
\longrightarrow D
$$

is an owned fit and a legitimate goodness-of-fit calculation. It is still not an independent prediction.

**Minimal honest reporting:** label the exact observations used for fitting, model selection, calibration decisions, and final validation. A held-out optical prediction is stronger than an optical residual from a joint fit, but both can be useful when named correctly.

## 6.2 Legitimate inputs

The following are legitimate inputs within their stated provenance:

- Independent astrometry, with its systematic treatment.
- Extracted observational counts, background counts, exposure metadata, calibrated ARFs/RMFs, and HST response curves.
- Measured photometry in a confirmed system, with calibration and aperture conventions.
- Timing measurements, with the correct estimator and instrument selection.
- Physical constants and explicitly adopted conventions.
- Independently measured sightline constraints, if actually obtained and documented.

The following are not independent measurements of the required physical parameters:

- Burwitz or Ho $N_{\rm H}$ as a prior on the same LETG data.
- Ho $T_\infty,R_\infty,z_g$ as fixed runtime truth.
- Kaplan’s predicted magnitudes.
- A fitted optical slope or excess treated as independent of the photometry used to derive it.
- A spin-down field treated as a direct measurement of the surface field.

The local review explicitly says the spin-down field depends on the “canonical model” and changes with assumed radius and moment of inertia.

Source: `refs/potekhin_ho_chabrier_2016_review.tex:1206`, anchors `the ``canonical model''` and `$\sqrt{I}/R^3$`.

## 6.3 Absorption and extinction ownership cannot be solved by renaming inputs

Wilms decomposes absorption into gas, molecules, and grains:

$$
\sigma_{\rm ISM}
=
\sigma_{\rm gas}+\sigma_{\rm molecules}+\sigma_{\rm grains}.
$$

Source: `refs/wilms_allen_mccray_2000_ism.tex:120`, anchor `molecules} + \sigma_{\rm grains}`.

Its gas contribution includes elemental abundances, ion fractions, depletion, and atomic cross-sections. Recomputing one atomic cross-section is therefore not equivalent to owning the entire attenuation model.

Moreover, the audit’s “Wilms plus Verner 1996” shorthand is incomplete for exact reproduction. The local Wilms text specifies helium from Yan et al. and metal cross-sections from Verner’s earlier work.

Source: `refs/wilms_allen_mccray_2000_ism.tex:223`, anchor `The cross section for helium is taken`; and `:241`, anchor `For all other elements we use the photoionization cross sections`.

Verner 1996 itself describes its inputs as theoretical cross-sections “interpolated and smoothed over resonances.”

Source: `refs/verner_1996_photoionization.tex:35`, anchor `on the Opacity Project theoretical cross sections interpolated and smoothed`.

Implementing its formula locally owns evaluation of the fit, not the underlying atomic calculation.

Likewise, a free $A_V$ does not make an imported CCM curve an owned physical calculation. CCM explicitly acknowledges:

> `real deviations of individual lines of sight from the mean law`

Source: `refs/cardelli_1989_extinction.txt:22`.

Computing grain extinction from first principles would still require a grain population, composition, and sightline normalization. Those quantities do not become known merely because the scattering solver is owned.

## 6.4 Minimal honest Milestone 2

Milestone 2 should do the following:

1. **Acquire and checksum actual measurements and responses.** Do not fit plots or published unfolded spectra as substitutes for counts.
2. **Fit an owned uniform-blackbody forward model to LETG counts.** Profile or marginalize background and declared calibration nuisances.
3. **Freeze that X-ray fit before optical prediction.** Render its Earth-distance spectrum and integrate the actual HST responses.
4. **Report the optical failure and pulse failure numerically.** A failure is the intended scientific result of this baseline.
5. **Fit a more flexible phenomenological model separately.** A two-temperature blackbody map can diagnose what additional area–temperature moment is needed; it cannot validate atmosphere microphysics.
6. **Carry unresolved compactness and geometry as an ensemble.** Do not manufacture a unique $(M,R)$ from the hot blackbody normalization.
7. **Retain yellow ownership status wherever absorption, extinction, or surface physics remains a stub.**

A labeled Wilms/Verner absorption implementation and CCM extinction implementation are defensible temporary dependencies for this diagnostic milestone. Their uncertainty must be explored by refitting under alternative assumptions; they cannot be assigned an invented universal percentage.

For illustration, Kaplan’s $\pm50\%$ perturbation around $A_V=0.05$ gives, using its tabulated coefficients:

| Quantity | Computed change |
|---|---:|
| F475W extinction | $\pm0.031$ mag |
| $2500\,\mathrm{\AA}$ extinction | $\pm0.058$ mag |
| Corresponding positive flux factors | 1.028964 and 1.054873 |

These are sensitivity calculations, not measured error bars. They are small compared with the blackbody deficit but material for a precise R5 claim.

*Execution note: the numerical calculations above ran in Python, and the quoted principal anchors were checked against local text. No archive-level fit was possible with the available inputs. The workspace is read-only, so this session could not update HANDOFF.md or the worklog.*

# 7. Concrete implementation plan

The optimizer should be the final composition layer. Neither fitted parameters nor photometric acceptance results should feed back into the forward-model dependency graph.

| Function | Inputs | Outputs and obligations |
|---|---|---|
| `load_observation_manifest` | Checksummed event/extraction products, responses, photometry, calibration versions | Source/background counts, scaling factors, energy grids, band definitions, covariance, fit/validation membership |
| `geometry_from_angular_radius` | $\rho_\infty,d,z_g$ | Consistent $R,GM,g_s$; reject unsupported compactness or inconsistent parameters |
| `surface_stokes` | Model family, patch conditions, energy, direction | Absolute spectral Stokes intensity; no hidden color correction or out-of-grid extrapolation |
| `earth_spectrum` | Geometry, surface map, phases | Intrinsic phase-resolved and averaged flux from the same ray/surface interface used by the renderer |
| `fold_letg` | Intrinsic flux, absorption model, ARF/RMF for each order, exposure | Expected source counts per observation/channel |
| `hst_stmag` | Intrinsic spectrum, extinction, exact throughput | Count rate, band-equivalent flux, STMAG |
| `pulse_observable` | Phase-resolved spectra, timing responses and selections | Predicted estimator matching the measured sinusoidal amplitude, plus harmonic diagnostics |
| `joint_objective` | Counts, magnitudes, timing, independent constraints | Separate likelihood terms, nuisance penalties, explicit invalid-domain diagnostics |
| `fit_and_diagnose` | Objective, model family, declared priors and starts | Conditional optimum, uncertainty samples or profiles, rank/degeneracy diagnostics, model-comparison record |
| `validate_frozen_fit` | Frozen fit, held-out observations, renderer output | Residuals with expected value, computed value, uncertainty/tolerance, and provenance |

## 7.1 Red–green verification targets

Write the failing tests before production implementation.

| Test | Expected value or limit | Tolerance/status |
|---|---|---|
| Uniform-disc flux | $\pi(R_\infty/d)^2B_\lambda(T_\infty)$ | Proposed numerical relative tolerance $10^{-8}$ |
| Zero-count Poisson term | $C(0\mid2)=4$ | Exact algebraic target; Python computed 4 |
| Saturated Poisson term | $C(3\mid3)=0$ | Exact algebraic target; Python computed 0 |
| Distance invariance | Same Earth flux when $R_\infty$ and $d$ scale together | Numerical precision |
| Blackbody redshift degeneracy | Same spectrum at fixed $T_\infty,R_\infty/d$ | Numerical precision |
| Isothermal absorbing slab | Analytic solution in Section 2.3 | Grid-convergence test |
| Optical benchmark | Excess $7.1\pm0.3$ at $2500\,\mathrm{\AA}$ | Current computation 7.100552; reproduction benchmark only |
| Full F475W integration | Kaplan prediction $27.13\pm0.05$ under matched assumptions | Pending throughput; proxy calculation is not this test |
| LETG blackbody recovery | $63.5\pm0.2\,\mathrm{eV}$; $(4.4\pm0.1)\,\mathrm{km}(d/120\,\mathrm{pc})$; $N_{\rm H}=(0.95\pm0.03)\times10^{20}$ | Pending counts/responses; compare under documented calibration and absorption assumptions |
| Pulse estimator | $0.0096\pm0.0006$ | Proposed numerical amplitude accuracy $0.00006$ before observational comparison |

The proposed photometric integration tolerance is **0.001 mag**, equivalent to a fractional flux change of **0.00092146**, computed in Python. It is a numerical engineering target, not an observational uncertainty.

For LETG recovery, do not force current calibrations and a Poisson likelihood to reproduce an older $\chi^2$ fit exactly. Print the displacement, then attribute it through controlled changes in calibration, extraction, absorption, and statistic.

Burwitz’s archive identifiers are **113, 3382, 3380, 3381, 3399**. The listed exposures sum in Python to **505.39 ks**.

Source: `refs/burwitz_2003_rxj1856_xray.tex:152`, anchors `113` and `LETGS`; total exposure appears at line 159.

Mutation tests should target omitted $\pi$, incorrect wavelength units, reversed redshift, missing photon weighting, incorrect background scaling, and dropped diffraction orders. Fuzzing should cover response-grid boundaries, empty bins, narrow spectral structure, and exact zero-count cases. Invalid model domains must produce explicit diagnostics rather than clamped spectra.

# 8. Missing sources and data to fetch

All entries below are **UNVERIFIED as newly acquired full sources or data products**. Bibliographic information is locally supported where stated; no unverified identifier should be treated as resolved.

| Priority | Missing item | Citation or retrieval lead |
|---|---|---|
| Blocking | Original RX J1856 optical/UV measurements, systems, aperture corrections, and uncertainties | **van Kerkwijk & Kulkarni 2001, A&A 378, 986**; arXiv ID not established here. Kaplan’s local bibliography explicitly gives `2001, \aap, 378, 986` at `refs/kaplan_2011_rxj1856_spectrum.tex:1371`. **GAPS §2.1 cites the other paper, A&A 380, 221; that should be corrected.** |
| Blocking | LETG source/background products and order-dependent responses | Chandra X-ray Center archive; Burwitz ObsIDs listed above. Fetch reduction/calibration documentation and pin versions. No arXiv identifier applies to the data products. |
| Blocking | HST throughput and calibration files | STScI calibration archive for the actual WFPC2/ACS/STIS measurements used. Record instrument, detector, epoch, aperture convention, and checksum. |
| Blocking for exact timing likelihood | Timing responses, event selections, and exposure weights | Original products corresponding to van Kerkwijk & Kaplan 2008 Table 1; a universal bolometric amplitude cannot replace them. |
| Statistical documentation | Cash 1979, ApJ 228, 939 | Full paper absent; local bibliographic anchor at `refs/vankerkwijk_kaplan_2008_rxj1856_timing.tex:518`. No arXiv ID established. |
| Ho reproduction | O’Donnell 1994, ApJ 422, 158 | Ho uses its extinction corrections; citation printed at `refs/ho_2007_rxj1856_thin_atm.tex:1152`. No arXiv ID established. |
| Absorption reproduction | Verner & Yakovlev 1995, A&AS 109, 125; Yan, Sadeghpour & Dalgarno 1998, ApJ 496, 1044 | Local Wilms bibliography at lines 1126 and 1154. These clarify dependencies not supplied by Verner 1996 alone. IDs not established. |
| Atomic ownership | Seaton et al. 1992, Rev. Mex. Astron. Astrofis. 23, 19, and the relevant atomic-method papers | Verner’s local bibliography gives this reference at `refs/verner_1996_photoionization.tex:382`. A review or table alone will not supply every implementable continuum calculation. |
| Geometry comparison | Ho 2007, MNRAS 380, 71, *Constraining the geometry of the neutron star RX J1856.5−3754* | arXiv:**0705.4543**, printed at `refs/potekhin_ho_chabrier_2016_review.tex:2144`; full paper unverified locally. Distinct from the thin-atmosphere paper. |
| Extinction ownership, if pursued | Weingartner & Draine 2001, ApJ 548, 296; Draine 2003, ARA&A 41, 241; Bohren & Huffman 1983, book | Candidate arXiv IDs **astro-ph/0008146** and **astro-ph/0304489**, respectively, **unverified** leads from the project gap list. Also requires observational constraints on the adopted grain population. |

Wilms 2000, Verner 1996, Predehl–Schmitt 1995, and CCM 1989 are already local. Their remaining problem is implementation and provenance, not acquisition.

# 9. Effect of stubbing A34 on the final render and R5

**Using the X-ray single-temperature blackbody:** the optical source is too faint by approximately the locally reported factor $7.3$. This corresponds to **2.15831 mag**, **86.3014% missing optical flux**, or **2.86790 exposure doublings**. These are Python conversions of Kaplan’s excess, not new measurements. Exposure or tone mapping can conceal the display discrepancy but cannot repair absolute radiance or R5.

**Using Ho’s parameters as a labeled preview:** the result is a conditional published-model reconstruction. Ho’s reported optical shortfall corresponds to approximately **0.1517–0.1980 mag** when the data are $15–20\%$ above the model. The F606W discrepancy is **five times** the quoted photometric error before the additional model and calibration terms are included. It therefore cannot be advertised as demonstrated agreement within the photometric errors alone.

**Failing to rescale Ho’s radius:** retaining $17\,\mathrm{km}$ at the adopted distance instead of the computed $14.9357\,\mathrm{km}$ increases flux by **29.5525%**, a magnitude shift of **−0.281115 mag**. Rescaling fixes that normalization error; a physically consistent atmosphere refit must also update surface gravity.

**Using an owned but flexible joint fit:** R5 may pass, including with a multitemperature blackbody. That establishes observational consistency at the fitted bands. It does not establish the atmosphere composition, condensed-surface state, or uniqueness of the surface map.

**For the resolved flyby:** no defensible universal percentage error follows from the current data. The explicit redshift family above changes radius and mass while preserving the distant blackbody spectrum. The pulse-amplitude expansion admits different contrasts and orientations with the same measured amplitude. Those alternatives affect lensing, angular emission, apparent surface structure, and the camera’s encounter with the radiation field. The uncertainty must therefore be rendered as a conditional ensemble.

Milestone 2 can honestly finish with an owned blackbody fitting and photometry harness, a quantified broadband failure, and clearly labeled absorption/extinction dependencies. A34 should remain conditional—and its recomputation branch yellow—until the required measurements are present and the remaining physical dependencies are owned.