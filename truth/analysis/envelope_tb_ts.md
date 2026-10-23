**Abstract.** The proposed replacement of `pipeline.envelope-tb-ts` is incomplete. A stationary one-dimensional envelope solve can replace the local vertical transport calculation, but it cannot reproduce the lateral heat transport already included in PPP15’s equatorial temperature fit. For the requested benchmark, interpreting $10^{13}\,\mathrm G$ as the **polar** field, direct evaluation gives $T_{\rm p}=1.244174\times10^6\,\mathrm K$, $T_{\rm eq}=2.831225\times10^5\,\mathrm K$, and $T_{\rm p}/T_{\rm eq}=4.394473$. These are fit evaluations, not predictions from an owned transport solver. The minimum honest replacement is a validated column solver followed by an axisymmetric envelope calculation, with explicit composition and microphysics assumptions. PPP15 is defensible as a labelled interim model; it cannot establish first-principles ownership or, by itself, compliance with R5.

# 1. What problem must actually be solved?

The database correctly classifies the relevant `surface.ppp15-*` entries as `provenance = "fit"` and `ownership = "stub"`. The repeated reimplementation route—“own 1-D magnetised envelope solve”—does not cover everything these entries represent.

PPP15 explicitly describes the underlying calculation as two-dimensional:

> `refs/potekhin_pons_page_2015_thermal.tex:3979–3987`: “with the 2D code”; “the non-radial heat flow increases the equatorial”.

Consequently, distinguish three objects:

1. **Local transport tensor:** $\boldsymbol\kappa(\rho,T,\boldsymbol B,\mathrm{composition})$.
2. **Independent-column relation:** $T_s(T_b,B,\Theta_B,g;\mathrm{composition})$, conditional on negligible lateral transport.
3. **Global envelope map:** $T_s(\theta,\phi)$, obtained from a spatial boundary-value problem with a specified magnetic geometry and bottom-temperature distribution.

The third object generally cannot be represented by a function of the four local arguments in the topic title. Two surface patches with identical local $T_b,B,\Theta_B,g$ can receive different lateral fluxes if their neighbouring field geometry differs. This is a consequence of the differential equations below, **DERIVED HERE**, rather than an additional empirical claim.

Throughout, $T_s$ means the local effective temperature, $\Theta_B$ the field–normal angle, and $\theta$ magnetic colatitude. They are different angles. Conductivity is denoted $\kappa$ and mass opacity $K$.

All numerical evaluations below were run in Python. Fit arithmetic was cross-checked using ordinary floating point and 60-digit arithmetic. No production envelope solver was implemented. The workspace is read-only, so the required handoff and worklog updates could not be written.

## 1.1 Tensor projection: exact geometry inside an approximate model

PY01 prints the required projection directly:

> `refs/potekhin_yakovlev_2001_cooling.tex:441–453`: “off-diagonal (Hall) component”; `\kappa=\kappa_\|\cos^2\theta + \kappa_\perp\sin^2\theta`.

**DERIVED HERE.** Let $\boldsymbol b=\boldsymbol B/B$, and define the conductivity tensor by

$$
\boldsymbol\kappa
=
\kappa_\perp\boldsymbol I
+(\kappa_\parallel-\kappa_\perp)\boldsymbol b\boldsymbol b
+\kappa_H\boldsymbol A,
\qquad
\boldsymbol A\boldsymbol v=\boldsymbol b\times\boldsymbol v.
$$

This states the adopted Hall-sign convention. Let $z$ increase outward along the surface normal $\boldsymbol n$. Under the independent-column assumption,

$$
\nabla T=\boldsymbol n\,\frac{dT}{dz},
\qquad
\boldsymbol F=-\boldsymbol\kappa\nabla T.
$$

Contracting with $\boldsymbol n$ gives

$$
F_n
=
-\left[
\kappa_\perp+
(\kappa_\parallel-\kappa_\perp)(\boldsymbol n\cdot\boldsymbol b)^2
\right]\frac{dT}{dz},
$$

because $\boldsymbol n\cdot(\boldsymbol b\times\boldsymbol n)=0$. Therefore

$$
\boxed{
\kappa_{\rm eff}
=
\kappa_\parallel\cos^2\Theta_B
+\kappa_\perp\sin^2\Theta_B,
\qquad
F_n=-\kappa_{\rm eff}\frac{dT}{dz}.
}
$$

The Hall component produces no normal flux for a purely normal temperature gradient. It also produces no local dissipation, since $\nabla T\cdot(\boldsymbol b\times\nabla T)=0$.

This projection is not the questionable part of GH83. The questionable step is prescribing the temperature field so that lateral derivatives can be neglected, then replacing the depth-dependent transport problem by a simple surface interpolation.

GH83 itself acknowledges the simplification:

> `refs/greenstein_hartke_1983_surface_T.txt:153–156`: “is beyond the scope of this paper”; “we also ignore the dependence of conductivity on temperature”.

GH83’s angle $\Phi$ is measured from the horizontal, so $\sin^2\Phi=\cos^2\Theta_B$. Confusing $\Phi$, $\Theta_B$, and colatitude changes the temperature map even if every conductivity is correct.

## 1.2 Hydrostatic structure and the temperature equation

PY01 supplies the thin-envelope hydrostatic reduction:

> `refs/potekhin_yakovlev_2001_cooling.tex:427–435`: `\mathrm{d} P/\mathrm{d} z=g\rho`; “surface gravity”.

Its $z$ increases inward. With the outward coordinate adopted here,

$$
\frac{dP}{dz}=-\rho g,
\qquad
g=\frac{GM}{R^2\sqrt{1-2GM/(Rc^2)}}.
$$

**DERIVED HERE**, eliminating $z$ between hydrostatics and conduction:

$$
\frac{dT}{dP}
=
\frac{F_n}{g\rho\kappa_{\rm eff}}.
$$

Define an equivalent transport opacity,

$$
K_{\rm eff}
=
\frac{16\sigma_{\rm SB}T^3}{3\rho\kappa_{\rm eff}},
\qquad
F_n=\sigma_{\rm SB}T_s^4.
$$

Then

$$
\boxed{
\frac{d\ln T}{d\ln P}
=
\frac{3PK_{\rm eff}}{16g}\frac{T_s^4}{T^4}.
}
$$

This reproduces PPP15 Eq. 38:

> `refs/potekhin_pons_page_2015_thermal.tex:2990–2995`: `\frac{3}{ 16}\,\frac{P\opac}{ g}\,\frac{\Ts^4}{ T^4}`.

The conductivity–opacity conversion is separately printed at:

> `refs/potekhin_pons_page_2015_thermal.tex:2589–2594`: `\frac{16\sSB T^3}{3\rho\opac}`.

The conversion must use the **total normal conductivity**. Electron and photon conductivities add:

$$
\kappa_{{\rm eff},{\rm total}}
=
\kappa_{{\rm eff},e}+\kappa_{{\rm eff},r},
\qquad
K_{\rm total}^{-1}=K_e^{-1}+K_r^{-1}.
$$

Source:

> `refs/potekhin_yakovlev_2001_cooling.tex:614–619`: “the two mechanisms work in parallel”; `K^{-1} = {K}_\mathrm{r}^{-1} + {K}_\mathrm{c}^{-1}`.

Adding the equivalent opacities would describe the wrong transport problem.

## 1.3 What the EOS and boundaries must provide

At each integration point the solver needs $\rho(P,T,B,\{X_i\})$, electron chemical potential, ionization state, phase, and thermodynamic derivatives.

PC13 decomposes the free energy as

$$
F=F_{\rm id}^{(i)}+F_{\rm id}^{(e)}+F_{ee}+F_{ii}+F_{ie},
\qquad
P=-\left(\frac{\partial F}{\partial V}\right)_{T,B}.
$$

Source:

> `refs/potekhin_chabrier_2013_mag_eos.tex:417–426`: `+ F_\mathrm{ee} + F_\mathrm{ii} + F_\mathrm{ie}`;  
> `refs/potekhin_chabrier_2013_mag_eos.tex:451–458`: `P=-(\partial F/\partial V)_T`.

For the prescribed-field hydrostatic model, use the thermodynamic pressure consistently, rather than substituting the transverse kinetic electron pressure:

> `refs/potekhin_chabrier_2013_mag_eos.tex:1520–1522`: “thermodynamic pressure $P$, which automatically includes”.

The ideal electron contribution can be recomputed from Landau-state sums and momentum integrals, solving charge neutrality for the chemical potential. PC13 provides both number density and grand potential:

> `refs/potekhin_chabrier_2013_mag_eos.tex:1052–1068`: `\sum_{n=0}^{\infty} \sum_\sigma`; “factor 2 for $n\geq1$”.

That does **not** make the nonideal terms owned. The decomposition identifies additional calculations.

For the boundary-value problem:

- **Bottom:** impose $T=T_b$ at a stated $\rho_b$, with composition and field specified there.
- **Top:** match to an atmosphere or condensed-surface calculation that returns the relation between outgoing flux and the thermodynamic state at an overlap depth.
- **Unknown:** normally the emergent flux, found by shooting or a coupled boundary solve.

PPP15 uses $\rho_b=10^{10}\,\mathrm{g\,cm^{-3}}$ for the Appendix-B fit:

> `refs/potekhin_pons_page_2015_thermal.tex:3975–3977`: `chosen $\rhob=10^{10}$`.

The condition $F=\sigma_{\rm SB}T_s^4$ defines effective temperature; it does not require a blackbody spectrum. Nor does it justify setting the material temperature equal to $T_s$ at an arbitrary atmosphere–envelope interface. PY01 explicitly distinguishes the Eddington boundary approximation from a radiative-transfer boundary:

> `refs/potekhin_yakovlev_2001_cooling.tex:418–425`: “A more accurate boundary condition”; “requires solution of the radiative transfer equation”.

**DERIVED HERE**, stationary energy conservation is

$$
\nabla\cdot\boldsymbol F=Q_{\rm heat}-Q_\nu.
$$

Constant flux follows only in a source-free one-dimensional column. In two dimensions, even without neutrino losses, the bottom and top fluxes of an individual column need not agree: lateral transport can enter or leave it. The surface-integrated energy balance must agree.

# 2. Nonmagnetic scaling: what is derivable, and what is fitted?

## 2.1 Rosseland transport and the sensitivity strip

PPP15 prints the Rosseland weighting:

> `refs/potekhin_pons_page_2015_thermal.tex:2634–2643`: “Rosseland mean”; `\frac{z^4\mathrm{e}^z}{(\mathrm{e}^z-1)^2}`.

Restoring kelvin units and using $x=h\nu/(k_BT)$,

$$
\frac{1}{K_R}
=
\int_0^\infty \frac{w(x)}{K_\nu}\,dx,
\qquad
w(x)=\frac{15}{4\pi^4}
\frac{x^4e^x}{(e^x-1)^2}.
$$

Together with the conductivity–opacity relation above,

$$
\kappa_r=\frac{16\sigma_{\rm SB}T^3}{3\rho K_R}.
$$

The thermal resistance per logarithmic pressure interval is, **DERIVED HERE**,

$$
\frac{dT}{d\ln P}
=
\frac{F}{g}\frac{P}{\rho(\kappa_r+\kappa_e)}.
$$

Layers in which this expression is large make a large contribution to the temperature drop. Neither a highly transparent radiative exterior nor a highly conducting deep interior necessarily controls that drop. The transition between them is the natural candidate.

PPP15 identifies it explicitly:

> `refs/potekhin_pons_page_2015_thermal.tex:2941–2946`: “sensitivity strip”; `$\kappa_\mathrm{r}=\kappa_e$`; `$\rho\sim10^5$\,--\,$10^7$`.

The equality is a locator for the important region, not a replacement for integrating through it.

## 2.2 Why the gravity exponent is robust

Write the structure equation as

$$
\frac{dT}{dP}
=
\frac{\sigma_{\rm SB}}{\rho\kappa}\frac{T_s^4}{g}.
$$

**DERIVED HERE.** If microphysics depends on $P,T$ and composition, with no independent $g$ dependence, then two envelopes with the same $T_s^4/g$ follow the same $T(P)$ differential equation. If the surface boundary contributes negligibly to $T_b$, their bottom temperatures are nearly identical. Hence

$$
\boxed{T_s^4\simeq g\,\mathcal F(T_b,\mathrm{composition}).}
$$

The qualification concerning the outer boundary is essential. PY01 makes precisely this argument:

> `refs/potekhin_yakovlev_2001_cooling.tex:1173–1181`: `scaling law $T_\mathrm{s} \propto g^{1/4}$`; “almost independent of the outer boundary”.

## 2.3 Why the temperature exponent is not universal

PY01 gives the free-free Rosseland scaling

$$
K_{\rm ff}\propto \rho T^{-7/2},
$$

when composition and the slowly varying Gaunt-factor contribution are held fixed:

> `refs/potekhin_yakovlev_2001_cooling.tex:697–709`: `{\rho\over T_6^{7/2}}`; the following integral defines $c_7$.

Thus

$$
\kappa_r\propto \frac{T^{13/2}}{\rho^2}.
$$

For strongly degenerate electrons, PPP15 gives

$$
\kappa_e=\frac{\pi^2k_B^2Tn_e}{3m_e^\ast\nu_{ei}},
\qquad
\nu_{ei}=\frac{4Zm_e^\ast e^4\Lambda}{3\pi\hbar^3}.
$$

Sources:

> `refs/potekhin_pons_page_2015_thermal.tex:1485–1493`: `\kappa = a \,{n T}/{m^\ast\nu}`; “squared”;  
> `refs/potekhin_pons_page_2015_thermal.tex:3808–3811`: `4 Z \mel^\ast e^4 \Lambda_{e\mathrm{i}}`.

The word “squared” refers to the Boltzmann constant when restoring kelvin units. Omitting it is a dimensional error.

**DERIVED HERE, restricted homology calculation.** Suppose that over the matching region

$$
P=C_P\rho^\gamma,\qquad
\kappa_e=C_eT\rho^a,
$$

with constant composition and Coulomb logarithm. Radiative transport then gives

$$
\frac{dT}{d\rho}\propto\frac{F}{g}\rho^\gamma T^{-13/2},
$$

so, neglecting the surface integration constant,

$$
T^{15/2}\propto\frac{F}{g}\rho^{\gamma+1}.
$$

At the radiative–conductive transition,

$$
T_t^{13/2}\rho_t^{-2}
\propto T_t\rho_t^a,
\qquad
\rho_t\propto T_t^{11/[2(a+2)]}.
$$

Therefore

$$
\frac{F}{g}\propto T_t^{\alpha_t},
\qquad
\alpha_t=
\frac{15}{2}
-\frac{11(\gamma+1)}{2(a+2)}.
$$

The limiting exponents follow from the electron gas, rather than being inserted:

- State counting gives $n_e\propto p_F^3$.
- For nonrelativistic electrons, integrating momentum flux gives $P\propto p_F^5$, hence $\gamma=5/3$. The conductivity above gives $a=1$.
- For ultrarelativistic electrons, $P\propto p_F^4$, hence $\gamma=4/3$. Since $m_e^\ast\propto p_F$, conductivity gives $a=1/3$.

The state-counting and dispersion inputs are local:

> `refs/potekhin_chabrier_2013_mag_eos.tex:260–264`: `$\pF = \hbar \,(3 \pi^2 \nel)^{1/3}$`; `$\EF = c\,\sqrt{(\mel c)^2 + (\pF)^2}-\mel c^2`.

Python evaluates the illustrative matching exponents as

| Idealized limit | $\alpha_t$ |
|---|---:|
| Nonrelativistic, degenerate | $47/18=2.611111$ |
| Ultrarelativistic, degenerate | $2.000000$ |

These are **not** competing precision predictions for the real envelope. They show why radiative–conductive matching does not force the exponent to equal 2.2. Partial degeneracy, changing relativity, scattering, ionization and the moving transition all matter.

A second derivation makes the near-quadratic dependence more transparent. Let $x=p_F/(m_ec)$ and keep $\Lambda$ constant. Using $dP=n_e\,d\mu$ for the zero-temperature gas and the conductivity above gives

$$
T_b^2-T_t^2
=
\frac{24Z^2e^4\Lambda}{\pi A m_u k_B^2c}
\frac{F}{g}
\left[J(x_b)-J(x_t)\right],
$$

where, **DERIVED HERE**,

$$
J(x)=\operatorname{arsinh}x-\frac{\sqrt{1+x^2}}{x},
\qquad
J'(x)=\frac{\sqrt{1+x^2}}{x^2}.
$$

The relation $dP=n_e\,d\mu$ follows by differentiating the zero-temperature grand-potential integral; PC13 supplies its finite-temperature parent at `refs/potekhin_chabrier_2013_mag_eos.tex:1057–1062`, quoted above.

If the bracket and $T_t/T_b$ varied weakly, $F/g$ would scale nearly as $T_b^2$. Their variation changes that slope. Determining the actual slope requires the envelope solution.

## 2.4 Where 2.2 actually comes from

PPP15 states that Gudmundsson and collaborators **fitted numerical solutions**:

> `refs/potekhin_pons_page_2015_thermal.tex:2930–2933`: “fitted the numerical solutions by a remarkably simple formula”.

The printed relation is

$$
T_b=1.288\times10^8
\left(\frac{T_{s6}^4}{g_{14}}\right)^{0.455}\,\mathrm K.
$$

Algebraically,

$$
T_{s6}^4
=
g_{14}
\left(\frac{T_b}{1.288\times10^8\,\mathrm K}\right)^{1/0.455}.
$$

Python gives

$$
\boxed{
\alpha=\frac1{0.455}=2.1978021978,
\qquad
T_s\propto T_b^{0.5494505495}.
}
$$

**The exponent 2.2 is the rounded reciprocal of a fitted exponent.** The local material supports the mechanism and approximate scaling, but not an exact first-principles derivation of that decimal.

PPP15 points to Ventura & Potekhin for an analytical derivation of a *similar* expression:

> `refs/potekhin_pons_page_2015_thermal.tex:2937–2939`: “An analytical derivation of a similar expression”.

That paper is a missing local source, not permission to attribute an unexamined derivation to it.

# 3. Magnetic transport and the failure near the equator

## 3.1 Classical suppression is only one part of the problem

For a nonquantizing field, Potekhin 1999 gives

$$
\tau_\parallel=\tau_0,\qquad
\tau_\perp^{\rm transport}=\frac{\tau_0}{1+(\omega_g\tau_0)^2},
\qquad
\tau_H=\frac{\omega_g\tau_0^2}{1+(\omega_g\tau_0)^2}.
$$

Source:

> `refs/potekhin_1999_magnetic_conductivity.tex:635–645`: `\tau_{xx}={\tau_0\over 1+(\omg\tau_0)^2}`.

In the degenerate, effectively elastic limit this gives

$$
\frac{\kappa_\perp}{\kappa_\parallel}
=
\frac1{1+(\omega_g\tau_0)^2}.
$$

The resulting tendency is clear: radial field lines permit efficient radial electron transport; tangential lines suppress it.

A deliberately simplified numerical diagnostic uses fully ionized Fe, $A=56$, $Z=26$, $\rho=10^6\,\mathrm{g\,cm^{-3}}$, local $T=10^8\,\mathrm K$, and the **assumption** $\Lambda=1$.

Using PY01’s relativity parameter and PPP15’s collision-time and magnetization estimates,

> `refs/potekhin_yakovlev_2001_cooling.tex:178–180`: `1.009 \,(\rho_6\,Z/ A )^{1/3}`;  
> `refs/potekhin_pons_page_2015_thermal.tex:2694–2695`: `1760\,\frac{B_{12}}{\sqrt{1+x_\mathrm{r}^2}}`;  
> `refs/potekhin_pons_page_2015_thermal.tex:3810–3811`: `5.7\times10^{-17}`,

Python gives

$$
x_r=0.781303,\qquad
\tau_0=1.72755\times10^{-18}\,\mathrm s,
$$

$$
\omega_g\tau_0=239.591,\qquad
\frac{\kappa_\perp}{\kappa_\parallel}=1.74201\times10^{-5}.
$$

This is an order-of-magnitude diagnostic, **not** an acceptable conductivity at that state: it intentionally omits quantization, realistic Coulomb logarithms and radiative transport.

In particular, it does not imply $T_{\rm eq}/T_{\rm p}=(1.74201\times10^{-5})^{1/4}$. The conductivities vary with depth, the sensitivity strip moves, photons also transport heat, and neighbouring columns exchange energy.

## 3.2 Quantization intersects the relevant layers

PY01 supplies

$$
\rho_B\simeq7045\,B_{12}^{3/2}\frac AZ
\ \mathrm{g\,cm^{-3}},
\qquad
T_B\simeq
\frac{1.343\times10^8B_{12}}{\sqrt{1+x_r^2}}\ \mathrm K.
$$

Source:

> `refs/potekhin_yakovlev_2001_cooling.tex:156–168`: `temperature $T\ll T_B$ and density $\rho<\rho_B$`; `\approx 7045`.

For the same illustrative Fe state and field, Python gives

$$
\rho_B=4.79839\times10^5\,\mathrm{g\,cm^{-3}},
\qquad
T_B=1.05829\times10^9\,\mathrm K.
$$

The quantization boundary lies within the nonmagnetic sensitivity-strip density range quoted above. Thus an unmagnetized conductivity with a transverse Hall suppression factor is not an adequate replacement.

PPP15 explicitly uses different collision times in quantizing fields:

> `refs/potekhin_pons_page_2015_thermal.tex:2825–2844`: “two different effective mean-free times”; “different Coulomb logarithms”.

Landau filling changes the electron density of states, EOS, scattering phase space and conductivity. PPP15 identifies the first major longitudinal-conductivity peak with occupation of the first excited Landau level:

> `refs/potekhin_pons_page_2015_thermal.tex:2871–2874`: “filling of the first”; “higher Landau levels”.

Photon transport is also anisotropic. PPP15’s angular and mode-dependent Rosseland construction appears at `refs/potekhin_pons_page_2015_thermal.tex:2603–2655`, with the phrase “effective opacities”. A scalar nonmagnetic radiative opacity is another approximation that requires testing.

## 3.3 Why a thin envelope can still require two dimensions

**DERIVED HERE.** In a meridional section,

$$
F_r=-\kappa_{rr}\frac{\partial T}{\partial r}
-\kappa_{r\theta}\frac1r\frac{\partial T}{\partial\theta}.
$$

The stationary temperature field obeys

$$
\frac1{r^2}\frac{\partial(r^2F_r)}{\partial r}
+
\frac1{r\sin\theta}
\frac{\partial(\sin\theta F_\theta)}{\partial\theta}
=Q_{\rm heat}-Q_\nu.
$$

A small geometric thickness does not establish that the second term is negligible. Near a tangential field, radial transport samples the small transverse conductivity while meridional transport samples the much larger longitudinal conductivity.

Exactly at a reflection-symmetric equator, $\partial_\theta T$ can vanish while $\partial_\theta(\sin\theta F_\theta)$ does not. Symmetry therefore does not rescue an isolated equatorial column.

PPP15 states the failure directly:

> `refs/potekhin_pons_page_2015_thermal.tex:3220–3227`: “the 1D approximation fails”; “overestimates the heat blanketing effect”; “width can be estimated as $\sim10$\% of the radius”.

The physical consequence is an equator warmer than an isolated-column calculation predicts. The Appendix-B fit includes this effect.

**Project correction:** implement the 1-D solver as a separately named approximation. Replacing the whole PPP15 map requires lateral transport or an explicit retained lateral-transport stub.

## 3.4 Quantitative limits of GH83

Two different error questions must not be conflated.

**Interpolation error within a 1-D model.** PY01 reports that the fourth-power GH interpolation reproduces its numerical angular dependence within approximately 30%, while replacing the power by 4.5 improves this to 10%:

> `refs/potekhin_yakovlev_2001_cooling.tex:1309–1323`: `(within $\approx30$\%)`; “yields better accuracy (within 10\%)”.

These are temperature errors. Python converts symmetric temperature deviations into flux factors:

| Temperature deviation | Bolometric flux factor |
|---|---:|
| $-30\%,+30\%$ | $0.2401,\ 2.8561$ |
| $-10\%,+10\%$ | $0.6561,\ 1.4641$ |

These published bounds concern interpolation against the PY01 calculations. They do not bound the error from omitting two-dimensional transport.

**Comparison against the PPP15 map at the requested benchmark.** To isolate angular shape, retain the same PPP15 pole and equator temperatures and compare

$$
T_{\rm GH}^4
=
T_{\rm p}^4\cos^2\Theta_B+
T_{\rm eq}^4\sin^2\Theta_B
$$

with PPP15 Eq. B.6. For this diagnostic only, adopt the flat-space dipole relation

$$
\cos^2\Theta_B
=
\frac{\cos^2\theta}
{\cos^2\theta+\tfrac14\sin^2\theta}.
$$

The dipole input is local:

> `refs/greenstein_hartke_1983_surface_T.txt:161–163`: `tan O = 2 cot /?`.

The OCR is poor; the corresponding field–normal relation is independently printed in PY01’s GR form:

> `refs/potekhin_yakovlev_2001_cooling.tex:484–490`: `\mathrm{tan}\,\theta = a\,\mathrm{tan}\,\chi`.

The flat-space diagnostic takes its $a\to1/2$ limit. It is not a GR-calibrated RX J1856 calculation.

| Colatitude | PPP15 $T_s$ [K] | GH interpolation [K] | GH/PPP15 $-1$ |
|---:|---:|---:|---:|
| $30^\circ$ | 1,224,706 | 1,219,593 | $-0.42\%$ |
| $60^\circ$ | 1,119,277 | 1,082,280 | $-3.31\%$ |
| $75^\circ$ | 912,495 | 857,079 | $-6.07\%$ |
| $80^\circ$ | 753,433 | 721,348 | $-4.26\%$ |
| $85^\circ$ | 501,068 | 527,489 | $+5.27\%$ |

Endpoints agree by construction. This comparison measures shape differences, not the missing equatorial heating.

As a further guard against false precision, evaluating PY01’s older transverse fit at the flat-dipole equatorial field gives $176{,}188\,\mathrm K$, below its stated applicability floor:

> `refs/potekhin_yakovlev_2001_cooling.tex:1227–1231`: `$T_\mathrm{s}>2\times10^5$`.

That extrapolated number must **not** be presented as a validated 1-D temperature or used to assign a precise percentage to PPP15’s lateral correction. The local sources establish its sign and importance, but do not provide a trustworthy pointwise error bar for this exact benchmark.

# 4. Microphysics inventory and the real recomputation tree

The audit’s transport ODE is a small part of the work. Its dependency tree must distinguish exact quadratures, physical approximations, empirical fits and numerical verification fixtures.

| Layer | Required calculation | Local evidence and ownership issue |
|---|---|---|
| **L0: global boundary problem** | Solve for $T(r,\theta)$ and surface flux, with atmosphere and bottom boundaries | PPP15 `refs/potekhin_pons_page_2015_thermal.tex:3984–3987`: “The use of the 2D code corrects”. Independent columns do not own this node. |
| **L1: transport assembly** | Rotate electron and photon tensors; solve hydrostatics; combine conductivities | PY01 `refs/potekhin_yakovlev_2001_cooling.tex:617–619`: `\kappa=\kappa_\mathrm{r}+\kappa_\mathrm{c}`. Algebraically straightforward, but only as accurate and owned as its children. |
| **L2: electron transport moments** | Energy-average collision operators; include thermoelectric correction | PPP15 `refs/potekhin_pons_page_2015_thermal.tex:2740–2742`: `\hat{\kappa} = \hat{\tilde\kappa} + T\hat{\alpha}\cdot \hat{Q}`. Using $\tilde\kappa$ alone outside strong degeneracy omits a term. |
| **L3: electron–ion scattering** | Screened collision integrals, relativistic factors, Landau transitions and non-Born corrections | Potekhin 1999 `refs/potekhin_1999_magnetic_conductivity.tex:682–691`: `q^3 |\phi_q|^2 S(q) |F(q)|^2`. Local equations exist, but the effective-potential implementation contains fits. |
| **L3: electron–phonon scattering** | Phonon modes, Debye–Waller factor, inelastic weighting, normal/Umklapp and multiphonon processes | BY95 `refs/baiko_yakovlev_1995_ephonon.tex:240–295`: “three branches of phonon vibrations”; “first Brillouin zone”. Its convenient final formulae are fits. |
| **L3: electron–electron and impurity scattering** | Partial-degeneracy collision calculation, screening, composition-dependent defects | PPP15 `refs/potekhin_pons_page_2015_thermal.tex:3899–3907`: “there is no”; “such averaging procedure”; “interpolation”. Quantizing $ee$ transport is not supplied by the classical formula. |
| **L4: EOS and screening** | Landau electron gas, ionic free energy, electron exchange/correlation and polarization | PC13 `refs/potekhin_chabrier_2013_mag_eos.tex:414–449`: “electron-electron, ion-ion, and ion-electron interactions”. Recomputing the ideal gas does not own the interaction terms. |
| **L5: liquid correlations** | OCP equilibrium sampling, finite-size control, $S(q)$, free-energy integration; extend beyond rigid-background OCP when needed | PC13 `refs/potekhin_chabrier_2013_mag_eos.tex:432–437`: “rigid (nonpolarizable) background”. Classical OCP is a model, not the complete magnetized plasma. |
| **L5: crystal dynamics** | Compute dynamical matrix and phonon eigenmodes; integrate transport-weighted spectra | BY95 `refs/baiko_yakovlev_1995_ephonon.tex:338–349`: “we have interpolated the tabulated data and then solved”. Those published tables would be another runtime dependency. |
| **Radiative branch** | Frequency-, angle- and mode-dependent free-free and scattering opacity; ionization and bound processes when relevant | PPP15 `refs/potekhin_pons_page_2015_thermal.tex:2646–2648`: “free-free absorption and Thomson”; `:2674–2677`: “presence of bound species”. Shares the atmosphere opacity tree. |
| **Model inputs** | Composition versus depth, impurity content, field geometry, $T_b$; neutrino/heating terms where needed | PPP15 `refs/potekhin_pons_page_2015_thermal.tex:3975–3978`: “non-accreted envelope”; “BSk21 EoS”. A pure-Fe column is not automatically identical to this benchmark. |

## 4.1 Electron collision quadrature is not ownership if its potential is fitted

Potekhin 1999 supplies a nonmagnetic Coulomb integral,

$$
\Lambda_0
=
\int_0^{2p/\hbar}
dq\,q^3|\phi_q|^2S(q)|F_{\rm nuc}(q)|^2
\left[1-\left(\frac{\hbar cq}{2\epsilon}\right)^2\right],
$$

and a Landau-state linear system for longitudinal transport. The relevant local anchors are:

> `refs/potekhin_1999_magnetic_conductivity.tex:682–685`: `\Lambda_0 = \int_0^{2p_0/\hbar}`;  
> `refs/potekhin_1999_magnetic_conductivity.tex:811–814`: `\left[ \varphi_{ns}(\epsilon)`.

However, the paper then replaces the scattering factors by an effective potential containing phenomenological corrections:

> `refs/potekhin_1999_magnetic_conductivity.tex:713–719`: “the expression”;  
> `refs/potekhin_1999_magnetic_conductivity.tex:752`: “phenomenological corrections to the Born approximation”.

Numerically integrating that fitted potential removes one fit layer, not the underlying fit dependency.

Moreover, the relaxation-time approximation is explicitly limited:

> `refs/potekhin_1999_magnetic_conductivity.tex:697–705`: “fails if inelastic processes”; `$T\ll T_{\rm p}$`.

The solver must distinguish electrical and thermal collision operators when energy transfers matter.

## 4.2 An OCP energy match does not validate transport

Slattery, Doolen & DeWitt Table I provides excess internal energies:

> `refs/slattery_doolen_dewitt_1982_ocp.txt:39`: “U/NkT is the excess internal energy per ion divided by kT.”

At $\Gamma=100$, $N=1024$, the printed liquid result is

$$
\frac{U_{\rm ex}}{Nk_BT}=-87.522,\qquad \sigma=0.001.
$$

Anchor:

> `refs/slattery_doolen_dewitt_1982_ocp.txt:93`: `100.0   —87.522    0.001`.

This is an excellent energy verification target. It is not a supplied $S(q)$ table and not a transport validation.

**DERIVED HERE.** An interaction energy constrains a weighted integral of pair correlations. Distinct correlation functions can share that integral while differing at the momentum transfers emphasized by the collision integral. Therefore energy agreement does not imply Coulomb-logarithm agreement.

Likewise, equilibrium classical Monte Carlo samples spatial configurations. It does not supply the time-dependent correlation spectrum needed for general inelastic scattering.

For free energy, PC00 provides the integration route:

$$
f_{ii}(\Gamma)=
\int_0^\Gamma \frac{u_{ii}(\Gamma')}{\Gamma'}\,d\Gamma',
$$

with the weak-coupling boundary

$$
u_{ii}\sim-\frac{\sqrt3}{2}\Gamma^{3/2}.
$$

Source:

> `refs/potekhin_chabrier_2000_eos.tex:343–360`: “Debye-H\"uckel energy”; `\int_0^\Gamma`.

Python gives $u_{ii}=-2.738613\times10^{-5}$ for the leading term at the chosen test point $\Gamma=10^{-3}$. This is an asymptotic test, not a finite-$\Gamma$ exact value.

## 4.3 The bottom is not necessarily a classical liquid

For an illustrative pure-Fe envelope at local $T=10^8\,\mathrm K$, calculate

$$
\Gamma=\frac{Z^2e^2}{a_i k_BT},
\qquad
T_{pi}=\frac{\hbar}{k_B}
\sqrt{\frac{4\pi Z^2e^2n_i}{A m_u}}.
$$

Source:

> `refs/potekhin_chabrier_2000_eos.tex:194–214`: `\Gamma = {(Ze)^2\over k_B T a}`; “ion plasma frequency”.

Using CODATA constants gives:

| $\rho$ [$\mathrm{g\,cm^{-3}}$] | $\Gamma$ | $T_{pi}$ [K] | $T/T_{pi}$ |
|---:|---:|---:|---:|
| $10^6$ | 40.1923 | $3.63635\times10^6$ | 27.5001 |
| $10^{10}$ | 865.917 | $3.63635\times10^8$ | 0.275001 |

These are diagnostic states, not an integrated temperature profile. They demonstrate why a classical-liquid OCP calculation cannot cover the whole depth range.

The constants used are local:

> `refs/tiesinga_2021_codata2018.txt:3403–3428`: “speed of light in vacuum”; “Planck”; “fine-structure constant”; “atomic mass unit”.

In Gaussian units the calculation used $e^2=\alpha\hbar c$, consistent with the project’s notation decision.

## 4.4 Two particularly important remaining approximations

**Electron–electron collisions.** PPP15 acknowledges that the quantizing-field conductivities were known for electron–ion transport but not the electron–electron mechanism, and applies a correction inferred at zero field:

> `refs/potekhin_pons_page_2015_thermal.tex:4020–4026`: “electron-ion but not electron-electron collision mechanism”; “collisions at $B=0$”.

Reproducing PPP15 faithfully therefore does not mean reproducing complete magnetic transport physics. That uncertainty survives even a perfect implementation.

**Heavy-element radiative opacity.** The audit overstates the need to solve all Fe bound processes before making useful envelope progress. PY01 explicitly found that fully ionized radiative opacity could adequately approximate the envelope relation for sufficiently warm cases:

> `refs/potekhin_yakovlev_2001_cooling.tex:648–653`: “not noticeably changed”; “fully ionized iron”.

But this is a conditional approximation. It must be tested by perturbing the upper boundary and omitted opacity contributions, especially for the cold equator. Switching to H/He changes the physical envelope; it does not validate a heavy-element PPP15 comparison or automatically remove partial ionization at strong field.

The “about 1200 LoC” estimate in the audit is consequently not a meaningful completion estimate. It might describe a transport driver using fitted microphysics. It does not describe ownership of all these leaves.

# 5. Numerical sanity chain for the requested benchmark

## 5.1 Interpret the field correctly

PPP15 defines

$$
T_9=\frac{T_b}{10^9\,\mathrm K},
\qquad
B_{12}=\frac{B_{\rm p}}{10^{12}\,\mathrm G}.
$$

Source:

> `refs/potekhin_pons_page_2015_thermal.tex:4049–4050`: `$B_{12}=B_\mathrm{p}/10^{12}$`.

Thus the requested calculation uses

$$
T_9=0.1,\qquad g_{14}=2.4,\qquad B_{12}=10.
$$

This is a global dipole-fit evaluation with polar field $10^{13}\,\mathrm G$. It is not two independent patches both having local field magnitude $10^{13}\,\mathrm G$.

## 5.2 Evaluate the actual fit

PPP15 Eqs. B.1–B.4 are

$$
a=\frac{0.337}{1+0.02\sqrt{B_{12}}},
$$

$$
T_0=(15.7T_9^{3/2}+1.36T_9)^{0.3796},
\qquad
T_1=1.13B_{12}^{0.119}T_9^a,
$$

$$
T_{\rm p}^{(0)}
=
10^6\,\mathrm K
\left\{
g_{14}\left[
T_1^4+(1+0.15\sqrt{B_{12}})T_0^4
\right]
\right\}^{1/4},
$$

$$
T_{\rm p}^{\max}
=
10^6\,\mathrm K
\left(5.2g_{14}^{0.65}+0.093\sqrt{g_{14}B_{12}}\right),
$$

$$
T_{\rm p}
=
T_{\rm p}^{(0)}
\left[1+\left(\frac{T_{\rm p}^{(0)}}{T_{\rm p}^{\max}}\right)^4\right]^{-1/4}.
$$

Source:

> `refs/potekhin_pons_page_2015_thermal.tex:4036–4064`: `T_\mathrm{p}^{(0)}=\left[g_{14}`; `a=0.337/(1+0.02\sqrt{B_{12}})`.

The ratio is

$$
q\equiv\frac{T_{\rm p}}{T_{\rm eq}}
=
1+
\frac{(1230T_9)^{3.35}B_{12}\sqrt{1+2B_{12}^2}}
{(B_{12}+450T_9+119B_{12}T_9)^4}
+
\frac{0.0066B_{12}^{5/2}}
{T_9^{1/2}+0.00258B_{12}^{5/2}}.
$$

Source:

> `refs/potekhin_pons_page_2015_thermal.tex:4072–4077`: `\frac{T_\mathrm{p}}{T_\mathrm{eq}}`; `0.0066\,B_{12}^{5/2}`.

Python results:

| Quantity | Computed value |
|---|---:|
| $a$ | 0.316954064830 |
| $T_0$ | 0.840381990133 |
| $T_1$ | 0.716348537639 |
| $T_{\rm p}^{(0)}$ | 1,244,260.270837 K |
| $T_{\rm p}^{\max}$ | 9,641,905.653476 K |
| First additive term in $q-1$ | 1.550897724955 |
| Second additive term in $q-1$ | 1.843575418994 |
| **$T_{\rm p}$** | **1,244,174.018792 K** |
| **$q$** | **4.394473143950** |
| **$T_{\rm eq}$** | **283,122.453599 K** |

The neutrino correction lowers the polar temperature by $0.006932\%$ and its flux by $0.027725\%$. This is consistent with PPP15’s statement:

> `refs/potekhin_pons_page_2015_thermal.tex:3123–3124`: `$T_\mathrm{b}\lesssim10^8$~K the neutrino emission does not`.

The fit correction is negligible for an initial calculation at this benchmark. That does not authorize omitting neutrinos throughout the much wider Appendix-B temperature domain.

## 5.3 Arithmetic verification is not physical validation

The expected values below are independent high-precision evaluations of the **same published formulae**.

| Quantity | 60-digit evaluation, rounded | Float64 evaluation, rounded | Relative difference | Arithmetic tolerance |
|---|---:|---:|---:|---:|
| $T_{\rm p}$ [K] | 1,244,174.01879183833 | 1,244,174.0187918385 | $1.87\times10^{-16}$ | $10^{-12}$ |
| $q$ | 4.39447314394953759 | 4.3944731439495381 | $2.02\times10^{-16}$ | $10^{-12}$ |
| $T_{\rm eq}$ [K] | 283,122.453599440036 | 283,122.45359944005 | below displayed Float64 resolution | $10^{-12}$ |

Computing $\sigma_{\rm SB}$ from the CODATA constants gives

$$
\sigma_{\rm SB}
=
5.6703744191844\times10^{-5}
\ \mathrm{erg\,cm^{-2}\,s^{-1}\,K^{-4}},
$$

consistent with the printed constant:

> `refs/tiesinga_2021_codata2018.txt:3424`: `5.670 374 419`.

The corresponding local fluxes are

$$
F_{\rm p}=1.358739650\times10^{20},
\qquad
F_{\rm eq}=3.643417919\times10^{17}
\quad
\mathrm{erg\,cm^{-2}\,s^{-1}},
$$

$$
\boxed{\frac{F_{\rm p}}{F_{\rm eq}}=q^4=372.929946.}
$$

For comparison, the Gudmundsson fit gives $T_s=1{,}083{,}077.36\,\mathrm K$ at the same $T_b,g$. PY01’s updated nonmagnetic expression gives $1{,}020{,}457.52\,\mathrm K$, lower by $5.78166\%$:

> `refs/potekhin_yakovlev_2001_cooling.tex:1200–1208`: `[(7\zeta)^{2.25}+(\zeta/3)^{1.25}]`.

This is a practical warning: requiring a modern solver to agree with every historical fit to a numerical tolerance would be scientifically wrong.

# 6. Concrete implementation and acceptance plan

The implementation should proceed breadth-first, but its provenance must remain explicit.

| Function or module | Inputs | Outputs and required checks |
|---|---|---|
| `electron_eos_landau` | $\rho,T,B,\{A_i,Z_i,X_i\}$ | $\mu_e,P_e$ and derivatives; controlled Landau and momentum truncation; nonmagnetic and classical limits; correct ground-level spin weight. |
| `ion_correlations` | $\Gamma$, particle number, phase, screening model | $U_{\rm ex}$, $S(q)$, uncertainty and finite-size diagnostics; OCP energy checks, then separate structure-factor checks. |
| `phonon_spectrum` | lattice, composition, density, screening | Frequencies and eigenvectors; translational modes, symmetry and quadrature convergence; no imported dynamical-matrix tables at runtime. |
| `collision_operator` | EOS state, correlations, phonons, $B$ | Energy-dependent longitudinal/transverse transport responses and validity flags; fitted potentials remain marked stubs. |
| `electron_transport` | collision responses and EOS | $\sigma,\alpha,\tilde\kappa,\kappa$; thermoelectric subtraction; positive symmetric conductivity; nonquantizing limit. |
| `radiative_transport` | spectral mode opacities, $\rho,T,B$ | Radiative conductivity tensor; isotropic limit must recover the ordinary Rosseland result without a factor-of-two error. |
| `envelope_column` | $T_b,g,B,\Theta_B,\rho_b$, composition, atmosphere boundary | $T(P),\rho(P),F,T_s$ and residuals; useful independently of the global map. |
| `envelope_axisymmetric` | $M,R,\boldsymbol B(r,\theta),T_b(\theta)$, composition and boundaries | $T(r,\theta),\boldsymbol F,T_s(\theta)$; equatorial regularity, lateral transport and integrated energy conservation. |
| `envelope_render_validation` | map, emergent-intensity model, geometry, measured passbands | Absolute Earth-distance photometry and resolved radiance comparisons, with envelope-model perturbations. |

### Required numerical gates

The following tolerances are **proposed engineering criteria**, not uncertainties claimed by the papers.

1. **Fit transcription fixture.** Reproduce the three arithmetic values above to relative $10^{-12}$. This tests formula transcription only.

2. **Tensor tests.** Test parallel and transverse limits, arbitrary rotations, and zero normal Hall contribution. Mutation tests must catch exchanged sine/cosine factors and use of colatitude in place of the field–normal angle.

3. **Column convergence.** Require relative $T_s$ changes below $10^{-3}$ under tighter integration and boundary matching, and relative flux conservation residual below $10^{-4}$ in the source-free problem.

4. **OCP fixture.** At $\Gamma=100,N=1024$, compare with $-87.522\pm0.001$. A provisional acceptance rule is agreement within three combined standard errors, including the simulation’s own uncertainty. The source-only three-sigma contribution is $0.003$. Do not compare a smaller particle number with this target without finite-size treatment.

5. **Nonmagnetic envelope comparison.** Use $1.083077\times10^6\,\mathrm K$ from Gudmundsson and $1.020458\times10^6\,\mathrm K$ from PY01 as distinct model targets. Record differences against the actual adopted microphysics. A discrepancy between them is not an integration error.

6. **Magnetic comparison.** Use PPP15’s $T_{\rm p}=1.244174\times10^6\,\mathrm K$ and $T_{\rm eq}=2.831225\times10^5\,\mathrm K$ as screening benchmarks. A provisional polar discrepancy threshold of 5% corresponds to $62{,}209\,\mathrm K$, but is not a published error bound. There is no defensible universal 5% equatorial tolerance at this field in the local text.

7. **Two-dimensional validation.** Compare matched assumptions against PPP15 Fig. 11 and, preferably, the underlying numerical data. Digitized curves require stated digitization errors. Demonstrate angular and depth convergence and conservation independently of agreement with a fit.

8. **Photometric stability.** Allocate at most 1% band-flux variation to numerical envelope error, equivalent to about $0.0108$ mag. Physical-model uncertainty is separate and may be larger.

PPP15 does **not** warrant a precision equatorial benchmark at $10^{13}\,\mathrm G$. Its ratio fit has errors up to 10% at smaller fields and deviations reaching 35% in the high-field oscillatory regime:

> `refs/potekhin_pons_page_2015_thermal.tex:4080–4090`: “deviations”; “reaching up to 35\%”; “up to”; `10\%`.

Neither bound is a calibrated uncertainty interval at the requested point. Before retiring the fit, the project needs a comparison with actual numerical solutions and a documented approximation budget, not merely agreement with Eq. B.5.

### Tests before production code

Follow the project’s red–green rule. Start with failing analytic-limit and conservation tests, then implement each kernel. Fuzz positive physical inputs across degeneracy, phase and Landau thresholds. Mutate:

- the hydrostatic/flux sign;
- opacity addition versus inverse-opacity addition;
- the Landau ground-state multiplicity;
- the fourth power in flux;
- field units and polar-versus-local field arguments;
- the thermoelectric subtraction;
- the angle used in tensor rotation.

Do not demand monotonicity in $B$: the local sources explicitly describe quantization oscillations.

The initial fully ionized column model can be a useful owned approximation. Its parent remains yellow while fitted collision physics, EOS corrections or lateral transport are still below it.

# 7. Missing sources and verification material to fetch

**All items below are UNVERIFIED as local scientific sources until the full text is stored, checked and read.** Bibliographic leads are not evidence for an equation. The three linked arXiv records were checked for identification only.

| Priority | Source to fetch | Purpose and local lead |
|---|---|---|
| Essential | **Gudmundsson, Pethick & Epstein (1983), ApJ 272, 286**, *Structure of neutron star envelopes*. arXiv ID not established. | Original nonmagnetic calculations, microphysics and sensitivity-strip tests. Local bibliography: `refs/potekhin_pons_page_2015_thermal.tex:4531–4533`, “Structure of neutron star envelopes”. |
| Essential | **Ventura & Potekhin (2001)**, in *The Neutron Star–Black Hole Connection*, NATO ASI, p. 393; [arXiv:astro-ph/0104003](https://arxiv.org/abs/astro-ph/0104003). | Analytical matching derivation. Local lead: `refs/potekhin_yakovlev_2001_cooling.tex:2208–2213`, `[astro-ph/0104003]`. |
| Essential | **Potekhin, Baiko, Haensel & Yakovlev (1999), A&A 346, 345**; [arXiv:astro-ph/9903127](https://arxiv.org/abs/astro-ph/9903127). | The unified collision model underlying Potekhin 1999’s effective potential. Local lead: `refs/potekhin_pons_page_2015_thermal.tex:4983–4987`, “Transport properties of degenerate electrons”. |
| Essential | **Baiko, Kaminker, Potekhin & Yakovlev (1998), PRL 81, 5556**; [arXiv:physics/9811052](https://arxiv.org/abs/physics/9811052). | Multiphonon and effective liquid structure factors. Local lead: `refs/potekhin_pons_page_2015_thermal.tex:4211–4215`, “Ion structure factors and electron transport”. |
| Essential for the global map | **Potekhin, Chabrier & Yakovlev (2007), Ap&SS 308, 353**, and **Viganò et al. (2013), MNRAS 434, 123**. arXiv IDs not established here. | Original column and multidimensional methods; obtain actual boundary conditions and numerical verification data. Local leads: `refs/potekhin_pons_page_2015_thermal.tex:5000–5003`, “Heat blanketing envelopes”; `:5260–5263`, `\MNRAS{434}{123}{2013}`. |
| Essential for OCP implementation | **Slattery, Doolen & DeWitt (1980), PRA 21, 2087**. arXiv ID not established. | The calculation method referenced by the local 1982 paper. Local lead: `refs/slattery_doolen_dewitt_1982_ocp.txt:33–34`, “same”; “Ref. 1”, with reference at `:256–257`. |
| Essential for owned screening/crystal dynamics | **Jancovici (1962), Nuovo Cimento 25, 428**; **Cohen & Keffer (1955), Phys. Rev. 99, 1128**; **Pollock & Hansen (1973), PRA 8, 3110**. Bibliographic spellings require verification; no arXiv IDs established. | Electron screening and dynamical-matrix derivations. Local bibliography: `refs/baiko_yakovlev_1995_ephonon.tex:916–928`, “Keffer”, “Jankovici”, “Pollock”. |
| Required for a controlled $ee$ treatment | **Shternin & Yakovlev (2006), PRD 74, 043004**; **Cassisi et al. (2007), ApJ 661, 1094**. arXiv IDs not established. | Dynamical screening and partial-degeneracy treatment. Local leads: `refs/potekhin_pons_page_2015_thermal.tex:5083–5086`, “collisions between degenerate electrons”; `:4336–4340`, “Updated electron-conduction opacities”. |
| Conditional on sensitivity tests | **Chugunov (2012), Astronomy Letters 38, 25**; **Chugunov & Haensel (2007), MNRAS 381, 1143**. arXiv IDs not established. | Low-temperature lattice transport and ion heat transport, particularly when transverse electron conduction is small. Local leads: `refs/potekhin_pons_page_2015_thermal.tex:4377–4385`, “low temperatures”; “Thermal conductivity of ions”. |
| Conditional on upper-layer sensitivity | **Seaton et al. (1994), MNRAS 266, 805**, Opacity Project material; magnetic heavy-element atomic calculations must be identified separately. All metadata and applicability unverified here. | Heavy-element bound opacity if omission affects $T_s$ or the atmosphere match. A nonmagnetic opacity compilation alone does not solve the magnetic atomic problem. |

A separate unresolved source requirement is a controlled quantizing-field electron–electron collision treatment. PPP15’s admission at lines 4020–4026 must not be silently replaced by a claim that the existing references solve it.

The verification data request should include PPP15’s temperature maps, conductivity curves and boundary assumptions. Those outputs would be **verification fixtures**, never runtime inputs.

# 8. Stub policy: effect on the hero render and R5

## 8.1 A labelled PPP15 model is acceptable for an interim hero render

The label must state:

- non-accreted heavy-element envelope;
- dipole geometry;
- fitted two-dimensional temperature map;
- assumed $T_b$, polar field and gravity;
- retained uncertainty in equatorial transport and microphysics.

This is more defensible than replacing PPP15 with an insufficiently validated 1-D model and calling the result first principles. It does not satisfy the project’s final ownership rule: labelling a fit does not turn it green.

The benchmark is also not a fitted RX J1856 solution. The supplied gravity does not determine both mass and radius, and the selected $T_b$ has not been constrained against the observations here.

## 8.2 The visible contrast is much smaller than the bolometric contrast

For a blackbody diagnostic, the local Planck law is available at

> `refs/marr_2011_planck_law.tex:42`: `B_\lambda =\frac{2hc^2}{\lambda^5}`.

**DERIVED HERE**, at fixed wavelength,

$$
\frac{\partial\ln B_\lambda}{\partial\ln T}
=
\frac{x}{1-e^{-x}},
\qquad
x=\frac{hc}{\lambda k_BT}.
$$

Thus $B_\lambda\propto T$ in the Rayleigh–Jeans limit, rather than $T^4$.

Using the benchmark temperatures and the CODATA-derived $hc/k_B=0.014387768775039\,\mathrm{m\,K}$:

| Diagnostic | Pole/equator contrast |
|---|---:|
| Temperature / Rayleigh–Jeans limit | 4.39447 |
| Local blackbody radiance at 475 nm | 4.58193 |
| Local blackbody radiance at 150 nm | 5.03153 |
| Bolometric flux | 372.92995 |

The values at 475 and 150 nm are illustrative monochromatic calculations, not HST passband predictions or predictions for a magnetized atmosphere.

The hero image therefore depends on the emitted spectrum, viewing angle and relativistic frequency transformation, not simply on a map of $\sigma T_s^4$.

## 8.3 “Negligible total flux” does not establish negligible optical error

PPP15 says of the cold equatorial regions:

> `refs/potekhin_pons_page_2015_thermal.tex:3229`: “their contribution to the total flux is negligible.”

For a symmetric map let $u=|\cos\theta|$. Its area average is

$$
\langle f\rangle=\int_0^1 f(u)\,du.
$$

PPP15’s angular profile is

$$
T_s(u)=T_{\rm eq}+(T_{\rm p}-T_{\rm eq})
\frac{(1+a_1+a_2)u^2}{1+a_1u+a_2u^2},
$$

$$
a_1=\frac{a_2T_9^{1/2}}3,
\qquad
a_2=\frac{10B_{12}}{T_9^{1/2}+0.1B_{12}T_9^{-1/4}}.
$$

Source:

> `refs/potekhin_pons_page_2015_thermal.tex:4101–4110`: `\frac{\Ts(\theta)-T_\mathrm{eq}}{T_\mathrm{p}-T_\mathrm{eq}}`.

For this benchmark, Python gives $a_1=5.032651907$ and $a_2=47.743928091$, with

$$
\langle T_s\rangle=1.005930095\times10^6\,\mathrm K,
\qquad
\langle T_s^4\rangle^{1/4}=1.081655616\times10^6\,\mathrm K.
$$

The normalized first and fourth moments were $0.808512378133$ and $0.571256987682$. Tightening quadrature changed them by less than $10^{-15}$ in absolute normalized units; estimated quadrature errors were below $8\times10^{-13}$.

For illustration, interpret the quoted band width as a **full** meridional width $0.1R$, so its half-width is $0.05$ radians. This convention gives:

| Equatorial-band diagnostic | Fraction of whole-surface quantity |
|---|---:|
| Area | 4.9979% |
| Rayleigh–Jeans luminosity | 1.5771% |
| Bolometric luminosity | 0.03895% |

The numerical contrast supports PPP15’s bolometric argument while showing why it does not automatically establish an optical or resolved-image error bound.

Replacing the PPP15 angular profile by the endpoint-matched GH diagnostic above changes the area-averaged Rayleigh–Jeans luminosity by $-2.1056\%$, equivalent to $+0.02310$ mag, and bolometric luminosity by $-7.4903\%$. These are whole-surface diagnostics, not a ray-traced observer flux.

## 8.4 An illustrative uncertainty ensemble

PPP15’s quoted 35% deviation is not a calibrated uncertainty at this benchmark. Nevertheless, it provides a useful **stress-test scale**, explicitly labelled as such.

Hold $T_{\rm p}$ and the B.6 shape coefficients fixed and perturb the ratio $q$ by $\pm35\%$:

| Ratio perturbation | $T_{\rm eq}$ [K] | $q$ | Rayleigh–Jeans luminosity change | Magnitude change |
|---|---:|---:|---:|---:|
| $q\rightarrow0.65q$ | 435,573 | 2.85641 | $+3.7570\%$ | $-0.04004$ |
| Nominal | 283,122 | 4.39447 | — | — |
| $q\rightarrow1.35q$ | 209,720 | 5.93254 | $-1.8089\%$ | $+0.01982$ |

At 475 nm the local blackbody pole/equator contrast spans 2.9222–6.3050 in this ensemble. That is a substantial uncertainty in a resolved temperature gradient even when integrated optical flux changes by only a few percent.

These variations are neither a confidence interval nor a bound on arbitrary field geometry, composition changes or atmospheric beaming. Those require separate model variants.

## 8.5 Consequence for R5

The magnitude response follows from the locally printed STMAG definition:

> `refs/kaplan_2011_rxj1856_spectrum.tex:488`: `${\rm STMAG}=-2.5\log_{10} F_{\lambda}-21.1`.

**DERIVED HERE**,

$$
\Delta m=-2.5\log_{10}\left(\frac{F'}F\right).
$$

For context, Walter reports

> `refs/walter_2010_rxj1856_distance.tex:862–866`: `25.36$\pm$0.15`; “standard deviation”.

The $0.15$ mag uncertainty corresponds to approximately $+14.82\%/-12.90\%$ in flux. Passing that one comparatively loose constraint would not validate the map.

Kaplan’s fitted optical/UV summary is much tighter:

> `refs/kaplan_2011_rxj1856_spectrum.tex:1035`: `$4.25\pm0.07$`;  
> `:1038–1041`: “absorbed flux density”; “for the fit”.

Its fractional normalization uncertainty is $1.6471\%$, approximately $0.01788$ mag. This is a **fitted summary**, not a substitute for the original per-filter measurements. The GH shape difference and the illustrative ratio ensemble are comparable to or larger than that scale.

The R5 consequence is therefore conditional:

- A PPP15 stub can support an honest Earth-distance forward-model test, provided the surface parameters are fitted by this project and the envelope assumption is disclosed.
- Passing photometry does not uniquely validate the pole-to-equator map. In the Rayleigh–Jeans regime, unresolved flux constrains a weighted temperature average; different maps can share that average.
- The quoted calculations suggest envelope-map differences of order a few percent in integrated optical flux, while allowing much larger changes in resolved equatorial brightness.
- No R5 pass can be claimed from this calculation: actual passbands, extinction, atmosphere spectra, GR visibility and viewing geometry have not been integrated.

**Disposition:** retain PPP15 as the labelled interim baseline, carry an explicit map uncertainty ensemble, implement and validate the column solver, and add axisymmetric lateral transport before claiming that the full `pipeline.envelope-tb-ts` path has been recomputed.