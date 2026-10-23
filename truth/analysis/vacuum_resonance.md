**Abstract.** Vacuum conversion must act at the resonance depth during radiative transfer. The v1 post-processing swap cannot reproduce its effect on the atmospheric spectrum or thermal structure. For a Ho-like thin hydrogen layer on RX J1856, however, the resonance lies above the continuum-forming region, making its direct effect on optical/UV intensity small while leaving potentially important changes in X-ray polarisation. Derivation from the local sources also exposes three discrepancies: PLCH04’s printed normalisation constant is slightly inaccurate; the early Lai–Ho paper has the wrong ion factor in its printed adiabatic energy; and the commonly quoted coefficient $2.52$ keV does not exactly follow from the local amplitude equations and CODATA constants.

## 1. Scope, conventions, and numerical provenance

This argument concerns the atmosphere contribution to `pipeline.vacuum-mode-conversion`. It does not establish the star’s magnetic-field map, solve its partially ionised atmosphere, or propagate polarisation through the exterior magnetosphere.

I distinguish:

- **SOURCED:** equations or physical assumptions printed locally, with a line reference and a greppable quotation.
- **DERIVED HERE:** deductions from those equations.
- **MODEL ESTIMATE:** calculations under explicitly simplified atmospheric assumptions.
- **UNVERIFIED:** missing references or conclusions requiring a calculation not performed here.

All numerical evaluations below were run in Python. Special-function calculations used `mpmath` at 60 decimal digits; the independent mode-evolution calculation used SciPy’s `DOP853`. These were analysis calculations, not an implementation of the production pipeline.

The numerical inputs are CODATA 2018 values in [refs/tiesinga_2021_codata2018.txt:3403](/home/tobias/Projects/neutronstars/refs/tiesinga_2021_codata2018.txt:3403), through line 3428. Greppable anchors include `299 792 458`, `6.626 070 15`, `9.109 383 7015`, `1.672 621 923 69`, `7.297 352 5693`, and `1.380 649`. Converted to CGS,

$$
\begin{aligned}
c&=2.99792458\times10^{10}\ {\rm cm\,s^{-1}},\\
h&=6.62607015\times10^{-27}\ {\rm erg\,s},\qquad \hbar=h/(2\pi),\\
m_e&=9.1093837015\times10^{-28}\ {\rm g},\\
m_p&=1.67262192369\times10^{-24}\ {\rm g},\\
\alpha_f&=7.2973525693\times10^{-3},\\
k_B&=1.380649\times10^{-16}\ {\rm erg\,K^{-1}}.
\end{aligned}
$$

I compute the Gaussian charge as $e=\sqrt{\alpha_f\hbar c}$, consistently with the local definition `$\alphaf=e^2/\hbar c` in [refs/potekhin_lai_chabrier_ho_2004_polarization.tex:1038](/home/tobias/Projects/neutronstars/refs/potekhin_lai_chabrier_ho_2004_polarization.tex:1038). This gives

$$
e=4.803204713884971\times10^{-10}\ {\rm esu},\qquad
B_Q=\frac{m_e^2c^3}{e\hbar}
=4.414005217486694\times10^{13}\ {\rm G}.
$$

All photon energies in the atmosphere calculations are **local surface-frame energies** unless explicitly labelled otherwise.

## 2. Euler–Heisenberg theory supplies the coefficients without a fit

### 2.1 The required tensor is the linear response to the wave

**SOURCED.** Heyl–Hernquist define the invariants and effective Lagrangian in [refs/heyl_hernquist_1997_qed.tex:109](/home/tobias/Projects/neutronstars/refs/heyl_hernquist_1997_qed.tex:109), through line 159. Relevant anchors are `I = F_{\mu\nu} F^{\mu\nu}` and `\label{eq:heweak}`. Their invariants are

$$
\mathcal I=2(B^2-E^2),\qquad
\mathcal K=-16(\mathbf E\cdot\mathbf B)^2.
$$

Normalise the Lagrangian so that its Maxwell term is $(E^2-B^2)/2$. The printed weak-field expansion becomes

$$
\mathcal L=
\frac{E^2-B^2}{2}
+
\frac{\alpha_f}{90\pi B_Q^2}
\left[(B^2-E^2)^2+7(\mathbf E\cdot\mathbf B)^2\right]
+\text{higher field powers}.
$$

Here and below the normalisation is equivalent to multiplying the Gaussian electromagnetic Lagrangian by $4\pi$; the resulting relative dielectric and permeability tensors are dimensionless.

The crucial instruction is printed at [refs/heyl_hernquist_1997_qed.tex:252](/home/tobias/Projects/neutronstars/refs/heyl_hernquist_1997_qed.tex:252):

> `To examine wave propagation, we must first linearize the relations`

The response tensors for a weak wave on a static magnetic background are therefore

$$
\epsilon^{\rm wave}_{ij}
=\frac{\partial^2\mathcal L}{\partial E_i\partial E_j},
\qquad
(\mu^{-1})^{\rm wave}_{ij}
=-\frac{\partial^2\mathcal L}{\partial B_i\partial B_j},
$$

as printed at lines 258–273.

**DERIVED HERE.** Put $C=\alpha_f/(90\pi B_Q^2)$ and differentiate before setting the background electric field to zero:

$$
\mathbf D
=
\mathbf E+C\left[-4(B^2-E^2)\mathbf E
+14(\mathbf E\cdot\mathbf B)\mathbf B\right],
$$

$$
\mathbf H
=
\mathbf B-C\left[4(B^2-E^2)\mathbf B
+14(\mathbf E\cdot\mathbf B)\mathbf E\right].
$$

Thus

$$
\epsilon^{\rm wave}_{ij}
=(1-4CB^2)\delta_{ij}+14CB_iB_j,
$$

$$
(\mu^{-1})^{\rm wave}_{ij}
=(1-4CB^2)\delta_{ij}-8CB_iB_j.
$$

Writing

$$
\epsilon^{\rm vac}=a\mathbf1+q\hat{\mathbf B}\hat{\mathbf B},
\qquad
(\mu^{-1})^{\rm vac}=a\mathbf1+m\hat{\mathbf B}\hat{\mathbf B},
\qquad
\hat a=a-1,
$$

gives

$$
\boxed{
\hat a=-2\delta,\qquad q=7\delta,\qquad m=-4\delta,\qquad
\delta=\frac{\alpha_f}{45\pi}\left(\frac{B}{B_Q}\right)^2.
}
$$

These reproduce PLCH04’s expressions at [refs/potekhin_lai_chabrier_ho_2004_polarization.tex:1031](/home/tobias/Projects/neutronstars/refs/potekhin_lai_chabrier_ho_2004_polarization.tex:1031), with anchor `\hat{a}= - \frac{2 \alphaf}{45\pi}`.

**Implementation consequence:** using the unlinearised relation between the background $\mathbf H$ and $\mathbf B$ would lose the longitudinal magnetic response $m$. That changes $q+m$, hence the resonance itself.

### 2.2 Exact one-loop coefficients and their reduction to PLCH04

Let $b=B/B_Q$ and $x=1/b$. This avoids a notation trap: Heyl–Hernquist use $\xi=b$ for a magnetic background, whereas PLCH04 use $\xi=b^{-1}$.

**SOURCED.** The closed form is printed in [refs/heyl_hernquist_1997_qed.tex:176](/home/tobias/Projects/neutronstars/refs/heyl_hernquist_1997_qed.tex:176), through line 199; the anchor `X_0(x) & = & 4 \int_0^{x/2-1}` locates it:

$$
\mathcal L_1(\mathcal I,0)
=\frac{\alpha_f}{2\pi}\frac{\mathcal I}{2}X_0(x),
$$

$$
\begin{aligned}
X_0(x)={}&4\int_1^{x/2}\ln\Gamma(t)\,dt
-\frac13\ln x+2\ln(4\pi)-4\ln A-\frac53\ln2\\
&-[\ln(4\pi)+1-\ln x]x
+\left[\frac34+\frac12\ln\frac2x\right]x^2,
\end{aligned}
$$

$$
X_1(x)=-2X_0(x)+xX_0'(x)+\frac23X_0''(x)-\frac{2}{9x^2},
\qquad
\ln A=\frac1{12}-\zeta'(-1).
$$

Their wave-response tensors at [refs/heyl_hernquist_1997_qed.tex:295](/home/tobias/Projects/neutronstars/refs/heyl_hernquist_1997_qed.tex:295), through line 307, imply

$$
\begin{aligned}
\hat a&=\frac{\alpha_f}{2\pi}[-2X_0+xX_0'],\\
q&=-\frac{\alpha_f}{2\pi}X_1,\\
m&=\frac{\alpha_f}{2\pi}[xX_0'-x^2X_0''].
\end{aligned}
$$

**DERIVED HERE.** Differentiating $X_0$ gives

$$
X(x)\equiv X_0'(x)
=
2\ln\Gamma(x/2)-\frac1{3x}
-\ln\frac{4\pi}{x}+x+x\ln\frac2x,
$$

$$
X'(x)=\psi(x/2)+\frac1{3x^2}+\frac1x+\ln\frac2x.
$$

Since

$$
X_0(x)=X_0(1)+\int_1^xX(t)\,dt,
$$

the coefficients become

$$
\boxed{
\hat a=\frac{\alpha_f}{2\pi}
\left[xX(x)-2\int_1^xX(t)\,dt-C_0\right],
\qquad C_0=2X_0(1),
}
$$

$$
\boxed{
\hat a+q=\frac{\alpha_f}{2\pi}
\left[\frac2{9x^2}-\frac23X'(x)\right],
\qquad
m=\frac{\alpha_f}{2\pi}[xX(x)-x^2X'(x)].
}
$$

These are PLCH04 A2–A6, at [refs/potekhin_lai_chabrier_ho_2004_polarization.tex:1046](/home/tobias/Projects/neutronstars/refs/potekhin_lai_chabrier_ho_2004_polarization.tex:1046). The printed anchor is

> `\xi X(\xi) - 2 \int_1^\xi X(\xi')\, \dd \xi' - 0.0329199`

### 2.3 The origin—and numerical defect—of $0.0329199$

**DERIVED HERE.** The constant is the integration normalisation inherited from the renormalised Euler–Heisenberg action:

$$
C_0=
-8\int_{1/2}^1\ln\Gamma(t)\,dt
+2\ln(4\pi)-8\ln A-\frac12-\frac73\ln2.
$$

It is not an empirical coefficient.

The magnetic-background limit of the proper-time integral printed at [refs/heyl_hernquist_1997_qed.tex:134](/home/tobias/Projects/neutronstars/refs/heyl_hernquist_1997_qed.tex:134), anchor `\int_0^\infty e^{-\zeta}`, gives an independent route:

$$
\mathcal L_1(B)
=
-\frac{\alpha_fB_Q^2}{2\pi}
\int_0^\infty\frac{e^{-t}}{t^3}
\left[bt\coth(bt)-1-\frac{(bt)^2}{3}\right]dt.
$$

Consequently,

$$
C_0=-2\int_0^\infty\frac{e^{-t}}{t^3}
\left[t\coth t-1-\frac{t^2}{3}\right]dt.
$$

The subtraction of the field-independent and quadratic terms fixes the vacuum and charge normalisations. Changing $C_0$ adds a quadratic-field term to the action and shifts the zero-field response.

The numerical comparison is:

| Quantity | Expected/source | Computed | Difference |
|---|---:|---:|---:|
| $\ln A$ | HH’s exact zeta expression | $0.2487544770337843$ | — |
| $C_0$ | PLCH04 prints **$0.0329199$** | **$0.03291978776184018$** | $-1.12238159823\times10^{-7}$ |
| $C_0$, independent proper-time route | Closed-form calculation | $0.03291978776184018$ | $-1.49\times10^{-57}$ in this numerical evaluation |

The PLCH04 number is close, but it is **not the correctly rounded seven-decimal value** of the local Heyl–Hernquist expression. Its discrepancy exceeds the half-last-digit tolerance $5\times10^{-8}$.

Using the printed constant literally leaves the spurious weak-field offset

$$
\Delta\hat a
=-\frac{\alpha_f}{2\pi}(0.0329199-C_0)
=-1.30354491220\times10^{-10}.
$$

The production implementation must calculate the normalisation or use a cancellation-safe equivalent. It must not copy the decimal and then clamp the result to zero.

### 2.4 Independent checks and numerical stability

**DERIVED HERE.** Expanding

$$
z\coth z-1-\frac{z^2}{3}
=-\frac{z^4}{45}+\frac{2z^6}{945}-\cdots
$$

inside the proper-time integral recovers the quartic action and the weak-field coefficients. For numerical quadrature, the expansion also removes catastrophic cancellation near the lower integration endpoint.

At $b=10^{-4}$, the exact closed form gives

$$
\frac{\hat a}{-2\delta}=0.9999999914285717,\quad
\frac{q}{7\delta}=0.9999999946938777,\quad
\frac{m}{-4\delta}=0.9999999828571435.
$$

All agree with the weak-field limit within a chosen relative tolerance $2\times10^{-8}$.

I also differentiated the proper-time integrand to calculate $\hat a$ and $m$ independently. At the two requested fields, the relative discrepancies from the closed form were below $5\times10^{-53}$ in the high-precision calculation. This is numerical agreement between representations, not a claim about QED’s physical accuracy.

The PLCH04 fits are useful verification targets. Their stated maxima appear at [refs/potekhin_lai_chabrier_ho_2004_polarization.tex:1099](/home/tobias/Projects/neutronstars/refs/potekhin_lai_chabrier_ho_2004_polarization.tex:1099):

> `The maximum errors`

| Test | Published approximate error | Computed absolute relative error |
|---|---:|---:|
| $\hat a$ fit at $b=0.07$ | $1.1\%$ | $1.11745\%$ |
| $q$ fit at $b=0.4$ | $2.3\%$ | $2.30001\%$ |
| $m$ fit at $b=0.3$ | $4.2\%$ | $4.21717\%$ |

These reproduce the stated accuracy to the precision with which the percentages were printed.

For Kohri–Yamada verification, use their refractive-index equations, not an assumed equality of notation. Comparing [refs/kohri_yamada_2002_vacuum.tex:329](/home/tobias/Projects/neutronstars/refs/kohri_yamada_2002_vacuum.tex:329), anchor `\mu^2_1 =`, with the vacuum tensors gives, in the low-energy limit,

$$
\chi_0=\hat a,\qquad \chi_1=-m,\qquad \chi_2=q.
$$

Their Fig. 2 caption at [refs/kohri_yamada_2002_vacuum.tex:867](/home/tobias/Projects/neutronstars/refs/kohri_yamada_2002_vacuum.tex:867) prints `$\chi_1$ reaches the limit value` and the approximate value $7.7\times10^{-4}$. The derived asymptote is

$$
-m\longrightarrow\frac{\alpha_f}{3\pi}
=7.74273155259\times10^{-4},
\qquad
q/b\longrightarrow\frac{\alpha_f}{3\pi}.
$$

A digitised figure comparison remains to be done; it is not replaced by this asymptotic check.

## 3. The resonance is a cancellation in the wave operator

### 3.1 Plasma plus vacuum

**SOURCED.** The cold-plasma tensor appears at [refs/ho_lai_2003_atmospheres_vacuum.tex:389](/home/tobias/Projects/neutronstars/refs/ho_lai_2003_atmospheres_vacuum.tex:389), anchor `For a cold electron-ion plasma`:

$$
\epsilon_{\rm p}=
\begin{pmatrix}
\varepsilon&ig&0\\
-ig&\varepsilon&0\\
0&0&\eta
\end{pmatrix}.
$$

With negligible damping, $u_e=(\omega_{Be}/\omega)^2\gg1$, and $v_e=(\omega_{pe}/\omega)^2\ll1$, their equations give approximately

$$
\eta\simeq1-v_e,\qquad
\varepsilon\simeq1,\qquad
g\simeq-\frac{v_e}{\sqrt{u_e}(1-u_i)}.
$$

The vacuum corrections add $\hat a$ to the diagonal, $q$ along $\mathbf B$, and $m$ to the inverse permeability along $\mathbf B$.

**DERIVED HERE.** Project onto the transverse wave plane, with one axis in the $\mathbf k$–$\mathbf B$ plane. To leading order in the small susceptibilities, the difference between the two diagonal wave-operator elements is

$$
D=(\eta-\varepsilon+q+m)\sin^2\theta_{kB}
\simeq(q+m-v_e)\sin^2\theta_{kB}.
$$

The off-diagonal coupling is proportional to $g\cos\theta_{kB}$. Therefore, for an isolated crossing at a nondegenerate angle,

$$
\boxed{v_e=q+m.}
$$

This is not cancellation of every component of the dielectric tensor. It is cancellation of the diagonal splitting relevant to the transverse modes.

The local resonance condition and its corrections are at [refs/ho_lai_2003_atmospheres_vacuum.tex:708](/home/tobias/Projects/neutronstars/refs/ho_lai_2003_atmospheres_vacuum.tex:708), anchor `\vel \approx \lp\vpq+\vpm\rp`. The printed corrections include finite $u_e^{-1}$ and vacuum-coefficient terms. A production root finder should use the full wave operator, rather than assume the leading cancellation is exact.

### 3.2 Density formula and the mass convention

For the comparison convention $n_e=Y_e\rho/m_p$,

$$
v_e=\frac{4\pi e^2Y_e\rho}{m_pm_e\omega^2}.
$$

Thus **DERIVED HERE**

$$
\boxed{
\rho_V=\frac{m_pm_e}{4\pi e^2Y_e}
\left(\frac{E}{\hbar}\right)^2(q+m).
}
$$

Define

$$
f_B^2=\frac{3\delta}{q+m}.
$$

This definition is printed at [refs/ho_lai_2003_atmospheres_vacuum.tex:733](/home/tobias/Projects/neutronstars/refs/ho_lai_2003_atmospheres_vacuum.tex:733), anchor `f(B) = \frac{\Evp}{\Evpo}`. Substitution gives

$$
\rho_V
=
0.9641435534\,
Y_e^{-1}E_1^2B_{14}^2f_B^{-2}
\ {\rm g\,cm^{-3}}.
$$

The local published formula is at [refs/ho_lai_2003_atmospheres_vacuum.tex:761](/home/tobias/Projects/neutronstars/refs/ho_lai_2003_atmospheres_vacuum.tex:761):

> `\densvp = 0.964\,Y_e^{-1}(B_{14}E_1)^2 f(B)^{-2}`

The calculated coefficient differs from $0.964$ by $1.44\times10^{-4}$, within its half-last-digit tolerance $5\times10^{-4}$.

The mass convention matters: using $m_u$ instead of $m_p$ gives $0.9571786748$. The production calculation should obtain $n_e$ from its EOS and composition, avoiding a hidden choice inside a decimal prefactor.

For $Y_e=1$:

| Local energy and field | Published formula with $f_B=1$ | Calculated $f_B$ | Calculated $\rho_V$ using exact one-loop $q+m$ |
|---|---:|---:|---:|
| $E=1$ keV, $B=10^{13}$ G | $9.6400\times10^{-3}$ | $0.9798393338$ | **$1.00422715986\times10^{-2}\ {\rm g\,cm^{-3}}$** |
| $E=0.3$ keV, $B=4\times10^{12}$ G | $1.38816\times10^{-4}$ | $0.9959245811$ | **$1.39975262504\times10^{-4}\ {\rm g\,cm^{-3}}$** |

These are leading cold-plasma resonance densities with exact one-loop vacuum coefficients, not roots of a complete partially ionised, dissipative tensor.

Replacing $f_B$ by unity underestimates these densities by about $3.99\%$ and $0.813\%$, respectively. That approximation is usable as a labelled low-field stub. It is not valid over the entire field interval claimed in v1: at $10^{15}$ G, the calculation gives $f_B=2.25560449$, so setting it to unity overestimates the resonance density by a factor $5.08775$.

## 4. Derivation of the jump probability

### 4.1 The coupled amplitude equations are already local

**SOURCED.** Lai–Ho expand the wave in local eigenmodes at [refs/lai_ho_2002_vacuum_resonance.tex:257](/home/tobias/Projects/neutronstars/refs/lai_ho_2002_vacuum_resonance.tex:257), anchor `\bE(z)=A_+\bE_+e^{i\phi_+}+A_-\bE_- e^{i\phi_-}`. After the WKB reduction and removal of a common phase, their equations at lines 282–295 are

$$
i\frac{d}{ds}
\begin{pmatrix}\mathcal A_+\\\mathcal A_-\end{pmatrix}
=
\begin{pmatrix}
-\Delta k/2&i\theta_m'\\
-i\theta_m'&\Delta k/2
\end{pmatrix}
\begin{pmatrix}\mathcal A_+\\\mathcal A_-\end{pmatrix},
$$

$$
\gamma=
\left|\frac{\Delta k}{2\theta_m'}\right|,
\qquad
\tan2\theta_m=\beta^{-1}.
$$

Here $s$ is distance **along the ray**, and primes denote $d/ds$.

The mixing parameter and refractive-index splitting are printed at lines 187–188 and 300–305:

$$
\beta\simeq
\frac{\sqrt{u_e}(1-u_i)\sin^2\theta_{kB}}
{2\cos\theta_{kB}}
\left(1-\frac{q+m}{v_e}\right),
$$

$$
n_+-n_-\simeq
-\frac{v_e\cos\theta_{kB}}
{\sqrt{u_e}(1-u_i)\sin2\theta_m}.
$$

### 4.2 Density gradient and the adiabaticity parameter

**DERIVED HERE.** For locally constant field, direction, and electron fraction,

$$
H_\rho=\left|\frac{d\ln\rho}{ds}\right|^{-1}.
$$

At resonance, $\sin2\theta_m=1$ and differentiation of $\tan2\theta_m=\beta^{-1}$ gives

$$
2|\theta_m'|_{\rm res}=|\beta'|_{\rm res}
=
\frac{\sqrt{u_e}|1-u_i|\sin^2\theta_{kB}}
{2|\cos\theta_{kB}|H_\rho}.
$$

Combining this with the mode splitting yields

$$
\boxed{
\gamma_{\rm res}
=
\frac{2\omega H_\rho}{c}
\frac{q+m}
{u_e(1-u_i)^2\tan^2\theta_{kB}}.
}
$$

Since

$$
q+m=\frac{\alpha_f}{15\pi}
\frac{B^2}{B_Q^2f_B^2},
\qquad
u_e=\left(\frac{\hbar eB}{m_ecE}\right)^2,
$$

the explicit field strength cancels:

$$
\boxed{
\gamma_{\rm res}
=
\frac{2\alpha_f}{15\pi}
\frac{H_\rho E^3}
{\hbar c(m_ec^2)^2}
\frac{1}
{f_B^2\tan^2\theta_{kB}(1-u_i)^2}
\equiv\left(\frac{E}{E_{\rm ad}}\right)^3.
}
$$

Therefore,

$$
E_{\rm ad}=
\left[
\frac{15\pi}{2\alpha_f}
\frac{\hbar c(m_ec^2)^2}{H_\rho}
f_B^2\tan^2\theta_{kB}(1-u_i)^2
\right]^{1/3}.
$$

Two source discrepancies follow directly.

**First: the ion factor.** The early source prints $|1-u_i|^{1/3}$ at [refs/lai_ho_2002_vacuum_resonance.tex:314](/home/tobias/Projects/neutronstars/refs/lai_ho_2002_vacuum_resonance.tex:314). Its own preceding equations require **$|1-u_i|^{2/3}$**. The latter appears in [refs/ho_lai_2003_atmospheres_vacuum.tex:846](/home/tobias/Projects/neutronstars/refs/ho_lai_2003_atmospheres_vacuum.tex:846), anchor `|1-\uion|^{2/3}`, and in VAL06.

**Second: the numerical coefficient.** The constants above give

$$
\boxed{
E_{\rm ad}
=
2.5528429467
\left[f_B|\tan\theta_{kB}|\,|1-u_i|\right]^{2/3}
\left(\frac{1\ {\rm cm}}{H_\rho}\right)^{1/3}
{\rm keV}.
}
$$

HL03 and VAL06 print $2.52$. The discrepancy is $1.3033\%$ in $E_{\rm ad}$, or $3.8101\%$ in $\gamma$ when comparing the derived result with use of the published coefficient.

This passes a deliberately loose $2\%$ comparison to the approximate published coefficient; it does **not** pass a last-printed-digit comparison. No local explanation for that numerical difference was established here. The owned implementation should calculate $\gamma$ from the tensor gap and gradient and retain this discrepancy as a verification issue.

The ion energy should likewise be calculated:

$$
E_{Bi}=\frac{\hbar ZeB}{Am_pc}
=0.6304902519\,(Z/A)B_{14}\ {\rm keV}.
$$

Compare the printed $0.6305B_{14}(Z/A)$ at [refs/ho_lai_2003_atmospheres_vacuum.tex:431](/home/tobias/Projects/neutronstars/refs/ho_lai_2003_atmospheres_vacuum.tex:431). VAL06’s prose at line 412 omits $B_{14}$; copying that prose literally would be wrong.

### 4.3 Linear crossing gives the Landau–Zener exponential

**DERIVED HERE.** Transform from instantaneous eigenmodes to a fixed transverse basis and remove the trace. Near one isolated crossing, write

$$
i\frac{d\psi}{ds}=
\begin{pmatrix}
\Lambda s/2&V\\
V&-\Lambda s/2
\end{pmatrix}\psi,
\qquad \Lambda>0,
$$

where the diagonal detuning is linear and $V$ is locally constant. The eigenvalue separation is

$$
\Delta k=\sqrt{\Lambda^2s^2+4V^2},
$$

and differentiation of the mixing angle gives

$$
\gamma_{\rm res}=\frac{4V^2}{\Lambda}.
$$

Eliminating the second component gives

$$
\psi_1''+
\left(V^2+\frac{\Lambda^2s^2}{4}+\frac{i\Lambda}{2}\right)\psi_1=0.
$$

Set

$$
z=e^{i\pi/4}\sqrt{\Lambda}\,s,\qquad
\delta_{\rm LZ}=\frac{V^2}{\Lambda}.
$$

Then

$$
\frac{d^2\psi_1}{dz^2}
+
\left(\frac12-i\delta_{\rm LZ}-\frac{z^2}{4}\right)\psi_1=0,
$$

the parabolic-cylinder equation with index $\nu=-i\delta_{\rm LZ}$.

For completeness, the connection coefficient can be obtained from the integral representation

$$
D_\nu(z)=
\frac{e^{-z^2/4}}{\Gamma(-\nu)}
\int_0^\infty e^{-zt-t^2/2}t^{-\nu-1}\,dt,
$$

initially for $\Re\nu<0$ and then by analytic continuation. Evaluating this integral and its derivative at zero, and matching the two independent solutions, gives

$$
D_\nu(z)=
e^{-i\pi\nu}D_\nu(-z)
+
\frac{\sqrt{2\pi}}{\Gamma(-\nu)}
e^{-i\pi(\nu+1)/2}D_{-\nu-1}(iz).
$$

The coefficient of the continued diabatic branch has modulus

$$
|e^{-i\pi\nu}|=e^{-\pi\delta_{\rm LZ}}.
$$

Squaring it gives

$$
\boxed{
P_{\rm jump}=e^{-2\pi V^2/\Lambda}
=e^{-\pi\gamma_{\rm res}/2}
=\exp\left[-\frac{\pi}{2}\left(\frac{E}{E_{\rm ad}}\right)^3\right].
}
$$

The local target is [refs/lai_ho_2002_vacuum_resonance.tex:330](/home/tobias/Projects/neutronstars/refs/lai_ho_2002_vacuum_resonance.tex:330), quoted exactly:

> `P_{\rm jump}=e^{-\pi\gamma_{\rm res}/2}`

**Zener 1932 is not local.** The historical source must be fetched; it is not being treated as evidence already on disk.

An independent numerical integration of the two-component equation, with $\Lambda=1$, $V=\sqrt{\gamma}/2$, initialised in an eigenstate at $s=-80$ and projected at $s=80$, gave:

| $\gamma_{\rm res}$ | Expected exponential | Computed jump probability | Absolute error | Chosen tolerance |
|---:|---:|---:|---:|---:|
| $0.1$ | $0.8546359992$ | $0.8546361121$ | $1.13\times10^{-7}$ | $10^{-5}$ |
| $1$ | $0.2078795764$ | $0.2078780448$ | $1.53\times10^{-6}$ | $10^{-5}$ |
| $5$ | $0.0003882032$ | $0.0003883682$ | $1.65\times10^{-7}$ | $10^{-5}$ |

Doubling the endpoint from $40$ to $80$ reduced each discrepancy. The largest norm error was $1.21\times10^{-9}$, below the chosen $10^{-8}$ tolerance.

### 4.4 Conditions that must accompany the formula

The derivation requires:

1. **Geometrical optics:** the background varies slowly compared with a wavelength.
2. **An isolated, approximately linear crossing:** the detuning changes sign once, with asymptotic mode labels identifiable on both sides.
3. **Slowly varying coupling:** $V$, field direction, and the constitutive coefficients change little across the conversion region.
4. **Negligible absorption and scattering during coherent conversion**, or an explicit dissipative intensity-matrix calculation replacing the lossless crossing.
5. **A nondegenerate geometry and valid plasma approximation.**

The WKB condition is stated at [refs/lai_ho_2002_vacuum_resonance.tex:260](/home/tobias/Projects/neutronstars/refs/lai_ho_2002_vacuum_resonance.tex:260), anchor `We shall adopt the WKB (or geometric optics) approximation`.

The dissipative limitation is explicit at [refs/van_adelsberg_lai_2006_magnetic_atm.tex:483](/home/tobias/Projects/neutronstars/refs/van_adelsberg_lai_2006_magnetic_atm.tex:483), anchor `the modal description of radiative`. Its criterion is

$$
|\Re(n_X-n_O)|\gg|\Im(n_X+n_O)|.
$$

Near a cyclotron pole, at mode collapse, or near exactly parallel propagation, substituting into the scalar formula is not a validity check. In particular, $\theta_{kB}=0$ removes the diagonal vacuum-resonance splitting used above; “perfect conversion of X into O” is then an ill-defined interpretation of a singular limit.

For a varying ionisation fraction, the relevant gradient is initially $d\ln n_e/ds$, and ultimately the derivative of the **full detuning**. Replacing it by $d\ln\rho/ds$ without checking the EOS is another approximation.

## 5. Conversion belongs inside the transfer solution

### 5.1 What VAL06 actually does

**SOURCED.** The interface transformation is printed at [refs/van_adelsberg_lai_2006_magnetic_atm.tex:773](/home/tobias/Projects/neutronstars/refs/van_adelsberg_lai_2006_magnetic_atm.tex:773):

> `I_X & \rightarrow & P_{\rm jump}\, I_X + (1-P_{\rm jump}) I_O`

Together with the following equation,

$$
\mathbf I_{\rm out}=\mathsf M(P)\mathbf I_{\rm in},
\qquad
\mathsf M(P)=
\begin{pmatrix}
P&1-P\\
1-P&P
\end{pmatrix}.
$$

The same subsection then says, at lines 777–782,

> `the standard Feautrier procedure`

and

> `transfer equation cannot be used here, as there is no simple way`

VAL06 instead integrates the directional transfer equations, stops at the resonance, mixes the intensities, and continues. The procedure is repeated for downward rays.

Therefore, “implement the VAL06 method inside the existing Feautrier solve” is not a literal description of the paper. It is a proposed new numerical formulation.

### 5.2 A defensible Feautrier formulation

**DERIVED HERE.** Let $y$ increase inward, let $\mu>0$, and define mode vectors

$$
\mathbf u=\frac{\mathbf I^\uparrow+\mathbf I^\downarrow}{2},
\qquad
\mathbf v=\frac{\mathbf I^\uparrow-\mathbf I^\downarrow}{2}.
$$

Away from the crossing, for a direction-even source function and diagonal extinction matrix $\mathsf K$,

$$
\mu\frac{d\mathbf u}{dy}=\mathsf K\mathbf v,\qquad
\mu\frac{d\mathbf v}{dy}=\mathsf K(\mathbf u-\mathbf S),
$$

so

$$
\mu^2\frac{d}{dy}
\left(\mathsf K^{-1}\frac{d\mathbf u}{dy}\right)
=\mathsf K(\mathbf u-\mathbf S).
$$

These follow from the directional transfer equation at [refs/van_adelsberg_lai_2006_magnetic_atm.tex:517](/home/tobias/Projects/neutronstars/refs/van_adelsberg_lai_2006_magnetic_atm.tex:517), anchor `\pm\mu{\partial I^j_{\nu}`. An asymmetric source requires retaining its symmetric and antisymmetric parts.

Let subscripts $d$ and $r$ denote the denser and rarer sides of the resonance. The correct interface conditions are

$$
\boxed{
\mathbf u_r+\mathbf v_r
=\mathsf M_\uparrow(\mathbf u_d+\mathbf v_d),
}
$$

$$
\boxed{
\mathbf u_d-\mathbf v_d
=\mathsf M_\downarrow(\mathbf u_r-\mathbf v_r).
}
$$

They couple the two modes **and** the two one-sided solutions. They should be inserted directly into a block system.

Do not invert $\mathsf M$: its determinant is $2P-1$, and it is singular at $P=1/2$. That singularity represents equal incoherent populations after the interface, not a physical failure of the transfer problem.

Each energy requires its own resonance location. Each ray requires its own $\theta_{kB}$ and scale height. In a plane-parallel column,

$$
\boxed{
H_{\rho,\rm ray}=
\frac{1}{|\mu|\,|d\rho/dy|}.
}
$$

A single $P_{\rm jump}(E)$ reused for every angular ordinate is generally wrong.

The scattering source must still couple directions and modes. SPW09 explicitly says photons `can be scattered in some other direction and in both modes` at [refs/suleimanov_potekhin_werner_2009_mag_atm.tex:330](/home/tobias/Projects/neutronstars/refs/suleimanov_potekhin_werner_2009_mag_atm.tex:330). Conversion is an additional interface coupling; it does not replace that redistribution.

### 5.3 Thermal balance and Stokes parameters

The atmosphere must be reconverged with conversion present. The constraints are printed at [refs/van_adelsberg_lai_2006_magnetic_atm.tex:534](/home/tobias/Projects/neutronstars/refs/van_adelsberg_lai_2006_magnetic_atm.tex:534), anchors `\kappa^{abs}_j` and `=\sigma_{\rm sb} T_{\rm eff}^4`:

$$
\int d\nu\sum_j\kappa_j^{\rm abs}
\left(\frac{B_\nu}{2}-J_\nu^j\right)=0,
\qquad
F(y)=\sigma_{\rm SB}T_{\rm eff}^4.
$$

Conservation of $I_X+I_O$ **at the crossing** does not imply an unchanged emergent spectrum: the two outgoing modes undergo different subsequent absorption, emission, and scattering.

Moreover, the population matrix $\mathsf M$ is an incoherent description. It does not provide the entire Stokes vector near resonance. VAL06 explicitly discusses phase averaging at [refs/van_adelsberg_lai_2006_magnetic_atm.tex:942](/home/tobias/Projects/neutronstars/refs/van_adelsberg_lai_2006_magnetic_atm.tex:942), with the anchor

> `from a random distribution, and average over the results.`

For PRD’s Stokes output, a short coherent or dissipative intensity-matrix bridge is required where the modal approximation fails. Phase averaging must have a documented physical basis. VAL06’s displayed Stokes transfer equation also **omits scattering**, as stated at [refs/van_adelsberg_lai_2006_magnetic_atm.tex:912](/home/tobias/Projects/neutronstars/refs/van_adelsberg_lai_2006_magnetic_atm.tex:912):

> `ignores scattering`

It is therefore not a complete substitute for the missing general polarised-transfer source.

## 6. Where the resonance lies in a thin RX J1856 atmosphere

### 6.1 Benchmark assumptions

The following is a **MODEL ESTIMATE**, not a computed partially ionised atmosphere.

Ho et al.’s comparison model supplies $B=4\times10^{12}$ G, $y_b=1.2\ {\rm g\,cm^{-2}}$, $T_\infty=4.34\times10^5$ K, and $z_g=0.22$ in [refs/ho_2007_rxj1856_thin_atm.tex:628](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:628), through line 641. Greppable anchors include `& 1.2 & 1.2` and `4.34(3) & 4.34(2)`.

The model gravity is stated at [refs/ho_2007_rxj1856_thin_atm.tex:252](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:252):

> `$g=1.1\times 10^{14}$`

The temperature transformation is printed at line 128 as `$\tinfty=\teff/(1+\zg)$`. Hence

$$
T=529480\ {\rm K},\qquad k_BT=0.0456270562\ {\rm keV}.
$$

For an isothermal, fully ionised hydrogen estimate,

$$
P=gy,\qquad P=\frac{2\rho k_BT}{m_p},
$$

so

$$
y=H_n\rho,\qquad
H_n=\frac{2k_BT}{m_pg}=0.794643561\ {\rm cm}.
$$

The corresponding scale-height formula is local at [refs/ho_lai_2003_atmospheres_vacuum.tex:853](/home/tobias/Projects/neutronstars/refs/ho_lai_2003_atmospheres_vacuum.tex:853), anchor `H_\rho\simeq 2kT/(\mpr g\cos\theta)`.

Choose a field normal to the surface and an illustrative ray at $45^\circ$:

$$
\mu=\cos\theta_{kB}=1/\sqrt2,\qquad
H_{\rho,\rm ray}=1.123795701\ {\rm cm}.
$$

The isothermal base density is then

$$
\rho_b=y_b/H_n=1.51011102\ {\rm g\,cm^{-3}}.
$$

The actual base is hotter and the gas is partially ionised. Ho et al. explicitly print `The atomic fraction is $<10\%$ throughout the atmosphere,` at [refs/ho_2007_rxj1856_thin_atm.tex:276](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:276). A small neutral fraction does not justify neglecting bound opacity near atomic features.

### 6.2 Deriving the continuum photospheric estimate

**SOURCED.** Lai–Ho give the simplified free-free opacity at [refs/lai_ho_2002_vacuum_resonance.tex:355](/home/tobias/Projects/neutronstars/refs/lai_ho_2002_vacuum_resonance.tex:355):

> `\kappa_0=9.3\,\rho_1T_6^{-1/2}E_1^{-3}G`

Here $G=1-e^{-E/k_BT}$, and they take $\kappa_O\simeq\kappa_0$, $\kappa_X\simeq\kappa_0/u_e$, neglect scattering and proton effects, and set the Gaunt factor to unity. Their photospheric criterion is the `decoupling layer (at which the optical depth is $2/3$)` at line 348.

**DERIVED HERE.** Write $\kappa_O=K(E,T)\rho$. Since $dy=H_n\,d\rho$,

$$
\tau_O(\rho,\mu)
=\frac1\mu\int_0^\rho K\rho'H_n\,d\rho'
=\frac{KH_n\rho^2}{2\mu}.
$$

Therefore,

$$
\rho_O=\left(\frac{4\mu}{3KH_n}\right)^{1/2},
\qquad
\rho_X\simeq\rho_O\sqrt{u_e}.
$$

The result is:

| $B$; local $E$ | $\rho_V$ | Estimated $\rho_O$ | Extrapolated thick-atmosphere $\rho_X$ |
|---|---:|---:|---:|
| $4\times10^{12}$ G; $0.3$ keV | $1.400\times10^{-4}$ | $5.010\times10^{-2}$ | $7.733$ |
| $4\times10^{12}$ G; $1$ keV | $1.555\times10^{-3}$ | $0.3047$ | $14.109$ |
| $10^{13}$ G; $1$ keV | $1.004\times10^{-2}$ | $0.3047$ | $35.272$ |

Densities are in ${\rm g\,cm^{-3}}$.

At the main $0.3$ keV benchmark,

$$
y_V=1.11230\times10^{-4},\quad
y_O=0.0398111,\quad
y_X^{\rm extrapolated}=6.14511
\quad{\rm g\,cm^{-2}}.
$$

Thus:

- The resonance density is about **358 times smaller** than the O-mode photospheric density.
- The O-mode photosphere lies inside the hydrogen layer.
- The extrapolated X-mode photosphere lies **below its base**.

The third point matters. This thin layer does not possess the two atmospheric photospheres assumed by a semi-infinite-atmosphere argument. X-mode radiation can originate at the condensed boundary. Nevertheless, both the O-mode continuum-forming region and that boundary lie below the vacuum resonance.

Ho et al.’s discussion of this finite-column behaviour is local at [refs/ho_2007_rxj1856_thin_atm.tex:330](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:330), anchor `optically thin at high`, and their boundary treatment includes reflections at lines 350–353.

### 6.3 Optical depth above resonance

For the same continuum estimate:

| Case | $\tau_O$ above $\rho_V$, along the ray | Thomson-depth estimate $0.4y_V/\mu$ |
|---|---:|---:|
| $4\times10^{12}$ G, $0.3$ keV | $5.20\times10^{-6}$ | $6.29\times10^{-5}$ |
| $4\times10^{12}$ G, $1$ keV | $1.74\times10^{-5}$ | $6.99\times10^{-4}$ |
| $10^{13}$ G, $1$ keV | $7.24\times10^{-4}$ | $4.51\times10^{-3}$ |

The rounded Thomson opacity is printed at [refs/van_adelsberg_lai_2006_magnetic_atm.tex:524](/home/tobias/Projects/neutronstars/refs/van_adelsberg_lai_2006_magnetic_atm.tex:524), anchor `\kappa_T=0.4`. Here it is a conservative continuum scattering scale, not a replacement for the magnetic scattering calculation.

These estimates support an almost transparent conversion layer for the Ho-like model. They do not include bound features or establish an error bar for the full atmosphere.

The published normal-field criterion agrees qualitatively. VAL06 give

$$
B_l\simeq6.6\times10^{13}
T_6^{-1/8}E_1^{-1/4}(1-e^{-E/k_BT})^{-1/4}\ {\rm G}
$$

at [refs/van_adelsberg_lai_2006_magnetic_atm.tex:432](/home/tobias/Projects/neutronstars/refs/van_adelsberg_lai_2006_magnetic_atm.tex:432), anchor `B \go B_l`. It evaluates to $9.659\times10^{13}$ G at $0.3$ keV and $7.146\times10^{13}$ G at $1$ keV for this temperature.

RX J1856’s Ho-like field is far below those values. The magnetar case, in which the resonance lies between optically thick formation layers and strongly changes the spectrum, must not be imported into this benchmark without checking the ordering.

## 7. Audit of v1

### 7.1 What was implemented

The relevant commits are `e5a0f62`, **“Depth-resolved P_jump (Eddington-Barbier weighting)”**, and `19540ba`, **“Document Avrett-Krook attempt and plan (D12)”**.

The surviving code does the following:

1. Solves the two Feautrier mode problems separately during temperature iteration.
2. Performs a final independent solve.
3. Extracts emergent intensities.
4. Finds a density-grid bracket for $\rho_V$.
5. Calculates one jump probability per frequency.
6. Mixes the emergent intensities with
   $$
   f_{\rm swap}=(1-P_{\rm jump})e^{-\bar\tau_V/\mu}.
   $$

See [archive/v1/src/atmosphere/magnetic_atmosphere.jl:145](/home/tobias/Projects/neutronstars/archive/v1/src/atmosphere/magnetic_atmosphere.jl:145), lines 209–224, and [line 269](/home/tobias/Projects/neutronstars/archive/v1/src/atmosphere/magnetic_atmosphere.jl:269), through line 312.

### 7.2 What is wrong or incomplete

**The total emergent spectrum cannot change.** The final swap preserves $I_X+I_O$ at every frequency and angle. Therefore the code cannot represent conversion-induced changes to spectral intensity, heating, or temperature structure. This is an algebraic limitation, independent of grid resolution.

**The exponential is not a general fraction of radiation emitted below resonance.** Integrating the transfer equation gives, for a purely absorbing semi-infinite atmosphere,

$$
I(0,\mu)=\int_0^\infty S(t)e^{-t/\mu}\frac{dt}{\mu}.
$$

Thus the actual below-resonance fraction is

$$
\frac{\int_{\tau_V}^\infty S(t)e^{-t/\mu}\,dt}
{\int_0^\infty S(t)e^{-t/\mu}\,dt}.
$$

It equals $e^{-\tau_V/\mu}$ only for constant $S$. Conversion also changes which mode’s extinction applies to the remainder of the path. A mode-averaged optical depth does not repair that.

A simple **DERIVED HERE** counterexample makes the defect measurable. Let a prescribed unit X-mode beam reach the crossing, with no emission in the upper layer. Choose upper optical depths $\tau_X=0.1$ and $\tau_O=10$:

$$
I_{\rm emergent}^{\rm no\ conversion}=e^{-0.1}=0.9048374180,
$$

$$
I_{\rm emergent}^{\rm complete\ conversion}=e^{-10}
=4.539992976\times10^{-5}.
$$

An emergent-intensity swap preserves the first total and cannot produce the second. Conservation at an interface is not conservation through a mode-dependent absorbing slab.

Other defects are visible in the implementation:

- **Wrong ray scale height:** `H_ρ = 1.0 / dρ_dy` omits $1/|\mu|$; the resulting single probability is reused across angular ordinates. See [dielectric_tensor.jl:508](/home/tobias/Projects/neutronstars/archive/v1/src/opacity/dielectric_tensor.jl:508).
- **Approximate resonance location:** quantities are evaluated at the first grid point beyond the crossing, and optical depth is integrated to that point rather than the interpolated resonance.
- **Missing-domain ambiguity:** an out-of-grid resonance becomes `NaN` and no swap. A resonance above a truncated numerical top is not equivalent to no physical crossing.
- **Invalid-input semantics:** a nonpositive gradient can return $P=0$, meaning complete conversion, while the docstring describes this as absence of a crossing. See [dielectric_tensor.jl:479](/home/tobias/Projects/neutronstars/archive/v1/src/opacity/dielectric_tensor.jl:479) and line 504.
- **Singular limits are asserted as physical answers:** exactly parallel propagation and exact ion-cyclotron coincidence return perfect conversion without establishing applicability.
- **The advertised broad validity of $f_B=1$ is false:** the high-field calculation in Section 3 disproves the comment at [dielectric_tensor.jl:414](/home/tobias/Projects/neutronstars/archive/v1/src/opacity/dielectric_tensor.jl:414).
- **Population conservation is called “unitarity” in the notes.** It tests a stochastic intensity map, not the coherent amplitude evolution.

### 7.3 The Avrett–Krook attempt

The D12 note explicitly records that no production implementation was completed: [archive/v1/notes/decisions.md:192](/home/tobias/Projects/neutronstars/archive/v1/notes/decisions.md:192), and line 252.

Its source diagnosis is supported by [refs/suleimanov_potekhin_werner_2009_mag_atm.tex:586](/home/tobias/Projects/neutronstars/refs/suleimanov_potekhin_werner_2009_mag_atm.tex:586):

> `second procedure is the Avrett-Krook flux correction`

The paper names the procedure and directs the reader to Kurucz for details; it does not supply the missing update formula there.

The surviving v1 solver instead adds a global grey flux scaling, explicitly labelled heuristic at [archive/v1/src/atmosphere/magnetic_atmosphere.jl:160](/home/tobias/Projects/neutronstars/archive/v1/src/atmosphere/magnetic_atmosphere.jl:160).

The central issue for this path is broader than the choice of temperature-correction algorithm: **conversion is absent from the operator whose radiative equilibrium is being solved.** Installing Avrett–Krook alone would not fix that.

### 7.4 What to salvage

Salvage the source anchors, coefficient-fit fixtures, dimensional checks, resonance bracketing concept, and stochastic-map conservation tests.

Replace the post-processing weighting, angularly shared probability, invalid-input defaults, and NaN diagnostics. Preserve v1 only as an archived comparison implementation; importing it at runtime would neither satisfy ownership nor establish correct transfer.

## 8. Expected effects on spectrum and polarisation

### 8.1 Local mode contrast

**DERIVED HERE.** For an incoherent pair define

$$
I=I_X+I_O,\qquad D=I_X-I_O.
$$

The interface gives

$$
I'=I,\qquad D'=(2P_{\rm jump}-1)D.
$$

Thus conversion can leave intensity unchanged while reducing, reversing, or nearly swapping the mode contrast.

For the benchmark $B=4\times10^{12}$ G, $T=529480$ K, $\mu=\cos\theta_{kB}=1/\sqrt2$, using the derived coefficient:

| Local energy, keV | Conversion probability $1-P$ | Mode-contrast multiplier $2P-1$ |
|---:|---:|---:|
| $0.1$ | $1.220\times10^{-4}$ | $0.999756$ |
| $0.3$ | $0.0029253$ | $0.994149$ |
| $0.5$ | $0.0133503$ | $0.973299$ |
| $1$ | $0.101574$ | $0.796852$ |
| $2$ | $0.575172$ | $-0.150345$ |
| $3$ | $0.944353$ | $-0.888707$ |

These are local-ray results, not disk-integrated observed polarisation fractions. Field geometry, phase averaging, and exterior propagation remain necessary.

The source’s distinction between intensity and polarisation is explicit at [refs/van_adelsberg_lai_2006_magnetic_atm.tex:441](/home/tobias/Projects/neutronstars/refs/van_adelsberg_lai_2006_magnetic_atm.tex:441), where the resonance `lies outside both photospheres` while polarisation remains affected.

### 8.2 Optical and UV conversion

Below the ion cyclotron energy, $|1-u_i|\simeq E_{Bi}^2/E^2$. Therefore the derived expression gives

$$
\gamma_{\rm res}\propto E^7
$$

at fixed field, geometry, and scale height. Dropping the ion factor would incorrectly replace this by $E^3$.

For the benchmark,

$$
E_{Bi}=0.0252196101\ {\rm keV},
\qquad
E_{Bi,\infty}=0.0206718115\ {\rm keV}.
$$

Using representative observed wavelengths and $E_{\rm surface}=(1+z_g)hc/\lambda_\infty$:

| Observed wavelength | Surface energy | $\rho_V$, ${\rm g\,cm^{-3}}$ | $1-P$ at the illustrative angle |
|---|---:|---:|---:|
| $6000$ Å | $0.002521012$ keV | $9.88459\times10^{-9}$ | $1.74613\times10^{-13}$ |
| $1500$ Å | $0.010084048$ keV | $1.58153\times10^{-7}$ | $3.97274\times10^{-9}$ |

These small probabilities are not uniform angular bounds. Nearly field-parallel rays can be much more adiabatic. For a normal field, the formal half-conversion cone estimates are approximately $2.42\times10^{-5}$ degrees and $0.00365$ degrees for those two wavelengths. In that limit, however, the assumptions defining an isolated X/O crossing must themselves be checked.

The correct conclusion is that ordinary optical/UV rays in this benchmark are strongly nonadiabatic—not that the angular endpoints may be silently assigned a universal answer.

### 8.3 Intensity and thermal feedback

For a transparent conversion region followed by a thin absorbing layer, expansion of the formal solution gives the rough scale

$$
\left|\frac{\Delta I}{I}\right|
\sim(1-P)|\tau_O-\tau_X|
$$

when the incident modes have substantial contrast. Using the continuum optical-depth scales above gives illustrative effects of about

$$
2.0\times10^{-7}\quad\text{at }0.3\ {\rm keV},
\qquad
7.3\times10^{-5}\quad\text{at }1\ {\rm keV}
$$

for the Ho-like field and ray.

These are estimates for the direct propagation effect, not validated spectral errors.

The highest conversion probabilities occur on a weak thermal tail in this model. Ho et al. describe the relevant base temperature as `$\tatm\approx 7\times 10^5$~K` at [refs/ho_2007_rxj1856_thin_atm.tex:464](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:464). A Planck-tail integration gives bolometric fractions above $1$ keV of

$$
5.65\times10^{-7}\quad(T=529480\ {\rm K}),
\qquad
5.33\times10^{-5}\quad(T=700000\ {\rm K}).
$$

These numbers explain why a large polarisation change at high energy need not cause a large change in optical thermal emission. They are not a substitute for reconverging the actual non-grey atmosphere.

The project plan should therefore distinguish two priorities:

- Vacuum conversion is a real requirement for the polarised spectral pipeline.
- For an RX J1856 R5 intensity fit, partially ionised opacity, the finite hydrogen column, the condensed boundary, and the fitted surface parameters are likely more consequential than this conversion correction.

## 9. Concrete implementation plan

The workspace was read-only in this session. The argument and numerical checks are delivered here; `HANDOFF.md` and the worklog were not modified.

### (a) Functions, inputs, outputs, and verification

Implement through red–green TDD, starting with tests that fail under the v1 treatment.

| Function | Inputs | Outputs and required behaviour |
|---|---|---|
| `vacuum_coefficients(B, constants)` | Field and fundamental constants | $\hat a,q,m$, numerical error estimate, low-energy validity metadata. Use the exact action or closed form with a controlled asymptotic branch. |
| `transverse_wave_operator(E, state, Bvec, khat)` | EOS state, damping, composition, geometry | Flux-normalised transverse operator, eigenmodes, splitting, and applicability diagnostics. Eliminate the longitudinal field consistently. |
| `locate_vacuum_crossings(E, ray, column)` | Energy, ray, atmospheric profiles | All isolated crossings and brackets, or explicit reasons for absence. Distinguish above-grid, below-base, and invalid constitutive states. |
| `crossing_parameters(operator, derivative)` | Local operator and derivative along the ray | Gap, detuning slope, $\gamma$, conversion width, and validity result. |
| `crossing_scattering_matrix(...)` | Valid isolated crossing | Coherent amplitude map and, when justified, its phase-averaged population map. Return $\log P$ as well as $P$. |
| `solve_two_mode_transfer(...)` | Column, angular/frequency grids, opacities, lower boundary | Internal and emergent intensities, depth-dependent flux and heating, with conversion in both ray directions. |
| `solve_atmosphere_with_conversion(...)` | Surface parameters and trial state | Reconverged structure plus equilibrium, flux, discretisation, and Stokes diagnostics. |

Use directional formal solutions first, because they directly implement VAL06’s procedure. A Feautrier block-interface implementation is an alternative to verify against that solver.

Required tests and numerical targets:

1. **Weak field:** at $b=10^{-4}$, recover the three ratios in Section 2 within $2\times10^{-8}$ of unity.
2. **Normalisation:** recover $C_0=0.03291978776184018$ from both closed-form and proper-time routes. Keep a regression exposing the discrepancy with printed $0.0329199$.
3. **Published fits:** reproduce the computed $1.11745\%$, $2.30001\%$, and $4.21717\%$ errors at their stated field points.
4. **Resonance:** recover $0.9641435534$ under the stated mass convention and the two requested densities. Separately converge full-tensor roots and quantify departure from the leading formula.
5. **Adiabaticity:** retain $2.5528429467$ as the derived constant-based result and $2.52$ as the approximate source target; do not hide their discrepancy.
6. **Amplitude evolution:** reproduce the three Landau–Zener probabilities above, including norm conservation and endpoint convergence.
7. **Transfer interface:** verify $P=0$, $P=1$, and $P=1/2$; positivity; conservation at the interface; and invariance of equal LTE mode populations.
8. **Physical placement:** the absorbing-slab test must change the total from $0.9048374180$ to $4.539992976\times10^{-5}$ under complete conversion. This explicitly rejects post-processing.
9. **Published transfer comparison:** reproduce VAL06’s comparison at $B=10^{14}$ G, $E=1$ keV, $\theta_{kB}=\pi/4$, and $T=5\times10^6$ K, specified at [refs/van_adelsberg_lai_2006_magnetic_atm.tex:974](/home/tobias/Projects/neutronstars/refs/van_adelsberg_lai_2006_magnetic_atm.tex:974). Its claim `the curves agree exactly except near the resonance` at line 988 is qualitative; numerical tolerances require figure digitisation or an independently specified benchmark.
10. **Atmosphere closure:** check radiative equilibrium and flux at every depth, then frequency, angle, resonance-position, and column convergence. Ho et al.’s local target is `deviations from constant flux` below $1\%$ at [refs/ho_2007_rxj1856_thin_atm.tex:238](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:238), through line 242.

Mutation tests should change the sign of $m$, remove $\mu$ from the ray gradient, exchange $P$ with $1-P$, remove the squared ion factor, and move conversion to the emergent boundary. Fuzzing should concentrate on weak fields, nearly degenerate angles, resonance-grid boundaries, and overlapping plasma features.

No pipeline node becomes owned merely because these analytic formulas have been documented.

### (b) Missing sources to fetch

The following are **UNVERIFIED LOCALLY** unless stated otherwise; the full sources are absent from the inspected local reference set.

- **Zener, C. (1932), Proceedings of the Royal Society A, 137, 696.** Historical nonadiabatic crossing derivation. No arXiv identifier expected. **Not local; bibliographic details unverified here.**
- **Lai, D. & Ho, W. C. G. (2003), ApJ, 588, 962; arXiv:astro-ph/0211315.** Needed for the complete intensity-matrix transfer, scattering terms, and mode-collapse treatment. The venue is cited locally at [refs/potekhin_lai_chabrier_ho_2004_polarization.tex:1179](/home/tobias/Projects/neutronstars/refs/potekhin_lai_chabrier_ho_2004_polarization.tex:1179), anchor `ApJ, 588, 962`. The [arXiv record](https://arxiv.org/abs/astro-ph/0211315) was checked for discovery metadata only; its physics was not used as local evidence.
- **Kurucz, R. L. (1970), SAO Special Report, 309.** Needed if the project specifically implements the Avrett–Krook procedure cited by SPW09. Local bibliography: [refs/suleimanov_potekhin_werner_2009_mag_atm.tex:1032](/home/tobias/Projects/neutronstars/refs/suleimanov_potekhin_werner_2009_mag_atm.tex:1032), anchor `SAO Special Report, 309`. **Contents unverified; no arXiv identifier known.**
- **Mihalas, D. (1978), *Stellar Atmospheres*, second edition, W. H. Freeman.** Transfer discretisation and temperature-correction reference. **Edition, publisher, and applicable equation locators unverified here.** Do not promote the v1 note’s proposed equation number into evidence.
- **Haxton, W. C. (1995), Annual Review of Astronomy and Astrophysics, 33, 459.** Optional supporting reference for the crossing connection used by LH02. Cited locally at [refs/lai_ho_2002_vacuum_resonance.tex:488](/home/tobias/Projects/neutronstars/refs/lai_ho_2002_vacuum_resonance.tex:488), anchor `Haxton, W.C. 1995, ARAA, 33, 459`. **Contents unverified; arXiv identifier not established here.**

Heyl–Hernquist, Kohri–Yamada, LH02, and HL03 are already local. Audit statements that list them as missing should be corrected rather than used to delay implementation.

### (c) Stub policy, final-render effect, and R5

**Defensible temporary stub:** omit atmospheric conversion for an explicitly labelled Ho-like, low-field **intensity** preview, provided the code reports resonance locations, checks optical depths through the crossing and overlying layer, and retains vacuum corrections in the local opacities. The continuum estimates above suggest a sub-$10^{-3}$ optical/UV flux effect is plausible. That is a target to verify, not an established uncertainty bound.

For a resonance demonstrably outside all interacting layers, a single external conversion map is also a defensible approximation for phase-averaged mode populations. Its conditions must be checked per energy and ray.

**Indefensible stub:** use the v1 exponential weighting and call the result an in-transfer solution; claim a validated Stokes spectrum from two populations alone; force scalar Landau–Zener answers at cyclotron poles or degenerate angles; or retain the same approximation when the resonance lies in an optically thick formation layer.

For the **final render**, the Ho-like optical intensity is expected to change little from this path alone. X-ray mode contrast can change substantially: the illustrative $1$ keV ray loses about one fifth of its initial contrast, while higher-energy contrast can reverse sign. Each camera ray must apply these criteria at its emitted energy. These local estimates do not determine the disk-integrated observed Stokes vector.

For **R5**, the path is probably a small correction in the Ho-like model. A $10^{-3}$ flux change corresponds to $0.001085$ mag. For scale, Kaplan et al.’s HST power-law fit gives $F_\lambda(2500\,{\rm \AA})=4.25\pm0.07$ in its tabulated units—a $1.647\%$ uncertainty—at [refs/kaplan_2011_rxj1856_spectrum.tex:1035](/home/tobias/Projects/neutronstars/refs/kaplan_2011_rxj1856_spectrum.tex:1035), anchor `$4.25\pm0.07$`. This is a fitted flux uncertainty, not a substitute for the individual-filter acceptance errors.

The minimum honest R5 decision is therefore to reconverge and refit the atmosphere with conversion enabled, compare the actual HST band integrals against the labelled stub, and retain the stub only if the difference lies within an explicitly allocated error budget. Passing R5 with a simplified intensity model would not certify this polarised-transfer path as complete.