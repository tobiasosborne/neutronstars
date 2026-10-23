**Abstract.** The zero-field Gaunt factor can be recomputed from the local sources without consuming a published opacity table. Its numerical difficulty is cancellation in the Coulomb-wave expression, followed by a thermal integral whose endpoint limits must be handled explicitly. Two distinctions are essential: the Born limit is not a unit Gaunt factor, and an exact zero-field calculation does not make PC03’s magnetic non-Born correction exact. Python calculations below reproduce three van Hoof thermal averages within their printed rounding intervals, expose a sign error in Hummer’s published fitting coordinate, and demonstrate failure at 53-bit precision. For representative RX J1856 parameters, replacing the exact factor by its Born limit changes free-free opacity substantially in the X-ray band, while replacing it by unity also produces a large optical error. Converting those opacity errors into photometric errors requires an atmosphere calculation.

# Exact Coulomb free-free Gaunt factors: the computation the project actually needs

## 1. Scope and normalization

“Exact” here means exact within the nonrelativistic, electric-dipole, isolated Coulomb-pair problem. It does not include magnetic quantization, partial ionization, plasma screening, or relativistic corrections.

The local van Hoof source explicitly considers “an unbound electron” absorbing a photon in the field of a “point-like charge,” and states that the theory is nonrelativistic: [refs/vanhoof_2014_gaunt.tex:115](/home/tobias/Projects/neutronstars/refs/vanhoof_2014_gaunt.tex:115), lines 115–127. These restrictions belong in the function’s physical contract.

Use the **absorption convention** throughout:

$$
E_f=E_i+h\nu,\qquad E_i>0,\qquad h\nu>0.
$$

The dimensionless variables printed in [refs/vanhoof_2014_gaunt.tex:137](/home/tobias/Projects/neutronstars/refs/vanhoof_2014_gaunt.tex:137), lines 137–171, are

$$
\epsilon_i=\frac{E_i}{Z^2{\rm Ry}},\qquad
\epsilon_f=\epsilon_i+w,\qquad
w=\frac{h\nu}{Z^2{\rm Ry}},
$$

$$
\eta_i=\epsilon_i^{-1/2},\qquad
\eta_f=(\epsilon_i+w)^{-1/2},
\qquad
\eta_i>\eta_f>0,
$$

$$
k_i=\eta_i^{-1},\qquad k_f=\eta_f^{-1},
$$

$$
x_H=-\frac{4\eta_i\eta_f}{(\eta_i-\eta_f)^2},
\qquad
\beta=\frac{\eta_i+\eta_f}{\eta_i-\eta_f}.
$$

Here $k_i,k_f$ are van Hoof’s dimensionless variables, not dimensional wave numbers. I write $x_H$ instead of the paper’s $x$ to distinguish it from the later thermal integration variable. A greppable defining anchor is `\eta_i = \frac{1}{\sqrt{\epsilon_i}}`. The energy unit is the infinite-nuclear-mass Rydberg; the same source prints `1~{\rm Ry} = \alpha^2 m_{\rm e} c^2 / 2` at line 142.

For fully ionized hydrogen, the normalization entering the atmosphere is McPHAC’s equation 12:

$$
\kappa_\nu^{\rm ff}
=
\underbrace{
\frac{4e^6}{3m_ehc}
\sqrt{\frac{2\pi}{3k_Bm_e}}\,
T^{-1/2}\frac{\rho}{(m_p+m_e)^2}
\nu^{-3}(1-e^{-u})
}_{\kappa_{\nu,K}}
\overline g_{\rm ff}(u,\gamma^2),
$$

$$
u=\frac{h\nu}{k_BT},\qquad
\gamma^2=\frac{Z^2{\rm Ry}}{k_BT}.
$$

The source is [refs/haakonsen_2012_mcphac.txt:233](/home/tobias/Projects/neutronstars/refs/haakonsen_2012_mcphac.txt:233), lines 233–248, anchored by “opacity is the free–free opacity (Rybicki & Lightman 1979),”. I checked the radical and mass denominator against PDF page 4.

Thus **Kramers opacity means setting $\overline g_{\rm ff}=1$**. It is a normalization, not the general Born result.

### Consistency with PC03 requires weak damping

PC03 defines

$$
\Lambda=\frac{\pi}{\sqrt3}\overline g
$$

and

$$
\nu_\alpha^{\rm ff}
=
\frac43\sqrt{\frac{2\pi}{m_ek_BT}}\,
\frac{n_e e^4}{\hbar\omega}
(1-e^{-u})\Lambda_\alpha^{\rm ff}.
$$

See [refs/potekhin_chabrier_2003_ff_opacity.tex:1390](/home/tobias/Projects/neutronstars/refs/potekhin_chabrier_2003_ff_opacity.tex:1390), lines 1390–1409; the exact anchor is `Coulomb logarithms $\Lambda=(\pi/\sqrt{3})\,\bar{g}$.`

Its fixed-ion absorption cross section contains

$$
\sigma_\alpha^{\rm ff}
=
\frac{4\pi e^2}{m_ec}
\frac{\nu_\alpha^{\rm ff}}
{(\omega+\alpha\omega_{ce})^2+\nu_{e,\alpha}^2},
$$

from [refs/potekhin_chabrier_2003_ff_opacity.tex:1344](/home/tobias/Projects/neutronstars/refs/potekhin_chabrier_2003_ff_opacity.tex:1344), anchored by “the free-free absorption by electrons is”.

**DERIVED HERE.** At $B=0$, substitute $\Lambda=(\pi/\sqrt3)\overline g$, $\omega=2\pi\nu$, $\hbar=h/(2\pi)$, and $n_e/\rho=1/(m_p+m_e)$. The result is

$$
\frac{n_e\sigma^{\rm ff}}{\rho}
=
\kappa_{\nu,K}\overline g_{\rm ff}
\frac{\omega^2}{\omega^2+\nu_e^2}.
$$

Consequently, the undamped Kramers expression follows only for $\nu_e/\omega\ll1$. The rejected database entry must be repaired by adding that qualification, not by declaring the two expressions unconditionally identical.

## 2. The exact Sommerfeld–Biedenharn expression

The locally printed expression is

$$
g_{\rm ff}=
\frac{2\sqrt3}{\pi\eta_i\eta_f}
\left[
(\eta_i^2+\eta_f^2+2\eta_i^2\eta_f^2)I_0
-
2\eta_i\eta_f
\sqrt{1+\eta_i^2}\sqrt{1+\eta_f^2}\,I_1
\right]I_0,
$$

where, for $l=0,1$,

$$
I_l=
\frac14(-x_H)^{l+1}
e^{\pi(\eta_i-\eta_f)/2}
\frac{
|\Gamma(l+1+i\eta_i)\Gamma(l+1+i\eta_f)|
}{
\Gamma(2l+2)
}
G_l,
$$

$$
G_l=
\beta^{-i(\eta_i+\eta_f)}
{}_2F_1
\left(l+1-i\eta_f,l+1-i\eta_i;2l+2;x_H\right).
$$

Source: [refs/vanhoof_2014_gaunt.tex:178](/home/tobias/Projects/neutronstars/refs/vanhoof_2014_gaunt.tex:178), lines 178–198; anchor: “The free-free Gaunt factor is given by Eq.” The source identifies the expression as KL61 equation 16, based on Biedenharn.

Although written using complex functions, $G_l$ is real for these physical arguments. Numerically, an imaginary residual is an error diagnostic. Taking an absolute value of $G_l$, or silently discarding an uncontrolled imaginary part, changes the calculation.

This closed form is an evaluation of Coulomb continuum matrix elements, not a fitted opacity. KL61 describes the underlying electrons as “normalized to plane waves of unit density at infinity” and states that the partial-wave sum “can be explicitly evaluated”: [refs/karzas_latter_1961_gaunt.txt:54](/home/tobias/Projects/neutronstars/refs/karzas_latter_1961_gaunt.txt:54), lines 54–55 and 100–102. The badly extracted KL61 equations need not be retranscribed when van Hoof supplies the usable expression.

### 2.1 Born limit: the factor is logarithmic, not unity

**DERIVED HERE from the preceding exact expression.**

Take $\eta_i,\eta_f\rightarrow0$ at fixed ratio, so $x_H$ and $\beta$ remain fixed. Write $L=\ln\beta$. The hypergeometric functions reduce to

$$
{}_2F_1(1,1;2;x)
=-\frac{\ln(1-x)}{x},
$$

$$
{}_2F_1(2,2;4;x)
=
6\,\frac{(x-2)\ln(1-x)-2x}{x^3}.
$$

These identities follow by integrating the geometric series once, and by evaluating
$6\int_0^1 t(1-t)(1-xt)^{-2}\,dt$, respectively. Since $1-x_H=\beta^2$, substitution gives

$$
I_0\longrightarrow\frac L2,
\qquad
I_1\longrightarrow
\frac{\eta_i^2+\eta_f^2}{4\eta_i\eta_f}L-\frac12.
$$

Therefore

$$
(\eta_i^2+\eta_f^2)I_0-2\eta_i\eta_f I_1
=\eta_i\eta_f,
$$

and

$$
\boxed{
g_{\rm Born}
=
\frac{\sqrt3}{\pi}
\ln\frac{k_f+k_i}{k_f-k_i}
=
\frac{\sqrt3}{\pi}
\ln\frac{p_f+p_i}{p_f-p_i}.
}
$$

This agrees with the logarithm in PC03’s Born cross section: [refs/potekhin_chabrier_2003_ff_opacity.tex:1420](/home/tobias/Projects/neutronstars/refs/potekhin_chabrier_2003_ff_opacity.tex:1420), lines 1420–1430, including the printed anchor `\ln\left|\frac{p_f+p_i}{p_f-p_i}\right|`.

A Python/mpmath check at $\eta_i=10^{-6}$, $\eta_f=5\times10^{-7}$ gives:

| Quantity | Value |
|---|---:|
| Born target, derived above: $(\sqrt3/\pi)\ln3$ | 0.6056966996081959 |
| Exact Coulomb expression | 0.6056976510343444 |
| Relative difference | $1.5707963\times10^{-6}$ |
| Proposed finite-parameter test tolerance | $2\times10^{-6}$ |

The test approaches the Born limit; it is not supposed to equal it exactly at finite $\eta$.

The requirement to make **both** the Born and classical limits give $g=1$ would be wrong. The Born limit reproduces **Kramers normalization multiplied by the Born logarithmic correction**.

### 2.2 The classical Kramers limit requires more than large $\eta$

Define

$$
r=\frac{\eta_f^2}{\eta_i^2},\qquad
q_C=(1-r)\eta_f
=\frac{w}{(\epsilon_i+w)^{3/2}}.
$$

Van Hoof’s corrected expansion requires $q_C\gg1$:

$$
\begin{aligned}
g_{\rm ff}={}&1+
c_1\frac{1+r}{q_C^{2/3}}
-c_2\frac{1-\frac43r+r^2}{q_C^{4/3}}\\
&-c_3\frac{1-\frac13r-\frac13r^2+r^3}{q_C^2}
+R,
\end{aligned}
$$

where

$$
c_1=\frac{\Gamma(1/3)}{5\,12^{1/3}\Gamma(2/3)},\quad
c_2=\frac{18\Gamma(2/3)}{35\,12^{2/3}\Gamma(1/3)},\quad
c_3=\frac3{175}.
$$

Source: [refs/vanhoof_2014_gaunt.tex:336](/home/tobias/Projects/neutronstars/refs/vanhoof_2014_gaunt.tex:336), lines 336–351; the validity anchor is `with $(1 - \eta_f^2/\eta_i^2) \, \eta_f \gg 1$`.

Python evaluation gives

$$
(c_1,c_2,c_3)=
(0.1728260369359927,\,
0.04959570168586959,\,
0.01714285714285714).
$$

The unit leading term can be derived rather than accepted as a normalization accident.

**DERIVED HERE using Appendix A’s leading terms.** Van Hoof writes two auxiliary functions $P,Q$ whose leading terms have a common factor

$$
\mathcal K=
-\frac{\sqrt3}{2\pi}
e^{\pi\eta_f}\beta^{-i(\eta_i+\eta_f)}.
$$

Their real and imaginary amplitudes have magnitudes

$$
A=
\frac{(1-r)^{1/3}\Gamma(2/3)}
{12^{1/3}\eta_f^{2/3}},
\qquad
B=
\frac{(\beta+1)(1-r)^{2/3}\Gamma(1/3)}
{12^{2/3}\beta\eta_f^{1/3}}.
$$

Thus $|Q^2-P^2|=4|\mathcal K|^2AB$. Using the identities printed at the end of Appendix A,

$$
\Gamma(1/3)\Gamma(2/3)=\frac{2\pi}{\sqrt3},
\qquad
\frac{(1-r)(\beta+1)}{\beta}
=\frac{2(\eta_i-\eta_f)}{\eta_i},
$$

one obtains

$$
|\Delta|\sim
\frac{e^{2\pi\eta_f}(\eta_i-\eta_f)}
{\pi\sqrt3\,\eta_i\eta_f}.
$$

Appendix A’s exact representation then gives

$$
g_{\rm ff}
=
\frac{\pi\sqrt3\eta_i\eta_f e^{-2\pi\eta_f}|\Delta|}
{(\eta_i-\eta_f)(1-e^{-2\pi\eta_i})(1-e^{-2\pi\eta_f})}
\longrightarrow1.
$$

The supporting equations are at [refs/vanhoof_2014_gaunt.tex:920](/home/tobias/Projects/neutronstars/refs/vanhoof_2014_gaunt.tex:920), lines 920–965; line 950 contains the greppable factor `|\Delta|`.

This proves recovery of the Kramers factor in the specified classical limit. It also exposes an invalid implementation shortcut: **large $\eta_i,\eta_f$ alone does not imply $g=1$**. If the photon becomes sufficiently soft, $q_C$ becomes small and the expansion is inapplicable.

## 3. Thermal averaging and its limits

Van Hoof §3 prints

$$
\boxed{
\overline g_{\rm ff}(u,\gamma^2)
=
\int_0^\infty e^{-t}
g_{\rm ff}
\left(\epsilon_i=\frac{t}{\gamma^2},
w=\frac{u}{\gamma^2}\right)\,dt.
}
$$

Equivalently,

$$
\eta_i=\frac{\gamma}{\sqrt t},
\qquad
\eta_f=\frac{\gamma}{\sqrt{t+u}}.
$$

Source: [refs/vanhoof_2014_gaunt.tex:579](/home/tobias/Projects/neutronstars/refs/vanhoof_2014_gaunt.tex:579), lines 579–589; the integration anchor is `\int_0^\infty {\rm e}^{-x}`. I renamed its $x$ to $t$.

### 3.1 Why the weight is $e^{-t}$

**DERIVED HERE.** The Maxwell momentum measure contributes

$$
p_i^2 e^{-p_i^2/(2m_ek_BT)}\,dp_i.
$$

The absorption cross section carries the $1/p_i$ factor printed in PC03 equation 42. Therefore

$$
p_i^2\,dp_i\;\frac1{p_i}
=p_i\,dp_i=m_e\,dE_i.
$$

After dividing by the thermally averaged Kramers prefactor and setting $t=E_i/k_BT$, the remaining weight is precisely $e^{-t}dt$.

Consequently:

- Do not insert an additional $\sqrt t$ Maxwell weight.
- Do not use a lower limit $u$ with these absorption variables.
- Do not put an extra $e^{-u}$ in this definition.
- The stimulated-emission factor $1-e^{-u}$ remains outside $\overline g$ in the absorption coefficient.

Changing to emission-state energies is possible, but requires transforming the integration variable and normalization consistently.

### 3.2 Thermal Born result

Substitution of the Born logarithm gives

$$
\overline g_{\rm Born}
=
\frac{\sqrt3}{\pi}
\int_0^\infty e^{-t}
\ln\frac{\sqrt{t+u}+\sqrt t}
{\sqrt{t+u}-\sqrt t}\,dt.
$$

**DERIVED HERE.** Set

$$
t=\frac u2(\cosh s-1).
$$

The logarithm becomes $s$. Integration by parts yields

$$
\int_0^\infty
s\,e^{-t}\frac u2\sinh s\,ds
=
e^{u/2}\int_0^\infty e^{-(u/2)\cosh s}\,ds.
$$

The final integral is the integral representation of $K_0(u/2)$, giving

$$
\boxed{
\overline g_{\rm Born}(u)
=
\frac{\sqrt3}{\pi}e^{u/2}K_0(u/2).
}
$$

This independently recovers PC03’s printed result, anchored by `\Lambda_\mathrm{cl}^\mathrm{ff} = \mathrm{e}^{u/2}K_0(u/2)` at [refs/potekhin_chabrier_2003_ff_opacity.tex:1436](/home/tobias/Projects/neutronstars/refs/potekhin_chabrier_2003_ff_opacity.tex:1436).

Python values are

| $u$ | $\overline g_{\rm Born}$ |
|---:|---:|
| $10^{-4}$ | 5.524271458100495 |
| 1 | 0.8402855441607155 |
| 10 | 0.3020221393366747 |
| 100 | 0.09747889394154515 |

These values alone rule out treating “Born” and “Kramers with $g=1$” as interchangeable.

### 3.3 Convergence for positive $u$

For fixed $u>0$ and $\gamma^2>0$:

1. At $t\rightarrow0$, $\eta_i\rightarrow\infty$, but $g_{\rm ff}$ approaches a finite endpoint value. Van Hoof explicitly confirms the existence of this limit and prints, for small $w$,

   $$
   g_{\rm ff}(0,w)
   =
   1+c_1w^{1/3}-c_2w^{2/3}-c_3w+O(w^{4/3}).
   $$

   Source: [refs/vanhoof_2014_gaunt.tex:512](/home/tobias/Projects/neutronstars/refs/vanhoof_2014_gaunt.tex:512), lines 512–522; anchor: `g_{\rm ff}(0,w) = 1 + c_1 w^{1/3}`.

2. At $t\rightarrow\infty$, both Coulomb parameters become small, and the kernel grows only logarithmically, as $\ln(4t/u)$ at leading order. The Maxwell factor makes the integral converge.

The raw hypergeometric evaluator should not be called at $t=0$. An open quadrature or the transformation $t=s^2$ avoids that singular parameterization without altering the integral.

### 3.4 Soft photons: finite integrals, divergent Gaunt factor

**DERIVED HERE from the transformed exact expression.** Put $\delta=\eta_i-\eta_f$ and let $\delta\rightarrow0$ at fixed $\eta=\eta_i$. In the connection formula, the relevant expansions are

$$
\Gamma(-i\delta)=\frac{i}{\delta}-\gamma_E+O(\delta),
$$

$$
\arg\Gamma(l+1+i\eta_i)-\arg\Gamma(l+1+i\eta_f)
=
\delta\,\Re\psi(l+1+i\eta)+O(\delta^2).
$$

Their product with $\beta^{-i\delta}$ gives

$$
I_l=
\frac12\left[
\ln\beta-\gamma_E-\Re\psi(l+1+i\eta)
\right]+o(1).
$$

Here $\psi=\Gamma'/\Gamma$, and $\gamma_E=-\Gamma'(1)$. The recurrence
$\psi(z+1)=\psi(z)+1/z$, obtained by differentiating $\Gamma(z+1)=z\Gamma(z)$, gives

$$
I_0-I_1=\frac{1}{2(1+\eta^2)}+o(1).
$$

Substitution in the exact Gaunt expression therefore yields

$$
g_{\rm ff}
=
\frac{\sqrt3}{\pi}
\left[
\ln\frac{4t}{u}
-\gamma_E
-\Re\psi\!\left(1+i\frac{\gamma}{\sqrt t}\right)
\right]+o(1).
$$

The transformation being expanded is printed at [refs/vanhoof_2014_gaunt.tex:251](/home/tobias/Projects/neutronstars/refs/vanhoof_2014_gaunt.tex:251), lines 251–260, with anchor `G_l = 2 \, \mathfrak{Re}`.

Thermal integration gives

$$
\overline g_{\rm ff}
=
\frac{\sqrt3}{\pi}
\left[
\ln\frac4u-2\gamma_E
-
\int_0^\infty e^{-t}
\Re\psi\!\left(1+i\frac{\gamma}{\sqrt t}\right)\,dt
\right]+o(1).
$$

Thus $\overline g$ diverges logarithmically as $u\rightarrow0$; it has no finite value to clamp to at $u=0$.

In the Born limit, $\psi(1)=-\gamma_E$, so

$$
\overline g_{\rm Born}
\sim
\frac{\sqrt3}{\pi}\left[\ln\frac4u-\gamma_E\right].
$$

In the large-$\gamma$, soft-photon regime, using the leading logarithmic behavior of $\psi$ gives

$$
\overline g_{\rm ff}
\sim
\frac{\sqrt3}{\pi}
\left[\ln\frac4u-\ln\gamma-\frac52\gamma_E\right].
$$

This agrees with Hummer’s printed long-wavelength expression, [refs/hummer_1988_gaunt.txt:241](/home/tobias/Projects/neutronstars/refs/hummer_1988_gaunt.txt:241), equation 3.1, anchored by “A long-wavelength approximation given by Scheuer”.

The distinction from the Kramers limit is now explicit: at thermally typical $t$, $q_C\sim\gamma u$ for small $u$. The limits $\gamma u\gg1$ and $\gamma u\ll1$ describe different classical regimes.

A numerical check of the derived soft expression at $\gamma^2=1$, $t=1$, $u=10^{-9}$ gives

$$
g_{\rm exact}=11.81922028954181,
\qquad
g_{\rm soft}=11.81922027059109,
$$

with absolute difference $1.8951\times10^{-8}$.

The logarithmic divergence is integrable in the frequency-integrated emissivity. Combining it with the Kramers absorption prefactor instead gives $\kappa_\nu^{\rm ff}\propto\nu^{-2}\ln(1/\nu)$ at fixed temperature in this isolated-pair model. These are mathematical limits, not permission to ignore the physical domain of the plasma calculation.

### 3.5 Hard photons: do not force the non-Born ratio to unity

**DERIVED HERE.** At fixed $\eta_i=a$ and $\eta_f=b\rightarrow0$, $x_H=-4b/a+O(b^2)$ and the leading hypergeometric term gives

$$
I_0\sim\frac ba e^{\pi a/2}|\Gamma(1+ia)|.
$$

Hence

$$
g_{\rm ff}
\sim
\frac{2\sqrt3}{\pi}\frac ba
e^{\pi a}|\Gamma(1+ia)|^2.
$$

The gamma normalization follows from the beta integral:

$$
|\Gamma(1+ia)|^2
=
\int_0^\infty\frac{s^{ia}}{(1+s)^2}\,ds
=
\frac{\pi a}{\sinh\pi a}.
$$

The second equality follows by a keyhole-contour evaluation around the positive-axis branch cut; the double pole at $s=-1$ supplies the residue. Therefore

$$
g_{\rm ff}\sim
\frac{4\sqrt3\,b}{1-e^{-2\pi a}},
$$

and

$$
\boxed{
\overline g_{\rm ff}(u,\gamma^2)
\sim\frac{A(\gamma)}{\sqrt u},
\quad
A(\gamma)=
4\sqrt3\,\gamma
\int_0^\infty
\frac{e^{-t}}{1-e^{-2\pi\gamma/\sqrt t}}\,dt.
}
$$

The coefficient integral is finite: its integrand tends to a constant at the lower endpoint and is Maxwell-suppressed at the upper endpoint.

For Born scattering, Laplace expansion of the $K_0$ integral gives

$$
\overline g_{\rm Born}
\sim
\sqrt{\frac3{\pi u}}
\left(1-\frac1{4u}+\frac9{32u^2}+\cdots\right).
$$

The absorption average thus decays algebraically. The thermal emissivity additionally contains its Boltzmann factor.

Most importantly,

$$
\lim_{u\to\infty}
\frac{\overline g_{\rm ff}}{\overline g_{\rm Born}}
=
\frac{A(\gamma)}{\sqrt{3/\pi}},
$$

which need not equal one. A hard outgoing electron does not make the incoming electron weakly Coulomb-distorted. The implementation must not use “large photon energy” alone as a switch to Born scattering.

For completeness, the thermal Kramers limit follows from the classical expansion at fixed positive $u$, $\gamma\rightarrow\infty$:

$$
\overline g_{\rm ff}
=
1+
c_1\frac{u+2}{(\gamma u)^{2/3}}
-
c_2\frac{u^2+\frac23u+\frac43}{(\gamma u)^{4/3}}
+\cdots
\longrightarrow1.
$$

This is obtained by substituting $r=t/(t+u)$ and integrating the resulting polynomials with $e^{-t}$.

## 4. Numerical evaluation: an error-controlled calculation, not a special-function call

### 4.1 Transform the argument before summation

Van Hoof supplies transformations keeping the hypergeometric series argument between zero and one half: [refs/vanhoof_2014_gaunt.tex:220](/home/tobias/Projects/neutronstars/refs/vanhoof_2014_gaunt.tex:220), lines 220–273. The explanatory anchor is “The last argument of the hypergeometric function is between 0 and”.

For $-1\le x_H<0$,

$$
G_l=
\beta^{-2l-2+i(\eta_f-\eta_i)}
{}_2F_1
\left(l+1+i\eta_i,l+1-i\eta_f;2l+2;
\frac{x_H}{x_H-1}\right).
$$

For $x_H<-1$,

$$
\begin{aligned}
G_l=2\Re\Bigg[
&\frac{\Gamma(2l+2)\Gamma(i\eta_f-i\eta_i)}
{\Gamma(l+1-i\eta_i)\Gamma(l+1+i\eta_f)}
\beta^{-2l-2+i(\eta_f-\eta_i)}
\\
&\times
{}_2F_1\left(
l+1+i\eta_i,l+1-i\eta_f;
1+i\eta_i-i\eta_f;
\frac1{1-x_H}
\right)
\Bigg].
\end{aligned}
$$

From the locally printed Taylor series, the summand recurrence is **DERIVED HERE**:

$$
T_0=1,\qquad
T_{n+1}
=
T_n\,z\frac{(a+n)(b+n)}{(c+n)(n+1)}.
$$

The source series is at [refs/vanhoof_2014_gaunt.tex:198](/home/tobias/Projects/neutronstars/refs/vanhoof_2014_gaunt.tex:198), anchored by `{}_2F_1(a,b;c;x) = \sum_{n=0}^{\infty}`.

A small $z$ does not eliminate cancellation when the complex parameters are large. Track both $\sum T_n$ and $\sum|T_n|$. Do not stop merely because one term is small. Once a bound $|T_{n+1}/T_n|\le q<1$ holds for the remaining tail, its absolute remainder is bounded by the next term divided by $1-q$.

### 4.2 Stabilize the elementary operations too

**DERIVED HERE.** Compute the Coulomb-parameter difference as

$$
\eta_i-\eta_f
=
\frac{w}
{\sqrt{\epsilon_i}\sqrt{\epsilon_i+w}
\left(\sqrt{\epsilon_i}+\sqrt{\epsilon_i+w}\right)},
$$

rather than subtracting nearly equal reciprocal square roots.

Combine the gamma factors and exponential through

$$
\begin{aligned}
\log P_l={}&-\log4+(l+1)\log(-x_H)
+\frac\pi2(\eta_i-\eta_f)\\
&+\Re\log\Gamma(l+1+i\eta_i)
+\Re\log\Gamma(l+1+i\eta_f)
-\log\Gamma(2l+2).
\end{aligned}
$$

Then $I_l=e^{\log P_l}G_l$, with further scaling if necessary. Reuse
$\Gamma(z+1)=z\Gamma(z)$ between $l=0$ and $l=1$.

Python gives

$$
\log_{10}|\Gamma(1+1000i)|=-680.2890869867416.
$$

That individual factor underflows in ordinary binary64 arithmetic although the final Gaunt factor is perfectly finite. Logarithmic assembly fixes this failure, but not the subsequent subtraction in the exact expression.

For that subtraction, define

$$
A=(\eta_i^2+\eta_f^2+2\eta_i^2\eta_f^2)I_0,
$$

$$
B=2\eta_i\eta_f\sqrt{1+\eta_i^2}\sqrt{1+\eta_f^2}I_1,
$$

and monitor

$$
\mathcal C_{\rm sub}=\frac{|A|+|B|}{|A-B|}.
$$

This measures only one source of cancellation; gamma phases and hypergeometric sums need their own accounting.

### 4.3 A computed precision failure

I evaluated the exact expression with mpmath at 53-, 160-, and 240-bit working precision. The 53-bit run retained mpmath’s wide exponent range, so this isolates mantissa loss rather than native floating-point underflow. For this diagnostic only, I took real parts to expose the resulting error; production code must reject uncontrolled imaginary residuals.

| $\epsilon_i$ | $w$ | 53-bit result | 240-bit reference | Relative error |
|---:|---:|---:|---:|---:|
| $10^{-6}$ | $10^{-8}$ | 1.073108191580731 | 1.073108259394692224 | $-6.3194\times10^{-8}$ |
| $10^{-8}$ | $10^{-12}$ | 1.305360399717817 | 1.304907699405679672 | $+3.4692\times10^{-4}$ |
| $10^{-10}$ | $10^{-15}$ | 2.898957412491485 | 1.304895799335132149 | $+1.2216$ |

At the last point, $\mathcal C_{\rm sub}=2.54807315323\times10^{10}$. The 160- and 240-bit answers agree to relative $5.5\times10^{-33}$.

A plausible, positive answer is therefore not evidence of accuracy.

Van Hoof explicitly reports “severe cancellation problems” and uses adaptive precision: [refs/vanhoof_2014_gaunt.tex:276](/home/tobias/Projects/neutronstars/refs/vanhoof_2014_gaunt.tex:276), lines 276–321. That source also identifies regions where even very high precision becomes impractical.

### 4.4 Use Appendix A with an honest remainder status

For $\epsilon_i=10^{-8}$, $w=10^{-6}$, Python gives

$$
q_C=985.185336841573,
$$

$$
g_{\rm series}=1.001757815008687975.
$$

The published remainder envelope gives

$$
1.4048\times10^{-11}<R<2.6016\times10^{-10}.
$$

However, van Hoof says the envelope “was derived by comparing the series expansion with exact” calculations: [refs/vanhoof_2014_gaunt.tex:353](/home/tobias/Projects/neutronstars/refs/vanhoof_2014_gaunt.tex:353).

That is a calibrated numerical envelope, not a supplied rigorous remainder proof. Accordingly:

- It is useful as a verification target.
- It must not be presented as an interval-certified error bound.
- An owned asymptotic branch should generate additional coefficients from Appendix A and verify overlap with the arbitrary-precision branch.
- A request for tighter accuracy than the supported remainder permits must trigger more computation or an explicit failure.

### 4.5 Precision and quadrature targets

The following are **proposed engineering budgets**, not physical uncertainties:

| Component | Initial target |
|---|---:|
| Reference Coulomb kernel | $10^{-12}$ relative |
| Reference thermal average | $10^{-8}$ relative |
| Atmosphere runtime average, including interpolation | $10^{-5}$ relative |
| Precision verification | Repeat with increased precision and require agreement within the allocated budget |

Use open adaptive quadrature, preferably with $t=s^2$ and intervals resolving both $t\sim u$ and $t\sim1$. Propagate kernel error separately from quadrature error.

Van Hoof’s integration estimates the remaining tail by holding $g$ constant, explicitly described at [refs/vanhoof_2014_gaunt.tex:608](/home/tobias/Projects/neutronstars/refs/vanhoof_2014_gaunt.tex:608), lines 608–620: “by assuming that $g_{\rm ff}$ is constant.” That is a reasonable estimator, but it is not an upper bound for a growing kernel.

**DERIVED HERE.** If an asymptotic analysis establishes

$$
g(t)\le C_0+C_1\ln(1+t),\qquad t\ge L,
$$

concavity of the logarithm gives

$$
\int_L^\infty e^{-t}g(t)\,dt
\le
e^{-L}
\left[
C_0+C_1\left(\ln(1+L)+\frac1{1+L}\right)
\right].
$$

The constants must be established from the kernel, not guessed from its last sampled value.

## 5. Verification against local tables and Hummer’s fit

### 5.1 Three thermal reference points

I computed the integral in two parameterizations:

- $t=s^2$, using 45-decimal-digit mpmath kernel evaluations and adaptive quadrature;
- direct $t$ integration on different subintervals, using 60-decimal-digit kernel evaluations.

The two results agree to the displayed precision. This is a quadrature and working-precision check, not an independent derivation of the Coulomb kernel.

All three selected points have $\gamma^2=1$, which also makes the Hummer fit independently reproducible from the even-index rows of its coefficient table.

| $u$ | Local van Hoof value and exact quoted token | Computed average | Computed minus printed | Allowed absolute difference¹ |
|---:|---|---:|---:|---:|
| 0.01 | [tex:558](/home/tobias/Projects/neutronstars/refs/vanhoof_2014_gaunt.tex:558), `2.6563$+0$` | 2.656310492668394 | $+1.04927\times10^{-5}$ | $5.00797\times10^{-5}$ |
| 1 | [tex:560](/home/tobias/Projects/neutronstars/refs/vanhoof_2014_gaunt.tex:560), `1.2937$+0$` | 1.293739531917590 | $+3.95319\times10^{-5}$ | $5.00388\times10^{-5}$ |
| 10 | [tex:561](/home/tobias/Projects/neutronstars/refs/vanhoof_2014_gaunt.tex:561), `9.7260$-1$` | 0.9725974478935024 | $-2.55211\times10^{-6}$ | $5.02918\times10^{-6}$ |

¹ Half a printed rounding unit plus the paper’s stated approximate relative error of $3\times10^{-8}$. The precision statement appears at [refs/vanhoof_2014_gaunt.tex:540](/home/tobias/Projects/neutronstars/refs/vanhoof_2014_gaunt.tex:540), anchored by “All entries have an approximate relative error of”.

The source cannot support a $10^{-8}$ comparison directly against its rounded printed values. The extra digits above are computed reference candidates; a tighter external test needs the electronic data and their error estimates.

### 5.2 The McPHAC file: interpret the columns correctly

The local file’s columns must be interpreted as

$$
(\gamma^2,\ u,\ \overline g_{\rm ff}).
$$

The warning is explicit: “the column labels u and γ 2 are reversed,” in [refs/haakonsen_2012_mcphac.txt:282](/home/tobias/Projects/neutronstars/refs/haakonsen_2012_mcphac.txt:282). Do not swap the actual data twice.

At the same points:

| Local fixture | Quoted value | Relative difference from the computed integral |
|---|---:|---:|
| [refs/code/McPHAC/gffgu.dat:1641](/home/tobias/Projects/neutronstars/refs/code/McPHAC/gffgu.dat:1641) | `2.655969` | $-1.28559\times10^{-4}$ |
| [refs/code/McPHAC/gffgu.dat:1661](/home/tobias/Projects/neutronstars/refs/code/McPHAC/gffgu.dat:1661) | `1.293904` | $+1.27126\times10^{-4}$ |
| [refs/code/McPHAC/gffgu.dat:1671](/home/tobias/Projects/neutronstars/refs/code/McPHAC/gffgu.dat:1671) | `0.9725416` | $-5.74214\times10^{-5}$ |

A proposed relative tolerance of $2\times10^{-3}$ is appropriate for a broad Sutherland-table compatibility test, with discrepancies reported rather than hidden. It is not the accuracy target for the new calculation.

Van Hoof reports a maximum thermal-table discrepancy of almost $0.13\%$, after correcting the transposition: [refs/vanhoof_2014_gaunt.tex:629](/home/tobias/Projects/neutronstars/refs/vanhoof_2014_gaunt.tex:629), lines 629–639, including “amounts to almost 0.13\%”.

### 5.3 Hummer is a Chebyshev fit, and its printed coordinate has a sign error

Hummer describes “A two-dimensional Chebyshev expansion” in [refs/hummer_1988_gaunt.txt:16](/home/tobias/Projects/neutronstars/refs/hummer_1988_gaunt.txt:16). PC03 calls it a Padé formula, but that is not the fit printed in the local Hummer paper. The database should record this source discrepancy.

Hummer’s equations 2.30–2.31 are

$$
\overline g_{\rm H}
=
\sum_{j=0}^{7}{}'C_j(\gamma^2)T_j(x_u),
\qquad
C_j=
\sum_{i=0}^{10}{}'D_{ij}T_i(x_g),
$$

$$
x_g=\frac{\log_{10}\gamma^2}{3}.
$$

Each prime halves the corresponding zeroth term.

The local text at [refs/hummer_1988_gaunt.txt:223](/home/tobias/Projects/neutronstars/refs/hummer_1988_gaunt.txt:223) contains `xu = (2 log « - 2.5)/5.5`. Inspection of PDF page 480 confirms that the minus sign is printed in the paper; it is not merely an OCR error.

**DERIVED HERE.** Mapping the stated interval $\log_{10}u\in[-4,1.5]$ onto $[-1,1]$ requires

$$
\boxed{x_u=\frac{2\log_{10}u+2.5}{5.5}.}
$$

The printed minus sign instead maps the endpoints to $-1.90909$ and $0.0909091$.

I visually transcribed the needed even-$i$ rows from Table 4 on PDF page 482. Its text anchor is “The coefficients Dij are given in Table 4” at [refs/hummer_1988_gaunt.txt:229](/home/tobias/Projects/neutronstars/refs/hummer_1988_gaunt.txt:229). At $\gamma^2=1$, the resulting ordinary, unprimed Chebyshev series has coefficients

$$
\begin{aligned}
\widetilde C={}&(
2.483311193546501,\,
-2.104762769191400,\,
0.449116066681900,\\
&-0.013829706548000,\,
-0.062878951909200,\,
-0.013918632415400,\\
&0.006498141356780,\,
0.004923773755750).
\end{aligned}
$$

Python evaluation gives:

| $u$ | Corrected-coordinate Hummer fit | Relative error against Coulomb integral | Using the printed coordinate |
|---:|---:|---:|---:|
| 0.01 | 2.655426719712282 | $-0.0332707\%$ | 5.570295652883200 |
| 1 | 1.293480754300608 | $-0.0200023\%$ | 3.197614155877816 |
| 10 | 0.9717526865876134 | $-0.0868562\%$ | 2.173970186644949 |

Hummer quotes maximum fitting error $0.007$ and primary-data accuracy about $0.003$: [refs/hummer_1988_gaunt.txt:213](/home/tobias/Projects/neutronstars/refs/hummer_1988_gaunt.txt:213), lines 213–230. A proposed $1.1\%$ compatibility tolerance accommodates both stated contributions. The fit must remain a verification fixture, with the coordinate correction explicitly documented.

## 6. Magnetic extension: what the zero-field result does and does not establish

PC03’s prescription is

$$
R_0(u,\gamma^2)
=
\frac{\overline g_{\rm ff}(u,\gamma^2)}
{\overline g_{\rm Born}(u)},
$$

$$
\boxed{
\Lambda_\alpha^{\rm approx}(B,T,\omega)
=
R_0(u,\gamma^2)
\Lambda_\alpha^{\rm Born}(B,T,\omega).
}
$$

Its wording is deliberately limited: “In order to partly correct” and “recover the non-Born Gaunt factor” at zero field. See [refs/potekhin_chabrier_2003_ff_opacity.tex:1853](/home/tobias/Projects/neutronstars/refs/potekhin_chabrier_2003_ff_opacity.tex:1853), lines 1853–1862.

**DERIVED HERE.** If the magnetic Born calculation has the correct limit,

$$
\lim_{B\to0}\Lambda_\alpha^{\rm Born}
=e^{u/2}K_0(u/2),
$$

then

$$
\lim_{B\to0}\Lambda_\alpha^{\rm approx}
=
\frac{\pi}{\sqrt3}\overline g_{\rm ff}.
$$

That proves zero-field consistency. It does not prove accuracy at finite field.

The rescaling assumes that the relative non-Born correction is independent of magnetic field and polarization. Nothing in the quoted prescription establishes that assumption or bounds its error.

It must also be applied to the **collision frequency and its occurrence in damping**, not merely multiplied onto an already assembled opacity. PC03’s finite-mass expression contains

$$
\sigma_\alpha^{\rm ff}\approx
\frac{4\pi e^2}{m_ec}
\frac{\omega^2\nu_\alpha^{\rm ff}}
{(\omega+\alpha\omega_{ce})^2
(\omega-\alpha\omega_{cp})^2
+\omega^2\widetilde\nu_\alpha^2}.
$$

The source is [refs/potekhin_chabrier_2003_ff_opacity.tex:1714](/home/tobias/Projects/neutronstars/refs/potekhin_chabrier_2003_ff_opacity.tex:1714), lines 1714–1740, with the equation anchor `\sigma_\alpha^\mathrm{ff} \approx`.

Near a resonance, changing $\nu^{\rm ff}$ changes both the numerator and the width. Opacity is then not simply proportional to $R_0$.

### A separate correction must not disappear into `gaunt_ratio`

PC03’s finite-proton-mass quantum correction, described as a “correction factor” “of the order of unity” at lines 1734–1738, is distinct from the zero-field non-Born rescaling.

Also, Appendix A is primarily the **proton–proton** calculation. The electron–proton finite-mass calculation is Appendix B, explicitly titled “Photoabsorption due to Electron-Proton Collisions” at [refs/potekhin_chabrier_2003_ff_opacity.tex:2395](/home/tobias/Projects/neutronstars/refs/potekhin_chabrier_2003_ff_opacity.tex:2395). The general transition formula in Appendix A is useful, but it does not make the Appendix B result non-Born.

### What a full magnetic calculation requires

PC03 itself states that an accurate treatment requires “a solution of the coupled-channel equations” and an $R$-matrix calculation, then says: “However, we will restrict to the first Born approximation.” See [refs/potekhin_chabrier_2003_ff_opacity.tex:1633](/home/tobias/Projects/neutronstars/refs/potekhin_chabrier_2003_ff_opacity.tex:1633).

The required calculation is therefore:

1. Solve the coupled longitudinal equations for Coulomb-distorted electron–proton continuum states at fixed field and transverse pseudomomentum.
2. Normalize incoming and outgoing channels consistently.
3. Evaluate the two-particle current matrix elements between those states.
4. Sum open final channels and thermally populated initial Landau states.
5. Integrate longitudinal momentum and transverse pseudomomentum.
6. Treat thresholds, resonances, and damping consistently.

The coupled equations and Coulomb coupling potential are local at [refs/potekhin_chabrier_2003_ff_opacity.tex:411](/home/tobias/Projects/neutronstars/refs/potekhin_chabrier_2003_ff_opacity.tex:411), anchored by “the set of the \emph{coupled channel equations}”. The transition formula

$$
d\sigma=
\frac{4\pi^2}{\omega c}
|\boldsymbol e\cdot\langle f|\boldsymbol j_{\rm eff}|i\rangle|^2
\delta(E_f-E_i-\hbar\omega)\,d\nu_f
$$

is printed at [refs/potekhin_chabrier_2003_ff_opacity.tex:2287](/home/tobias/Projects/neutronstars/refs/potekhin_chabrier_2003_ff_opacity.tex:2287), where $d\nu_f$ is explicitly the “density of final states”. The thermal state sums are at lines 2631–2667.

This is a separate research path. Owning $R_0$ retires the zero-field table dependency; it does not retire the uncertainty in magnetic non-Born physics.

## 7. Concrete implementation and verification plan

The minimal honest implementation is a reusable Coulomb kernel, a controlled thermal integral, and a separately labelled magnetic prescription. Reproducing the entire published parameter grid is unnecessary before solving the atmosphere’s actual domain.

| Function | Inputs | Outputs and required checks |
|---|---|---|
| `coulomb_parameters` | $\epsilon_i,w>0$ | $\eta_i,\eta_f,\delta,x_H,\beta$; stable differences; explicit invalid-input errors |
| `hyp2f1_transformed` | Complex parameters, transformed $z$, tolerance | Value, absolute error estimate, term count, cancellation estimate |
| `gaunt_ff_exact` | $\epsilon_i,w$, tolerance | Real positive $g$, error estimate, precision used, method; controlled imaginary residual |
| `gaunt_ff_asymptotic` | $\epsilon_i,w$, tolerance | Series value, remainder status, validity diagnostics; explicit refusal outside supported accuracy |
| `gaunt_ff_thermal` | $u,\gamma^2>0$, tolerance | $\overline g$, separate kernel/quadrature/tail errors |
| `gaunt_ff_born_thermal` | $u>0$ | Scaled-Bessel or asymptotic evaluation, avoiding exponential overflow |
| `gaunt_ratio` | $u,\gamma^2$ | $R_0$ with propagated error |
| `ff_collision_frequency` | Local plasma state, polarization, selected magnetic approximation | Collision frequency and provenance/approximation metadata |
| `ff_cross_section` | Frequencies and collision rates | Cross section with the same corrected rates in numerator and damping |
| `gaunt_fixture_reader` | Checksummed verification file | Explicitly named columns; never imported by the runtime atmosphere path |

Follow red–green development:

1. **Parameter and parser tests first.** Use unequal $u,\gamma^2$ to detect column transposition. Reject nonpositive, nonfinite, or out-of-contract inputs with their values in the error.
2. **Born reduction.** Recover the logarithmic expression and the numerical target $0.6056966996081959$ in its limiting test.
3. **Cancellation regression.** Recover $1.304895799335132149$ at $(\epsilon_i,w)=(10^{-10},10^{-15})$ to the allocated kernel tolerance. Verify that insufficient precision is detected.
4. **Classical and soft limits.** Check the derived powers, the unit Kramers limit, and the soft digamma expression. Test overlap between exact and asymptotic branches.
5. **Thermal fixtures.** Meet the three rounding-aware van Hoof comparisons above; separately match the computed reference values to relative $10^{-8}$.
6. **Legacy comparisons.** Print the Sutherland and Hummer discrepancies; do not tune the new answer toward the older tables.
7. **Magnetic limit.** Require the corrected magnetic collision rate to approach the exact zero-field rate as $B\rightarrow0$. Test weak damping separately from resonance behavior.
8. **Atmosphere response.** Perturb free-free opacity while reconverging the thermal structure; measure changes in spectral intensity, angular dependence, and synthetic HST measurements.

Mutation testing should change the hypergeometric phases, omit the final $I_0$, reverse the energy ordering, exchange $u$ and $\gamma^2$, insert an erroneous $\sqrt t$ weight, and remove a stimulated-emission factor. Each mutation should be killed by a specific numerical test.

Fuzz logarithmic input coordinates and concentrate samples near $x_H=-1$, $\eta_i\simeq\eta_f$, the thermal endpoint, and the asymptotic transition. Agreement between two precision levels is necessary but does not replace branch-overlap and analytic-limit tests.

The current `pipeline.gaunt-factor` scope is appropriate: magnetic correction is downstream. Its missing exact-expression and averaging dependencies should be added. `opacity.approx-ff-non-born-correction` needs an explicit dependency on the owned zero-field ratio, while retaining its status as a magnetic approximation.

All calculations in this argument ran in memory. No production implementation was added. The read-only filesystem prevented the required `HANDOFF.md` and worklog updates.

## 8. Missing sources and data to fetch

The principal zero-field calculation is already locally sourced. These acquisitions improve independent verification or address adjacent physical limits; they should not delay the first owned kernel.

| Acquisition | Purpose | Bibliographic status |
|---|---|---|
| Biedenharn, 1956, *Physical Review* **102**, 262 | Original exact reduction and normalization audit | **UNVERIFIED acquisition**; citation printed in [refs/vanhoof_2014_gaunt.txt:726](/home/tobias/Projects/neutronstars/refs/vanhoof_2014_gaunt.txt:726): `Biedenharn L. C., 1956, Physical Review, 102, 262`. No arXiv identifier verified. |
| Menzel & Pekeris, 1935, *MNRAS* **96**, 77 | Complete contour derivation underlying Appendix A; independent remainder work | **UNVERIFIED acquisition**; citation at [refs/vanhoof_2014_gaunt.txt:734](/home/tobias/Projects/neutronstars/refs/vanhoof_2014_gaunt.txt:734): `Menzel D. H., Pekeris C. L., 1935, MNRAS, 96, 77`. No arXiv identifier verified. |
| Bethe & Salpeter, 1957, *Quantum Mechanics of One- and Two-Electron Atoms*, Springer | Independently normalized continuum-wave and matrix-element calculation | **UNVERIFIED acquisition**; title and publisher printed at [refs/potekhin_chabrier_2003_ff_opacity.tex:2723](/home/tobias/Projects/neutronstars/refs/potekhin_chabrier_2003_ff_opacity.tex:2723). No arXiv identifier verified. |
| van Hoof et al., 2014, *MNRAS* **444**, 420, electronic tables with error columns | More precise and broader verification fixtures | **UNVERIFIED data acquisition**. Correct paper identifier: **arXiv:1407.5048**, printed in local TeX line 1 and corroborated by the [arXiv record](https://arxiv.org/abs/1407.5048). The audit’s other identifier is wrong. |
| Nozawa, Itoh & Kohyama, 1998, *ApJ* **507**, 530 | Quantify departure from the nonrelativistic calculation at the hard-energy end | **UNVERIFIED acquisition**; citation printed at [refs/vanhoof_2014_gaunt.txt:735](/home/tobias/Projects/neutronstars/refs/vanhoof_2014_gaunt.txt:735). No arXiv identifier verified. |

The quoted author strings are bibliographic leads, not substitutes for reading the acquired papers. Hummer, KL61, van Hoof’s paper, PC03, and the local Potekhin 1994 source should not remain on a “missing paper” list.

## 9. Stub policy, spectrum impact, and R5

### 9.1 Quantified opacity effects at representative RX J1856 parameters

Use Ho’s published model parameters only as a labelled reference case, not as newly inferred measurements. The local source gives $T_\infty\simeq4.3\times10^5$ K, $z\simeq0.22$, and $B\simeq(3\!-\!4)\times10^{12}$ G: [refs/ho_2007_rxj1856_thin_atm.tex:70](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:70), lines 70–78. Its description is explicitly “thin, magnetic, partially ionized hydrogen” at line 67.

**DERIVED HERE.** Applying the redshift to both photon energy and thermal temperature gives

$$
T_{\rm loc}=(1+z)T_\infty=524600\ {\rm K},
$$

$$
k_BT_{\rm loc}=45.2065302932\ {\rm eV},
\qquad
\gamma^2=0.300967427377.
$$

The Python calculation used CODATA inputs from [refs/tiesinga_2021_codata2018.txt:3403](/home/tobias/Projects/neutronstars/refs/tiesinga_2021_codata2018.txt:3403), lines 3403–3426: the greppable values include `299 792 458`, `6.626 070 15`, `9.109 383 7015`, `7.297 352 5693`, `1.380 649`, and `1.602 176 634`. The Rydberg was computed from $\alpha^2m_ec^2/2$, giving $13.60569312296775$ eV.

The following photon energies are at infinity. Surface energies are larger by $1+z$, while $u=E_\infty/(k_BT_\infty)$.

| $E_\infty$ | $u$ | Exact zero-field $\overline g$ | Born $\overline g$ | Exact/Born |
|---:|---:|---:|---:|---:|
| 2 eV | 0.0539745 | 2.11051012 | 2.11224789 | 0.999177 |
| 0.1 keV | 2.698725 | 1.05153092 | 0.55370992 | 1.899065 |
| 0.3 keV | 8.096175 | 0.81156936 | 0.33402347 | 2.429678 |
| 1 keV | 26.987251 | 0.55582748 | 0.18643264 | 2.981385 |
| 10 keV | 269.872514 | 0.21694866 | 0.05942997 | 3.650492 |

At this temperature, Born-only free-free opacity is lower by approximately $47\%$, $59\%$, $66\%$, and $73\%$ at the four X-ray energies. These are **free-free opacity differences at fixed state**, not emergent-flux differences.

The optical agreement at 2 eV is accidental. Across $E_\infty=1.5$–$3.5$ eV, the computed exact/Born ratio runs from $0.984109$ to $1.039617$. At 8.3 eV it is $1.142301$. In a deliberately cooler test layer at $T_{\rm loc}=200000$ K, the 8.3 eV ratio becomes $1.342751$. The atmosphere must evaluate the factor at the local temperature, not once at $T_{\rm eff}$.

Setting $g=1$ is much worse in the optical: at 2 eV it lowers the free-free opacity by $52.6\%$.

At the reference temperature, the high-$u$ ratio derived earlier tends to $4.04476538$, rather than unity. This rules out a hard-energy Born fallback even when the outgoing electron energy greatly exceeds the Rydberg.

The mathematical accuracy should also not be confused with physical accuracy at the upper band edge: a 10 keV photon at infinity has surface energy $12.2$ keV here, or $0.0238748\,m_ec^2$. That dimensionless ratio is not an error bound, but it makes a relativistic comparison necessary before claiming uniformly tiny physical errors throughout the band.

### 9.2 Magnetic uncertainty is a different budget

For $B=4\times10^{12}$ G, PC03’s printed approximation

$$
\beta_e\simeq134.3\,B_{12}/T_6
$$

gives $\beta_e\simeq1024$ at the reference temperature. Source: [refs/potekhin_chabrier_2003_ff_opacity.tex:144](/home/tobias/Projects/neutronstars/refs/potekhin_chabrier_2003_ff_opacity.tex:144), anchored by `\approx 134.3\,B_{12}/T_6`.

The zero-field factor is therefore not directly applicable as the magnetic coefficient. It supplies only PC03’s approximate rescaling.

Using the same rounded cyclotron relation and the CODATA mass ratio gives a proton cyclotron energy about $25.21$ eV at the surface. Optical photons lie below it. Finite-ion-motion effects therefore enter the optical calculation through a separate part of the cross section, irrespective of how accurately $R_0$ is evaluated.

Nor does solving free-free scattering remove bound opacity. Ho’s local source explicitly reports bound-bound and bound-free features even with atomic fraction below $10\%$: [refs/ho_2007_rxj1856_thin_atm.tex:274](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:274), lines 274–292, including “the other features are due to bound-bound and bound-free transitions.”

There is no locally established percentage uncertainty for PC03’s non-Born rescaling at this field. Assigning it a precise error bar would be invented.

### 9.3 Translating opacity error into a render estimate

**DERIVED HERE, as a sensitivity model rather than an atmosphere prediction.**

Suppose a fixed background has $T(y)\propto y^a$, and multiplying the relevant opacity by $f$ moves its formation column approximately as $y_\nu\propto f^{-1}$. Then

$$
\delta\ln T_\nu\simeq-a\,\delta\ln f.
$$

For a thermal source, summing oscillator occupations gives
$\langle n\rangle=\sum n e^{-nu}/\sum e^{-nu}=1/(e^u-1)$. Differentiation therefore gives

$$
\frac{\partial\ln B_\nu}{\partial\ln T}
=
\frac{u}{1-e^{-u}},
$$

so that

$$
\delta\ln I_\nu
\simeq
-a\,\frac{u}{1-e^{-u}}\,\delta\ln f.
$$

In the optical Rayleigh–Jeans regime this is approximately $-a\,\delta\ln f$. Taking $a=1/4$ as an explicit illustrative choice:

- Replacing the exact optical factor $2.11051$ by unity gives a flux multiplier $2.11051^{1/4}=1.20530$, or magnitude shift $-0.20274$.
- A $0.13\%$ opacity perturbation gives about $0.000353$ mag.
- A $1\%$ opacity perturbation gives about $0.00271$ mag.
- The computed Born-ratio variation across 1.5–3.5 eV corresponds to about $0.0149$ mag of differential response in this toy model.

The magnitude conversion follows the locally printed definition
`${\rm STMAG}=-2.5\log_{10} F_{\lambda}-21.1$` at [refs/kaplan_2011_rxj1856_spectrum.tex:488](/home/tobias/Projects/neutronstars/refs/kaplan_2011_rxj1856_spectrum.tex:488).

These estimates are not bounds. Scattering, depth-dependent opacity, temperature readjustment, and the thin atmosphere’s lower boundary can change the response. Where other opacity dominates, the perturbation to total opacity is also reduced by the free-free fraction.

### 9.4 What can be stubbed honestly

| Treatment | Defensible use | Consequence |
|---|---|---|
| Published thermal table or corrected Hummer fit | Verification; explicitly labelled development comparisons | Small numerical discrepancy within a validated domain, but no ownership claim. Not an acceptable final runtime dependency under this task. |
| Born factor | Analytic tests and labelled comparison spectra | Useful diagnostic; inadequate as an unqualified X-ray opacity substitute at the reference state. |
| $g=1$ | Kramers-limit tests | Not a credible general optical atmosphere approximation here. |
| Exact zero-field factor with PC03 magnetic rescaling | A labelled magnetic baseline with sensitivity comparisons | Correct zero-field limit; unquantified finite-field non-Born uncertainty remains. |
| Silent table clamp or extrapolation | None | Destroys both the numerical error statement and traceability. |

For R5, a validated thermal-table approximation would probably be a smaller error source than the magnetic and surface-model assumptions, but that judgment must be checked by reconverged atmosphere runs. It does not waive the ownership requirement.

The local photometric scale is concrete: Walter reports `25.36$\pm$0.15` in F475W at [refs/walter_2010_rxj1856_distance.tex:862](/home/tobias/Projects/neutronstars/refs/walter_2010_rxj1856_distance.tex:862). Kaplan’s fitted optical/UV summary gives `$4.25\pm0.07$` at 2500 Å, a computed fractional uncertainty of $1.65\%$, at [refs/kaplan_2011_rxj1856_spectrum.tex:1035](/home/tobias/Projects/neutronstars/refs/kaplan_2011_rxj1856_spectrum.tex:1035). The latter is a fitted summary, not a replacement for the individual HST measurements.

The illustrative $0.20$ mag error from a unit Gaunt factor is large enough to threaten R5 directly. A validated subpercent thermal fit would plausibly produce much smaller direct optical changes. Born-only X-ray opacity can still bias the inferred temperature, emitting area, and hence optical normalization even where its direct optical correction happens to be small.

**The minimal honest completion is therefore an owned zero-field Coulomb calculation, an explicitly approximate magnetic extension, and an atmosphere sensitivity calculation carried through the actual HST bandpasses. Until that final calculation is done, the remaining effect on R5 is a measured sensitivity to be obtained—not an error bar that this analytic argument can certify.**