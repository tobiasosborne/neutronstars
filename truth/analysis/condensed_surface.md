## Abstract

The reflectivity fits can be replaced by a direct Maxwell boundary calculation, but this does **not** complete the condensed-surface physics. It produces a conditional answer given a dielectric tensor, density, ion-response model, and material temperature. The unresolved ion response has a large optical consequence: in one explicitly calculated example below, free-ion and fixed-ion emissivities differ by a factor of $4.65$ at $2$ eV. Condensation itself is plausible for iron under some RX J1856 parameter choices, but is not established by its fitted effective temperature. A thin atmosphere can hide the substrate’s direct optical emission; whether it also suppresses the uncertainty in the emergent spectrum requires a coupled transfer calculation.

All numerical results below were calculated with Python. These are exploratory calculations, not a committed production implementation. The workspace is read-only, so I could not update `HANDOFF.md` or the worklog.

## 1. What this critical path can establish

There are two bookkeeping corrections. In the current `truth/DEPENDENCY_AUDIT.md`, A24–A27 concern the BSk EOS, TOV, opacity-table reader, and nuclear-density placeholder. Condensed emission appears under Part C’s new Milestone-3 dependencies. The current spine node is `pipeline.surface-radiance-condensed`, and the relevant entries reside mainly in `truth/entries/condensed.toml`.

The substantive correction is more important: **solving Fresnel reflection from an assumed dielectric tensor owns the electromagnetic boundary calculation, not the underlying material response.**

Van Adelsberg et al. explicitly introduce their closure as “consider the free electron gas model for the condensed matter” and later call it “only a first approximation to the true” dielectric tensor. [refs/van_adelsberg_2005_condensed_surface.tex:395](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:395), [refs/van_adelsberg_2005_condensed_surface.tex:513](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:513).

The minimal honest calculation therefore has distinct layers:

1. **Electromagnetics:** compute the complex reflection matrix from a supplied passive dielectric tensor.
2. **Material closure:** supply density, damping, and free-ion, fixed-ion, or ultimately lattice response.
3. **Phase consistency:** determine whether condensed matter and the assumed overlying column can exist at the computed interface temperature.
4. **Atmosphere and observation:** propagate the boundary radiation through the atmosphere, then fit the observed spectrum and photometry.

Replacing the first layer’s fit does not make the other layers green.

## 2. Derivation of the interface problem

### 2.1 Assumptions and conventions

The calculation below assumes:

- A plane, sharp interface; a homogeneous condensed half-space at $z<0$; vacuum at $z>0$.
- A constant magnetic field and composition within the patch.
- A local dielectric response, with spatial dispersion neglected.
- Unit magnetic permeability for the wave calculation.
- No separately imposed sheet current or free surface-charge layer.
- A smooth surface, hence specular reflection.
- Local thermodynamic equilibrium in the emitting material.
- A material temperature $T_s$ approximately constant over the emitting penetration depth.
- No radiation incident from an independently illuminated rear boundary.

These are model assumptions, not inferred properties of RX J1856. The plane-patch approximation is explicit in “We treat this small patch as plane, neglecting” its curvature and roughness. The bare-surface calculation assumes temperature “at least a few times lower than $T_{\rm crit}$”. [refs/potekhin_2012_condensed_surface.tex:503](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:503), [refs/van_adelsberg_2005_condensed_surface.tex:607](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:607).

Adopt fields proportional to $\exp[i(\boldsymbol{k}\cdot\boldsymbol{x}-\omega t)]$, matching the local source’s `\mathrm{e}^{\mathrm{i}(\bm{k}\cdot\bm{r}-\omega t)}`. [refs/potekhin_2012_condensed_surface.tex:433](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:433).

Choose the incidence plane to be $xz$, and define

$$
s=\sin\theta,\qquad c_\theta=\cos\theta,
$$

$$
\boldsymbol{n}_i=(s,0,-c_\theta),\qquad
\boldsymbol{n}_r=(s,0,c_\theta),\qquad
\boldsymbol{n}=\frac{c_{\rm light}\boldsymbol{k}}{\omega}.
$$

The field direction is an arbitrary unit vector $\boldsymbol b$. Its azimuth in this coordinate system replaces the paper’s choice of an arbitrary incidence azimuth relative to a fixed field plane.

Use explicit linear polarizations:

$$
\boldsymbol e_s^{\,i}=\boldsymbol e_s^{\,r}=(0,1,0),
$$

$$
\boldsymbol e_p^{\,i}=(c_\theta,0,s),\qquad
\boldsymbol e_p^{\,r}=(c_\theta,0,-s).
$$

These agree with the vector definitions in Potekhin’s Appendix B after setting its incidence azimuth to zero; the anchor is `\bm{e}_1^\mathrm{(i)} = \bm{e}_1^\mathrm{(r)}`. [refs/potekhin_2012_condensed_surface.tex:2096](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:2096).

**Do not implement the basis from the prose labels alone.** Lines 531–532 describe the first polarization as parallel to the incidence plane, whereas Appendix B explicitly makes it perpendicular. Also, Potekhin reverses the amplitude-index order relative to van Adelsberg: “order of the subscripts in $r_{mj}$ and $t_{mj}$”. [refs/potekhin_2012_condensed_surface.tex:531](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:531), [refs/potekhin_2012_condensed_surface.tex:1787](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:1787).

Hereafter the **first matrix index is outgoing polarization**.

### 2.2 Cold-plasma dielectric tensor, with ions and damping

In coordinates aligned with the field, write

$$
\boldsymbol\epsilon_B=
\begin{pmatrix}
\epsilon_\perp&ig&0\\
-ig&\epsilon_\perp&0\\
0&0&\epsilon_\parallel
\end{pmatrix},
$$

$$
\epsilon_\perp\pm g
\simeq
1-\frac{\omega_{pe}^2}
{(\omega\pm\Omega_e)(\omega\mp\Omega_i)
+i\omega\nu_\perp},
\qquad
\epsilon_\parallel
\simeq
1-\frac{\omega_{pe}^2}{\omega(\omega+i\nu_\parallel)}.
$$

These are the locally printed equations, with renamed tensor components. Verbatim anchors are `\epsilon & i\,g & 0` and `i\,\gamma_{ei}^{(tr)}`; the source defines $\gamma=\nu/\omega$. [refs/van_adelsberg_2005_condensed_surface.tex:404](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:404), [refs/van_adelsberg_2005_condensed_surface.tex:411](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:411), [refs/van_adelsberg_2005_condensed_surface.tex:464](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:464).

Here

$$
\Omega_e=\frac{eB}{m_ec_{\rm light}},\qquad
\Omega_i=\frac{ZeB}{m_ic_{\rm light}},\qquad
\omega_{pe}^2=\frac{4\pi n_e e^2}{m_e},
\qquad m_i=A m_u.
$$

The same definitions occur with the anchor “the electron and ion cyclotron frequencies”. [refs/potekhin_lai_chabrier_ho_2004_polarization.tex:519](/home/tobias/Projects/neutronstars/refs/potekhin_lai_chabrier_ho_2004_polarization.tex:519).

**DERIVED HERE: origin of the ion denominator.** Consider the circular component whose dielectric eigenvalue is $\epsilon_\perp+g$. A momentum-conserving relaxation model gives

$$
[\nu_\perp-i(\omega+\Omega_e)]v_e-\nu_\perp v_i
=-\frac{eE}{m_e},
$$

$$
-\delta\nu_\perp v_e+
[\delta\nu_\perp-i(\omega-\Omega_i)]v_i
=\frac{ZeE}{m_i},
\qquad
\delta=\frac{Zm_e}{m_i}.
$$

Since $\Omega_i=\delta\Omega_e$, the determinant reduces to

$$
-\left[(\omega+\Omega_e)(\omega-\Omega_i)
+i\omega\nu_\perp(1+\delta)\right].
$$

Using $j=-n_e e(v_e-v_i)$ gives

$$
\epsilon_+
=
1-\frac{\omega_{pe}^2(1+\delta)}
{(\omega+\Omega_e)(\omega-\Omega_i)
+i\omega\nu_\perp(1+\delta)}.
$$

Neglecting the small mass-ratio correction yields the printed approximation. The source specifically corrects an earlier damping expression and identifies the omitted factor as `1 + Z m_e/A m_p`. [refs/van_adelsberg_2005_condensed_surface.tex:399](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:399).

The two ion prescriptions are therefore:

- **Free ions:** retain $\Omega_i$ in the denominator.
- **Fixed ions:** remove ion dynamics while holding the physical electron density fixed, giving

$$
\epsilon_\perp\pm g
=
1-\frac{\omega_{pe}^2}
{\omega(\omega\pm\Omega_e+i\nu_\perp)}.
$$

Do not implement “$m_i\rightarrow\infty$” by also changing the composition used to compute $n_e$. It is a dynamical limit of a specified material. The source describes it as “the motion of the ions is neglected”. [refs/van_adelsberg_2005_condensed_surface.tex:521](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:521).

**DERIVED HERE: rotation without angle-specific formulas.** Define $[\boldsymbol b]_\times\boldsymbol v=\boldsymbol b\times\boldsymbol v$. Then

$$
\boxed{
\boldsymbol\epsilon
=
\epsilon_\perp\mathbf I
+(\epsilon_\parallel-\epsilon_\perp)\boldsymbol b\boldsymbol b^T
-ig[\boldsymbol b]_\times .
}
$$

For $\boldsymbol b=\hat{\boldsymbol z}$ this reproduces the printed tensor. It also exposes the sign convention and supports arbitrary orientation without copying a specialized rotated matrix.

### 2.3 Deriving the two transmitted modes

Tangential phase continuity fixes $n_x=s$ and $n_y=0$. Only the normal component $q=n_z$ is unknown:

$$
\boldsymbol n_t=(s,0,q).
$$

The local source explicitly rewrites the transmitted wavevector using its normal component, `+ n_{z,j}\,{\hat\mathbf{z}}`. [refs/potekhin_2012_condensed_surface.tex:1843](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:1843).

The wave equation is

$$
\left[
\boldsymbol\epsilon+\boldsymbol n_t\boldsymbol n_t^T
-(s^2+q^2)\mathbf I
\right]\boldsymbol E=0,
$$

with source anchor `{\bm\epsilon} + \bm{n}_j\otimes\bm{n}_j`. [refs/potekhin_2012_condensed_surface.tex:1859](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:1859).

**DERIVED HERE:** expanding this equation gives

$$
\begin{pmatrix}
\epsilon_{xx}-q^2&\epsilon_{xy}&\epsilon_{xz}+sq\\
\epsilon_{yx}&\epsilon_{yy}-s^2-q^2&\epsilon_{yz}\\
\epsilon_{zx}+sq&\epsilon_{zy}&\epsilon_{zz}-s^2
\end{pmatrix}
\boldsymbol E=0.
$$

Its determinant is quartic in $q$. This formulation does not introduce the extra branches generated by squaring expressions involving a complex refraction angle.

For implementation, a direct four-component Maxwell eigenproblem is preferable to naively expanding that determinant. From

$$
\boldsymbol H=\boldsymbol n_t\times\boldsymbol E,
\qquad
\boldsymbol n_t\times\boldsymbol H=-\boldsymbol\epsilon\boldsymbol E,
$$

one obtains

$$
H_z=sE_y,\qquad
E_z=-\frac{\epsilon_{zx}E_x+\epsilon_{zy}E_y+sH_y}{\epsilon_{zz}}.
$$

Define $d=\epsilon_{zz}$ and the transverse Schur complement

$$
\bar\epsilon_{ab}
=\epsilon_{ab}-\frac{\epsilon_{az}\epsilon_{zb}}{d},
\qquad a,b\in\{x,y\}.
$$

Then, for $\boldsymbol f=(E_x,E_y,H_x,H_y)^T$,

$$
\boxed{
q\boldsymbol f=
\begin{pmatrix}
-s\epsilon_{zx}/d&-s\epsilon_{zy}/d&0&1-s^2/d\\
0&0&-1&0\\
-\bar\epsilon_{yx}&s^2-\bar\epsilon_{yy}&0&s\epsilon_{yz}/d\\
\bar\epsilon_{xx}&\bar\epsilon_{xy}&0&-s\epsilon_{xz}/d
\end{pmatrix}
\boldsymbol f .
}
$$

This is a derivation from the same Maxwell equations, not a new material approximation. If $d$ is singular or the elimination is ill-conditioned, use another elimination or the unreduced generalized system; do not replace $d$ by a floor.

For the adopted time convention and material occupying $z<0$,

$$
|\boldsymbol E|\propto
\exp[-(\omega/c_{\rm light})\operatorname{Im}(q)z].
$$

Hence the decaying modes satisfy

$$
\boxed{\operatorname{Im}q<0.}
$$

The source motivates its selection by decay in the “transmitted wave region, leading to the condition:”. Its displayed notation uses $n_j$, but the relevant quantity in this geometry is the **normal component** $n_{z,j}$. [refs/potekhin_2012_condensed_surface.tex:1945](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:1945).

In the lossless limit, select propagating modes with inward normal energy flux,

$$
S_z\propto\operatorname{Re}(\boldsymbol E\times\boldsymbol H^*)_z<0,
$$

and decaying evanescent modes. Equivalently, continue from positive dissipation toward zero. Selecting by $\operatorname{Re}q$ alone is not the physical criterion.

At repeated eigenvalues, retain the complete two-dimensional physical subspace. Computing one nullvector twice produces two identical transmitted columns and a singular boundary system.

### 2.4 Fresnel matrix for arbitrary incidence and field direction

The locally printed boundary conditions include

$$
\Delta\boldsymbol E\times\hat{\boldsymbol z}=0,
\qquad
\Delta\boldsymbol H\times\hat{\boldsymbol z}=0,
$$

with anchor `\Delta\bE\times\hatz=0`. Normal $\boldsymbol D$ and $\boldsymbol B$ continuity are not additional independent equations for these plane-wave solutions. [refs/van_adelsberg_2005_condensed_surface.tex:784](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:784), [refs/van_adelsberg_2005_condensed_surface.tex:796](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:796).

**DERIVED HERE:** in the chosen $(s,p)$ basis, the incident and reflected tangential-field matrices are

$$
W_i=
\begin{pmatrix}
0&c_\theta\\
1&0\\
c_\theta&0\\
0&-1
\end{pmatrix},
\qquad
W_r=
\begin{pmatrix}
0&c_\theta\\
1&0\\
-c_\theta&0\\
0&1
\end{pmatrix}.
$$

Let $W_t$ contain the tangential fields of the two physical transmitted modes. For incident amplitude vector $\boldsymbol a$,

$$
W_i\boldsymbol a+W_r r\boldsymbol a=W_t t\boldsymbol a.
$$

Thus

$$
\boxed{
\begin{pmatrix}r\\t\end{pmatrix}
=
[\,W_r\; -W_t\,]^{-1}(-W_i),
\qquad
r=
\begin{pmatrix}
r_{ss}&r_{sp}\\
r_{ps}&r_{pp}
\end{pmatrix}.
}
$$

The inverse denotes a linear solve, not an instruction to form an inverse numerically. The free-ion and fixed-ion matrices follow by inserting their respective dielectric tensors into this same calculation.

Because both exterior channels are in vacuum at the same incidence angle, their flux-normalization factors agree. Consequently,

$$
\mathcal R_{ab}=|r_{ab}|^2
$$

is the reflected power fraction from incident channel $b$ into outgoing channel $a$.

A useful analytic reduction occurs for field and incidence normal to the interface. The transmitted modes are circular:

$$
n_\pm^2=\epsilon_\perp\pm g,
\qquad
r_\pm=\frac{1-n_\pm}{1+n_\pm},
$$

up to the chosen reflected-basis phase convention. The mean reflectivity is

$$
R=\frac12\left|
\frac{n_+-1}{n_++1}
\right|^2+
\frac12\left|
\frac{n_--1}{n_-+1}
\right|^2.
$$

This is the printed normal-incidence result; the source gives `n_{1}^2=\epsilon+g`, `n_{2}^2=\epsilon-g`. [refs/van_adelsberg_2005_condensed_surface.tex:862](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:862), [refs/van_adelsberg_2005_condensed_surface.tex:865](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:865).

### 2.5 Kirchhoff’s law requires row sums—and full coherency

The source derives emission by putting the body in a blackbody cavity and requiring the outgoing radiation to remain blackbody. Its scalar result is `I_\nu^{(e)}=(1-R)B_\nu(T)`. [refs/van_adelsberg_2005_condensed_surface.tex:614](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:614), [refs/van_adelsberg_2005_condensed_surface.tex:674](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:674).

**DERIVED HERE: matrix form.** Let $C$ be the electric-field coherency matrix, normalized so that $\operatorname{tr}C=I_E$. An unpolarized equilibrium incident field has

$$
C_{\rm eq}=\frac{B_E(T_s)}2\mathbf I.
$$

Reflection produces $rC_{\rm eq}r^\dagger$. Equilibrium therefore requires

$$
\boxed{
C_{\rm em}=
\frac{B_E(T_s)}2\left(\mathbf I-r r^\dagger\right).
}
$$

Taking its diagonal gives

$$
J_a=1-\sum_b|r_{ab}|^2,\qquad
I_{E,a}=\frac12J_aB_E,
$$

and

$$
J=1-\frac12\operatorname{tr}(rr^\dagger).
$$

This reproduces Potekhin’s explicit **row-sum** definition, `R_j = |r_{j1}|^2+|r_{j2}|^2`. [refs/potekhin_2012_condensed_surface.tex:1802](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:1802).

By contrast, absorption of a specified incident channel involves the corresponding **column** of $r$, or $\mathbf I-r^\dagger r$. Do not interchange these matrices.

For the convention

$$
I=C_{ss}+C_{pp},\quad
Q=C_{ss}-C_{pp},\quad
U=2\operatorname{Re}C_{sp},\quad
V=-2\operatorname{Im}C_{sp},
$$

the full emitted Stokes vector follows immediately. These Stokes signs are an explicitly chosen convention.

This exposes a shortcoming in the spine’s promise to return only two-mode intensity and linear polarization. Equal $J_s$ and $J_p$ do not imply unpolarized emission. In the calculated normal-incidence, fixed-ion example at $B=1.5\times10^{13}$ G, $\eta_\rho=0.55$, and $E=2$ eV,

$$
I/B_E=0.1896811,\qquad Q/I=0,\qquad V/I=-1
$$

in the zero-damping limit and the convention above. One circular channel is completely reflected. A renderer that stores only the two linear diagonal intensities loses this information.

Finally, these lossless calculations mean the limit of an absorbing **semi-infinite** material as damping tends to zero. A finite, genuinely transparent slab would require a rear transmission channel; its emission cannot be inferred from $1-R$ alone.

## 3. Material parameters and their sensitivity

### 3.1 Frequencies and density

Charge neutrality gives $n_e=Zn_i$ and $n_i=\rho/(Am_u)$. Therefore, **DERIVED HERE**,

$$
\omega_{pi}^2
=\frac{4\pi Z^2e^2n_i}{Am_u}
=\omega_{pe}^2\frac{Zm_e}{Am_u}.
$$

For numerical comparisons I used the printed rounded energy coefficients:

$$
E_{ce}=115.8B_{13}\ {\rm keV},
$$

$$
E_{ci}=0.0635(Z/A)B_{13}\ {\rm keV},
$$

$$
E_{pe}=0.0288\sqrt{\rho Z/A}\ {\rm keV},
\qquad
E_{pi}=6.75\times10^{-4}(Z/A)\sqrt{\rho}\ {\rm keV}.
$$

Here $\rho$ is in $\mathrm{g\,cm^{-3}}$. Anchors are `1158\, B_{14}`, `0.0288\,\left ( {Z \over A} \right )`, and `=6.75\times 10^{-4}`. [refs/van_adelsberg_2005_condensed_surface.tex:429](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:429), [refs/van_adelsberg_2005_condensed_surface.tex:432](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:432), [refs/van_adelsberg_2005_condensed_surface.tex:498](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:498).

These rounded coefficients explain the numerical precision of the following tables; they are not a substitute for the project’s eventual CODATA implementation.

Potekhin’s density estimate is

$$
\boxed{
\rho_s=
8.9\times10^3\eta_\rho A Z^{-0.6}B_{13}^{1.2}
\ {\rm g\,cm^{-3}}.
}
$$

The verbatim anchor is `8.9\times10^3\,\eta\,AZ^{-0.6}\,B_{13}^{1.2}`. The paper identifies $\eta_\rho=1$ with the ion-sphere model. [refs/potekhin_2012_condensed_surface.tex:320](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:320), [refs/potekhin_2012_condensed_surface.tex:327](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:327).

For $^{56}$Fe, locally quoted alternatives include:

- Zero-temperature Thomas–Fermi:
  $\eta_\rho\simeq0.2+0.0028/B_{13}^{0.56}$, over its stated field range. Anchor: `$\eta\approx0.2+0.0028/B_{13}^{0.56}$`. [refs/potekhin_2012_condensed_surface.tex:330](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:330).
- A representation of Medin–Lai results:
  $\eta_\rho\simeq0.55$, with the later review explicitly giving $0.55\pm0.11$ for $5\leq B_{12}\leq1000$. Anchor: `$\xi=0.55\pm0.11$`. [refs/potekhin_chabrier_2013_mag_eos.tex:2220](/home/tobias/Projects/neutronstars/refs/potekhin_chabrier_2013_mag_eos.tex:2220).

The latter spread is not an observational confidence interval. The same sources report that cohesive energies vary by a factor of two between approximations, and that a finite-temperature Thomas–Fermi model predicts no magnetic condensation. Anchors: “factor of two” and “does not predict magnetic”. [refs/potekhin_2012_condensed_surface.tex:346](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:346), [refs/potekhin_2012_condensed_surface.tex:336](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:336).

### 3.2 Density changes both the emissivity and the spectral break

**DERIVED HERE:** at fixed composition and field,

$$
E_{pe}\propto\rho_s^{1/2},\qquad
E_{pi}\propto\rho_s^{1/2}.
$$

Define

$$
C=\frac{E_{pe}^2}{E_{ce}}\propto\rho_s.
$$

For normal incidence, negligible damping, and photon energy small compared with $E_{ce}$, the circular refractive indices reduce to

$$
n_a^2\simeq1+\frac{C}{E+E_{ci}},
\qquad
n_b^2\simeq1-\frac{C}{E-E_{ci}}.
$$

This is the printed approximation with anchor `1 \pm \frac{\Epe^2}{\Ece(E\pm\Eci)}`. [refs/potekhin_2012_condensed_surface.tex:637](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:637).

The second mode becomes evanescent between $E_{ci}$ and approximately

$$
E_C=E_{ci}+C.
$$

The source prints `\EC=\Eci+\Epe^2/\Ece`. [refs/potekhin_2012_condensed_surface.tex:590](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:590).

For $E\ll E_{ci}$, both free-ion indices approach

$$
n_0=\sqrt{1+\Xi},\qquad
\Xi=\frac{E_{pe}^2}{E_{ce}E_{ci}},
$$

so **DERIVED HERE**

$$
\boxed{J_0=\frac{4n_0}{(1+n_0)^2}.}
$$

The density sensitivity is

$$
\frac{\partial\ln J_0}{\partial\ln\rho_s}
=
\frac{1-n_0}{1+n_0}
\frac{\Xi}{2(1+\Xi)}.
$$

Python evaluation at $B=10^{13}$ G, $Z=26$, $A=56$ gives:

| Density prescription | $\eta_\rho$ | $\rho_s$ [$\mathrm{g\,cm^{-3}}$] | $E_{pe}$ [keV] | $E_C$ [keV] | Free-ion $J_0$ |
|---|---:|---:|---:|---:|---:|
| Thomas–Fermi representation | 0.2028 | 14,310.7 | 2.34756 | 0.077073 | 0.944434 |
| Lower comparison density | 0.44 | 31,048.9 | 3.45787 | 0.132737 | 0.870864 |
| Medin–Lai representation | 0.55 | 38,811.2 | 3.86601 | 0.158550 | 0.842064 |
| Upper comparison density | 0.66 | 46,573.4 | 4.23501 | 0.184364 | 0.816231 |
| Ion-sphere model | 1.00 | 70,565.8 | 5.21293 | 0.264151 | 0.750841 |

Thus a density variation of $\pm20\%$ around $\eta_\rho=0.55$ changes this low-energy free-ion emissivity by approximately $+3.42\%$ and $-3.07\%$. Changing from the published fit’s $\eta_\rho=1$ to $0.55$ changes it by $12.15\%$, and substantially moves the suppression-band endpoint.

Near a moving edge, the emissivity at a fixed photon energy can change much more sharply than this low-energy estimate. Therefore **a density nuisance parameter cannot be represented by a single multiplicative uncertainty on radiance**.

### 3.3 Ion response dominates the optical ambiguity

For fixed ions, **DERIVED HERE**,

$$
n_a^2\simeq1+\frac CE,\qquad
n_b^2\simeq1-\frac CE.
$$

At $E\ll C$, one mode is evanescent and

$$
J_{\rm fixed}
\simeq\frac{2n_a}{(1+n_a)^2}
\longrightarrow2\sqrt{\frac EC}.
$$

The free-ion emissivity instead approaches the finite $J_0$ above.

A representative calculation using the full cold tensor, normal field and incidence, negligible damping, $B=1.5\times10^{13}$ G, and $\eta_\rho=0.55$ gives:

$$
\rho_s=6.31344\times10^4\ {\rm g\,cm^{-3}},
$$

$$
E_{ci}=44.2232\ {\rm eV},\quad
E_{pe}=4.93081\ {\rm keV},\quad
E_{pi}=78.7449\ {\rm eV}.
$$

| Surface-frame photon energy | Free-ion $J$ | Fixed-ion $J$ |
|---|---:|---:|
| $2$ eV | 0.882649 | 0.189681 |
| $8$ eV | 0.879772 | 0.306122 |

The free/fixed ratio at $2$ eV is $4.65333$. These are conditional patch radiances, not a disk-integrated RX J1856 prediction.

Nor is one ion prescription physically established. The source estimates the lattice scale by $\omega_D\simeq0.4E_{pi}/\hbar$ and says magnetic forces are “not completely” negligible compared with lattice forces. [refs/van_adelsberg_2005_condensed_surface.tex:497](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:497), [refs/van_adelsberg_2005_condensed_surface.tex:508](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:508).

For this example,

$$
\frac{\hbar\omega_D}{E_{ci}}=0.71225.
$$

Neither neglecting lattice binding nor freezing all ionic motion follows from a large separation of these scales.

The papers expect the actual spectrum to lie between the two prescriptions, but that expectation is **not a mathematical bound** on a dielectric tensor with additional lattice resonances. The explicit qualification “future work is needed” remains relevant to the proposed ownership route. [refs/van_adelsberg_2005_condensed_surface.tex:523](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:523), [refs/van_adelsberg_2005_condensed_surface.tex:527](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:527).

### 3.4 Damping is another approximation, not a solved leaf

The adopted low-frequency prescription is

$$
\nu_\parallel=\tau_\parallel^{-1},\qquad
\nu_\perp=\tau_\perp^{-1},
$$

with the exact anchor `\nu^{(l)}_{ei}=1/\tau_\|,~\nu^{(tr)}_{ei}=1/\tau_\perp`. It substitutes zero-frequency conductivity calculations for the unknown general finite-frequency response. [refs/van_adelsberg_2005_condensed_surface.tex:557](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:557), [refs/van_adelsberg_2005_condensed_surface.tex:565](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:565).

Potekhin 1999 is now on disk. Its $\tau_\parallel$ depends on a kinetic equation, and its $\tau_\perp$ on Landau-state sums. Anchors are “a kinetic equation for the electron spin density matrix” and `\Psi(\epsilon) = \sum_{nn'\gamma}`. [refs/potekhin_1999_magnetic_conductivity.tex:789](/home/tobias/Projects/neutronstars/refs/potekhin_1999_magnetic_conductivity.tex:789), [refs/potekhin_1999_magnetic_conductivity.tex:947](/home/tobias/Projects/neutronstars/refs/potekhin_1999_magnetic_conductivity.tex:947).

Two qualifications must survive implementation:

- $\tau_\perp$ is the effective collision time, **not** the transverse conductivity component $\tau_{xx}$; the source explicitly relates them through the Hall denominator. [refs/potekhin_1999_magnetic_conductivity.tex:968](/home/tobias/Projects/neutronstars/refs/potekhin_1999_magnetic_conductivity.tex:968), anchor `\tau_{xx}={\tau_\perp\over 1+(\omg\tau_\perp)^2}`.
- Inserting an arbitrary frequency-dependent damping rate into both real and imaginary parts of a Drude expression need not satisfy the Kramers–Kronig relation. PLCH04 explicitly says the transformed real part “does not coincide with the real” part obtained using the same damping function. [refs/potekhin_lai_chabrier_ho_2004_polarization.tex:539](/home/tobias/Projects/neutronstars/refs/potekhin_lai_chabrier_ho_2004_polarization.tex:539).

A direct DC transport calculation would remove one fit dependency. It would not establish an exact optical conductivity.

Vacuum polarization is a separate approximation. The local justification is that condensed density greatly exceeds the vacuum-resonance estimate

$$
\rho_V\simeq0.96(A/Z)B_{14}^2(E/{\rm keV})^2
\ {\rm g\,cm^{-3}}.
$$

Anchor: `\rho_V\simeq 0.96 (A/Z) B_{14}^2`. [refs/van_adelsberg_2005_condensed_surface.tex:751](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:751).

At the example above and $E=10$ keV, Python gives $\rho_V=4.65231\ {\rm g\,cm^{-3}}$ and $\rho_s/\rho_V=1.35706\times10^4$. This supports neglecting vacuum corrections **inside that condensate calculation**; it does not justify dropping them throughout an atmosphere.

## 4. Condensation: the relevant temperature is the interface temperature

### 4.1 Why the dielectric tensor cannot determine the phase

**DERIVED HERE:** transfer a small number of particles between vapor and condensate at fixed temperature and volume. The free-energy variation is

$$
\delta F=(\mu_{\rm vapor}-\mu_{\rm condensed})\delta N.
$$

Coexistence requires equal chemical potentials, together with mechanical equilibrium.

For an illustrative ideal, neutral, monatomic vapor,

$$
F_{\rm vapor}
=Nk_BT\left[
\ln\left(\frac{n\lambda_T^3}{g_{\rm int}}\right)-1
\right],
\qquad
\lambda_T=\sqrt{\frac{2\pi\hbar^2}{m_i k_BT}}.
$$

Taking $\mu_{\rm condensed}\simeq-Q_s$ gives

$$
n_{\rm vapor}\simeq
\frac{g_{\rm int}}{\lambda_T^3}
\exp\left(-\frac{Q_s}{k_BT}\right).
$$

This derives the strong sensitivity to cohesive energy. It does **not** calculate a magnetized vapor containing atoms, ions, and molecular chains.

Potekhin’s local account says Medin–Lai obtained the critical temperature by “equating the vapor density to $\rho_\mathrm{s}$” and found

$$
T_{\rm crit}\simeq0.08\,Q_s/k_B.
$$

The factor $0.08$ corresponds to $Q_s/(k_BT_{\rm crit})=12.5$; it is a reported numerical result, not a consequence of the cold-plasma tensor. [refs/potekhin_2012_condensed_surface.tex:351](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:351), [refs/potekhin_2012_condensed_surface.tex:366](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:366).

**Missing primary evidence:** the Medin–Lai papers themselves. The local account is sufficient to define a labelled comparison model, but not to reconstruct their cohesive-energy calculation.

For iron, Potekhin 2012 reports

$$
T_{\rm crit}\simeq(5+2B_{13})10^5\ {\rm K},
\qquad 0.5\leq B_{13}\leq100.
$$

The anchor is `$T_\mathrm{crit} \approx (5+2\,B_{13})\times10^5$~K`. [refs/potekhin_2012_condensed_surface.tex:368](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:368).

This is a representation of published calculations. Making $T_{\rm crit}$ a free parameter does not turn the cohesive-energy dependency into a measured or conventional leaf.

### 4.2 What RX J1856 actually constrains

The temperatures in the local papers are model dependent.

Ho’s atmosphere fit gives

$$
T_\infty=4.34(2)\times10^5\ {\rm K},\qquad z=0.22(2)
$$

for its $B=4\times10^{12}$ G model. The table also fixes $y_H=1.2\ {\rm g\,cm^{-2}}$. Anchors are `4.34(3) & 4.34(2)` and `$\ycol$ (g cm$^{-2}$) & 1.2 & 1.2`. [refs/ho_2007_rxj1856_thin_atm.tex:628](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:628), [refs/ho_2007_rxj1856_thin_atm.tex:635](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:635).

Using its printed relation $T_\infty=T_{\rm eff}/(1+z)$ gives

$$
T_{\rm eff}=529\,480\ {\rm K},
\qquad k_BT_{\rm eff}=45.6271\ {\rm eV}.
$$

Source anchor: `$\tinfty=\teff/(1+\zg)$`. [refs/ho_2007_rxj1856_thin_atm.tex:128](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:128).

The conversion uses the locally printed CODATA values `1.380 649 × 10−23` and `1.602 176 634 × 10−19`. [refs/tiesinga_2021_codata2018.txt:3731](/home/tobias/Projects/neutronstars/refs/tiesinga_2021_codata2018.txt:3731), [refs/tiesinga_2021_codata2018.txt:3426](/home/tobias/Projects/neutronstars/refs/tiesinga_2021_codata2018.txt:3426).

But Ho’s **condensed layer** is hotter:

$$
T_s\simeq7\times10^5\ {\rm K},
\qquad k_BT_s=60.3213\ {\rm eV}.
$$

The paper explicitly states `\tatm\approx 7\times 10^5`. [refs/ho_2007_rxj1856_thin_atm.tex:464](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:464).

This matters immediately. Ho quotes an older condensation estimate near $5.5\times10^5$ K for that field, with substantial uncertainty. Re-evaluating its printed expression gives $550\,584.5$ K. The nominal substrate temperature exceeds it. [refs/ho_2007_rxj1856_thin_atm.tex:438](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:438), anchor `$T\la 10^{5.5}B_{12}^{2/5}`.

The later Potekhin fit starts at $5\times10^{12}$ G; applying it at Ho’s $4\times10^{12}$ G would be an extrapolation.

Timing supplies a different field estimate:

$$
B_{\rm dip}\simeq1.5\times10^{13}\ {\rm G},
$$

explicitly inferred by “Assuming the star is spinning down by magnetic dipole radiation”. [refs/vankerkwijk_kaplan_2008_rxj1856_timing.tex:425](/home/tobias/Projects/neutronstars/refs/vankerkwijk_kaplan_2008_rxj1856_timing.tex:425), [refs/vankerkwijk_kaplan_2008_rxj1856_timing.tex:434](/home/tobias/Projects/neutronstars/refs/vankerkwijk_kaplan_2008_rxj1856_timing.tex:434).

At that field, the Potekhin representation gives

$$
T_{\rm crit}=800\,000\ {\rm K},
\qquad
k_BT_{\rm crit}=68.9387\ {\rm eV}.
$$

This agrees at the intended approximate level with the timing paper’s statement `$kT_{\rm cond}\simeq70\,$eV`. [refs/vankerkwijk_kaplan_2008_rxj1856_timing.tex:457](/home/tobias/Projects/neutronstars/refs/vankerkwijk_kaplan_2008_rxj1856_timing.tex:457).

As a **conditional comparison**, carrying Ho’s temperatures to this field gives

$$
T_{\rm eff}/T_{\rm crit}=0.66185,\qquad
T_s/T_{\rm crit}=0.875.
$$

The second ratio is the relevant one for the substrate. Even the quoted $20\%$ variation around the temperature representation spans $640\,000$–$960\,000$ K, straddling the assumed substrate temperature.

That comparison is not a self-consistent fit: Ho’s temperature structure was computed at a different field. It demonstrates why a condensation label cannot be assigned before solving the atmosphere.

The alternative X-ray blackbody fit gives $k_BT_\infty=63.5\pm0.2$ eV. [refs/burwitz_2003_rxj1856_xray.tex:320](/home/tobias/Projects/neutronstars/refs/burwitz_2003_rxj1856_xray.tex:320), anchor `63.5\pm 0.2`.

Python converts this to $736\,887\pm2321$ K. Applying the illustrative redshift $z=0.22$ would give a local color temperature $899\,002$ K, or $77.47$ eV—above the nominal condensation threshold. This is another model-dependent comparison, not a determination of the material temperature.

**Conclusion:** iron condensation is possible, but neither established nor securely in the negligible-vapor regime. The cold-plasma calculation cannot settle this.

Also distinguish material temperature from effective temperature even for a bare surface. **DERIVED HERE** from the emissivity law,

$$
F=\int dE\int_{\rm outward}d\Omega\,
\mu J(E,\Omega)B_E(T_s).
$$

Defining $F=\sigma T_{\rm eff}^4$, a surface with suppressed emissivity generally has $T_s>T_{\rm eff}$. Equating the two silently biases the phase test.

## 5. A thin atmosphere changes the boundary problem

### 5.1 Column depth and transfer

Define downward mass column

$$
y(z)=\int_z^\infty \rho(z')\,dz',
\qquad 0\leq y\leq y_H.
$$

**DERIVED HERE:** in plane-parallel hydrostatic balance with approximately constant gravity, $dP/dy=g$. For an opacity per mass $\kappa_E$,

$$
d\tau_E=\kappa_E\,dy.
$$

For an uncoupled ray with source function $S_E$, the outward transfer equation and formal solution are

$$
\mu\frac{\partial I_E}{\partial y}
=\kappa_E(I_E-S_E),
$$

$$
I_E(0,\mu)
=
I_E(y_H,\mu)e^{-\tau_E(y_H)/\mu}
+\int_0^{\tau_E(y_H)}
S_E(t,\mu)e^{-t/\mu}\frac{dt}{\mu}.
$$

This is the direct integration of the first-order transfer equation. Ho’s local optical-depth definition is `\taulam=\int\chi_\lambda\,ds`, and its thin-atmosphere discussion explicitly distinguishes low-energy optically thick from high-energy optically thin radiation. [refs/ho_2007_rxj1856_thin_atm.tex:316](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:316), [refs/ho_2007_rxj1856_thin_atm.tex:330](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:330).

Thus $y_H$ sets the lower integration boundary and affects both the optical depth and the temperature structure. It is not a grey transmission coefficient.

### 5.2 Correct lower boundary

For incoherent atmospheric modes $M\in\{X,O\}$,

$$
\boxed{
I_M^+(E,\Omega)
=
\frac12J_M(E,\Omega)B_E(T_s)
+\sum_{M'}\mathcal R_{MM'}(E,\Omega)I_{M'}^-(E,\Omega_{\rm spec}).
}
$$

Here

$$
J_M=1-\sum_{M'}\mathcal R_{MM'}.
$$

These are Potekhin’s printed inner-boundary equations, including the cross-mode term `+ R_\mathrm{XO}`. [refs/potekhin_2012_condensed_surface.tex:1067](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:1067).

**DERIVED HERE:** if $U_i$ and $U_r$ contain the incident and reflected atmospheric polarization vectors expressed in the vacuum linear bases,

$$
r_{\rm mode}=U_r^{-1}r_{\rm linear}U_i,
\qquad
\mathcal R_{MM'}=|(r_{\rm mode})_{MM'}|^2,
$$

provided the exterior-mode normalization and nearly transverse approximation justify this conversion.

For a general refractive atmosphere, solve the boundary problem using that atmosphere’s own dielectric tensor and flux-normalized modes. The substrate density and the gas density immediately above it are distinct inputs.

Potekhin explicitly says the atmosphere interface differs from the vacuum interface and permits replacement only when the atmosphere is “sufficiently rarefied”. [refs/potekhin_2012_condensed_surface.tex:2063](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:2063).

This approximation must be checked at the optical/UV energies being rendered, not inferred from its usefulness for X-rays.

The intensity-only modal description also has a domain restriction. Ho–Lai give

$$
|\operatorname{Re}(n_+-n_-)|
\gg
|\operatorname{Im}(n_++n_-)|
$$

as the Faraday-depolarization condition. [refs/ho_lai_2003_atmospheres_vacuum.tex:958](/home/tobias/Projects/neutronstars/refs/ho_lai_2003_atmospheres_vacuum.tex:958), anchor “valid only in the limit of Faraday depolarization”.

Where it fails, use coupled polarized transfer or coherency evolution; do not silently treat nonorthogonal, mutually coherent modes as two independent intensities.

### 5.3 What reproducing Ho 2007 would—and would not—verify

Potekhin 2012 states that Ho 2007 “retained only the emission terms” in the boundary equations. [refs/potekhin_2012_condensed_surface.tex:1104](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:1104).

Therefore separate two tests:

- **Historical reproduction:** recover Ho-like spectra under Ho-like assumptions.
- **Corrected computation:** include reflected incident radiation and mode conversion at the interface.

Agreement with the first does not validate the second, and the corrected calculation need not return identical fitted parameters.

Partial ionization is also substantive. Ho states that its opacities differ enough from “fully ionized opacities” to change the structure and continuum. [refs/ho_2007_rxj1856_thin_atm.tex:167](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:167).

The boundary Maxwell solve and transfer integration are derivable. Ho’s fitted $(T_\infty,R_\infty,z,B,y_H)$, the condensed-density representation, cohesive-energy criterion, and imported opacity tables are not independently owned merely because the surrounding solver is new.

## 6. Numerical verification: what the local papers actually print

### 6.1 Two explicitly printed emissivity targets

The requested papers do not provide two high-precision tabulated emissivities in the examined passages. They do provide two explicit approximate values, which can be tested without inventing precision.

All calculations below use free ions, zero damping interpreted as the limiting-absorption result, $B=10^{13}$ G, $\eta_\rho=1$, $A=56$, and $Z=26$.

| Target | Printed expected value and anchor | Computed value | Absolute error | Proposed tolerance |
|---|---|---:|---:|---:|
| Normal incidence near the upper ion-related suppression edge | $J\simeq0.5$; `J\rightarrow 0.5` | 0.4893078435 | 0.0106921565 | 0.02 |
| Oblique incidence near the plasma feature, without damping | $J\simeq0.5$; “in the absence of collisional damping” | 0.4999999675 | $3.25\times10^{-8}$ | 0.005 |

The first target is printed at [refs/van_adelsberg_2005_condensed_surface.tex:890](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:890); the second at [refs/van_adelsberg_2005_condensed_surface.tex:907](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:907).

The evaluation points are specified, rather than selected by eye:

- **First target:** solve the full zero-damping circular-mode cutoff,

  $$
  (E+E_{ce})(E-E_{ci})=E_{pe}^2.
  $$

  The positive root is $E_*=0.2636181546$ keV. Evaluate at $E=(1-10^{-6})E_*$ and $\theta=0$.

- **Second target:** $E=(1+10^{-6})E_{pe}=5.212938424$ keV and $\theta=\pi/3$.

The tolerances reflect the approximate source statements. The much smaller numerical residual of the second calculation does not make the published statement an eight-decimal datum.

An additional figure target is available. Visual inspection of Figure 3 gives approximately $J=0.75$ at $E=0.001$ keV for normal incidence. Using V05’s own density normalization and $A=55.9$, Python gives $J=0.7514823$, an absolute difference of $0.0014823$. A figure-reading tolerance of $0.02$ is reasonable. The caption anchor is “Same as Fig.~2, except for”. [refs/van_adelsberg_2005_condensed_surface.tex:849](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:849). This is a **figure-derived** target, not a tabulated datum.

### 6.2 A source sign check worth preserving

**DERIVED HERE:** set $\epsilon_\perp=1$ and $g=0$, retaining $\epsilon_\parallel=1-v_e$. The extraordinary-mode determinant gives

$$
q_p^2=1-\frac{s^2}{1-v_e},
\qquad
n_p^2=1-\frac{v_es^2}{1-v_e}.
$$

V05’s approximate equation instead prints `1+{v_e\over {(1-v_e)}}`. [refs/van_adelsberg_2005_condensed_surface.tex:901](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:901).

The sign in that displayed approximation is inconsistent with the reduction of the stated wave equation. Use the dielectric tensor and Maxwell determinant as the numerical oracle, not that approximate index formula. The approximate half-emissivity target remains usable.

### 6.3 Checks performed and implementation hazards found

The exploratory calculation compared the four-component boundary solution against the independent circular-index formula at normal incidence:

| Energy [keV] | Boundary calculation $J$ | Circular formula $J$ | Absolute difference |
|---:|---:|---:|---:|
| 0.002 | 0.8826492024850205 | 0.8826492024850208 | $2.22\times10^{-16}$ |
| 0.008 | 0.8797719101213375 | 0.8797719101213375 | 0 |
| 0.1 | 0.4858820267029624 | 0.4858820267029624 | $5.55\times10^{-17}$ |
| 1.0 | 0.9987273971172512 | 0.9987273971172512 | 0 |

These use $B=1.5\times10^{13}$ G and $\eta_\rho=0.55$. A synthetic isotropic dielectric $\epsilon=4$ gives the analytic result $J=8/9=0.8888888889$.

A reproducible exploratory sweep of 300 cases varied energy, field, orientation, density factor, and ion prescription. With an explicitly synthetic damping energy $\hbar\nu=0.1$ eV, it gave:

- Maximum relative boundary residual: $4.62\times10^{-16}$.
- Maximum normalized eigenproblem residual: $1.63\times10^{-15}$.
- Minimum eigenvalue of $\mathbf I-rr^\dagger$: $2.37\times10^{-9}$.

This is a numerical consistency check, not validation of that damping rate.

Two failures are particularly instructive:

1. Naively expanding the three-dimensional determinant produced $J=-5.9996$ in one case despite a boundary residual of $3.31\times10^{-16}$. The mode calculation was wrong; the boundary solver faithfully matched the wrong modes.
2. Using $\mathbf I-r^\dagger r$ as the emission matrix produced a blackbody-cavity residual of $0.00672846$ in an oblique example. A row/column mutation is therefore detectable.

Production assertions must cover the dielectric’s passivity, mode residuals, propagation direction, boundary residuals, and

$$
0\preceq\mathbf I-rr^\dagger\preceq\mathbf I.
$$

A small linear-system residual alone is inadequate evidence.

## 7. Concrete implementation plan

Implement this as a conditional physical model with explicit provenance on every material input.

| Function | Inputs | Outputs and verification |
|---|---|---|
| `condensed_material_state` | $B,T_s,Z,A$, density/phase prescription or computed EOS | $\rho_s$, phase assessment, provenance, validity interval. Reject extrapolation. At $B=4\times10^{12}$ G, $\eta_\rho=1$, compute $23\,499.90\ {\rm g\,cm^{-3}}$ versus Ho’s printed $2.35\times10^4$; relative discrepancy $4.25\times10^{-6}$. |
| `condensed_collision_rates` | $B,T_s,\rho_s,Z,A$, transport model | $\nu_\parallel,\nu_\perp$, approximation flags. Separate effective collision times from tensor conductivity components. |
| `cold_condensed_dielectric` | $E,B,\rho_s,Z,A,\nu_\parallel,\nu_\perp$, ion prescription | Complex field-frame tensor. Verify zero-field, circular-eigenvalue, fixed-ion, and positive-dissipation limits. |
| `rotate_dielectric` | Tensor, unit field vector | Surface-frame tensor. Verify rotational covariance and field-reversal behavior, including Stokes signs. |
| `transmitted_mode_subspace` | Tensor, tangential wavevector | Physical invariant subspace, normal wavevectors, residuals, inward-flux diagnostics. Handle repeated modes and singular eliminations explicitly. |
| `interface_reflection` | Incident/reflected exterior modes, transmitted subspace | Complex $r$, transmission amplitudes, flux diagnostics. Compare with the circular formula and both printed half-emissivity targets. |
| `condensed_emission_stokes` | $r,E,T_s$, basis convention | Absolute $(I,Q,U,V)$ from $\mathbf I-rr^\dagger$. Verify equilibrium-cavity closure and positivity. |
| `atmosphere_bottom_boundary` | $r$, atmospheric modes, downward radiation, $T_s$ | Emitted-plus-reflected upward radiation, including cross-mode terms. |
| `solve_thin_atmosphere` | $y_H,g,B$, composition, flux, boundary operator | Converged temperature structure, mode/coherency radiation, substrate temperature, emergent Stokes radiance. Recheck phase consistency after convergence. |
| `earth_photometry` | Emergent radiance, surface map, GR transport, distance, extinction, measured passbands | Predicted band fluxes and residuals against measured data; material-uncertainty propagation. |

The Ho density target is printed with `\rho=2.35\times 10^4`. [refs/ho_2007_rxj1856_thin_atm.tex:473](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:473).

A further dimensional target is Potekhin’s printed iron cyclotron energy at $B=4\times10^{13}$ G: `\Eci =0.118`. The calculated value is $0.117928571$ keV, error $7.14\times10^{-5}$ keV; tolerance $5\times10^{-4}$ keV respects the printed precision. [refs/potekhin_2012_condensed_surface.tex:1277](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:1277).

Follow red-green TDD in this order:

1. Scalar Fresnel and repeated-mode tests.
2. Circular-mode analytic comparisons.
3. Arbitrary-orientation residual, passivity, and basis-covariance tests.
4. Blackbody-cavity and full-Stokes tests.
5. Reflection-capable atmosphere boundary tests.
6. Published atmosphere comparisons and photometry.

Mutation tests should reverse the damping sign, choose the wrong normal-wavevector branch, interchange row and column emissivities, omit cross-polarization amplitudes, and duplicate a repeated-mode eigenvector. Fuzzing should concentrate around cyclotron resonances, plasma features, grazing incidence, nearly tangent fields, and eigenvalue coalescence.

Use a numerical residual target such as $10^{-10}$ in well-conditioned cases, separately from the much looser physical comparison tolerances. Report conditioning and fail with the offending parameters when that numerical target is unattainable.

## 8. Missing sources to fetch

The following bibliography is **unverified against the original works**. Identifiers are retained only where local papers or the project’s gap list supply them.

| Missing source | Bibliographic lead, marked unverified | Purpose |
|---|---|---|
| Medin & Lai 2006a | **Unverified:** *Physical Review A* 74, 062507; astro-ph/0607166 | Atomic/molecular calculations needed to reconstruct cohesive-energy comparisons. |
| Medin & Lai 2006b | **Unverified:** *Physical Review A* 74, 062508; astro-ph/0607277 | Condensed/infinite-chain calculations and equilibrium structural parameters. Verify the pairing of article number and identifier during acquisition. |
| Medin & Lai 2007 | **Unverified:** *MNRAS* 382, 1833; arXiv:0708.3863 | Saturated vapor, cohesive energies, and condensation criterion. |
| Chabrier 1993 | **Unverified:** *ApJ* 414, 695; arXiv identifier unknown | Coulomb-crystal vibrational physics underlying the quoted lattice scale. |
| Usov, Grebenshchikov & Ulinich 1980 | **Unverified:** *Soviet Physics JETP* 51, 148; arXiv identifier unknown | Magnetic modification of ionic lattice dynamics. |
| Ziman 1979 | **Unverified:** *Principles of the Theory of Solids*, Cambridge University Press; no arXiv identifier supplied | Bound-ion dielectric response and the additional restoring-force terms. |
| Potekhin & Pavlov 1997 | **Unverified:** *ApJ* 483, 414; arXiv identifier unknown | Bound-free atomic opacity needed for the partially ionized atmosphere branch. |

The two 2006 arXiv leads are printed as `astro-ph/0607166` and `astro-ph/0607277` in Ho’s bibliography. [refs/ho_2007_rxj1856_thin_atm.tex:1141](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:1141). Potekhin’s bibliography gives `\pra, 74, 062508` and `\mnras, 382, 1833`. [refs/potekhin_2012_condensed_surface.tex:1646](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:1646). The 2007 identifier appears in the timing paper’s bibliography. [refs/vankerkwijk_kaplan_2008_rxj1856_timing.tex:540](/home/tobias/Projects/neutronstars/refs/vankerkwijk_kaplan_2008_rxj1856_timing.tex:540).

The lattice leads are locally listed as “Chabrier, G. 1993, ApJ, 414, 695”, “Sov.\ Phys.--JETP, 51, 148”, and “Ziman, J.M. 1979, Principles of the Theory of Solids”. [refs/van_adelsberg_2005_condensed_surface.tex:1459](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:1459), [refs/van_adelsberg_2005_condensed_surface.tex:1671](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:1671), [refs/van_adelsberg_2005_condensed_surface.tex:1688](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:1688). The atomic-opacity lead is “ApJ, 483, 414”. [refs/potekhin_lai_chabrier_ho_2004_polarization.tex:1264](/home/tobias/Projects/neutronstars/refs/potekhin_lai_chabrier_ho_2004_polarization.tex:1264).

Fetching these sources does not guarantee a ready-made, validated optical dielectric tensor for the actual material. That remains a possible research task. Potekhin 1999, PLCH04, and Ho–Lai 2003 are already local and should not remain listed as missing.

## 9. Stub policy, final-render impact, and R5

### 9.1 What is defensible for a hero render

A labelled development render may use the Potekhin emissivity fits within their stated iron, density, energy, field, and geometry assumptions. The source gives a fit range of approximately $10^{12}$–$10^{14}$ G and $1$ eV–$10$ keV. [refs/potekhin_2012_condensed_surface.tex:685](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:685), anchor `\mbox{1~eV}\lesssim E`.

Its reported error is within $10\%$ over most of the logarithmic energy interval, with larger errors for some inclined geometries. Anchor: “geometric settings, the fit error lies within 10\%”. [refs/potekhin_2012_condensed_surface.tex:896](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:896).

That is **fit-to-model error**, not total physical uncertainty.

A defensible ensemble should distinguish:

- Free-ion and fixed-ion material closures.
- Density prescriptions and their validity ranges.
- Bare and atmosphere-covered surfaces.
- Phase-consistent cases from hypothetical condensed cases.
- Intensity-only previews from full-Stokes outputs.

The fit’s prescription to truncate unphysical mode reflectivities cannot be copied into a fail-loud production path. The paper says “one should truncate the mode-specific” reflectivities. [refs/potekhin_2012_condensed_surface.tex:2217](/home/tobias/Projects/neutronstars/refs/potekhin_2012_condensed_surface.tex:2217). Such a point should instead be flagged and recomputed with the direct boundary solver.

Surface roughness is an additional, independent uncertainty. The source’s geometrical multiple-reflection example makes emission closer to blackbody when roughness scales exceed wavelength. [refs/van_adelsberg_2005_condensed_surface.tex:1164](/home/tobias/Projects/neutronstars/refs/van_adelsberg_2005_condensed_surface.tex:1164), anchor “Effect of surface roughness on the reflectivity.” There is no basis here for turning that example into an arbitrary visual texture or an adjustable grey emissivity.

### 9.2 Estimated optical/UV effect

For the bare-patch example above, choosing fixed rather than free ions reduces the $2$ eV radiance by a factor $4.65333$, equivalent to $1.66941$ magnitudes.

The conversion follows by subtracting the locally printed magnitude definition, `STMAG = −2.5 log10 Fλ − 21.1`. [refs/kaplan_2011_rxj1856_spectrum.txt:345](/home/tobias/Projects/neutronstars/refs/kaplan_2011_rxj1856_spectrum.txt:345).

Between $2$ and $8$ eV, the ratio of fixed-ion to free-ion spectral-shape changes is $1.61916$, equivalent to a color difference of $0.52322$ magnitudes at those two sampled energies. This is not grey dimming.

**DERIVED HERE:** in the Rayleigh–Jeans limit, $B_\nu\propto\nu^2$. The low-energy free-ion limit therefore retains that slope, whereas the fixed-ion asymptote $J\propto E^{1/2}$ gives $I_\nu\propto\nu^{5/2}$. The renderer must use the full Planck function, but the asymptote explains why ionic response can alter visible/UV color.

For an absorbing atmosphere at fixed source function, the direct boundary contribution is suppressed by $e^{-\tau/\mu}$; a factor below $0.01$ requires $\tau/\mu>4.60517$. This is only a direct-transmission estimate. Scattering, cross-mode reflection, and changes in the temperature structure can transmit substrate dependence indirectly.

Ho’s particular model states that the atmosphere is optically thick at the longer wavelengths where condensed emission differs from blackbody, so the substrate features are hidden. [refs/ho_2007_rxj1856_thin_atm.tex:484](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:484), anchor “the atmosphere is optically thick at longer wavelengths”.

That result supports testing a screened-substrate approximation. It does not authorize assuming screening for every field, column, mode, or phase-consistent atmosphere.

### 9.3 Consequence for R5

A typical $10\%$ radiance error corresponds to approximately $0.10$–$0.11$ magnitudes; a $20\%$ error corresponds to approximately $0.20$–$0.24$ magnitudes. The bare-surface ion ambiguity in the worked example is much larger.

For scale, Walter’s locally printed F475W result is `25.36$\pm$0.15`. [refs/walter_2010_rxj1856_distance.tex:862](/home/tobias/Projects/neutronstars/refs/walter_2010_rxj1856_distance.tex:862). Its magnitude uncertainty corresponds to approximately $+14.8\%/-12.9\%$ in flux. This single measurement is not the complete R5 dataset, and its passband and zero-point handling must be reproduced.

Consequently:

- **Bare condensed render:** the current material uncertainty can change optical brightness by factors and UV/optical color appreciably. A precise R5 claim is unjustified without propagating it.
- **Thin atmosphere:** a substrate stub may be photometrically adequate if a coupled recalculation shows that varying the substrate assumptions changes every required band by less than its allocated error budget. This must be demonstrated for the converged atmosphere, including its interface temperature.
- **Project-wide R5:** this path need not block R5 if another permitted surface model passes the measured optical/UV constraints. R5 requires at least one model; it does not require proving that RX J1856 is condensed.

The minimal honest endpoint is therefore an owned Maxwell solver with full complex amplitudes and Stokes emission, coupled to a reflection-capable atmosphere, while retaining explicit material and phase uncertainty. A hero image may be shown earlier as a labelled ensemble member. It cannot yet be presented as a uniquely determined condensed surface or as an R5-certified prediction.