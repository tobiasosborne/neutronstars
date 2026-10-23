## Abstract

Interstellar extinction is small enough that it cannot explain RX J1856’s optical excess, but large enough to change a marginal R5 verdict. The calculations below give a change of approximately $0.0665$ mag in WFPC2/F606W between $A_V=0.05$ and $0.12$, compared with a locally reported photometric uncertainty equivalent to $0.0326$ mag. For a model initially centred on that measurement, an extinction uncertainty of approximately $0.034$ mag is therefore consequential.

Two premises of the requested calculation need correction. Kaplan tabulates **effective wavelengths, not pivot wavelengths**, and Verner (1996) does **not** supply all the inner-shell data needed to calculate the complete soft-X-ray ISM opacity. I calculate the supported optical estimates, the Verner H/He contribution, Morrison–McCammon transmissions, and a separately labelled numerical benchmark extracted from Wilms’s local figure. A complete Wilms/Verner elemental-sum evaluation remains blocked by missing atomic inputs.

For M2, the defensible policy is an explicitly yellow extinction/absorption branch, with independently varied $A_V$ and $N_H$, actual passband integration, and separate reporting of observational residuals and model uncertainty. Recomputing grain scattering or atomic cross-sections does not determine the material composition along this sightline.

*All numerical evaluations below were run in Python. Relevant scanned equations and tables were visually checked. The workspace is read-only, so this response does not update `HANDOFF.md`, the worklog, or the truth database.*

# Extinction and absorption as a critical path to R5

## 1. What the local evidence actually supports

The relevant observational bands are identified explicitly in [refs/ho_2007_rxj1856_thin_atm.tex:518](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:518):

> `F170W, F300W, F450W, and F606W photometry`

These are WFPC2 measurements. Ho additionally uses a STIS far-UV spectrum; that is not interchangeable with an integrated ACS/F140LP measurement.

Kaplan gives the following representative wavelengths:

- WFPC2/F170W: $1769$ Å, in [refs/kaplan_2011_rxj1856_spectrum.tex:706](/home/tobias/Projects/neutronstars/refs/kaplan_2011_rxj1856_spectrum.tex:706), whose row contains `1355 & 1442 & 1517 & 1769 & 2283`.
- WFPC2/F300W, F450W and F606W: $2955$, $4520$ and $5932$ Å, respectively, in [refs/kaplan_2011_rxj1856_spectrum.tex:751](/home/tobias/Projects/neutronstars/refs/kaplan_2011_rxj1856_spectrum.tex:751), with the continuation `4520 & 4709 & 5066 & 5844 & 5932 & 7025`.

However, Kaplan explicitly calls them:

> `approximate effective`

at [refs/kaplan_2011_rxj1856_spectrum.tex:716](/home/tobias/Projects/neutronstars/refs/kaplan_2011_rxj1856_spectrum.tex:716), and specifies that they depend on a typical neutron-star spectrum. **Calling these pivot wavelengths would introduce a false factual claim.** The requested pivot-wavelength calculation cannot be completed from these sources as stated. The calculations below use the printed effective wavelengths and label the approximation.

A second trap is more serious for R5. Kaplan’s tables containing the RX J1856 entries with $\pm0.05$ mag errors are headed:

> `Predicted \hst\ Optical Photometry from Extrapolated X-ray Blackbodies`

[refs/kaplan_2011_rxj1856_spectrum.tex:730](/home/tobias/Projects/neutronstars/refs/kaplan_2011_rxj1856_spectrum.tex:730). Their errors are described as:

> `uncertainties in the`

> `X-ray spectra`

at lines 759–760. These are **not measured HST photometric errors** and must not become R5 tolerances.

Two usable local error statements are:

1. Ho’s WFPC2/F606W discussion gives `3\% (photometric uncertainty)` at [refs/ho_2007_rxj1856_thin_atm.tex:711](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:711).
2. Walter reports `25.36$\pm$0.15` for F475W at [refs/walter_2010_rxj1856_distance.tex:862](/home/tobias/Projects/neutronstars/refs/walter_2010_rxj1856_distance.tex:862). This is ACS/HRC photometry; Kaplan’s ACS/WFC effective wavelength and predictions cannot silently be substituted for its calibration.

The individual measured WFPC2 magnitudes and uncertainties still need harvesting from the original photometric analysis. Moreover, `truth/GAPS.md` §2.1 names the wrong van Kerkwijk–Kulkarni paper: Kaplan’s bibliography gives `2001, \aap, 378, 986` at [refs/kaplan_2011_rxj1856_spectrum.tex:1371](/home/tobias/Projects/neutronstars/refs/kaplan_2011_rxj1856_spectrum.tex:1371), and Ho independently gives `2001a, A\&A, 378, 986` at [refs/ho_2007_rxj1856_thin_atm.tex:1204](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:1204).

## 2. Optical/UV attenuation: derive the observable before choosing a dust model

### 2.1 The CCM89 prescription

CCM89 prints

$$
\left\langle \frac{A(\lambda)}{A(V)}\right\rangle
=a(x)+\frac{b(x)}{R_V},
\qquad x=\lambda^{-1},
$$

with $x$ expressed numerically in $\mu{\rm m}^{-1}$. Its Eq. 1 is at [refs/cardelli_1989_extinction.txt:336](/home/tobias/Projects/neutronstars/refs/cardelli_1989_extinction.txt:336), anchored by:

> `<A(X)/A(V)> = a(x) + b(x)/Rv`

The angle brackets matter: this is a fitted mean relation. The applied model is consequently

$$
A_\lambda=A_V k(\lambda;R_V),
\qquad
k=a+b/R_V.
$$

For the optical segment, set $y=x-1.82$:

$$
\begin{aligned}
a={}&1+0.17699y-0.50447y^2-0.02427y^3
 +0.72085y^4+0.01979y^5-0.77530y^6+0.32999y^7,\\
b={}&1.41338y+2.28305y^2+1.07233y^3
 -5.38434y^4-0.62251y^5+5.30260y^6-2.09002y^7.
\end{aligned}
$$

These are CCM89 Eqs. 3a–3b, visually checked on PDF page 5; the text anchor is `a(x) = 1 + 0.17699y` at [refs/cardelli_1989_extinction.txt:311](/home/tobias/Projects/neutronstars/refs/cardelli_1989_extinction.txt:311). The segment is $1.1\leq x\leq3.3$.

For $3.3\leq x\leq8$,

$$
\begin{aligned}
a(x)&=1.752-0.316x-\frac{0.104}{(x-4.67)^2+0.341}+F_a(x),\\
b(x)&=-3.090+1.825x+\frac{1.206}{(x-4.62)^2+0.263}+F_b(x).
\end{aligned}
$$

For $x<5.9$, $F_a=F_b=0$. Otherwise, writing $u=x-5.9$,

$$
F_a=-0.04473u^2-0.009779u^3,
\qquad
F_b=0.2130u^2+0.1207u^3.
$$

Sources: [refs/cardelli_1989_extinction.txt:431](/home/tobias/Projects/neutronstars/refs/cardelli_1989_extinction.txt:431), anchor `a(x) = 1.752 - 0.316x`; line 435, anchor `b(x) = -3.090 + 1.825x`; and lines 371–379, anchor `0.04473`. These were checked against PDF page 6 rather than trusting the OCR’s superscripts.

At the adopted wavelengths, F170W and F300W use the UV branch; F450W and F606W use the optical branch.

### 2.2 Numerical corrections at the locally printed wavelengths

Using $R_V=3.1$, the Python evaluations are:

| WFPC2 band | Printed effective wavelength, Å | Computed $k(\lambda)$ | $A_\lambda$ for $A_V=0.05$ | $A_\lambda$ for $A_V=0.12$ | Difference |
|---|---:|---:|---:|---:|---:|
| F170W | 1769 | 2.509791 | 0.125490 mag | 0.301175 mag | 0.175685 mag |
| F300W | 2955 | 1.843240 | 0.092162 mag | 0.221189 mag | 0.129027 mag |
| F450W | 4520 | 1.280654 | 0.064033 mag | 0.153678 mag | 0.089646 mag |
| F606W | 5932 | 0.918130 | 0.045907 mag | 0.110176 mag | 0.064269 mag |

These are **monochromatic estimates**, not synthetic band extinctions.

The sign follows directly from the magnitude definition. Kaplan prints `${\rm STMAG}=-2.5\log_{10} F_{\lambda}-21.1` at [refs/kaplan_2011_rxj1856_spectrum.tex:488](/home/tobias/Projects/neutronstars/refs/kaplan_2011_rxj1856_spectrum.tex:488). Subtracting the intrinsic and attenuated magnitudes gives, **DERIVED HERE**,

$$
A_\lambda=m_{\lambda,\rm obs}-m_{\lambda,0}
=-2.5\log_{10}\frac{F_{\lambda,\rm obs}}{F_{\lambda,0}},
$$

hence

$$
F_{\lambda,\rm obs}=F_{\lambda,0}\,10^{-0.4A_\lambda}.
$$

Thus extinction is added when predicting observed magnitudes and subtracted when dereddening observations. At F170W, the computed transmissions are $0.890849$ and $0.757757$ for the two $A_V$ choices; at F606W they are $0.958600$ and $0.903503$.

A normalization check gives:

| Expected value and source | Computed | Error | Numerical tolerance |
|---|---:|---:|---:|
| $k(x=1.82)=1$, CCM89 Table 3, [refs/cardelli_1989_extinction.txt:347](/home/tobias/Projects/neutronstars/refs/cardelli_1989_extinction.txt:347), anchor `1.0000` | 1.000000 | 0 | $10^{-12}$ |

This checks arithmetic and normalization, not the accuracy of CCM89 for this sightline.

### 2.3 Why a wavelength substitution is insufficient for R5

Let $W_i(\lambda)$ be the measured detector response kernel, including conversion from incident spectral energy flux to the recorded signal. Assume a linear detector and a foreground attenuation screen. Then the band signal is, **DERIVED HERE**,

$$
C_i(A_V)=\int W_i(\lambda)F_{\lambda,0}
\exp[-qA_Vk(\lambda)]\,d\lambda,
\qquad q=0.4\ln10.
$$

Define a normalized unextinguished weight,

$$
w_i(\lambda)=
\frac{W_i(\lambda)F_{\lambda,0}}
{\int W_i(\lambda)F_{\lambda,0}\,d\lambda}.
$$

Taking the ratio of signals eliminates the zero point:

$$
A_i=-2.5\log_{10}\frac{C_i(A_V)}{C_i(0)}
=-\frac{2.5}{\ln10}
\ln\left\langle e^{-qA_Vk}\right\rangle_{w_i}.
$$

Expanding the exponential and then the logarithm gives

$$
A_i
=A_V\langle k\rangle_{w_i}
-\frac{q}{2}A_V^2\operatorname{Var}_{w_i}(k)
+O(A_V^3).
$$

The computed coefficient is $q/2=0.460517$. Even at first order, the relevant coefficient is a spectrum-weighted average, not $k$ evaluated at a nominal wavelength.

Differentiating the unexpanded expression gives the useful exact sensitivity,

$$
\frac{\partial A_i}{\partial A_V}
=\langle k\rangle_{w_i,A_V},
$$

where the weight now includes extinction. This derivative, evaluated through the actual throughput, is what the R5 uncertainty calculation should use.

Kaplan supplies approximate band-integrated extinction coefficients as verification targets: `2.65` for F170W at [refs/kaplan_2011_rxj1856_spectrum.tex:707](/home/tobias/Projects/neutronstars/refs/kaplan_2011_rxj1856_spectrum.tex:707), and `1.93 & 1.30` together with `0.95` for the other WFPC2 bands at [refs/kaplan_2011_rxj1856_spectrum.tex:754](/home/tobias/Projects/neutronstars/refs/kaplan_2011_rxj1856_spectrum.tex:754).

Using those published coefficients gives this separate, approximately band-integrated estimate:

| WFPC2 band | Published $A_i/A_V$ | $A_i(0.05)$ | $A_i(0.12)$, linearized | Difference |
|---|---:|---:|---:|---:|
| F170W | 2.65 | 0.1325 | 0.3180 | 0.1855 |
| F300W | 1.93 | 0.0965 | 0.2316 | 0.1351 |
| F450W | 1.30 | 0.0650 | 0.1560 | 0.0910 |
| F606W | 0.95 | 0.0475 | 0.1140 | 0.0665 |

All entries are magnitudes except the coefficient column. The $0.12$ column is an extrapolation of the *coefficient*, explicitly labelled here, not a replacement for computing the band integral. Kaplan reports coefficient variations below `1\%` over the extinction range discussed at [refs/kaplan_2011_rxj1856_spectrum.tex:664](/home/tobias/Projects/neutronstars/refs/kaplan_2011_rxj1856_spectrum.tex:664).

This distinction already matters numerically: at F170W and $A_V=0.12$, the monochromatic and band estimates differ by $0.016825$ mag.

Ho explicitly acknowledges the same limitation:

> `convolving the models with the filter bandpasses`

[refs/ho_2007_rxj1856_thin_atm.tex:704](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:704).

The minimal honest M2 computation therefore integrates the extincted spectrum through the throughput. It does not require a Mie solver before it requires correct photometry.

### 2.4 What extinction uncertainty can flip R5?

There is no model-independent single answer: it depends on the residual before changing extinction.

For fixed surface parameters, define

$$
r_i=m_{i,\rm model}-m_{i,\rm obs},
$$

and let $\epsilon_i$ be the published magnitude uncertainty adopted by the R5 test. Locally,

$$
r_i(A_V+\delta A_V)
\simeq r_i(A_V)+k_{i,\rm eff}\delta A_V.
$$

The accepted interval is consequently, **DERIVED HERE**,

$$
\frac{-\epsilon_i-r_i}{k_{i,\rm eff}}
\leq\delta A_V\leq
\frac{\epsilon_i-r_i}{k_{i,\rm eff}}.
$$

For an initially passing model, the distance to its nearest boundary is

$$
|\delta A_V|_{\rm nearest}
=\frac{\epsilon_i-|r_i|}{k_{i,\rm eff}}.
$$

From Ho’s quoted F606W photometric error of $3\%$,

$$
\epsilon_{606}
\simeq\frac{2.5}{\ln10}(0.03)
=0.032572\ {\rm mag}.
$$

The exact bright/faint magnitude excursions are $0.032093$ and $0.033071$ mag. Using Kaplan’s approximate band coefficient $k_{606}=0.95$:

- A model centred on the measurement reaches the boundary after $|\delta A_V|\simeq0.034286$ mag.
- A model already displaced by $0.020$ mag reaches its nearest boundary after only $0.013234$ mag.
- The $0.05\rightarrow0.12$ scenario change produces $0.0665$ mag, or approximately $2.04$ times the quoted photometric uncertainty.

The last two residuals/scenarios are calculations, not additional observations.

For UV band $i$, the corresponding answer is $\epsilon_i/k_{i,\rm eff}$ for a centred model. A numerical RX J1856 UV verdict cannot honestly be supplied until its measured per-band uncertainties are local. Kaplan’s predicted-magnitude errors cannot fill that gap.

The assertion “small extinction cannot change R5” is therefore wrong. The weaker assertion “small extinction cannot explain the optical excess” survives. Kaplan’s RX J1856 row gives an excess `$7.1\pm0.3$` at $2500$ Å at [refs/kaplan_2011_rxj1856_spectrum.tex:1035](/home/tobias/Projects/neutronstars/refs/kaplan_2011_rxj1856_spectrum.tex:1035). A factor $7.1$ is a computed $2.128146$ mag discrepancy, far larger than the attenuation changes calculated above.

Nor does Ho’s fit automatically establish R5. Ho states that F606W is `15\% above our model spectra` at [refs/ho_2007_rxj1856_thin_atm.tex:710](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:710), then invokes model and normalization uncertainties in addition to photometric error. Under a fixed-spectrum, $k_{606}=0.95$ illustration, reducing $A_V$ from $0.12$ to $0.05$ increases the prediction by $6.316\%$ but leaves the observation approximately $8.17\%$ above it. This is not a formal reanalysis of Ho; it demonstrates why a broad physical error budget cannot be substituted for R5’s “within published errors.”

### 2.5 What owning the dust calculation would actually mean

**DERIVED HERE, conditional on independent grains and extinction out of the measured beam.** If species $s$ has grain number density per size $n_s(a,\ell)$ and extinction cross-section $C_{{\rm ext},s}(a,\lambda)$, an infinitesimal path removes the fraction

$$
\frac{-dF_\lambda}{F_\lambda}
=d\ell\sum_s\int n_s(a,\ell)C_{{\rm ext},s}(a,\lambda)\,da.
$$

Integrating and converting optical depth to magnitudes gives

$$
A_\lambda=\frac{2.5}{\ln10}
\int d\ell\sum_s\int
n_s(a,\ell)\,C_{{\rm ext},s}(a,\lambda)\,da.
$$

Define the efficiency by $C_{\rm ext}=\pi a^2Q_{\rm ext}$. Then an owned route needs:

1. A locally sourced electromagnetic boundary-value calculation for $Q_{\rm ext}$.
2. Grain material optical properties over the entire instrument band.
3. Size distributions, material fractions, and column normalizations.
4. An explicit treatment of the spherical-grain assumption and any light scattered back into the measurement aperture.
5. Constraints appropriate to the RX J1856 sightline.

**UNSOURCED locally:** the Mie coefficient equations and their stable numerical evaluation. Bohren–Huffman must be fetched before implementing them. Weingartner–Draine and Draine are candidate sources for the material and population model, not presently local authorities for equations in this argument.

The ownership distinction is fundamental. Solving Maxwell’s equations for a prescribed grain population does not determine that population. Importing a fitted grain-size distribution and fitted optical constants can merely replace one empirical dependency with several.

A measured extinction-curve dataset could instead be an allowed observational input, provided the underlying measurements, uncertainties and sightline applicability are documented. Renaming the CCM polynomial “measured data” would not accomplish that.

For a change $\delta k_i$ at fixed $A_V$,

$$
\delta m_i\simeq A_V\delta k_i.
$$

At $A_V=0.12$, moving a centred F606W prediction by its quoted photometric error requires $|\delta k_{606}|\simeq0.2714$; at $A_V=0.05$, it requires $0.6514$. These are useful comparison thresholds, **not bounds on the difference between CCM and an uncomputed grain model**.

An owned grain calculation might change a marginal verdict. It cannot be claimed to do so—or not to do so—before constraining its inputs and evaluating its band integrals.

## 3. X-ray absorption: the elemental sum and its missing leaves

### 3.1 Derivation and the difference between a neutral-atom sum and Wilms

For an independent mixture, define $\sigma_Z(E)$ as the removal cross-section per atom and $A_Z=n_Z/n_{\rm H}$. An infinitesimal hydrogen column gives

$$
dI_E=-I_E\sum_Z A_Z\sigma_Z(E)\,dN_{\rm H}.
$$

If the composition is constant along the path, integration yields, **DERIVED HERE**,

$$
\sigma_{\rm atom}(E)=\sum_Z A_Z\sigma_Z(E),
\qquad
T(E)=\frac{I_E}{I_{E,0}}
=\exp[-N_{\rm H}\sigma_{\rm atom}(E)].
$$

This is the requested sum under explicit neutral, atomic, undepleted assumptions. It is not the complete Wilms model.

Wilms prints

$$
\sigma_{\rm ISM}
=\sigma_{\rm gas}+\sigma_{\rm molecules}+\sigma_{\rm grains},
$$

at [refs/wilms_allen_mccray_2000_ism.tex:121](/home/tobias/Projects/neutronstars/refs/wilms_allen_mccray_2000_ism.tex:121), anchor `\sigism = \sigma_{\rm gas}`, and

$$
\sigma_{\rm gas}
=\sum_{Z,i}A_Za_{Z,i}(1-\beta_{Z,i})\sigma_{\rm bf}(Z,i)
$$

at [refs/wilms_allen_mccray_2000_ism.tex:142](/home/tobias/Projects/neutronstars/refs/wilms_allen_mccray_2000_ism.tex:142), anchor `\sum_{Z,i} A_Z`. Its exponential attenuation is printed at line 129, anchor `I_{\rm obs}(E) =`.

Thus “Wilms abundances plus Verner cross-sections” needs a named phase and ionization prescription before it defines an ISM model.

### 3.2 Abundances

Wilms’s table uses $12+\log_{10}A_Z$, not the abundance itself. The heading is `12+\log A_Z` at [refs/wilms_allen_mccray_2000_ism.tex:1325](/home/tobias/Projects/neutronstars/refs/wilms_allen_mccray_2000_ism.tex:1325). Exponentiating its ISM column gives:

| Element | Printed logarithmic abundance | Computed abundance per H |
|---|---:|---:|
| He | 10.99 | $9.7723722\times10^{-2}$ |
| C | 8.38 | $2.3988329\times10^{-4}$ |
| N | 7.88 | $7.5857758\times10^{-5}$ |
| O | 8.69 | $4.8977882\times10^{-4}$ |
| Ne | 7.94 | $8.7096359\times10^{-5}$ |
| Fe | 7.43 | $2.6915348\times10^{-5}$ |

Sources: [refs/wilms_allen_mccray_2000_ism.tex:1336](/home/tobias/Projects/neutronstars/refs/wilms_allen_mccray_2000_ism.tex:1336)–1340 and 1353; short row anchors include `10.99`, `8.38`, `8.69`, and `7.43`.

These are adopted environmental estimates. Wilms explicitly describes uncertainties as:

> `0.1\,dex or higher`

[refs/wilms_allen_mccray_2000_ism.tex:382](/home/tobias/Projects/neutronstars/refs/wilms_allen_mccray_2000_ism.tex:382). A $0.1$ dex increase is a computed factor $1.258925$. This is not automatically an uncertainty of that size in the *total* opacity: each element must be weighted by its contribution.

### 3.3 The Verner fit, including its validity boundary

Verner prints

$$
\begin{aligned}
\sigma(E)&=\sigma_0F(y)\ {\rm Mb},\\
x&=\frac{E}{E_0}-y_0,\qquad y=\sqrt{x^2+y_1^2},\\
F(y)&=\left[(x-1)^2+y_w^2\right]
y^{P/2-5.5}
\left(1+\sqrt{y/y_a}\right)^{-P}.
\end{aligned}
$$

Source: [refs/verner_1996_photoionization.tex:179](/home/tobias/Projects/neutronstars/refs/verner_1996_photoionization.tex:179)–191, with greppable anchors `\sigma(E)=\sigma_0` and `F(y) = \left[(x-1)^2`. The paper specifies energy in eV and `$10^{-18}$~cm$^2$` per Mb.

As an algebraic check, at sufficiently high energy the three factors scale as $E^2$, $E^{P/2-5.5}$ and $E^{-P/2}$. Their product is therefore $E^{-3.5}$, **DERIVED HERE**, matching the printed asymptote at [refs/verner_1996_photoionization.tex:213](/home/tobias/Projects/neutronstars/refs/verner_1996_photoionization.tex:213), anchor `E^{-3.5}`.

The visually verified neutral H and He rows are:

| Parameter | H I | He I |
|---|---:|---:|
| $E_{\rm th}$, eV | 13.60 | 24.59 |
| $E_{\max}$, eV | 50000 | 50000 |
| $E_0$, eV | 0.4298 | 13.61 |
| $\sigma_0$, Mb | 54750 | 949.2 |
| $y_a$ | 32.88 | 1.469 |
| $P$ | 2.963 | 3.188 |
| $y_w$ | 0 | 2.039 |
| $y_0$ | 0 | 0.4434 |
| $y_1$ | 0 | 2.136 |

Sources: [refs/verner_1996_photoionization.txt:402](/home/tobias/Projects/neutronstars/refs/verner_1996_photoionization.txt:402)–403, anchors `5.475+4   3.288+1` and `9.492+2   1.469+0`. The scan was essential because the extracted text loses some minus signs in exponents.

The critical qualification is:

> `the present fits should be supplemented by the fits to the inner shell`

[refs/verner_1996_photoionization.tex:244](/home/tobias/Projects/neutronstars/refs/verner_1996_photoionization.tex:244).

For example, the local C I row has $E_{\max}=291.0$ eV, [refs/verner_1996_photoionization.txt:417](/home/tobias/Projects/neutronstars/refs/verner_1996_photoionization.txt:417), anchor `2.910+2`; O I has $E_{\max}=538.0$ eV at line 430, anchor `5.380+2`. Evaluating these outer-shell fits alone at the requested higher energies and calling the result the total atomic opacity would omit absorption channels.

Furthermore, Wilms’s actual prescription is not uniformly Verner (1996):

- H uses Band et al.; [refs/wilms_allen_mccray_2000_ism.tex:216](/home/tobias/Projects/neutronstars/refs/wilms_allen_mccray_2000_ism.tex:216): `use the fitting formula of \citet{band:90a}`.
- He uses Yan et al.; line 224: `evaluation of \citet{yan:98a}`.
- Other elements use Verner’s earlier subshell calculations; line 243: `by \citet{verner:95a}`.

The audit’s shorthand therefore conflates distinct implementations.

### 3.4 Numerical transmissions supported by the local inputs

First, evaluate the complete neutral H and He contributions using the local Verner parameters and Wilms’s He abundance:

$$
\sigma_{\rm HHe}
=\sigma_{\rm H}+A_{\rm He}\sigma_{\rm He}.
$$

For the requested $N_{\rm H}=10^{20}\ {\rm cm^{-2}}$:

| $E$, keV | $\sigma_{\rm H}$, cm² | $A_{\rm He}\sigma_{\rm He}$, cm²/H | $\sigma_{\rm HHe}$, cm²/H | $T_{\rm HHe}$ |
|---|---:|---:|---:|---:|
| 0.2 | $2.2047104\times10^{-21}$ | $5.3934979\times10^{-21}$ | $7.5982083\times10^{-21}$ | 0.467750 |
| 0.5 | $1.1351310\times10^{-22}$ | $3.1161063\times10^{-22}$ | $4.2512373\times10^{-22}$ | 0.958379 |
| 1.0 | $1.1408673\times10^{-23}$ | $3.2738184\times10^{-23}$ | $4.4146858\times10^{-23}$ | 0.995595 |

**These are not full ISM transmissions.** Within this stipulated neutral-atom model, adding positive metal opacity lowers the transmission.

An independent 40-digit decimal evaluation of the H fit at $200$ eV gives $2.2047104085638850\times10^{-21}\ {\rm cm^2}$; binary floating-point evaluation gives $2.2047104085638926\times10^{-21}$. Relative arithmetic disagreement is $3.43\times10^{-15}$, below the chosen $10^{-12}$ numerical tolerance. This checks evaluation, not atomic physics.

For a complete-ISM comparison available *locally*, I extracted the dotted MRN curve from Wilms Fig. 1. Its source is [refs/wilms_allen_mccray_2000_ism.tex:486](/home/tobias/Projects/neutronstars/refs/wilms_allen_mccray_2000_ism.tex:486), anchor `cosab.eps`; the caption specifies:

> `absorptivity including grains with a MRN distribution`

at line 489.

The extraction used the vector coordinates in `refs/tex/wilms_allen_mccray_2000_ism/cosab.eps`, not a screenshot estimate. For reproducibility, the main-panel mapping is

$$
X=10647+6599\log_{10}E_{\rm keV},
\qquad
Q=\frac{400(Y-1015)}{8261-1015},
\qquad
\sigma=Q\,10^{-24}E_{\rm keV}^{-3}.
$$

These equations are **DERIVED HERE from the plotted tick coordinates**. Interpolation was restricted to adjacent curve vertices at the requested energies, away from edges. The extracted $Q$ values are $68.15296$, $78.77450$ and $173.88352$.

This yields:

| $E$, keV | Wilms Fig. 1 MRN $\sigma$, cm²/H | Computed $\tau=N_{\rm H}\sigma$ | Computed $T=e^{-\tau}$ |
|---|---:|---:|---:|
| 0.2 | $\approx8.5191\times10^{-21}$ | 0.85191 | $\approx0.4266$ |
| 0.5 | $\approx6.3020\times10^{-22}$ | 0.06302 | $\approx0.9389$ |
| 1.0 | $\approx1.7388\times10^{-22}$ | 0.01739 | $\approx0.9828$ |

**Provenance: published-figure verification target, not a recomputed elemental sum.** I propose a conservative $2\%$ cross-section tolerance when using these extracted values as coarse verification fixtures. This is an engineering tolerance for figure extraction and matching assumptions, not a confidence interval for the physical model. Its first-order transmission equivalents are approximately $0.0073$, $0.0012$ and $0.00034$.

The complete Wilms/Verner sum requested in the brief cannot be honestly claimed until the missing subshell and species prescriptions are local. The figure-derived table supplies an explicitly limited comparison rather than disguising that omission.

### 3.5 Morrison–McCammon comparison

Morrison–McCammon Table 2 gives

$$
\sigma_{\rm MM83}(E)
=(c_0+c_1E+c_2E^2)E^{-3}\,10^{-24}\ {\rm cm^2/H},
$$

with $E$ in keV. Source: [refs/morrison_mccammon_1983_ism.txt:136](/home/tobias/Projects/neutronstars/refs/morrison_mccammon_1983_ism.txt:136), anchor `Cross section per hydrogen atom =`; the formula and coefficients were visually checked on PDF page 3.

The applicable rows are:

| Energy interval, keV | $c_0$ | $c_1$ | $c_2$ | Local line |
|---|---:|---:|---:|---:|
| 0.100–0.284 | 34.6 | 267.9 | −476.1 | 122 |
| 0.400–0.532 | 71.4 | 66.8 | −51.4 | 124 |
| 0.867–1.303 | 120.6 | 169.3 | −47.7 | 127 |

For example, the greppable first-row anchor is `34.6    267.9     -476.1`.

The Python results are:

| $E$, keV | $\sigma_{\rm MM83}$, cm²/H | $T_{\rm MM83}$ | Wilms figure $T$ | Absolute transmission difference |
|---|---:|---:|---:|---:|
| 0.2 | $8.6420\times10^{-21}$ | 0.421389 | $\approx0.426599$ | $\approx0.005210$ |
| 0.5 | $7.3560\times10^{-22}$ | 0.929080 | $\approx0.938925$ | $\approx0.009845$ |
| 1.0 | $2.4220\times10^{-22}$ | 0.976071 | $\approx0.982762$ | $\approx0.006691$ |

Decimal arithmetic reproduces the three cross-sections from the printed coefficients with zero discrepancy at the displayed precision.

The MM83 cross-sections exceed the Wilms figure values by approximately $1.44\%$, $16.73\%$ and $39.29\%$. Small absolute transmission differences at the higher energies do not imply small differences in fitted column.

For a fixed inferred optical depth,

$$
N_{{\rm H},2}
=N_{{\rm H},1}\frac{\sigma_1}{\sigma_2},
$$

**DERIVED HERE** from $\tau=N_{\rm H}\sigma$. The energy dependence of the above ratios means there is no single exact conversion between MM83 and Wilms $N_{\rm H}$.

More generally,

$$
\delta\ln T=-\sigma\,\delta N_{\rm H}-N_{\rm H}\delta\sigma.
$$

A deliberately imposed $10\%$ increase in the Wilms figure cross-section changes transmission by approximately $-8.17\%$, $-0.628\%$ and $-0.174\%$ at the three energies. This is a sensitivity experiment, not an assigned physical uncertainty.

### 3.6 The owned atomic route is deeper than the audit suggests

The audit’s small wrapper estimate must not be interpreted as the effort needed to own photoionization physics.

An actual route requires bound-state structure, continuum states, transition matrix elements, shell thresholds and occupations, followed by sums over channels and a controlled treatment of resonances. Verner itself distinguishes calculations with an `insufficient number of` target states in the close-coupling expansion at [refs/verner_1996_photoionization.tex:267](/home/tobias/Projects/neutronstars/refs/verner_1996_photoionization.tex:267). That is a convergence problem in the atomic calculation, not in the final abundance sum.

The breadth-first route should be:

1. Own analytic hydrogen photoionization.
2. Own or independently validate helium, which is important in the lowest-energy calculation above.
3. Add C, N, O and Ne with all relevant subshells.
4. Add Fe and the remaining species, proving any proposed omission by its contribution to the folded spectrum.
5. Add molecular, depletion and grain treatments under explicit assumptions.

**UNSOURCED locally:** the implementable atomic Hamiltonian, continuum normalization, dipole cross-section prefactor, and R-matrix or Hartree–Fock–Slater algorithm. These need dedicated sources before production code.

Fetching Opacity Project tables would supply verification targets. Using them at runtime would remain a table dependency under Rule 7. Likewise, reproducing an approximate Hartree–Fock–Slater model would establish ownership of that approximation, not erase its physical uncertainty.

## 4. The $N_{\rm H}$–$A_V$ relation is a consistency diagnostic

Predehl–Schmitt Fig. 3 gives the regression

$$
N_{\rm H}\simeq
(1.79\pm0.03)\times10^{21}A_V
\ {\rm cm^{-2}}.
$$

Source: [refs/predehl_schmitt_1995_nh_av.txt:1048](/home/tobias/Projects/neutronstars/refs/predehl_schmitt_1995_nh_av.txt:1048)–1050, anchored by:

> `Visual extinction vs. equivalent hydrogen column density`

The numerical relation was checked against PDF page 9 because the OCR damages both the uncertainty sign and exponent.

Two source facts prevent treating this as a precise conversion for RX J1856.

First, the X-ray columns in that regression were fitted using MM83:

> `by Morrison & McCammon (1983)`

[refs/predehl_schmitt_1995_nh_av.txt:249](/home/tobias/Projects/neutronstars/refs/predehl_schmitt_1995_nh_av.txt:249).

Second, the authors discuss a zero-point offset corresponding to $2.3\times10^{20}\ {\rm cm^{-2}}$, visually verified on the same page. The local anchor is `our zero point off-` at [refs/predehl_schmitt_1995_nh_av.txt:1073](/home/tobias/Projects/neutronstars/refs/predehl_schmitt_1995_nh_av.txt:1073). This corresponds to a computed $0.12849$ mag when divided by the regression slope—larger than either extinction scenario here. It is not a measured uncertainty to attach mechanically to RX J1856; it demonstrates the danger of treating a population regression as an exact low-column identity.

### 4.1 Implied extinctions, kept separate by fit

Applying the nominal slope only:

| Source/model | Published $N_{\rm H}/10^{20}\ {\rm cm^{-2}}$ | Implied $A_V$, mag | Propagated quoted $N_{\rm H}$ error only |
|---|---:|---:|---:|
| Burwitz, Chandra LETGS | $0.95\pm0.03$ | 0.053073 | 0.001676 |
| Burwitz, XMM RGS | $0.87\pm0.08$ | 0.048603 | 0.004469 |
| Ho, atmosphere with $B=3\times10^{12}$ G | $1.18(2)$ | 0.065922 | 0.001117 |
| Ho, atmosphere with $B=4\times10^{12}$ G | $1.30(2)$ | 0.072626 | 0.001117 |
| Ho, blackbody | $0.91(1)$ | 0.050838 | 0.000559 |

Burwitz sources: [refs/burwitz_2003_rxj1856_xray.tex:248](/home/tobias/Projects/neutronstars/refs/burwitz_2003_rxj1856_xray.tex:248), anchor `0.95&$\pm$0.03`; line 251, anchor `0.87&$\pm$0.08`. Ho source: [refs/ho_2007_rxj1856_thin_atm.tex:634](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:634), anchor `1.18(2) & 1.30(2) & 0.91(1)`.

These small propagated errors exclude regression scatter, sightline variation, absorption-model differences and surface-model covariance. They are not legitimate complete errors on $A_V$.

Conversely,

$$
A_V=0.05\Rightarrow N_{\rm H}=8.95\times10^{19}\ {\rm cm^{-2}},
$$

$$
A_V=0.12\Rightarrow N_{\rm H}=2.148\times10^{20}\ {\rm cm^{-2}}.
$$

Ho explicitly reports trying the relation and free extinction, finding the fits `too unconstrained`, and then fixing $A_V$; [refs/ho_2007_rxj1856_thin_atm.tex:690](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:690)–695. Kaplan likewise says reddening is `highly degenerate with the power-law index` at [refs/kaplan_2011_rxj1856_spectrum.tex:597](/home/tobias/Projects/neutronstars/refs/kaplan_2011_rxj1856_spectrum.tex:597).

The proposed inference should therefore use independently parameterized attenuation:

$$
\mathcal L=
\mathcal L_{\rm HST}(D_{\rm HST}\mid\theta,A_V,R_V,\mathcal D)
\,
\mathcal L_{\rm X}(D_{\rm X}\mid\theta,N_{\rm H},\mathcal A),
$$

where $\theta$ denotes surface parameters, $\mathcal D$ the extinction model, and $\mathcal A$ the absorption model. This factorization assumes conditionally independent datasets and is **DERIVED HERE as the proposed statistical model**.

Only afterward should one inspect

$$
\eta=
\frac{N_{\rm H}}
{1.79\times10^{21}A_V}.
$$

Here “consistency prior” means an external population expectation used to assess plausibility. A Bayesian prior multiplied into the likelihood would itself be a fit input; that is deliberately excluded from the baseline proposed here. No additional deterministic $A_V(N_{\rm H})$ constraint should be imposed.

## 5. M2 stub policy

The defensible temporary policy is:

| Component | M2 treatment | Required uncertainty/provenance label |
|---|---|---|
| CCM89 | Implement the printed function and integrate through measured throughputs | `fit`, `stub`; mean extinction law, sightline applicability unquantified |
| $A_V$ | Keep $0.05$ and $0.12$ as separately named scenarios; profile or fit independently when data permit | Neither scenario is an independent direct measurement of dust column |
| $R_V$ | Use $3.1$ for the requested reference calculation | Fixed assumption; sensitivity requires an explicitly declared range |
| Ho reproduction | Include its extinction-law variant | Ho used `corrections from O'Donnell~(1994)`, [refs/ho_2007_rxj1856_thin_atm.tex:698](/home/tobias/Projects/neutronstars/refs/ho_2007_rxj1856_thin_atm.tex:698); plain CCM is not exact reproduction |
| X-ray absorption | Use a complete, versioned prescription after fetching missing inputs | Cross-sections, abundances, ionization, molecules and depletion recorded separately |
| MM83 | Alternative-model and historical verification target | Model spread is not a calibrated confidence interval |
| Predehl–Schmitt | Post-fit consistency diagnostic | No deterministic conversion and no narrow regression-slope-only prior |
| R5 | Report residuals against observational errors | Report physical-model uncertainty alongside the test; do not enlarge its tolerance silently |

The branch remains yellow even when its numerical evaluator is owned: evaluating a published fit correctly does not change the fit’s provenance.

An incomplete metal sum, missing passbands, or missing measured photometry cannot be made defensible by attaching a generic percentage. Those are missing computations or inputs, not characterized uncertainties.

## 6. Concrete implementation plan

All production work should follow red-green TDD. The scratch calculations in this argument are not proposed runtime code.

| Function | Inputs | Outputs and failure conditions | Verification |
|---|---|---|---|
| `ccm89_k(lambda_angstrom, Rv)` | Wavelength, positive $R_V$ | Dimensionless $k$; reject non-finite values and unsupported wavelength ranges | $k(1.82\ \mu{\rm m}^{-1})=1$; four monochromatic values above; numerical tolerance $10^{-12}$ for normalization |
| `attenuate_optical(F_lambda, lambda, Av, curve)` | Spectrum with explicit units, nonnegative $A_V$ | Attenuated spectrum and provenance | $A_V=0$ identity; positive extinction reduces flux; inverse recovers input |
| `synthetic_band_signal(F_lambda, response)` | Spectrum, checksummed detector response | Signal, magnitude, coverage diagnostics | Recover the actual published measured photometry; compare extinction coefficients with Kaplan’s rounded band targets under matching assumptions |
| `band_extinction_and_derivative(...)` | Same spectrum/response, $A_V,R_V$ | $A_i$ and $\partial A_i/\partial A_V$ | Finite-difference derivative against the weighted-mean derivation; small-$A_V$ limit |
| `verner_outer_sigma(E_eV, species)` | Energy, explicitly identified parameter row | Cross-section and shell coverage metadata | H at $200$ eV: $2.20471040856\times10^{-21}\ {\rm cm^2}$; fail on requests for a total atomic cross-section with missing inner shells |
| `atomic_sigma(E, species, atomic_model)` | Complete channel definitions | Total atomic cross-section and channel contributions | Verner/experimental targets; threshold coverage and high-energy asymptotes |
| `ism_sigma(E, composition, phase_model)` | Atomic cross-sections, abundances, ion fractions, molecular/grain assumptions | Cross-section per total H nucleus | Full Wilms figure targets above with proposed $2\%$ coarse comparison tolerance |
| `mm83_sigma(E_keV)` | Energy within printed intervals | Historical fitted cross-section | $8.6420\times10^{-21}$, $7.3560\times10^{-22}$, $2.4220\times10^{-22}\ {\rm cm^2/H}$ |
| `xray_transmission(NH, sigma)` | Nonnegative column and opacity | $e^{-N_{\rm H}\sigma}$ | MM83 transmissions above; zero-column identity; monotonicity |
| `fold_xray(photon_spectrum, ARF, RMF, exposure)` | Photon spectrum and measured response files | Expected channel counts | Recover published fits only with matched data, responses and absorption assumptions |
| `r5_report(model, photometry, attenuation_scenarios)` | Predictions, actual measured errors, named scenarios | Per-band residuals, pass/fail, extinction margins | F606W centred-model sensitivity approximately $0.0343$ mag in $A_V$ |
| `nh_av_consistency(NH, Av, absorption_model)` | Independently inferred values | $\eta$ and applicability warnings | Numerical conversions in §4; never returns a forced fit parameter |

Mutation testing should target the extinction sign, wavelength units, $b/R_V$, the eV/keV conversion, the Mb conversion, the factor $10^{20}$ in columns, and omission of shell contributions. Fuzzing should exercise validity boundaries, malformed parameter rows, non-finite inputs and energy points on either side of absorption edges.

Do not require global energy monotonicity of the cross-section: absorption edges invalidate that test. Require monotonic transmission with column at fixed energy and fixed composition.

The order of work should be: correct the photometric source identification; fetch the measurements and responses; implement band integration; complete the absorption stub; fit and report the attenuation ensemble. Start the deep dust or atomic route after this establishes which uncertainty actually controls R5.

## 7. Missing sources to fetch

Every item below is **UNVERIFIED FOR LOCAL IMPLEMENTATION** until its contents are stored under `refs/`, checked, and anchored. Online links identify acquisition targets only; no additional numerical physics above depends on them.

| Priority | Source | Purpose |
|---|---|---|
| Immediate | **van Kerkwijk, M. H., & Kulkarni, S. R. (2001), A&A 378, 986–995**, *Optical spectroscopy and photometry of the neutron star RX J1856.5−3754*, **astro-ph/0106265** | Actual WFPC2 photometry, errors and calibration. This corrects GAPS §2.1; identity confirmed by the [authors’ repository](https://authors.library.caltech.edu/records/td93m-msp87). |
| Immediate | **STScI instrument teams, versioned WFPC2 and ACS/HRC calibration datasets and documentation**; no single year/arXiv identifier | Exact throughput, magnitude system, detector response and true pivot wavelengths; checksum each file |
| Immediate | **Verner, D. A., & Yakovlev, D. G. (1995), A&AS 109, 125**; arXiv ID not established | Missing inner-shell fits and complete parameter files; locally identified as Paper II at `refs/verner_1996_photoionization.tex:390` |
| Immediate | **Yan, M., Sadeghpour, H. R., & Dalgarno, A. (1998), ApJ 496, 1044**; arXiv ID not established | Wilms H/He/molecular prescription; bibliography at `refs/wilms_allen_mccray_2000_ism.tex:1154` |
| Immediate | **Band, I. M., Trzhaskovskaya, M. B., Verner, D. A., & Yakovlev, D. G. (1990), A&A 237, 267**; arXiv ID not established | Exact historical Wilms H-fit reproduction; bibliography at `refs/wilms_allen_mccray_2000_ism.tex:914` |
| Immediate | **O’Donnell, J. E. (1994), ApJ 422, 158**; arXiv ID not established | Ho’s optical extinction-law correction; locally listed at `refs/ho_2007_rxj1856_thin_atm.tex:1152` |
| Immediate | **Chandra archive and calibration teams, observation-specific LETG products and responses** | Joint inference; preserve extraction and calibration versions |
| Dust ownership | **Weingartner, J. C., & Draine, B. T. (2001), ApJ 548, 296**, **astro-ph/0008146** | Grain populations, size distributions and the distinction between their fitted inputs and recomputed outputs; [author preprint](https://arxiv.org/abs/astro-ph/0008146) |
| Dust ownership | **Draine, B. T. (2003), ARA&A 41, 241–289**, **astro-ph/0304489** | Dust materials, observational constraints and limitations; [author preprint](https://arxiv.org/abs/astro-ph/0304489) |
| Dust ownership | **Bohren, C. F., & Huffman, D. R. (1983), Wiley**, *Absorption and Scattering of Light by Small Particles*; no arXiv ID | Mie derivation, normalization and numerical method. Verify edition: the [publisher’s online record](https://onlinelibrary.wiley.com/doi/book/10.1002/9783527618156) describes a later edition. |
| Atomic ownership | **Seaton, M. J., Yan, Y., Mihalas, D., & Pradhan, A. K. (1994), MNRAS 266, 805–828**, *Opacities for stellar envelopes*; arXiv ID not established | Opacity Project overview and routes to species-specific calculations; [publisher record](https://academic.oup.com/mnras/article/266/4/805/982644) |
| Atomic ownership | **Bethe, H. A., & Salpeter, E. E. (1957), Springer**, *Quantum Mechanics of One- and Two-Electron Atoms*; no arXiv ID | Owned hydrogen calculation and normalization; locally referenced at `refs/wilms_allen_mccray_2000_ism.tex:917` |
| Atomic validation | **Samson, J. A. R., He, Z. X., Yin, L., & Haddad, G. N. (1994), J. Phys. B 27, 887**; arXiv ID not established | Experimental He cross-sections and uncertainties; local bibliography at `refs/verner_1996_photoionization.tex:375` |
| Environmental inputs | **Cardelli et al. (1996), ApJ 467, 334; Meyer et al. (1997), ApJ 490, L103; Meyer et al. (1998), ApJ 493, 222; Snow & Witt (1996), ApJ 468, L65**; arXiv IDs not established | Underlying abundance evidence cited by Wilms, plus separately selected laboratory grain optical constants |

Seaton’s overview alone is not an implementation specification. The next atomic source selection must identify the actual continuum method and species calculations. Likewise, a dust review does not supply a measured grain distribution for RX J1856.

## 8. Estimated effect on the final render and R5 if this path stays stubbed

**Earth-distance photometry.** Keeping a complete CCM stub while varying $A_V$ between the two published scenarios changes the approximate WFPC2 magnitudes by $0.0665$–$0.1855$ mag. That is enough to change a marginal R5 result. Omitting extinction entirely is not defensible against the quoted F606W precision. Using a complete labelled stub is defensible for M2, provided the result remains conditional on that attenuation model.

**Blackbody rejection.** The extinction ambiguity does not remove the factor-of-several optical excess. It therefore does not overturn the qualitative reason for proceeding beyond the X-ray-normalized single blackbody.

**Inferred surface normalization.** At fixed spectral shape, changing the inferred intrinsic F606W signal between the two extinction scenarios changes it by approximately $6.32\%$. If that change is absorbed entirely into an emitting-area normalization, the inferred linear radius changes by approximately $3.11\%$, since area scales as radius squared. The corresponding approximate F170W–F606W colour change is $0.1190$ mag, so one normalization adjustment cannot absorb the whole effect.

**X-ray-constrained surface parameters.** At fixed intrinsic spectrum, changing the column from Burwitz’s $0.95$ to Ho’s $1.30$ in units of $10^{20}\ {\rm cm^{-2}}$ changes the Wilms-benchmark transmitted flux by factors approximately $0.7422$, $0.9782$ and $0.9939$ at $0.2$, $0.5$ and $1$ keV. This strongly energy-dependent change must be propagated through a refit. It cannot be converted into a defensible temperature, radius or rendered-colour error without the spectrum and response matrices.

**Close flyby.** PRD excludes interstellar extinction at the camera. At fixed surface parameters, these foreground-screen stubs therefore have zero direct effect on the close-camera radiance calculation. Their effect is indirect: they alter the parameters and surface models inferred from Earth observations. The optical calculation above establishes a several-percent normalization sensitivity, while the X-ray propagation remains unquantified until refitting.

**R5 verdict.** A stubbed path can support a statement such as “this surface model meets the measured photometric tolerances under this explicitly named attenuation prescription.” It cannot support an unconditional R5 claim when the extinction scenario changes the result, and it cannot turn the provenance tree green. The decisive next computation is the measured-passband, jointly refitted attenuation ensemble—not an uncalibrated first-principles dust model.