| `L_{sc}` | scaled luminance entering the tone operator, L_sc = t_exp k_exp Y_pre (project convention); same symbol as the Reinhard scaled luminance. The synonym L_in is retired | 1 | `L_{sc}` (colour.display-exposure-convention) | colour.display-exposure-convention | reuse existing L_sc row; no new symbol |
| `Y_{pre}` | unexposed Y of the rendered spectrum (Y of pipeline.display-colour before exposure); bare L is retired | same as Y | `Y_{pre}` (colour.display-exposure-convention) | colour.display-exposure-convention | new |
| `t_{exp}` | exposure time factor of the stated camera model (logged per frame); t_exp > 0 | s | `t_{exp}` (colour.display-exposure-convention) | colour.display-exposure-convention | new |
| `k_{exp}` | exposure scale factor such that L_sc = t_exp k_exp Y_pre (logged per frame); k_exp > 0 | [Y_pre]^-1 s^-1 | `k_{exp}` (colour.display-exposure-convention) | colour.display-exposure-convention | new |
| `L_d` | existing row (colour.notation-reinhard-luminance); no change | 1 | `L_d` | colour.notation-reinhard-luminance | unchanged |
| `I_{psf}` | specific intensity of the ideal Airy PSF at off-axis angle theta_psf | same as I | `I_{psf}` (rt.airy-intensity-pattern) | rt.airy-intensity-pattern | new; replaces bare I(theta) |
| `I_{peak}` | central (theta_psf = 0) intensity of the Airy PSF; replaces I_0, which is emitted surface intensity | same as I | `I_{peak}` (rt.airy-intensity-pattern) | rt.airy-intensity-pattern | new; I_0 collision avoided |
| `\theta_{psf}` | small off-axis angle in the camera frame, theta_psf = r_fp/f_foc (paraxial); theta alone is the surface-normal ray angle | rad | `\theta_{psf}` (rt.airy-intensity-pattern) | rt.airy-intensity-pattern | new; theta collision avoided |
| `x_{psf}` | Airy argument pi D_ap theta_psf/lambda; bare x is CIE chromaticity | 1 | `x_{psf}`, `r'` (r' = pi D_ap r_fp/(f_foc lambda)) | rt.airy-intensity-pattern | new; x collision avoided |
| `r_{fp}` | radial distance in the focal plane; bare r is stellar areal radius | cm | `r` (source) | rt.airy-intensity-pattern | new; r collision avoided |
| `\mathrm{J}_1` | Bessel function of the first kind, order 1 (distinct from J_l emissivity and J = 1 - R) | 1 | `J_1` (rt.airy-intensity-pattern) | rt.airy-intensity-pattern | upright J spelling avoids the J_l source spelling |
| `D_{ap}` | circular aperture diameter of the camera model (logged parameter); bare D is observer distance | cm | `D` (source) | rt.airy-intensity-pattern | new; D collision avoided |
| `\lambda` | wavelength | cm | `\lambda` | colour.notation-cie-cmf | unchanged, reuse |
| `f_{foc}` | camera focal length (source symbol f) | cm | `f` | rt.airy-intensity-pattern | new |
| `\theta_{1,psf}` | angular position of the first Airy zero, (x_1/pi) lambda/D_ap ~ 1.22 lambda/D_ap | rad | `\theta_{1,psf}` | rt.airy-first-zero | new |
| `x_1` | first positive root of J_1, x_1 = 3.8317059702 | 1 | `x_1` (derivation from Eqs. 4-5) | rt.airy-first-zero | new |
| `K_\lambda`, `A_\lambda` | unit-integral Airy kernel K_lambda = A_lambda / integral of A_lambda over the support; A_lambda = [2 J_1(x_psf)/x_psf]^2 peak-normalised | sr^-1, 1 | `K_\lambda`, `A_\lambda` | rt.psf-aperture-convention | new |
