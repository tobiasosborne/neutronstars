# Failure Log

Record of what went wrong, why, and what we learned.

## F1: Wiley book chapter downloads blocked by Cloudflare Turnstile
**Date:** 2026-03-24
**What:** wget, curl, and headless Playwright all fail to download chapter PDFs from Wiley Online Library, despite TIB VPN access confirmed (book landing pages load fine).
**Why:** Cloudflare Turnstile CAPTCHA on individual chapter/PDF URLs. Even with webdriver spoofing, the Turnstile challenge blocks automated download.
**Workaround:** Book landing pages load in stealth headless mode. PDFs must be downloaded manually via browser. 4 books affected: Rybicki & Lightman, Shapiro & Teukolsky, Mészáros, Mihalas.
**Impact:** The key equations from these books are also available in the arXiv papers we successfully downloaded (Potekhin papers contain the opacity formulae; Haakonsen describes the Feautrier method; Potekhin 2013 has the TOV/EOS). Books needed primarily for additional verification and the original derivations.

## F2: Magnetic B=10^12 G atmosphere converged to excessive flux
**Date:** 2026-05-13
**What:** After splitting magnetic true absorption from scattering, `solve_magnetic_atmosphere(1e6, 2e14, 1e12, π/4, gaunt; K=50, M=8, N=100)` converges by `max|ΔT/T| < 1e-4` but reports `F/σT⁴ = 16.1600`, not ≈1.
**Why:** The magnetic solver used a non-magnetic grey initial temperature profile after extending the column to magnetic optical depth, which imposed a too-hot LTE bottom boundary. The convergence criterion also ignored the surface flux error, so it could accept a locally balanced but wrongly normalized atmosphere.
**Fix:** Added a magnetic directional grey starting profile, a damped target-flux correction, a flux-aware convergence criterion, and a thermalization-depth guard before applying the LTE bottom boundary. The same `B=10^12 G` case now converges with `F/σT⁴ ≈ 1.000`.
**Remaining impact:** Stronger fields, e.g. RX J1856-like `B≈3×10^12 G`, are now flux-normalized in test runs but still show local temperature-correction oscillations. Treat those spectra as not fully converged until the magnetic correction scheme is improved further.

## F3: Suleimanov Fig. 2 free-free opacity mismatch
**Date:** 2026-05-14
**What:** Digitized Suleimanov, Potekhin & Werner (2009) Fig. 2 and overplotted local magnetic opacity components for the paper parameters `T=7e6 K`, `ρ=30 g cm^-3`, `B=1e14 G`, and `θ_B={5,45,90} deg`. Electron scattering branches are often close at the figure-digitization level, but free-free/true absorption is systematically below the digitized magenta branches by about 1-2.4 dex outside resonance-marker bands.
**Why:** Two implementation errors in the Potekhin & Chabrier (2003) magnetic Coulomb logarithm suppressed the magnetic free-free opacity. Eq. (44e) was transcribed with division by `sqrt(1/4+y/β_e)` instead of multiplication, and Eq. (44a) was truncated at `y=20` although the published integral is over `0..∞`. At the Fig. 2 low-energy longitudinal point this gave `Λ_0≈0.026` instead of `Λ_0≈7.415`.
**Fix:** Corrected Eq. (44e), integrated Eq. (44a) over the infinite domain, and made the Eq. (52) damping split source-faithful by removing proton-proton damping from the electron damping channel. The Fig. 2 metric gate now passes: the clean free-free branches at `θ_B=5°` and `90°` have median absolute residuals below `0.15 dex`.
**Current evidence:** `verification/figures/suleimanov_2009_fig2/opacity_comparison.png` and `verification/data/suleimanov_2009_fig2/opacity_comparison_metrics.json`.
**Remaining impact:** The comparison excludes the annotated proton-cyclotron and vacuum-resonance bands. Vacuum-polarized local mode vectors are now included, but partial/adiabatic mode conversion through the resonance is still not implemented, so spectra across resonance regions are not yet fully physical.

## F4: Suleimanov Fig. 2 oblique-angle branch mismatch remains
**Date:** 2026-05-14
**What:** A cleaner branch-level validation view was added after the basic nearest-branch Fig. 2 metric passed. With resonance-marker bands and tall annotation segments masked, `theta_B=90°` is good and `theta_B=5°` is mixed, but `theta_B=45°` remains badly mismatched. The clean branch metrics show median absolute residuals of about `2.7-3.0 dex` for the upper magenta/blue branches at `theta_B=45°`.
**Why:** Not resolved. This is not explained solely by the vacuum resonance, because the clean validation masks the vacuum-resonance band. Likely suspects are oblique-angle mode/branch convention, B-frame cyclic weight sign/order at intermediate angles, missing scattering redistribution, or a real formula difference between van Adelsberg & Lai (2006) and the local P&C component implementation.
**Current evidence:** `verification/figures/suleimanov_2009_fig2/clean_branch_validation.png` and `verification/data/suleimanov_2009_fig2/clean_branch_metrics.json`.
**Next diagnostic:** For `T=7e6 K`, `ρ=30 g cm^-3`, `B=1e14 G`, `theta_B=45°`, compare fixed-energy intermediate quantities (`β`, `K_j`, `K_z,j`, cyclic weights, `σ_α`, and final mode opacities) directly against the published equations before changing the RT layer.
