# Design Decisions

## D1: Start with Phase 2 tracer bullet using blackbody placeholder
**Date:** 2026-03-24
**Decision:** Implement end-to-end pipeline with modified blackbody atmosphere before self-consistent RT.
**Rationale:** Working pipeline with known approximations beats perfect opacity module with no pipeline. The blackbody placeholder (f_col correction from Suleimanov+ 2009) is a standard observational approximation.

## D2: CGS units throughout physics code
**Date:** 2026-03-24
**Decision:** All physics functions use CGS (erg, cm, g, s, Gauss).
**Rationale:** Standard in astrophysics literature. All reference papers use CGS. Avoids unit conversion errors when comparing against published values.

## D3: Schwarzschild first, Kerr deferred to Phase 4
**Date:** 2026-03-24
**Decision:** Ray tracing starts with Schwarzschild metric (no rotation).
**Rationale:** Most NS spin periods >> R/c ≈ 40 μs, so rotation effects are small. Schwarzschild has exact elliptic integral solutions (Pechenick+ 1983). Kerr adds significant complexity (Dexter & Agol 2009).

## D4: Fully ionised H atmosphere only in Phase 2-3
**Date:** 2026-03-24
**Decision:** Defer partial ionisation (Potekhin & Chabrier 2004) to Phase 4.
**Rationale:** Partial ionisation matters for T_eff < 5×10⁵ K and B > 10¹² G. The canonical tracer bullet target (T_pole = 10⁶ K) is fully ionised.

## D5: Julia native — no runtime dependency on external Fortran/C
**Date:** 2026-03-24
**Decision:** All physics implemented in Julia. External codes (McPHAC, Potekhin Fortran) used only for verification.
**Rationale:** Reproducibility and transparency. Every equation must be traceable in our code.
