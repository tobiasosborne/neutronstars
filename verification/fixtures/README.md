# Verification fixtures

Published numerical tables used as **verification ground truth only**. Nothing
in `src/` may read from this directory at runtime; see CLAUDE.md §5 for the
explicit exception that scopes these files to `verification/` and `test/`
consumers.

If you are adding a new fixture, add an entry below describing its provenance,
the published reference it derives from, and the verification script that
consumes it.

## `potekhin_tables/` — magnetised hydrogen EOS + Rosseland opacity tables

- **Source:** A. Y. Potekhin, magnetised hydrogen EOS + opacity tables.
  Underlying physics: Potekhin & Chabrier, A&A 374, 213 (2001) (EOS) and
  Potekhin & Chabrier, A&A 399, 1007 (2003) (opacities). See
  `refs/potekhin_chabrier_2001_eos.pdf` and
  `refs/potekhin_chabrier_2003_ff_opacity.pdf`.
- **Consumed by:** `verification/potekhin_table_comparison.jl`, via the parser
  in `src/eos/potekhin_table_reader.jl` (the parser module is in `src/` for
  organisational reasons but is *not* included from `NeutronStar.jl` and is
  loaded only by the verification script — see the header note in that file).
- **Format:** plain-text `.dat` files, one per magnetic-field strength.
  Each file concatenates isotherm blocks; every data row has 14 columns
  (lg R, lg P, PV/(NkT), U/(NkT), S/(Nk), Cv/(Nk), chi_T, chi_rho, x_H,
  x_H0, x_H2, x_pert, lg K_parallel, lg K_perp). Full column definitions and
  the unit conventions live in `potekhin_tables/hmagtab.txt`.
- **File-name convention** (encodes lg B in Gauss):
  - `hmm*.dat`  — lg B = 10.5 … 12.0 (2-line header per isotherm).
    `hmm10a8.dat` = lg B = 10.845, `hmm11a8.dat` = lg B = 11.845.
  - `hmag*.dat` — lg B = 11.9 … 13.5 (5-line header per isotherm).
  - `hmn*.dat`  — lg B = 13.5 … 15.0 (5-line header per isotherm).
- **Grid:** lg R from -7.4 to 3.6 in steps of 0.2 (upper limit 3.4 for
  lg B > 13.5). Temperature grid depends on B band; see `hmagtab.txt`.
- **Size:** ~7.4 MB across 47 `.dat` files + 1 documentation file.
