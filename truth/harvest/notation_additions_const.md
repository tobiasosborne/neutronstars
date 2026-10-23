# Notation additions: const lane (IAU)

Dictionary row format: `| sym | meaning | units | source entry |`

| sym | meaning | units | source |
|---|---|---|---|
| `\mathrm{au}` | astronomical unit, exactly 149 597 870 700 m (IAU 2012 B2) | m | const.au |
| `\mathrm{pc}` | parsec, exactly 648000/\pi au (IAU 2015 B2) | m | const.parsec |
| `(GM)^{N}_\odot` | nominal solar mass parameter, 1.3271244e20 m^3 s^-2 (conversion constant, IAU 2015 B3); sym_source `\mathcal{(GM)}^{\mathrm N}_\odot` | m^3 s^-2 | const.gm-sun-nominal |
| `R^{N}_\odot` | nominal solar radius, 6.957e8 m (conversion constant); sym_source `\mathcal{R}^{\mathrm N}_\odot` | m | const.r-sun-nominal |
| `L^{N}_\odot` | nominal solar luminosity, 3.828e26 W (conversion constant); sym_source `\mathcal{L}^{\mathrm N}_\odot` | W | const.l-sun-nominal |
| `q_{esu}`, `q_{C}` | charge in Gaussian esu / SI coulomb; q_esu = 10 c[m/s] q_C | esu, C | const.esu-from-coulomb |

Note: `L^{N}_\odot` and `R^{N}_\odot` do not clash with the retired bare `L` (they are decorated). `\odot` rows must not be confused with the observed `R^\infty` / `L_{ph}^\infty`.
