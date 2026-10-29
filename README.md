# QMThermoProperties.jl

[![Build Status](https://github.com/se-schmitt/QMThermoProperties.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/se-schmitt/QMThermoProperties.jl/actions/workflows/CI.yml?query=branch%3Amain)

Thermodynamic models built from quantum-mechanical (QM) molecular data in the
[CHAOS database](https://doi.org/10.48550/arXiv.2511.19002), for use with [Clapeyron.jl](https://github.com/ClapeyronThermo/Clapeyron.jl).

| Model | Type | What it gives you |
|---|---|---|
| [`QMIdeal`](#qmideal-ideal-gas-model) | `Clapeyron.IdealModel` | Ideal-gas properties (`cₚ`, `h`, `s`, `g`, …) from statistical mechanics |
| [`COSMOSAC02_CHAOS`](#cosmo-sac-activity-coefficient-models) | `Clapeyron.COSMOSAC02` | COSMO-SAC (2002) activity coefficients with CHAOS sigma profiles |
| [`COSMOSAC10_CHAOS`](#cosmo-sac-activity-coefficient-models) | `Clapeyron.COSMOSAC10` | COSMO-SAC (2010) activity coefficients with CHAOS sigma profiles |

All models are regular Clapeyron models: every Clapeyron property and equilibrium
function works with them, and `QMIdeal` can be used as the ideal part of any Clapeyron
equation of state.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/se-schmitt/QMThermoProperties.jl")
Pkg.add(["Clapeyron", "RDKitMinimalLib"])   # RDKitMinimalLib is optional, see below
```

Before building the first model, download the CHAOS data once (~60 MB):

```julia
using QMThermoProperties
QMThermoProperties.download_chaos()
```

The data is stored in a Julia [scratch space](https://github.com/JuliaPackaging/Scratch.jl)
and persists across sessions and package updates. See [CHAOS database](#chaos-database)
for details.

### RDKit (optional)

The CHAOS database is keyed by RDKit-canonical SMILES. When
[RDKitMinimalLib.jl](https://github.com/eloyfelix/RDKitMinimalLib.jl) is loaded, the
package canonicalizes the SMILES of your components automatically and can compute
missing molar masses:

```julia
using QMThermoProperties, Clapeyron, RDKitMinimalLib
```

Without it, SMILES are used as given (with a one-time warning), so they must already be
RDKit-canonical. RDKitMinimalLib.jl does not work on Windows.

## Quick start

```julia
using QMThermoProperties, Clapeyron, RDKitMinimalLib

# ideal gas heat capacity of an equimolar water/ethanol mixture
ideal = QMIdeal(["water", "ethanol"])
isobaric_heat_capacity(ideal, 1e5, 298.15, [0.5, 0.5])    # J mol⁻¹ K⁻¹

# QM ideal-gas part inside a cubic equation of state
model = PR(["water", "ethanol"]; idealmodel = QMIdeal)
enthalpy(model, 1e5, 400.0, [0.5, 0.5])

# activity coefficients from COSMO-SAC with CHAOS sigma profiles
act = COSMOSAC10_CHAOS(["water", "ethanol"])
activity_coefficient(act, 1e5, 298.15, [0.5, 0.5])        # ≈ [1.379, 1.115]
bubble_pressure(act, 298.15, [0.5, 0.5])
```

## Specifying components

Every model looks up its components in CHAOS by canonical SMILES. There are three ways to
supply them:

1. **By name.** Names are resolved through Clapeyron's identifier database
   (`properties/identifiers.csv`), which provides the SMILES:
   ```julia
   QMIdeal(["water", "carbon dioxide"])
   ```
2. **By SMILES via Clapeyron's `smiles` string macro**, which also searches the identifier
   database (only works for SMILES listed there):
   ```julia
   QMIdeal(smiles"CCCO")
   COSMOSAC02_CHAOS([smiles"CCO", smiles"O"])
   ```
3. **Via `userlocations`**, for any molecule in CHAOS and with any name:
   ```julia
   QMIdeal("my propanol"; userlocations = (; SMILES = ["CCCO"]))

   # skip canonicalization (no RDKit needed): give the canonical SMILES directly
   QMIdeal("my propanol"; userlocations = (; canonicalsmiles = ["CCCO"], Mw = [60.096]))
   ```
   `userlocations` can also be a list of paths to CSV files in Clapeyron's format.

Recognized parameters:

| Parameter | Used by | Meaning |
|---|---|---|
| `SMILES` | all | SMILES string, canonicalized with RDKit before lookup |
| `canonicalsmiles` | all | RDKit-canonical SMILES, used as is (takes precedence over `SMILES`) |
| `Mw` | `QMIdeal` | molar mass in g mol⁻¹; if missing it is computed with RDKit |

A component that is not in CHAOS, or whose QM calculation is flagged as not converged in
CHAOS, raises an error naming its SMILES.

## `QMIdeal`: ideal gas model

```julia
QMIdeal(components;
    userlocations = String[],
    reference_state = nothing,
    verbose = false,
    vib_correction = 0.95461)
```

Ideal-gas Helmholtz energy from the molecular partition function in the rigid-rotor /
harmonic-oscillator approximation. All molecular data (vibrational frequencies, rotational
constants, symmetry numbers, spin multiplicities) come from the QM calculations in CHAOS,
so no fitted parameters are involved.

For a mixture with amounts $n_i$ in volume $V$, the reduced Helmholtz energy is

$$
\frac{A^\mathrm{ig}}{nRT} = \sum_i x_i \left[ \ln\left(\frac{n_i N_\mathrm{A}}{V}\Lambda_i^3\right) - \ln q_{\mathrm{rot},i} - \ln q_{\mathrm{vib},i} - \ln q_{\mathrm{el},i} \right] - 1,
\qquad
\Lambda_i = \frac{h}{\sqrt{2\pi m_i k_\mathrm{B} T}}
$$

with the molecular partition functions

| Contribution | Expression | Data from CHAOS |
|---|---|---|
| Translation | via $\Lambda_i$, $m_i = M_{w,i}/N_\mathrm{A}$ | — (`Mw` from Clapeyron/RDKit/user) |
| Rotation, nonlinear | $q_\mathrm{rot} = \dfrac{\sqrt{\pi}}{\sigma}\left(\dfrac{8\pi^2 k_\mathrm{B}T}{h^2}\right)^{3/2}\sqrt{I_A I_B I_C}$ | rotational constants $B_k$, $I_k = h/(8\pi^2 B_k)$; symmetry number $\sigma$ |
| Rotation, linear | $q_\mathrm{rot} = \dfrac{8\pi^2 I k_\mathrm{B}T}{\sigma h^2}$ | as above (one constant) |
| Vibration | $q_\mathrm{vib} = \prod_k \dfrac{\exp(-h\nu_k / 2k_\mathrm{B}T)}{1-\exp(-h\nu_k / k_\mathrm{B}T)}$ | harmonic frequencies $\tilde\nu_k$, $\nu_k = s\,c_0\tilde\nu_k$ |
| Electronic | $q_\mathrm{el} = \omega_0$ | spin multiplicity of the ground state |

$s$ is the frequency scaling factor `vib_correction` (default `0.95461`, see [Gond et al. (2025)](https://doi.org/10.48550/arXiv.2511.19002)). Use `vib_correction = 1` for unscaled frequencies.

The vibrational partition function is referenced to the bottom of the potential well, so
the internal energy contains the zero-point energy. The electronic energy itself is not
included: absolute energies are relative to the electronic ground state of each molecule,
not to the elements. Use a reference state (below) to set a common zero.

### Usage

```julia
using QMThermoProperties, Clapeyron

ideal = QMIdeal("ethanol")

isobaric_heat_capacity(ideal, 1e5, 300.0)   # J mol⁻¹ K⁻¹
enthalpy(ideal, 1e5, 300.0)                 # J mol⁻¹
entropy(ideal, 1e5, 300.0)                  # J mol⁻¹ K⁻¹
gibbs_free_energy(ideal, 1e5, 300.0)        # J mol⁻¹
```

As the ideal part of an equation of state, pass the type (default settings) or a
constructed model (custom settings):

```julia
PCSAFT("ethanol"; idealmodel = QMIdeal)
PCSAFT("ethanol"; idealmodel = QMIdeal("ethanol"; vib_correction = 1.0))
```

**Reference states.** Any Clapeyron reference state can be applied. Reference states that
need a liquid phase (e.g. `:ashrae`, `:nbp`, `:iir`) only make sense on the full equation of
state, not on `QMIdeal` alone:

```julia
PR("ethanol"; idealmodel = QMIdeal, reference_state = :ashrae)

# for the ideal gas alone: fix h and s at a chosen (T0, P0)
QMIdeal("ethanol"; reference_state = ReferenceState(:volume; T0 = 298.15, P0 = 1e5, H0 = 0.0, S0 = 0.0, phase = :vapour))
```

### Zero-point energy

```julia
zpe(model, z = [1.0])
```

Returns the zero-point energy in J of an ideal gas model for the amounts `z` in mol,
computed as $A^\mathrm{ig}(T \to 0)$:

```julia
zpe(QMIdeal("water"))                              # J mol⁻¹
zpe(QMIdeal(["water", "ethanol"]), [0.5, 0.5])     # J for 0.5 mol each
```

### Validation

The test suite compares zero-point energy, heat capacity, entropy and thermal energy
against Gaussian thermochemistry output (`vib_correction = 1`) for molecules from O₂ to
pentacyclic aromatics, with agreement better than 0.05 %.

## COSMO-SAC activity coefficient models

```julia
COSMOSAC02_CHAOS(components;
    puremodel = PR,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)

COSMOSAC10_CHAOS(components; ...)     # same keywords
```

These construct Clapeyron's `COSMOSAC02` and `COSMOSAC10` activity coefficient models
with parameters taken from CHAOS instead of the default sigma-profile database. This makes
COSMO-SAC available for every molecule in CHAOS.

| Model | Method | Parameters from CHAOS |
|---|---|---|
| `COSMOSAC02_CHAOS` | Lin & Sandler (2002) | total sigma profile `Pi`, cavity area `A` (Å²), cavity volume `V` (Å³) |
| `COSMOSAC10_CHAOS` | Hsieh, Sandler & Lin (2010) | sigma profiles split into non-hydrogen-bonding (`Pnhb`), hydroxyl (`POH`) and other hydrogen-bonding (`POT`) segments, `A`, `V` |

Sigma profiles use 51 bins from −0.025 to 0.025 e Å⁻². The models use the standard
COSMO-SAC parameters implemented in Clapeyron.

`puremodel` is the equation of state for the pure components, used for vapour pressures
and fugacities in phase equilibrium calculations (`pure_userlocations` goes to it).

```julia
using QMThermoProperties, Clapeyron, RDKitMinimalLib

model = COSMOSAC10_CHAOS(["toluene", "ethanol"])

activity_coefficient(model, 1e5, 300.0, [0.3, 0.7])
excess(model, 1e5, 300.0, [0.3, 0.7], gibbs_free_energy)   # J mol⁻¹
bubble_pressure(model, 300.0, [0.3, 0.7])

# different pure-component model, components given by SMILES
COSMOSAC02_CHAOS(["toluene", "ethanol"];
    userlocations = (; SMILES = ["CC1=CC=CC=C1", "CCO"]),
    puremodel = PCSAFT)
```

## CHAOS database

```julia
QMThermoProperties.download_chaos(; force = false, raw = false)
```

Makes the CHAOS data available. Three files are kept in the package scratch space: QM
data for `QMIdeal`, and sigma profiles for the 2002 and 2010 COSMO-SAC variants. They are
loaded lazily, at the first model construction in a session.

- Default: downloads the converted databases (~60 MB) from the GitHub release and verifies
  their SHA256 checksum.
- `force = true`: replaces an existing database.
- `raw = true`: instead downloads the full raw CHAOS data from Zenodo (~4.3 GB, needs
  ~16 GB of temporary disk space) and converts it locally. This requires RDKitMinimalLib.jl.
  It is only needed to rebuild the converted databases (maintainers can package the
  result with `QMThermoProperties.bundle_chaos()`).

The data is removed together with the package by `Pkg.gc()`.

## References

1. D. Gond, J. Arweiler, T. Specht, H. Hasse, and F. Jirasek: CHAOS – A Consistent Large-Scale Database for Sigma-Profiles and Other Molecular Descriptors (2025), DOI: [10.48550/arXiv.2511.19002](https://doi.org/10.48550/arXiv.2511.19002).
2. A. Klamt: Conductor-like Screening Model for Real Solvents: A New Approach to the Quantitative Calculation of Solvation Phenomena, J. Phys. Chem. 99 (1995) 2224–2235, DOI: [10.1021/j100007a062](https://doi.org/10.1021/j100007a062).
3. S.-T. Lin and S. I. Sandler: A Priori Phase Equilibrium Prediction from a Segment Contribution Solvation Model, Ind. Eng. Chem. Res. 41 (2002) 899–913, DOI: [10.1021/ie001047w](https://doi.org/10.1021/ie001047w).
4. C.-M. Hsieh, S. I. Sandler, and S.-T. Lin: Improvements of COSMO-SAC for Vapor–Liquid and Liquid–Liquid Equilibrium Predictions, Fluid Phase Equilib. 297 (2010) 90–97, DOI: [10.1016/j.fluid.2010.06.011](https://doi.org/10.1016/j.fluid.2010.06.011).

If you use this package, please cite the CHAOS paper (1) and, for the COSMO-SAC models, the
corresponding method papers.

## License

MIT, see [LICENSE](LICENSE).
