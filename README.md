# EXPLORAE — Automated Scoring and Integration of Protein–Protein Interaction Metrics

EXPLORAE computes, aggregates and exports multiple structural‑confidence, biophysical, and energy‑based metrics from AlphaFold‑Multimer predictions into an Excel file. It scans interaction model folders produced by AlphaFold‑Multimer, computes interface confidence and energy metrics (ipSAE, pDockQ2, ipTM, pTM, PRODIGY predictions, Rosetta interface analysis), and writes the results into a master Excel/CSV sheet for downstream analysis.

---

## METRICS TABLE (MANDATORY)

The table below lists the metrics produced by EXPLORAE, their units, source, and short meaning. The Excel output uses the column names shown in the table; units are reported in the table and are the canonical units for each metric.

| Metric | Unit | Source | Meaning |
|--------|------|--------|---------|
| ipSAE | dimensionless | AlphaFold (PAE → transform) | Interface confidence |
| pDockQ2 | dimensionless | AF + logistic regression | Interface correctness score |
| ipTM | dimensionless | AlphaFold | Interface TM-score |
| pTM | dimensionless | AlphaFold | Global TM-score |
| ipTM+pTM | dimensionless | AlphaFold | Combined ranking score |
| PRODIGY Kd | M | PRODIGY | Predicted dissociation constant |
| PRODIGY ΔG internal | kcal/mol | PRODIGY | Binding free energy |
| Rosetta dG_cross | REU | Rosetta | Interface interaction energy |
| Rosetta dSASA_int | Å² | Rosetta | Buried SASA |
| dG_SASA_ratio | REU/Å² | Rosetta | Normalized interface energy |

Excel output column names (examples used in the code):
- `ipsae` (ipSAE, dimensionless)
- `pdockq2` (pDockQ2, dimensionless)
- `prodigy_kd` (PRODIGY Kd, M)
- `prodigy_dg_internal` (PRODIGY ΔG internal, kcal/mol)
- `dG_rosetta` (Rosetta dG_cross, REU)
- `dG_SASA_ratio` (dG_cross / dSASA_int, REU/Å²)
- `ipTM+pTM` (combined score)
- `pTM` (pTM/ipTM entries)

---

## EXPLANATION OF METRICS

### 1. AlphaFold confidence metrics
- What they measure:
  - ipSAE: an interface-centric transform of AlphaFold Predicted Aligned Error (PAE). It uses a PTM‑like transform to convert PAE values into a 0–1 confidence per residue pair and aggregates them to give an interface confidence score.
  - pDockQ2: a data‑driven interface correctness score combining local confidence (pLDDT) and transformed PAE; mapped through an empirically fitted logistic function.
  - ipTM: an interface TM‑score proxy derived from PAE using the ptm transformation used by ipSAE/pTM proxies.
  - pTM: AlphaFold’s predicted TM (global confidence) as produced by the AlphaFold pipeline.
- Statistical/physical/confidence-based:
  - ipSAE: confidence‑based (derived from predicted errors, empirical transform).
  - pDockQ2: statistical/regression (empirically fitted logistic model combining confidence and PAE features).
  - ipTM & pTM: confidence‑based, proxies of structural correctness derived from predicted error matrices or AlphaFold outputs.

Notes on calculation (implementation in `src/ipsae.py`):
- ptm transform: ptm_func(x, d0) = 1 / (1 + (x/d0)^2), where x is PAE and d0 is an empirically chosen scale (often a function of chain length or a fixed value).
- ipSAE: per‑residue means of ptm_func aggregated across interface pairs, then the maximum residue value for a pair is selected and symmetrical max is used to report a pair value.
- pDockQ2: mean Cβ‑pLDDT × mean_ptm (ptm computed with d0=10.0), then mapped by:
  - pDockQ2 = 1.31 / (1 + exp(-0.075*(x - 84.733))) + 0.005
  where x = mean_plddt * mean_ptm (constants from the implemented logistic fit).

---

### 2. PRODIGY metrics
- What they measure:
  - PRODIGY ΔG internal (`ba_val`): an empirically predicted binding free energy (reported in kcal/mol) produced by a linear regression over intermolecular contact counts and interface amino‑acid composition metrics.
  - PRODIGY Kd (`kd_val`): a predicted dissociation constant (M) derived from the PRODIGY ΔG internal using the thermodynamic relation used in the code.
- Statistical/physical/confidence-based:
  - PRODIGY ΔG internal: statistical/regression (IC_NIS linear model).
  - Kd: numeric transform of the regression output using standard thermodynamic relation.
- Exact formulas (from the local implementation `src/modules/models.py` and `src/modules/utils.py`):
  - IC_NIS linear model:
    - ba_val = -0.09459*CC + -0.10007*AC + 0.19577*PP + -0.22671*AP + 0.18681*nis_a + 0.13810*nis_c - 15.9433
      - CC, AC, PP, AP are counts of contact types (charged/ apolar/ polar combinations)
      - nis_a, nis_c are % apolar and % charged NIS residues (0–100)
      - ba_val unit: kcal/mol (reported as PRODIGY ΔG internal)
  - Conversion to Kd (function `dg_to_kd`):
    - temp_in_K = T(°C) + 273.15
    - RT (kcal·mol⁻1·K⁻1) constant used: 0.0019858775
    - kd = exp( ba_val / (RT) )
      - Thus kd (M) = exp(ΔG / (R·T)) according to the local code.
    - Inverse relation (thermodynamic):
      - ΔG = RT * ln(Kd)
    - Conversions to kJ: 1 kcal = 4.184 kJ (if needed).
- Notes:
  - PRODIGY ΔG is an empirical energy; kd is derived via the exponential transform. Check sign conventions (the code uses kd = exp(ΔG/RT) consistent with the implementation of dg_to_kd).

---

### 3. Rosetta metrics
- What they measure:
  - Rosetta `dG_cross`: Rosetta’s interface energy (REU), typically computed as Score(complex) - Score(partnerA) - Score(partnerB) (empirical Rosetta energy units).
  - Rosetta `dSASA_int`: buried solvent‑accessible surface area (Å²) at the interface.
  - `dG_SASA_ratio`: dG_cross divided by dSASA_int (units: REU/Å²); a normalized indicator of energy per buried surface area.
- Statistical/physical/confidence-based:
  - Rosetta metrics are physics‑inspired empirical energy terms (force‑field based) and surface calculations; they are not strictly thermodynamic ΔG in kcal/mol and are best interpreted relatively (REU).
- Calculation notes (as used by `explorae.py`):
  - The code uses PyRosetta InterfaceAnalyzerMover:
    - `InterfaceAnalyzerMover(interface_str, False)` applied to a `pose`.
    - After application, `pose.scores` contains `"dG_cross"` and `"dSASA_int"`.
  - `dG_SASA_ratio` is computed as `dG_cross / dSASA_int` if `dSASA_int` is nonzero.
- Units:
  - dG_cross: REU (Rosetta Energy Units) — not directly kcal/mol.
  - dSASA_int: Å².
  - dG_SASA_ratio: REU/Å².

---

## INSTALLATION

Clone the repository, create a Python virtual environment, install Python dependencies, and optionally install PyRosetta.

```bash
# Clone the repo
git clone https://github.com/adboussif/explorae.git
cd explorae

# Create and activate a Python virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install Python requirements
pip install --upgrade pip
pip install -r requirements.txt

# Optional: install PyRosetta (if you plan to compute Rosetta metrics)
# Example using conda channels (adjust depending on your provider / licenses):
conda create -n explorae_pyro python=3.10 -y
conda activate explorae_pyro
# Example install (user must have access to PyRosetta channel/package):
conda install -c rosetta -c conda-forge pyrosetta -y