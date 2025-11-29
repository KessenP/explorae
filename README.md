# EXPLORAE — Automated Scoring and Integration of Protein–Protein Interaction Metrics

EXPLORAE computes, aggregates and exports multiple structural‑confidence, biophysical, and energy‑based metrics from AlphaFold‑Multimer predictions into an Excel file. It scans interaction model folders produced by AlphaFold‑Multimer, computes interface confidence and energy metrics (ipSAE, pDockQ2, ipTM, pTM, PRODIGY predictions, Rosetta interface analysis), and writes the results into a master Excel/CSV sheet for downstream analysis.

---

## METRICS TABLE

The table below lists the metrics produced by EXPLORAE, their units, source, and short meaning. The Excel output uses the column names shown in the table; units are reported in the table and are the canonical units for each metric.

---

## EXPLANATION OF METRICS

### 1. AlphaFold confidence metrics
- What they measure:
  - ipSAE: an interface-centric transform of AlphaFold Predicted Aligned Error (PAE). It uses a PTM‑like transform to convert PAE values into a 0–1 confidence per residue pair and aggregates them to give an interface confidence score.
  - pDockQ2: a data‑driven interface correctness score combining local confidence (pLDDT) and transformed PAE; mapped through an empirically fitted logistic function.
- Statistical/physical/confidence-based:
  - ipSAE: confidence‑based (derived from predicted errors, empirical transform).
  

  ---

  ## Table des métriques (rapide)

  | Metric | Unité | Source | Sens |
  |--------|-------|--------|------|
  | ipSAE |  | AlphaFold (PAE → transform) | confiance sur l'interface |
  | pDockQ2 |  | AF + régression | score de qualité d'interface |
  | ipTM+pTM |  | AlphaFold | score combiné pour le classement |
  | PRODIGY Kd | M | PRODIGY | constante de dissociation prédite |
  | PRODIGY ΔG internal | kcal/mol | PRODIGY | énergie libre prédite (empirique) |
  | Rosetta dG_cross | REU | Rosetta (PyRosetta) | énergie d'interface (score Rosetta) |
  | Rosetta dSASA_int | Å² | Rosetta | surface enfouie (SASA) |
  | dG_SASA_ratio | REU/Å² | Rosetta | énergie normalisée par surface |

  Les noms de colonnes écrits dans l'Excel sont : `ipsae`, `pdockq2`, `prodigy_kd`, `prodigy_dg_internal`, `dG_rosetta`, `dG_SASA_ratio`, `ipTM+pTM`.

  ---

  ## Explication des groupes de métriques

  ### AlphaFold (confiance)
  - ipSAE : transforme la matrice PAE pour donner une confiance locale sur l'interface. C'est adimensionnel et entre 0 et 1 (plus grand = meilleur).
  - pDockQ2 : combine pLDDT local (confiance) et PAE transformé.
  - ipTM / pTM : mesures de type TM issues des sorties AlphaFold (ici on lit `ranking_debug.json` pour `iptm` et `iptm+ptm`).

  Ces métriques sont basées sur la confiance / erreurs prédites par AlphaFold.

  ### PRODIGY
  - PRODIGY calcule d'abord un `ΔG` empirique (appelé `ba_val`) via un modèle linéaire sur les contacts et la composition d'interface. Cette valeur est en kcal/mol.
  - Ensuite, `Kd` est dérivé par la relation utilisée dans le code : `kd = exp(ΔG / (R*T))` (R en kcal·mol⁻1·K⁻1). Inversement `ΔG = R*T*ln(Kd)`.

  ### Rosetta (PyRosetta)
  - `dG_cross` : score d'interface renvoyé par `InterfaceAnalyzerMover` (unités REU = Rosetta Energy Units).
  - `dSASA_int` : surface enfouie (Å²) calculée par Rosetta.
  - `dG_SASA_ratio` : `dG_cross / dSASA_int` (REU/Å²), utile pour normaliser l'énergie par surface.

  Ces valeurs viennent de Rosetta (score empirique, utile en comparaisons internes) — elles ne sont pas directement en kcal/mol.

  ---

  ## Installation rapide

  ```bash
  git clone https://github.com/adboussif/explorae.git
  cd explorae
  python3 -m venv .venv
  source .venv/bin/activate
  pip install -r requirements.txt


  ## Usage de base

  ```bash
  python ./src/explorae.py ./test.xlsx ./interactions
  ```

  Ce que fait le script :
  - parcourt `./interactions` (un dossier par interaction),
  - lit `ranking_debug.json` pour retrouver le modèle top et ipTM/pTM,
  - calcule ipSAE / pDockQ2 via `ipsae.py` (fichiers de résumé),
  - lance PRODIGY ( CLI + parsing),
  - lance PyRosetta InterfaceAnalyzer si PyRosetta est installée,
  - met à jour l'Excel/CSV en ajoutant les colonnes (ou les crée si manquantes).

  ---

  ## Options (principales)

  | Option | Usage |
  |--------|-------|
  | `--pae` | seuil PAE (Å) pour ipSAE (par défaut `10`) |
  | `--dist` | cutoff distance (Å) pour ipSAE (par défaut `10`) |
  | `--id-col` | nom de la colonne ID dans l'Excel (par défaut `jobs`) |
  | `--sheet` | onglet Excel (index ou nom, par défaut `0`) |

  Exemple :
  ```bash
  python ./src/explorae.py master.xlsx ./interactions --pae 12 --dist 8 --id-col jobs --sheet 0
  ```

  ---

  ## Structure du dépôt

  ```
  explorae/
  ├─ src/
  │  ├─ explorae.py    # script principal
  │  ├─ ipsae.py       # calcule ipSAE, pDockQ2, etc.
  │  ├─ cli.py         # wrapper PRODIGY
  │  └─ modules/       # code utilitaire (prodigy, parsers, utils...)
  ├─ interactions/     # y mettre un dossier par interaction
  └─ requirements.txt
  ```

  ---

  ## Exemple de sortie console

  ```
  === INTERACTION_XYZ ===
  ipSAE   : 0.78
  pDockQ2 : 0.41
  PRODIGY Kd (M): 2.3e-07
  PRODIGY ΔG internal : -8.9 kcal/mol
  iptm+ptm: 1.10
  pTM : 0.55
  dG_SASA_ratio : -0.011 (REU/Å²)
  dG Rosetta (dG_cross) : -38.5 REU
  ```

  ---

