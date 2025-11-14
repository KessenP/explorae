# 🧬 EXPLORAE — Automated Scoring and Excel Integration of Interaction Metrics

**EXPLORAE** is a command-line tool that automatically computes and integrates **interaction quality metrics** (`ipSAE`, `pDockQ2`, `PRODIGY Kd`) for AlphaFold-Multimer predictions.  
It scans interaction subdirectories, selects the top-ranked model from each `ranking_debug.json`, performs metric calculations, and updates an Excel summary file with the results.

---

## 📦 Installation

Clone the repository and create a virtual environment:

```bash
git clone https://github.com/adboussif/explorae.git
cd explorae
python -m venv venv
source venv/bin/activate      # (or: venv\Scripts\activate on Windows)
pip install -r requirements.txt
```

---

## ⚙️ Repository Structure

```
explorae/
│
├── src/
│   ├── explorae.py           # Main CLI script
│   ├── cli.py                # Local PRODIGY implementation
│   ├── ipsae.py              # ipSAE + pDockQ2 metrics computation
│   ├── modules/              # Local dependencies
│   │   ├── prodigy.py
│   │   ├── freesasa_tools.py
│   │   ├── aa_properties.py
│   │   ├── models.py
│   │   ├── utils.py
│   │   └── data/
│   │       ├── naccess.config
│   │       └── __init__.py
│   └── __init__.py
│
├── interactions/             # Folder with all interaction subfolders
│   ├── A0A3Q7E6L3_SOLLC_and_RdRPL/
│   │   ├── ranking_debug.json
│   │   ├── unrelaxed_model_X_multimer_v3_pred_0.pdb
│   │   └── result_model_X_multimer_v3_pred_0.pkl
│   └── ...
│
└── test.xlsx                 # Excel summary file to update
```

---

## 🚀 Usage

### Basic usage (recommended)
By default, the tool looks for the top-ranked model in each interaction folder and adds the computed metrics to the Excel file:

```bash
python ./src/explorae.py ./test.xlsx ./interactions
```

This will:
- Detect each interaction subfolder under `./interactions/`
- Select the best-ranked model from `ranking_debug.json`
- Compute `ipSAE`, `pDockQ2`, and `PRODIGY Kd` for that model
- Update the Excel file (`test.xlsx`) in place, adding or updating the columns:
  - `ipsae`
  - `pdockq2`
  - `prodigy_kd`

---

### Command-line options

| Option        | Default   | Description |
|---------------|-----------|-------------|
| `--pae FLOAT` | `10`      | PAE cutoff (Å) used for ipSAE computation |
| `--dist FLOAT`| `10`      | Distance cutoff (Å) for defining contacts in ipSAE |
| `--id-col STR`| `"jobs"`  | Column name in the Excel file used to match interaction IDs |
| `--sheet S`   | `0`       | Excel sheet index (`0`, `1`, …) or sheet name (e.g. `"Sheet1"`) |

Examples:

```bash
# Change only PAE cutoff
python ./src/explorae.py ./test.xlsx ./interactions --pae 8

# Change PAE and distance cutoffs
python ./src/explorae.py ./test.xlsx ./interactions --pae 12 --dist 6

# Use a different Excel sheet and ID column
python ./src/explorae.py ./test.xlsx ./interactions --sheet AlphaPulldown --id-col Interaction_ID

### Example output

Console output example:

```
=== A0A3Q7E6L3_SOLLC_and_RdRPL ===
ipSAE   : 0.035683
pDockQ2 : 0.147800
PRODIGY Kd (M): 1.23e-06
[OK] Excel mis à jour (2 lignes, 0 manquantes)
```

---

## 🧩 Notes

- `ranking_debug.json` **must exist** in each interaction folder. The first entry in `"order"` determines which model (`unrelaxed_...pdb`, `result_...pkl`) is used.
- The script will skip folders missing either `.pdb`, `.pkl`, or the `ranking_debug.json` file.
- If the Excel file is open in Excel, the update will fail with `PermissionError`; close it and rerun the script.

---
