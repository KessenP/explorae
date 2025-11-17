#!/usr/bin/env python3
# Usage: python explorae.py path/to/master.xlsx path/to/interactions

from __future__ import annotations
import sys, json, re, shlex, subprocess, shutil, os
from pathlib import Path
from typing import Optional, Dict, Tuple, List
import pandas as pd
import argparse

import pyrosetta
from pyrosetta import pose_from_pdb
from pyrosetta import *
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

#from pyrosetta import pose_from_pdb, get_fa_scorefxn


# ---------- CONFIG MODIFIABLE ----------
ID_COL = "jobs"         # colonne d'ID dans l'Excel (doit matcher le nom des dossiers)
SHEET  = 0              # 0 = premier onglet ; ou "Sheet1"
PAE_CUTOFF  = 10
DIST_CUTOFF = 10
IPSae_SCRIPT = Path(__file__).resolve().parent / "ipsae.py"
PRODIGY_CMD = [sys.executable, str((Path(__file__).parent / "cli.py").resolve()), "--showall"]
# --------------------------------------


COL_IPSAE   = "ipsae"
COL_PDOCKQ2 = "pdockq2"
COL_PRODIGY = "prodigy_kd"
COL_IPTM_PTM = "ipTM+pTM"
COL_PTM = "pTM"
COL_REF2015 = "REF2015"
COL_dG_SASA_ratio = "dG_SASA_ratio"

def log(m: str): print(m, flush=True)

def parse_debug_file(dirpath: Path) -> Optional[str]:
    f = dirpath / "ranking_debug.json"
    debug_data={}
    if not f.exists():
        return None
    data = json.loads(f.read_text())
    order = data.get("order") or []
    if order:
        best_model=order[0]
        iptm_ptm=data.get("iptm+ptm")[best_model]
        iptm=data.get("iptm")[best_model]
        debug_data={'model':best_model, 'iptm+ptm':iptm_ptm, 'iptm':iptm}
    return debug_data if order else None

def expected_files(dirpath: Path, top_model: str) -> Tuple[Path, Path]:
    
    pdb = dirpath / f"unrelaxed_{top_model}.pdb"
    pkl = dirpath / f"result_{top_model}.pkl"
    #pdb = os.path.join(dirpath, f"unrelaxed_{top_model}.pdb")
    #pkl = os.path.join(dirpath, f"result_{top_model}.pkl")
    return pdb, pkl

def run_ipsae(pae_pkl: Path, pdb: Path, workdir: Path) -> Tuple[Optional[float], Optional[float]]:
    if not IPSae_SCRIPT.exists():
        log(f"[ERROR] ipsae.py introuvable: {IPSae_SCRIPT}")
        return None, None
    cmd = [sys.executable, str(IPSae_SCRIPT), str(pae_pkl), str(pdb), str(PAE_CUTOFF), str(DIST_CUTOFF)]
    try:
        subprocess.run(cmd, cwd=str(workdir), check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        log(f"[ERROR] ipsae.py a échoué dans {workdir.name}\nSTDERR: {e.stderr.decode(errors='ignore')}")
        return None, None

    # Résumé généré par ipsae.py: <pdb_stem>_<pae>_<dist>.txt
    stem = str(pdb)
    base = stem[:-4] if stem.endswith((".pdb", ".cif")) else str(pdb.with_suffix(""))
    summary = Path(f"{base}_{int(PAE_CUTOFF):02d}_{int(DIST_CUTOFF):02d}.txt")
    if not summary.exists():
        summary = workdir / f"{Path(base).name}_{int(PAE_CUTOFF):02d}_{int(DIST_CUTOFF):02d}.txt"
    if not summary.exists():
        log(f"[ERROR] Résumé ipSAE introuvable: {summary}")
        return None, None

    ipsae_vals: List[float] = []
    pdq2_vals: List[float] = []
    txt = summary.read_text().splitlines()
    for line in txt:
        if " max " in line:
            parts = line.split()
            try:
                if parts[4] == "max":
                    ipsae_vals.append(float(parts[5]))
                    pdq2_vals.append(float(parts[12]))
            except Exception:
                pass
    if not ipsae_vals or not pdq2_vals:
        for line in txt:
            if " asym " in line:
                parts = line.split()
                try:
                    ipsae_vals.append(float(parts[5]))
                    pdq2_vals.append(float(parts[12]))
                except Exception:
                    pass

    ipSAE   = max(ipsae_vals) if ipsae_vals else None
    pDockQ2 = max(pdq2_vals) if pdq2_vals else None
    return ipSAE, pDockQ2

def run_prodigy(pdb_path: Path) -> Optional[float]:
    """
    Exécute le CLI local sur LE FICHIER PDB sélectionné.
    On passe --showall pour garantir l'impression du Kd.
    """
    if not pdb_path.exists():
        log(f"[WARN] PDB introuvable pour PRODIGY: {pdb_path}")
        return None

    cmd = PRODIGY_CMD + [str(pdb_path)]
    try:
        p = subprocess.run(
            cmd, check=True, capture_output=True, text=True,
            cwd=str(pdb_path.parent)
        )
        out = p.stdout + "\n" + p.stderr
    except subprocess.CalledProcessError as e:
        log(f"[WARN] PRODIGY a échoué sur {pdb_path.name}: {e}\nSTDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}")
        return None

    # Parsers tolérants : "dissociation constant (M) at 25.0˚C: 1.23e-06"
    # variantes d'espaces/degrés/majuscule
    patterns = [
        r"dissociation\s+constant.*?:\s*([0-9.+\-eE]+)",
        r"\bkd[_\s]*val\b[^0-9eE+\-]*([0-9.+\-eE]+)",
        r"\bKd\s*\(M\)\s*[:=]\s*([0-9.+\-eE]+)",
    ]
    for pat in patterns:
        m = re.search(pat, out, flags=re.IGNORECASE)
        if m:
            try:
                return float(m.group(1))
            except Exception:
                pass

    # Si on n'a rien matché, loggons la sortie pour diagnostic
    head = "\n".join(out.splitlines()[:40])
    log(f"[INFO] Kd non détecté dans la sortie PRODIGY pour {pdb_path.name}. "
        f"Extrait de sortie:\n{head}\n--- fin extrait ---")
    return None




def update_excel_write_safe(excel: Path, sheet, id_col: str, updates: Dict[str, Dict[str, Optional[float]]]):
    df = pd.read_excel(excel, sheet_name=sheet)

    # Colonnes métriques: crée si absentes
    for c in (COL_IPSAE, COL_PDOCKQ2, COL_PRODIGY):
        if c not in df.columns:
            df[c] = pd.NA

    ok, miss = 0, 0
    for inter_id, vals in updates.items():
        mask = df[id_col].astype(str) == str(inter_id)
        if mask.any():
            idx = df.index[mask][0]
            for k, v in vals.items():
                if v is not None:
                    df.at[idx, k] = v
            ok += 1
        else:
            miss += 1
            log(f"[WARN] ID '{inter_id}' absent de l’Excel")

    # Écriture sûre: fichier temporaire puis remplacement
    tmp_path = excel.with_suffix(".tmp.xlsx")
    try:
        with pd.ExcelWriter(tmp_path, engine="openpyxl", mode="w") as w:
            df.to_excel(w, sheet_name=sheet if isinstance(sheet, str) else "Sheet1", index=False)
        os.replace(tmp_path, excel)  # atomique sous Windows 10+
        log(f"[OK] Excel mis à jour ({ok} lignes, {miss} manquantes)")
    except PermissionError:
        log(
            "[ERROR] Fichier Excel verrouillé (ouvert dans Excel). "
            f"Ferme '{excel.name}' et relance. Un export a été écrit ici: {tmp_path}"
        )

def parse_args():
    parser = argparse.ArgumentParser(
        description="EXPLORAE: enrich AlphaFold-Multimer interaction models with ipSAE, pDockQ2, and PRODIGY Kd, then update a master Excel file."
    )
    parser.add_argument(
        "excel",
        help="Path to the master Excel file (e.g. master.xlsx)"
    )
    parser.add_argument(
        "interactions_root",
        help="Root directory containing one folder per interaction (each with ranking_debug.json, unrelaxed_*.pdb, result_*.pkl, etc.)"
    )
    parser.add_argument(
        "--pae",
        type=float,
        default=PAE_CUTOFF,
        help=f"PAE cutoff (Å) for ipSAE (default: {PAE_CUTOFF})"
    )
    parser.add_argument(
        "--dist",
        type=float,
        default=DIST_CUTOFF,
        help=f"Distance cutoff (Å) for contact definition in ipSAE (default: {DIST_CUTOFF})"
    )
    parser.add_argument(
        "--id-col",
        default=ID_COL,
        help=f"Column name in Excel used as interaction ID (default: '{ID_COL}')"
    )
    parser.add_argument(
        "--sheet",
        default=str(SHEET),
        help=(
            f"Excel sheet index or name (default: {SHEET}). "
            "Use an integer index (0,1,2,...) or a sheet name like 'Sheet1'."
        )
    )
    return parser.parse_args()

def dico_to_csv(dataset, filepath):

   
    # Convertir en DataFrame
    df = pd.DataFrame.from_dict(dataset, orient='index')

    # Ajouter une colonne pour les clés si nécessaire
    df.index.name = 'id'

    print("dataset is:", df)

    # Sauvegarder en CSV
    df.to_csv(filepath, encoding='utf-8')

def calculate_ref2015_energy_terms(pdb_path):
    """
    Calculate REF2015 energy terms for a given PDB file.
    
    Args:
        pdb_path (str): Path to the PDB file.
    
    Returns:
        float total score value
    """
    # Initialize PyRosetta if not already done
    pyrosetta.init()
    
    # Load the pose from PDB
    pose = pose_from_pdb(pdb_path)
    
    # Set up the REF2015 scoring function
    scorefxn = get_fa_scorefxn()  # Default is REF2015 in PyRosetta
    
    # Score the pose
    total_score = scorefxn(pose)
    
    return total_score



def analyze_interface(pdb_file, chains_partner1, chains_partner2):
    """
    Analyse l'interface entre deux ensembles de chaînes dans un complexe.
    
    Arguments :
        pdb_file (str) : chemin vers le fichier PDB
        chains_partner1 (str) : ex. "A" ou "AHL"
        chains_partner2 (str) : ex. "B" ou "CD"
    
    Retourne :
        dict : métriques calculées (dSASA, dG, dG_norm)
    """

    # Initialize PyRosetta if not already done
    #pyrosetta.init()

    # Charger la pose
    pose = pose_from_pdb(pdb_file)

    # Format Rosetta attendu : "A_B" ou "AHL_CD"
    interface_str = f"{chains_partner1}_{chains_partner2}"

    # Lancer l'analyse
    iam = InterfaceAnalyzerMover(interface_str, False)
    iam.apply(pose)

    # Extraction des valeurs stockées dans la pose
    scores = pose.scores

     # Les noms corrects dans ta version de Rosetta
    dG = scores.get("dG_cross", None) # énergie de liaison
    dSASA = scores.get("dSASA_int", None) # dSASA enterrée

    # Normalisation
    dG_norm = dG / dSASA if (dG is not None and dSASA not in (0, None)) else None

    return {
        "interface": interface_str,
        "dG": dG,
        "dSASA": dSASA,
        "dG_norm": dG_norm
    }


def main():
    # Initialize PyRosetta if not already done
    pyrosetta.init()

    args = parse_args()

    # On résout les chemins
    excel = Path(args.excel).resolve()
    root  = Path(args.interactions_root).resolve()

    #if not excel.exists():
    #    sys.exit(f"[ERROR] Excel introuvable: {excel}")
    if not root.exists():
        sys.exit(f"[ERROR] Dossier interactions introuvable: {root}")

    # On met à jour les globals avec les valeurs éventuellement overridées
    global PAE_CUTOFF, DIST_CUTOFF, ID_COL, SHEET
    PAE_CUTOFF = args.pae
    DIST_CUTOFF = args.dist
    ID_COL = args.id_col

    # sheet : si l'utilisateur donne un entier -> index, sinon nom
    sheet_arg = args.sheet
    if isinstance(sheet_arg, str) and sheet_arg.isdigit():
        SHEET = int(sheet_arg)
    else:
        SHEET = sheet_arg  # peut être str (nom d’onglet) ou déjà int

    log(f"[INFO] Excel        : {excel}")
    log(f"[INFO] Interactions : {root}")
    log(f"[INFO] SHEET        : {SHEET}")
    log(f"[INFO] ID_COL       : {ID_COL}")
    log(f"[INFO] PAE_CUTOFF   : {PAE_CUTOFF}")
    log(f"[INFO] DIST_CUTOFF  : {DIST_CUTOFF}")

    updates: Dict[str, Dict[str, Optional[float]]] = {}

    

    for inter_dir in sorted([p for p in root.iterdir() if p.is_dir()]):
        inter_id = inter_dir.name
        log(f"\n=== {inter_id} ===")
        debug_data = parse_debug_file(inter_dir)
        #debug_data={'model':best_model, 'iptm+ptm':iptm_ptm, 'iptm':iptm}

        top=debug_data['model']
        if not top:
            log("[WARN] ranking_debug.json absent ou sans 'order' -> skip")
            continue

        pdb, pkl = expected_files(inter_dir, top)
        if not pdb.exists() or not pkl.exists():
        #if (not pdb or not pkl):
            log(f"[WARN] Fichiers manquants pour {top}: pdb={pdb.exists()} pkl={pkl.exists()} -> skip")
            #log(f"[WARN] Fichiers manquants pour {top}: pdb={pdb} pkl={pkl} -> skip")
            continue

        ipsae_val, pdq2_val = run_ipsae(pkl, pdb, inter_dir)
        kd_val = run_prodigy(pdb)

        # calculate_interface_dG_dSASA_ratio attend une string
        # in our case, chains interacting are B and C
        dG_SASA_results=analyze_interface(str(pdb), "B", "C")

        if ipsae_val is not None: log(f"ipSAE   : {ipsae_val:.6f}")
        if pdq2_val  is not None: log(f"pDockQ2 : {pdq2_val:.6f}")
        if kd_val    is not None: log(f"PRODIGY Kd (M): {kd_val:.3e}")
        if debug_data['iptm+ptm']    is not None: log(f"iptm+ptm: {debug_data['iptm+ptm']:.3e}")
        if debug_data['iptm']     is not None: log(f"iptm : {debug_data['iptm'] :.3e}")
        if dG_SASA_results["dG_norm"] is not None: log(f"dG_SASA_ratio (kJ/A): {dG_SASA_results["dG_norm"]:.3e}")

        
        updates[inter_id] = {
            COL_IPSAE:   ipsae_val,
            COL_PDOCKQ2: pdq2_val,
            COL_PRODIGY: kd_val,
            COL_IPTM_PTM: debug_data['iptm+ptm'], 
            COL_PTM : debug_data['iptm'],
            COL_dG_SASA_ratio : dG_SASA_results["dG_norm"]
        }

    #update_excel_write_safe(excel, SHEET, ID_COL, updates)
    #with open("data/interactions_négatives.json", "w", encoding="utf-8") as f:
    #    json.dump(updates, f)

    # write to csv file
    dico_to_csv(updates, excel)

if __name__ == "__main__":
    main()