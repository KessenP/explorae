#!/usr/bin/env python3
# Usage: python explorae.py path/to/master.xlsx path/to/interactions

from __future__ import annotations
import sys, json, re, shlex, subprocess, shutil, os
from pathlib import Path
from typing import Optional, Dict, Tuple, List, Any
import pandas as pd
import argparse
import io
import contextlib

# Import PyRosetta optionnel
try:
    import pyrosetta
    from pyrosetta import pose_from_pdb, get_fa_scorefxn
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
    HAS_PYROSETTA = True
except ImportError:
    HAS_PYROSETTA = False
    log_buffer = []  # Pour stocker les messages avant que log() soit défini
    log_buffer.append("[WARN] PyRosetta non installé. Les métriques Rosetta (dG, dSASA) ne seront pas calculées.")


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
COL_PTM = "ipTM"
# COL_REF2015 = "REF2015"
COL_dG_SASA_ratio = "dG_SASA_ratio"
COL_PRODIGY_DG_INTERNAL = "prodigy_dg_internal"
COL_DG_ROSETTA = "dG_rosetta"


def log(m: str) -> None:
    print(m, flush=True)


def parse_debug_file(dirpath: Path) -> Optional[Dict[str, Any]]:
    """
    Lit ranking_debug.json et renvoie un dict:
    {
        'model': <nom_du_top_model>,
        'iptm+ptm': float,
        'iptm': float
    }
    ou None si absent / invalide.
    """
    f = dirpath / "ranking_debug.json"
    debug_data: Dict[str, Any] = {}
    if not f.exists():
        return None

    data = json.loads(f.read_text())
    order = data.get("order") or []
    if not order:
        return None

    best_model = order[0]
    iptm_ptm = data.get("iptm+ptm", {}).get(best_model)
    iptm = data.get("iptm", {}).get(best_model)
    debug_data = {"model": best_model, "iptm+ptm": iptm_ptm, "iptm": iptm}
    return debug_data


def expected_files(dirpath: Path, top_model: str) -> Tuple[Path, Path]:
    pdb = dirpath / f"unrelaxed_{top_model}.pdb"
    pkl = dirpath / f"result_{top_model}.pkl"
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
    if stem.endswith((".pdb", ".cif")):
        base = stem[:-4]
    else:
        base = str(pdb.with_suffix(""))

    summary = Path(f"{base}_{int(PAE_CUTOFF):02d}_{int(DIST_CUTOFF):02d}.txt")
    if not summary.exists():
        summary = workdir / f"{Path(base).name}_{int(PAE_CUTOFF):02d}_{int(DIST_CUTOFF):02d}.txt"
    if not summary.exists():
        log(f"[ERROR] Résumé ipSAE introuvable: {summary}")
        return None, None

    ipsae_vals: List[float] = []
    pdq2_vals: List[float] = []
    txt = summary.read_text().splitlines()

    # Première tentative : lignes " max "
    for line in txt:
        if " max " in line:
            parts = line.split()
            try:
                if parts[4] == "max":
                    ipsae_vals.append(float(parts[5]))
                    pdq2_vals.append(float(parts[12]))
            except Exception:
                pass

    # Fallback : lignes " asym "
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


def run_prodigy(pdb_path: Path) -> Tuple[Optional[float], Optional[float]]:
    """
    Exécute le CLI local sur LE FICHIER PDB sélectionné.
    On passe --showall pour garantir l'impression du Kd.

    Retourne (kd, dg_internal) :
      - kd (float) en M ou None,
      - dg_internal (float) ou None.
    """
    if not pdb_path.exists():
        log(f"[WARN] PDB introuvable pour PRODIGY: {pdb_path}")
        return None, None


    try:
        from modules.parsers import parse_structure
        from modules.prodigy import Prodigy

        models, _, _ = parse_structure(str(pdb_path))
        if models:
            model = models[0]
            prodigy = Prodigy(model=model, name=pdb_path.stem, temp=25.0)
            prodigy.predict()
            # prodigy.ba_val : predicted binding affinity (kcal/mol as per modules.prodigy)
            # prodigy.kd_val : predicted dissociation constant (M)
            kd = float(prodigy.kd_val) if getattr(prodigy, "kd_val", None) is not None else None
            ba = float(prodigy.ba_val) if getattr(prodigy, "ba_val", None) is not None else None
            return kd, ba
    except Exception as e:
        # Import or processing may fail if dependencies missing (freesasa, biopython...),
        # on tombe alors sur le fallback qui appelle la CLI et parse la sortie.
        log(f"[INFO] run_prodigy direct failed ({e}), falling back to CLI parsing")

    cmd = PRODIGY_CMD + [str(pdb_path)]
    try:
        p = subprocess.run(
            cmd, check=True, capture_output=True, text=True,
            cwd=str(pdb_path.parent)
        )
        out = p.stdout + "\n" + p.stderr
    except subprocess.CalledProcessError as e:
        log(f"[WARN] PRODIGY a échoué sur {pdb_path.name}: {e}\nSTDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}")
        return None, None

    # Parsers tolérants pour Kd
    patterns_kd = [
        r"dissociation\s+constant.*?:\s*([0-9.+\-eE]+)",
        r"\bkd[_\s]*val\b[^0-9eE+\-]*([0-9.+\-eE]+)",
        r"\bKd\s*\(M\)\s*[:=]\s*([0-9.+\-eE]+)",
    ]
    kd_val = None
    for pat in patterns_kd:
        m = re.search(pat, out, flags=re.IGNORECASE)
        if m:
            try:
                kd_val = float(m.group(1))
                break
            except Exception:
                pass

    # Parsers tolérants pour le ΔG interne (reporté par PRODIGY)
    dg_internal = None
    # Cherche des motifs comme 'Delta G', 'ΔG', 'Binding energy', 'Binding affinity' avec un nombre et optionnellement une unité
    dg_patterns = [
        r"(?:ΔG|Delta\s*G|Binding\s+energy|Binding\s+affinity)[^0-9\-\+\.eE\n\r]{0,40}([+\-]?[0-9]*\.?[0-9]+(?:[eE][+\-]?\d+)?)(?:\s*(kcal|kJ|kj|kcal/mol|kJ/mol))?",
        r"(?:binding\s+affinity)[^0-9\-\+\.eE\n\r]{0,40}([+\-]?[0-9]*\.?[0-9]+(?:[eE][+\-]?\d+)?)",
    ]
    for pat in dg_patterns:
        m = re.search(pat, out, flags=re.IGNORECASE)
        if m:
            try:
                dg_internal = float(m.group(1))
                break
            except Exception:
                pass

    if kd_val is None:
        head = "\n".join(out.splitlines()[:40])
        log(
            f"[INFO] Kd non détecté dans la sortie PRODIGY pour {pdb_path.name}. "
            f"Extrait de sortie:\n{head}\n--- fin extrait ---"
        )

    return kd_val, dg_internal

def init_pyrosetta_quiet():
    """
    Initialise PyRosetta en silence (sans bannière ni logs C++).
    Retourne True si succès, False sinon.
    """
    if not HAS_PYROSETTA:
        log("[INFO] PyRosetta non disponible, calculs Rosetta désactivés")
        return False
    
    try:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            pyrosetta.init("-mute all")
        log("[OK] PyRosetta initialisé")
        return True
    except Exception as e:
        log(f"[ERROR] Échec initialisation PyRosetta: {e}")
        return False

def load_table_auto(path: Path, sheet):
    """
    Charge automatiquement un fichier Excel (.xlsx/.xls/.ods)
    ou un CSV. Si l'extension dit 'xlsx' mais que le contenu
    est un CSV (très courant quand un fichier a été écrasé),
    tente un fallback en CSV.
    """
    suffix = path.suffix.lower()

    # --- Excel cas nominal ---
    if suffix in [".xlsx", ".xls", ".ods"]:
        try:
            return pd.read_excel(path, sheet_name=sheet)
        except Exception:
            # Fallback CSV si le fichier n'est pas un vrai Excel
            try:
                return pd.read_csv(path)
            except Exception as e:
                raise RuntimeError(
                    f"Impossible de lire '{path}'. Ni Excel ni CSV ne fonctionne.\n{e}"
                )

    # --- CSV simple ---
    if suffix == ".csv":
        return pd.read_csv(path)

    # --- Autre extension : tente Excel -> CSV ---
    try:
        return pd.read_excel(path, sheet_name=sheet)
    except Exception:
        return pd.read_csv(path)


def save_table_auto(df: pd.DataFrame, path: Path, sheet):
    """
    Sauvegarde automatiquement dans le même format que l'entrée :
      - XLSX / XLS → Excel
      - CSV → CSV
    """
    suffix = path.suffix.lower()

    if suffix in [".xlsx", ".xls"]:
        with pd.ExcelWriter(path, engine="openpyxl", mode="w") as w:
            df.to_excel(w, sheet_name=sheet if isinstance(sheet, str) else "Sheet1", index=False)
        return

    if suffix == ".csv":
        df.to_csv(path, index=False)
        return

    # fallback par défaut → Excel
    with pd.ExcelWriter(path, engine="openpyxl", mode="w") as w:
        df.to_excel(w, sheet_name=sheet if isinstance(sheet, str) else "Sheet1", index=False)


def update_excel_write_safe(
    excel: Path,
    sheet,
    id_col: str,
    updates: Dict[str, Dict[str, Optional[float]]]
) -> None:
    df = load_table_auto(excel, sheet)

    # Déterminer la colonne d'ID :
    # - si id_col (ex: 'jobs') existe, on l'utilise
    # - sinon, on prend la première colonne du tableau
    if id_col in df.columns:
        effective_id_col = id_col
    else:
        effective_id_col = df.columns[0]
        log(
            f"[INFO] Colonne '{id_col}' absente, "
            f"utilisation de la première colonne '{effective_id_col}' comme ID."
        )

    # Colonnes métriques: crée si absentes
    for c in (
        COL_IPSAE,
        COL_PDOCKQ2,
        COL_PRODIGY,
        COL_PRODIGY_DG_INTERNAL,
        COL_IPTM_PTM,
        COL_PTM,
        COL_dG_SASA_ratio,
    ):
        if c not in df.columns:
            df[c] = pd.NA


    ok, miss = 0, 0
    for inter_id, vals in updates.items():
        mask = df[effective_id_col].astype(str) == str(inter_id)
        if mask.any():
            idx = df.index[mask][0]
            for k, v in vals.items():
                if v is not None:
                    df.at[idx, k] = v
            ok += 1
        else:
            miss += 1
            log(f"[WARN] ID '{inter_id}' absent du tableau")

    # Sauvegarde sûr (tmp puis remplace)
    tmp_path = excel.with_suffix(excel.suffix + ".tmp")
    try:
        save_table_auto(df, tmp_path, sheet)
        os.replace(tmp_path, excel)
        log(f"[OK] Tableau mis à jour ({ok} lignes, {miss} manquantes)")
    except PermissionError:
        log(
            "[ERROR] Fichier Excel/CSV verrouillé (ouvert ?). "
            f"Fermes '{excel.name}' puis relance. Export temporaire: {tmp_path}"
        )


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "EXPLORAE: enrich AlphaFold-Multimer interaction models with "
            "ipSAE, pDockQ2, PRODIGY Kd, ipTM/pTM and dG/dSASA, then update a master Excel file."
        )
    )
    parser.add_argument(
    "excel",
    help="Path to the master table file (Excel or CSV, e.g. results.xlsx or results.csv)"
    )

    parser.add_argument(
        "interactions_root",
        help=(
            "Root directory containing one folder per interaction "
            "(each with ranking_debug.json, unrelaxed_*.pdb, result_*.pkl, etc.)"
        )
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
    parser.add_argument(
        "--af3",
        action="store_true",
        help="Active le mode AF3 (mmCIF + confidences.json + summary_confidences.json)"
    )
    return parser.parse_args()


def dico_to_csv(dataset: Dict[str, Dict[str, Any]], filepath: Path) -> None:
    df = pd.DataFrame.from_dict(dataset, orient="index")
    df.index.name = "id"
    df.to_csv(filepath, encoding="utf-8")
    log(f"[OK] Export CSV écrit: {filepath}")



def detect_chains_from_pdb(pdb_file: str) -> Tuple[str, str]:
    """
    Détecte automatiquement les deux premières chaînes d'un fichier PDB.
    
    Arguments:
        pdb_file (str): chemin vers le fichier PDB
    
    Retourne:
        Tuple[str, str]: (chain1, chain2) - les deux premières chaînes trouvées
    """
    chains = set()
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    chain_id = line[21].strip()
                    if chain_id:
                        chains.add(chain_id)
                    if len(chains) >= 2:
                        break
        
        if len(chains) < 2:
            log(f"[WARN] Moins de 2 chaînes trouvées dans {pdb_file}, utilisation de B et C par défaut")
            return "B", "C"
        
        chain_list = sorted(list(chains))
        return chain_list[0], chain_list[1]
    
    except Exception as e:
        log(f"[WARN] Erreur détection chaînes dans {pdb_file}: {e}, utilisation de B et C par défaut")
        return "B", "C"


def analyze_interface(pdb_file: str, chains_partner1: str, chains_partner2: str) -> Dict[str, Optional[float]]:
    """
    Analyse l'interface entre deux ensembles de chaînes dans un complexe.

    Arguments :
        pdb_file (str) : chemin vers le fichier PDB
        chains_partner1 (str) : ex. "A" ou "AHL"
        chains_partner2 (str) : ex. "B" ou "CD"

    Retourne :
        dict : métriques calculées (dSASA, dG, dG_norm)
    """
    if not HAS_PYROSETTA:
        return {
            "interface": f"{chains_partner1}_{chains_partner2}",
            "dG": None,
            "dSASA": None,
            "dG_norm": None,
        }
    
    try:
        pose = pose_from_pdb(pdb_file)

        # Format Rosetta attendu : "A_B" ou "AHL_CD"
        interface_str = f"{chains_partner1}_{chains_partner2}"

        iam = InterfaceAnalyzerMover(interface_str, False)
        iam.apply(pose)

        scores = pose.scores

        dG = scores.get("dG_cross", None)
        dSASA = scores.get("dSASA_int", None)

        dG_norm = dG / dSASA if (dG is not None and dSASA not in (0, None)) else None

        return {
            "interface": interface_str,
            "dG": dG,
            "dSASA": dSASA,
            "dG_norm": dG_norm,
        }
    except Exception as e:
        log(f"[ERROR] Analyse Rosetta échouée pour {pdb_file}: {e}")
        return {
            "interface": f"{chains_partner1}_{chains_partner2}",
            "dG": None,
            "dSASA": None,
            "dG_norm": None,
        }


def main():
    args = parse_args()

    excel = Path(args.excel).resolve()
    root = Path(args.interactions_root).resolve()

    if not root.exists():
        sys.exit(f"[ERROR] Dossier interactions introuvable: {root}")

    global PAE_CUTOFF, DIST_CUTOFF, ID_COL, SHEET
    PAE_CUTOFF = args.pae
    DIST_CUTOFF = args.dist
    ID_COL = args.id_col

    sheet_arg = args.sheet
    if isinstance(sheet_arg, str) and sheet_arg.isdigit():
        SHEET = int(sheet_arg)
    else:
        SHEET = sheet_arg

    if getattr(args, "af3", False):
        from af3_handle import run_af3
        log("[INFO] Bascule vers le workflow AF3")
        run_af3(excel, root, SHEET, ID_COL, PAE_CUTOFF, DIST_CUTOFF)
        return

    # Collect directories
    dirs = sorted([p for p in root.iterdir() if p.is_dir()])
    total = len(dirs)
    
    # Initialisation PyRosetta (optionnel)
    has_rosetta = init_pyrosetta_quiet()
    
    log("="*60)
    log("Format: AF2")
    log(f"PyRosetta: {'initialisé' if has_rosetta else 'non disponible'}")
    log(f"Structures à traiter: {total}")
    log("="*60)

    updates: Dict[str, Dict[str, Optional[float]]] = {}

    for idx, inter_dir in enumerate(dirs, 1):
        inter_id = inter_dir.name
        # Progress bar
        progress = int(50 * idx / total)
        bar = "█" * progress + "░" * (50 - progress)
        print(f"\r[{bar}] {idx}/{total}", end="", flush=True)

        debug_data = parse_debug_file(inter_dir)
        if not debug_data:
            continue

        top = debug_data["model"]
        if not top:
            continue

        pdb, pkl = expected_files(inter_dir, top)
        if not pdb.exists() or not pkl.exists():
            continue

        ipsae_val, pdq2_val = run_ipsae(pkl, pdb, inter_dir)
        kd_val, dg_internal = run_prodigy(pdb)

        # Détection automatique des chaînes pour Rosetta
        chain1, chain2 = detect_chains_from_pdb(str(pdb))
        dG_SASA_results = analyze_interface(str(pdb), chain1, chain2)

        # Récupération du dG Rosetta (dG_cross) déjà calculé par analyze_interface
        dG_rosetta = dG_SASA_results.get("dG") if isinstance(dG_SASA_results, dict) else None

        updates[inter_id] = {
            COL_IPSAE: ipsae_val,
            COL_PDOCKQ2: pdq2_val,
            COL_PRODIGY: kd_val,
            COL_PRODIGY_DG_INTERNAL: dg_internal,
            COL_DG_ROSETTA: dG_rosetta,
            COL_IPTM_PTM: debug_data.get("iptm+ptm"),
            COL_PTM: debug_data.get("iptm"),
            COL_dG_SASA_ratio: dG_SASA_results["dG_norm"],
        }
    
    print()  # New line after progress bar
    if not excel.exists():
        log("="*60)
        # Si l'utilisateur a donné un .csv → on convertit en Excel avec même nom de base
        if excel.suffix.lower() == ".csv":
            excel_output = excel.with_suffix(".xlsx")
            df = pd.DataFrame.from_dict(updates, orient="index")
            df.index.name = ID_COL
            with pd.ExcelWriter(excel_output, engine="openpyxl", mode="w") as w:
                df.to_excel(w, sheet_name=str(SHEET) if isinstance(SHEET, int) else SHEET, index=True)
            log(f"✓ Fichier Excel créé: {excel_output.resolve()}")
            log("="*60)
            return

        # Sinon → on crée un Excel
        df = pd.DataFrame.from_dict(updates, orient="index")
        df.index.name = ID_COL
        with pd.ExcelWriter(excel, engine="openpyxl", mode="w") as w:
            df.to_excel(w, sheet_name=str(SHEET) if isinstance(SHEET, int) else SHEET, index=True)
        log(f"✓ Fichier Excel créé: {excel.resolve()}")
        log("="*60)
        return

    # 2) Si le fichier existe déjà → on le met à jour
    # Si c'est un CSV, on le convertit en Excel
    log("="*60)
    if excel.suffix.lower() == ".csv":
        excel_output = excel.with_suffix(".xlsx")
        update_excel_write_safe(excel, SHEET, ID_COL, updates)
        # Après mise à jour, sauvegarder en Excel
        df = load_table_auto(excel, SHEET)
        with pd.ExcelWriter(excel_output, engine="openpyxl", mode="w") as w:
            df.to_excel(w, sheet_name=str(SHEET) if isinstance(SHEET, int) else SHEET, index=False)
        log(f"✓ Fichier Excel mis à jour: {excel_output.resolve()}")
    else:
        update_excel_write_safe(excel, SHEET, ID_COL, updates)
        log(f"✓ Fichier Excel mis à jour: {excel.resolve()}")
    log("="*60)

if __name__ == "__main__":
    main()
