#!/usr/bin/env python3
"""
AF3 handling workflow:
- consumes folders with model.cif + confidences.json + summary_confidences.json
- computes ipSAE / pDockQ2 via ipsae.py (AF3 mode)
- computes PRODIGY Kd via modules.prodigy (fallback to cli.py)
- updates Excel/CSV with same columns as AF2 pipeline
"""
from __future__ import annotations
import json, sys, subprocess, re, os
from pathlib import Path
from typing import Optional, Dict, Any, Tuple, List
import pandas as pd

# Column names aligned with explorae.py
COL_IPSAE   = "ipsae"
COL_PDOCKQ2 = "pdockq2"
COL_PRODIGY = "prodigy_kd"
COL_PRODIGY_DG_INTERNAL = "prodigy_dg_internal"
COL_IPTM_PTM = "ipTM+pTM"
COL_PTM = "ipTM"
COL_dG_SASA_ratio = "dG_SASA_ratio"
COL_DG_ROSETTA = "dG_rosetta"

IPSae_SCRIPT = Path(__file__).resolve().parent / "ipsae.py"
PRODIGY_CMD = [sys.executable, str((Path(__file__).parent / "cli.py").resolve()), "--showall"]


def log(msg: str) -> None:
    print(msg, flush=True)


def load_table_auto(path: Path, sheet):
    suffix = path.suffix.lower()
    if suffix in [".xlsx", ".xls", ".ods"]:
        try:
            return pd.read_excel(path, sheet_name=sheet)
        except Exception:
            return pd.read_csv(path)
    if suffix == ".csv":
        return pd.read_csv(path)
    try:
        return pd.read_excel(path, sheet_name=sheet)
    except Exception:
        return pd.read_csv(path)


def save_table_auto(df: pd.DataFrame, path: Path, sheet):
    suffix = path.suffix.lower()
    if suffix in [".xlsx", ".xls"]:
        with pd.ExcelWriter(path, engine="openpyxl", mode="w") as w:
            df.to_excel(w, sheet_name=sheet if isinstance(sheet, str) else "Sheet1", index=False)
        return
    if suffix == ".csv":
        df.to_csv(path, index=False)
        return
    with pd.ExcelWriter(path, engine="openpyxl", mode="w") as w:
        df.to_excel(w, sheet_name=sheet if isinstance(sheet, str) else "Sheet1", index=False)


def update_excel_write_safe(
    excel: Path,
    sheet,
    id_col: str,
    updates: Dict[str, Dict[str, Optional[float]]]
) -> None:
    df = load_table_auto(excel, sheet)
    if id_col in df.columns:
        effective_id_col = id_col
    else:
        effective_id_col = df.columns[0]
        log(f"[INFO] Colonne '{id_col}' absente, utilisation de '{effective_id_col}' comme ID.")

    for c in (
        COL_IPSAE,
        COL_PDOCKQ2,
        COL_PRODIGY,
        COL_PRODIGY_DG_INTERNAL,
        COL_IPTM_PTM,
        COL_PTM,
        COL_dG_SASA_ratio,
        COL_DG_ROSETTA,
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

    tmp_path = excel.with_suffix(excel.suffix + ".tmp")
    try:
        save_table_auto(df, tmp_path, sheet)
        os.replace(tmp_path, excel)
        log(f"[OK] Tableau mis à jour ({ok} lignes, {miss} manquantes)")
    except PermissionError:
        log(
            "[ERROR] Fichier Excel/CSV verrouillé (ouvert ?). "
            f"Ferme '{excel.name}' puis relance. Export temporaire: {tmp_path}"
        )


def parse_summary_confidences(dirpath: Path) -> Optional[Dict[str, Any]]:
    f = dirpath / "summary_confidences.json"
    if not f.exists():
        return None
    try:
        data = json.loads(f.read_text())
        iptm = data.get("iptm")
        ptm = data.get("ptm")
        # iptm+ptm absent côté AF3 : on le reconstruit 0.8*iptm + 0.2*ptm si possible
        iptm_ptm = 0.8 * iptm + 0.2 * ptm
        return {
            "iptm+ptm": iptm_ptm,
            "iptm": iptm,
            "ptm": ptm,
        }
    except Exception as e:
        log(f"[WARN] summary_confidences.json illisible: {e}")
        return None


def find_af3_files(inter_dir: Path) -> Optional[Tuple[Path, Path, Path]]:
    cif = inter_dir / "model.cif"
    pae = inter_dir / "confidences.json"
    summary = inter_dir / "summary_confidences.json"
    if cif.exists() and pae.exists() and summary.exists():
        return cif, pae, summary
    log(f"[WARN] Fichiers AF3 manquants dans {inter_dir.name}")
    return None


def run_ipsae_af3(pae_json: Path, cif: Path, workdir: Path, pae_cutoff: float, dist_cutoff: float) -> Tuple[Optional[float], Optional[float]]:
    if not IPSae_SCRIPT.exists():
        log(f"[ERROR] ipsae.py introuvable: {IPSae_SCRIPT}")
        return None, None
    cmd = [sys.executable, str(IPSae_SCRIPT), str(pae_json), str(cif), str(pae_cutoff), str(dist_cutoff)]
    try:
        subprocess.run(cmd, cwd=str(workdir), check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        log(f"[ERROR] ipsae.py a échoué dans {workdir.name}\nSTDERR: {e.stderr.decode(errors='ignore')}")
        return None, None

    stem = str(cif)
    base = stem[:-4] if stem.endswith(".cif") else str(cif.with_suffix(""))
    summary = Path(f"{base}_{int(pae_cutoff):02d}_{int(dist_cutoff):02d}.txt")
    if not summary.exists():
        summary = workdir / f"{Path(base).name}_{int(pae_cutoff):02d}_{int(dist_cutoff):02d}.txt"
    if not summary.exists():
        log(f"[ERROR] Résumé ipSAE introuvable: {summary}")
        return None, None

    ipsae_vals: List[float] = []
    pdq2_vals: List[float] = []
    for line in summary.read_text().splitlines():
        if " max " in line:
            parts = line.split()
            try:
                if parts[4] == "max":
                    ipsae_vals.append(float(parts[5]))
                    pdq2_vals.append(float(parts[12]))
            except Exception:
                pass
    if not ipsae_vals or not pdq2_vals:
        for line in summary.read_text().splitlines():
            if " asym " in line:
                parts = line.split()
                try:
                    ipsae_vals.append(float(parts[5]))
                    pdq2_vals.append(float(parts[12]))
                except Exception:
                    pass
    ipSAE = max(ipsae_vals) if ipsae_vals else None
    pDockQ2 = max(pdq2_vals) if pdq2_vals else None
    return ipSAE, pDockQ2


def run_prodigy_structure(struct_path: Path) -> Tuple[Optional[float], Optional[float]]:
    if not struct_path.exists():
        log(f"[WARN] Structure introuvable pour PRODIGY: {struct_path}")
        return None, None
    try:
        from modules.parsers import parse_structure
        from modules.prodigy import Prodigy
        models, _, _ = parse_structure(str(struct_path))
        if models:
            model = models[0]
            prodigy = Prodigy(model=model, name=struct_path.stem, temp=25.0)
            prodigy.predict()
            kd = float(getattr(prodigy, "kd_val", None)) if getattr(prodigy, "kd_val", None) is not None else None
            ba = float(getattr(prodigy, "ba_val", None)) if getattr(prodigy, "ba_val", None) is not None else None
            return kd, ba
    except Exception as e:
        log(f"[INFO] run_prodigy direct failed ({e}), fallback CLI")

    cmd = PRODIGY_CMD + [str(struct_path)]
    try:
        p = subprocess.run(cmd, check=True, capture_output=True, text=True, cwd=str(struct_path.parent))
        out = p.stdout + "\n" + p.stderr
    except subprocess.CalledProcessError as e:
        log(f"[WARN] PRODIGY a échoué sur {struct_path.name}: {e}\nSTDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}")
        return None, None

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

    dg_internal = None
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
        log(f"[INFO] Kd non détecté dans PRODIGY pour {struct_path.name}. Extrait:\n{head}\n--- fin extrait ---")

    return kd_val, dg_internal


def detect_chains_from_mmcif(cif_file: Path) -> Tuple[str, str]:
    try:
        from Bio.PDB.MMCIFParser import MMCIFParser
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(cif_file.stem, str(cif_file))
        model = structure[0]
        chains = [c.id for c in model.get_chains()]
        if len(chains) >= 2:
            return chains[0], chains[1]
    except Exception as e:
        log(f"[WARN] Détection chaînes mmCIF échouée ({e}), fallback B/C")
    return "B", "C"


def run_af3(excel: Path, interactions_root: Path, sheet, id_col: str, pae_cutoff: float, dist_cutoff: float) -> None:
    updates: Dict[str, Dict[str, Optional[float]]] = {}
    
    # Collect directories
    dirs = sorted([p for p in interactions_root.iterdir() if p.is_dir()])
    total = len(dirs)
    
    log("="*60)
    log("Format: AF3")
    log("PyRosetta: désactivé (mmCIF non supporté)")
    log(f"Structures à traiter: {total}")
    log("="*60)
    
    for idx, inter_dir in enumerate(dirs, 1):
        inter_id = inter_dir.name
        # Progress bar
        progress = int(50 * idx / total)
        bar = "█" * progress + "░" * (50 - progress)
        print(f"\r[{bar}] {idx}/{total}", end="", flush=True)

        files = find_af3_files(inter_dir)
        if not files:
            continue
        cif, pae_json, _ = files

        conf_data = parse_summary_confidences(inter_dir) or {}
        ipsae_val, pdq2_val = run_ipsae_af3(pae_json, cif, inter_dir, pae_cutoff, dist_cutoff)
        kd_val, dg_internal = run_prodigy_structure(cif)

        chain1, chain2 = detect_chains_from_mmcif(cif)
        dG_SASA_ratio = None
        dG_rosetta = None

        updates[inter_id] = {
            COL_IPSAE: ipsae_val,
            COL_PDOCKQ2: pdq2_val,
            COL_PRODIGY: kd_val,
            COL_PRODIGY_DG_INTERNAL: dg_internal,
            COL_DG_ROSETTA: dG_rosetta,
            COL_IPTM_PTM: conf_data.get("iptm+ptm"),
            COL_PTM: conf_data.get("iptm"),
            COL_dG_SASA_ratio: dG_SASA_ratio,
        }
    
    print()  # New line after progress bar
    log("="*60)
    if not excel.exists():
        df = pd.DataFrame.from_dict(updates, orient="index")
        df.index.name = id_col
        with pd.ExcelWriter(excel, engine="openpyxl", mode="w") as w:
            df.to_excel(w, sheet_name=str(sheet) if isinstance(sheet, int) else sheet, index=True)
        log(f"✓ Fichier Excel créé: {excel.resolve()}")
        log("="*60)
        return

    if excel.suffix.lower() == ".csv":
        excel_output = excel.with_suffix(".xlsx")
        update_excel_write_safe(excel, sheet, id_col, updates)
        df = load_table_auto(excel, sheet)
        with pd.ExcelWriter(excel_output, engine="openpyxl", mode="w") as w:
            df.to_excel(w, sheet_name=str(sheet) if isinstance(sheet, int) else sheet, index=False)
        log(f"✓ Fichier Excel mis à jour: {excel_output.resolve()}")
    else:
        update_excel_write_safe(excel, sheet, id_col, updates)
        log(f"✓ Fichier Excel mis à jour: {excel.resolve()}")
    log("="*60)


__all__ = ["run_af3"]
