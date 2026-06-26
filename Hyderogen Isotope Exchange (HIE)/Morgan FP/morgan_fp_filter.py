"""
Morgan ECFP fingerprints (RDKit) for HIE — drop redundant near-constant bits.

Reads:
  - RDKit Fragments - Published/Data/Training_Data.csv (NAMES, class)
  - RDKit Fragments - Published/Names_smiles.csv
  - RDKit Fragments - Published/Names_smiles_Predictions.csv

Writes:
  - Training_morgan_filtered.csv, Prediction_morgan_filtered.csv

Requires: rdkit, pandas, numpy
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

SCRIPT_DIR = Path(__file__).resolve().parent
PUBLISHED = SCRIPT_DIR.parent / "RDKit Fragments - Published"
TRAINING_CSV = PUBLISHED / "Data" / "Training_Data.csv"
NAMES_SMILES_CSV = PUBLISHED / "Names_smiles.csv"
NAMES_SMILES_PRED_CSV = PUBLISHED / "Names_smiles_Predictions.csv"

OUTPUT_TRAIN = SCRIPT_DIR / "Training_morgan_filtered.csv"
OUTPUT_PRED = SCRIPT_DIR / "Prediction_morgan_filtered.csv"

MORGAN_RADIUS = 2
MORGAN_NBITS = 2048
DEFAULT_SIMILAR_FRACTION = 0.9


def _validate_smiles_name_columns(df: pd.DataFrame, path: Path) -> None:
    cols = [c.strip() for c in df.columns]
    if cols != ["SMILES", "Name"]:
        raise ValueError(f"Expected columns [SMILES, Name] in {path}, got {list(df.columns)}")


def smiles_to_fp_array(smiles: str) -> np.ndarray:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, MORGAN_RADIUS, nBits=MORGAN_NBITS)
    arr = np.zeros((MORGAN_NBITS,), dtype=np.float64)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def build_smiles_lookup() -> dict[str, str]:
    ns = pd.read_csv(NAMES_SMILES_CSV)
    pr = pd.read_csv(NAMES_SMILES_PRED_CSV)
    _validate_smiles_name_columns(ns, NAMES_SMILES_CSV)
    _validate_smiles_name_columns(pr, NAMES_SMILES_PRED_CSV)
    smap: dict[str, str] = {}
    for df, path in (ns, NAMES_SMILES_CSV), (pr, NAMES_SMILES_PRED_CSV):
        for _, row in df.iterrows():
            name = str(row["Name"])
            smi = str(row["SMILES"])
            if name in smap and smap[name] != smi:
                raise ValueError(f"Conflicting SMILES for Name={name!r} in {path}")
            smap[name] = smi
    return smap


def load_molecule_tables() -> tuple[list[str], list[str], dict[str, object], np.ndarray]:
    train = pd.read_csv(TRAINING_CSV)
    if "NAMES" not in train.columns or "class" not in train.columns:
        raise ValueError(f"Training_Data.csv must have NAMES and class; got {list(train.columns)}")
    pr = pd.read_csv(NAMES_SMILES_PRED_CSV)
    _validate_smiles_name_columns(pr, NAMES_SMILES_PRED_CSV)

    smap = build_smiles_lookup()
    class_map = dict(zip(train["NAMES"].astype(str), train["class"], strict=True))

    train_names = train["NAMES"].astype(str).tolist()
    pred_names = pr["Name"].astype(str).tolist()
    overlap = set(train_names) & set(pred_names)
    if overlap:
        raise ValueError(f"Name(s) appear in both training and predictions: {sorted(overlap)}")

    ordered_names = train_names + pred_names
    missing = [n for n in ordered_names if n not in smap]
    if missing:
        raise ValueError(f"No SMILES for names: {missing}")

    X = np.vstack([smiles_to_fp_array(smap[n]) for n in ordered_names])
    return train_names, pred_names, class_map, X


def column_dominant_fractions(X: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    frac_zero = X.mean(axis=0)
    frac_one = 1.0 - frac_zero
    dominant_frac = np.maximum(frac_zero, frac_one)
    dominant_value = np.where(frac_zero >= frac_one, 0, 1)
    return dominant_frac, dominant_value, frac_zero


def redundant_bit_mask(X_train: np.ndarray, *, similar_fraction: float) -> np.ndarray:
    if not 0.0 < similar_fraction <= 1.0:
        raise ValueError(f"similar_fraction must be in (0, 1], got {similar_fraction}")
    dominant_frac, _, _ = column_dominant_fractions(X_train)
    return dominant_frac >= similar_fraction


def filter_redundant_columns(
    X: np.ndarray,
    *,
    similar_fraction: float,
    fit_mask_on: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    drop_mask = redundant_bit_mask(fit_mask_on, similar_fraction=similar_fraction)
    keep_idx = np.flatnonzero(~drop_mask)
    return X[:, keep_idx], keep_idx


def mfp_column_names(bit_indices: np.ndarray) -> list[str]:
    return [f"MFP_{int(i)}" for i in bit_indices]


def write_feature_csvs(
    *,
    ordered_names: list[str],
    class_map: dict[str, object],
    n_train: int,
    feature_cols: list[str],
    X_features: np.ndarray,
    train_path: Path,
    pred_path: Path,
) -> None:
    out = pd.DataFrame({"Name": ordered_names})
    out[feature_cols] = X_features
    out["class"] = [class_map.get(n, "") for n in ordered_names]
    out.iloc[:n_train].to_csv(train_path, index=False)
    out.iloc[n_train:].to_csv(pred_path, index=False)


def run_filter(
    train_names: list[str],
    pred_names: list[str],
    class_map: dict[str, object],
    X: np.ndarray,
    *,
    similar_fraction: float,
) -> None:
    n_tr = len(train_names)
    X_train = X[:n_tr, :]
    X_f, keep_idx = filter_redundant_columns(
        X,
        similar_fraction=similar_fraction,
        fit_mask_on=X_train,
    )
    n_drop = int(MORGAN_NBITS - len(keep_idx))
    fp_cols = mfp_column_names(keep_idx)

    write_feature_csvs(
        ordered_names=train_names + pred_names,
        class_map=class_map,
        n_train=n_tr,
        feature_cols=fp_cols,
        X_features=X_f,
        train_path=OUTPUT_TRAIN,
        pred_path=OUTPUT_PRED,
    )
    print(
        f"Wrote {OUTPUT_TRAIN} ({n_tr} rows) and {OUTPUT_PRED} ({len(pred_names)} rows); "
        f"kept {len(keep_idx)} / {MORGAN_NBITS} Morgan bits "
        f"(dropped {n_drop} with >= {similar_fraction:.0%} identical values on training set)"
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Morgan fingerprints for HIE (redundant-bit filter).")
    p.add_argument(
        "--similar-fraction",
        type=float,
        default=DEFAULT_SIMILAR_FRACTION,
        metavar="F",
        help=(
            "Drop a bit if the same value (0 or 1) appears in >= F of training molecules "
            f"(default {DEFAULT_SIMILAR_FRACTION})"
        ),
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    train_names, pred_names, class_map, X = load_molecule_tables()
    run_filter(
        train_names,
        pred_names,
        class_map,
        X,
        similar_fraction=args.similar_fraction,
    )


if __name__ == "__main__":
    main()
