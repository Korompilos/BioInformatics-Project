import os
from pathlib import Path

# project root (parent of the folder that contains config.py)
ROOT = Path(os.getenv("BIO_ROOT_DIR", Path(__file__).resolve().parents[1]))

# base results directory (overridable by env vars)
RESULTS_DIR = Path(os.getenv("BIO_RESULTS_DIR", ROOT / "results"))
DATASET_DIR = RESULTS_DIR / "dataset"
JSON_DIR = RESULTS_DIR / "json"
VISUALIZATION_DIR = RESULTS_DIR / "visualization"

# source directory (Q1–Q4)
SRC_DIR = Path(os.getenv("BIO_SRC_DIR", ROOT / "source2025"))

# --- dataset paths ---
DATASET_A_PATH = DATASET_DIR / "datasetA.txt"
DATASET_B_PATH = DATASET_DIR / "datasetB.txt"
DATASET_C_PATH = DATASET_DIR / "datasetC.txt"
Q2_ALIGNMENT_PATH = DATASET_DIR / "q2_alignment.txt"

# --- results paths (JSONs) ---
HMM_PROFILE_BEFORE = JSON_DIR / "hmm_profile_before.json"
HMM_PROFILE_AFTER  = JSON_DIR / "hmm_profile_after.json"
Q4_RESULTS_PATH    = JSON_DIR / "q4_results.json"

# --- visualize output paths ---
VIZ_Q3 = VISUALIZATION_DIR / "Q3"
VIZ_Q4 = VISUALIZATION_DIR / "Q4"

def ensure_dirs():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    DATASET_DIR.mkdir(parents=True, exist_ok=True)
    JSON_DIR.mkdir(parents=True, exist_ok=True)
    VISUALIZATION_DIR.mkdir(parents=True, exist_ok=True)
