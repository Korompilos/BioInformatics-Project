import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

AUX_DIR = Path(os.getenv("BIO_AUX_DIR", ROOT / "auxiliary2025"))
RESULTS_DIR = Path(os.getenv("BIO_RESULTS_DIR", ROOT / "results"))
SRC_DIR = Path(os.getenv("BIO_SRC_DIR", ROOT / "source2025"))

# --- input data paths ---
DATASET_A_PATH = AUX_DIR / "datasetA.txt"
DATASET_B_PATH = AUX_DIR / "datasetB.txt"
DATASET_C_PATH = AUX_DIR / "datasetC.txt"
Q2_ALIGNMENT_PATH = RESULTS_DIR / "q2_alignment.txt"

# --- visualize output paths ---
VISUALIZATION_DIR = Path(os.getenv("BIO_VIZ_DIR", RESULTS_DIR / "visualization"))
VIZ_Q4 = VISUALIZATION_DIR / "Q4"

# -- results paths ---
HMM_PROFILE_BEFORE = RESULTS_DIR / "hmm_profile_before.json"
HMM_PROFILE_AFTER = RESULTS_DIR / "hmm_profile_after.json"
Q4_RESULTS_PATH = RESULTS_DIR / "q4_results.json"


def ensure_dirs():
    AUX_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

