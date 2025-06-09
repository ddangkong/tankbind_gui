#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TankBind GUI with Pipeline Integration and PDBQT Support
- GUI: PyQt5, runs in tankbind_gpu kernel
- Pipeline: Runs in pymol kernel via subprocess
- PDBQT: Added support for run_pdbqt.py with QProgressBar for AutoDock Vina docking
- Affinity Prediction: Unchanged from original code
- Prerequisite: 'pip install ninja' in the 'tankbind_gpu' conda environment
"""

import site, sys
site.ENABLE_USER_SITE = False
user_site = site.getusersitepackages()
sys.path[:] = [p for p in sys.path if not p.startswith(user_site)]

import types, importlib.metadata as md
dist = md.distribution("torchmetrics")
tm = types.ModuleType("torchmetrics")
tm.__version__ = dist.version
tm.__file__ = str(dist.locate_file("torchmetrics/__init__.py"))
sys.modules["torchmetrics"] = tm
for sub in [
    "torchmetrics.functional",
    "torchmetrics.functional.audio",
    "torchmetrics.functional.audio.pit",
    "torchmetrics.functional.image",
    "torchmetrics.functional.image._deprecated",
    "torchmetrics.utilities",
    "torchmetrics.utilities.imports",
    "torchmetrics.utilities.data",
    "torchmetrics.utilities.checks",
    "torchmetrics.aggregation",
    "torchmetrics.metric",
    "torchmetrics.classification",
]:
    stub = types.ModuleType(sub)
    if sub.endswith("audio.pit"):
        stub.permutation_invariant_training = lambda *a, **k: None
        stub.pit_permutate = lambda *a, **k: None
    sys.modules[sub] = stub

import os, csv, gc, shutil, logging, subprocess, hashlib, time
from typing import Optional
from tqdm import tqdm
import torch
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from Bio.PDB import PDBParser
from torch.utils import cpp_extension
import torchdrug.utils.torch as td_torch
from PyQt5.QtWidgets import (QMessageBox, QApplication, QMainWindow, QWidget, 
                             QVBoxLayout, QHBoxLayout, QLabel, QPushButton, 
                             QLineEdit, QTableWidget, QTableWidgetItem, 
                             QScrollArea, QProgressBar, QFileDialog, QHeaderView)
from PyQt5.QtCore import Qt, QSize, QTimer, QThread, pyqtSignal, QProcess
from PyQt5.QtGui import QIcon
import nest_asyncio

# XDG_RUNTIME_DIR 설정
if 'XDG_RUNTIME_DIR' not in os.environ:
    os.environ['XDG_RUNTIME_DIR'] = '/tmp/runtime-jovyan'
    os.makedirs(os.environ['XDG_RUNTIME_DIR'], exist_ok=True)
    os.chmod(os.environ['XDG_RUNTIME_DIR'], 0o700)

# Jupyter 이벤트 루프 통합
nest_asyncio.apply()
os.environ["QT_LOGGING_RULES"] = "qt5ct.debug=false"

# --- 경로 설정 ---
PDB_ROOT = "/home/jovyan/Chembl(효능데이터베이스)/TestPDBs/"
CENTER_PATH = "/home/jovyan/Chembl(효능데이터베이스)/center_coords.txt"
TANKBIND_DIR = "/home/jovyan/TankBind"
OUT_DIR = "/home/jovyan/TankBindTest/interface/pdb_peptide_run"
CSV_PATH = os.path.join(OUT_DIR, "affinity_all.csv")
FINAL_OUT_FILE = os.path.join(OUT_DIR, "high_affinity.csv")
TMP_BASE = os.path.join(OUT_DIR, "tmp_ds")
PIPELINE_SCRIPT = "/home/jovyan/TankBindTest/interface/pipeline.py"
PYMOL_PYTHON = "/home/jovyan/miniconda/envs/pymol/bin/python"
# PDBQT 관련 경로
RUN_PDBQT_SCRIPT = "/home/jovyan/TankBindTest/interface/run_pdbqt.py"
OBABEL_EXECUTABLE = "/home/jovyan/miniconda/envs/pymol/bin/obabel"
AUTODOCK_VINA_EXECUTABLE = "/home/jovyan/miniconda/envs/autodock/bin/vina"
PDBQT_OUTPUT_BASE_DIR = os.path.join(OUT_DIR, "vina_pymol_results")
HISTORY_FILE_PATH = os.path.join(OUT_DIR, "search_history.txt")

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(TMP_BASE, exist_ok=True)
os.makedirs(PDBQT_OUTPUT_BASE_DIR, exist_ok=True)

# 로깅 설정
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(os.path.join(OUT_DIR, "tankbind.log")),
        logging.StreamHandler(sys.stdout)
    ]
)
log = logging.getLogger(__name__)

# --- GPU / torchdrug 확장 로드 ---
device = "cuda" if torch.cuda.is_available() else "cpu"
log.info(f"device = {device}")

td_torch.sparse_coo_tensor = lambda idx, val, size: torch.sparse_coo_tensor(idx, val, size)
os.environ["TORCH_EXTENSIONS_DIR"] = TANKBIND_DIR
cpp_extension.VERBOSE = False
build_dir = os.path.join(TANKBIND_DIR, "torch_ext")
os.makedirs(build_dir, exist_ok=True)
td_torch.torch_ext = td_torch.load_extension(
    "torch_ext",
    [os.path.join(sys.prefix, "lib", f"python{sys.version_info.major}.{sys.version_info.minor}",
                  "site-packages", "torchdrug", "utils", "extension", "torch_ext.cpp")],
    build_directory=build_dir,
    with_cuda=True, verbose=False
)

# --- TankBind 패키지 import & 모델 로드 ---
sys.path.insert(0, os.path.join(TANKBIND_DIR, "tankbind"))
from feature_utils import get_protein_feature, extract_torchdrug_feature_from_mol, three_to_one
from data import TankBind_prediction
from model import get_model
from torch_geometric.loader import DataLoader

try:
    model = get_model(0, log, device)
    model.load_state_dict(torch.load(
        os.path.join(TANKBIND_DIR, "saved_models", "self_dock.pt"),
        map_location=device
    ))
    model.eval()
except Exception as e:
    log.error(f"Model loading failed: {str(e)}")
    raise

EDGE_PAD_DIM = next(
    m.edge_linear.in_features
    for m in model.modules() if hasattr(m, "edge_linear")
)
log.info(f"edge_feature dim = {EDGE_PAD_DIM}")

# --- pocket center 로드 ---
try:
    center_df = pd.read_csv(
        CENTER_PATH, header=None, sep=": ", engine="python",
        names=["file", "coords"]
    )
    center_df["pdb"] = center_df["file"].str.replace(".pdb", "", regex=False)
    CENTER_MAP = {
        row.pdb: row.coords.strip("[]").replace(" ", "")
        for _, row in center_df.iterrows()
    }
    log.debug(f"Loaded {len(CENTER_MAP)} PDBs in CENTER_MAP: {list(CENTER_MAP.keys())[:5]}...")
except Exception as e:
    log.error(f"Failed to load center coords: {str(e)}")
    raise

# --- 결과 파일 초기화 ---
done = set()
if os.path.exists(CSV_PATH):
    with open(CSV_PATH) as f:
        reader = csv.reader(f)
        header = next(reader, None)
        for row in reader:
            if len(row) >= 2:
                done.add(f"{row[0]}\t{row[1]}")
            elif row:
                log.warning(f"Malformed row in {CSV_PATH}: {row}")
else:
    with open(CSV_PATH, "w", newline="", encoding='utf-8') as f:
        csv.writer(f).writerow(["pdb", "peptide", "affinity"])

# --- 단백질 특징 캐시 ---
prot_cache = {}
def get_prot_feat(pdbid: str):
    if pdbid not in prot_cache:
        path = os.path.join(PDB_ROOT, f"{pdbid}.pdb")
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure(pdbid, path)
            residues = list(structure.get_residues())
            filtered = [r for r in residues if three_to_one.get(r.get_resname())]
            if len(filtered) < len(residues):
                log.warning(f"[{pdbid}] 비표준 잔기 {len(residues)-len(filtered)}개 제외 (Original: {len(residues)}, Filtered: {len(filtered)})")
            if not filtered:
                log.error(f"[{pdbid}] No standard residues found after filtering. Cannot get protein feature.")
                prot_cache[pdbid] = None
                return None
            prot_cache[pdbid] = get_protein_feature(filtered)
            log.info(f"Computed protein features for {pdbid}")
        except Exception as e:
            log.error(f"Failed to process PDB {pdbid}: {str(e)}")
            prot_cache[pdbid] = None
            return None
    return prot_cache[pdbid]

# --- peptide → Mol 함수 ---
def peptide_to_mol(seq: str, max_try: int = 3) -> Optional[Chem.Mol]:
    seq = seq.strip().upper()
    AA_MAP = {"O": "P"}
    seq2 = "".join(AA_MAP.get(a, a) for a in seq)
    log.info(f"[Mol] Processed sequence: {seq2}")
    if any(ch not in "ACDEFGHIKLMNPQRSTVWY" for ch in seq2):
        log.error(f"[Mol] Invalid sequence: {seq2}")
        return None
    
    mol = Chem.MolFromFASTA(seq2)
    if mol is None:
        log.error(f"[Mol] MolFromFASTA failed: {seq2}")
        return None
    log.info(f"[Mol] Molecule created with {mol.GetNumAtoms()} atoms")

    params = AllChem.ETKDGv3()
    params.useBasicKnowledge = True
    params.useExpTorsionAnglePrefs = True
    params.maxIterations = 2000
    params.randomSeed = 42
    
    for i in range(max_try):
        log.info(f"[Mol] Attempting 3D embedding {i+1}/{max_try} for {seq2}")
        if AllChem.EmbedMolecule(mol, params) == 0:
            log.info(f"[Mol] EmbedMolecule successful on attempt {i+1}")
            try:
                if AllChem.UFFOptimizeMolecule(mol, maxIters=1000) == 0:
                    log.info(f"[Mol] 3D embedding and optimization successful for {seq2}")
                    return mol
                else:
                    log.warning(f"[Mol] UFFOptimizeMolecule failed on attempt {i+1}")
            except RuntimeError as e:
                log.warning(f"[Mol] UFFOptimizeMolecule runtime error on attempt {i+1} for {seq2}: {e}")
        else:
            log.warning(f"[Mol] EmbedMolecule failed on attempt {i+1}")
        
        mol.RemoveAllConformers()

    log.error(f"[Mol] 3D embedding failed after {max_try} tries: {seq2}")
    return None

# --- 친화도 예측 함수 ---
def run_affinity_for_seq(seq: str):
    results = []
    try:
        pdb_files_all = os.listdir(PDB_ROOT)
        pdb_files = [f for f in pdb_files_all if f.endswith(".pdb")]
        if not pdb_files:
            log.error(f"No PDB files found in {PDB_ROOT}. Cannot run affinity predictions.")
            return results
    except FileNotFoundError:
        log.error(f"PDB_ROOT directory not found: {PDB_ROOT}")
        return results
    except Exception as e:
        log.error(f"Error listing PDB files in {PDB_ROOT}: {e}")
        return results

    total_pdbs = len(pdb_files)
    log.debug(f"Processing {total_pdbs} PDB files from {PDB_ROOT} for sequence: {seq}")

    mol = peptide_to_mol(seq)
    
    if mol is None:
        log.error(f"Molecule generation failed for sequence {seq}. Marking all PDBs as MolFail for this sequence.")
        for pdb_file in pdb_files:
            pdbid = pdb_file[:-4]
            key = f"{pdbid}\t{seq}"
            if key not in done:
                with open(CSV_PATH, "a", newline="", encoding='utf-8') as f:
                    csv.writer(f).writerow([pdbid, seq, "MolFail"])
                done.add(key)
            results.append((pdbid, seq, "MolFail"))
        return results

    try:
        raw_ligand_features = extract_torchdrug_feature_from_mol(mol, has_LAS_mask=True)
        coords, node_f, edge_l, edge_f_orig, pair_dis = raw_ligand_features[:5]
        
        d = edge_f_orig.size(1)
        if d < EDGE_PAD_DIM:
            edge_f = torch.cat([edge_f_orig, edge_f_orig.new_zeros(edge_f_orig.size(0), EDGE_PAD_DIM - d)], 1)
        elif d > EDGE_PAD_DIM:
            edge_f = edge_f_orig[:, :EDGE_PAD_DIM]
        else:
            edge_f = edge_f_orig
        
        ligand_feature_tuple = (coords, node_f, edge_l, edge_f, pair_dis)
        log.info(f"[Mol] Ligand features extracted successfully for sequence {seq}")
    except Exception as e:
        log.error(f"Ligand feature extraction failed for sequence {seq}: {str(e)}. Marking as FeatExtFail for all PDBs.")
        for pdb_file in pdb_files:
            pdbid = pdb_file[:-4]
            key = f"{pdbid}\t{seq}"
            if key not in done:
                with open(CSV_PATH, "a", newline="", encoding='utf-8') as f:
                    csv.writer(f).writerow([pdbid, seq, "FeatExtFail"])
                done.add(key)
            results.append((pdbid, seq, "FeatExtFail"))
        return results

    for i, pdb_file in enumerate(tqdm(pdb_files, desc=f"PDB screening for {seq}", leave=False)):
        pdbid = pdb_file[:-4]
        if pdbid not in CENTER_MAP:
            log.warning(f"[{pdbid}] center coordinates missing — skipping this PDB.")
            continue

        key = f"{pdbid}\t{seq}"
        if key in done:
            log.info(f"[{pdbid}|{seq}] Already processed (found in 'done' set), skipping affinity prediction.")
            try:
                df_temp = pd.read_csv(CSV_PATH, dtype={'pdb': str, 'peptide': str, 'affinity': str})
                existing_row = df_temp[(df_temp['pdb'] == pdbid) & (df_temp['peptide'] == seq)]
                if not existing_row.empty:
                    results.append((pdbid, seq, str(existing_row.iloc[0]['affinity'])))
            except Exception as read_err:
                log.error(f"Error reading {CSV_PATH} for existing 'done' entry {key}: {read_err}")
            continue

        tmp_dir = os.path.join(TMP_BASE, f"tmp_{pdbid}_{seq.replace('/','_')}")
        shutil.rmtree(tmp_dir, ignore_errors=True)
        os.makedirs(tmp_dir, exist_ok=True)
        
        score = "PredErr"
        try:
            protein_feature = get_prot_feat(pdbid)
            if protein_feature is None:
                log.error(f"[{pdbid}|{seq}] Failed to get protein features. Setting score to ProtErr.")
                score = "ProtErr"
            else:
                ds = TankBind_prediction(
                    root=tmp_dir,
                    data=pd.DataFrame([[pdbid, "lig", "center", CENTER_MAP[pdbid]]],
                                      columns=["protein_name", "compound_name", "pocket_name", "pocket_com"]),
                    protein_dict={pdbid: protein_feature},
                    compound_dict={"lig": ligand_feature_tuple}
                )
                if not ds or len(ds) == 0:
                    log.error(f"[{pdbid}|{seq}] Dataset creation resulted in an empty dataset.")
                    score = "DataErr"
                else:
                    data_loader = DataLoader(
                        ds, batch_size=1,
                        follow_batch=["x", "y", "compound_pair"],
                        shuffle=False
                    )
                    batch = next(iter(data_loader), None)
                    if batch is None:
                        log.error(f"[{pdbid}|{seq}] DataLoader returned no data.")
                        score = "DataLoadErr"
                    else:
                        batch = batch.to(device)
                        with torch.no_grad():
                            _, aff = model(batch)
                        score = f"{aff.item():.4f}"
                        log.info(f"[{pdbid}|{seq}] Prediction successful: Affinity = {score}")

        except Exception as e:
            log.error(f"[{pdbid}|{seq}] Prediction error during TankBind model call: {str(e)}", exc_info=True)
            score = "PredErr" 

        shutil.rmtree(tmp_dir, ignore_errors=True)
        with open(CSV_PATH, "a", newline="", encoding='utf-8') as f:
            csv.writer(f).writerow([pdbid, seq, score])
        done.add(key)
        results.append((pdbid, seq, score))
        
        if (i + 1) % 20 == 0:
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
            gc.collect()

    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    gc.collect()

    log.debug(f"Predicted {len(results)} results for sequence {seq}: {[(r[0], r[2]) for r in results][:5]}...")
    return results

# --- 파이프라인 실행 함수 ---
def run_pipeline_script():
    if not os.path.exists(PIPELINE_SCRIPT):
        log.error(f"Pipeline script not found: {PIPELINE_SCRIPT}")
        return False
    if not os.path.exists(PYMOL_PYTHON):
        log.error(f"PyMOL Python not found: {PYMOL_PYTHON}")
        return False
    try:
        log.info(f"Attempting to run pipeline: {PYMOL_PYTHON} {PIPELINE_SCRIPT}")
        process = subprocess.Popen(
            [PYMOL_PYTHON, PIPELINE_SCRIPT],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True, 
            bufsize=1, 
            universal_newlines=True 
        )
        log.info(f"Started pipeline script: {PIPELINE_SCRIPT} (PID: {process.pid})")

        def log_output(pipe, log_func):
            try:
                with pipe:
                    for line in iter(pipe.readline, ''):
                        log_func(f"[Pipeline PID:{process.pid}] {line.strip()}")
            except Exception as e:
                log.error(f"Error in logging pipeline output: {e}")

        from threading import Thread
        stdout_thread = Thread(target=log_output, args=(process.stdout, log.info), daemon=True)
        stderr_thread = Thread(target=log_output, args=(process.stderr, log.error), daemon=True)
        stdout_thread.start()
        stderr_thread.start()
        
        return True
    except Exception as e:
        log.error(f"Failed to run pipeline: {str(e)}")
        return False

# --- PDB ID와 단백질 이름 매핑 ---
PDB_TO_PROTEIN = {
    "5F19": "Cyclooxygenase-2 (COX-2)", "5F1A": "Cyclooxygenase-2 (COX-2)", "5IKQ": "Cyclooxygenase-2 (COX-2)",
    "5IKR": "Cyclooxygenase-2 (COX-2)", "5IKT": "Cyclooxygenase-2 (COX-2)", "5IKV": "Cyclooxygenase-2 (COX-2)",
    "5KIR": "Cyclooxygenase-2 (COX-2)", "2AZ5": "Tumor Necrosis Factor-α (TNF-α)", "7JRA": "Tumor Necrosis Factor-α (TNF-α)",
    "7KP9": "Tumor Necrosis Factor-α (TNF-α)", "7KPA": "Tumor Necrosis Factor-α (TNF-α)",
    "6O0P": "Apoptosis regulator BCL-2 (BCL-2)", "8GCP": "Prostaglandin E₂ receptor EP4 (PTGER4)",
    "8GDB": "Prostaglandin E₂ receptor EP4 (PTGER4)", "5YWY": "Prostaglandin E₂ receptor EP4 (PTGER4)",
    "8GD9": "Prostaglandin E₂ receptor EP4 (PTGER4)", "7D7M": "Prostaglandin E₂ receptor EP4 (PTGER4)",
    "8GDA": "Prostaglandin E₂ receptor EP4 (PTGER4)", "8GCM": "Prostaglandin E₂ receptor EP4 (PTGER4)",
    "5YHL": "Prostaglandin E₂ receptor EP4 (PTGER4)", "2FST": "p38α Mitogen-Activated Protein Kinase (MAPK14)",
    "3LFF": "p38α Mitogen-Activated Protein Kinase (MAPK14)", "3OEF": "p38α Mitogen-Activated Protein Kinase (MAPK14)",
    "3ZS5": "p38α Mitogen-Activated Protein Kinase (MAPK14)", "4EHV": "p38α Mitogen-Activated Protein Kinase (MAPK14)",
    "5WJJ": "p38α Mitogen-Activated Protein Kinase (MAPK14)", "6SFI": "p38α Mitogen-Activated Protein Kinase (MAPK14)",
    "4GEO": "p38α Mitogen-Activated Protein Kinase (MAPK14)", "6QYX": "p38α Mitogen-Activated Protein Kinase (MAPK14)",
    "2FSL": "p38α Mitogen-Activated Protein Kinase (MAPK14)", "2QD9": "p38α Mitogen-Activated Protein Kinase (MAPK14)",
    "3FMK": "p38α Mitogen-Activated Protein Kinase (MAPK14)", "4A7U": "Superoxide dismutase [Cu-Zn] (SOD1)",
    "4A7V": "Superoxide dismutase [Cu-Zn] (SOD1)", "6Z3V": "Superoxide dismutase [Cu-Zn] (SOD1)",
    "5O3Y": "Superoxide dismutase [Cu-Zn] (SOD1)", "7T8G": "Superoxide dismutase [Cu-Zn] (SOD1)",
    "2WZ6": "Superoxide dismutase [Cu-Zn] (SOD1)", "8CCX": "Superoxide dismutase [Cu-Zn] (SOD1)",
    "1DGF": "Catalase (CAT)", "1DGB": "Catalase (CAT)", "1DGG": "Catalase (CAT)", "1DGH": "Catalase (CAT)",
    "2OBI": "Glutathione peroxidase 4 (GPX4)", "5H5Q": "Glutathione peroxidase 4 (GPX4)",
    "5H5R": "Glutathione peroxidase 4 (GPX4)", "5H5S": "Glutathione peroxidase 4 (GPX4)",
    "6HKQ": "Glutathione peroxidase 4 (GPX4)", "7P8W": "Catalase (CAT)",
    "7U4J": "Glutathione peroxidase 4 (GPX4)", "7U4K": "Glutathione peroxidase 4 (GPX4)",
    "7VD9": "Catalase (CAT)", "8EL9": "Catalase (CAT)", "8HID": "Catalase (CAT)"
}

# --- PyQt5 UI ---
class RecentTransactionsWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("TankBind Screening")
        self.setMinimumSize(1280, 800)
        self.setStyleSheet("background-color:#FFFFFF;")

        self.pipeline_running = False
        self.pipeline_start_time = None
        self.pdbqt_running = False  # PDBQT 프로세스 상태 추적

        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        # --- 사이드바 ---
        side = QWidget()
        side.setFixedWidth(300)
        side.setStyleSheet("""
            QWidget { background-color: #080325; color: #E2E8F0; }
            QLabel { color: #FFFFFF; font-size: 16px; font-weight: bold; }
            QPushButton {
                background-color: #191A38; color: #FFFFFF; text-align: left;
                padding-left: 12px; border: none; border-radius: 4px; font-size: 14px;
            }
            QPushButton:hover { background-color: #2A2C4D; }
        """)
        sl = QVBoxLayout(side)
        sl.setContentsMargins(18, 24, 18, 32)
        sl.setSpacing(24)
        logo = QLabel("🚀 Nongshim")
        logo.setStyleSheet("color:#FFF;font-size:26px;font-weight:bold; margin-bottom: 12px;")
        sl.addWidget(logo)
        history_label = QLabel("Search History")
        history_label.setStyleSheet("margin-top:10px; font-size:16px; font-weight:bold;")
        sl.addWidget(history_label, alignment=Qt.AlignLeft)
        
        self.hist_cont = QWidget()
        self.hist_lyt = QVBoxLayout(self.hist_cont)
        self.hist_lyt.setContentsMargins(0, 0, 0, 0)
        self.hist_lyt.setSpacing(4)
        self.hist_lyt.setAlignment(Qt.AlignTop)
        
        hist_scroll = QScrollArea()
        hist_scroll.setWidgetResizable(True)
        hist_scroll.setWidget(self.hist_cont)
        hist_scroll.setStyleSheet("""
            QScrollArea { border: none; background-color: transparent; }
            QScrollBar:vertical { border: none; background: #1A1C3A; width: 8px; margin: 0px 0px 0px 0px; }
            QScrollBar::handle:vertical { background: #4A507A; min-height: 25px; border-radius: 4px; }
            QScrollBar::handle:vertical:hover { background: #5A608A; }
            QScrollBar::handle:vertical:pressed { background: #3A406A; }
            QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0px; background: none; }
            QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical { background: none; }
            QScrollBar:horizontal { border: none; background: #1A1C3A; height: 8px; margin: 0px 0px 0px 0px; }
            QScrollBar::handle:horizontal { background: #4A507A; min-width: 25px; border-radius: 4px; }
            QScrollBar::handle:horizontal:hover { background: #5A608A; }
            QScrollBar::handle:horizontal:pressed { background: #3A406A; }
            QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal { width: 0px; background: none; }
            QScrollBar::add-page:horizontal, QScrollBar::sub-page:horizontal { background: none; }
        """)
        sl.addWidget(hist_scroll, 1) 

        export_buttons_container = QWidget()
        export_buttons_layout = QHBoxLayout(export_buttons_container)
        export_buttons_layout.setContentsMargins(0,5,0,0)
        export_buttons_layout.setSpacing(10)

        self.export_csv_button = QPushButton("Export CSV")
        self.export_csv_button.setFixedHeight(35)
        self.export_csv_button.setStyleSheet("""
            background-color: #1E88E5; color: #FFFFFF; text-align: center;
            padding: 5px 10px; border: none; border-radius: 4px; font-size: 13px; font-weight: bold;
        """)
        self.export_csv_button.clicked.connect(self.export_to_csv)
        export_buttons_layout.addWidget(self.export_csv_button)

        self.export_excel_button = QPushButton("Export Excel")
        self.export_excel_button.setFixedHeight(35)
        self.export_excel_button.setStyleSheet("""
            background-color: #43A047; color: #FFFFFF; text-align: center;
            padding: 5px 10px; border: none; border-radius: 4px; font-size: 13px; font-weight: bold;
        """)
        self.export_excel_button.clicked.connect(self.export_to_excel)
        export_buttons_layout.addWidget(self.export_excel_button)
        sl.addWidget(export_buttons_container)
        root.addWidget(side)
        
        # --- 메인영역 ---
        mc = QWidget()
        mc.setStyleSheet("background-color: #FFFFFF;")
        root.addWidget(mc, 1)
        ml = QVBoxLayout(mc)
        ml.setContentsMargins(30, 30, 30, 30)
        ml.setSpacing(20)
        
        hr = QWidget()
        hrl = QHBoxLayout(hr)
        hrl.setContentsMargins(0, 0, 0, 0)
        hrl.setSpacing(15)
        title = QLabel("Protein Affinity System")
        title.setStyleSheet("font-size:26px;font-weight:bold;color:#080325;")
        hrl.addWidget(title)
        hrl.addStretch()
        self.input = QLineEdit()
        self.input.setPlaceholderText("Peptide 서열 입력: 최대 15개 아미노산 (예: AAG, GPP)")
        self.input.setFixedHeight(40)
        self.input.setStyleSheet(
            "border:1px solid #CBD5E0;border-radius:20px;padding:5px 15px;font-size:15px;"
        )
        self.input.returnPressed.connect(self.on_search)
        hrl.addWidget(self.input)
        sb = QPushButton("Search")
        sb.setFixedHeight(40)
        sb.setStyleSheet(
            "background-color:#080325;color:white;border-radius:20px;padding: 0 20px;font-size:15px;font-weight:bold; text-align: center;"
        )
        sb.clicked.connect(self.on_search)
        hrl.addWidget(sb)
        ml.addWidget(hr)

        self.progress_bar = QProgressBar()
        self.progress_bar.setFixedHeight(12)
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setValue(0)
        self.progress_bar.setVisible(False)
        self.progress_bar.setStyleSheet("""
            QProgressBar {
                border: 1px solid #B0B0B0; border-radius: 6px; text-align: center;
                background-color: #E8E8E8; color: #333333; font-size: 10px;
            }
            QProgressBar::chunk { background-color: #080325; border-radius: 5px; }
        """)
        ml.addWidget(self.progress_bar)

        self.table = QTableWidget(0, 14)
        header_labels = [
            "PDB ID", "Protein", "Peptide", "Affinity", "MolWt", "LogP", "HBD", "HBA",
            "RotBonds", "TPSA", "Pass_Lipinski", "Pass_Veber", "Pass_All", "Reason"
        ]
        self.table.setSortingEnabled(True)
        self.table.setAlternatingRowColors(True)
        self.table.verticalHeader().setDefaultSectionSize(50)
        self.table.verticalHeader().setVisible(False)
        self.table.setStyleSheet("""
            QTableWidget {
                border: 1px solid #D0D0D0; selection-background-color: #AED6F1;
                selection-color: #000000; background-color: #FFFFFF;
                alternate-background-color: #F9F9F9; gridline-color: #E0E0E0;
            }
            QTableWidget QHeaderView::section {
                background-color: #F0F4F8; color: #080325; padding: 8px;
                font-size: 13px; font-weight: bold; border-style: none;
                border-bottom: 1px solid #D0D0D0;
            }
            QTableWidget QHeaderView::section:horizontal { border-right: 1px solid #D0D0D0; }
            QTableWidget QHeaderView::section:horizontal:last { border-right: none; }
            QTableWidget QScrollBar:vertical { border: none; background: #F0F0F0; width: 12px; margin: 0px 0px 0px 0px; }
            QTableWidget QScrollBar::handle:vertical { background: #BCC6CC; min-height: 30px; border-radius: 6px; }
            QTableWidget QScrollBar::handle:vertical:hover { background: #A8B3B9; }
            QTableWidget QScrollBar::handle:vertical:pressed { background: #93A0A6; }
            QTableWidget QScrollBar::add-line:vertical, QTableWidget QScrollBar::sub-line:vertical { height: 0px; background: none; }
            QTableWidget QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical { background: none; }
            QTableWidget QScrollBar:horizontal { border: none; background: #F0F0F0; height: 12px; margin: 0px 0px 0px 0px; }
            QTableWidget QScrollBar::handle:horizontal { background: #BCC6CC; min-width: 30px; border-radius: 6px; }
            QTableWidget QScrollBar::handle:horizontal:hover { background: #A8B3B9; }
            QTableWidget QScrollBar::handle:horizontal:pressed { background: #93A0A6; }
            QTableWidget QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal { width: 0px; background: none; }
            QTableWidget QScrollBar::add-page:horizontal, QScrollBar::sub-page:horizontal { background: none; }
        """)
        TOOLTIP_DESCRIPTIONS = {
            "PDB ID": "단백질 구조 데이터베이스(PDB)의 4~6자리 고유코드", "Protein": "PDB ID에 해당하는 단백질(효소) 이름",
            "Peptide": "스크리닝 대상 펩타이드 서열 (입력된 아미노산 문자 코드)",
            "Affinity": "TankBind 예측 결합 친화도 점수\n └ 7 이상: 강결합, 5~7: 중간, 5 미만: 약결합",
            "MolWt": "분자량 (g/mol)\n └ 500 이하가 경구 약물에 유리, 800 이상이면 흡수성 저하",
            "LogP": "소수성 지표 cLogP\n └ 1~3 이상적, 5↑ 용해도↓/독성↑, 0↓ 막투과성↓",
            "HBD": "수소 결합 공여자 수 (H-donor)\n └ 0~5 권장, 많을수록 막투과성·생체이용률↓",
            "HBA": "수소 결합 수용자 수 (H-acceptor)\n └ 0~10 권장, 과다 시 투과성↓",
            "RotBonds": "회전 가능한 결합 수\n └ 10 이하가 흡수성·결합 엔트로피에 유리",
            "TPSA": "위상 극성 표면적 (Å²)\n └ 140 이하: 경구 흡수 가능, 90 이하: CNS 투과 가능",
            "Pass_Lipinski": "경구투여 후보 물질의 흡수 가능성 여부\n └ MolWt≤500, LogP≤5, HBD≤5, HBA≤10 일때 흡수율↑",
            "Pass_Veber": "Veber 규칙 충족 여부\n └ TPSA≤140 & RotBonds≤10 → 장내 흡수 예측",
            "Pass_All": "모든 약물 유사성 필터 통합 결과\n └ True면 초기 ADME 기준 통과",
            "Reason": "필터 탈락 시 원인 상세 (예: MolWt>500)"
        }
        
        for i, label_text in enumerate(header_labels):
            tooltip_text = TOOLTIP_DESCRIPTIONS.get(label_text, label_text)
            header_item = QTableWidgetItem(label_text)
            header_item.setToolTip(tooltip_text)
            self.table.setHorizontalHeaderItem(i, header_item)
            self.table.horizontalHeader().setSectionResizeMode(i, QHeaderView.Stretch)
            log.debug(f"Set tooltip for header {label_text}: {tooltip_text if tooltip_text else 'No tooltip text'}")

        ml.addWidget(self.table, 1)

        self.last_modified = 0
        self.last_hash = ""
        self.current_sequence = None
        self.processed_pipeline_pdbs_for_sequence = set()

        self.timer = QTimer(self)
        self.timer.timeout.connect(self.check_pipeline_results)
        self.timer.start(5000)

        self.load_search_history()
        self.table.cellClicked.connect(self.on_pdb_cell_clicked)

        # QProcess 초기화
        self.pdbqt_process = QProcess(self)
        self.pdbqt_process.started.connect(self.on_pdbqt_started)
        self.pdbqt_process.finished.connect(self.on_pdbqt_finished)
        self.pdbqt_process.errorOccurred.connect(self.on_pdbqt_error)

    def on_pdb_cell_clicked(self, row, column):
        if column != 0:
            return
        pdb_id = self.table.item(row, column).text() if self.table.item(row, column) else None
        if not pdb_id or not self.current_sequence:
            self.show_error("No PDB ID or peptide sequence selected.")
            return
        if self.pdbqt_running or self.pipeline_running:
            self.show_error("Another process (PDBQT or pipeline) is running. Please wait.")
            return

        protein_pdb_path = os.path.join(PDB_ROOT, f"{pdb_id}.pdb")
        paths = {
            "PDBQT Script": RUN_PDBQT_SCRIPT,
            "Open Babel": OBABEL_EXECUTABLE,
            "Vina": AUTODOCK_VINA_EXECUTABLE,
            "PyMOL Python": PYMOL_PYTHON,
            "Protein PDB": protein_pdb_path
        }
        for name, path in paths.items():
            if not os.path.exists(path):
                self.show_error(f"{name} not found at:\n{path}")
                return
        if pdb_id not in CENTER_MAP:
            self.show_error(f"Center coordinates for {pdb_id} not found.")
            return
        try:
            coords = [float(c.strip()) for c in CENTER_MAP[pdb_id].strip("[] \n").split(',')]
            if len(coords) != 3:
                raise ValueError()
        except:
            self.show_error(f"Invalid center coordinates for {pdb_id}: {CENTER_MAP[pdb_id]}")
            return

        # GPU 메모리 정리
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
            gc.collect()

        # QProcess로 PDBQT 실행
        cmd = [
            paths["PDBQT Script"],
            "--peptide_seq", self.current_sequence,
            "--pdb_id", pdb_id,
            "--protein_pdb_file", paths["Protein PDB"],
            "--base_out_dir", PDBQT_OUTPUT_BASE_DIR,
            "--obabel_executable", paths["Open Babel"],
            "--vina_executable", paths["Vina"],
            "--center_x", str(coords[0]),
            "--center_y", str(coords[1]),
            "--center_z", str(coords[2])
        ]
        log.info(f"Starting Vina/PyMOL: {PYMOL_PYTHON} {' '.join(cmd)}")
        self.pdbqt_process.setProgram(PYMOL_PYTHON)
        self.pdbqt_process.setArguments(cmd)
        self.pdbqt_process.setWorkingDirectory(os.path.dirname(RUN_PDBQT_SCRIPT))
        self.pdbqt_process.start()

    def on_pdbqt_started(self):
        self.pdbqt_running = True
        self.progress_bar.setVisible(True)
        self.progress_bar.setMaximum(0)  # 불확정 상태
        self.progress_bar.setFormat("Running Vina/PyMOL...")
        log.info("PDBQT process started.")
        QApplication.processEvents()

    def on_pdbqt_finished(self, exit_code, exit_status):
        self.pdbqt_running = False
        self.progress_bar.setVisible(False)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setFormat("%p%")
        if exit_code == 0:
            QMessageBox.information(self, "Success", "Vina docking and PyMOL session generation completed.")
            log.info("PDBQT process completed successfully.")
        else:
            self.show_error(f"PDBQT process failed with exit code {exit_code}.")
        QApplication.processEvents()

    def on_pdbqt_error(self, error):
        self.pdbqt_running = False
        self.progress_bar.setVisible(False)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setFormat("%p%")
        self.show_error(f"PDBQT process error: {error}")
        log.error(f"PDBQT process error: {error}")
        QApplication.processEvents()

    def _create_history_button(self, seq_text):
        btn = QPushButton(seq_text)
        btn.setFixedHeight(30)
        btn.setStyleSheet(
            "background:transparent;color:#E2E8F0;font-size:13px;"
            "text-align:left;padding-left:12px;border:none;"
            "border-bottom: 1px solid #2A2C4D;"
        ) 
        btn.clicked.connect(lambda _, s=seq_text: self.fill_and_search(s, from_history=True))
        return btn

    def load_search_history(self):
        if os.path.exists(HISTORY_FILE_PATH):
            try:
                with open(HISTORY_FILE_PATH, "r", encoding="utf-8") as f:
                    history_sequences = [line.strip() for line in f if line.strip()]
                
                for i in range(self.hist_lyt.count()):
                    item = self.hist_lyt.takeAt(0)
                    if item and item.widget():
                        item.widget().deleteLater()

                for seq in reversed(history_sequences): 
                    btn = self._create_history_button(seq)
                    self.hist_lyt.insertWidget(0, btn) 
                log.info(f"Loaded {len(history_sequences)} sequences from search history: {HISTORY_FILE_PATH}")
            except Exception as e:
                log.error(f"Failed to load search history: {e}")

    def save_search_history(self):
        history_sequences = []
        for i in range(self.hist_lyt.count()):
            widget = self.hist_lyt.itemAt(i).widget()
            if isinstance(widget, QPushButton):
                history_sequences.append(widget.text())
        try:
            with open(HISTORY_FILE_PATH, "w", encoding="utf-8") as f:
                for seq in history_sequences: 
                    f.write(seq + "\n")
            log.info(f"Saved {len(history_sequences)} sequences to search history: {HISTORY_FILE_PATH}")
        except Exception as e:
            log.error(f"Failed to save search history: {e}")

    def closeEvent(self, event):
        self.save_search_history()
        super().closeEvent(event)

    def on_search(self):
        seq = self.input.text().strip().upper()
        if not seq:
            QMessageBox.warning(self, "Input Error", "펩타이드 서열을 입력하세요. 최대 15개 아미노산")
            return
        
        valid_chars = set("ACDEFGHIKLMNPQRSTVWYO")
        if not all(char in valid_chars for char in seq):
            QMessageBox.warning(self, "Input Error", "유효한 아미노산 코드를 포함한 펩타이드 서열을 입력하세요 (ACDEFGHIKLMNPQRSTVWYO).")
            return

        processed_seq = seq.replace('O', 'P')
        
        for i in range(self.hist_lyt.count()):
            widget = self.hist_lyt.itemAt(i).widget()
            if isinstance(widget, QPushButton) and widget.text() == processed_seq:
                widget.deleteLater() 
                log.info(f"Sequence {processed_seq} found in history, moved to top.")
                break
        
        new_btn = self._create_history_button(processed_seq)
        self.hist_lyt.insertWidget(0, new_btn) 

        self.last_modified = 0 
        self.last_hash = ""    
        self.processed_pipeline_pdbs_for_sequence.clear()

        affinity_fully_done = False
        num_expected_pdbs = 0
        try:
            num_expected_pdbs = len([f for f in os.listdir(PDB_ROOT) if f.endswith(".pdb") and f[:-4] in CENTER_MAP])
        except Exception as e:
            log.error(f"Could not count expected PDBs: {e}")

        if os.path.exists(CSV_PATH) and num_expected_pdbs > 0:
            try:
                df_affinity = pd.read_csv(CSV_PATH, dtype=str, usecols=['pdb', 'peptide', 'affinity'])
                df_seq_affinity = df_affinity[df_affinity['peptide'] == processed_seq]
                
                valid_affinity_count = 0; error_affinity_present = False
                for _, row_aff in df_seq_affinity.iterrows():
                    try: float(row_aff['affinity']); valid_affinity_count +=1
                    except ValueError: error_affinity_present = True; break
                
                if not error_affinity_present and valid_affinity_count >= num_expected_pdbs:
                    affinity_fully_done = True
                    log.info(f"All {num_expected_pdbs} PDB affinity results for {processed_seq} already exist and are valid in {CSV_PATH}.")
            except Exception as e:
                log.error(f"Error checking {CSV_PATH} for existing results: {str(e)}")

        if affinity_fully_done:
            self.fill_and_search(processed_seq, from_history=True) 
            if run_pipeline_script():
                log.info(f"Pipeline (re)started for fully processed affinity of {processed_seq}")
                self.pipeline_running = True; self.pipeline_start_time = time.time()
                self.progress_bar.setVisible(True); self.progress_bar.setMaximum(0) 
                self.progress_bar.setFormat("Running ADMET & Visualization Pipeline...")
            else:
                log.error(f"Failed to (re)start pipeline for {processed_seq}")
            return
        
        self.fill_and_search(processed_seq, from_history=False) 

    def fill_and_search(self, seq, from_history=False):
        self.table.setSortingEnabled(False) 
        self.table.clearContents(); self.table.setRowCount(0)
        
        self.current_sequence = seq
        self.input.setText(seq) 
        self.processed_pipeline_pdbs_for_sequence.clear()

        log.info(f"{'Loading results from history' if from_history else 'Starting new affinity prediction'} for sequence: {seq}")

        if from_history:
            self.progress_bar.setVisible(True); self.progress_bar.setMaximum(0) 
            self.progress_bar.setFormat("Loading existing results..."); QApplication.processEvents()

            results_to_display = []
            if os.path.exists(CSV_PATH):
                try:
                    df_affinity = pd.read_csv(CSV_PATH, dtype=str)
                    df_affinity_seq = df_affinity[df_affinity['peptide'] == seq]
                    for _, r in df_affinity_seq.iterrows():
                        protein_name = PDB_TO_PROTEIN.get(r['pdb'], "")
                        results_to_display.append((r['pdb'], protein_name, r['peptide'], r['affinity'], "", "", "", "", "", "", "", "", "", ""))
                    log.debug(f"Loaded {len(df_affinity_seq)} affinity results from {CSV_PATH} for {seq}")
                except Exception as e: log.error(f"Error reading {CSV_PATH} for history: {str(e)}")
            
            self.display_results(results_to_display)
            self.check_pipeline_results(force_check=True) 

            QApplication.processEvents() 
            if not self.pipeline_running: 
                self.progress_bar.setVisible(False); self.progress_bar.setMaximum(100)
                self.progress_bar.setFormat("%p%")
        else: 
            self.progress_bar.setVisible(True); self.progress_bar.setMaximum(100); self.progress_bar.setValue(0)
            self.progress_bar.setFormat("Predicting Affinity: %p%")
            self.affinity_thread = AffinityThread(seq)
            self.affinity_thread.started.connect(lambda: self.progress_bar.setValue(10))
            self.affinity_thread.finished.connect(lambda: self.progress_bar.setValue(100))
            self.affinity_thread.result.connect(self.display_affinity_results)
            self.affinity_thread.error.connect(self.show_error)
            self.affinity_thread.start()
        
        QApplication.processEvents()
        self.table.setSortingEnabled(True)

    def display_affinity_results(self, affinity_results_list):
        self.progress_bar.setValue(100)
        log.debug(f"Displaying {len(affinity_results_list)} new affinity prediction results for {self.current_sequence}")
        
        formatted_results = []
        for r_pdb, r_seq, r_aff_val in affinity_results_list:
            protein_name = PDB_TO_PROTEIN.get(r_pdb, "")
            formatted_results.append((r_pdb, protein_name, r_seq, r_aff_val, "", "", "", "", "", "", "", "", "", ""))
        
        self.display_results(formatted_results)

        if run_pipeline_script():
            log.info(f"Pipeline script started automatically for sequence: {self.current_sequence}")
            self.pipeline_running = True; self.pipeline_start_time = time.time()
            self.progress_bar.setVisible(True); self.progress_bar.setMaximum(0) 
            self.progress_bar.setFormat("Running ADMET & Visualization Pipeline...")
        else:
            log.error(f"Failed to start pipeline script automatically for {self.current_sequence}")
            self.show_error("Failed to start ADMET & Visualization pipeline.")
            self.progress_bar.setVisible(False); self.progress_bar.setMaximum(100)
            self.progress_bar.setFormat("%p%")

    def display_results(self, results_list_of_tuples):
        self.table.setSortingEnabled(False)
        current_table_data = {} 
        for i in range(self.table.rowCount()):
            pdb_id_item = self.table.item(i, 0)
            if pdb_id_item:
                pdb_id = pdb_id_item.text()
                row_data = [self.table.item(i, j).text() if self.table.item(i,j) else "" for j in range(self.table.columnCount())]
                current_table_data[pdb_id] = row_data
        
        for r_tuple in results_list_of_tuples:
            if len(r_tuple) != 14: 
                log.warning(f"Skipping malformed result tuple (length {len(r_tuple)} instead of 14): {r_tuple}")
                continue
            pdb_id, protein_name, peptide_val, affinity_val = r_tuple[0], r_tuple[1], r_tuple[2], r_tuple[3]
            if peptide_val != self.current_sequence: continue

            is_pipeline_data = any(val != "" for val in r_tuple[4:])
            if pdb_id in current_table_data:
                if is_pipeline_data: current_table_data[pdb_id] = list(r_tuple)
                else: 
                    current_is_error = str(current_table_data[pdb_id][3]).lower() in ["err", "molfail", "proterr", "dataerr", "prederr", "featextfail", "dataloaderr"]
                    new_is_not_error = str(affinity_val).lower() not in ["err", "molfail", "proterr", "dataerr", "prederr", "featextfail", "dataloaderr"]
                    has_existing_pipeline_data = any(val != "" for val in current_table_data[pdb_id][4:])
                    if (current_is_error and new_is_not_error) or not has_existing_pipeline_data:
                         current_table_data[pdb_id] = list(r_tuple)
            else: current_table_data[pdb_id] = list(r_tuple)
        
        final_display_list = list(current_table_data.values())
        self.table.clearContents(); self.table.setRowCount(len(final_display_list))

        for i, row_data in enumerate(final_display_list):
            for j, value in enumerate(row_data):
                item = QTableWidgetItem(str(value))
                if j == 3:  
                    try:
                        val_float = float(value)
                        if val_float > 7.0: item.setForeground(Qt.darkGreen)
                        elif val_float > 5.0: item.setForeground(Qt.darkYellow) 
                        else: item.setForeground(Qt.black)
                    except ValueError: item.setForeground(Qt.red) 
                elif j >= 10 and j <= 12:  
                    if str(value).lower() == "true": item.setForeground(Qt.darkGreen)
                    elif str(value).lower() == "false": item.setForeground(Qt.red)
                self.table.setItem(i, j, item)
        
        self.table.setSortingEnabled(True)
        if self.table.columnCount() > 3: self.table.sortByColumn(3, Qt.DescendingOrder) 
        QApplication.processEvents()

    def check_pipeline_results(self, force_check=False):
        if not self.current_sequence: return
        if not self.pipeline_running and not force_check: return
        if not os.path.exists(FINAL_OUT_FILE):
            if self.pipeline_running: log.debug(f"Pipeline output file {FINAL_OUT_FILE} not found yet.")
            return

        if not force_check and self.pipeline_start_time and (time.time() - self.pipeline_start_time < 2):
            log.debug("Pipeline recently started, short grace period before stability check.")
            return

        try:
            last_known_size = os.path.getsize(FINAL_OUT_FILE)
            last_known_mod_time = os.path.getmtime(FINAL_OUT_FILE)
            required_stable_duration_secs = 2
            consecutive_stable_checks = 0

            for _ in range(required_stable_duration_secs): 
                QApplication.processEvents(); time.sleep(0.5)
                if not os.path.exists(FINAL_OUT_FILE): 
                    log.warning(f"File {FINAL_OUT_FILE} disappeared during stability check."); return
                current_size = os.path.getsize(FINAL_OUT_FILE); current_mod_time = os.path.getmtime(FINAL_OUT_FILE)
                if current_size == last_known_size and current_mod_time == last_known_mod_time:
                    consecutive_stable_checks += 1
                else: 
                    last_known_size = current_size; last_known_mod_time = current_mod_time
                    consecutive_stable_checks = 0
                    log.debug(f"Pipeline file {FINAL_OUT_FILE} changed. Resetting stability count.")
            
            if consecutive_stable_checks >= required_stable_duration_secs: 
                log.info(f"Pipeline file {FINAL_OUT_FILE} has been stable for {required_stable_duration_secs * 0.5} seconds.")
                with open(FINAL_OUT_FILE, 'rb') as f_hash: current_file_hash = hashlib.md5(f_hash.read()).hexdigest()
                
                if force_check or last_known_mod_time > self.last_modified or current_file_hash != self.last_hash:
                    log.info(f"Processing updated pipeline file {FINAL_OUT_FILE}. Stable ModTime: {last_known_mod_time}, Hash: {current_file_hash}")
                    self.last_modified = last_known_mod_time; self.last_hash = current_file_hash
                    df_pipeline = pd.read_csv(FINAL_OUT_FILE, dtype=str)
                    df_pipeline['affinity'] = pd.to_numeric(df_pipeline['affinity'], errors='coerce').fillna(0.0) 
                    df_pipeline_seq = df_pipeline[df_pipeline['peptide'] == self.current_sequence]
                    
                    pipeline_results_for_display = []; newly_processed_count = 0
                    for _, row in df_pipeline_seq.iterrows():
                        result_key = (row.get('pdb',""), row.get('peptide',""), f"{row.get('affinity',0.0):.4f}") 
                        protein_name = PDB_TO_PROTEIN.get(row.get('pdb', ""), "")
                        pipeline_row_tuple = (
                            row.get('pdb', ""), protein_name, row.get('peptide', ""), str(row.get('affinity', "")),
                            row.get('MolWt', ""), row.get('LogP', ""), row.get('HBD', ""), row.get('HBA', ""),
                            row.get('RotBonds', ""), row.get('TPSA', ""),
                            str(row.get('Pass_Lipinski', "")).lower(), str(row.get('Pass_Veber', "")).lower(),
                            str(row.get('Pass_All', "")).lower(), row.get('Reason', "")
                        )
                        pipeline_results_for_display.append(pipeline_row_tuple)
                        if result_key not in self.processed_pipeline_pdbs_for_sequence:
                            self.processed_pipeline_pdbs_for_sequence.add(result_key); newly_processed_count +=1

                    if newly_processed_count > 0 or force_check: 
                        log.info(f"Found {len(pipeline_results_for_display)} pipeline results for {self.current_sequence}. {newly_processed_count} new/forced.")
                        self.display_results(pipeline_results_for_display)
                    else:
                        log.info(f"Pipeline file {FINAL_OUT_FILE} for {self.current_sequence} processed, no new PDB entries this cycle.")

                    if self.pipeline_running:
                        self.pipeline_running = False; self.progress_bar.setVisible(False)
                        self.progress_bar.setMaximum(100); self.progress_bar.setFormat("%p%")
                else: 
                    log.info(f"Pipeline file {FINAL_OUT_FILE} stable but unchanged. Pipeline finished for {self.current_sequence}.")
                    if self.pipeline_running: 
                        self.pipeline_running = False; self.progress_bar.setVisible(False)
                        self.progress_bar.setMaximum(100); self.progress_bar.setFormat("%p%")
            else: 
                log.debug(f"Pipeline file {FINAL_OUT_FILE} not stable.")
                if self.pipeline_running and not self.progress_bar.isVisible(): 
                    self.progress_bar.setVisible(True); self.progress_bar.setMaximum(0) 
                    self.progress_bar.setFormat("Running ADMET & Visualization Pipeline...")
        except FileNotFoundError:
            if self.pipeline_running and not self.progress_bar.isVisible():
                self.progress_bar.setVisible(True); self.progress_bar.setMaximum(0)
                self.progress_bar.setFormat("Running ADMET & Visualization Pipeline...")
        except Exception as e:
            log.error(f"Error during pipeline result check for {FINAL_OUT_FILE}: {str(e)}", exc_info=True)
            if self.pipeline_running: 
                self.progress_bar.setVisible(True); self.progress_bar.setMaximum(0) 
                self.progress_bar.setFormat("Pipeline Error. Checking...")
        finally: QApplication.processEvents()

    def show_error(self, error_msg):
        self.progress_bar.setVisible(False); self.progress_bar.setMaximum(100)
        self.progress_bar.setFormat("%p%"); self.pipeline_running = False; self.pdbqt_running = False
        QMessageBox.critical(self, "Error", f"An error occurred: {error_msg}")
        log.error(f"GUI Error Displayed: {error_msg}")

    def export_to_csv(self):
        if self.table.rowCount() == 0:
            QMessageBox.information(self, "Export CSV", "내보낼 데이터가 없습니다."); return
        default_filename = f"{self.current_sequence or 'export'}_results.csv"
        path, _ = QFileDialog.getSaveFileName(self, "CSV 파일 저장", default_filename, "CSV Files (*.csv)")
        if path:
            try:
                with open(path, 'w', newline='', encoding='utf-8-sig') as csvfile:
                    writer = csv.writer(csvfile)
                    headers = [self.table.horizontalHeaderItem(i).text() for i in range(self.table.columnCount())]
                    writer.writerow(headers)
                    for row in range(self.table.rowCount()):
                        row_data = [self.table.item(row, col).text() if self.table.item(row, col) else ""
                                    for col in range(self.table.columnCount())]
                        writer.writerow(row_data)
                QMessageBox.information(self, "Export CSV", f"데이터가 성공적으로 저장되었습니다:\n{path}")
                log.info(f"Exported table data to CSV: {path}")
            except Exception as e:
                QMessageBox.critical(self, "내보내기 오류", f"CSV 파일 저장 실패: {e}"); log.error(f"Error exporting to CSV: {e}")

    def export_to_excel(self):
        if self.table.rowCount() == 0:
            QMessageBox.information(self, "Export Excel", "내보낼 데이터가 없습니다."); return
        default_filename = f"{self.current_sequence or 'export'}_results.xlsx"
        path, _ = QFileDialog.getSaveFileName(self, "Excel 파일 저장", default_filename, "Excel Files (*.xlsx)")
        if path:
            try:
                headers = [self.table.horizontalHeaderItem(i).text() for i in range(self.table.columnCount())]
                data = []
                for row in range(self.table.rowCount()):
                    row_data = [self.table.item(row, col).text() if self.table.item(row, col) else ""
                                for col in range(self.table.columnCount())]
                    data.append(row_data)
                df = pd.DataFrame(data, columns=headers)
                df.to_excel(path, index=False, engine='openpyxl') 
                QMessageBox.information(self, "Export Excel", f"데이터가 성공적으로 저장되었습니다:\n{path}")
                log.info(f"Exported table data to Excel: {path}")
            except ImportError:
                QMessageBox.critical(self, "내보내기 오류", "Excel 파일 저장 실패: `openpyxl` 라이브러리가 필요합니다. (예: pip install openpyxl)")
                log.error("Error exporting to Excel: openpyxl not found.")
            except Exception as e:
                QMessageBox.critical(self, "내보내기 오류", f"Excel 파일 저장 실패: {e}"); log.error(f"Error exporting to Excel: {e}")

# --- 친화도 예측 스레드 ---
class AffinityThread(QThread):
    result = pyqtSignal(list) 
    error = pyqtSignal(str)
    def __init__(self, seq): super().__init__(); self.seq = seq.strip().upper()
    def run(self):
        try:
            log.info(f"AffinityThread started for sequence: {self.seq}")
            results_list = run_affinity_for_seq(self.seq)
            self.result.emit(results_list)
            log.info(f"AffinityThread completed for sequence: {self.seq}")
        except Exception as e:
            log.error(f"AffinityThread error for sequence {self.seq}: {str(e)}", exc_info=True)
            self.error.emit(str(e))

if __name__ == "__main__":
    try:
        log_file_path = os.path.join(OUT_DIR, "tankbind.log")
        if os.path.exists(log_file_path) and os.path.getsize(log_file_path) > 10 * 1024 * 1024: 
            os.remove(log_file_path); logging.info("Old log file was large and has been removed.")
    except Exception as e: logging.warning(f"Could not check/clear old log file: {e}")

    try:
        from IPython import get_ipython
        ipython = get_ipython()
        if ipython and hasattr(ipython, 'run_line_magic'):
             ipython.run_line_magic('gui', 'qt5')
    except ImportError: pass 
    except Exception as e: log.warning(f"Could not set IPython GUI magic: {e}")

    app = QApplication.instance() 
    if not app: app = QApplication(sys.argv)
    
    app.setStyleSheet("""
        QToolTip {
            background-color: #FFFFFF;
            color: #333740;
            border: 1px solid #1E1F22;
            border-radius: 4px;
            padding: 4px 6px;
            font-size: 12px;
            opacity: 230;
        }
    """)

    main_window = RecentTransactionsWindow()
    main_window.show()
    
    if not (hasattr(sys, 'ps1') and sys.ps1) and not os.getenv('PYCHARM_HOSTED'):
        sys.exit(app.exec_())
    else:
        log.info("Running in an interactive environment or PyCharm, app.exec_() might be managed externally.")
