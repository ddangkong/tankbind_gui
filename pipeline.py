#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pipeline.py: Runs in pymol kernel. (Simplified Logic Version)

import sys
import os
import re
import pandas as pd
import gc
import logging
try:
    from rdkit.Chem import AllChem, MolFromSequence, Descriptors, Lipinski
    from rdkit.ML.Cluster import Butina
    from rdkit.DataStructs import BulkTanimotoSimilarity
except ImportError as e:
    print(f"Error: RDKit not installed: {e}")
    sys.exit(1)
from tqdm import tqdm

# --- 로깅 설정 (기존과 동일) ---
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("/home/jovyan/TankBindTest/interface/pdb_peptide_run/tankbind_pipeline.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
log = logging.getLogger(__name__)

# --- 설정 (기존과 동일) ---
CSV_PATH = "/home/jovyan/TankBindTest/interface/pdb_peptide_run/affinity_all.csv"
SEQ_COL = "peptide"
PASS_COL = "Pass_All" # 이 변수는 더 이상 필터링에 사용되지 않지만, 다른 곳에서 참조할 수 있어 유지합니다.
STEP1_OUT_FILE = "/home/jovyan/TankBindTest/interface/pdb_peptide_run/step1_filtered.csv"
# 더 이상 사용되지 않는 파일 경로들이지만, 혹시 모를 참조를 위해 남겨둡니다.
STAGE1_OUT_FILE = "/home/jovyan/TankBindTest/interface/pdb_peptide_run/stage1_reps.csv"
STAGE1_DONE_FILE = "/home/jovyan/TankBindTest/interface/pdb_peptide_run/stage1_done_batches.txt"
STAGE2_OUT_FILE = "/home/jovyan/TankBindTest/interface/pdb_peptide_run/step3_clustered.csv"
FINAL_OUT_FILE = "/home/jovyan/TankBindTest/interface/pdb_peptide_run/high_affinity.csv"
CHUNK_SIZE = 500
NBITS = 1024
RADIUS = 2
TANIMOTO_CUTOFF = 0.4

os.makedirs("/home/jovyan/TankBindTest/interface/pdb_peptide_run", exist_ok=True)


# --- 함수 정의 (기존과 동일) ---
# 아래 함수들은 main()에서 직접 호출되지 않더라도, 다른 함수에 의해 사용될 수 있으므로 그대로 둡니다.

def load_data(path, required_cols=None):
    try:
        log.debug(f"Loading CSV: {path}")
        df = pd.read_csv(path, dtype=str)
        if required_cols:
            missing = [col for col in required_cols if col not in df.columns]
            if missing:
                raise ValueError(f"Missing columns: {missing}")
        log.debug(f"Loaded {len(df)} rows with columns: {df.columns.tolist()}")
        # PDB ID 로깅은 데이터가 많을 경우 너무 길어질 수 있으므로, 필요시 주석 해제합니다.
        # log.debug(f"PDB IDs: {df['pdb'].unique().tolist()}")
        return df
    except Exception as e:
        log.error(f"Error loading '{path}': {e}")
        sys.exit(1)

def gen_fps(seqs):
    pat = re.compile(r'^[A-Za-z]+$')
    fps, idxs = [], []
    for i, seq in enumerate(tqdm(seqs, desc="Generating fingerprints", unit="seq")):
        if isinstance(seq, str) and pat.match(seq):
            try:
                m = MolFromSequence(seq)
                if m:
                    fp = AllChem.GetMorganFingerprintAsBitVect(m, RADIUS, nBits=NBITS)
                    fps.append(fp)
                    idxs.append(i)
                del m
            except Exception as e:
                log.warning(f"Failed to process sequence {seq}: {e}")
        if i % 500 == 0:
            gc.collect()
    return fps, idxs

def cluster_and_extract(fps_list, base_idxs, desc="clustering"):
    n = len(fps_list)
    dists = []
    for i in tqdm(range(n), desc=f"Distance calc {desc}", unit="row"):
        sims = BulkTanimotoSimilarity(fps_list[i], fps_list[i+1:])
        dists.extend([1.0 - s for s in sims])
    clusters = Butina.ClusterData(dists, n, 1.0 - TANIMOTO_CUTOFF, isDistData=True)
    reps, reasons = [], []
    for cid, cluster in enumerate(clusters):
        rep_local = cluster[0]
        rep_glob = base_idxs[rep_local]
        reps.append(rep_glob)
        reasons.append(f"Cluster {cid}: size={len(cluster)}, rep idx={rep_glob}")
    return reps, reasons

def step1_physchem_filter(input_path, output_path):
    log.info("Step 1: Physicochemical Properties Calculation")
    df = load_data(input_path, required_cols=[SEQ_COL, 'pdb', 'affinity'])
    
    # 'err', 'FeatExtFail' 등 텍스트가 포함된 행을 제거하고 숫자형 친화도 값만 남깁니다.
    df['affinity_numeric'] = pd.to_numeric(df['affinity'], errors='coerce')
    df_clean = df.dropna(subset=['affinity_numeric']).copy()
    
    log.info(f"Loaded {len(df)} rows, {len(df_clean)} have valid affinity values.")

    if df_clean.empty:
        log.warning("No valid affinity data to process.")
        return df_clean

    df_clean['MolWt'] = 0.0
    df_clean['LogP'] = 0.0
    df_clean['HBD'] = 0
    df_clean['HBA'] = 0
    df_clean['RotBonds'] = 0
    df_clean['TPSA'] = 0.0
    df_clean['Pass_Lipinski'] = False
    df_clean['Pass_Veber'] = False
    df_clean['Pass_All'] = False
    df_clean['Reason'] = ''

    # .at 접근자는 인덱스가 고유해야 하므로, 인덱스를 리셋합니다.
    df_clean = df_clean.reset_index(drop=True)

    for idx, row in tqdm(df_clean.iterrows(), total=len(df_clean), desc="Calculating ADMET properties", unit="entry"):
        seq = row[SEQ_COL]
        reasons = []
        if not isinstance(seq, str) or len(seq) == 0:
            reasons = ['Invalid sequence']
        else:
            try:
                m = MolFromSequence(seq)
                if m is None:
                    reasons = ['Parse error']
                else:
                    mw = Descriptors.MolWt(m)
                    logp = Descriptors.MolLogP(m)
                    hbd = Lipinski.NumHDonors(m)
                    hba = Lipinski.NumHAcceptors(m)
                    rot = Descriptors.NumRotatableBonds(m)
                    tpsa = Descriptors.TPSA(m)
                    
                    df_clean.at[idx, 'MolWt'] = round(mw, 2)
                    df_clean.at[idx, 'LogP'] = round(logp, 2)
                    df_clean.at[idx, 'HBD'] = hbd
                    df_clean.at[idx, 'HBA'] = hba
                    df_clean.at[idx, 'RotBonds'] = rot
                    df_clean.at[idx, 'TPSA'] = round(tpsa, 2)
                    
                    lip_ok = (mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10)
                    veb_ok = (rot <= 10 and tpsa <= 140)
                    
                    df_clean.at[idx, 'Pass_Lipinski'] = lip_ok
                    df_clean.at[idx, 'Pass_Veber'] = veb_ok
                    df_clean.at[idx, 'Pass_All'] = (lip_ok and veb_ok)
                    
                    if not lip_ok: reasons.append('Lipinski fail')
                    if not veb_ok: reasons.append('Veber fail')
                    if lip_ok and veb_ok: reasons.append('Pass all')
                    del m
            except Exception as e:
                reasons = [f"RDKit error: {str(e)}"]
        df_clean.at[idx, 'Reason'] = '; '.join(reasons)

    log.info(f"ADMET calculation complete: {len(df_clean)} entries processed")
    df_clean = df_clean.drop(columns=['affinity_numeric']) # 임시 컬럼 제거
    
    # STEP1_OUT_FILE 저장은 이제 선택사항이지만, 디버깅을 위해 유지할 수 있습니다.
    df_clean.to_csv(output_path, index=False)
    log.info(f"Saved intermediate file to '{output_path}'")
    return df_clean

# 아래 step2, step3 함수들은 이제 main에서 호출되지 않습니다.
# def step2_cluster_stage1(...):
# def step3_cluster_stage2(...):

def step4_save_final_results(input_df, output_file):
    log.info("Saving final results file")
    df = input_df.copy()
    
    # affinity 값을 소수점 4자리까지 포맷팅합니다.
    df['affinity'] = pd.to_numeric(df['affinity']).map('{:.4f}'.format)
    
    log.info(f"Processed {len(df)} candidates for final output.")
    df.to_csv(output_file, index=False)
    log.info(f"Saved final results to '{output_file}'")
    return df

# --- ❗ [수정] 메인 실행 로직 단순화 ---
def main():
    try:
        log.info("Starting simplified pipeline")

        # 1단계: 모든 유효한 친화도 결과에 대해 ADMET 특성을 계산합니다.
        df_with_admet = step1_physchem_filter(CSV_PATH, STEP1_OUT_FILE)
        
        gc.collect()

        # 2단계: ADMET 특성이 추가된 전체 데이터프레임을 최종 결과 파일로 저장합니다.
        # (기존의 step2, step3 클러스터링 로직은 건너뜁니다.)
        final_df = step4_save_final_results(df_with_admet, FINAL_OUT_FILE)
        
        log.info("Pipeline completed successfully")
        print("\nFinal Results:")
        # 최종 데이터프레임이 비어있지 않은 경우에만 출력합니다.
        if not final_df.empty:
            print(final_df.to_string())
        else:
            print("No valid data to display.")

    except Exception as e:
        log.error(f"Pipeline error: {str(e)}", exc_info=True)
        print(f"Pipeline error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
