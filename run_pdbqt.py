#!/usr/bin/env python
# -*- coding: utf-8 -*-
# run_pdbqt.py (수정 완료 버전)

import os
import sys
import subprocess
import argparse
import logging
from typing import Optional
import time

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("오류: RDKit 라이브러리를 찾을 수 없습니다.", file=sys.stderr)
    sys.exit(1)

# --- 로거 설정 (기존과 동일) ---
log_formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(filename)s(%(lineno)d): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setFormatter(log_formatter)
if not logger.hasHandlers():
    logger.addHandler(console_handler)

# --- ❗ [수정 1] PyMOL GUI 충돌을 유발하는 코드 블록 전체 삭제 ---
# 이 스크립트 내에서 더 이상 'import pymol'을 직접 사용하지 않습니다.
# 대신, 외부 명령어로 PyMOL을 헤드리스(Headless) 모드로 호출합니다.

def peptide_to_mol_for_pdbqt(seq: str, max_try: int = 5) -> Optional[Chem.Mol]:
    """펩타이드 서열로부터 3D 구조의 RDKit Mol 객체를 생성합니다. (기존 코드와 동일)"""
    seq = seq.strip().upper()
    AA_MAP = {"O": "P"}
    seq2 = "".join(AA_MAP.get(a, a) for a in seq)
    logger.info(f"[MolGen-PDBQT] 처리된 서열: {seq2}")
    
    valid_aas = "ACDEFGHIKLMNPQRSTVWY"
    if not all(char in valid_aas for char in seq2):
        logger.error(f"[MolGen-PDBQT] 잘못된 아미노산 코드 포함: {seq2}. 유효한 코드: {valid_aas}")
        return None
    
    mol = Chem.MolFromFASTA(seq2)
    if mol is None:
        mol = Chem.MolFromSequence(seq2)
        if mol is None:
            logger.error(f"[MolGen-PDBQT] MolFromFASTA 및 MolFromSequence 모두 실패: {seq2}")
            return None
    
    num_atoms_initial = mol.GetNumAtoms()
    if num_atoms_initial == 0:
        logger.error(f"[MolGen-PDBQT] 원자 수가 0인 분자가 생성됨: {seq2}")
        return None
    mol = Chem.AddHs(mol)
    logger.info(f"[MolGen-PDBQT] 수소 추가 후 원자 수: {mol.GetNumAtoms()}")

    params = AllChem.ETKDGv3()
    params.useRandomCoords = True
    params.maxIterations = 2000
    
    embedding_successful = False
    for i in range(max_try):
        params.randomSeed = (os.getpid() + i * 42 + int(time.time() * 1000)) % (2**32)
        logger.info(f"[MolGen-PDBQT] 3D 임베딩 시도 {i+1}/{max_try} (서열: {seq2}, 시드: {params.randomSeed})")
        mol.RemoveAllConformers()
        conf_id = AllChem.EmbedMolecule(mol, params)

        if conf_id == 0:
            logger.info(f"[MolGen-PDBQT] EmbedMolecule 성공 (시도 {i+1})")
            actual_conf_id = mol.GetConformers()[0].GetId() if mol.GetNumConformers() > 0 else -1
            if actual_conf_id == -1:
                logger.warning(f"[MolGen-PDBQT] EmbedMolecule 성공했으나 컨포머 ID 얻기 실패.")
                continue

            try:
                opt_result = AllChem.UFFOptimizeMolecule(mol, maxIters=1000, confId=actual_conf_id)
                if opt_result == 0:
                    logger.info(f"[MolGen-PDBQT] 3D 임베딩 및 UFF 최적화 성공: {seq2}")
                    embedding_successful = True
                    break
                else:
                    logger.warning(f"[MolGen-PDBQT] UFFOptimizeMolecule 최적화 실패 (결과: {opt_result}, 시도 {i+1})")
            except RuntimeError as e:
                logger.warning(f"[MolGen-PDBQT] UFFOptimizeMolecule 런타임 오류 (시도 {i+1}): {e}")
            except Exception as e_ff:
                logger.warning(f"[MolGen-PDBQT] UFF ForceField 관련 오류 (시도 {i+1}): {e_ff}")
        else:
            logger.warning(f"[MolGen-PDBQT] EmbedMolecule 실패 (결과 코드: {conf_id}, 시도 {i+1})")

    if not embedding_successful or mol.GetNumConformers() == 0:
        logger.error(f"[MolGen-PDBQT] {max_try}번의 시도 후 3D 임베딩 최종 실패: {seq2}")
        return None
        
    logger.info(f"[MolGen-PDBQT] {seq2}에 대한 Mol 객체 생성 완료.")
    return mol

def run_shell_command(command_list, working_dir=None) -> bool:
    """쉘 명령어를 실행하고 결과를 로깅합니다. (기존 코드와 동일)"""
    logger.info(f"명령어 실행: {' '.join(command_list)}")
    try:
        process = subprocess.run(command_list, capture_output=True, text=True, check=False, timeout=300) # 타임아웃 5분
        if process.stdout: logger.info(f"STDOUT:\n{process.stdout.strip()}")
        # Vina는 정상 작동 시에도 stderr로 진행률을 출력하므로 warning으로 처리합니다.
        if process.stderr: logger.warning(f"STDERR:\n{process.stderr.strip()}")
        
        # Vina는 성공해도 0이 아닌 값을 반환할 때가 있으므로, stderr에 'Error'가 없을 때도 성공으로 간주할 수 있습니다.
        # 여기서는 단순하게 returncode로 확인합니다.
        if process.returncode != 0:
             logger.error(f"명령어 실행이 0이 아닌 코드로 종료되었습니다: {process.returncode}")
             return False
        return True
    except subprocess.TimeoutExpired:
        logger.error(f"명령어 실행 시간 초과: {' '.join(command_list)}")
        return False
    except FileNotFoundError:
        logger.error(f"명령어 실행 파일({command_list[0]})을 찾을 수 없습니다.")
        return False
    except Exception as e:
        logger.error(f"명령어 실행 중 예외 발생: {e}", exc_info=True)
        return False

def convert_to_pdbqt(input_file: str, output_file: str, obabel_exec: str, is_receptor: bool = False) -> bool:
    """Open Babel을 사용하여 PDB/SDF를 PDBQT로 변환합니다. (기존 코드와 동일)"""
    if not os.path.exists(input_file):
        logger.error(f"[PDBQT Convert] 입력 파일 없음: {input_file}")
        return False
    
    # 리간드는 PDB에서, 수용체는 SDF에서 시작하는 등의 차이가 있을 수 있으므로 input 포맷을 명시해주는 것이 더 안정적일 수 있습니다.
    # 예: obabel_cmd = [obabel_exec, "-ipdb", input_file, "-opdbqt", "-O", output_file]
    obabel_cmd = [obabel_exec, input_file, "-O", output_file]
    if is_receptor:
        # 수용체 준비: 비극성 수소 병합 및 원자 유형 할당
        obabel_cmd.extend(["-xr", "-A", "hydrogens"])
    else: # Ligand
        # 리간드 준비: 모든 수소 추가 및 7.4 pH에서의 양성자화 상태 설정
        obabel_cmd.extend(["-xh", "-p", "7.4"])
    
    logger.info(f"[PDBQT Convert] {'수용체' if is_receptor else '리간드'} PDBQT 변환: {input_file} -> {output_file}")
    if run_shell_command(obabel_cmd):
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            logger.info(f"[PDBQT Convert] PDBQT 변환 성공: {output_file}")
            return True
        else:
            logger.error(f"[PDBQT Convert] PDBQT 변환 후 파일이 없거나 비어있음: {output_file}")
            return False
    else:
        logger.error(f"[PDBQT Convert] PDBQT 변환 실패 (obabel 실행 오류).")
        return False

def create_and_run_pymol_session(protein_file: str, poses_file: str, session_output_file: str) -> None:
    """
    ❗ [수정 2] PyMOL 스크립트(.pml) 파일을 생성하고 헤드리스 모드로 실행하여 세션(.pse) 파일을 저장하는 새 함수.
    """
    pml_script_content = f"""
reinitialize
load {os.path.abspath(protein_file)}, receptor
load {os.path.abspath(poses_file)}, ligand_poses

color gray80, receptor
show_as cartoon, receptor
remove resn HOH # 물 분자 제거
remove hydro # 수소 숨기기

# 모든 리간드 포즈를 다른 색으로 표시하고 스타일 적용
split_states ligand_poses
util.cbag ligand_poses_0*

show_as sticks, ligand_poses
set stick_radius, 0.2, ligand_poses

# 결합 부위 주변만 표시 (예시)
select binding_site, receptor within 5 of ligand_poses
show lines, binding_site
show sticks, binding_site

zoom ligand_poses, 5
center ligand_poses and state 1

# PyMOL 세션 파일 저장
save {os.path.abspath(session_output_file)}
quit
"""
    pml_script_path = os.path.join(os.path.dirname(session_output_file), "visualize.pml")
    try:
        with open(pml_script_path, "w") as f:
            f.write(pml_script_content)
        logger.info(f"PyMOL 스크립트 생성: {pml_script_path}")

        # PyMOL 실행 파일은 환경 PATH에 있거나 절대 경로여야 합니다.
        pymol_exec = "pymol"
        # -c: 커맨드라인, -q: 조용한 시작, -k: 스크립트 후 종료, -r: 스크립트 실행
        pymol_cmd = [pymol_exec, "-c", "-q", "-k", "-r", pml_script_path]
        
        logger.info(f"PyMOL 헤드리스 실행: {' '.join(pymol_cmd)}")
        if run_shell_command(pymol_cmd):
            logger.info(f"PyMOL 세션 파일 저장 성공: {session_output_file}")
        else:
            logger.error("PyMOL 세션 파일 저장 실패.")

    except Exception as e:
        logger.error(f"PyMOL 세션 생성 중 오류: {e}")


def main():
    parser = argparse.ArgumentParser(description="펩타이드 PDBQT 생성, Vina 도킹 및 PyMOL 세션 생성 스크립트")
    # 인자 파싱 부분은 기존과 동일합니다. 메인 GUI에서 모든 인자를 잘 전달해 줄 것입니다.
    parser.add_argument("--peptide_seq", required=True, help="리간드 펩타이드 서열")
    parser.add_argument("--pdb_id", required=True, help="단백질 PDB ID (파일 네이밍 및 원본 PDB 파일명)")
    parser.add_argument("--protein_pdb_file", required=True, help="원본 단백질 PDB 파일 전체 경로")
    parser.add_argument("--base_out_dir", required=True, help="결과 파일을 저장할 기본 디렉토리")
    parser.add_argument("--obabel_executable", required=True, help="Open Babel 실행 파일 경로")
    parser.add_argument("--vina_executable", required=True, help="AutoDock Vina 실행 파일 경로")
    parser.add_argument("--center_x", type=float, required=True, help="도킹 박스 중심 X")
    parser.add_argument("--center_y", type=float, required=True, help="도킹 박스 중심 Y")
    parser.add_argument("--center_z", type=float, required=True, help="도킹 박스 중심 Z")
    parser.add_argument("--box_size_x", type=float, default=25.0, help="Vina 검색 박스 X 크기")
    parser.add_argument("--box_size_y", type=float, default=25.0, help="Vina 검색 박스 Y 크기")
    parser.add_argument("--box_size_z", type=float, default=25.0, help="Vina 검색 박스 Z 크기")
    parser.add_argument("--num_modes", type=int, default=9, help="Vina 결합 모드 수")
    parser.add_argument("--exhaustiveness", type=int, default=8, help="Vina 검색 철저도")
    
    args = parser.parse_args()
    logger.info(f"run_pdbqt_vina.py 실행 시작. 인자: {args}")

    specific_out_dir = os.path.join(args.base_out_dir, f"{args.pdb_id}_{args.peptide_seq.replace('/','_')}")
    os.makedirs(specific_out_dir, exist_ok=True)
    logger.info(f"결과 저장 경로: {specific_out_dir}")

    # 1. 리간드 준비
    mol_ligand = peptide_to_mol_for_pdbqt(args.peptide_seq)
    if not mol_ligand:
        logger.error(f"리간드 Mol 객체 생성 실패: {args.peptide_seq}")
        print(f"FAILURE: Ligand Mol generation failed for {args.peptide_seq}", file=sys.stderr)
        sys.exit(1)

    ligand_temp_pdb_file = os.path.join(specific_out_dir, f"{args.peptide_seq}_temp_ligand.pdb")
    ligand_pdbqt_file = os.path.join(specific_out_dir, f"{args.peptide_seq}_ligand.pdbqt")
    try:
        Chem.MolToPDBFile(mol_ligand, ligand_temp_pdb_file)
        logger.info(f"리간드 임시 PDB 파일 저장: {ligand_temp_pdb_file}")
    except Exception as e:
        logger.error(f"리간드 임시 PDB 파일 저장 실패: {e}")
        print(f"FAILURE: Ligand PDB saving failed for {args.peptide_seq}", file=sys.stderr)
        sys.exit(1)
    
    if not convert_to_pdbqt(ligand_temp_pdb_file, ligand_pdbqt_file, args.obabel_executable, is_receptor=False):
        logger.error(f"리간드 PDBQT 변환 실패: {ligand_temp_pdb_file}")
        print(f"FAILURE: Ligand PDBQT conversion failed for {args.peptide_seq}", file=sys.stderr)
        sys.exit(1)

    # 2. 수용체 준비
    receptor_pdbqt_file = os.path.join(specific_out_dir, f"{args.pdb_id}_receptor.pdbqt")
    if not convert_to_pdbqt(args.protein_pdb_file, receptor_pdbqt_file, args.obabel_executable, is_receptor=True):
        logger.error(f"수용체 PDBQT 변환 실패: {args.protein_pdb_file}")
        print(f"FAILURE: Receptor PDBQT conversion failed for {args.pdb_id}", file=sys.stderr)
        sys.exit(1)

    # 3. Vina 실행
    vina_out_poses_pdbqt = os.path.join(specific_out_dir, f"{args.pdb_id}_{args.peptide_seq}_vina_poses.pdbqt")
    vina_log_file = os.path.join(specific_out_dir, f"{args.pdb_id}_{args.peptide_seq}_vina.log")

    vina_cmd = [
        args.vina_executable,
        "--receptor", receptor_pdbqt_file,
        "--ligand", ligand_pdbqt_file,
        "--center_x", str(args.center_x),
        "--center_y", str(args.center_y),
        "--center_z", str(args.center_z),
        "--size_x", str(args.box_size_x),
        "--size_y", str(args.box_size_y),
        "--size_z", str(args.box_size_z),
        "--out", vina_out_poses_pdbqt,
        "--log", vina_log_file,
        "--num_modes", str(args.num_modes),
        "--exhaustiveness", str(args.exhaustiveness)
    ]
    logger.info(f"Vina 실행 준비. 명령어: {' '.join(vina_cmd)}")
    if not run_shell_command(vina_cmd, working_dir=specific_out_dir):
        logger.error("AutoDock Vina 실행 실패.")
        print(f"FAILURE: Vina execution failed for {args.pdb_id}-{args.peptide_seq}", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(vina_out_poses_pdbqt) or os.path.getsize(vina_out_poses_pdbqt) == 0:
        logger.error(f"Vina 출력 파일({vina_out_poses_pdbqt})이 생성되지 않았거나 비어있습니다. 로그: {vina_log_file}")
        print(f"FAILURE: Vina output file missing or empty for {args.pdb_id}-{args.peptide_seq}", file=sys.stderr)
        sys.exit(1)
    
    logger.info(f"Vina 도킹 완료. 결과: {vina_out_poses_pdbqt}")

    # ❗ [수정 3] PyMOL 시각화 부분을 새로운 함수 호출로 변경
    pymol_session_file = os.path.join(specific_out_dir, f"{args.pdb_id}_{args.peptide_seq}_session.pse")
    
    logger.info("PyMOL 세션 파일 생성 시작...")
    create_and_run_pymol_session(
        protein_file=args.protein_pdb_file,
        poses_file=vina_out_poses_pdbqt,
        session_output_file=pymol_session_file
    )

    print(f"SUCCESS: Vina docking and PyMOL session generation complete. Output dir: {specific_out_dir}")
    sys.exit(0)


if __name__ == "__main__":
    main()
