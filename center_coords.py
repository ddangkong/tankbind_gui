import sys

# Python 3.9 이상에서 sys.setcheckinterval가 제거되었으므로, setswitchinterval을 대신 사용하도록 monkey-patch합니다.
if not hasattr(sys, "setcheckinterval"):
    sys.setcheckinterval = sys.setswitchinterval

import os
import numpy as np

# pymol2 인터페이스 사용. 이는 apt-get으로 설치된 PyMOL 환경에서 사용하기 좋습니다.
import pymol2

pdb_dir = "/home/jovyan/Chembl(효능데이터베이스)/PDBs/"
output_file = "/home/jovyan/Chembl(효능데이터베이스)/center_coords.txt"

def process_pdb(pdb_file, cmd):
    try:
        # 각 파일마다 상태를 초기화합니다.
        cmd.reinitialize()
        pdb_path = os.path.join(pdb_dir, pdb_file)
        cmd.load(pdb_path, "protein")
        
        # 리간드(organic or het) 선택
        cmd.select("ligands", "organic or het")
        if cmd.count_atoms("ligands") > 0:
            cmd.select("pockets", "(polymer) and byres (ligands around 8)")
        else:
            cmd.select("pockets", "(polymer) and byres (polymer around 8)")
        
        # 선택된 영역의 좌표를 추출하여 중심값 계산
        coords = cmd.get_coords("pockets")
        if coords is None or len(coords) == 0:
            raise ValueError("No coordinates found")
        center = np.mean(coords, axis=0)
        return f"{pdb_file}: {center.tolist()}"
    except Exception as e:
        return f"{pdb_file}: Failed ({str(e)})"
    finally:
        cmd.delete("all")

if __name__ == "__main__":
    results = []
    
    # pymol2 인터페이스로 PyMOL 인스턴스를 생성합니다.
    with pymol2.PyMOL() as pymol:
        cmd = pymol.cmd
        for f in os.listdir(pdb_dir):
            if f.endswith(".pdb"):
                result = process_pdb(f, cmd)
                print(result)  # 진행 상황 출력
                results.append(result)
                
    with open(output_file, 'w') as f:
        f.write("\n".join(results))
        
    print(f"Saved coordinates to {output_file}")
