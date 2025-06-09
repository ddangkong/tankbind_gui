TankBind GUI: 단백질-펩타이드 친화도 예측 및 도킹 시뮬레이션
TankBind의 빠른 결합 친화도 예측과 AutoDock Vina의 정밀한 도킹 시뮬레이션 기능을 통합한 연구용 GUI 애플리케이션

이 프로그램은 복잡한 명령어 없이, 그래픽 인터페이스를 통해 하나의 펩타이드 서열을 다수의 단백질 구조에 대해 스크리닝하고, 높은 친화도를 보이는 후보에 대해 즉시 도킹 시뮬레이션 및 3D 시각화를 수행할 수 있도록 돕습니다.

✨ 주요 기능
결합 친화도 예측: TankBind (https://github.com/luwei0917/TankBind) 모델을 사용하여 단백질-펩타이드 간의 결합 친화도(pK_d)를 빠르게 예측합니다.
다중 타겟 스크리닝: 단일 펩타이드 서열을 폴더 내 모든 PDB 파일에 대해 일괄적으로 스크리닝합니다.
ADME/T 예측: pipeline.py 스크립트 연동을 통해 Lipinski's Rule of Five, Veber's Rule 등 약물 유사성(ADME/T) 프로퍼티를 계산합니다.
원클릭 도킹 & 시각화:
AutoDock Vina (https://vina.scripps.edu/): 예측 결과 테이블에서 버튼 하나만 누르면 Vina를 이용한 도킹 시뮬레이션이 자동으로 실행됩니다.
PyMOL (https://pymol.org/2/): 도킹 완료 후 단백질과 리간드의 결합 모습을 즉시 확인할 수 있는 .pse 세션 파일이 자동 생성됩니다.
결과 관리:
모든 예측 및 분석 결과는 표 형태로 정렬하여 볼 수 있습니다.
결과를 CSV 또는 Excel 파일로 손쉽게 내보낼 수 있습니다.
검색했던 펩타이드 서열은 히스토리에 자동 저장됩니다.
⚙️ 사전 요구사항
Git: 코드를 내려받기 위해 필요합니다.
Anaconda 또는 Miniconda: Conda 환경 관리를 위해 반드시 설치되어 있어야 합니다.
NVIDIA GPU & CUDA Toolkit: TankBind 모델은 GPU 환경에서 가장 효율적으로 작동합니다. (CUDA 11.8 기반으로 환경이 구성됩니다.)
🚀 설치 및 실행 방법
1단계: 프로젝트 복제 (Clone)
Bash

git clone https://github.com/ddangkong/tankbind_gui.git
cd tankbind_gui
2단계: Conda 환경 생성 및 활성화
참고: 아래 명령어로 PyTorch, RDKit, PyMOL, Vina 등 모든 필요한 라이브러리가 포함된 tankbind-suite 환경이 생성됩니다. 인터넷 속도에 따라 10~20분 소요될 수 있습니다.

Bash

conda env create -f environment.yml
conda activate tankbind-suite
3단계: 데이터 준비 (PDB 파일 다운로드)
이 애플리케이션을 사용하려면 스크리닝할 단백질의 PDB 파일들이 로컬 컴퓨터에 저장되어 있어야 합니다.

단백질 구조 파일(.pdb 형식)은 공식 RCSB Protein Data Bank (PDB) (https://www.rcsb.org/) 웹사이트에서 다운로드할 수 있습니다.

원하는 단백질을 검색하여 PDB 파일을 다운로드한 후, 다음 단계에서 PDB_ROOT 변수에 지정할 폴더 안에 저장해 주세요. 또는, PDBe (유럽) (https://www.pdbe.org/)나나) PDBj (일본) (https://pdbj.org/) 같은 미러 사이트를 이용할 수도 있습니다.

4단계: 스크립트 내부 경로 설정 (실행 전 필수!)
중요: 실행 전 <your_main_script_name>.py의 상단 경로 설정을 로컬 환경에 맞게 수정해야 합니다.

Python

# --- 경로 설정 ---
PDB_ROOT = "C:/Users/YourUser/Documents/TestPDBs/"      # 3단계에서 PDB 파일들을 저장한 폴더 경로
CENTER_PATH = "C:/Users/YourUser/Documents/center_coords.txt" # 중심 좌표 파일 경로
TANKBIND_DIR = "C:/Users/YourUser/Documents/TankBind"     # TankBind 모델 폴더 경로
OUT_DIR = "C:/Users/YourUser/Documents/TankBind_Results" # 결과 파일이 저장될 폴더 경로
# ... (다른 경로들도 필요시 수정) ...
참고: 통합 환경을 사용하므로 PYMOL_PYTHON이나 AUTODOCK_VINA_EXECUTABLE 등은 수정할 필요가 없습니다.

5단계: 프로그램 실행
모든 설정이 완료되었습니다. 아래 명령어로 프로그램을 실행합니다.

Bash

python <your_main_script_name>.py
📜 라이선스
이 프로젝트는 MIT 라이선스를 따릅니다. 자세한 내용은 LICENSE 파일을 참고하세요.

🙏 감사의 말
훌륭한 예측 모델을 제공해 준 TankBind (https://github.com/luwei0917/TankBind) 연구팀에 감사드립니다.
이 프로젝트는 RDKit, PyMOL, AutoDock Vina 등 여러 오픈소스 도구에 크게 의존하고 있습니다.
&lt;br>

&lt;br>

🇬🇧 English Version
TankBind GUI: Protein-Peptide Affinity Prediction and Docking Simulation
A research-oriented GUI application that integrates the rapid binding affinity prediction of TankBind with the precise docking simulation of AutoDock Vina.

This application enables users to screen a peptide sequence against multiple protein structures through a graphical interface, eliminating the need for complex command-line operations. It allows for immediate docking simulations and 3D visualization for high-affinity candidates.

✨ Key Features
Binding Affinity Prediction: Rapidly predicts the binding affinity (pK_d) between proteins and peptides using the TankBind (https://github.com/luwei0917/TankBind) model.
Multi-Target Screening: Screens a single peptide sequence against all PDB files in a specified folder in one batch.
ADME/T Prediction: Calculates drug-likeness properties, such as Lipinski's Rule of Five and Veber's Rule, through integration with the pipeline.py script.
One-Click Docking & Visualization:
AutoDock Vina (https://vina.scripps.edu/): Automatically runs a docking simulation with Vina at the click of a button in the results table.
PyMOL (https://pymol.org/2/): Automatically generates a .pse session file to instantly visualize the protein-ligand binding pose after docking is complete.
Results Management:
Displays all prediction and analysis results in a sortable table.
Easily exports results to CSV or Excel files.
Automatically saves searched peptide sequences to a history panel.
⚙️ Prerequisites
Git: Required to clone the repository.
Anaconda or Miniconda: Must be installed for Conda environment management.
NVIDIA GPU & CUDA Toolkit: The TankBind model operates most efficiently in a GPU environment. The provided environment is configured based on CUDA 11.8.
🚀 Installation and Usage
Step 1: Clone the Project
Bash

git clone https://github.com/ddangkong/tankbind_gui.git
cd tankbind_gui
Step 2: Create and Activate the Conda Environment
Note: The following command creates a Conda environment named tankbind-suite that includes all necessary libraries and programs (PyTorch, RDKit, PyMOL, Vina, etc.). This process may take 10-20 minutes depending on your internet connection.

Bash

# Create the Conda environment from the .yml file
conda env create -f environment.yml

# Activate the newly created environment
conda activate tankbind-suite
Step 3: Data Preparation (Download PDB Files)
This application requires local PDB files of the target proteins for screening.

Protein structure files (in .pdb format) can be downloaded from the official RCSB Protein Data Bank (PDB) (https://www.rcsb.org/).

Search for your proteins of interest, download their PDB files, and place them into the folder you will specify for the PDB_ROOT variable in the next step. Alternatively, you can use mirror sites such as PDBe (Europe) (https://www.pdbe.org/) or PDBj (Japan) (https://pdbj.org/).

Step 4: Configure Paths in the Script (Required!)
Important: Before running the application for the first time, you must edit the data paths in the main script to match your local system.

Open the main Python script (<your_main_script_name>.py) and modify the # --- 경로 설정 --- (Path Settings) section at the top.

Python

# --- Path Settings ---
PDB_ROOT = "C:/Users/YourUser/Documents/TestPDBs/"      # Path to the folder where you saved PDB files in step 3
CENTER_PATH = "C:/Users/YourUser/Documents/center_coords.txt" # Path to the center coordinates file
TANKBIND_DIR = "C:/Users/YourUser/Documents/TankBind"     # Path to the TankBind model folder
OUT_DIR = "C:/Users/YourUser/Documents/TankBind_Results" # Path where result files will be saved
# ... (Modify other paths as needed) ...
Note: Since we are using a unified environment, you do not need to modify paths like PYMOL_PYTHON or AUTODOCK_VINA_EXECUTABLE.

Step 5: Run the Application
Once all settings are configured, run the application with the following command:

Bash

python <your_main_script_name>.py
📜 License
Distributed under the MIT License. See the LICENSE file for more information.

🙏 Acknowledgments
Special thanks to the TankBind (https://github.com/luwei0917/TankBind) research team for providing their excellent prediction model.
This project heavily relies on open-source tools like RDKit, PyMOL, and AutoDock Vina.
