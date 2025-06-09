# TankBind GUI: 단백질-펩타이드 친화도 예측 및 도킹 시뮬레이션

TankBind의 빠른 결합 친화도 예측과 AutoDock Vina의 정밀한 도킹 시뮬레이션을 통합한 연구용 GUI 애플리케이션입니다. 복잡한 명령어 없이 그래픽 인터페이스로 펩타이드 서열을 여러 단백질 구조에 대해 스크리닝하고, 고친화도 후보에 대해 도킹 시뮬레이션과 3D 시각화를 즉시 수행할 수 있습니다.

---

## 주요 기능

* **결합 친화도 예측**: TankBind 모델로 단백질-펩타이드 결합 친화도(pKd)를 빠르게 예측합니다.
* **다중 타겟 스크리닝**: 한 펩타이드 서열을 폴더 내 모든 PDB 파일에 대해 일괄 스크리닝할 수 있습니다.
* **ADME/T 예측**: `pipeline.py`를 통해 Lipinski's Rule of Five, Veber's Rule 등 약물 유사성을 계산합니다.
* **원클릭 도킹 & 시각화**:
    * **AutoDock Vina**: 결과 테이블에서 버튼 클릭으로 도킹 시뮬레이션을 자동 실행합니다.
    * **PyMOL**: 도킹 후 단백질-리간드 결합 모습을 `.pse` 파일로 즉시 시각화합니다.
* **결과 관리**:
    * 예측 및 분석 결과를 정렬 가능한 표로 표시합니다.
    * CSV 또는 Excel 파일로 결과를 내보낼 수 있습니다.
    * 검색한 펩타이드 서열을 히스토리에 자동 저장합니다.

---

## 사전 요구사항

* **Git**: 코드 복제에 필요합니다.
* **Anaconda 또는 Miniconda**: Conda 환경 관리에 필요합니다.
* **NVIDIA GPU & CUDA Toolkit**: TankBind 모델 최적화 (CUDA 11.8 기반)에 필요합니다.

---

## 설치 및 실행 방법

1.  **프로젝트 복제**:
    ```bash
    git clone [https://github.com/ddangkong/tankbind_gui.git](https://github.com/ddangkong/tankbind_gui.git)
    cd tankbind_gui
    ```

2.  **Conda 환경 생성 및 활성화**:
    **참고**: 다음 명령어로 PyTorch, RDKit, PyMOL, Vina 등 필요한 라이브러리가 포함된 `tankbind-suite` 환경이 생성됩니다. 인터넷 속도에 따라 10~20분 소요될 수 있습니다.
    ```bash
    conda env create -f environment.yml
    conda activate tankbind-suite
    ```

3.  **데이터 준비 (PDB 파일 다운로드)**:
    * 스크리닝할 단백질의 PDB 파일을 로컬에 저장합니다.
    * RCSB Protein Data Bank ([https://www.rcsb.org/](https://www.rcsb.org/))에서 다운로드하거나, PDBe ([https://www.ebi.ac.uk/pdbe/](https://www.ebi.ac.uk/pdbe/))나 PDBj ([https://pdbj.org/](https://pdbj.org/)) 같은 미러 사이트를 사용할 수 있습니다.
    * 다운로드한 PDB 파일을 다음 단계에서 지정할 폴더에 저장합니다.

4.  **스크립트 내부 경로 설정 (필수)**:
    **중요**: 실행 전 `tankbind_gui.py`의 상단 경로 설정을 로컬 환경에 맞게 수정해야 합니다.
    ```python
    PDB_ROOT = "C:/Users/YourUser/Documents/TestPDBs/"  # PDB 파일 저장 폴더
    CENTER_PATH = "C:/Users/YourUser/Documents/center_coords.txt"  # 중심 좌표 파일
    TANKBIND_DIR = "C:/Users/YourUser/Documents/TankBind"  # TankBind 모델 폴더
    OUT_DIR = "C:/Users/YourUser/Documents/TankBind_Results"  # 결과 저장 폴더
    ```
    `tankbind-suite` 환경을 사용하므로 `PYMOL_PYTHON`, `AUTODOCK_VINA_EXECUTABLE` 등은 수정할 필요 없습니다.

5.  **프로그램 실행**:
    ```bash
    python tankbind_gui.py
    ```

---

## 라이선스

MIT 라이선스. 자세한 내용은 `LICENSE` 파일을 참조하세요.

---

## 감사의 말

* TankBind 연구팀 ([https://github.com/luwei0917/TankBind](https://github.com/luwei0917/TankBind))
* RDKit ([https://www.rdkit.org/](https://www.rdkit.org/)), PyMOL ([https://pymol.org/2/](https://pymol.org/2/)), AutoDock Vina ([https://vina.scripps.edu/](https://vina.scripps.edu/)) 등 오픈소스 도구


---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# TankBind GUI: Protein-Peptide Affinity Prediction and Docking Simulation

This research-oriented GUI application combines **TankBind's rapid binding affinity prediction** with **AutoDock Vina's precise docking simulation**. It offers a user-friendly graphical interface to screen a peptide sequence against multiple protein structures, instantly performing docking simulations and 3D visualizations for high-affinity candidates, all without complex command-line operations.

---

## Key Features

* **Binding Affinity Prediction**: Quickly predicts protein-peptide binding affinity (pKd) using the TankBind model.
* **Multi-Target Screening**: Screen a single peptide sequence against all PDB files within a specified folder.
* **ADME/T Prediction**: Calculates drug-likeness properties (e.g., Lipinski's Rule of Five, Veber's Rule) via `pipeline.py`.
* **One-Click Docking & Visualization**:
    * **AutoDock Vina**: Automatically runs docking simulations with a single click directly from the results table.
    * **PyMOL**: Generates `.pse` files for immediate 3D visualization of protein-ligand binding poses post-docking.
* **Results Management**:
    * Displays prediction and analysis results in a sortable table.
    * Allows export of results to **CSV or Excel files**.
    * Automatically saves searched peptide sequences to a history panel.

---

## Prerequisites

* **Git**: Essential for cloning the project repository.
* **Anaconda or Miniconda**: Required for managing the Conda environment.
* **NVIDIA GPU & CUDA Toolkit**: Optimized for TankBind model performance (CUDA 11.8-based).

---

## Installation and Usage

1.  **Clone the Project**:
    ```bash
    git clone [https://github.com/ddangkong/tankbind_gui.git](https://github.com/ddangkong/tankbind_gui.git)
    cd tankbind_gui
    ```

2.  **Create and Activate the Conda Environment**:
    * **Note**: This command will create a `tankbind-suite` environment, installing all necessary libraries including PyTorch, RDKit, PyMOL, and Vina. This process may take **10-20 minutes** depending on your internet speed.
    ```bash
    conda env create -f environment.yml
    conda activate tankbind-suite
    ```

3.  **Data Preparation (Download PDB Files)**:
    * Download your target protein PDB files and save them locally.
    * You can obtain PDB files from the **RCSB Protein Data Bank** (https://www.rcsb.org/) or mirror sites like **PDBe** (https://www.ebi.ac.uk/pdbe/) or **PDBj** (https://pdbj.org/).
    * Place the downloaded PDB files into the folder that you will specify for `PDB_ROOT` in the next step.

4.  **Configure Paths in the Script (Mandatory)**:
    * **Important**: Before running, you **must edit** the path settings at the top of `tankbind_gui.py` to match your local system's directory structure:
    ```python
    PDB_ROOT = "C:/Users/YourUser/Documents/TestPDBs/"  # Folder where PDB files are stored
    CENTER_PATH = "C:/Users/YourUser/Documents/center_coords.txt"  # File for center coordinates
    TANKBIND_DIR = "C:/Users/YourUser/Documents/TankBind"  # TankBind model folder
    OUT_DIR = "C:/Users/YourUser/Documents/TankBind_Results"  # Folder for saving results
    ```
    * You do **not** need to modify `PYMOL_PYTHON` or `AUTODOCK_VINA_EXECUTABLE` if you are using the `tankbind-suite` environment.

5.  **Run the Application**:
    ```bash
    python tankbind_gui.py
    ```

---

## License

This project is released under the **MIT License**. Please refer to the `LICENSE` file for more details.

---

## Acknowledgments

We extend our gratitude to:

* The **TankBind research team** (https://github.com/luwei0917/TankBind) for their foundational work.
* The developers of essential open-source tools: **RDKit** (https://www.rdkit.org/), **PyMOL** (https://pymol.org/2/), and **AutoDock Vina** (https://vina.scripps.edu/).
