# 🔭 ThyroScope: Clinical-Grade WES Pipeline for Thyroid and Parathyroid Disorders

[![Docker Image Version](https://img.shields.io/docker/v/hanyunseo01/thyroid_pipeline/v23?color=blue&label=Docker%20Image)](https://hub.docker.com/r/hanyunseo01/thyroid_pipeline)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20macOS%20%7C%20Linux-lightgrey)](https://www.docker.com/)

**ThyroScope** is a streamlined, containerized bioinformatics pipeline designed to "scope out" and analyze thyroid-specific variants from **Whole Exome Sequencing (WES)** data.

By leveraging a **"Lightweight Containerization Strategy"** and a **"Virtual Panel approach"**, ThyroScope allows clinicians to bypass complex command-line interfaces. It automates the entire workflow from raw FASTQ to a clinical-grade Excel report with a single click, focusing specifically on 25 high-priority genes associated with thyroid and parathyroid disorders.

![ThyroScope Pipeline Diagram](canvas-image-1-1768740173287.png)

---

## 🚀 Key Features

* **🔭 Targeted Precision:** Implements a "Virtual Panel" to filter WES data for 25 key genes (e.g., *BRAF, RET, TSHR*), drastically reducing incidental findings.
* **🐳 Dockerized & Reproducible:** Encapsulates GATK 4.4, BWA, Samtools, and Python dependencies in a portable Docker image, ensuring identical results on any OS.
* **📉 Lightweight Architecture:** Heavy reference databases (hg38 Ref, SnpEff DB) are externalized, keeping the Docker image light (~2GB) and easy to distribute.
* **📊 Comprehensive QC Reports:**
    * **MultiQC:** Aggregates quality metrics (Trimming, Alignment, Duplication) into a single interactive HTML report.
    * **Mosdepth:** Provides rapid coverage statistics to verify diagnostic depth (e.g., >20x) for target regions.
* **🌍 Population & Clinical Annotation:**
    * **ClinVar & dbSNP:** Auto-fetches the latest pathogenicity data and rsIDs via `MyVariant.info`.
    * **gnomAD:** Includes Global Allele Frequency (AF) to help differentiate rare variants from common polymorphisms.
    * **Functional Prediction:** SIFT and PolyPhen-2 scores included.
* **⚡ One-Click Automation:** Includes a Windows Batch script (`.bat`) for "Drag-and-Drop" style execution.

---

## 🔬 Pipeline Workflow (Methods)
ThyroScope automates the following bioinformatics steps in a sequential manner:

1.  **QC & Trimming** ✂️
    * **Tools:** `FastQC`, `Trimmomatic`
    * **Details:** Adapter removal, Quality trimming (`SlidingWindow:4:15`).
2.  **Alignment** 🧬
    * **Tool:** `BWA-MEM`
    * **Reference:** Aligned to the **GRCh38 (hg38)** reference genome.
3.  **Post-Processing** 🧹
    * **Tools:** `Samtools` (Sort/Index), `GATK MarkDuplicates`.
    * **Details:** Sorting BAM files and marking PCR duplicates to ensure accurate variant calling.
4.  **Coverage Analysis** 📉
    * **Tool:** `Mosdepth`
    * **Details:** Rapid depth-of-coverage check specifically for the targeted regions.
5.  **Variant Calling** 🔍
    * **Tool:** `GATK HaplotypeCaller`
    * **Details:** Calling germline variants with `--interval-padding 100` to capture essential splice site regions.
6.  **Annotation** 📝
    * **Tools:** `SnpEff` (HGVS notation), `MyVariant.info` API.
    * **Databases:** ClinVar, gnomAD (Global AF), dbSNP (rsID), dbNSFP (SIFT/PolyPhen).
7.  **Reporting** 📊
    * **Output:** A custom Python script aggregates all data into a **Hybrid Clinical Excel Report** and generates a **MultiQC** HTML summary.

---

## 🛠️ System Requirements

* **OS:** Windows 10/11 (Pro/Home), macOS, or Linux.
* **Software:** [Docker Desktop](https://www.docker.com/products/docker-desktop) (Must be installed and running).
* **Hardware:** Minimum 16GB RAM recommended (for GATK & Java heap).
* **Storage:** At least 50GB of free space (for reference data and WES output).

---

## 📥 Installation & Setup

Since ThyroScope uses a **lightweight strategy**, you must download the reference bundles separately.

1.  **Clone this Repository** (or download the ZIP):
    ```bash
    git clone [https://github.com/hanyunseo01/ThyroScope.git](https://github.com/hanyunseo01/ThyroScope.git)
    cd ThyroScope
    ```

2.  **Download Reference Data:**
    * Download the `ref` (BWA Indices) and `snpEff_db` folders from our repository storage.
    * [**📂 Download Link (Google Drive / Dropbox)**](#) *(Link to be updated)*
    * Place them in your project folder.

3.  **Directory Structure:**
    Ensure your folder looks exactly like this before running:
    ```text
    ThyroScope/
    ├── 📂 data/             <-- Put your FASTQ files here (e.g., Patient_1.fq.gz)
    ├── 📂 ref/              <-- Contains hg38.fasta, .bwt, .pac, etc.
    ├── 📂 snpEff_db/        <-- Contains 'hg38' folder
    ├── 📄 Run_Analysis.bat  <-- For Windows Users
    ├── 📄 pipeline_script.py
    └── 📄 Dockerfile
    ```

---

## 🏃‍♀️ How to Run

### 🖥️ Option A: Windows (Recommended for Clinicians)
We provide a batch script for a seamless experience.

1.  Place your paired-end FASTQ files (e.g., `PatientA_1.fq.gz`, `PatientA_2.fq.gz`) into the `data/` folder.
2.  Double-click **`Run_Analysis.bat`**.
3.  A terminal window will open, and the pipeline will start automatically.
4.  Wait for the message: `[INFO] Pipeline Completed Successfully!`

### 🐧 Option B: Linux / macOS
Run the Docker container manually using the following command:

```bash
# Move to your project directory
cd /path/to/ThyroScope

# Run the pipeline (v23)
docker run --rm \
  -v "${PWD}/data:/data" \
  -v "${PWD}/ref:/data/ref" \
  -v "${PWD}/snpEff_db:/pipeline/snpEff/data" \
  hanyunseo01/thyroid_pipeline:v23
```
---

## 📄 Output Files

After the analysis completes, check the `data/` folder for these key files:

| File Name | Description |
| :--- | :--- |
| **`*_Clinical_Report.xlsx`** | **The Final Report.** A comprehensive Excel file containing filtered variants annotated with Variant ID (rsID), ClinVar, gnomAD AF, SIFT/PolyPhen, and Impact. |
| **`*_MultiQC_Report.html`** | **Quality Control.** Interactive graphs showing read quality, mapping rates, and duplicate levels. Open with any web browser. |
| **`*_coverage.mosdepth.summary.txt`** | **Depth Statistics.** Shows how well the target genes were covered (e.g., mean depth, % bases > 20x). |
| `*.bam` / `*.vcf` | Intermediate alignment and variant calling files for further manual inspection (IGV). |

### 🩺 Inside the Clinical Report (.xlsx)

The **Clinical Report** is designed for immediate clinical interpretation. It aggregates data from **SnpEff**, **ClinVar**, **gnomAD**, and **dbSNP** into a single view.

| Column Category | Columns Included | Description |
| :--- | :--- | :--- |
| **Target Info** | `Gene`, `Transcript ID` | The gene symbol (e.g., *BRAF*) and the specific transcript used for annotation. |
| **Variant Identity** | `Variant ID (rsID)` | The dbSNP reference ID (e.g., *rs113488022*), crucial for cross-referencing with literature. |
| **Genomic Location** | `Chromosome`, `Position`, `Ref`, `Alt` | Exact genomic coordinates (GRCh38) and the specific base change. |
| **Mutation Detail** | `DNA Change`, `Protein Change` | HGVS notation describing the change at the DNA (c.) and Protein (p.) level. |
| **Clinical Significance** | `ClinVar`, `gnomAD AF` | **ClinVar:** Clinical interpretation (e.g., *Pathogenic, Benign*).<br>**gnomAD AF:** Global Allele Frequency to identify rare variants vs. common polymorphisms. |
| **Impact Prediction** | `Effect`, `Impact`, `Feature Type` | **Effect:** Type of mutation (e.g., *missense_variant*).<br>**Impact:** Predicted severity (*HIGH, MODERATE, LOW, MODIFIER*).<br>**Feature Type:** Affected feature (e.g., *transcript*). |
| **In Silico Scores** | `SIFT`, `PolyPhen` | Computational predictions of how the variant affects protein function. |

> **💡 Tip for Clinicians:** Start by filtering the **`Impact`** column for **HIGH** or **MODERATE**, and check the **`ClinVar`** column for known pathogenic variants.

---

## 🛠️ Troubleshooting

If you encounter issues while running **ThyroScope**, check the solutions below for common problems.


#### 1. `is not a valid Windows path` Error ⚠️

**Problem:** This usually happens when running the `docker run` command in **PowerShell** using CMD-style syntax (`%cd%`).

**Solution:** * **If using PowerShell:** Use `${PWD}` instead of `%cd%`.
  * *Example:* `-v "${PWD}\data:/data"`
* **If using Command Prompt (CMD):** Keep using `%cd%`.

#### 2. `SnpEff database not found!` Error 📂

**Problem:** The pipeline cannot find the `hg38` database files inside the container.

**Solution:** * **Check Folder Name:** Verify if your folder name is `snpEff_db` (all lowercase) or `SnpEff_db`.
* **Match Docker Mount:** The Docker mount command `-v .../snpEff_db:/pipeline/snpEff/data` must match your actual host folder name **exactly** (case-sensitive).
* **Verify Structure:** Ensure the `hg38` folder is located directly inside the `snpEff_db` folder on your host machine.

#### 3. `No input FASTQ files found` Error 🔍

**Problem:** The pipeline script cannot detect your raw sequencing data.

**Solution:** * **Check Naming Convention:** Ensure your files are named with the suffix `_1.fq.gz` and `_2.fq.gz`.
  * *Example:* `Patient01_1.fq.gz`, `Patient01_2.fq.gz`
* **Check Mounting:** Verify if the data folder is correctly mounted to the `/data` path in the Docker container.

#### 4. Memory/RAM Crash (GATK Errors) 💻

**Problem:** `GATK HaplotypeCaller` or `SnpEff` may crash if Docker is not allocated enough RAM.

**Solution 1: Docker Desktop Settings**
* Open **Docker Desktop Settings** > **Resources**.
* Increase the allocated **Memory** to at least **16GB** (8GB minimum).

**Solution 2: WSL2 Configuration (Advanced)**
If you are using Docker with the WSL2 backend, you can manually limit or allocate memory by editing the `.wslconfig` file:
1. Press `Win + R`, type `%UserProfile%`, and press **Enter**.
2. Create or edit a file named `.wslconfig` (ensure it has no `.txt` extension).
3. Paste the following configuration:
   ```ini
   [wsl2]
   memory=16GB    # Limits VM memory in WSL2
   processors=8   # Optional: Limits number of CPU cores
   ```
4. Restart WSL by running `wsl --shutdown` in PowerShell, then restart Docker Desktop.

> [!TIP]
> **Recommended Allocation:** It is generally recommended to set the WSL2 memory to **50–75%** of your total system RAM to maintain overall system stability.

---
## 📜 Citation & Contact

If you use **ThyroScope** in your research, please cite the following paper:

> 
### ✉️ Contact

For technical support, bug reports, or collaboration inquiries, please contact:

* **Developer:** Yun-seo Han ([yunseo21c@korea.ac.kr](mailto:yunseo21c@korea.ac.kr))
* **Lab:** H-Lee Lab, Department of Life Sciences
