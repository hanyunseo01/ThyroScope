# 🔭 ThyroScope: Clinical-Grade WES Pipeline for Thyroid and Parathyroid Disorders

[![Docker Image Version](https://img.shields.io/docker/v/hanyunseo01/thyroid_pipeline/v23?color=blue&label=Docker%20Image)](https://hub.docker.com/r/hanyunseo01/thyroid_pipeline)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20macOS%20%7C%20Linux-lightgrey)](https://www.docker.com/)

**ThyroScope** is a streamlined, containerized bioinformatics pipeline designed to "scope out" and analyze thyroid-specific variants from **Whole Exome Sequencing (WES)** data.

By leveraging a **"Lightweight Containerization Strategy"** and a **Virtual Panel approach**, ThyroScope allows clinicians to bypass complex command-line interfaces. It automates the entire workflow—from raw FASTQ to a clinical-grade Excel report—with a single click, focusing specifically on 25 high-priority genes associated with thyroid and parathyroid disorders.
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
