import os
import subprocess
import sys
import glob
import pandas as pd

# ==========================================
# 1. 환경 설정 (Configuration) - Docker 내부 경로 기준
# ==========================================

# 환자 ID (나중에 입력 파일명에서 자동 추출하도록 변경 가능)
SAMPLE_ID = "Patient_Test"

# 파일 경로 설정 (도커 마운트 경로)
DATA_DIR = "/data"
REF_DIR = "/data/ref"
PIPELINE_DIR = "/pipeline"

# 참조 유전체 및 타겟 파일
# 주의: 사용자가 data/ref 폴더에 해당 파일들을 넣어둬야 함
REF_GENOME = f"{REF_DIR}/Homo_sapiens_assembly38.fasta"
TARGET_BED = f"{PIPELINE_DIR}/targets.bed"

# 툴 경로 설정
# Trimmomatic (Docker 설치 시 기본 경로)
TRIMMOMATIC_JAR = "/usr/share/java/trimmomatic.jar"
# 어댑터 파일 (Illumina 장비 표준)
ADAPTERS = "/usr/share/trimmomatic/adapters/TruSeq3-PE.fa"

# 기타 툴 명령어
BWA = "bwa"
GATK = "gatk"
SAMTOOLS = "samtools"
ANNOVAR_CMD = "table_annovar.pl" # ANNOVAR는 라이선스 문제로 별도 설치/설정 필요

# ==========================================
# 2. 실행 함수 (Wrapper Function)
# ==========================================
def run_command(cmd, step_name):
    print(f"\n[INFO] Starting Step: {step_name}")
    print(f"Command: {cmd}")
    try:
        subprocess.check_call(cmd, shell=True)
        print(f"[SUCCESS] {step_name} completed.\n")
    except subprocess.CalledProcessError:
        print(f"[ERROR] {step_name} failed.")
        # 파이프라인 중단 (에러 발생 시)
        sys.exit(1)

# ==========================================
# 3. 메인 파이프라인 (Main Workflow)
# ==========================================
def main():
    print(">>> Starting Thyroid/Parathyroid WES Pipeline (STAR Protocol Compliant) <<<")
    
    # 1) 입력 FASTQ 파일 찾기
    fastqs = sorted(glob.glob(f"{DATA_DIR}/*.fastq.gz"))
    if len(fastqs) < 2:
        print(f"[ERROR] No FASTQ files found in {DATA_DIR}. Please check your data.")
        sys.exit(1)
        
    r1_raw = fastqs[0]
    r2_raw = fastqs[1]
    print(f"Found Input Files: {r1_raw}, {r2_raw}")

    # ---------------------------------------------------------
    # Step 0: Trimming (Trimmomatic) - STAR Protocol 필수 단계
    # ---------------------------------------------------------
    # 어댑터 제거 및 퀄리티 낮은 염기서열 절단
    r1_paired = f"{DATA_DIR}/{SAMPLE_ID}_R1_paired.fq.gz"
    r1_unpaired = f"{DATA_DIR}/{SAMPLE_ID}_R1_unpaired.fq.gz"
    r2_paired = f"{DATA_DIR}/{SAMPLE_ID}_R2_paired.fq.gz"
    r2_unpaired = f"{DATA_DIR}/{SAMPLE_ID}_R2_unpaired.fq.gz"

    cmd_trim = (
        f"java -jar {TRIMMOMATIC_JAR} PE -threads 4 "
        f"{r1_raw} {r2_raw} "
        f"{r1_paired} {r1_unpaired} "
        f"{r2_paired} {r2_unpaired} "
        f"ILLUMINACLIP:{ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
    )
    run_command(cmd_trim, "0. Trimming (Trimmomatic)")

    # ---------------------------------------------------------
    # Step 1: Alignment (BWA-MEM)
    # ---------------------------------------------------------
    # ★ 중요: Trimming이 완료된 파일(r1_paired)을 입력으로 사용
    bam_output = f"{DATA_DIR}/{SAMPLE_ID}.bam"
    cmd_bwa = (
        f"{BWA} mem -t 4 "
        f"-R '@RG\\tID:{SAMPLE_ID}\\tSM:{SAMPLE_ID}\\tPL:ILLUMINA' "
        f"{REF_GENOME} {r1_paired} {r2_paired} | "
        f"{SAMTOOLS} view -bS - > {bam_output}"
    )
    run_command(cmd_bwa, "1. Alignment (BWA)")

    # ---------------------------------------------------------
    # Step 2: Sorting & Indexing
    # ---------------------------------------------------------
    sorted_bam = f"{DATA_DIR}/{SAMPLE_ID}.sorted.bam"
    cmd_sort = f"{SAMTOOLS} sort -o {sorted_bam} {bam_output}"
    run_command(cmd_sort, "2. Sorting (Samtools)")
    
    cmd_index = f"{SAMTOOLS} index {sorted_bam}"
    run_command(cmd_index, "2-1. Indexing")

    # ---------------------------------------------------------
    # Step 3: Mark Duplicates (GATK)
    # ---------------------------------------------------------
    dedup_bam = f"{DATA_DIR}/{SAMPLE_ID}.dedup.bam"
    metrics_file = f"{DATA_DIR}/{SAMPLE_ID}.metrics.txt"
    
    cmd_dedup = (
        f"{GATK} MarkDuplicates "
        f"-I {sorted_bam} -O {dedup_bam} -M {metrics_file}"
    )
    run_command(cmd_dedup, "3. Mark Duplicates (GATK)")
    
    # 인덱싱 (GATK 다음 단계를 위해 필요)
    run_command(f"{SAMTOOLS} index {dedup_bam}", "3-1. Dedup Indexing")

    # ---------------------------------------------------------
    # Step 4: Variant Calling (Targeted)
    # ---------------------------------------------------------
    vcf_output = f"{DATA_DIR}/{SAMPLE_ID}.vcf"
    
    cmd_hc = (
        f"{GATK} HaplotypeCaller "
        f"-R {REF_GENOME} "
        f"-I {dedup_bam} "
        f"-O {vcf_output} "
        f"-L {TARGET_BED} "       # 선생님의 타겟 유전자 리스트 사용
        f"--interval-padding 100" # 유전자 주변 100bp까지 포함
    )
    run_command(cmd_hc, "4. Variant Calling (GATK HaplotypeCaller)")

    # ---------------------------------------------------------
    # Step 5: Annotation (ANNOVAR) - 예시 코드
    # ---------------------------------------------------------
    # 주의: ANNOVAR DB가 설치되어 있어야 작동합니다. (Docker 이미지에 DB 포함 권장)
    # 만약 ANNOVAR가 없다면 이 단계에서 에러가 날 수 있으므로 try-except 처리하거나 주석 처리
    anno_output = f"{DATA_DIR}/{SAMPLE_ID}.anno"
    try:
        cmd_anno = (
            f"{ANNOVAR_CMD} {vcf_output} {DATA_DIR}/humandb/ -buildver hg38 "
            f"-out {anno_output} -remove -protocol refGene,cytoBand,gnomad211_exome,clinvar_20221231 "
            f"-operation g,r,f,f -nastring . -vcfinput"
        )
        # run_command(cmd_anno, "5. Annotation (ANNOVAR)") 
        # ▲ 실제 실행하려면 위 주석(#)을 해제하고 ANNOVAR DB를 연결해야 합니다.
        print("[INFO] Step 5 (Annotation) skipped in this demo script. (Requires ANNOVAR DB)")
        
    except Exception as e:
        print(f"[WARNING] Annotation step failed or skipped: {e}")

    # ---------------------------------------------------------
    # Step 6: Final Reporting (Excel)
    # ---------------------------------------------------------
    print(">>> Generating Clinical Report... <<<")
    
    # ANNOVAR 결과 파일이 있다고 가정하고 엑셀 변환 (없으면 VCF만 변환)
    # 여기서는 간단히 VCF를 파싱하여 엑셀로 만드는 예시를 보여드립니다.
    
    try:
        # VCF 파일 읽기 (주석 줄 제외)
        vcf_cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
        df = pd.read_csv(vcf_output, sep="\t", comment='#', names=vcf_cols)
        
        # 엑셀 저장
        excel_output = f"{DATA_DIR}/{SAMPLE_ID}_Clinical_Report.xlsx"
        df.to_excel(excel_output, index=False)
        print(f"[SUCCESS] Excel Report Saved: {excel_output}")
        
    except Exception as e:
        print(f"[WARNING] Could not generate Excel report: {e}")

    print("\n>>> All Pipeline Steps Completed Successfully! <<<")

if __name__ == "__main__":
    main()
