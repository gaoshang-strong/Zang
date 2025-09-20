# DNA-seq Somatic Variant Calling Pipeline

本仓库包含一个基于 **GATK4 Mutect2** 的 DNA-seq 分析流程，用于从原始 FASTQ 文件到 **体细胞突变检测**（somatic variant calling）。  
整个流程通过 Bash 脚本实现，依赖一系列主流工具（BWA、samtools、GATK 等），并按照标准顺序完成比对、去重、BQSR 及变异检测。  

---

## 📦 环境与依赖

在运行前，请确保安装以下工具并加入 `PATH`：

- [BWA](http://bio-bwa.sourceforge.net/) ≥ 0.7.17  
- [samtools](http://www.htslib.org/) ≥ 1.10  
- [GATK4](https://gatk.broadinstitute.org/) ≥ 4.2  
- GNU coreutils（bash, mkdir, cp 等）

### 参考文件与资源
- 人类参考基因组（`REF_FASTA`）
- gnomAD Germline resource（`GERMLINE_RESOURCE`）
- Panel of Normals（可选，`PON_VCF`）
- Known sites VCF（用于 BQSR）

所有路径和资源文件通过 `config.sh` 配置。

---

## ⚙️ 配置文件

### `config.sh`
该文件定义了整个流程运行所需的全局变量，包括：

- **输入数据**
  - `FASTQ_BASE_DIR`：FASTQ 文件目录
  - `SAMPLE_SHEET`：样本清单（CSV，包含 Run title, File1, File2, type, patient）

- **输出目录**
  - `WORKDIR`：分析结果根目录，自动创建子目录：
    - `00_pairs`：样本配对信息
    - `01_align`：比对结果
    - `02_markdup`：去重结果
    - `03_bqsr`：BQSR 结果
    - `04_mutect2`：变异检测结果
    - `05_filter`：后续过滤结果
    - `logs`：日志文件

- **参考资源**
  - `REF_FASTA`
  - `GERMLINE_RESOURCE`
  - `PON_VCF`
  - `KNOWN_SITES_1/2/3`

- **工具与线程数**
  - `THREADS`
  - `GATK`
  - `BWA`
  - `SAMTOOLS`

---

## 🚀 运行方式

整个流程由 `run_pipeline.sh` 管理，按顺序执行各个步骤并记录日志：

```bash
bash run_pipeline.sh
```

日志存放在 `analysis_mutect2/logs/`。

---

## 🔬 流程步骤

### **Step 1: 比对与排序**
脚本：`01_align_and_sort.sh`

- 使用 **BWA-MEM** 将 FASTQ 比对到参考基因组  
- 添加 read group 信息  
- 用 samtools sort 排序并生成 BAM 文件  
- 输出目录：`patient/Tumor|Normal/01_align/`  
- 产物：
  - `*.sorted.bam`
  - `*.sorted.bam.bai`

---

### **Step 2: 去重**
脚本：`02_markduplicates.sh`

- 使用 **GATK MarkDuplicates** 去除 PCR duplicates  
- 自动生成 metrics 文件和 BAM index  
- 输出目录：`patient/Tumor|Normal/02_markdup/`  
- 产物：
  - `*.markdup.bam`
  - `*.metrics.txt`

---

### **Step 3: BQSR（碱基质量分数重校正）**
脚本：`03_bqsr.sh`

- 使用 **BaseRecalibrator** 和 **ApplyBQSR** 对 BAM 文件进行质量校正  
- 输入：已去重 BAM + known sites VCF  
- 输出目录：`patient/Tumor|Normal/03_bqsr/`  
- 产物：
  - `*.bqsr.bam`
  - `*.bqsr.bai`
  - `*.recal.table`

如果未设置 known sites，会跳过 BQSR 并直接复制 BAM。

---

### **Step 4: 体细胞突变检测**
脚本：`04_mutect2.sh`

- 使用 **GATK Mutect2** 进行体细胞突变检测  
- 输入：
  - `Tumor/03_bqsr/*.bqsr.bam`
  - `Normal/03_bqsr/*.bqsr.bam`
- 参数：
  - `--germline-resource`：gnomAD
  - `--panel-of-normals`：可选
  - `--native-pair-hmm-threads`：线程数
- 输出目录：`patient/Tumor/04_mutect2/`  
- 产物：
  - `*.unfiltered.vcf.gz`
  - `*.f1r2.tar.gz`

---

## 📊 输出结果

最终的主要结果文件为：
- **VCF**（变异检测结果，未过滤）：`patient/Tumor/04_mutect2/patient.unfiltered.vcf.gz`  
- 后续可通过 GATK `FilterMutectCalls` 进行过滤，并进一步注释。

---

## 📖 日志管理

每一步都会在 `analysis_mutect2/logs/` 下生成独立的 log 文件，方便调试和溯源：
- `01_align_and_sort.log`
- `02_markduplicates.log`
- `03_bqsr.log`
- `04_mutect2.log`
- `pipeline.log`（整体流程记录）

---

## ✅ 总结

本 pipeline 实现了一个 **标准的 DNA-seq 体细胞突变检测流程**，包含：
1. FASTQ → 比对 → BAM  
2. 去重  
3. BQSR  
4. Mutect2 检测  

用户可根据实际需求修改 `config.sh` 中的路径与资源文件。
