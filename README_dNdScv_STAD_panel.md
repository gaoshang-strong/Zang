# 基于 dNdScv 的面板突变选择信号分析（STAD，panel + TCGA 同义稳背景）

> 本 README 以你提供的最终 R 代码为蓝本，系统性记录了从数据准备、坐标统一、与 TCGA 合并、构建 dNdScv 输入、质量控制、两轮 dNdScv 拟合，到基因层筛选与扩展分析（免疫治疗/幽门螺杆菌）的完整流程。**目标读者**：熟悉 MAF 与 R 的分子流行病学/肿瘤生信工程师。

---

## 目录
- [一、背景与总体思路](#一背景与总体思路)
- [二、数据与路径](#二数据与路径)
- [三、依赖与环境](#三依赖与环境)
- [四、TCGA-STAD 下载（hg38）](#四tcga-stad-下载hg38)
- [五、坐标统一：hg19 → hg38（liftOver）](#五坐标统一hg19--hg38liftover)
- [六、清洗并构建 dNdScv 5 列输入](#六清洗并构建-dndscv-5-列输入)
- [七、限定基因集合到面板基因](#七限定基因集合到面板基因)
- [八、第一次 dNdScv 与典型告警](#八第一次-dndscv-与典型告警)
- [九、参考碱基校正 + 相邻位点过滤 + 第二次拟合](#九参考碱基校正--相邻位点过滤--第二次拟合)
- [十、面板专属结果：借 TCGA 仅同义突变稳背景](#十面板专属结果借-tcga-仅同义突变稳背景)
- [十一、基因层结果抽取与筛选](#十一基因层结果抽取与筛选)
- [十二、可视化示例](#十二可视化示例)
- [十三、与免疫治疗/幽门螺杆菌的扩展分析](#十三与免疫治疗幽门螺杆菌的扩展分析)
- [十四、限制、常见坑与改进建议](#十四限制常见坑与改进建议)
- [十五、可复现性与记录](#十五可复现性与记录)

---

## 一、背景与总体思路
- **问题**：面板 MAF（~500+ 病人）**同义（Silent）很少**，直接用 dNdScv 拟合中性背景会不稳。  
- **策略**：从 **TCGA-STAD（hg38）** 下载大队列 MAF，**仅借用其同义突变**稳住背景，而 **非同义信号完全来自面板 cohort**。  
- **一致性**：两边坐标统一到 **hg38**；dNdScv 的 **RefCDS** 与默认 **协变量**原生支持 hg38。  
- **范围**：避免机会空间不一致，**限定到面板基因**（`gene_list = panel_genes_clean`）。

> 术语提醒：**Silent=同义突变**，**Nonsense=无义突变**。

---

## 二、数据与路径

- 面板 MAF（hg19）：`/ShangGaoAIProjects/Zang/MAF_data/gc_vcfs/somt.maf`  
- TCGA-STAD MAF（hg38）：用 `TCGAbiolinks::GDCprepare()` 获取  
- 输出前缀：`dNdscv/Sep072025_dNdscv`

---

## 三、依赖与环境

```r
library(data.table)
library(dplyr)
library(TCGAbiolinks)
library(GenomicRanges); library(GenomeInfoDb); library(rtracklayer)
library(dndscv)
library(BSgenome.Hsapiens.UCSC.hg38)
```

> 建议在末尾记录 `sessionInfo()`。

---

## 四、TCGA-STAD 下载（hg38）

```r
library(TCGAbiolinks)

query <- GDCquery(
  project = "TCGA-STAD",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
  access = "open"
)

# 首次需要取消注释：GDCdownload(query)
stad_maf <- GDCprepare(query)

write.table(
  stad_maf,
  "/ShangGaoAIProjects/Zang/MAF_data/TCGA_STAD.maf",
  sep = "\t",
  row.names = FALSE
)
```

---

## 五、坐标统一：hg19 → hg38（liftOver）

将面板 MAF 从 hg19 抬升到 hg38，保留一对一映射：

```r
library(data.table)
library(GenomicRanges)
library(rtracklayer)

maf <- fread("/ShangGaoAIProjects/Zang/MAF_data/gc_vcfs/somt.maf", data.table = FALSE)

chr <- ifelse(grepl("^chr", maf$Chromosome), maf$Chromosome, paste0("chr", maf$Chromosome))
gr19 <- GRanges(seqnames = chr, ranges = IRanges(start = maf$Start_Position, end = maf$End_Position))

chain <- rtracklayer::import("hg19ToHg38.over.chain", format = "chain")
lo <- liftOver(gr19, chain)

keep <- elementNROWS(lo) == 1
gr38 <- unlist(lo[keep])

maf38 <- maf[keep, ]
maf38$Chromosome     <- sub("^chr", "", as.character(seqnames(gr38)))
maf38$Start_Position <- start(gr38)
maf38$End_Position   <- end(gr38)
maf38$NCBI_Build     <- "GRCh38"

write.table(
  maf38,
  "/ShangGaoAIProjects/Zang/MAF_data/gc_vcfs/somt_hg38.maf",
  sep = "\t",
  row.names = FALSE
)
```

---

## 六、清洗并构建 dNdScv 5 列输入

合并（TCGA + panel）到统一字段，并导出 5 列：

```r
need <- c("Tumor_Sample_Barcode","Chromosome","Start_Position","End_Position",
          "Reference_Allele","Tumor_Seq_Allele2","Variant_Classification","Hugo_Symbol")

strip_chr <- function(x) sub("^chr", "", x, ignore.case = TRUE)

stad_maf$Chromosome <- strip_chr(as.character(stad_maf$Chromosome))
maf38$Chromosome    <- strip_chr(as.character(maf38$Chromosome))

stad_maf$Tumor_Sample_Barcode <- paste0("TCGA_",  stad_maf$Tumor_Sample_Barcode)
maf38$Tumor_Sample_Barcode    <- paste0("PANEL_", maf38$Tumor_Sample_Barcode)

combined_maf <- dplyr::bind_rows(stad_maf[, need], maf38[, need])
data.table::fwrite(combined_maf, "/ShangGaoAIProjects/Zang/MAF_data/gc_vcfs/combined_TCGA_PANEL_hg38_MAF.csv.gz")

keep_vc <- c("Silent","Missense_Mutation","Nonsense_Mutation","Splice_Site",
             "Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins")

mutations <- combined_maf %>%
  dplyr::filter(Variant_Classification %in% keep_vc) %>%
  dplyr::transmute(
    sampleID = as.character(Tumor_Sample_Barcode),
    chr      = as.character(Chromosome),
    pos      = as.integer(Start_Position),
    ref      = toupper(Reference_Allele),
    mut      = toupper(Tumor_Seq_Allele2)
  ) %>%
  dplyr::distinct(sampleID, chr, pos, ref, mut, .keep_all = TRUE)

data.table::fwrite(mutations, "/ShangGaoAIProjects/Zang/dNdscv/combined_for_dndscv_5col_input.csv.gz")
```

---

## 七、限定基因集合到面板基因

为避免机会空间不一致，限定分析到面板基因：

```r
panel_genes <- sort(unique(maf38$Hugo_Symbol))

tokens <- unique(unlist(strsplit(panel_genes, "[^A-Za-z0-9]+")))
tokens <- tokens[nchar(tokens) > 0]

drop_pat <- paste(c("^LOC\\d+", "^MIR\\d", "^LINC\\d", "^RP\\d", "^RNU\\d", "^MT-",
                    "-AS\\d*$", "AS\\d*$"), collapse="|")
tokens2 <- tokens[!grepl(drop_pat, tokens, ignore.case = TRUE)]

## 从包内 refcds_GRCh38_hg38.rda 获取 hg38 RefCDS 的合法蛋白编码基因名 ref_genes（不同版本对象名略有差异）
## panel_genes_clean <- intersect(tokens2, ref_genes)
```

> 注：不同版本的 `dndscv` 数据对象名可能不同，推荐 `load(system.file("data","refcds_GRCh38_hg38.rda", package="dndscv"))` 后查看对象结构再提取基因名。

---

## 八、第一次 dNdScv 与典型告警

直接在 `mutations` 跑 dNdScv 会出现：
- **Mutations observed in contiguous sites**：相邻位点（疑似把 MNV 拆成 2 个 SNV）；  
- **Same mutations in different sampleIDs**：跨样本复发（可能真实也可能技术/同患者多样本）；  
- **wrong reference base**：参考碱基错误，常见于**反向互补**。

---

## 九、参考碱基校正 + 相邻位点过滤 + 第二次拟合

**用 hg38 基因组批量校正 ref/mut 为正义链（仅 SNP）**：

```r
library(GenomicRanges); library(Biostrings); library(BSgenome.Hsapiens.UCSC.hg38)

gr <- GRanges(seqnames = paste0("chr", mutations$chr),
              ranges   = IRanges(start = mutations$pos, end = mutations$pos))
ref_true_vec <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr))

comp <- function(x) chartr("ACGT","TGCA", toupper(x))
is_snp    <- nchar(mutations$ref)==1 & nchar(mutations$mut)==1
needs_fix <- is_snp & (mutations$ref == comp(ref_true_vec))

mutations_fix <- mutations
mutations_fix$ref[needs_fix] <- ref_true_vec[needs_fix]
mutations_fix$mut[needs_fix] <- comp(mutations$mut[needs_fix])
```

**去掉相邻位点（疑似 MNV/复杂替换）**：

```r
library(dplyr)

mutations_fix <- mutations_fix %>%
  arrange(sampleID, chr, pos) %>%
  group_by(sampleID, chr) %>%
  mutate(
    prev_pos = dplyr::lag(pos),
    next_pos = dplyr::lead(pos),
    adj      = (pos - prev_pos == 1) | (next_pos - pos == 1)
  ) %>%
  ungroup() %>%
  filter(!adj) %>%
  select(sampleID, chr, pos, ref, mut)
```

**第二次拟合**：

```r
library(dndscv)

dnds2 <- dndscv(
  mutations_fix,
  refdb = "hg38",
  gene_list = panel_genes_clean,
  max_muts_per_gene_per_sample = Inf,
  max_coding_muts_per_sample    = Inf
)
```

---

## 十、面板专属结果：借 TCGA 仅同义突变稳背景

**思路**：TCGA 只保留 **Silent** 合并到 5 列输入，用于稳住背景；**非同义全部来自面板**。

```r
keep_vc_panel <- c("Silent","Missense_Mutation","Nonsense_Mutation","Splice_Site",
                   "Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins")
keep_vc_tcga  <- c("Silent")

panel_5 <- maf38 %>%
  dplyr::filter(Variant_Classification %in% keep_vc_panel) %>%
  dplyr::transmute(sampleID=paste0("PANEL_", Tumor_Sample_Barcode),
                   chr=Chromosome, pos=Start_Position,
                   ref=toupper(Reference_Allele), mut=toupper(Tumor_Seq_Allele2))

tcga_syn <- stad_maf %>%
  dplyr::filter(Variant_Classification %in% keep_vc_tcga) %>%
  dplyr::transmute(sampleID=paste0("TCGA_SYN_", Tumor_Sample_Barcode),
                   chr=Chromosome, pos=Start_Position,
                   ref=toupper(Reference_Allele), mut=toupper(Tumor_Seq_Allele2))

mutations_hybrid <- dplyr::bind_rows(panel_5, tcga_syn) %>%
  dplyr::distinct(sampleID, chr, pos, ref, mut)

dnds_panel_only <- dndscv(
  mutations_hybrid,
  refdb = "hg38",
  gene_list = panel_genes_clean,
  max_muts_per_gene_per_sample = Inf,
  max_coding_muts_per_sample   = Inf
)
```

**本次（示例）全局 dN/dS：**
```
wmis = 1.6769 (95% CI: 1.5899–1.7685)
wnon = 2.0137 (95% CI: 1.8008–2.2516)
wspl = 1.9398 (95% CI: 1.6533–2.2760)
wtru = 1.9905 (95% CI: 1.8074–2.1922)
wall = 1.6997 (95% CI: 1.6125–1.7917)
```
> 结论：面板 cohort 展现明确正向选择信号（尤其截短类）。

---

## 十一、基因层结果抽取与筛选

你的 `sel_cv` 列结构：  
- 观测：`n_syn, n_mis, n_non, n_spl, n_ind`  
- 基因级 dN/dS：`wmis_cv, wnon_cv, wspl_cv, wind_cv`  
- p/q 值：`pmis_cv, ptrunc_cv, pallsubs_cv, pind_cv, pglobal_cv` 与对应 q 值（`q*`）

**抽取 + 三类筛选：**

```r
sig_genes <- dnds_panel_only$sel_cv %>%
  arrange(qglobal_cv) %>%
  select(gene_name, qglobal_cv,
         wmis_cv, wnon_cv, wspl_cv, wind_cv,
         n_mis, n_non, n_spl, n_ind,
         qmis_cv, qtrunc_cv, qallsubs_cv, qind_cv)

res <- dnds_panel_only$sel_cv %>%
  mutate(
    n_nonsyn = n_mis + n_non + n_spl,
    trunc_n  = n_non + n_spl + n_ind,
    any_w_gt1 = (wmis_cv>1) | (wnon_cv>1) | (wspl_cv>1) | (wind_cv>1)
  )

# 1) 核心驱动集（均衡）
core_drivers <- res %>%
  filter(qglobal_cv <= 0.10, any_w_gt1, n_nonsyn >= 3) %>%
  arrange(qglobal_cv) %>%
  select(gene_name, qglobal_cv, wmis_cv, wnon_cv, wspl_cv, wind_cv,
         n_mis, n_non, n_spl, n_ind)

# 2) TSG 候选（失活富集）
tsg_candidates <- res %>%
  filter(qtrunc_cv <= 0.10, (wnon_cv>1 | wspl_cv>1 | wind_cv>1), trunc_n >= 2) %>%
  arrange(qtrunc_cv) %>%
  select(gene_name, qtrunc_cv, wnon_cv, wspl_cv, wind_cv,
         n_non, n_spl, n_ind, qglobal_cv)

# 3) Oncogene 候选（错义富集）
onc_candidates <- res %>%
  filter(qglobal_cv <= 0.10, wmis_cv > 1, n_mis >= 2) %>%
  arrange(qglobal_cv) %>%
  select(gene_name, qglobal_cv, wmis_cv, n_mis, wnon_cv, wspl_cv)
```

---

## 十二、可视化示例

**全局 dN/dS 误差棒图、基因散点图、Top 基因条形图**（可选 Oncoplot）：

```r
# 需: ggplot2 dplyr scales maftools data.table
# 图1：全局 dN/dS
g <- dnds_panel_only$globaldnds %>%
  dplyr::mutate(class = factor(name, levels = c("wmis","wnon","wspl","wtru","wall"),
                        labels = c("Missense","Nonsense","Splice","Truncating","All nonsyn")))

ggplot2::ggplot(g, ggplot2::aes(x = class, y = mle)) +
  ggplot2::geom_hline(yintercept = 1, linetype = 2) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = cilow, ymax = cihigh), width = 0.15) +
  ggplot2::labs(x = NULL, y = "dN/dS (MLE ± 95% CI)", title = "Global dN/dS (panel-only, syn from TCGA)") +
  ggplot2::theme_bw(base_size = 12)

# 准备基因层表与散点/条形图略（参见你现有脚本或聊天记录中的模板）
```

> Oncoplot 需要标准 MAF 列；若已有 `panel_maf_hg38`，可用 `maftools::read.maf()` 后 `oncoplot()`。

---

## 十三、与免疫治疗/幽门螺杆菌的扩展分析

### A. 免疫治疗（IO）
- **TMB（面板近似）**：`TMB ≈ 非同义数 / 捕获大小(Mb)`；联合 `wall/wtru` 分层疗效（ORR/PFS/OS）。  
- **MSI/MMR**：MLH1/MSH2/MSH6/PMS2 截短显著（`qtrunc_cv`）+ 工具（MSIsensor 等）；MSI-H 往往 `wall` 更高。  
- **抗原呈递/免疫逃逸基因**：B2M、HLA-A/B/C、JAK1/2、STAT1、TAP1/2、CIITA 在 `qtrunc_cv/qglobal_cv` 的显著性。  
- **免疫编辑**：可呈递肽段 vs 不可呈递区域的 dN/dS（需 HLA 分型与 epitope 预测）。

### B. 幽门螺杆菌（H. pylori）
- **HP+ vs HP− 分组**：比较 `wall/wtru` 与基因层 q 值分布；可用 bootstrap/置换。  
- **APOBEC 近似特征**：在面板上粗评 TCW→T/K 富集，比较 HP 分组。

### C. 通用扩展
- **通路富集**（Reactome/KEGG）：PI3K、WNT、TGF-β、Cell cycle、DNA repair、JAK/STAT。  
- **共突变/互斥**：如 TP53 vs ARID1A，PIK3CA vs RHOA 等。  
- **临床关联**：分期分级、年龄性别、治疗线数、吸烟饮酒，与 dN/dS 或特定基因显著性的关系（多变量回归）。

---

## 十四、限制、常见坑与改进建议

- **面板机会空间**：`gene_list` 按**基因**限制，未按捕获区间精确限制 CDS——实务常用；若需严格，需自建区域化 RefCDS。  
- **版本差异**：`sel_cv` 列名（是否带 `_cv`）随版本不同；建议先 `colnames()`。  
- **分隔符**：写表请用 `sep = "\t"`，避免误写成 `"/t"`。  
- **相邻位点**：MNV 拆分会抬高假阳性；已做相邻位点过滤。  
- **参考碱基错误**：反向互补较常见；已用 hg38 批量校正（仅 SNP）。  
- **重复突变**：同一患者多样本需按患者去重；不同患者复发可保留。  
- **子模型复杂度**：面板同义稀疏时，可考虑 `sm="12r_3w"` 以稳估计。  
- **LiftOver**：个别 indel/重复区映射不稳，必要时剔除。

---

## 十五、可复现性与记录

建议保存关键对象与环境：

```r
saveRDS(list(mutations_fix = mutations_fix,
             panel_genes_clean = panel_genes_clean,
             dnds_panel_only = dnds_panel_only),
        file = "dNdscv/objects_for_reproducibility.rds")

sessionInfo()
```

---

### 致谢与备注
- 本流程遵循 dNdScv 的统计建模：以 **三核苷酸子模型 + 过度离散回归** 估计中性背景率 λ，并在 **同义/非同义类别**上做比率检验（`w = dN/dS`）。  
- “借 TCGA 同义稳背景”是针对 **面板同义稀疏** 的实用办法，**不引入 TCGA 非同义**以保持 cohort 效应可解释性。
