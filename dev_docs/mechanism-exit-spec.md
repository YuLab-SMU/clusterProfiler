# clusterProfiler 机制出口：`evidence()` 分期实现规范

版本：v2 · 2026-10-25 · 状态：已过独立审核，待采纳

> **v2 修订记录**：独立审核（含实测复现）发现 3 处硬错误、若干建议。本版已修正，并在相应位置标记「v2 改」。v1 的两处错误出自我自己的判断失误（GSON 契约、`pscore` 通道语义），不是笔误。

---

## 0. 这份文档要解决什么

把 clusterProfiler 从「富集到通路」推进到「多源证据裁决出的机制候选」，并让 LLM 只在被裁定的范围内说话。

**不在本规范范围内**：新的统计检验、新的富集引擎。ORA/GSEA 都已经在包里，这里只加「证据从哪来、怎么定级、怎么并线、怎么交给模型」。

**用户只需要记住两个函数**：`evidence()` 整理证据，`interpret()` 说人话。

---

## 1. 范围与分期

### 1.1 上游基因集的完整清单

| 来源 | 集合 | 上游实体 | 方向性 | 数据 | 默认等级 |
|---|---|---|---|---|---|
| `tf` | MSigDB C3 `TFT` | TF 蛋白 | 有 | GMT 2.83 MB / 1114 集 | B（motif-only 降 C） |
| `mirna` | MSigDB C3 `MIR` | miRNA | 负调控 | GMT 2.78 MB / 2598 集 | C |
| `perturbagen` | MSigDB C2 `CGP` | 化学物、基因扰动 | 有 | GMT 2.97 MB / 3555 集 | B（与实测扰动一致可升 A） |
| `cellstate` | MSigDB C8 | 细胞类型 | 无 | GMT 1.13 MB / 866 集 | C |
| `tissue` | MSigDB **C4** `CGN` | 组织（来源签名） | 无 | GMT 0.29 MB / 427 集 | C |
| `complex` | CORUM / hu.MAP | 蛋白复合体 | 无 | **无 GMT，需自建** | B |

**文件名实测（2026-10-25）**：

```text
c3.tft.v2026.2.Hs.symbols.gmt   200  2821791
c3.mir.v2026.2.Hs.symbols.gmt   200  2779338
c2.cgp.v2026.2.Hs.symbols.gmt   200  2966321
c8.all.v2026.2.Hs.symbols.gmt   200  1129789
c4.cgn.v2026.2.Hs.symbols.gmt   200   291423
c2.cgn.v2026.2.Hs.symbols.gmt   404      <-- C2 下没有 CGN
```

> **v2 改**：v1 曾写「书里 `regulator-enrichment.qmd` 写的是 C2 CGN，需要同步修正」。**这句是错的**——书里本来就写 C4，无需修正。只保留「C2 下没有 CGN」这一实测事实。

### 1.2 分期

| 期 | 内容 | 交付判据 |
|---|---|---|
| **P1** | `mechanism_sets()` + `evidence()` 骨架 + `xRef` + `diagnostics` | §7.1 |
| **P2** | `consensus()` + `writable` + 网络证据定级 | §7.2 |
| **P3** | `interpret()` 消费 bundle | §7.3 |
| **P4** | MCP 通道：`perturbagen` / `coexpression` / `genetic` | §7.4 |
| **P5** | `cellstate` / `tissue` / `complex` | 等级恒 C/B 的回归测试 |

---

## 2. 关键设计决策

### D1：`GSON` 只能走 `gson =`，不能塞进 `TERM2GENE =`

`enricher()` 与 `GSEA()` 的签名里都有 `gson = NULL`。2026-10-25 实测：

```r
enricher(gene, TERM2GENE = G)
#> ERROR: attempt to set 'colnames' on an object with less than two dimensions

enricher(gene, gson = G)
#> OK；且 @organism = "Homo sapiens", @keytype = "SYMBOL", @ontology = "MSigDB C3:TFT"
```

**只有 `gson =` 会把 `species`/`keytype`/`gsname` 带进结果的三个槽**；`TERM2GENE =`（数据框路径）会把它们写成 `unknown`。`$regulators` 的命名空间判断依赖这三个槽，所以这是契约不是风格。

因此 `mechanism_sets()` 返回 `GSON`，由 `enricher(gson = ...)` / `GSEA(gson = ...)` 消费。**写进测试 T0。**

### D2：各类分别做 BH 校正，不合池

TF 1114 集、miRNA 2598 集、CGP 3555 集，靶定义方式与集合大小分布不同。合池做 BH 会稀释或抬高各自的显著性：

```text
$tf     在 1114 个集合内校正
$mirna  在 2598 个集合内校正
$regulators  仅合并结果，不合并校正（故列名 padj_within_source）
```

### D3：等级是「规则表」，不是「每个 source 一行」

`$xRef` 是规则表（每行一条规则 + 匹配条件），按序匹配，第一条命中即用：

| 序 | source | match | default_level |
|---|---|---|---|
| 1 | `tf` | term 匹配 `*_UNKNOWN` 或无因子名 | C |
| 2 | `tf` | term 匹配 `_TARGET_GENES$` | B |
| 3 | `tf` | 其它（legacy motif） | B |
| 4 | `mirna` | 全部 | C |
| 5 | `perturbagen` | 全部 | B |
| 6 | `cellstate` | 全部 | C |
| 7 | `tissue` | 全部 | C |
| 8 | `complex` | 全部 | B |

**顺序有意义**：第 1 条必须在前，否则 `AACTTT_UNKNOWN` 会因不含 `_TARGET_GENES` 落到第 3 条拿 B。

### D4：缓存——P1 必须一并修 `.onAttach`

`R/zzz.R:8` 是**无条件赋值** `options(clusterProfiler_cache_dir = tempdir())`：会覆盖用户在 `.Rprofile` 里设的值，且 `tempdir()` 每个会话都不同 → **`mechanism_sets()` 缓存永不命中，每次重下 2.8 MB**。

v2 决定 **P1 里修**：

```r
if (is.null(getOption("clusterProfiler_cache_dir"))) {
    options(clusterProfiler_cache_dir = yulab.utils::user_dir("clusterProfiler"))
}
```

影响面已核对：全包仅 `R/ppi.R:196` 与本方案使用该 option。

### D5：GMT 版本是字面量

`release` 是函数参数，默认常量 `"2026.2.Hs"`，写进 GSON `version` 槽。`accessed_date` 记首次下载日期（命中不重写）。不用 `Sys.Date()` 当版本号。

### D6：输入形态收敛为「基因列表」

`enrichResult` 展开 `geneID`、data.frame 取 `design$gene_col`、或直接给向量。只支持一种内部表示。

### D7：不做自动降级到「猜」

证据缺席时返回空表 + `diagnostics` 写明原因，不返回「最像的」替代结果。

---

## 3. 用户入口

| 用户 | 带什么 | 真正想问 | 调用 |
|---|---|---|---|
| A 差异表达 | DE 表 / 显著基因 | 谁在管这些基因 | `evidence(genes, design = ...)` |
| B 扰动实验 | 敲低/过表达 DE | 驱动者是谁 | `evidence(genes, design = list(treatment = "knockdown", gene = "HIF1A"))` |
| C 已有结果 | `enrichResult` | 结论靠不靠得住 | `evidence(enrich_result)` |
| D 写文章 | 任一 | 这句话能说到什么程度 | `evidence(...)$writable` |

### 3.1 `design` 的完整字段（v2 补全）

v1 的参数表引用了 `design$gene_col` 与 `design$coverage_floor`，但这里没定义。现补全：

```r
design <- list(
  contrast      = "HIF1A knockdown vs control",  # 进 prompt 与 provenance
  system        = "HeLa",
  species       = "Homo sapiens",                # 决定集合物种
  treatment     = "knockdown",                   # none|knockdown|knockout|overexpression|drug
  gene          = "HIF1A",                       # 被扰动基因（可选）
  direction     = NULL,                          # up|down|NULL
  gene_col      = "gene",                        # genes 为 data.frame 时取哪列
  coverage_floor = 0.30                          # 低于此值判 untested 而非 absent
)
```

`treatment` 决定 A 级证据是否可及，**不让模型从自由文本里猜**。

`coverage_floor` 用 **`n_annotated / n_input`** 计算（不是 `n_sig_covered / n_input`），见 §4.4。

---

## 4. API 规范

### 4.1 `mechanism_sets()`（P1）

```r
mechanism_sets(source = c("tf", "mirna", "perturbagen", "cellstate", "tissue"),
               species = "Homo sapiens",
               release = "2026.2.Hs")
```

| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `source` | character(1) | 必填 | 见 §1.1 |
| `species` | character(1) | `"Homo sapiens"` | 目前只有 `*.Hs.*` |
| `release` | character(1) | `"2026.2.Hs"` | 字面常量（D5） |

**返回**：`GSON`

| 槽 | 值 |
|---|---|
| `gsid2gene` | `data.frame(gsid, gene)`，`gene` 为 SYMBOL |
| `gsid2name` | `data.frame(gsid, name)` |
| `species` | 参数值 |
| `gsname` | `"MSigDB C3:TFT"` / `"C3:MIR"` / `"C2:CGP"` / `"C8"` / `"C4:CGN"` |
| `keytype` | `"SYMBOL"` |
| `version` | `release` |
| `accessed_date` | 首次下载日期 |

**实现细节**

1. 文件名映射按 §1.1 实测表（`c8.all`、`c4.cgn`、`c2.cgp` 不能按 subcollection 名推）。
2. URL：`https://data.broadinstitute.org/gsea-msigdb/msigdb/release/<release>/<file>`
3. 缓存：`file.path(getOption("clusterProfiler_cache_dir"), file)` + `.meta`（记 release 与下载日期）。命中判据：文件存在 **且** `.meta` release 一致；不一致则重下并 **warn**。
4. 解析：`name<TAB>url<TAB>gene1<TAB>...`；**必须处理 CRLF**。
5. `gsid` 保持原名。归一到调控者名是 `evidence()` 的职责。
6. `unique()` 作用于 `gsid2gene`。

**错误**：`source` 非法 / 下载非零 / 下到 HTML / 解析为空 → `rlang::abort()`，消息带 URL 与目标路径。

---

### 4.2 `evidence()`（P1 骨架）

**v2 改**：删掉四个 logical 开关，统一走 `upstream`；`pvalueCutoff`/`qvalueCutoff` 改名并降级为打印参数；`seed` 推迟到 P2（P1 无随机过程，不生效的参数比没有更糟）。

```r
evidence(genes,
         design        = list(),
         upstream      = c("tf", "mirna"),
         universe      = NULL,
         minGSSize     = 5,
         maxGSSize     = Inf,
         print_pvalue  = 1,
         print_qvalue  = 1)
```

| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `genes` | character / `enrichResult` / data.frame | 必填 | data.frame 取 `design$gene_col` |
| `design` | list | `list()` | §3.1 |
| `upstream` | character | `c("tf","mirna")` | 唯一来源入口 |
| `universe` | character / NULL | `NULL` | `NULL` 时用 `org.*.eg.db` 全部 SYMBOL |
| `minGSSize` | numeric | **5** | 与书一致（见 §7.1 的 1023/1030 差异） |
| `maxGSSize` | numeric | `Inf` | TF 集最大 2000 靶，默认 500 会截掉最重要的行 |
| `print_pvalue` / `print_qvalue` | numeric | `1` | **只影响打印，不影响检验与裁决** |

> **v2 改**：v1 的 `pvalueCutoff = 0.05` 会让人以为它过滤了检验。裁决必须用全表，否则会把 `single_source` 误判成 `absent`。

**返回** `evidenceBundle`（`class = c("evidenceBundle","list")`）

| 元素 | P1 产出 | 说明 |
|---|---|---|
| `$input` | ✅ | genes / universe / design / namespace |
| `$sets` | ✅ | 每 source：gsname、version、n_sets（来自 GSON 槽，D1） |
| `$tf`、`$mirna`… | ✅ | 每类一份 `enrichResult`，各自类内 BH |
| `$regulators` | ✅ | 裁决长表（§4.3） |
| `$xRef` | ✅ | 规则表（D3） |
| `$diagnostics` | ✅ | 覆盖率四列（§4.4） |
| `$provenance` | ✅ | release / URL / 时间 / 包版本 |
| `$writable` | P2 | §4.6 |
| `$perturbagen` / `$coexpression` / `$genetic` | P4 | 与 `upstream` 枚举同名 |
| `$ppi` | P2 | 网络证据（不是证据源，故独立命名） |

> **v2 改**：v1 写 `$perturbation`，与 `upstream = "perturbagen"` 不一致。统一用 source 枚举值；`$network` → `$ppi`。

---

### 4.3 `$regulators` 的列（v2 改）

列名对齐既有约定，不再发明第三套命名：

```text
regulator           归一后的名字（FOXO3、AP1）
source              tf / mirna / perturbagen / cellstate / tissue / complex
term                原始集合 ID
named               TRUE / FALSE（motif-only 且无因子名者 FALSE）
level               A / B / C（按 $xRef 规则匹配）
padj_within_source  该类内 BH 后的值（名字本身是 D2 的教学点）
setSize             集合大小（对齐既有 setSize）
Count               与输入列表的交集（对齐既有 Count）
```

> `padj` → `padj_within_source`，`n_members` → `setSize`，`n_selected` → `Count`。

**冻结策略（v2 改）**：**已有列名冻结，允许追加**。`direction`（P2）、`kind`（P2+）、分通道计数（P4）以新增列出现，不改名。v1 说「一次性定死」对 P4 太紧。

**调控者归一化**（书里已验证，逐字搬运）：

```r
mechanism_regulator <- function(term) {
    tf <- sub("_TARGET_GENES$", "", term)
    tf <- sub("_[CQ][0-9A-Z]*_[0-9]+$", "", tf)   # SRF_C, HNF3_Q6, AP1_Q4_01
    tf <- sub("_[CQ][0-9]+$", "", tf)             # AP1_Q4
    sub("^[ACGTUN]+_", "", tf)                    # TGANTCA_AP1 -> AP1
}
```

`shortcut: motif-only 归一后仍可能是 UNKNOWN / TATA_01，这类标 named = FALSE，consensus() 里不进 A/B`

---

### 4.4 `$diagnostics` 的四列（v2 新增）

v1 只有 `uncovered` 一个数，审核指出它会误导：`enricher()` 的 `GeneRatio` 分母是**被注释到的基因数**（850），不是选中基因数（897）。200 个「未覆盖」里有一部分根本没进入检验。

四个数，**两个不同的分母，不能混用**：

```text
n_input            897   进入 evidence() 的基因数（selected）
n_annotated        850   n_input 中被任一所选集合注释到的（= GeneRatio 分母）
n_sig_covered      697   n_input 中被【显著】调控者覆盖的
uncovered_sig      200   n_input - n_sig_covered

coverage_annotated = n_annotated   / n_input   = 850/897 = 0.948
coverage_sig       = n_sig_covered / n_input   = 697/897 = 0.777
```

**v2.1 校准（2026-10-25，实测）**：审核曾指出 850「不是 GeneRatio 分母」。实测**两个方向都核对过**，结论是 850 **确实就是** `GeneRatio` 分母：

```r
# 合成例：selected = {g1,g2,g9,g8}，只有 g1,g2 被任何集合注释
#> S1  GeneRatio 2/2      # 分母 = 2 = 被注释的数
#> S2  GeneRatio 1/2

# airway：selected 897，被任一 c3.tft 集合注释 850
#> 首行 GeneRatio: 157/850   # 157 是 Count，850 是分母
#> 850 == length(intersect(selected, unique(t2g$gene))) -> TRUE
```

**但审核的实质担忧成立**：`coverage_annotated = 0.948` 看起来很高，却只说明「基因出现在某个集合里」，不说明「该来源能对这批基因做出判断」。用 0.948 当 `absent`/`untested` 的判据过于宽松。

**v2.1 决定**：`coverage_floor` 用 **`coverage_annotated`**（注释覆盖），因为 `absent` 的定义就是「注释到了但没显著」——注释不到的基因本来就不能被判为 `absent`，只能判 `untested`。所以：

- `coverage_annotated` 决定**该来源整体**是否 `untested`（低于 floor 时，该来源全部降级）；
- 单个基因的两个数分开用：`n_annotated` 之外的是 `untested`，`n_annotated` 之内但不在 `n_sig_covered` 里的是 `absent`。

**`coverage_floor` 的默认值不能拍 0.30**（v2 的 0.30 无数据支撑，审核也已指出）。v2.1 改为：**默认 `NULL`，即不设阈值**；只有用户显式给出 `design$coverage_floor` 时才做整体降级。理由：在 22899 个基因的 universe 上，TF 集合只覆盖 18759，任何固定阈值都可能在某个集合上误判，而「不设阈值 + 逐基因三分类」是保守但不会错的选择。

---

### 4.5 `consensus()`（P2）

```r
consensus(x, entity = c("regulator", "gene"), min_independent = 2)
```

| verdict | 条件 |
|---|---|
| `supported` | 同一 entity 被 ≥ `min_independent` 类不同来源支持，且至少一个 `level` 为 B 及以上 |
| `single_source` | 只有一类来源支持 |
| `absent` | 来源覆盖充足（`coverage_annotated >= coverage_floor`）但该实体未出现 |
| `untested` | 该来源未运行，或覆盖不足 |

**`absent` 与 `untested` 的区分是本规范最有价值的观点**：前者是证据，后者是缺证据。

### 4.6 `writable`（P2）

| level | pattern | example（填此行数据后） |
|---|---|---|
| A | 受 `<扰动>` 调控（实测敲低/过表达） | "HIF1A 敲低后该模块整体下调" |
| B | `<实体>` 的靶集合富集 / 文献支持的调控 | "FOXO3 的靶基因在本列表中富集" |
| C | 表达相关 / 功能关联，**不可写成调控** | "HIF1A 与 PGK1 在组织中正相关" |
| absent | 未检测到该来源证据（**不等于否定**） | "未发现 miRTarBase 支持的直接证据" |

---

### 4.7 `interpret(bundle)`（P3，v2 改）

v1 说「只新增一个分支」低估了改动量。审核实测：`process_enrichment_input()` 对非四类合法对象**直接 abort**（`R/interpret.R:207-210`），而 `interpret(list(tf = enrichResult, mirna = enrichResult))` 能跑。

**v2 方案**：在 `.interpret_run()` 的**分派点**判断：

```r
if (inherits(x, "evidenceBundle")) {
    x <- process_evidence_input(x)                 # 拆成有名 list
}
x <- process_enrichment_input(x, n_pathways)       # 既有逻辑一行不改
```

`process_evidence_input()` 负责：取 `$regulators` 里 `level %in% c("A","B")` 且 `verdict %in% c("supported","single_source")` 的实体 → 组装 `List(tf = enrichResult, mirna = enrichResult)`。

**P3 的第一版可以不引入 `evidenceBundle` 专属逻辑**：既然 `interpret(list(tf=..., mirna=...))` 已通，先交付「拆包」，`verdict`/`unsupported_claims` 后加。

**prompt 预算（v2 新增）**：225 行显著 TF 的 `geneID` 列共 40,585 字符 ≈ 10.1k token，已超 `max_tokens = 8192`。§4.7 必须显式给 `n_pathways` 默认值；§7.3 的「有 key」判据要加前置断言**prompt 未触顶**，否则「`unsupported_claims` 为空」什么都证明不了。

**`unsupported_claims` 的实现（v2 澄清）**：`evidence_status()` 的第二个参数是 `enrichResult`/`data.frame`，**不能直接把 bundle 传进去**（会在 `gcol` 检查处 abort）。实现对 bundle 里每个 `$tf`/`$mirna` 分别调用再合并。

---

## 5. 依赖与兼容

只列与本期相关的变化：

| 依赖 | 现状 | 变化 |
|---|---|---|
| `gson` | `Imports (>= 0.0.7)` | 不变 |
| `aisdk` | `Imports` | 不变 |
| `aisdk.mcp` | 不在 `Imports` 也不在 `Suggests` | **P4 进 `Suggests`** |
| `yulab.utils` | `Imports (>= 0.2.3)` | 不变（D4 用到 `user_dir()`） |

**向后兼容**：`evidence()` / `consensus()` / `mechanism_sets()` 是新增导出（`NAMESPACE` 已核对无同名）。`interpret()` 只在 `.interpret_run()` 加一个前置分派，既有四种输入的字节输出由现有 `test-interpret-*.R` 快照守护。

---

## 6. 实现顺序

| 期 | 新文件 | 测试 |
|---|---|---|
| P1 | `R/mechanism_sets.R`、`R/evidence.R`、`R/xref.R`；改 `R/zzz.R` | `test-mechanism-sets.R`、`test-evidence.R`（含 T0） |
| P2 | `R/consensus.R`、`R/writable.R`、`R/ppi_evidence.R` | `test-consensus.R`、**`test-ppi-evidence.R`**（v2 补：现有 `getPPI` 无测试） |
| P3 | `R/interpret-evidence.R` | `test-interpret-bundle.R` |
| P4 | `R/evidence-mcp.R` | `test-evidence-mcp.R`（全 mock） |

**冻结项收窄（v2 改）**：只冻结两处——①顶层元素名（§4.2 表），②`$regulators` 最小组：`regulator / source / term / named / level / padj_within_source / setSize / Count`。其余列名可追加不可改名。

---

## 7. 验证判据

### 7.1 P1 验收（全部可离线）

**T0（v2 新增，防 D1 回归）**

```r
expect_error(enricher(genes, TERM2GENE = G), "colnames")
x <- enricher(genes, gson = G, ...)
expect_equal(x@keytype, "SYMBOL")
expect_equal(x@ontology, "MSigDB C3:TFT")
```

**T1（airway 复现）**——`minGSSize` 必须显式，两个值都断言：

| 断言 | 期望 |
|---|---|
| `$tf` 被测集合数，`minGSSize = 5` | **1030** |
| `$tf` 被测集合数，`minGSSize = 10` | **1023** |
| `$tf` 中 `padj < 0.05`（两种参数下） | **225** |
| `AACTTT_UNKNOWN` | `named = FALSE`，rank 1 |
| `$diagnostics$n_input` | **897** |
| `$diagnostics$n_annotated` | **850** |
| `$diagnostics$n_sig_covered` | **697** |
| `$diagnostics$uncovered_sig` | **200** |
| `$mirna` 前三 | `MIR23C` / `MIR23A_3P_MIR23B_3P` / `MIR3910` |

> **v2 改**：v1 只写 1030 未标参数（按默认 10 跑会得 1023，测试会挂）；`uncovered = 200` 未区分口径。现在四列分开断言。

**T2（测试数据可执行性，v2 新增）**

审核指出 §7.1 依赖的 `datasets/de_table.tsv` 与 `cache/*.gmt` **都不在包里**，包自身 testthat 访问不到。方案：

- 897 个基因向量存为 `tests/testthat/fixtures/airway_selected_897.txt`（纯文本，约 7 KB）；
- 从 1030 个 set 中抽 **20 个**（含 `AACTTT_UNKNOWN` 与 ≥3 个 GTRD 集）做 `fixtures/tiny_tft.gmt`，对这批数据重算期望值；
- 1030/225/200 那组完整数字**只在书仓库验证**，作为 P1 人工验收，不作 testthat 断言。

### 7.2 P2 验收

**网络定级规则（v2 重写）**

v1 用 `pscore >= 0.4 → B`，审核实测该判据不可达且通道语义用错。2026-10-25 实测：

```text
8 个缺氧基因 24 条边：pscore 最大 0
NFKB1-NFKBIA（文献 physical 0.999）：score=0.999 pscore=0.000 escore=0.999 dscore=0.900
```

**`pscore` 恒为 0，连真物理互作都是 0**；携带证据的是 **`escore`（experiments）** 与 `dscore`（databases）。

v2 规则：

| 条件 | level | 依据 |
|---|---|---|
| `escore >= 0.4` | **B** | 有实验证据 |
| `dscore >= 0.4` | **B** | 有数据库收录 |
| 其它（纯 `tscore`/`ascore`） | C | 只有文本/共表 |

`pscore` 保留为布尔标注 `channel_experiment`，**不参与定级**。

| 断言 | 期望 |
|---|---|
| HIF1A 缺氧模块（8 基因） | 实测 **22–24 条边，5–7 条 B / 15–19 条 C**（见 v2.2） |
| NFKB1–NFKBIA | **B**（escore 0.999） |
| 任何网络证据 | `level != "A"` |

**M3 的实现路径（v2 改）**：v1 说「`network_type` 加 `"all"`」是错的——`getPPI()` 的 `igraph` 分支硬编码 `res[, c(3,4,6)]`（`R/ppi.R:164-168`），**只有 `score` 一列**，加不加 `"all"` 都拿不到 `escore`。

v2 方案：P2 **不用 `getPPI()`**，改用 `getPPI(output = "data.frame")`（已返回全 13 列，实测可用）+ 新代码里的定级逻辑。`getPPI()` **一行不改**。

### 7.3 P3 验收

| 场景 | 期望 |
|---|---|
| 无 API key | `eval: FALSE`，构建不失败 |
| 有 key，且 prompt 未触顶 | `unsupported_claims` 为空 |
| **prompt 触顶** | 显式 warning，不静默截断 |
| 输入含 C 级实体 | 该实体不出现在 `key_mechanisms` |

### 7.4 P4 验收

端点不可达时：`$diagnostics` 记录 source、错误类别、时间；`evidence()` **不失败**，其余来源照常返回。全部 mock。

---

## 8. 风险

| 风险 | 影响 | 处理 |
|---|---|---|
| 覆盖率假象 | 无结果可能是没注释 | `diagnostics` 四列；`absent`/`untested` 分开 |
| 过滤顺序 | 先过滤会把 `single_source` 误判 `absent` | 裁决用全表；`print_*` 只影响打印 |
| motif-only 假调控者 | `AACTTT_UNKNOWN` 会排第一 | `named = FALSE`，不进 A/B |
| MSigDB 版本漂移 | 换 release 数字全变 | `release` 是参数；`.meta` 不匹配时 warn |
| `pscore` 语义误用 | 连真互作都判 C | 用 `escore`/`dscore`；`pscore` 只做标注 |
| `.onAttach` 覆盖用户设置 | 缓存永不命中 | P1 一并修（D4） |
| prompt 超出上下文 | 静默截断导致假通过 | §7.3 加触顶断言 |

---

## 9. 明确不做的

- ❌ 不改 `getPPI()` 的返回值（P2 用 `output = "data.frame"` 绕过）
- ❌ 不删 `interpret()` 的 `regulatory_drivers`（破坏性变更）
- ❌ 不新增相似度/聚类（`reduce_redundancy()` 已够）
- ❌ 不做自动「猜」缺席证据
- ❌ P1 不引入任何网络依赖

---

## 10. 审核已确认正确的部分（保留备查）

独立审核实测确认：1030/225 可复现、miRNA 前三可复现、`AACTTT_UNKNOWN` 确为 rank 1、`c4.cgn` 200 + `c2.cgn` 404、`gson=` 路径保留三个元数据槽、`interpret(list(tf=, mirna=))` 可跑、D2 分源 BH 理由成立、D7「不猜」原则正确、新增导出无命名冲突。

**审核未能验证的**（不阻塞）：`aisdk.mcp` 端点可达性（P4）、`gsname` 超长是否被 `enrichResult@ontology` 拒绝、`coverage_floor = 0.30` 的合理默认值。

---

## 附：v2.2 实测校正（2026-10-25，最终审核后）

规范中有三条断言与实现实测不符，以此为准：

1. **§7.2「HIF1A 缺氧模块全部 C」不成立。** 实测 8 基因查询返回 22–24 条边，其中 5–7 条为 B（`EGLN3–HIF1A` escore 0.999、`PGK1–HIF1A` 0.463、`CA9–HIF1A` 0.720 等）。`escore >= 0.4 -> B` 的规则本身正确，错的是「全部 C」这条历史断言。

2. **`coverage_floor` 默认 `NULL`**（不是 0.30）。默认下逐基因三分类仍给出 `absent`/`untested`；`coverage_floor` 只在调用者显式给出时做整体降级。

3. **§4.2 参数表**：`evidence()` 实际有 9 个参数（`print_pvalue`/`print_qvalue` 只影响打印），`design` 有 10 个字段但只有 `species`/`treatment`/`gene_col`/`coverage_floor`/`background` 被消费。

最终审核另指出并被采纳的设计修订：`$writable` 不再独立计算，直接取 `consensus()` 的结果（逐行判 verdict 会让同一实体得到互相矛盾的句子）；`unsupported_claims` 的家族前缀判定改为单向（`startsWith(vocab, tok)`），反向会让 1018 个真实符号免检。
