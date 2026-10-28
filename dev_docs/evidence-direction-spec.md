# `evidence()` 重设计：方向为主轴，四层证据

状态：待实现。日期 2026-10-28。
本文件与 `dev_docs/adjudicate-migration-spec.md` **配套**，见 §9 的一致性约束。
本文件不改动 `mechanism-exit-spec.md` 的其余各节。

背景一句话：`evidence()` 现在收集的是「谁的靶集在基因列表里富集」，而一个没有方向的基因列表问不出机制——同一批基因可能是某因子被激活的结果，也可能是它被抑制的结果。本次把**方向**变成横贯各层的输入约束，并补上一直缺的通路层。

---

## 1. 决策

| # | 决策 | 理由 |
|---|---|---|
| **D1** | 规范输入是**带方向的差异结果**：DE 表（`gene` + `logFC`）或 named numeric | 无方向的列表只能回答「像不像某个因子的靶」，回答不了「这个因子在起什么作用」 |
| **D2** | 无方向输入**不报错**，方向层全部报 `untested` 并在 `print()` 说明 | `untested` 的既有语义正是「这个来源问不了这个问题」。硬报错会挡住「手上只有一列命中基因」的真实情况。若日后改为硬报错，是一行 `rlang::abort` |
| **D3** | 富集结果从**输入**降为**一层可选证据**：`pathway = TRUE \| <enrichResult> \| FALSE` | 现在的 `evidence(enrichResult)` 从 `@result` 反推基因，被 `pvalueCutoff` 静默截断（airway 897→483）。基因来自 DE 表才是完整的，富集结果只负责「看哪几条通路」 |
| **D4** | 删除 `evidence(enrichResult)` 入口，报错指路 `evidence(genes, pathway = result)`；**保留** `evidence(gseaResult)` | 前者有更好的替代且本身有损；后者的 `@geneList` 是排序本身，无法用 `pathway=` 表达 |
| **D5** | 四层证据，**L3 单独成表**，不进 `regulators` | L3 由同一批 TF 集合派生，若作为行进入 `regulators`，`adjudicate()` 会把它当成独立见证者——与当初删掉 biotrigger 是同一个错误 |
| **D6** | L2 的方向由「交集限定 + 符号检验」决定，取代「所有符号相同才算 up」 | 现在的规则在混合列表上必然给出 `mixed`，那是尺子的问题，不是数据的问题 |
| **D7** | 通道增减按**实测**：biotrigger / similarity / network 各跑一次对比，只有指标变好才回加 | 当初删除的理由是实测；重新评估也必须用实测，不能因为「设计变了」就假定它们现在有用 |
| **D7'** | **已测（2026-10-28），结果见 `mechanism-accuracy.md` §R3**：biotrigger 与基线逐位相同（MRR 0.083 → 0.083），不回复；network 的 `$ppi` 从未被任何代码读取（`mechanism_pair_evidence()` 零调用点），保持默认关闭并改标为「参考信息」；similarity 函数已删除，未测 | 三声明的实测落地 |

---

## 2. 输入契约

```r
evidence(genes,
         pathway  = FALSE,     # FALSE | TRUE | <enrichResult/gseaResult>
         upstream = "tf",
         network  = FALSE,
         design   = list(),
         universe = NULL,
         minGSSize = 5, maxGSSize = Inf,
         tft = NULL, print_pvalue = 1, print_qvalue = 1)
```

`genes` 接受四种形态，前两种是规范形态：

| 形态 | 方向 | 说明 |
|---|---|---|
| `data.frame` + `design$logfc_col` | **有** | 推荐。DE 表的原样形态 |
| named numeric | **有** | `setNames(de$log2FoldChange, de$gene)` |
| `gseaResult` | 有（排序） | 走 GSEA 路径，见现有 `.mechanism_from_gsearesult()` |
| `character` | **无** | 允许；L1/L2/L3 的方向层全部 `untested`，`print()` 说明原因 |
| `enrichResult` | — | **报错**（D4），指路 `evidence(genes, pathway = result)` |

`design` 新增一个字段：

| 字段 | 默认 | 用途 |
|---|---|---|
| `logfc_col` | `"log2FoldChange"` | `genes` 为 data.frame 时，取哪一列作方向 |

不加自动猜测列名的逻辑：猜错会静默改变方向层的结果。

---

## 3. 各层定义

方向以「调用方对比里的方向」为准：`up` = 在这批数据上升高。

| 层 | 问题 | 数据源 | 方向的含义 | 输出位置 |
|---|---|---|---|---|
| **L1** 上游调控者 | 谁的靶基因在**上调**里富集 | `upstream` 里的基因集来源 | 上调富集 ⇒ 该因子的程序在增强 | `$regulators` |
| **L2** 实测扰动 | 敲低该因子，**我的上调基因**是否随之下降 | GPSAdb | 与 L1 方向一致 ⇒ 支持它驱动这些基因 | `$regulators` |
| **L3** 通路层 | 这条**富集通路**由谁调控 | 同一批基因集 + 通路基因 | 通路整体上调 + 靶集在通路里富集 | **`$pathways`（独立表）** |
| **L4** 网络 | 互作是否支持 | STRING | 无方向，仅旁证 | `$ppi`（已有，独立） |

**`n_sources` 只数 L1/L2。** L3 与 L4 是旁证与细化，不是投票者。

---

## 4. L1：已有，仅确认

方向拆分已在 4.21.4.004 实现（`.mechanism_directional_collection()`）：带符号输入时上调、下调两半分别检验，统一 BH。本次不改数值逻辑。

需要改的只有一处：`direction` 列的取值必须与 §3 的定义一致（`up` = 该集合在**上调**半显著）。已如此。

---

## 5. L2：重做方向判定

### 5.1 现状的问题

`R/evidence_perturbation.R` 现在这样定方向：

```r
signs <- sign(moved$log2fc)   # moved = 输入基因里显著移动的
direction <- if (all(signs > 0)) "up" else if (all(signs < 0)) "down" else "mixed"
```

`moved` 已经是「输入基因的子集」，但输入本身是混合方向的一批基因。一个因子的敲低必然让一部分基因升、一部分降，所以只要列表够大，`all(signs > 0)` 几乎不可能成立，**`mixed` 是这把尺子必然的读数**。实测在气管数据上，所有 supported 的 `measured_direction` 全是 `mixed`。

### 5.2 新规则

对每个数据集，分别对上调半和下调半计算：

1. 取 `up` 与该数据集显著移动基因的交集，记 `n_rise` / `n_fall`
2. `binom.test(n_rise, n_rise + n_fall, 0.5)` 得 `sign_padj`
3. `measured_direction`：
   - `sign_padj < 0.05` 且 `n_rise > n_fall` → `"up"`
   - `sign_padj < 0.05` 且 `n_fall > n_rise` → `"down"`
   - 否则 → `"mixed"`（**方向不成立**，而不是没测）
4. 下调半同样处理，结果落在 `down_*` 列

### 5.3 新增列（`$regulators`，仅 L2 行）

| 列 | 含义 |
|---|---|
| `n_agree` / `n_disagree` | 上调半里，与该因子被激活的预期一致 / 不一致的基因数 |
| `sign_padj` | 上述 2×N 二项检验的 BH 校正 p |
| `measured_direction` | 语义不变（数据集让这些基因往哪走），取值规则按 5.2 |

`measured_direction` 的名字与语义**不改**，只把判定规则换掉：`consensus()` 已按此列拆分，不受影响。

### 5.4 与 L1 的一致性（D5 的落地）

`adjudicate()` 的 verdicts 表新增一列 `perturbation_effect`（三值：`"activator"` / `"repressor"` / `NA`）：

- L1 有方向、L2 有方向、且两者符合该 `method` 的预期 → `"activator"`（该因子促进这个程序）
- 两者都有但方向相反 → `"repressor"`（该因子压制这个程序）
- 任一缺失，或符号检验判不出方向 → `NA`

**不是布尔值，这是实现时改掉的一个设计错误。** 原计划叫 `direction_concordant` 并返回 `TRUE/FALSE`，但在缺氧数据上，「敲低后我的上调基因继续上调」会得到 `FALSE`，读起来像「证据矛盾」——其实它说明该因子**抑制**这批基因。激活与抑制都是机制，压成一个真假值会丢掉更有信息量的一半：实测 25 个 supported 里有 7 个是这种情况。

**它不改变 verdict。** `supported` 的判定条件保持 `adjudicate-migration-spec.md` §5.2 原样。理由：并入必要条件会让现有 25 个 supported 全部重判，而那批结论本身没有错——它们只是缺一个加强标记。

`method` 的方向预期写在一处（`mechanism_method_sign()`），不在别处复制。

**取哪个字段**：GPSAdb 的 `method` 是混合词表（实测 `siRNA`、`knockout`），而 `method_label` 是规范化的（`Knock_down`、`Knock_out`）。所以优先用 `method_label`，缺失时退回 `method`，两者都认不出 → `NA`（**不猜**）。

| 类别 | 包含 | 因子的预期 | 与 `direction` 的关系 |
|---|---|---|---|
| 去除 | `Knock_down` / `Knock_out` / `siRNA` / `shRNA` / `CRISPR` / `knockdown` / `knockout` | 因子减少 ⇒ 靶基因应下降 | **相反** = 一致 |
| 增加 | `Over_expression` / `overexpression` / `activation` | 因子增加 ⇒ 靶基因应上升 | **相同** = 一致 |
| 认不出 | 其余 | 无法判定 | `NA` |

---

## 6. L3：通路层（新增）

### 6.1 通路从哪来

| `pathway` 取值 | 行为 |
|---|---|
| `FALSE`（默认） | 不做 |
| `TRUE` | 在 `evidence()` 内部跑 `enrichGO(ont = "BP")`（`design$species` 决定 OrgDb） |
| `<enrichResult>` | 直接用调用方给的富集结果，不重算 |
| `<gseaResult>` | 同上，另外可带每条通路的方向（NES 符号） |

默认 `FALSE`：`enrichGO` 要 OrgDb 且需要数秒，默认打开会让最简单的调用变慢且可能因缺包失败。

### 6.2 统计量与成本

对每条通路，把**通路自己的基因**作为输入，用同一批上游集合跑 `enricher()`——不新增统计量。

已实测（缺氧上调 351 基因 → GO BP 125 条）：

```
[response to hypoxia] 33 基因 → HIF1A p=3.9e-14, SOX4 6.2e-13, RELA 7.4e-12, IRF3, IRF7
```

第一名即正确答案。成本控制两条，都已量过：

1. **按基因集合去重**：125 条术语只有 84 个不同基因集合（省 33%）。GO 的近义术语常常字面不同而成员完全相同。
2. **只算 `pathway_top` 条**（默认 20），按 `padj` 排序后去重再截断。

### 6.3 `$pathways` 表

```
pathway, pathway_description, pathway_direction, n_pathway_genes,
regulator, source, term, Count, setSize, padj_within_pathway, geneID
```

`padj_within_pathway` = 该通路内所有被测集合的 BH 校正值。「每条通路自己校正」与 `padj_within_source` 是同一套口径，不跨通路合并。

`pathway_direction` = 该通路基因在调用方对比里的整体方向（`up` / `down` / `NA`）。

### 6.4 不参与裁决

`$pathways` 不进 `$regulators`，`adjudicate()` 不读它（D5）。它在 `interpret(adjudication)` 里作为**上下文**出现，用于写「X 调控了 Y 通路」，但不得据此产生新的 `supported` 实体。

已核实不会被现有代码误读：`.process_evidence_input()` 只遍历 `x$provenance$upstream` 里出现的名字（`R/interpret_evidence.R:2-4`），而 `pathways` 不在 `upstream` 里，所以新增的 `$pathways` 不会漏进旧的 prompt 构造路径。

---

## 7. 被暂缓的通道（D7）

不改代码，只测量。用 `dev_docs/mechanism-accuracy-benchmark.R`（30 个 GPSAdb 真值数据集）跑四组：

| 组 | 配置 |
|---|---|
| baseline | `upstream = c("tf", "perturbagen")` |
| +biotrigger | 加 `evidence_biotrigger()` |
| +similarity | 加 `evidence_similarity()`（现已删除，需先恢复） |
| +network | `network = TRUE` |

指标：`recall@1`、`recall@10`、`MRR`、新增的 supported 实体是否包含真值调控者。

**准入判据**：新增通道使 MRR 提升或使真值调控者进入 recall@10，且不引入仅由同一数据集重复计数得来的 supported。任一不满足，保持删除。

四组结果写进 `dev_docs/mechanism-accuracy.md`，作为去留依据。

**已执行（2026-10-28）**：结论见 `mechanism-accuracy.md` 的 §R3。biotrigger 不回复，network 保持关闭并改标参考，similarity 未测（实现已删）。

---

## 8. 改动清单

| 文件 | 改动 |
|---|---|
| `R/evidence.R` | 输入契约（§2）、`pathway` 参数与 L3 调用、`$pathways` 组装、`print.evidenceBundle` 报告层可用性与无法提问的原因 |
| `R/mechanism_design.R`（或 `mechanism_design()`） | 新增 `logfc_col`、`pathway_ontology`、`pathway_top` |
| `R/evidence_pathway.R`（新） | L3：`.mechanism_pathway_layer()` |
| `R/evidence_perturbation.R` | §5 的方向判定；新增 `sign_padj` / `n_agree` / `n_disagree` |
| `R/mechanism_xref.R` 或新文件 | `mechanism_method_sign()`：一处定义 method → 方向预期 |
| `R/consensus.R` | 新增 `perturbation_effect` 列（不参与 verdict） |
| `R/interpret_evidence.R` | `.process_evidence_input()` 适配新 bundle；L3 作为上下文进入 prompt |
| `NAMESPACE` / `man/` | 由 roxygen 生成 |
| `dev_docs/mechanism-accuracy.md` | §7 的四组结果 |

不改：`mechanism_sets()`、`mechanism_tfunion()`、集结级定级规则、`padj_within_source` 的口径。

---

## 9. 与 `adjudicate-migration-spec.md` 的一致性约束

两份 spec 同时生效，必须满足：

1. **顺序**：本文件先落地（R1–R3），`adjudicate` 改名最后做（R4）。`adjudicate` 的冻结项依赖 `$regulators` 的列，先改 evidence 再冻结，只冻结一次。
2. **冻结列更新**：`adjudicate-migration-spec.md` §5.1 的冻结列需追加 `perturbation_effect`。这是本文件 §5.4 带来的唯一变更，其余列名与顺序不动。
3. **verdict 判定条件不变**：`supported` 仍需 ≥ `min_independent` 个不同来源且至少一条 measured。`perturbation_effect` 是报告，不是条件（D5）。
4. **`interpret(evidenceBundle)` 仍然要删**：`adjudicate-migration-spec.md` D2 不变。本文件新增的 `$pathways` 不构成保留该分支的理由。
5. **prompt 对照表需扩一行**：`adjudicate-migration-spec.md` §4.2 的措辞表要补「L3 通路上下文不得被写成新的机制实体」。
6. **同一概念只有一个名字**：`direction` = 调用方对比里的方向；`measured_direction` = 扰动数据集让基因走的方向。两份 spec 及各层输出都不得再引入第三种说法。

---

## 10. 测试清单

新增 `tests/testthat/test-evidence-direction.R`：

1. DE 表输入按 `logfc_col` 分出上下调，且两半独立检验。
2. 未给方向时，L1/L2/L3 方向层报 `untested`，`print()` 含说明。
3. `evidence(enrichResult)` 报错并指路 `pathway=`。
4. `gseaResult` 输入仍走 GSEA 路径（排序未被丢弃）。
5. `pathway = FALSE` 时 `$pathways` 为 `NULL` 或空表，且不触发 `enrichGO`。
6. `pathway = <enrichResult>` 时复用给定的通路，不重算。
7. L3 按基因集合去重：两条成员完全相同的术语只算一次。
8. **L3 不进 `$regulators`，且 `n_sources` 不受 L3 影响**（D5 的核心断言）。

新增 `tests/testthat/test-evidence-concordance.R`：

9. `measured_direction` 在符号检验不显著时为 `mixed`，显著时取多数方。
10. `perturbation_effect`：knockdown 下 L1=up、L2=down → `"activator"`；L2=up → `"repressor"`；缺一方 → `NA`。overexpression 下预期方向相反。
11. `perturbation_effect` 为 `"repressor"` 时 verdict 不变（证明它不是判定条件）。

修改：`test-mechanism.R` 中依赖 `evidence(enrichResult)` 的用例，改走 `pathway=`。

---

## 11. 验收判据

1. `grep -rn "evidence(enrichResult\|inherits(genes, \"enrichResult\")" R/` 只剩报错分支。
2. 缺氧数据端到端：`$pathways` 里 `response to hypoxia` 的 top regulator 是 `HIF1A`。
3. 同一份数据，`upstream = c("tf","perturbagen")` 下 `measured_direction` 不再**全部**是 `mixed`（与当前全 `mixed` 形成对照，值记入 `mechanism-accuracy.md`）。
4. `n_sources` 在 `pathway = TRUE` 与 `FALSE` 两种配置下对同一实体一致。
5. §7 四组对比结果写入 `mechanism-accuracy.md`，每个通道有明确的去留结论。
6. `make check` 通过（先读 `Makefile`）。

---

## 12. 风险

- **风险**：`binom.test` 的通路在每次 L2 行上多一次检验，`sign_padj` 需在同一数据集内做 BH，否则等于多重比较未校正。缓解：按 `dataset` 分组校正，写进 §5.2。
- **风险**：L3 的通路数量由调用方的富集阈值决定，一个宽松阈值会带来上百条通路。缓解：`pathway_top` 截断 + 按集合去重，两者都可配置。
- **风险**：删 `evidence(enrichResult)` 会打断现有用法。未发布，报错已指路，接受（与 `adjudicate` spec D2 同一处理方式）。
- **回退**：R1–R3 一组提交，R4 单独提交。
