# `adjudicate()` 落地规范：把裁决从 `interpret()` 里拿出来

状态：待实现。日期 2026-10-28。
本文件**取代** `dev_docs/mechanism-exit-spec.md` 的 §4.5（`consensus()`）与 §4.7（`interpret(bundle)`）；其余各节仍然有效。

**执行顺序**：`dev_docs/evidence-direction-spec.md` 先落地（它改变 `$regulators` 的列），本文件的改名与冻结在其后。冻结只做一次。

背景一句话：`consensus()` 做的是裁决（决定哪些主张站得住），但名字说的是共识；而裁决目前藏在 `interpret()` 内部，用户在花模型调用之前看不到那张 verdict 表。本次改动把裁决变成一等公民。

---

## 1. 三条决策

| # | 决策 | 理由 |
|---|---|---|
| **D1** | `consensus()` 改名为 `adjudicate()`，**不留别名** | 该函数 2026-10-25 才加入（`git log --diff-filter=A`），未进任何发布版，改名的外部成本为零。"consensus" 在生物信息学里专指聚类共识（ConsensusClusterPlus），与"判定"无关 |
| **D2** | 删除 `interpret(evidenceBundle)` 分支 | 保留了它，用户就可以绕过裁决直接花模型调用，D1 的意义会被抵消。裁决必须是**显式一步** |
| **D3** | `interpret(adjudication)` 的措辞受 verdict 约束 | 这是"可查"从口号变成约束的地方：`single_source` 不得被写成机制，`untested` 不得被写成"未发现" |

## 2. 目标 API

```
evidence(genes)  →  adjudicate(bundle)  →  interpret(adjudication)
  收集 + 分级           裁决（四类 verdict）       语言化（唯一调模型的一步）
```

`adjudicate()` 与 `evidence()` 都不调用模型；只有 `interpret()` 调。

## 3. `adjudicate()`

签名与现有 `consensus()` 完全一致（`R/consensus.R:71`）：

```r
adjudicate(x, min_independent = 2, coverage_floor = NULL)
```

**输入**：必须是 `evidenceBundle`（`evidence()` 的返回）。其他输入报错并指路：

> `` `x` must be an `evidenceBundle` from `evidence()`. ``
> i: "Collect the evidence first, then adjudicate it."

**输出**：新的 S3 对象 `adjudication`，结构：

```r
structure(
  list(
    verdicts    = <现有 data.frame，列名与顺序冻结>,   # 见 §5
    evidence    = <输入 evidenceBundle 原样>,
    diagnostics = <evidenceBundle$diagnostics>,
    parameters  = list(min_independent = , coverage_floor = ,
                       species = , design = ),
    provenance  = list(evidence_sha = , call = , time = )
  ),
  class = c("adjudication", "list")
)
```

`print.adjudication()` 的摘要必须让读者**不用跑模型就知道该不该继续**：四类 verdict 各多少条、`supported` 里 top 3 的 `regulator` 与 `max_level`、以及 `diagnostics` 里 `untested` 的来源数。

**不允许**：把 verdict 计算逻辑搬进新函数或改数值。`adjudicate()` 的 body 就是原 `consensus()` 的 body，唯一改动是返回包装。

## 4. `interpret()` 的分派

### 4.1 两个入口，第三类明确报错

| 输入 | 行为 |
|---|---|
| `enrichResult` / `gseaResult` / `compareClusterResult` / `data.frame` | 老路径，行为与 prompt 文本**逐字节不变**（由现有 `test-interpret-*.R` 快照守护） |
| `adjudication` | 新路径：按 verdict 构造 prompt |
| `evidenceBundle` | **报错** |

报错原文（照抄，测试会断言这段）：

```
`interpret()` no longer accepts an `evidenceBundle`.
i: Call `adjudicate()` first, then pass the `adjudication` to `interpret()`.
i: The verdict decides which claims may be written; skipping it would let a
   single-source finding be narrated as a mechanism.
```

### 4.2 `interpret(adjudication)` 的措辞约束

prompt 里必须带上 verdict 与 level 两列，并在指令中写死这张对照表——**模型不得升级措辞**：

| verdict | level | 允许写成 | 禁止写成 |
|---|---|---|---|
| `supported` | A | 受 X 调控（实测敲低/过表达） | —— |
| `supported` | B | X 的靶集合富集 / 有 ChIP 或文献支持的调控 | "被证实""因果" |
| `single_source` | 任意 | 有单一来源提示 / 待独立验证 | 调控、机制、驱动 |
| `absent` | 任意 | 未检测到该来源的证据（**不等于否定**） | "没有""不存在" |
| `untested` | 任意 | **不写成任何结论**；必要时写"未测" | "未发现" |

实现要点：

- 复用 `mechanism_writable()`（`R/consensus.R:263`）已给的句子，不要把措辞再写一遍——两处措辞必然漂移。
- 只有 `supported` 与 `single_source` 进 prompt 的主张段；`absent` / `untested` 只作为"边界说明"出现。
- prompt 顶部加一行固定声明：*"Only `supported` findings may be stated as mechanisms; `single_source` is a hypothesis; `untested` is not an absence."*
- **通路层是上下文，不是投票者**：`adjudication$pathways`（见 `evidence-direction-spec.md` §6）只用于写「X 调控了 Y 通路」这类句子，不得据它产生新的机制实体，也不得把它当成第三条独立证据。它与 `biotrigger` 当初被删掉是同一个陷阱：由同一批集合派生出来的东西不能给自己作证。

## 5. 冻结项（不许改）

1. `verdicts` 的列，**顺序以 `consensus()` 的实际输出为准**（下面这一行是从代码里跑出来的，不是回忆的）：
   `regulator / sources / n_sources / max_level / verdict / padj / family / evidence_kind / direction / measured_direction / writable`。
   本文档 2026-10-28 的早期草稿把 `verdict / padj / family / max_level` 写成了另一种次序，与代码不符，已按代码订正。
   末尾**追加** `perturbation_effect`（由 `dev_docs/evidence-direction-spec.md` §5.4 引入；三值 `activator`/`repressor`/`NA`，不是布尔）。
   列名可追加，不可改名，不可插在中间。
2. 四类 verdict 的判定条件（`supported` 需 ≥ `min_independent` 类不同来源、同一实体、且至少一个 B 及以上）。
3. `absent` 与 `untested` 的区分，以及 `coverage_floor` 的语义。
4. 老路径 prompt 的字节输出。

## 6. 改动清单（按文件）

| 文件 | 改动 |
|---|---|
| `R/consensus.R` → `R/adjudicate.R` | `git mv`；函数改名；返回 `adjudication`；新增 `print.adjudication()`；roxygen 里所有 `consensus()` 引用改为 `adjudicate()`（含 `@seealso`、`@param verdict A adjudicate() verdict`） |
| `R/evidence.R:736` | `cv <- consensus(bundle)` → `adjudicate(bundle)` |
| `R/collapse.R:42` | `consensus_findings(consensus(x))` → `consensus_findings(adjudicate(x))`；`consensus_findings()` 本身不对外，保持原名（或改为内部 `.consensus_findings()`，二选一，不要两个都留） |
| `R/mechanism_graph.R:56` | 同上 |
| `R/interpret_dialogue.R:73` | `x$writable %||% consensus(x)` → `adjudicate(x)` |
| `R/interpret.R:548-559` | 删除 `evidenceBundle` 分支，替换为 §4.1 的报错；新增 `adjudication` 分支走 `.process_adjudication_input()` |
| `R/interpret.R:557` | 原报错里 "Inspect `consensus(x)`" 一句随分支一起消失（新报错已给出指引） |
| `NAMESPACE` | `export(consensus)` → `export(adjudicate)`；新增 `S3method(print, adjudication)` |
| `man/` | 由 roxygen 重新生成（`consensus.Rd` → `adjudicate.Rd`） |
| `tests/testthat/test-mechanism.R` | 412–467 行的 `consensus(` 全部改 `adjudicate(` |
| `NEWS.md` | 顶部小节的下一版条目（见 §8） |
| `DESCRIPTION` | 版本号递增（见 §8） |

## 7. 测试清单

新增 `tests/testthat/test-adjudicate.R`：

1. 返回对象继承 `adjudication`，且 `$verdicts` 列名与顺序等于冻结项。
2. 四类 verdict 各出现至少一次（用现有 `test-mechanism.R:412-467` 的固件迁移）。
3. `coverage_floor` 生效时，覆盖不足的来源整体报 `untested`。
4. `print()` 输出包含四类计数，且不含逐条 regulator（避免刷屏）。
5. 非 `evidenceBundle` 输入报错，消息含 "must be an `evidenceBundle`"。

新增 `tests/testthat/test-interpret-adjudication.R`：

6. `interpret(evidenceBundle)` 报错，消息含 "no longer accepts"（断言 §4.1 三层原文）。
7. `interpret(adjudication)` 的 prompt 只包含 `supported` 与 `single_source`；断言 prompt 字符串里出现对照声明那一句。
8. 同一 `adjudication` 里把某条 `supported` 人为改成 `single_source`，prompt 不再把该 regulator 放进机制段（用字符串断言，不跑模型）。

修改 `test-interpret-*.R`：确认四个老入口的快照**未变**（这一条是 D2 的安全网：删分支不能顺手改 prompt）。

## 8. 文档与发布

- roxygen：`adjudicate()` 的 `@details` 里写清它与 `interpret()` 的分工，并用一句话说明改名理由（"决定主张是否成立，而不只是找出各方一致之处"）。
- 类文档：`@name adjudication-class`，列出五个字段。
- `?interpret`：更新"输入类型"一节，明确三类输入与一句指引（"要断言机制，先 `adjudicate()`"）。
- `NEWS.md` 顶部小节（当前 `# clusterProfiler 4.21.4.007`）之上新增一条，一句话说清"改名 + 删兼容 + 新入口"三件事，并注明这是**未发布 API 的重命名**，不涉及兼容负担。
- `DESCRIPTION`：`4.21.4.007` → `4.21.4.008`（按本仓规则，开发期版本号每次提交 +1，正式号由作者决定）。

`NEWS.md` 条目草稿：

```markdown
+ `consensus()` is now `adjudicate()`, and `interpret()` no longer accepts an
  `evidenceBundle`. The function decides which claims survive a change of
  evidence type -- `supported`, `single_source`, `absent`, `untested` -- so
  "consensus" named the wrong operation, and a verdict object is something a
  caller can read, save and compare before paying for a model call. Adjudication
  is now an explicit step: `evidence()` collects, `adjudicate()` decides,
  `interpret()` writes. Both entry points take one step each and neither calls a
  model until the last. Verdict semantics, column order and the prompt of the
  four legacy inputs are unchanged (2026-10-28, Tue)
```

## 9. 验收判据

1. `grep -rn "consensus(" R/` 只剩 `consensus_findings(`（若保留该名）与历史注释；`interpret_dialogue.R`、`evidence.R`、`collapse.R`、`mechanism_graph.R` 四处调用全部改为 `adjudicate(`。
2. `interpret(evidenceBundle)` 报错，消息逐字匹配 §4.1。
3. `interpret(adjudication)` 的 prompt 不含 `untested` 实体的结论句（测试 7 覆盖）。
4. 四个老入口的快照测试全绿且**未修改快照文件**。
5. `make check` 通过（先读 `Makefile` 目标，不要自拼命令）。
6. 没有任何文件同时存在 `consensus(` 与 `adjudicate(` 两套用户入口。

## 10. 风险与回退

- **风险**：删 `interpret(evidenceBundle)` 会打断任何按 `mechanism-exit-spec.md` §4.7 写好的用法。该用法未发布，且报错已给出替代路径，接受。
- **风险**：把 verdict 措辞写在两处（`mechanism_writable()` 与 prompt 指令）会漂移。缓解：prompt 只引用 `writable` 列，指令表只写"不得升级"的规则。
- **回退**：全部改动放在一个提交里（`git mv` + 改名 + 删分支 + 测试），回退即 `git revert` 单一提交。

---

## 附：为什么值得为一个函数名做这次改动

报告里要讲的核心区分是"解说 vs 裁决"：解说所有人都在做，裁决是把多源证据折成一个带等级的判定。如果代码里两者共用一个入口、裁决藏在调用栈内部，这个区分在实现层就不存在；改名并让它返回一个对象之后，第 15–19 页那三页在代码里各有一个对应物，听众问"这个区别只是你们的说法吗"，答案是一个可运行的调用链。
