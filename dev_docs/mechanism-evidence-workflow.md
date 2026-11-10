# 机制证据工作流：发现、检验、披露，不自动认定因果

适用于正式基线4.21.5之后的开发修复。这里的输出是**受证据约束的候选解释**，不是唯一上游驱动者的识别结果；不加入靶基因数量权重。

## 1. 候选发现与通路检验

下面变量由调用者提供：`signed_all`是按gene symbol命名的完整测量log2FC向量，`ego`是富集结果，`tf_network`含`source`、`target`、可选`weight`。

```r
# 在看排名之前按实验设计确定排序方式，不按预期标签挑选较好结果。
shortlist <- mechanism_discover(signed_all, tf_network, top = 5,
                               rank_by = "score") # KO/KD可预先选择"abs"
pf <- mechanism_preflight(ego, signed = signed_all,
                          drivers = shortlist$source)
pe <- pathway_evidence(ego, drivers = shortlist$source, signed = signed_all)
r <- mechanism_report(pe, preflight = pf, network = tf_network)

r$verdicts
r$candidate_evidence
r$perturbation_tests
r$gaps
```

- 预检是当前输入和覆盖范围的启发式检查，不是因果充分条件。`description`不等于未来实验无用。
- 数据库查询失败是覆盖未知，不是没有数据。若传`inventory`，报告只对该表范围作陈述；保存表与原始输入才能脱离实时数据库复现。
- A级要求同一候选有足够方向一致、显著的细胞系证据，并通过成员支持度要求。缺失/空白细胞系身份不能计作独立复现；两个不同候选各自一条证据不能拼成一个A级候选。
- 为兼容保留的`coherent-refuted`名称只表示所测数据没有显著支持，不是证明候选为假。报告分类为`supported`、`conflicting`、`undetermined`，不再使用绝对的`excluded`。
- 原有最佳结果行不是证据全集。一个A级行仍可能在别的数据集有反向结果；报告保留等级以避免增加未经验证的新规则，同时在`candidate_evidence$assessment == "mixed"`及逐条记录中明确披露冲突。不能把A级直接改写成“无反证”或“已确证驱动”。

`perturbation_tests`的状态含义：

| 状态 | 所测数据允许的解释 |
|---|---|
| `supports` | 显著且扰动ES与该通路方向相反，符合当前方向解释；不单独保证成员支持度或因果身份 |
| `conflicts` | 显著但方向不符合当前解释；不能据此排除该因子的所有生物学作用 |
| `not-significant` | 未达到当前阈值；不是等效性检验或无效证明 |
| `failed` | 调用失败，未取得可解释统计量 |
| `not-computable` | 缺方向、缺统计量、非有限值或零ES，不能作该方向判定 |

未能解析到数据集的候选没有可虚构的dataset行；通过保留的请求候选及`gaps`披露覆盖限制。保留属性随普通R对象保存；将结果导出CSV会丢失属性，应同时保存完整RDS或分别保存证据表。

`module`只是请求候选之间的无符号、两两Jaccard邻域，不是互斥分区，不证明生物学等价、无法区分或上下游顺序。`global_effect_reason`只是所提供/查得扰动的DE数量中位数过阈值，不证明必需性或该通路的非特异性；两者均不重排候选或自动改变等级。

## 2. 先裁决，再交给语言模型

这是另一条证据接口，不应把通路级报告冒充`adjudication`。`signed_de`是调用者确定的有符号DE基因，`measured_genes`是测量背景。

```r
bundle <- evidence(signed_de, universe = measured_genes,
                   upstream = c("tf", "perturbagen"))
a <- adjudicate(bundle)
a$verdicts
# 模型与服务配置就绪后才运行；会调用所配置的模型。
interpret(a)
```

裁决对象记录实际参数与证据指纹。显式`coverage_floor`覆盖默认值，并影响实际筛选；改变证据p值应改变指纹，而仅改变行列顺序或运行时间不应改变。共享同一dataset的不同来源不自动构成两个独立见证。

传给解释器的结构化授权记录精确候选身份、裁决、方向、扰动效应与允许措辞。`single_source`不能被写成多来源确证；`supported`也不能写成唯一因果上游。只有扰动证据而没有富集对象的合法输入仍应保留。

**这是输入约束和可审查记录，不是自然语言输出的形式化验证器。** 本次没有用实时候选全集重跑原历史验证，也没有评估不同模型的遵从率；不能据测试通过宣称模型不会越界。

## 3. 重现实证基准与验收

[固定17份基准的方法、结果和边界](<mechanism-benchmark-4215.md>)与[运行脚本](<mechanism-benchmark-4215.R>)明确区分完整发现排名、历史扰动验证、集合覆盖及大小匹配基线。全部17份必须保留在分母；重复study/driver分组披露。候选标签来自实验条件，不是排他的因果金标准。

包验证遵循项目Makefile：先修改R中的roxygen，再`make rd`生成文档，最终`make check`。离线回归检验的是确定性逻辑边界；真实17份排名基准检验的是这一固定探索性语料上的观察，二者都不替代外部独立因果验证。
