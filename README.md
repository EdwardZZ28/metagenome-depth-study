# How deep is deep enough? 污水宏基因组测序深度对 AMR 指标的影响

**Working report / 阶段性工作记录 — 2026-10-11**
Status: core set complete (32 samples, 5 matrices, 320 datasets, read-based and assembly-based
metrics). Bin quality (CheckM2, GTDB-Tk) and the two larger validation layers are still running.
状态:核心组已完成(32 个样本、5 种基质、320 个数据集,读长与组装两端指标)。bin 质量
(CheckM2、GTDB-Tk)和两个更大的验证层仍在运行。

---

## 1. The question / 研究问题

**EN.** Shotgun metagenomics is now routine in environmental AMR surveillance, but there is no
agreed sequencing depth for wastewater. Published studies use anywhere from 5 to 30 Gbp per
sample and almost none report whether their conclusions would change at a different depth. The
gap matters most for assembly-derived results (plasmids, mobile genetic elements, metagenome
assembled genomes), which are the fastest-growing part of the field.

This study answers the question the only way it can be answered cleanly: take one sample, split
its own reads into a depth ladder, run an identical pipeline at every rung, and see which numbers
move. Because every comparison is within a sample, library preparation, extraction chemistry,
sequencing batch and community composition are all held constant. Depth is the only variable.

**中文.** 鸟枪法宏基因组已经是环境 AMR 监测的常规手段,但污水样本该测多深没有共识。已发表研究
的深度从 5 到 30 Gbp 不等,几乎没有人报告"换个深度结论会不会变"。这个缺口对组装类结果
(质粒、可移动遗传元件、宏基因组组装基因组)影响最大,而这恰恰是近年增长最快的方向。

本研究用唯一干净的方式回答这个问题:取同一个样本,把它自己的 reads 抽成一个深度梯度,每一档
跑完全相同的流程,看哪些数字会变。因为所有比较都在样本内部进行,建库、提取、测序批次和群落组成
全部固定,唯一变化的就是深度。

---

## 2. Design / 设计

| | |
|---|---|
| Samples / 样本 | 32 metagenomes with clean depth >= 25 Gbp: 17 wastewater influent, 1 effluent, 6 preserved wastewater, 4 river water, 4 soil |
| Origin / 来源 | public SRA/ENA runs plus in-house Australian samples; influent draws on 17 BioProjects |
| Depth rungs / 深度档位 | 2, 5, 10, 15, 20, 25 Gbp of clean paired reads |
| Seeds / 随机种子 | 3 independent seqtk seeds at 2 and 5 Gbp, 1 elsewhere |
| Datasets / 数据集总数 | 320 |
| Reference / 参照 | the 25 Gbp rung of the same sample |

Sections 3.1 to 3.4 report the first 11 influent samples (66 datasets), where no more than two
samples came from any one BioProject. Sections 3.5 and 3.6 report the full core set.
第 3.1 到 3.4 节基于最初的 11 个进水样本(66 个数据集),当时每个 BioProject 最多取 2 个样本;
第 3.5 和 3.6 节基于完整的核心组。

The full pipeline at every rung: seqtk subsampling, ARGs-OAP (read-based ARG abundance),
Kraken2/Bracken (taxonomy), MEGAHIT (assembly), RGI/CARD (contig-level ARGs),
geNomad (plasmids, viruses, conjugative elements), Bowtie2 + MetaBAT2 (binning).
See [METHODS.md](METHODS.md) for versions, parameters and exactly how each number is counted.

---

## 3. Headline results / 主要结果

### 3.1 Abundance is immune to depth; everything else is not
### 3.1 丰度对深度免疫,其余指标都不是

![recovery curves](figures/fig1_recovery.png)

Median recovery relative to the same sample at 25 Gbp:
相对于同一样本 25 Gbp 的中位回收率:

| Depth (Gbp) | ARG abundance (copies/16S) | ARG subtype richness | RGI hits (id >=80%) | Contigs | Plasmid contigs | Conjugative | Virus contigs | Bins |
|---|---|---|---|---|---|---|---|---|
| 2  | **100%** | 51% | 21% | 11% | 11% | 11% | 5%  | 14% |
| 5  | **100%** | 66% | 43% | 28% | 26% | 29% | 19% | 29% |
| 10 | **100%** | 80% | 61% | 51% | 49% | 50% | 43% | 51% |
| 15 | **100%** | 89% | 79% | 68% | 68% | 68% | 62% | 72% |
| 20 | **100%** | 95% | 93% | 87% | 85% | 84% | 84% | 91% |

**EN.** ARG abundance normalised to 16S is flat to two decimal places across a 12.5-fold range of
depth. Every count-based metric, by contrast, is still climbing at the top of the ladder. At the
depth most published wastewater studies actually use, about 10 to 14 Gbp, roughly **half** of the
plasmid contigs, conjugative elements and genome bins that 25 Gbp would have found are simply
absent. Note that these percentages are relative to 25 Gbp, which is a budget ceiling rather than
a true value; section 3.7 restates the same data as marginal return per Gbp, which does not depend
on that denominator.

**中文.** 以 16S 归一化的 ARG 丰度在 12.5 倍的深度范围内小数点后两位都不变。相比之下,所有
计数型指标到梯度顶端仍在上升。在已发表污水研究实际常用的深度(约 10 到 14 Gbp)下,25 Gbp
本可以找到的质粒 contig、接合性元件和基因组 bin,**大约有一半根本不会出现**。注意这些百分比
是相对 25 Gbp 而言的,而 25 Gbp 是预算上限而非真值;3.7 节用每 Gbp 边际收益重新表述同一批
数据,那个指标不依赖这个分母。

### 3.2 How much 10 Gbp recovers depends on the sample
### 3.2 10 Gbp 能回收多少,取决于样本本身

![complexity](figures/fig2_complexity.png)

| Sample | Contigs at 25 Gbp | ARG richness recovered at 10 Gbp | Contigs recovered at 10 Gbp |
|---|---|---|---|
| SRR12147797 | 37,029 | 90% | **85%** |
| ERR13731085 | 106,121 | 83% | 64% |
| ERR12510812 | 173,031 | 78% | 52% |
| ERR13731241 | 261,059 | 83% | 59% |
| SRR35672372 | 269,944 | 80% | 41% |
| SRR27197501 | 235,491 | 82% | **33%** |
| SRR35242016 | 552,726 | 78% | 39% |

**EN.** ARG subtype recovery at 10 Gbp is tightly clustered between 75% and 90% in all 11 samples,
so that part of the answer generalises. Contig recovery is not: it ranges from 33% to 85% and
tracks community complexity. A simple, cheap proxy for complexity (contig count from a shallow
run) predicts how much deeper a sample needs to go. That is more useful than a single
recommended number, because it allows a two-pass design: sequence everything shallow, then
top up only the complex samples.

**中文.** 11 个样本在 10 Gbp 的 ARG 种类回收率全部落在 75% 到 90% 之间,这一条可以推广。
contig 回收率则不然,从 33% 到 85%,并且和群落复杂度相关。复杂度可以用浅测序的 contig 数
便宜地估计出来,从而预测某个样本还需要测多深。这比给一个固定推荐值更有用,因为它支持两段式
设计:先全部浅测,再只对复杂样本补测。

### 3.3 Low depth misses the genes that matter most
### 3.3 低深度漏掉的恰恰是最要紧的基因

![high-risk](figures/fig3_highrisk.png)

| Depth (Gbp) | High-risk gene detections | Samples with >=1 detection |
|---|---|---|
| 2  | 3  | 3 / 11 |
| 5  | 8  | 6 / 11 |
| 10 | 16 | 7 / 11 |
| 25 | **22** | 7 / 11 |

Genes recovered only at higher depth include **NDM-1, IMP-1, KPC-75, CTX-M-15, CTX-M-30,
MCR-4.3 and mcr-10.1**.

**EN.** Carbapenemases and mobile colistin resistance genes are exactly what a surveillance
programme exists to detect, and they are rare enough in wastewater that shallow sequencing
misses them. A programme that reports only normalised ARG abundance would conclude its depth
was adequate, while silently under-reporting carbapenemase detections by roughly 30% at 10 Gbp.

**中文.** 碳青霉烯酶和可移动黏菌素耐药基因正是监测项目最想抓到的目标,而它们在污水里足够稀有,
浅测序会漏掉。一个只报告归一化 ARG 丰度的项目会认为自己的深度够用,而实际上在 10 Gbp 下
碳青霉烯酶的检出会被低估约三成。

---

### 3.4 What is lost is not a random sample of what was there
### 3.4 丢掉的不是随机的一部分

Plasmid contigs called at 25 Gbp were used as a reference set and searched for in the assemblies
of the same sample at lower depth (minimap2, alignment identity >= 95%). Two thresholds are
reported: **strict** (alignment covers >= 80% of the reference contig, i.e. the plasmid was
rebuilt) and **lenient** (>= 30%, i.e. the plasmid was seen at all). Full output in
[`data/plasmid_fate_summary.txt`](data/plasmid_fate_summary.txt).

以 25 Gbp 检出的质粒 contig 为参照集,在同一样本较低深度的组装结果中搜索(minimap2,比对
相似度 >= 95%)。报告两个判定标准:**严格**(比对覆盖参照 contig 的 80% 以上,即质粒被完整
重建)和**宽松**(覆盖 30% 以上,即质粒被看到)。完整输出见
[`data/plasmid_fate_summary.txt`](data/plasmid_fate_summary.txt)。

Recovery at 10 Gbp, median across the 11 samples / 10 Gbp 时 11 个样本的中位回收率:

| Group / 分组 | n contigs | Strict / 严格 | Lenient / 宽松 |
|---|---|---|---|
| 1-2 kb | 8,645 | 6% | 9% |
| 2-5 kb | 21,908 | 7% | 15% |
| 5-10 kb | 9,633 | 13% | 26% |
| > 10 kb | 4,282 | **17%** | **39%** |
| carries an ARG / 携带 ARG | 290 | **19%** | **40%** |
| no ARG / 不携带 ARG | 44,178 | 10% | 18% |
| conjugative / 接合性 | 8,362 | 8% | 21% |
| non-conjugative / 非接合性 | 36,106 | 10% | 18% |

**EN.** Two things stand out. First, long plasmids survive shallow sequencing far better than
short ones: at 10 Gbp a plasmid over 10 kb is four times more likely to be detected than one of
1 to 2 kb. Second, and against expectation, **plasmids carrying an ARG are recovered about twice
as often as plasmids without one** (40% vs 18% lenient, 19% vs 10% strict), consistently at every
rung. The likely reason is not the ARG itself but its host: acquired ARGs sit disproportionately
on plasmids of abundant organisms, which are covered more deeply and assemble more easily.

This has a specific and awkward consequence. The numerator of "fraction of plasmids carrying an
ARG" falls more slowly with depth than the denominator, so **a shallow library inflates that
fraction**. At 10 Gbp the denominator retains 18% of reference plasmids while the numerator
retains 40%, so the estimated ARG-carriage rate is roughly doubled. Mobilome studies that report
plasmid-borne ARG proportions from shallow libraries are therefore biased upwards, not downwards,
even though they see fewer plasmids overall. Conjugation genes, by contrast, make almost no
difference to recovery (21% vs 18%), so conjugative elements are not unusually hard to assemble.

**中文.** 有两点值得注意。第一,长质粒在浅测序下的存活率远高于短质粒:10 Gbp 时,大于 10 kb
的质粒被检出的概率是 1 到 2 kb 质粒的四倍。第二,与预期相反,**携带 ARG 的质粒回收率约为不
携带者的两倍**(宽松标准 40% 对 18%,严格标准 19% 对 10%),且在每一档都保持这个方向。
原因很可能不是 ARG 本身,而是它的宿主:获得性 ARG 不成比例地位于高丰度菌的质粒上,这些质粒
覆盖度更高,也更容易组装。

这带来一个具体而别扭的后果。"携带 ARG 的质粒占比"这个指标,分子随深度下降的速度慢于分母,
所以**浅测序会高估这个占比**。10 Gbp 时分母保留 18%,分子保留 40%,估算出的携带率大约被
放大一倍。也就是说,基于浅测序文库报告质粒携带 ARG 比例的 mobilome 研究,偏差方向是偏高
而非偏低,尽管它们看到的质粒总数更少。相比之下,接合基因对回收率几乎没有影响(21% 对 18%),
说明接合性元件并不特别难组装。

**Caveat / 注意.** Only 290 of 44,468 reference plasmid contigs carry an ARG at identity >= 80%,
so the ARG-carriage comparison rests on a small numerator and has not yet been tested for
per-sample consistency.
44,468 条参照质粒 contig 中只有 290 条在 identity >= 80% 下携带 ARG,因此这一组比较的分子
很小,尚未检验其在各样本间是否稳定。

---

### 3.5 The depth a sample needs depends on what it is
### 3.5 一个样本需要多深,取决于它是什么

The ladder was extended to 32 samples across five matrices, with three subsampling seeds at the
2 and 5 Gbp rungs. Median recovery relative to the same sample at 25 Gbp:

梯度扩展到五种基质共 32 个样本,2 和 5 Gbp 档各用三个随机种子。相对同一样本 25 Gbp 的中位
回收率:

![by matrix](figures/fig4_by_matrix.png)

| Matrix / 基质 | n | Contigs @2 | @5 | @10 | @15 | Plasmid @10 |
|---|---|---|---|---|---|---|
| Wastewater influent / 污水进水 | 17 | 10% | 26% | 50% | 67% | 48% |
| Wastewater effluent / 污水出水 | 1 | 6% | 18% | 41% | 62% | 53% |
| Preserved wastewater / 加保存剂污水 | 6 | 6% | 18% | 43% | 65% | 39% |
| River water / 河水 | 4 | 5% | 17% | 39% | 62% | 26% |
| **Soil / 土壤** | 4 | **0%** | **1%** | **12%** | **35%** | **8%** |

**EN.** ARG abundance normalised to 16S holds at 95% to 104% of the full-depth value in every
matrix at every rung, so that result generalises beyond wastewater. Nothing else does. Soil is a
different problem entirely: at 5 Gbp its assembly recovers 1% of the contigs and 1% of the plasmid
contigs that 25 Gbp finds, and produces no genome bins at all. River water sits between soil and
wastewater, and is noticeably worse for plasmids (26% at 10 Gbp against 48% for influent).
Taxonomic profiling is far more robust than assembly in all matrices: Bracken species counts
reach 78% to 91% at 10 Gbp where contigs reach 12% to 50%.

The three wastewater forms (influent, effluent, and influent held with preservatives) give almost
identical curves, which is a useful negative result: adding a preservative does not change how
deep the sample needs to be sequenced.

**中文.** 以 16S 归一化的 ARG 丰度在所有基质、所有档位都维持在满深度值的 95% 到 104%,所以
这一条可以推广到污水以外。其余指标都不行。土壤完全是另一回事:5 Gbp 时它的组装只能回收
25 Gbp 所得 contig 的 1%、质粒 contig 的 1%,而且一个基因组 bin 都产不出。河水介于土壤和
污水之间,质粒回收尤其差(10 Gbp 时 26%,进水为 48%)。在所有基质中,物种分类都比组装稳健
得多:10 Gbp 时 Bracken 物种数达到 78% 到 91%,而 contig 只有 12% 到 50%。

污水的三种形态(进水、出水、加保存剂的进水)曲线几乎重合,这是一个有用的阴性结果:加保存剂
不改变该样本需要测多深。

### 3.6 Subsampling noise is small / 抽样随机性很小

![seed noise](figures/fig5_seed_noise.png)

Three independent seqtk seeds were run at the 2 and 5 Gbp rungs of every core sample. Median
coefficient of variation across seeds: ARG abundance 1.9%, ARG richness 2.7%, contigs 0.4%,
plasmid contigs 3.5%. The curves above are therefore not an artefact of which reads happened to
be drawn. The only large coefficients of variation come from counts near zero, which is itself the
soil and river-water result restated: at 2 Gbp those samples yield zero, one or two plasmid
contigs, so a 173% coefficient of variation describes the values 0, 0, 2.

每个核心样本的 2 和 5 Gbp 档都用三个独立的 seqtk 随机种子跑过。种子间变异系数中位数:
ARG 丰度 1.9%,ARG 种类 2.7%,contig 0.4%,质粒 contig 3.5%。所以上面的曲线不是"碰巧抽到
哪些 reads"造成的。唯一出现的大变异系数都来自接近零的计数,而那本身就是土壤和河水结果的另
一种说法:2 Gbp 时这些样本只产出 0、1 或 2 条质粒 contig,173% 的变异系数描述的是 0、0、2
这组数。

---

### 3.7 The useful question is not "is it saturated" but "what does the next Gbp buy"
### 3.7 有用的问题不是"饱和了没有",而是"下一个 Gbp 买到什么"

Recovery percentages depend on the 25 Gbp denominator, and 25 Gbp is a ceiling set by what is
affordable, not a true value. Marginal return does not have that problem: for each sample the gain
between two adjacent rungs is divided by the Gbp between them. Full output in
[`data/marginal_returns.txt`](data/marginal_returns.txt).

回收率百分比依赖 25 Gbp 这个分母,而 25 Gbp 是按可负担程度设定的上限,不是真值。边际收益没有
这个问题:对每个样本,取相邻两档之间的增量除以两档之间的 Gbp 数。完整输出见
[`data/marginal_returns.txt`](data/marginal_returns.txt)。

Median marginal gain per extra Gbp, wastewater influent (17 samples):
污水进水(17 个样本)每多测 1 Gbp 的中位边际收益:

| Metric / 指标 | 2→5 | 5→10 | 10→15 | 15→20 | 20→25 |
|---|---|---|---|---|---|
| ARG abundance (copies/16S) | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| ARG subtypes | 40.7 | 18.6 | 12.6 | 9.6 | 8.0 |
| Bracken species | 9,733 | 4,722 | 2,120 | 1,126 | 701 |
| Contigs | 10,249 | 10,340 | 9,325 | 8,305 | 7,515 |
| Plasmid contigs | 230 | 179 | 147 | 156 | 157 |
| Virus contigs | 187 | 182 | 192 | 202 | 194 |
| Bins | 8.3 | 6.6 | 5.4 | 5.2 | 5.0 |

Deceleration, defined as the 20→25 gain divided by the 2→5 gain. A value near zero means the
metric has levelled off inside the window; near one means it is still climbing at the same rate;
above one means it is accelerating.
减速比定义为 20→25 的增益除以 2→5 的增益。接近 0 表示该指标在窗口内已经趋平;接近 1 表示
仍以同样速率上升;大于 1 表示还在加速。

| Matrix / 基质 | copies/16S | ARG subtypes | species | contigs | plasmid | virus | bins |
|---|---|---|---|---|---|---|---|
| Influent / 进水 | 0.03 | 0.20 | 0.07 | 0.58 | 0.58 | 0.77 | 0.60 |
| Effluent / 出水 | -0.14 | 0.19 | 0.03 | 0.88 | 0.44 | 1.01 | 1.10 |
| Preserved / 加保存剂 | 0.12 | 0.23 | 0.16 | 0.86 | 1.02 | 0.42 | 1.01 |
| River water / 河水 | -0.03 | 0.22 | 0.02 | 1.02 | 1.69 | 1.20 | 0.45 |
| **Soil / 土壤** | 0.08 | 0.22 | 0.10 | **16.69** | **14.64** | **3.81** | **3.60** |

**EN.** Three behaviours, not one. Abundance, ARG richness and taxonomic richness flatten inside
the window: in influent the species count buys 9,733 species per Gbp at the 2→5 step and 701 at
20→25, a deceleration of 0.07. Assembly-derived counts do not flatten: contigs, plasmid contigs
and bins still return 58% to 60% of their initial rate at the top of the ladder, and virus contigs
return 77% in influent and 101% in effluent, meaning no deceleration at all. In that regime every
extra Gbp buys about the same amount, so the decision is purely a budget question: roughly 7,500
contigs, 157 plasmid contigs and 5 bins per Gbp in influent.

Soil is a third shape. Its deceleration ratios are far above one (contigs 16.7, plasmid contigs
14.6), which means the curve is convex: the 20→25 step returns more than ten times what the 2→5
step returned. Soil assembly has not reached its productive range by 25 Gbp, and spending
anything below about 10 Gbp on soil assembly returns almost nothing. This is a different statement
from "soil needs more depth": the shape of the curve, not only its height, differs from water.

**中文.** 是三种行为,不是一种。丰度、ARG 种类和物种丰富度在窗口内趋平:进水的物种数在 2→5
段每 Gbp 买到 9,733 个物种,到 20→25 段只剩 701 个,减速比 0.07。组装类计数不趋平:contig、
质粒 contig 和 bin 在梯度顶端仍保持初始速率的 58% 到 60%,病毒 contig 在进水是 77%、出水
101%,也就是完全没有减速。在这种情形下,每多测 1 Gbp 买到的东西大致相同,所以决策纯粹是预算
问题:进水每 Gbp 约 7,500 个 contig、157 个质粒 contig、5 个 bin。

土壤是第三种形状。它的减速比远大于 1(contig 16.7,质粒 contig 14.6),意味着曲线是凸的:
20→25 这一段的收益是 2→5 段的十倍以上。土壤组装到 25 Gbp 还没有进入高产区间,而在约 10 Gbp
以下做土壤组装几乎得不到任何东西。这和"土壤需要更深"是不同的陈述:和水样相比,不同的不只是
曲线的高度,还有它的形状。

---

### 3.8 Shallow data does not degrade every conclusion equally. One kind of conclusion is wrong, not merely incomplete
### 3.8 浅测序并非均匀地降低所有结论的质量,有一类结论是错的,而不只是不完整

Everything above measures how much is missing. This section asks a different question: taking each
sample on its own, does the **answer** change with depth? Three checks against the same sample at
25 Gbp. Full output in
[`data/conclusion_stability.txt`](data/conclusion_stability.txt).

上面的内容衡量的都是"缺了多少"。这一节问的是另一个问题:就每个样本自身而言,**答案**会不会
随深度改变?三项检验,都与同一样本的 25 Gbp 结果比较。完整输出见
[`data/conclusion_stability.txt`](data/conclusion_stability.txt)。

**Check 1, species composition.** Spearman correlation with the 25 Gbp profile at 2 Gbp: 0.84 to
0.94 across the five matrices, rising to 0.96 to 0.98 by 10 Gbp. More directly: in every matrix,
**all ten of the top ten species at 25 Gbp are already the top ten at 2 Gbp**, with no exceptions
among the 32 samples.

**检验一,物种组成。** 2 Gbp 与 25 Gbp 的 Spearman 相关在五种基质中为 0.84 到 0.94,到 10 Gbp
升至 0.96 到 0.98。更直接的是:在每一种基质中,**25 Gbp 的前十名物种在 2 Gbp 时已经全部位列
前十**,32 个样本无一例外。

**Check 2, ARG drug-class composition.** Correlation at 2 Gbp is 0.95 to 0.99; the top five drug
classes are the same five, except occasionally in river water where contig-level ARGs are scarce
to begin with.

**检验二,ARG 药物类别组成。** 2 Gbp 时相关系数 0.95 到 0.99;前五名药物类别就是那五个,
只有河水偶尔差一个,而河水本身 contig 层面的 ARG 就极少。

**Check 3, the fraction of plasmid contigs that carry an ARG.** This one does change, and in a
consistent direction:

**检验三,携带 ARG 的质粒 contig 占全部质粒 contig 的比例。** 这一项确实变了,而且方向一致:

| Matrix / 基质 | 2 Gbp | 5 | 10 | 15 | 20 | 25 Gbp | Inflation at 2 Gbp |
|---|---|---|---|---|---|---|---|
| Influent / 进水 | 0.98% | 1.01% | 0.70% | 0.63% | 0.67% | 0.62% | **1.6x** |
| Effluent / 出水 | 0.78% | 0.17% | 0.35% | 0.32% | 0.31% | 0.37% | **2.1x** |
| Preserved / 加保存剂 | 0.46% | 0.32% | 0.22% | 0.31% | 0.30% | 0.21% | **2.2x** |

**EN.** The explanation is in section 3.4: plasmids carrying an ARG are recovered about twice as
often at low depth as plasmids without one, because acquired ARGs sit disproportionately on
plasmids of abundant organisms. The numerator of this ratio therefore falls more slowly with depth
than the denominator, and a shallow library inflates the estimate by roughly a factor of two.

This matters more than any of the count results. A missing contig announces itself as missing. A
ratio does not: 0.98% looks like a measurement, carries no warning that it is twice the value a
deeper library would give, and would be reported, plotted and compared across studies exactly as
if it were correct. Mobilome studies that compare plasmid-borne ARG fractions between sites,
seasons or countries are comparing a quantity that moves with sequencing depth, and depth is
rarely matched across such comparisons.

Taken together, the three checks give a simple rule. **Composition survives shallow sequencing;
counts are incomplete but honest about it; ratios between two depth-sensitive quantities are
biased and give no sign of it.**

**中文.** 解释在 3.4 节:携带 ARG 的质粒在低深度下的回收率约为不携带者的两倍,因为获得性
ARG 不成比例地位于高丰度菌的质粒上。于是这个比值的分子随深度下降得比分母慢,浅测序文库会把
估计值放大约一倍。

这比前面任何计数结果都重要。缺失的 contig 会自己暴露出来,比例不会:0.98% 看起来就是一个
测量值,不带任何提示说明它是更深文库所得数值的两倍,而且会被原样报告、作图,并与其他研究
比较。那些在站点、季节或国家之间比较质粒携带 ARG 比例的 mobilome 研究,比较的是一个随测序
深度变化的量,而这类比较很少对深度做过匹配。

三项检验合起来给出一条简单规则:**组成能经受浅测序;计数不完整,但它对此是诚实的;两个
depth-sensitive 量之间的比值则有偏,且不给出任何提示。**

---

## 4. What this means in practice / 实践含义

**EN.**

1. If the deliverable is ARG abundance or drug-class composition, 5 to 10 Gbp is sufficient and
   spending more is waste.
2. If the deliverable is a list of which ARGs are present, 10 Gbp returns about 80% of the
   subtypes and the curve is still rising.
3. If the deliverable involves plasmids, mobility or MAGs, there is no depth inside this window at
   which the return stops. Every extra Gbp of influent buys about 7,500 contigs, 157 plasmid
   contigs and 5 bins, at a rate that barely falls between 2 and 25 Gbp, so the choice is a budget
   decision rather than a saturation point. Published mobilome and MAG results from 10 to 15 Gbp
   libraries are systematic underestimates of what a deeper library would have found.
4. For surveillance of clinically high-risk genes, depth is not a tuning parameter. It decides
   whether a carbapenemase is reported at all.
5. Do not carry a wastewater depth across to soil. At 5 Gbp a soil assembly returns 1% of what
   25 Gbp returns and no genome bins, so assembly-based work on soil below 10 Gbp produces
   almost nothing. River water needs more depth than wastewater for plasmids specifically.
6. Where a project only needs who is present rather than what they carry, taxonomic profiling
   survives shallow sequencing in every matrix tested.

**中文.**

1. 如果产出是 ARG 丰度或药物类别组成,5 到 10 Gbp 就够,再多是浪费。
2. 如果产出是"检出了哪些 ARG"的清单,10 Gbp 能拿到约八成 subtype,曲线仍在上升。
3. 如果产出涉及质粒、移动性或 MAG,在这个窗口内不存在"收益停止"的深度。进水每多测 1 Gbp
   约买到 7,500 个 contig、157 个质粒 contig 和 5 个 bin,而且这个速率在 2 到 25 Gbp 之间
   几乎不下降,所以这是预算决策而非饱和点。基于 10 到 15 Gbp 的已发表 mobilome 和 MAG 结果
   是对更深文库所能发现内容的系统性低估。
4. 对高风险基因的监测而言,深度不是一个可调参数,它直接决定碳青霉烯酶会不会被报出来。
5. 不要把污水的深度标准搬到土壤。5 Gbp 时土壤组装只能得到 25 Gbp 所得的 1%,且没有任何
   基因组 bin,所以低于 10 Gbp 做土壤的组装类分析几乎得不到东西。河水在质粒这一项上需要
   比污水更深。
6. 如果项目只需要知道有哪些物种、而不是它们携带什么,物种分类在所有测试过的基质里都能
   经受浅测序。

---

## 5. Known limitations / 已知问题

**EN.**

- **25 Gbp is a budget ceiling, not a true value.** Every recovery percentage in sections 3.1 to
  3.5 is relative to the 25 Gbp rung of the same sample, and section 3.7 shows that assembly
  metrics are still climbing there. True recovery against the full community is therefore lower
  than every percentage quoted. Marginal return (section 3.7) is the denominator-free alternative.
- **Effluent rests on one sample.** Only one effluent metagenome reached 25 Gbp, so the effluent
  row of section 3.5 is a single curve, not a median. Eight more effluent samples are queued at a
  15 Gbp ceiling in the next layer.
- **Soil and river water have four samples each**, all from one flood study, so matrix and study
  are confounded for those two rows.
- **MetaBAT2 returned zero bins** for two datasets (R8188_S1 at 25 Gbp, SRR23730275 at 10 Gbp)
  and those are excluded from bin statistics.
- **Bin quality is not yet assessed.** CheckM2 and GTDB-Tk have not been run, so "bins" here
  means raw MetaBAT2 output, not high-quality MAGs. The genuinely interesting number, how many
  near-complete genomes each depth yields, is still missing.
- **Influent only.** Whether the same curves hold for effluent, sludge, river water, soil and air
  is the subject of Layer 2.
- **geNomad plasmid and virus counts are contig counts**, not de-replicated element counts. They
  scale with assembly size by construction, which is part of what is being measured, but it means
  the absolute values are not comparable to studies that count elements differently.

**中文.**

- **25 Gbp 是预算上限,不是真值。** 3.1 到 3.5 节的所有回收率都相对于同一样本的 25 Gbp 档,
  而 3.7 节表明组装类指标在那里仍在上升。因此相对真实群落的回收率比文中引用的每一个百分比
  都更低。边际收益(3.7 节)是不依赖分母的替代指标。
- **出水只有 1 个样本。** 只有一个出水宏基因组达到 25 Gbp,所以 3.5 节出水那一行是单条曲线
  而非中位数。下一层会补 8 个出水样本,上限 15 Gbp。
- **土壤和河水各只有 4 个样本**,且全部来自同一项洪水研究,因此这两行里基质和研究是混淆的。
- **两个数据集 MetaBAT2 产出 0 个 bin**(R8188_S1 的 25 Gbp、SRR23730275 的 10 Gbp),
  已从 bin 统计中排除。
- **尚未评估 bin 质量。** CheckM2 和 GTDB-Tk 还没跑,所以这里的 bins 是 MetaBAT2 的原始
  输出,不是高质量 MAG。真正有意思的数字(每个深度能得到多少接近完整的基因组)还没有。
- **只做了进水。** 出水、污泥、河水、土壤和空气是否遵循同样的曲线,属于第二层的内容。
- **geNomad 的质粒和病毒计数是 contig 计数**,不是去冗余后的元件计数。它按定义会随组装规模
  增长,这本身正是要测量的对象,但这意味着绝对值不能和采用其他计数方式的研究直接比较。

---

## 6. Next steps / 下一步

1. **Finish the running jobs.** 168 datasets take effluent from 1 to 8 samples and add 4 air
   samples, both capped at 15 Gbp; GTDB-Tk is classifying the 1,627 high-quality bins.
   完成正在跑的任务:168 个数据集把出水从 1 个扩到 8 个并加入 4 个空气样本(上限 15 Gbp);
   GTDB-Tk 正在对 1,627 个高质量 bin 做分类。
2. **Extend soil to 30 and 35 Gbp.** Section 3.7 shows the soil curve is convex and still
   accelerating at 25 Gbp, so the productive range begins somewhere above the current ceiling.
   The four soil samples have 31 to 40 Gbp of clean data, enough for two more rungs.
   把土壤扩展到 30 和 35 Gbp:3.7 节显示土壤曲线是凸的,到 25 Gbp 仍在加速,高产区间开始于
   当前上限之上。四个土壤样本有 31 到 40 Gbp 的 clean 数据,够再加两档。
3. **Report MAG yield as marginal return**, i.e. how many additional classifiable high-quality
   genomes each extra Gbp buys, rather than recovery against 25 Gbp.
   把 MAG 产出改用边际收益报告,即每多测 1 Gbp 多得几个可分类的高质量基因组。
4. **Repurpose the third layer.** It was designed to test recovery ratios at scale, but recovery
   needs a 25 Gbp denominator that those 171 samples do not have. Better use: sequence each at
   2 and 10 Gbp and test whether the contig count at 2 Gbp predicts what 10 Gbp will yield, which
   is the practical form of the complexity result in section 3.2.
   第三层改变用途:它原本用来在更大规模上检验回收率,但回收率需要 25 Gbp 分母,而那 171 个
   样本没有。更好的用法是每个样本只跑 2 和 10 Gbp 两档,检验 2 Gbp 的 contig 数能否预测
   10 Gbp 的产出,这正是 3.2 节复杂度结果的可操作形式。
5. **Per-sample curves as the main figure.** Each of the 32 samples is an independent within-sample
   experiment, so plotting all 32 curves is stronger evidence than a median.
   把每个样本的曲线作为主图:32 个样本各自都是一次独立的样本内实验,画出 32 条曲线比中位数
   更有说服力。

## 7. Repository layout / 仓库结构
data/L1_summary.tsv one row per dataset, all extracted metrics / 每个数据集一行的全部指标
figures/ the three figures in this report / 报告中的三张图
scripts/summarise_L1.py builds data/L1_summary.tsv from pipeline output / 从流程输出生成汇总表
scripts/plot_L1.py builds figures/ from data/L1_summary.tsv / 从汇总表生成图
METHODS.md versions, parameters, and how every number is counted / 版本、参数与计数方式


Raw pipeline outputs (assemblies, RGI tables, geNomad summaries, bins) are held on institutional
storage and are not in this repository because of their size.
原始流程输出(组装、RGI 表、geNomad 汇总、bin)因体积原因保存在机构存储上,不在本仓库内。
