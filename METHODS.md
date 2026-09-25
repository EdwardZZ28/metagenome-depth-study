# Methods / 方法

Everything below is what was actually run, not a plan. Every tool ran inside a pinned
BioContainers image under Apptainer, so the version is fixed by the image tag.
下面写的都是实际运行过的内容,不是计划。所有工具都在 Apptainer 下以固定的 BioContainers
镜像运行,版本由镜像标签锁定。

---

## 1. Input data / 输入数据

Starting material is **quality-filtered (clean) paired reads**, not raw reads. For the public
datasets these were produced with fastp 0.23.2 during earlier work; the in-house samples used the
same version and settings. Depth is therefore always quoted as clean bases, read from the fastp
JSON field `summary.after_filtering.total_bases`.

起点是**质控后的 clean 双端 reads**,不是原始 reads。公开数据集的 clean reads 由 fastp 0.23.2
在此前工作中产生,自有样本使用相同版本和参数。因此文中深度一律指 clean 碱基数,取自 fastp JSON
的 `summary.after_filtering.total_bases` 字段。

Sample selection rules / 样本选择规则:

- Illumina or MGI short-read, paired end / Illumina 或 MGI 短读长双端
- clean depth >= 25 Gbp, so that every rung can be cut from the sample itself
  clean 深度 >= 25 Gbp,保证每一档都能从样本自身切出
- fastp retention >= 80% / fastp 保留率 >= 80%
- at most 2 samples per BioProject / 每个 BioProject 最多 2 个样本

---

## 2. Subsampling / 抽样
seqtk sample -s<SEED> <R1> <FRAC> | gzip -1 > <ID>_1.fq.gz
seqtk sample -s<SEED> <R2> <FRAC> | gzip -1 > <ID>_2.fq.gz
FRAC = min(1.0, target_Gbp / clean_Gbp)


R1 and R2 use the **same seed**, which is what keeps pairs synchronised. Fraction is computed on
bases, not reads, because read length differs between datasets (104 to 150 bp).
R1 和 R2 使用**同一个随机种子**,这是保证配对同步的关键。抽样比例按碱基数计算而非 reads 数,
因为不同数据集的读长不同(104 到 150 bp)。

Every dataset carries a unique id `<sample>__t<tier>__s<seed>`, and every contig is renamed to
`<id>__k141_<n>` immediately after assembly. MEGAHIT contig names repeat across samples, and
mixing them up is a real and previously encountered failure mode.
每个数据集有唯一编号 `<样本>__t<档位>__s<种子>`,组装完成后立刻把每条 contig 改名为
`<编号>__k141_<n>`。MEGAHIT 的 contig 名在不同样本间会重复,混淆是真实发生过的错误。

Verification: at the 2 Gbp rung of R8188_S1, seqtk returned 7,308,794 read pairs totalling
2.002 Gbp, with identical read counts in R1 and R2.
验证:R8188_S1 的 2 Gbp 档,seqtk 返回 7,308,794 个 read pair,合计 2.002 Gbp,R1 与 R2
reads 数完全相等。

---

## 3. Tools and databases / 工具与数据库

| Step | Tool | Image tag | Database |
|---|---|---|---|
| subsampling | seqtk | `1.5--h577a1d6_1` | |
| quality control | fastp | `0.23.2--hadf994f_5` | |
| read-based ARGs | ARGs-OAP | `3.2.4--pyhdfd78af_0` | SARG (bundled) |
| taxonomy | Kraken2 | `2.17.1--pl5321h077b44d_0` | GTDB R220 |
| abundance re-estimation | Bracken | `3.1p1--hc52dbad_0` | GTDB R220, `-r 100 -l S -t 10` |
| assembly | MEGAHIT | `1.2.9--haf24da9_8` | |
| contig ARGs | RGI | `6.0.8--pyh05cac1d_0` | CARD 4.0.1 |
| plasmids, viruses | geNomad | `1.12.0--pyhdfd78af_0` | geNomad DB v1.9 |
| mapping | Bowtie2 | `2.5.5--h63e9258_1` | |
| BAM handling | samtools | `1.24--h9dcdb79_1` | |
| binning | MetaBAT2 | `2.18--h38e344b_2` | |
| bin quality (pending) | CheckM2 | `1.1.0--pyh7e72e81_1` | CheckM2 v3 |
| bin taxonomy (pending) | GTDB-Tk | `2.4.1--pyhdfd78af_1` | GTDB R220 |

GTDB-Tk 2.4.1 is paired with R220 deliberately: it is the officially matched combination and it is
the same database used in the related thesis chapters, so results stay comparable.
GTDB-Tk 2.4.1 与 R220 是有意搭配的:这是官方匹配组合,也和相关论文章节用的是同一个数据库,
结果可以互相比较。

Kraken2 uses the GTDB R220 index rather than an NCBI index so that read-level and bin-level
taxonomy share one naming system.
Kraken2 使用 GTDB R220 索引而非 NCBI 索引,使读长层面和 bin 层面的分类命名统一。

One practical note: ARGs-OAP builds its own BLAST and DIAMOND indices on first run, which fails
inside a read-only container image. The bundled database directory is copied out, indexed once,
and bind-mounted back over its original path.
一个实践细节:ARGs-OAP 首次运行时要自建 BLAST 和 DIAMOND 索引,而容器镜像是只读的,会失败。
解决办法是把自带数据库目录拷出来,建一次索引,再挂载回原路径。

---

## 4. Commands / 命令

```bash
# read-based ARGs
args_oap stage_one -i reads -o argoap -f fq.gz -t $T
args_oap stage_two -i argoap -t $T

# taxonomy
kraken2 --db $KRAKEN_DB --threads $T --paired --gzip-compressed \
        --report <id>.kreport --output /dev/null <id>_1.fq.gz <id>_2.fq.gz
bracken -d $KRAKEN_DB -i <id>.kreport -o <id>.bracken_S.tsv -r 100 -l S -t 10

# assembly
megahit -1 <id>_1.fq.gz -2 <id>_2.fq.gz -o megahit -t $T --min-contig-len 1000

# contig ARGs
rgi main -i <id>.contigs.fa -o rgi/<id> -t contig -a DIAMOND -n $T --local --clean

# mobile genetic elements
genomad end-to-end --cleanup --splits 16 -t $T <id>.contigs.fa genomad $GENOMAD_DB

# mapping and binning
bowtie2-build --threads $T <id>.contigs.fa idx/<id>
bowtie2 -x idx/<id> -1 <id>_1.fq.gz -2 <id>_2.fq.gz -p $T --no-unal \
  | samtools sort -@ 12 -m 2G -o <id>.bam -
jgi_summarize_bam_contig_depths --outputDepth depth.txt <id>.bam
metabat2 -i <id>.contigs.fa -a depth.txt -o bins/<id>.bin -m 1500 -t 16 --seed 1
```

`--min-contig-len 1000` matters: every contig-derived count in this report is over contigs of at
least 1 kb, which is also the threshold used in the related chapters.
`--min-contig-len 1000` 很重要:本报告中所有基于 contig 的计数都建立在 >= 1 kb 的 contig 上,
这也是相关章节使用的阈值。

---

## 5. How each number is counted / 每个数字怎么算

| Column in `L1_summary.tsv` | Definition / 定义 |
|---|---|
| `arg_subtypes` | Number of data rows in `argoap/normalized_16S.subtype.txt`, i.e. distinct ARG subtypes detected. `normalized_16S.subtype.txt` 的数据行数,即检出的不同 ARG subtype 数 |
| `copies_16S` | Sum of the last column of the same file: total ARG abundance in ARG copies per 16S rRNA gene copy. 同一文件最后一列之和,即每个 16S 拷贝对应的 ARG 拷贝总数 |
| `contigs` | Contigs >= 1 kb in the renamed assembly. 改名后组装文件中 >= 1 kb 的 contig 数 |
| `assembly_bp` | Total assembled bases over those contigs. 上述 contig 的碱基总数 |
| `bins` | `.fa` files produced by MetaBAT2, no quality filter applied. MetaBAT2 产出的 `.fa` 文件数,未做质量过滤 |
| `plasmid` | Data rows in geNomad `*_plasmid_summary.tsv`. geNomad 质粒汇总表的数据行数 |
| `conjugative` | Plasmid rows whose `conjugation_genes` field is neither empty nor `NA`. 质粒行中 `conjugation_genes` 字段非空且非 `NA` 的行数 |
| `virus` | Data rows in geNomad `*_virus_summary.tsv`. geNomad 病毒汇总表的数据行数 |
| `rgi_ps` | RGI hits with `Cut_Off` of `Perfect` or `Strict`. `Cut_Off` 为 `Perfect` 或 `Strict` 的 RGI 命中数 |
| `rgi_id80` | Of those, hits with `Best_Identities` >= 80. 其中 `Best_Identities` >= 80 的命中数 |
| `highrisk` | Distinct `Best_Hit_ARO` values among `rgi_id80` hits matching the high-risk pattern below. `rgi_id80` 命中中匹配下列高风险模式的不同 `Best_Hit_ARO` 数 |
| `highrisk_genes` | The gene names behind that count, semicolon separated. 上述计数对应的基因名,分号分隔 |

The **identity >= 80% filter is not optional**. RGI default output at Perfect/Strict is dominated by
low-identity hits to accessory and efflux genes, and in related work that filter changed the
glycopeptide share of a dataset from 44% to 1%. Anything reported here as an ARG hit has passed it.

**identity >= 80% 这个过滤不是可选项**。RGI 默认的 Perfect/Strict 输出被低相似度的辅助基因和
外排泵命中主导,在相关工作中这个过滤把某数据集的糖肽类占比从 44% 变成 1%。本报告中称为
ARG 命中的都通过了这个过滤。

High-risk gene pattern (case-insensitive regular expression) / 高风险基因模式(不区分大小写的正则):

\b(NDM|KPC|OXA-48|OXA-23|IMP-|VIM-|MCR-|mcr-|VanA|vanA|VanB|vanB|CTX-M|TEM-|SHV-)


This is a deliberately broad net that includes carbapenemases, extended-spectrum
beta-lactamases, mobile colistin resistance and vancomycin resistance. It also catches TEM-1 and
SHV-1, which are common and not themselves alarming; they are kept because excluding them would
require a judgement call that belongs in the discussion, not in the counting.
这是有意放宽的范围,包含碳青霉烯酶、超广谱 beta 内酰胺酶、可移动黏菌素耐药和万古霉素耐药。
它也会捕捉到常见且本身并不可怕的 TEM-1 和 SHV-1;保留它们是因为排除需要一个判断,
而那个判断属于讨论部分,不属于计数环节。

---

## 6. Recovery statistic / 回收率统计

For a metric m, sample s and rung t:

recovery(m, s, t) = 100 x value(m, s, t) / value(m, s, 25 Gbp)


Reported values are the **median across samples**, and a sample contributes only if it is complete
at both t and 25 Gbp. This is why the number of samples per row differs in the recovery table.
Taking per-sample ratios first and then the median is deliberate: taking the median first would let
one very complex sample dominate the numerator and denominator separately.

报告值是**样本间的中位数**,且只有在 t 和 25 Gbp 两档都完整的样本才参与计算,这就是回收率
表格中每行样本数不同的原因。先算每个样本的比值再取中位数是有意的:先取中位数会让某个特别复杂的
样本分别主导分子和分母。

---

## 7. Reproducing the tables and figures / 复现表格与图

```bash
python3 scripts/summarise_L1.py   # pipeline output -> data/L1_summary.tsv
python3 scripts/plot_L1.py        # data/L1_summary.tsv -> figures/
```

`summarise_L1.py` needs the pipeline output tree; `plot_L1.py` needs only the committed TSV, so
every figure in this report can be regenerated from what is in this repository.
`summarise_L1.py` 需要流程输出目录;`plot_L1.py` 只需要仓库内的 TSV 文件,所以本报告中的每
一张图都可以仅凭本仓库内容重新生成。

---

## 8. Plasmid fate analysis / 质粒去向分析

The question is not how many plasmid contigs each depth yields, but **which** plasmids are lost.
问题不是每个深度产出多少质粒 contig,而是**哪些**质粒被丢掉了。

Contig names differ between depth rungs because each rung is assembled independently, so
plasmids cannot be matched by name. They are matched by sequence alignment instead.
不同深度档位是独立组装的,contig 名字不同,因此质粒无法按名字匹配,只能按序列比对匹配。

```bash
# reference: plasmid contigs called at 25 Gbp for this sample
# query:     plasmid contigs called at a lower rung of the same sample
minimap2 -x asm10 -t 8 --secondary=no <25Gbp_plasmid.fna> <tier_plasmid.fna> > <sample>__t<tier>.paf
```

Recovery rule / 判定规则:

1. Keep alignments with `nmatch / blocklen >= 0.95` (identity). 保留相似度 >= 0.95 的比对。
2. Merge the kept alignment intervals on each reference contig. 合并参照 contig 上的比对区间。
3. A reference plasmid counts as recovered if merged coverage / reference length >= `MIN_COV`.
   合并后覆盖度除以参照长度 >= `MIN_COV` 即判为回收。

Two values of `MIN_COV` are reported because they answer different questions:
报告两个 `MIN_COV` 值,因为它们回答不同的问题:

| `MIN_COV` | Meaning / 含义 |
|---|---|
| 0.80 (strict) | the plasmid was **rebuilt** at this depth / 该深度下质粒被**完整重建** |
| 0.30 (lenient) | the plasmid was **seen** at this depth, possibly in fragments / 该深度下质粒**被看到**,可能是碎片 |

The gap between the two matters: at low depth a plasmid is often split across several short
contigs, each covering too little of the reference to pass the strict threshold. The strict number
therefore measures assembly completeness, the lenient number measures detection.
两者的差距本身有意义:低深度下一个质粒常被拆成若干条短 contig,每条覆盖参照的比例都不足以
通过严格阈值。因此严格标准衡量组装完整性,宽松标准衡量检出。

Grouping variables come from the 25 Gbp results of the same sample:
分组变量来自同一样本 25 Gbp 的结果:

- **length** from the `length` column of the geNomad plasmid summary / 取自 geNomad 质粒汇总表的 `length` 列
- **conjugative** if `conjugation_genes` is neither empty nor `NA` / `conjugation_genes` 非空且非 `NA`
- **carries an ARG** if the contig name appears in the RGI output with `Cut_Off` of `Perfect` or
  `Strict` and `Best_Identities >= 80`. The RGI `Contig` column holds the plain contig name, so
  the two tables join directly.
  若 contig 名出现在 RGI 输出中且 `Cut_Off` 为 `Perfect` 或 `Strict`、`Best_Identities >= 80`,
  则判为携带 ARG。RGI 的 `Contig` 列就是纯 contig 名,两张表可以直接连接。

Reported rates are the **median of the per-sample recovery rate**, and a sample contributes to a
group only if it has at least 20 reference contigs in that group. Pooling all contigs first would
let the two or three most complex samples dominate every group.
报告值是**各样本回收率的中位数**,且某样本在某分组中至少要有 20 条参照 contig 才参与计算。
若先把所有 contig 合并再计算,两三个最复杂的样本会主导每一个分组。

Run it with / 运行方式:

```bash
python3 scripts/plasmid_fate.py                 # strict, MIN_COV=0.80
MIN_COV=0.3 python3 scripts/plasmid_fate.py     # lenient
```

This script needs the geNomad and RGI outputs plus the PAF files, none of which are in this
repository; the summary it produces is committed as `data/plasmid_fate_summary.txt`.
该脚本需要 geNomad 和 RGI 的输出以及 PAF 文件,这些都不在本仓库内;它产生的汇总结果以
`data/plasmid_fate_summary.txt` 提交。
