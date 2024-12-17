# 简介

![HapPhase](static/HapPhase_logo.png)

HapPhase是一款用于扩增子长读长测序数据的单体型分型和拷贝数计算的软件。它主要针对像SMN基因这样的复杂基因，通过利用SNP位点的连锁信息和等位基因频率比例，准确地进行单体型分型计算，从而推断基因的拷贝数。

## 功能特点

- **连锁分析**：基于SNP的连锁关系，对扩增子长读长数据进行单体型连锁识别。
- **拷贝数计算**：利用单体型的比例信息和点突变的ratio信息，推断目标基因的拷贝数。
- **高效解析**：支持VCF和BAM文件的快速解析，适用于大规模数据分析。

## 安装

```shell
# 克隆项目仓库
git clone https://github.com/[your-repo]/HapPhase.git

# 进入项目目录
cd HapPhase

# 构建项目
cargo build --release
```

## 使用方法

```shell
./happhase --input <input.bam> --vcf <variants.vcf> [其他选项]
```

### 参数说明

- `--input` / `-i`：输入的BAM格式的比对文件。
- `--vcf` / `-v`：包含SNP信息的VCF文件。
- `--region` / `-r`：指定分析的基因区域（可选）。
- `--snp_num` / `-s`：用于分型的最小SNP数量，默认值为4。
- `--min_ratio` / `-m`：SNP的最小等位基因频率，默认值为0.2。
- `--max_ratio` / `-M`：SNP的最大等位基因频率，默认值为0.8。
- `--min_proportion`：单体型的最小比例，默认值为0.1。
- `--min_count`：单体型的最小支持读数数目，默认值为2。
- `--method`：相似性计算方法，默认值为`cosine`

### 示例

```shell
./happhase --input sample.bam --vcf variants.vcf --region "chr5:68900000-69000000" --snp_num 4 --min_ratio 0.1 --max_ratio 0.7 --method cosine
```

## 结果解读

运行完成后，程序将输出识别到的单体型数量以及每个单体型的详细信息，包括对应的SNP位点和等位基因序列。

## 许可证

MIT License. 详见 LICENSE。

## 联系方式

如有任何问题或建议，请联系：[jiancghen2@genomics.cn](mailto:jiancghen2@genomics.cn) 或 [cherryamme@qq.com](mailto:cherryamme@qq.com)。
```
