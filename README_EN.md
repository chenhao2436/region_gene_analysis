# Region Gene Analysis

**简体中文：** [README.md](README.md) | **English**

A collection of research scripts for plant BSA and genetic-mapping studies. The repository covers SNP-index / ED calculation and plotting, candidate-gene extraction from BSA signal regions, cross-population shared-peak detection, regional VCF extraction and variant annotation, candidate-gene functional integration, and supporting GVCF/VCF utilities.

> **Status: research toolbox / alpha.** This repository has grown across multiple research stages and is not a single, uniformly configured pipeline. Several legacy Bash drivers still contain absolute paths from the original server. Prefer the Python module entry points with explicit CLI arguments. Annovar, reference genomes, GFF3/GTF annotations, VCFs, and other external resources must be provided by the user.

## Main modules

| Module | Purpose | Main inputs | Main outputs |
|---|---|---|---|
| BSA calculation and plotting | F2 SNP-index / ED calculation and static plots | VCF, parent and bulk samples | SNP-index / ED tables, PDFs |
| `gene_extract` | Expand significant BSA windows and extract genes | Multiple population plot.tsv files, GFF3 | Windows, regions, candidate-gene tables |
| `gene_common` | Summarize candidate genes by population support | Per-sample candidate-gene table | BED/TSV sets supported by at least N or all populations |
| `gene_region_vcf` | Extract regional variants and run Annovar annotation | BED, joint VCF, FASTA, GTF | Regional VCF, avinput, functional annotations |
| `gene_out` | Integrate pepper genes with Arabidopsis RBH/TAIR annotations and variant annotations | Gene intersections, RBH, TAIR, Annovar outputs | Integrated candidate-gene table |
| `gene_extract_SharePeak` | Smooth BSA signals, detect peaks, and cluster peaks across samples | Multiple plot.tsv files, GFF3 | Shared-peak clusters, support intervals, candidate genes |
| `interact_html` | Generate interactive BSA HTML figures | Multiple population plot.tsv files | Interactive HTML |
| GVCF/VCF utilities | Prepare GVCF lists, joint genotyping, and file audits | GVCFs, reference genome | Joint VCF and audit reports |

Each module has its own inputs and parameters. They do not need to be run as one mandatory pipeline.

## Repository layout

```text
region_gene_analysis/
├── F2_BSA.py
├── bsa.fig.py
├── auto_bsa_*.sh
├── audit_bsa_sample_chain.sh
├── gene_extract/
│   ├── scripts/              # gene_extract / gene_common / gene_region_vcf / gene_out
│   └── md/                   # Module notes
├── gene_extract_SharePeak/
│   ├── scripts/              # Shared-peak detection and clustering
│   └── md/
├── interact_html/
│   ├── scripts/              # Interactive BSA figures
│   └── md/
├── ploidy30_bsa_bundle_20260318_upload_v2/
│   └── gvcf_plidy_30_BSA.py
└── gvcf2vcf_genome_gpt_3m6d.sh
```

The repository also retains research notes, backups, historical scripts, and archives. Use the module entry points documented here rather than treating every historical copy as a separate product.

## Software environment

The Conda environment file covers common Python dependencies used by the main analysis modules:

```bash
conda env create -f environment.yml
conda activate region-gene-analysis
```

Core Python dependencies:

- Python 3.10+
- NumPy
- SciPy (peak detection in Shared Peak)
- Matplotlib (static plots)
- Plotly (interactive plots)

Some modules also require external executables or databases:

| Tool/resource | Purpose |
|---|---|
| bcftools | VCF region extraction, sample and genotype queries |
| tabix / bgzip | VCF compression, indexing, and region access |
| GATK 4 | GVCF joint-genotyping helper scripts |
| Annovar | Variant functional annotation; install separately and follow its license |
| UCSC `gtfToGenePred` | GTF conversion |
| Reference FASTA | Coordinate and sequence reference |
| GFF3 / GTF | Gene coordinates and annotations |
| TAIR annotation and RBH file | Map pepper candidate genes to Arabidopsis functional annotations |

These external tools are not all installed by the Conda environment file. Confirm their availability and versions before running the corresponding module.

## Recommended analysis path: from BSA signals to candidate genes

### 1. Prepare BSA plot.tsv files

Generate one `plot.tsv` per population using an upstream BSA workflow. For a valid comparison, use the same reference assembly, consistent chromosome naming, and compatible window settings.

### 2. Extract genes around significant BSA windows

The main entry point `gene_extract/scripts/gene_extract` accepts:

- `--plot-tsv`: one or more population plot.tsv files.
- `--gff3`: genome GFF3 annotation.
- `--output-dir`: output directory.
- `--delta-threshold`: DELTA_SNP_INDEX threshold (0.4 is used in the historical example).
- `--flank`: flanking expansion in bp (1,000,000 bp in the example).

Example:

```bash
python gene_extract/scripts/gene_extract \
  --plot-tsv sample_36.plot.tsv sample_70.plot.tsv sample_168.plot.tsv \
  --gff3 reference.gff3 \
  --output-dir ./results/gene_extract \
  --delta-threshold 0.4 \
  --flank 1000000
```

Outputs include significant windows, sample regions, genes per window, sample candidate genes, and BED/TSV files. See `gene_extract/README.md` and `gene_extract/md/` for fields and interpretation.

### 3. Summarize cross-population support

```bash
python gene_extract/scripts/gene_common \
  --sample-candidate-tsv ./results/gene_extract/gene_extract.sample_candidate_genes.tsv \
  --output-dir ./results/gene_common \
  --min-support 2 \
  --all-support 5
```

`--min-support` defines the minimum number of populations required for inclusion; `--all-support` defines the full-support count. Set these values according to the actual number of input populations rather than assuming every analysis has five.

### 4. Alternative strategy: Shared Peak

`gene_extract_SharePeak/scripts/gene_share_peak` smooths the BSA signal, detects positive and negative peaks, then clusters peaks by chromosome, direction, and distance.

```bash
python gene_extract_SharePeak/scripts/gene_share_peak \
  --plot-tsv sample_36.plot.tsv sample_70.plot.tsv sample_168.plot.tsv \
  --gff3 reference.gff3 \
  --output-dir ./results/shared_peak \
  --smooth-window 10 \
  --min-width 5 \
  --prominence 0.05 \
  --distance-window 3000000 \
  --flank 1000000 \
  --min-support 2
```

This is a different signal-selection strategy from the fixed DELTA-threshold approach in `gene_extract`. Do not combine their outputs without documenting the methods.

### 5. Interactive BSA figures

```bash
python interact_html/scripts/interact_html \
  --plot-tsv sample_36.plot.tsv sample_70.plot.tsv \
  --output-dir ./results/interactive
```

Use `--sample-labels` to set display labels. `--cluster-tsv` and `--min-cluster-support` can overlay Shared Peak clusters in script versions that support those options. Outputs generally include combined and subplot HTML files. Run `--help` to confirm the arguments for the exact script version.

## Regional VCF extraction and Annovar annotation

`gene_extract/scripts/gene_region_vcf` extracts joint-VCF records overlapping candidate-gene BED intervals with bcftools and invokes the external Annovar toolchain. Required inputs include:

- Candidate-gene BED and matching TSV.
- bgzip-compressed and indexed joint VCF.
- Reference FASTA matching the VCF.
- GTF matching the genome assembly.
- Annovar installation and executable scripts.
- `bcftools`, `tabix`, and `gtfToGenePred`.

Example:

```bash
python gene_extract/scripts/gene_region_vcf \
  --gene-bed ./results/gene_common/gene_common.at_least_2.bed \
  --gene-table ./results/gene_common/gene_common.at_least_2.tsv \
  --vcf joint.vcf.gz \
  --reference-fasta reference.fa \
  --gtf reference.gtf \
  --annovar-dir /path/to/annovar \
  --bcftools-bin bcftools \
  --tabix-bin tabix \
  --gtfToGenePred-bin gtfToGenePred \
  --label at_least_2 \
  --output-dir ./results/variant_annotation
```

The current source retains default GTF, Annovar, and bcftools paths from the original server. Always provide `--gtf`, `--annovar-dir`, and `--bcftools-bin` explicitly instead of relying on those defaults.

## Important limitations

- This is a multi-module research toolbox, not one universal end-to-end workflow. Run and validate modules independently.
- Some root-level and subdirectory Bash drivers still contain absolute paths such as `/data2/chenh/`; they may not work on other machines.
- `gene_region_vcf` requires Annovar and other external tools. Annovar is not distributed here.
- TAIR annotations, pepper-to-Arabidopsis RBH mappings, reference genomes, and genome annotations must be supplied and version-matched by the user.
- BSA thresholds, window sizes, peak width, prominence, clustering distance, and population-support count all influence results. These are analysis parameters, not universal biological conclusions.
- Historical scripts, backups, research notes, and older configurations remain in the repository and may differ from the current entry points.
- This delivery does not claim end-to-end execution of every module on real data. The documented capabilities and commands are based on static source inspection and recorded usage.

## Pre-run checklist

1. Reference FASTA, VCF, GFF3, and GTF must match the same genome assembly.
2. Sample names must match the VCF header exactly.
3. Confirm VCF indexes (`.tbi` / `.csi`) before regional access.
4. Verify plot.tsv column order, especially DELTA_SNP_INDEX and SNP-count columns.
5. Test with a small region or a few samples before genome-wide or large-cohort runs.
6. Record commands, parameters, software versions, and logs.

## Roadmap

- [ ] Replace hard-coded Bash paths with reusable CLI parameters or configuration files.
- [ ] Provide shared example configurations and input templates for modules.
- [ ] Add synthetic fixtures and automated tests for core Python modules.
- [ ] Clean generated caches and Python bytecode from version control.
- [ ] Record software versions, defaults, and methodological boundaries for each module.

## Support

When opening an issue, specify the module, command, software versions, input-column definitions, and relevant logs. Do not commit large sequencing datasets, unauthorized data, or private server configuration to the public repository.
