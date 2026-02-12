# Analyzing ppmSeq Output CRAMs - Guide for Preprocessing and SingleReadSNV

## Overview

This document addresses how to determine whether ppmSeq CRAMs have been preprocessed (trimmed, aligned, sorted), how to run the `ppmSeq_preprocess` workflow if needed, and how to connect the output to the `SingleReadSNV` workflow for somatic variant detection.

---

## 1. Determining If CRAMs Are Already Preprocessed

The key question is whether the CRAM files have already gone through ppmSeq-specific adapter trimming. You can check this by inspecting the CRAM tags present in the reads:

### Already trimmed (ppmSeq v1 adapters)
If the CRAM contains the **`st`** and **`et`** tags (Start_loop and End_loop pattern match results), the data has already been trimmed with the ppmSeq trimmer. These tags indicate ppmSeq adapter trimming was performed.

### Already trimmed (legacy v5 adapters)
If the CRAM contains **`as`** and **`ts`** tags instead, it was trimmed with the older legacy v5 ppmSeq adapter format.

### Native adapter trimming only (partially trimmed)
If the CRAM contains the **`a3`** tag (native adapter start position) but **not** `st`/`et` tags, then only UG native adapters were trimmed, but the ppmSeq adapters and loop structures were not. This happens when `application_type` was set to `'native'` instead of `'ppmSeq'` on the sequencer.

### Not trimmed at all
If none of the above tags are present, the CRAM has not been trimmed.

**To inspect tags in a CRAM file**, you can use `samtools view` to check a few reads:
```bash
samtools view <your_file>.cram | head -5 | tr '\t' '\n' | grep -E '^(st|et|as|ts|a3):'
```

---

## 2. Running the ppmSeq_preprocess Workflow

### Choosing the Right Template

Based on the state of your input CRAMs, select the appropriate input template from
`workflows/ppmSeq_preprocess/input_templates/`:

| Input Data State | Template | adapter_version |
|---|---|---|
| Untrimmed ppmSeq v1 data | `ppmSeq_preprocess_template-ppmSeq.json` | `"v1"` |
| Untrimmed legacy v5 data | `ppmSeq_preprocess_template-ppmSeq_legacy_v5.json` | `"legacy_v5"` |
| Native adapters trimmed, ppmSeq adapters NOT trimmed | `ppmSeq_preprocess_template-ppmSeq_post_native_adapter_trimming.json` | `"v1"` |

### Correct Parameters for Each Template

#### Template 1: Standard ppmSeq (v1) - Most Common

```json
{
    "ppmSeqPreprocess.adapter_version": "v1",
    "ppmSeqPreprocess.trimmer_parameters": {
        "formats_description": "s3://ultimagen-workflow-resources-us-east-1/trimmer-formats/2.30/public/ppmSeq/ppmSeq.json",
        "format": "ppmSeq",
        "pattern_files": [
            "s3://ultimagen-workflow-resources-us-east-1/trimmer-formats/2.30/public/ppmSeq/ppmSeq_start_loop.csv",
            "s3://ultimagen-workflow-resources-us-east-1/trimmer-formats/2.30/public/ppmSeq/ppmSeq_end_loop.csv"
        ],
        "extra_args": "--histogram %1.%2.%4.%5.histogram.csv",
        "failure_read_group": "unmatched",
        "cram_reference": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta"
    },
    "ppmSeqPreprocess.ua_parameters": {
        "ua_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/UA/b38-v45-79372c0.uai",
        "ref_alt": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt",
        "ua_extra_args": "--seed-score-ratio 0.5 --vector --soft-clipping",
        "v_aware_alignment_flag": false
    },
    "ppmSeqPreprocess.sorter_params": {
        "mark_duplicates": true,
        "aligned": true,
        "mark_duplicates_ends_read_uncertainty": 0,
        "coverage_intervals": "s3://ultimagen-workflow-resources-us-east-1/interval_lists/coverage_intervals.hg38.tar.gz"
    },
    "ppmSeqPreprocess.steps": {
        "trim": true,
        "align": true,
        "sort": true
    }
}
```

#### Template 2: Legacy v5

The main differences from v1:
- `adapter_version`: `"legacy_v5"` instead of `"v1"`
- `trimmer_parameters.formats_description`: uses `ppmSeq_legacy_v5.json`
- `trimmer_parameters.format`: `"ppmSeq_legacy_v5"` instead of `"ppmSeq"`
- No `pattern_files` needed
- Histogram extra_args uses `%1.%2.%3.%4.%5` (additional index `%3`)

UA and sorter parameters remain identical.

#### Template 3: Post Native Adapter Trimming

The main differences from v1:
- `trimmer_parameters.formats_description`: uses `ppmSeq_post_native_adapter_trimming.json`
- `trimmer_parameters.format`: `"ppmSeq_post_native_adapter_trimming"`
- Still uses `adapter_version: "v1"` and the same pattern_files

UA and sorter parameters remain identical.

### Common Required Inputs (All Templates)

```json
{
    "ppmSeqPreprocess.input_cram_bam_list": ["<your_file_1>.cram", "<your_file_2>.cram"],
    "ppmSeqPreprocess.base_file_name": "<sample_name>",
    "ppmSeqPreprocess.references": {
        "ref_fasta": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta",
        "ref_fasta_index": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
        "ref_dict": "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.dict"
    },
    "ppmSeqPreprocess.ref_fastas_cram": [
        "s3://ultimagen-workflow-resources-us-east-1/hg38/v0/Homo_sapiens_assembly38.fasta",
        "s3://ultimagen-workflow-resources-us-east-1/hg38/methyl_seq_ref/hg38_Lambda_pUC19.fa",
        "s3://ultimagen-workflow-resources-us-east-1/hg19/v0/Homo_sapiens_assembly19.fasta",
        "s3://ultimagen-workflow-resources-us-east-1/hg38/rna-seq/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta"
    ],
    "ppmSeqPreprocess.cpu": 32,
    "ppmSeqPreprocess.preemptible_tries": 1,
    "ppmSeqPreprocess.no_address": true,
    "ppmSeqPreprocess.create_md5_checksum_outputs": false
}
```

### Valid adapter_version Values

The workflow validates `adapter_version` against this enum:
- `"v1"` - Current standard ppmSeq adapter
- `"legacy_v5"` - Older ppmSeq adapter (generally unavailable since 2024)
- `"legacy_v5_start"` - Legacy v5 variant (start-specific)
- `"legacy_v5_end"` - Legacy v5 variant (end-specific)
- `"dmbl"` - DMBL adapter variant

### Preprocessing Outputs

The workflow produces:
- **`output_cram_bam`** - Trimmed, aligned, sorted CRAM file
- **`output_cram_bam_index`** - CRAM index (.crai)
- **`sorter_stats_json`** - Sorter statistics in JSON (needed by SingleReadSNV)
- **`sorter_stats_csv`** - Sorter statistics in CSV
- **`report_html`** - ppmSeq QC report (HTML) - **this is the QC report file**
- **`application_qc_h5`** - QC metrics in HDF5 format
- **`aggregated_metrics_json`** - QC metrics in JSON format
- **`trimmer_histogram`** / **`trimmer_stats`** - Trimmer statistics
- **`bedgraph_mapq0`** / **`bedgraph_mapq1`** - Coverage bedgraphs

---

## 3. Connecting to SingleReadSNV Workflow

Once preprocessing completes, the SingleReadSNV workflow requires the following outputs from the preprocessing step:

| SingleReadSNV Input | Source from ppmSeq_preprocess |
|---|---|
| `input_cram_bam_list` | `output_cram_bam` |
| `input_cram_bam_index_list` | `output_cram_bam_index` |
| `sorter_json_stats_file_list` | `sorter_stats_json` |

Use the **`single_read_snv_template-ppmSeq.json`** template (for v1 adapter data) which includes:
- ppmSeq-specific CRAM tags to copy (`st`, `et`, `tm`, `a3`, `rq`, `MI`, `DS`, `sd`, `ed`, `l1`-`l7`, `q2`-`q6`)
- ppmSeq-specific features list for the XGBoost model
- Appropriate training parameters (`tp_train_set_size: 3200000`, `fp_train_set_size: 3200000`, `num_CV_folds: 3`)

The SingleReadSNV workflow will **not work correctly** on untrimmed CRAMs because:
1. It expects ppmSeq-specific tags (`st`, `et`, etc.) in the CRAM reads, which are added by the trimmer
2. The sorter_json_stats_file is required for computing downsampling rates and coverage statistics

---

## 4. Troubleshooting Common Issues

### "Unable to get the workflow to successfully run"

Common failure causes with the ppmSeq_preprocess workflow:

1. **Wrong adapter_version**: Ensure the adapter version matches the actual adapter chemistry used during sequencing. If unsure, try `"v1"` first (the current standard).

2. **Missing pattern_files**: The v1 template requires two pattern files (`ppmSeq_start_loop.csv` and `ppmSeq_end_loop.csv`). Legacy v5 does not require pattern files.

3. **Wrong formats_description**: Each adapter version uses a different trimmer formats JSON file. Ensure these match:
   - v1: `ppmSeq.json`
   - legacy_v5: `ppmSeq_legacy_v5.json`
   - post native: `ppmSeq_post_native_adapter_trimming.json`

4. **Input already trimmed**: If the CRAM was already trimmed (check for `st`/`et` tags), running the trimmer again will fail or produce incorrect results. In that case, preprocessing is not needed - proceed directly to SingleReadSNV.

5. **S3 resource access**: All resource files reference `s3://ultimagen-workflow-resources-us-east-1/`. Ensure your HealthOmics run role has read access to this bucket.

6. **Cloud provider**: If running on AWS HealthOmics, you may need to add:
   ```json
   "ppmSeqPreprocess.cloud_provider_override": "aws"
   ```

### SingleReadSNV failures on unpreprocessed data

The SingleReadSNV workflow expects:
- A **trimmed, aligned, sorted** CRAM with ppmSeq tags
- A **sorter_stats.json** file from the sorter step
- Sufficient coverage (default minimum: 8x) to train the model; if below this, a pre-trained model must be supplied

If running on CRAMs that came directly from the sequencer without preprocessing, the workflow will fail because the required tags and statistics files are missing.

---

## 5. Summary of Recommended Workflow

```
Raw CRAM from sequencer
    │
    ▼
[Check tags: st/et present?]
    │
    ├── YES → CRAMs are already preprocessed
    │         → Transfer QC report to bucket
    │         → Proceed to SingleReadSNV with ppmSeq template
    │
    └── NO → Run ppmSeq_preprocess workflow
              │
              ├── Check for a3 tag (native adapter trimmed)?
              │   ├── YES → Use post_native_adapter_trimming template
              │   └── NO  → Use standard ppmSeq (v1) template
              │
              ▼
         Preprocessed CRAM + sorter_stats.json + QC report
              │
              ▼
         Run SingleReadSNV with single_read_snv_template-ppmSeq.json
              │
              ▼
         Somatic variant calls (FeatureMap VCF with SNVQ scores)
```
