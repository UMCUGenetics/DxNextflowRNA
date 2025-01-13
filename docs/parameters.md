# umcugenetics/dxnextflowrna pipeline parameters
UMCU Genetics RNA seq Workflow

## Input/output options
Define where the pipeline should find input data and save output data.


| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `analysis_id` |  | `string` |  | True |  |
| `input` | Path to input data. | `string` |  | True |  |
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. | `string` |  | True |  |
| `email` | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your
user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.</small></details>| `string` |  |  |  |
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. | `string` |  |  |  |


## Reference genome options


Reference genome related files and options required for the workflow.


| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `fasta` | Path to FASTA genome file. | `string` |  |  |  |
| `fai` | Path to FASTA genome index file. | `string` |  |  |  |
| `gencode_version_name` |  | `string` | GRCh38_gencode_v22_CTAT_lib_Mar012021 |  |  |
| `genome_base` |  | `string` |  |  |  |
| `gene_bed` |  | `string` |  |  |  |
| `gtf` | Path to GTF annotation file. | `string` |  |  |  |
| `ref_flat` |  | `string` |  |  |  |
| `rrna_database_manifest` |  | `string` |  |  |  |
| `rrna_intervals` |  | `string` |  |  |  |
| `sortmerna_index` |  | `string` |  |  |  |
| `sortmerna_index_versions` |  | `string` |  |  |  |
| `star_index` | Path to directory or tar.gz archive for pre-built STAR index. | `string` |  |  |  |


## Institutional config options


Parameters used to describe centralised config profiles. These should not be edited.


| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |  | True |
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't
need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.</small></details>| `string` |
https://raw.githubusercontent.com/nf-core/configs/master |  | True |
| `config_profile_name` | Institutional config name. | `string` |  |  | True |
| `config_profile_description` | Institutional config description. | `string` |  |  | True |
| `config_profile_contact` | Institutional config contact information. | `string` |  |  | True |
| `config_profile_url` | Institutional config URL link. | `string` |  |  | True |


## Generic options


Less common options for the pipeline, typically set in a config file.


| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `cluster_options` |  | `string` | --mail-user null --mail-type FAIL --account=diaggen |  |  |
| `version` | Display version and exit. | `boolean` |  |  | True |
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory.
This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | link |  | True |
| `email_on_fail` | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does
not exit successfully.</small></details>| `string` |  |  | True |
| `plaintext_email` | Send plain-text email instead of HTML. | `boolean` |  |  | True |
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | `string` | 25.MB |  | True |
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |  | True |
| `hook_url` | Incoming hook URL for messaging service <details><summary>Help</summary><small>Incoming hook URL for messaging service. Currently, MS Teams and Slack are supported.</small></details>| `string` |  |  | True |
| `multiqc_config` | Custom config file to supply to MultiQC. | `string` | /home/cog/edejong2/repos/DxNextflowRNA_create_pipeline/assets/multiqc_config.yml |  | True |
| `multiqc_logo` | Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file | `string` |  |  | True |
| `multiqc_methods_description` | Custom MultiQC yaml file containing HTML including a methods description. | `string` |  |  |  |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | True |  | True |
| `pipelines_testdata_base_path` | Base URL or local path to location of pipeline test dataset files | `string` | https://raw.githubusercontent.com/nf-core/test-datasets/ |  | True |
| `rseqc_modules` |  | `string` | bam_stat,infer_experiment,inner_distance,junction_annotation,junction_saturation,read_distribution,read_duplication |  |  |
| `save_non_ribo_reads` |  | `boolean` |  |  |  |
| `seq_platform` | Sequencing platform | `string` | Illumina |  |  |
| `seq_center` | Sequencing center | `string` | UMCU Genetics |  |  |

## Other parameters

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `dx_tracks_path` |  | `string` | /hpc/diaggen/software/development/Dx_tracks/ |  |  |
