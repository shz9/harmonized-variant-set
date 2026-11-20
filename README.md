# Harmonized variant sets for ancestry-specific LD references

This standalone workflow builds a versioned candidate panel from **HapMap3+ ∪
genome-wide-significant curated GWAS Catalog variants**, then selects variants
within each ancestry using **gnomAD v4.1 genomes MAF > 0.5% AND MAC > 20**.
GRCh38 is the canonical build; validated GRCh37 exports retain the same stable
variant identities. It produces variant manifests for use by downstream LD
pipelines.

This directory is the complete project root and can be copied or moved into a
new repository under any name. It contains all code, tests, configuration and
Python requirements; no files or Python modules from a parent repository are
needed. Missing reference datasets are downloaded automatically to the configured
paths, or existing shared copies can be reused.

```text
.
├── README.md
├── .gitignore
├── config.json             # paths, ancestry mappings and selection thresholds
├── requirements.txt        # Python dependencies for processing and tests
├── run.sh                  # public setup/check/download/submit entry point
├── workflow.py             # offline processing and submission implementation
├── job.sbatch              # internal SLURM worker
└── tests/
    └── test_workflow.py
```

## Selection procedure

1. **Freeze discovery sources.** Start with the published LDpred2 HapMap3+ map,
   which contains 1,444,196 SNPs. Source rows without GRCh38 coordinates,
   valid rsIDs or valid SNP alleles are excluded and recorded with reasons in
   `candidates/hapmap3_unresolved.tsv.gz`; counts appear in the candidate-stage
   manifest. Add **curated
   GWAS Catalog association records with p < 5 × 10⁻⁸**, across all traits and
   study ancestries. This does not search full summary-statistics datasets.
   Preserve the complete significant association records, including study
   accessions and trait annotations, for provenance. See the
   [HapMap3+ documentation](https://privefl.github.io/bigsnpr/articles/LDpred2.html)
   and [Catalog download formats/releases](https://www.ebi.ac.uk/gwas/docs/file-downloads/).
2. **Resolve biological alleles on GRCh38.** HapMap3+ supplies an allele pair;
   orient it against the reference FASTA, complementing where needed. For
   Catalog records, require a unique biallelic explicit DNA allele pair from
   local GRCh38 dbSNP at the Catalog position and current rsID. SNPs, short
   indels and other explicit sequence alleles are supported; symbolic alleles,
   ambiguous multiallelic dbSNP records, and unresolved mappings are audited.
   Association rows representing interactions or multiple variants are not
   expanded into independent hits. No positional guesses or online lookups run
   on compute nodes.
3. **Normalize and deduplicate.** Use `bcftools norm -f ... -c e`, followed by
   sorting. REF must match the FASTA. Left-align and trim indels. Use
   `chr:position:REF:ALT` on GRCh38 as the stable `VARIANT_KEY`; rsIDs are aliases,
   not primary keys. Deduplicate by the normalized allele key while retaining
   both source labels and all resolved rsID aliases. Coordinates are **1-based
   VCF positions**, never BED coordinates. See
   [bcftools normalization](https://samtools.github.io/bcftools/bcftools.html#norm).
4. **Annotate ancestry counts from gnomAD v4.1 genomes.** Require a PASS,
   biallelic record with the exact normalized allele key. Use ancestry-specific
   adjusted `AC_<group>` and `AN_<group>`, rather than global AF, raw counts or
   filtering allele frequency (FAF). Use genomes for consistent genome-wide
   ascertainment; do not substitute exome or joint exome/genome frequencies.
   v4.1 includes an allele-number correction relative to v4.0; see the
   [release notes](https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1).
5. **Select independently within each ancestry.** For a biallelic variant:

   ```text
   MAC = min(AC, AN - AC)
   MAF = MAC / AN
   keep iff PASS AND AN > 0 AND MAF > 0.005 AND MAC > 20
   ```

   AC counts the ALT allele; it is not necessarily the minor allele. The
   inequalities are strict: MAC must be at least 21; exactly 0.5% fails. Missing
   counts, no called alleles, absent frequency records, duplicate records and
   failed reference QC have distinct audit outcomes. An absent gnomAD variant
   is not automatically assigned frequency zero: distinguishing nonvariation
   from lack of calls requires the additional
   [all-called-sites AN resource](https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1),
   which this initial implementation does not ingest.
6. **Build and verify a GRCh37 mapping.** Lift normalized GRCh38 variants using
   Picard `LiftoverVcf`, normalize against the GRCh37 FASTA, then lift back and
   normalize again against GRCh38. Require exact recovery of the original
   chromosome, position, REF and ALT; one-to-one forward and reverse mappings;
   no target allele-key collisions; and agreement with HapMap3+'s native
   GRCh37 position where available. Reject cross-chromosome mappings. Picard
   handles reverse complements and can recover SNP REF/ALT swaps. See its
   [allele-aware liftover documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360036449552-LiftoverVcf-Picard).
   A failed mapping removes a variant from the GRCh37 export, not the GRCh38
   master or ancestry export. REF/ALT swaps for complex indels are conservatively
   rejected if Picard cannot recover them; they are not forced through.

The unfiltered master retains all successfully canonicalized candidates, even
if they fail every ancestry's frequency criteria. This preserves rare GWAS hits
and allows a future panel version to change thresholds without losing provenance.

## Ancestry definitions and interpretation

Mappings are explicit in [config.json](config.json), using the
[gnomAD group definitions](https://gnomad.broadinstitute.org/news/2023-11-genetic-ancestry):

| Output ancestry | Default gnomAD genomes group | Interpretation |
| --- | --- | --- |
| AFR | `afr` | African/African American proxy |
| AMR | `amr` | Admixed American proxy |
| CSA | `sas` | South Asian proxy; does not cover all Central Asian populations |
| EAS | `eas` | East Asian proxy |
| EUR | `nfe` | Non-Finnish European proxy |
| MID | `mid` | Middle Eastern proxy |

These broad reference groups do not reproduce every biobank's ancestry
composition. `EUR: ["nfe", "fin"]` is supported if a broader European definition
is desired: sum the disjoint groups' AC and AN first, then compute MAF/MAC.
Do not average group frequencies. Group-specific counts must exist in the VCF
header; a missing group is an error, never a fallback to another ancestry.

The MAC threshold is evaluated **in gnomAD**, not in the target biobank. Small
reference groups will have a higher effective frequency cutoff because MAC > 20
can dominate MAF > 0.5%. Target biobanks must still report and check their own
ancestry-specific counts, missingness, imputation quality, sample quality and
relatedness before computing LD. A variant common in gnomAD may be rare or
unavailable in a particular cohort. This reference selection does not establish
that LD can be estimated accurately from a small ancestry sample.

## Running on the cluster

Run from this project's root (the directory containing `run.sh`):

```bash
# One command: environment, downloads, reference preparation and analysis.
bash run.sh --account def-sgravel --chromosomes 22

# Production: choose a fresh output_dir/version after the pilot.
bash run.sh --chromosomes 1-22 --max-concurrent 2 --download-workers 4 --scratch-dir /scratch/szabad
```

No manual resource download or index creation is required. The launcher starts
a detached **login-node controller** and prints its log/status paths. It survives
terminal disconnection; rerunning the same command reuses the controller or
restarts a stopped one. On a new run, omitting `--chromosomes` selects 1–22; on a
retry it reuses the saved chromosome list.

**Compute nodes do not need internet access. All network operations run on the
login node.** The controller downloads at most two chromosome VCFs at a time
(default), with four parallel HTTP connections per VCF. It submits a chromosome
SLURM task only after that VCF and its index are available on shared scratch.
Scientific extraction, allele normalization, ancestry filtering and liftover run
in strictly offline SLURM jobs. All jobs request 1 CPU, 8 GiB and one hour by
default; `--download-workers` controls login-node transfer concurrency only.

The controller polls SLURM, reuses completed stages, and advances the download
window as workers finish and remove raw files. `--max-concurrent` bounds both
outstanding chromosome jobs and the managed scratch staging window across all
array submissions. It never queues unstaged chromosomes or the final merge
before every required chromosome is complete or has a valid worker dependency.

Raw gnomAD files go under `/scratch/szabad/harmonized-variant-set/<output-id>/chr_N/`,
where the output ID isolates projects/panel versions. `--scratch-dir` or the
`scratch_dir` config key changes the root. At two jobs, budget up to roughly
94 GB for raw gnomAD VCFs, plus the shared dbSNP reference and other working files.
The default dbSNP path is also on `/scratch/szabad/`; it is retained for reuse.
Existing configured local VCFs are read in place and never deleted.

Downloads use parallel HTTP ranges against the public GCS endpoint, without
requiring Google credentials or `gcloud`. They stream into a single preallocated
scratch file rather than keeping a second set of chunk files. A completed-slice
checkpoint enables resume, conditional ETag requests keep slices on one object
version, and the completed GCS file is verified against its supplied MD5.
Progress is logged in the chromosome job's SLURM output. Common-reference
transfers continue to use resumable curl downloads.

Each job extracts and normalizes only candidate positions, retains exact
candidate alleles and writes `output/<version>/gnomad_counts/chr_N/counts.tsv.gz`.
Columns are `VARIANT_KEY,CHROM,POS,REF,ALT,FILTER` plus `AC_<group>,AN_<group>` for
all source ancestry groups needed by the configuration. Per-ancestry pooled
counts, MAF/MAC and selection outcomes remain in the existing ancestry audit
and export files. Duplicate reference records, missing counts, failed filters
and absent variants retain their separate QC outcomes.

The count table and manifest are published atomically and checked before the
job deletes its managed raw gnomAD VCF and index. Failed downloads/extractions
retain scratch data for retry. A matching count cache skips both downloading
and extraction on a worker retry; input/ancestry-group signatures and a table
SHA256 prevent stale or corrupt cache reuse. Rerunning the same launcher command
automatically resumes missing/failed stages and reuses completed stages and
active jobs in the same output directory.
Download helpers explicitly refuse to run when `SLURM_JOB_ID` is set. An offline
worker with a missing local input fails immediately instead of trying HTTP. On
upgrade, the controller cancels this run's pending chromosome/merge jobs created
by the incompatible network-on-compute worker and replaces them after staging.
Already running jobs are allowed to finish before their slots are reused.

All stages default to `01:00:00`. Override with `--time 00:30:00` or the
`wall_time` configuration key. On rerun, pending pipeline jobs with longer limits
are shortened in place; running jobs are left alone. Timeouts retain partial
downloads and count caches for the next run.

The default account is read from `job.sbatch` (currently `def-sgravel`). Supply
`--account def-sgravel` explicitly, or set `slurm_account` in `config.json`, to
override it. `SBATCH_ACCOUNT` is also supported when the configuration does not
specify an account. Every submitted stage receives the selected `--account`.
An account rejection leaves the run recoverable: rerun the same command after
correcting the account. No output-directory change or provenance deletion is
needed.

The login controller checks SLURM state itself with `squeue` and `sacct`. It reuses
pending/running jobs, submits only missing or failed chromosomes, and replaces
pending downstream jobs whose dependencies failed or were replaced. A durable
submission journal and run-specific SLURM comments recover jobs accepted before
the launcher could save their IDs. Submission and worker-stage locks prevent
duplicate launchers or requeued workers from repeating completed work. If the
scheduler cannot be queried reliably, submission stops before duplicating jobs;
rerun when the scheduler is available.

If a newly submitted worker fails, the controller records the error and a bounded
worker-log excerpt and stops rather than repeatedly resubmitting the same failure.
Rerunning the launcher authorizes another attempt using retained downloads and
caches. Account rejections and transfer failures are also saved in controller
status. Check progress or the last failure without submitting anything:

```bash
bash run.sh --status
```

Status is stored in `output/<version>/provenance/controller.json`; the printed
`logs/hvs-controller-*.out` path contains transfer progress and submission history.
A lock on the shared output directory ensures only one login controller runs,
even when the launcher is invoked from multiple login nodes. Rerunning while it
is active updates operational settings such as concurrency and wall time.

Completed runs are a no-op on rerun. Incomplete stage directories are retained
under `output/<version>/recovery/` before being rebuilt; count caches and partial
scratch downloads are preserved. Account, wall time, scratch and download-worker settings
can change on retry. Changes to scientific selection settings or an explicitly
different chromosome list require a separate output version to avoid mixing
panels.

Optional commands:

```bash
bash run.sh --setup-only                         # software only; no data or jobs
bash run.sh --download-sources                   # HapMap3+ and Catalog only
bash run.sh --download-references --chromosomes 22 # common references only; no jobs
bash run.sh --check-only --chromosomes 22          # offline check; requires prepared indexes
```

`--download-references` also stages discovery sources. It downloads the common
FASTA, chain and dbSNP resources and the matching dbSNP index. It does not fetch
gnomAD. FASTA unpacking, `.fai`/`.dict` generation, and indexing of existing local
VCFs happen in the preparation job on a normal run. After download-only, invoke the normal run to complete preparation.

Paths in the configuration, including `--config`, are relative to this directory
unless absolute. This also holds when invoking `run.sh` by its path from another
working directory. The default dbSNP location is
`data/reference/GCF_000001405.40.gz` within this project; set an absolute path
to reuse an independently staged shared copy. The project-local `.gitignore`
excludes reference data, environments, intermediate outputs, logs and Python
caches. Final merged variant sets and their manifest remain available to Git.

A normal `bash run.sh` fetches all missing discovery and reference resources
before submission. Shared reference paths in `config.json` are reused when
present. The optional `reference_urls` object can override individual download
URLs using the same keys as the resource paths; `gnomad_vcf` URLs support
`{chrom}`. Local paths remain the processing inputs.

HapMap3+ downloads use Figshare's API endpoint, with an alternate download host
as fallback. The website download URL can return an empty bot challenge, which
is detected and rejected. Downloads retry up to three times, validate gzip/ZIP
integrity, and replace the destination only after validation succeeds. Reruns
restart incomplete `.partial` files and repair empty destination files; completed
sources are reused. No manual removal of partial files is needed.

Environment creation and all downloads run on login; scientific extraction runs
on offline compute nodes. `--setup-only` installs the
Python dependencies and checks that bcftools, the configured Picard executable
and SLURM sbatch/curl are available, without requiring reference data or submitting
jobs. It can be rerun to reuse the environment and install missing Python
dependencies. `--check-only`
does not install, download or submit. The download command reuses existing files;
to update sources, stage new dated files and select a new panel version. The
Catalog snapshot is identified automatically from the staged file when
`gwas_catalog_release` is `"auto"` (the default). Existing `REPLACE_...`
placeholders, empty values and null also select automatic identification, so
previously downloaded files work without editing that field or downloading again.
When the filename or ZIP member includes the official `rYYYY-MM-DD` and `eNNN`
tokens, the workflow records the release date and Ensembl version. Otherwise it
records an exact `local-sha256:...` snapshot identity and marks the upstream
release as unknown. It never infers a release from download time or today's
`latest` release. This metadata is frozen into submission provenance and stage
manifests; it does not rewrite your working configuration. An explicitly supplied
release string is preserved.

Identification runs offline during preflight. It reads ZIP directory metadata
or, for generic filenames, hashes only the staged Catalog file; it does not scan
the large reference VCFs/FASTAs. Download receipts also record source SHA256
checksums, and empty or incomplete new downloads are rejected. Prefer a
[dated Catalog release](https://ftp.ebi.ac.uk/pub/databases/gwas/releases/)
when assembling a reproducible release.

The entry point loads `python/3.13`, `bcftools` and `picard` through modules when
available, and creates its own `env/harmonized` environment. Override
`PYTHON_MODULE` and `HVS_TOOL_MODULES` consistently for your cluster.
The default `picard_command` is `["java", "-jar", "$EBROOTPICARD/picard.jar"]`,
matching the EasyBuild Picard module. Both the launcher and SLURM jobs load the
module before reading the configuration. The workflow expands `EBROOTPICARD`,
checks that the JAR exists, and records the resolved command in provenance.
No separate Picard download is needed.

For a different installation, set `picard_command` to `["picard"]` for a wrapper
on PATH, or `["java", "-jar", "/shared/software/picard.jar"]` for a fixed JAR path.
Set `HVS_TOOL_MODULES=""` when the tools and Java are already available on PATH.
Python 3.11+, its venv support, curl, SLURM and the configured tool modules must be
provided by the cluster; setup does not install system software. Picard JARs
configured with `java -jar` must also be staged at the configured path.

Setup updates pip before installing dependencies, prefers prebuilt wheels, and
allows `pyreadr` to build from source when no compatible wheel is available.
A source build needs a C compiler, Python headers and zlib, bzip2 and lzma
development headers supplied by the cluster (load the appropriate compiler and
library modules). NumPy, pandas and pysam still require wheels.
If `pyreadr` reports "No matching distribution found", also check whether your
configured package index contains it; setup uses your existing pip index settings.
After updating the repository, rerun `bash run.sh --setup-only`; an existing
`env/harmonized` can be reused.

The submission creates a **prepare → candidates → chromosome array → merge**
chain using `afterok` dependencies. The prepare job is omitted if all common
references and indexes already pass preflight. It requests 1 CPU, 8 GiB and
one hour; chromosome jobs use the resources described above, while candidates
and merge jobs request 1 CPU, 8 GiB and one hour. Java is capped at 4 GiB.
Preparation unpacks FASTAs, creates FASTA indexes and Picard dictionaries, and
builds absent indexes for existing local VCFs. Preparation writes the final
preflight inventory before downstream jobs can start. `SLURM_TMPDIR` holds
small working subsets; full gnomAD files are staged by the login controller on
the configured shared scratch. No worker stage performs network requests.

Automatically staged resources (or provide existing paths in `config.json`):

| Resource | Requirements |
| --- | --- |
| HapMap3+ | Published `map_hm3_plus.rds`, or an equivalent TSV with `chr,pos,pos_hg38,rsid,a0,a1` |
| GWAS Catalog | Full v1.0.2 curated associations ZIP containing a single TSV, or unpacked TSV; snapshot identified automatically; verify GRCh38 build |
| gnomAD | Official **v4.1 genomes** per-chromosome sites VCFs plus `.tbi`/`.csi` indexes; no sample genotypes required |
| dbSNP | GRCh38 sites VCF plus `.tbi`/`.csi`, staged locally or at an explicitly configured shared path |
| GRCh38 / GRCh37 | Matching FASTAs with `chr1`…`chr22`, `.fai` and Picard `.dict` files |
| Chains | Matching UCSC `hg38ToHg19.over.chain.gz` and `hg19ToHg38.over.chain.gz` |

Use the [official gnomAD downloads](https://gnomad.broadinstitute.org/downloads)
for current public mirrors and checksums. The release resource path is
`gs://gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz`,
as defined by the
[official resource code](https://github.com/broadinstitute/gnomad_qc/blob/main/gnomad_qc/v4/resources/release.py).
These are large files: check available shared storage or reuse an institutional
copy. The downloader uses the public HTTPS mirror of this pinned v4.1 genomes
path. dbSNP downloads are pinned to [b156](https://ftp.ncbi.nih.gov/snp/archive/b156/VCF/).
FASTA and chain downloads come from UCSC hg38/hg19. Each downloaded reference
gets a `.download.json` receipt containing its URL, size and timestamp. No
full-file hashes or scans of these large references run in login-node preflight.

Reference assemblies, contig naming and chains must agree. The defaults expect
UCSC-style hg38/hg19 primary autosomes; do not silently mix `1` and `chr1` or use
GRCh37/hg19 references with differing sequence conventions. FAI lengths provide
a preliminary build check; allele checks and round trips provide variant-level
validation. Record upstream checksums for the FASTAs, chains and VCFs alongside
the staged resources. Preflight records their paths, sizes and modification times;
these are not cryptographic verification. Use only trusted, pre-normalized
gnomAD/dbSNP sources: indexed position queries precede normalization, so an
arbitrarily padded, nonnormalized indel could otherwise be missed upstream.

Configuration, workflow code and the worker script are frozen into immutable
revision files under `output_dir/provenance` before submission. Queued/running
jobs keep their original revision when the launcher is updated. `jobs.json`
records current job IDs, `submissions.json` records submission attempts and
recovered jobs, and `run.json` records the requested chromosomes and account.
Source hashes, the resolved Catalog snapshot, configuration, code hashes,
exclusion counts and ancestry selection counts are recorded in stage manifests.
Outputs appear atomically after each stage completes. Retry the same launcher
command to recover; there is no need to inspect job IDs, cancel failed dependency
chains manually, delete provenance, or discard cached data.

## Outputs and build conventions

```text
output/v0.1.0/
  provenance/                 # submitted config/code, input metadata, job IDs
  candidates/                 # per-chromosome candidates, association evidence, rejection audits
    hapmap3_unresolved.tsv.gz  # source rows lacking usable coordinates/IDs/alleles
  chr_22/
    master.GRCh38.tsv.gz       # all canonical variants, source labels, rsID aliases, mapping status
    EUR.audit.tsv.gz          # every master variant plus counts and selection outcome
    EUR.GRCh38.tsv.gz         # selected ancestry variants
    EUR.GRCh37.tsv.gz         # selected variants with validated GRCh37 mapping
    unresolved_candidates.tsv.gz
    rejected37.vcf            # Picard rejects, with FILTER reasons
    rejected_roundtrip.vcf
    manifest.json
  gnomad_counts/
    chr_22/
      counts.tsv.gz           # candidate AC/AN records for configured source ancestry groups
      manifest.json           # source metadata, candidate signature, table checksum
  merge/                      # concatenated master and ancestry/build manifests
```

The default output directory is inside this Git repository. After a successful
run, `output/v0.1.0/merge/` contains the final compressed master and ancestry/build
variant sets (`*.tsv.gz`) and `manifest.json`. These files are not ignored, so
they can be committed and uploaded to GitHub:

```bash
git add output/v0.1.0/merge/
git commit -m "Add v0.1.0 variant sets"
git push
```

The ignore rules also allow final exports under `output/<version>/merge/` for
future versions. Keep `output_dir` in that layout to use these rules.

Every ancestry is exported. `VARIANT_KEY` always identifies the normalized
**GRCh38** allele, providing the cross-build join key. `CHROM,POS,REF,ALT` describe
the file's named build. GRCh38 exports also contain `CHROM37,POS37,REF37,ALT37`
where valid. GRCh37 exports contain `CHROM38,POS38,REF38,ALT38` and a
`TARGET_VARIANT_KEY` for the GRCh37 allele. Their frequency columns are explicitly
named `AC_GRCH38_ALT`, `AN_GRCH38`, `MAF_GRCH38`, `MAC_GRCH38`: liftover may swap
REF/ALT, so the original AC must not be interpreted as a target-build ALT count.
MAF/MAC remain allele-symmetric for these biallelic mappings.

## Using the panel across biobanks

Freeze the **same version and ancestry manifest** for every biobank. Match each
biobank's variants by normalized, reference-checked build-specific allele keys,
then propagate the stable GRCh38 key into LD metadata. Resolve reversed allele
coding explicitly; signed Pearson correlations change sign when dosage coding
is flipped for one variant. Never force the same rsID onto different alleles or
match indels using position alone. Palindromic SNPs are safe when REF/ALT and the
reference build are known; allele-only or strand-uncertain matches need additional
evidence and otherwise remain unresolved.

For every ancestry and chromosome, report: requested panel size, exact allele
matches, variants with unresolved representation, imputation-quality failures,
missingness failures, frequency/count failures, block-coverage failures, and
variants actually exported into LD. Report pairwise shared variants relative to
each bank's eligible set and to their union. Investigate unexplained losses
before release; agree on acceptance thresholds after the chromosome-22 pilot.
Do not conceal failures by redefining the panel as an intersection of available
biobank variants. A shared reference selection improves comparability; source
coverage and cohort QC still determine achievable overlap.

Use broadly ascertained WGS genotypes when building LD references. For example,
an HGDP+1KG ancestry-classifier training panel applied a pooled HWE selection
upstream; see the
[source preparation notebook](https://zenodo.org/records/14286454/files/GIANT_Extract_HGDP_1KG_Data.ipynb?download=1).
Harmonizing identifiers cannot restore variants absent from that source. This
workflow does not reconstruct missing genotypes.

HapMap3+ plus significant Catalog hits is a useful reproducible baseline, but
does not guarantee 1–3 million variants **after ancestry filtering**, nor complete
fine-mapping coverage. Many Catalog hits overlap HapMap3+; many others are rare
in a particular ancestry. Fine mapping generally needs additional well-imputed
variants within each locus, including variants without prior associations. For
a later release, evaluate ancestry-specific tagging/coverage and add a fixed,
documented backbone of common WGS variants or locus-specific dense panels rather
than padding to a target count. Publish a shared ancestry-specific core plus
optional dense extensions, retaining stable identities across versions.

## Validation

```bash
python3 -m venv env/harmonized
env/harmonized/bin/python -m pip install --upgrade pip
env/harmonized/bin/python -m pip install --prefer-binary --only-binary=numpy,pandas,pysam -r requirements.txt
env/harmonized/bin/python -m unittest discover -s tests -v
bash -n run.sh job.sbatch
```

Use Python 3.11 or newer (the cluster default is `module load python/3.13`). If
`run.sh` has already created the environment, only the test and syntax-check
commands are needed. Tests use small synthetic references and do not require
SLURM, downloaded datasets, external bcftools/Picard executables or another
repository. The reference integration tests use system curl and a local HTTP
server, with synthetic data and simulated Picard/SLURM commands. The same `requirements.txt` supplies all their dependencies.

Tests cover strict MAF/MAC boundaries, count pooling, significant association
parsing, native-build checks, reciprocal allele mapping, collision rejection,
and build-specific exports. A tiny offline integration fixture exercises real
FASTA/VCF IO and bcftools through pysam, with Picard isolated by a mock and a
command-contract test. A real chromosome-22 run with the staged reference files
is still required to validate the external Picard/chain/reference combination,
measure retention and confirm cluster runtime before production.


## Functional expansion and AlphaGenome 1 (optional v0.2 workflow)

The original `config.json` still reproduces the HapMap3+/GWAS workflow.
Use a separate configuration/output version for functional discovery:

```bash
# Functional annotations are downloaded automatically on login.
bash run.sh --config config.functional.json --chromosomes 22
# Use a fresh output_dir before switching a chromosome pilot to production.
bash run.sh --config config.functional.json --chromosomes 1-22
```

The example selects up to **5 million master variants**, with a desired minimum
of **3 million**. These are master counts, not guaranteed counts per ancestry.
It preserves the complete canonicalized HapMap3+/Catalog backbone, including its
rare hits. New variants must pass the original ancestry-specific gnomAD genomes
`MAF > 0.005 AND MAC > 20` rule in at least one configured ancestry.
If fewer than 3 million eligible variants are available, the report explicitly
says `below_minimum`; unrelated variants are not added to fill a quota.

### Functional regions and supported alleles

`config.functional.json` freezes GENCODE v50 basic and ENCODE Registry V4.
The worker discovers variants directly from gnomAD, independently of rsIDs or
GWAS significance. Regions are:

- CDS, including annotated start/stop codons;
- 20 intronic bases adjoining each exon–intron junction;
- UTRs and lncRNA exons;
- promoters from 2 kb upstream through 500 bp downstream of each basic
  protein-coding/lncRNA transcript's TSS, respecting strand (both endpoints included);
- all eight ENCODE cCRE classes, preserving their individual element identifiers.

GTF coordinates are converted from 1-based inclusive to BED's 0-based half-open
convention. BED sources use exact boundaries. Overlaps retain every category,
gene and element label. Splice flanks are selection windows, not predicted
splice consequences; CDS membership does not distinguish missense from synonymous.
Optional conservation or BaselineLD-derived intervals can be added under `beds`
with category `CONSERVED` or another supported category. BaselineLD's intercept,
MAF and LD-related annotations must not be exported as functional intervals.

New alleles are biallelic SNPs and normalized simple insertions/deletions of
1–49 bases excluding the shared VCF anchor. Other replacements, symbolic and
multiallelic records are excluded. Common PASS alleles are normalized **before**
intersecting functional intervals; deletion spans and insertion breakpoints
are considered, avoiding anchor-only overlap errors. Indels are retained without
an AlphaGenome score. Existing build-validation rules still apply: GRCh37 failure
does not remove a GRCh38 variant.

### AlphaGenome 1 / Atlas AVI

Use `config.functional-avi.json` after placing the downloaded, unzipped
`alphagenome_variant_impact_score_snvs.tsv.gz` **and its `.tbi` index** at the
configured shared path. Obtain the original file through the
[official AlphaGenome downloads page](https://deepmind.google.com/science/alphagenome/downloads).
That page currently requires sign-in; this repository does not bundle the data,
credentials or a guessed direct download URL. A versioned direct URL may be set
in the resource configuration for automatic download of the data file; its index
must also be staged. No model inference/API calls occur on workers.

The native tabix input has these columns:

```text
#CHROM  POS  REF  ALT  raw_score  PHRED
```

Workers use indexed `bcftools annotate` matching chromosome, position, REF and
ALT against the common-variant candidates. They do not load or enumerate all
possible substitutions in Python. The example uses **AVI PHRED >=30** as a
configurable high-score inclusion rule, including variants outside the functional
intervals. AVI combines AlphaGenome 1 regulatory predictions with other evidence;
it is distinct from an individual AlphaGenome track score and from AlphaGenome 2.
The threshold is a prioritization choice, not a pathogenicity classification or
a percentile calibrated within common variants. See the
[AVI VEP plugin documentation](https://github.com/Ensembl/VEP_plugins/blob/release/116/AVI.pm).

Smaller score tables are also supported. For example, replace the `alphagenome1`
resource list with:

```json
[
  {
    "path": "data/functional/ag1_scores.tsv.gz",
    "release": "your verified AlphaGenome 1 score snapshot",
    "build": "GRCh38",
    "model": "AlphaGenome1",
    "format": "table",
    "kind": "quantile",
    "score_column": "quantile_score",
    "scorer": "your recorded scorer and calibration",
    "threshold": 0.99
  }
]
```

Tables need explicit `CHROM POS REF ALT` columns, a `VARIANT_KEY`, or native tidy
`variant_id` values such as `chr22:250:A>C`. Alleles are normalized/reference-checked
before matching. `kind` can be `avi_phred`, `quantile`, or `raw`; thresholds and
release/model provenance are required. Quantile/raw selection uses absolute
magnitude, so the table must contain comparable scorer/track outputs. Use
`filters` (column names mapped to lists of accepted strings) to restrict an
export to particular scorers, tracks or tissues, or provide separate resources
with separate thresholds. Raw values are never compared across score sources.
Optional `columns` maps input column names to canonical names, and `separator`
can be changed from tab to comma. Missing/nonfinite scores fail explicitly.

### Comparison with standard PRS panels

To prioritize high-score variants missing from specified panels, supply their
**normalized GRCh38 allele manifests**, for example:

```json
"prs_panels": [
  {
    "name": "SBayesRC_7M",
    "path": "data/functional/sbayesrc7m.grch38.tsv.gz",
    "release": "your verified SBayesRC panel snapshot",
    "build": "GRCh38",
    "normalized": true
  }
]
```

The manifest needs `VARIANT_KEY` or `CHROM POS REF ALT`. Native SBayesRC
`snp.info` A1/A2 columns are not REF/ALT, and GRCh37 positions must not be relabeled
as GRCh38; orient, normalize and validate any converted panel first. Membership
uses exact alleles, not rsID aliases. No SBayesRC manifest is bundled.
`PRS_STATUS=absent_from_supplied_panels` means absent from **all supplied lists**,
not proven to lack LD proxies. With no comparison lists it is `not_assessed`;
the workflow never silently claims novelty against unspecified standard panels.

### Size selection, audit and operational behavior

The merge applies the global budget after collecting the chromosome candidate
pools. Priority is: preserved backbone; qualifying AlphaGenome variants absent
from the supplied PRS panels; other qualifying AlphaGenome variants; CDS; splice
flanks; UTR; promoters; conserved regions; enhancers; lncRNA exons; other regulatory
regions. Within an oversubscribed tier, each candidate is assigned deterministically
to one eligible ancestry using an allele-key hash. Equal budgets are allocated
across nonempty ancestry strata and unused capacity is redistributed. Remaining
ties use the stable hash. This avoids a global-European-frequency ranking, but
does not guarantee equal ancestry counts or optimal LD tagging.

The final merge contains:

- the existing master and ancestry/build manifests, restricted to the same selected identities;
- `panel_selection.tsv.gz`: every candidate's selection status, priority, ancestry
  stratum and reason, including budget exclusions;
- `functional_expansion.json`: pool/backbone/final counts, size status, category
  and allele-type counts, AlphaGenome selection counts and PRS comparison coverage.

Category counts overlap and must not be summed. Full candidate pools and ancestry
audits remain under the chromosome directories, and direct-discovery counts are
cached under `functional_counts/chr_N/`. Sources/releases, annotation hashes and
large indexed-resource identities are recorded. Changed scientific settings or
inputs require a fresh output version. Both discovery and backbone count caches
must be available before a retry can skip staging the raw gnomAD VCF.
Managed gnomAD files are deleted only after both caches have been published;
existing candidate-only caches cannot supply the functional expansion.

This implementation expands functional coverage; it does not measure LD coverage
or add a genetic-map-based genome-wide scaffold. The local regression suite uses
real FASTA/VCF/tabix IO and bcftools through pysam; production still requires the
normal cluster tools and staged references.
