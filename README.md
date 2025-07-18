# BARCODE: BActeRial COre DEfensome

## Description

BARCODE is a pipeline for analyzing bacterial genomes to identify core and quasi-core defense systems. It processes genomes through the following steps:
- **Quality Control** using [CheckM](https://pubmed.ncbi.nlm.nih.gov/25977477/)
- **Genome Dereplication** using [dRep](https://pubmed.ncbi.nlm.nih.gov/28742071/)
- **Annotation** using [Prokka](https://pubmed.ncbi.nlm.nih.gov/24642063/)
- **Defense-system Detection** using [DefenseFinder](https://pubmed.ncbi.nlm.nih.gov/35538097/)
- **Core and Quasi-core Analysis** to identify defense systems present across genomes

## 1. Environment Preparation

To run this pipeline use the pre-built **Docker image** that includes both required environments:
- **core_quasicore_pipeline**: For quality control, dereplication, annotation, and defense system detection.
- **core_quasicore_analysis**: For core and quasi-core analysis.

First, ensure that Docker is installed and running on your system. Then, pull the Docker image:

```bash
docker pull angibvg/barcode:latest
```

Next, launch the container and mount your working directory:

```bash
docker run --platform your_platform -it -v /path/to/your/folder:/data -w /data angibvg/barcode
```

- `--platform`: Select your platform (e.g., linux/amd64).
- `-v`: Mounts the host directory /path/to/your/folder into the container at /data, so all files you put in that host folder are accessible (and writable) inside the container.
- `-w`: Sets the working directory inside the container to /data, so the container’s shell or entry‑point will start in that folder.

## 2. Pipeline Usage

### 2.1. Input Directory Structure

Before running the pipeline, prepare a parent folder (`INPUT_DIR`) containing one subdirectory per species. Each species folder must contain its raw genome FASTA files (`*.fasta` or `*.fna`).

**Example Directory Structure:**

```bash
INPUT_DIR/
├── Species__A/
│   ├── GCA_XXXXXX1.fna
│   └── GCA_XXXXXX2.fna
├── Species__B/
│   ├── GCA_YYYYYY1.fna
│   └── GCA_YYYYYY2.fna
└── ...
```

No additional structure is required. The pipeline will create a `results/` directory inside `INPUT_DIR` to store all outputs.

### 2.2. Launch the Pipeline: `BARCODE.sh`

To run the entire pipeline, which includes quality control, dereplication, annotation, defense system detection, and core/quasi-core analysis, use the following command:

```bash
./BARCODE.sh -i INPUT_DIR/ -t 12 -g 10
```

**Options:**
- `-i INPUT_DIR`: Parent folder holding one subfolder per species.
- `-t THREADS`: (Optional) Number of threads to use (default is 10).
- `-g MIN_GENOMES`: (Optional) Minimum number of genomes to retain per species after quality control and dereplication (default is 2).



**Pipeline Steps:**

1. **Quality Control (CheckM)**
   - Runs `CheckM lineage_wf` on each species’ genome set.
   - Retains genomes with completeness ≥ 90% and contamination ≤ 5%.
   - Rejects low-quality genomes and logs them under `results/rejected/low_quality/`.
   - If fewer than `MIN_GENOMES` genomes remain, the species is moved to `results/rejected/dRep_too_few/too_few_genomes/`.

2. **Dereplication (dRep)**
   - Runs `dRep dereplicate` on QC-passed genomes.
   - Keeps representative genomes; moves duplicates to `results/rejected/duplicates/`.
   - If fewer than `MIN_GENOMES` genomes remain, the species is moved to `results/rejected/dRep_too_few/`.

3. **Annotation (Prokka)**
   - Annotates each retained genome with Prokka.
   - Outputs results under `results/prokka/SPECIES/`.

4. **Defense-system Detection (DefenseFinder)**
   - Runs DefenseFinder on each Prokka-derived `.faa` file (protein sequences).
   - Outputs results under `results/defense/SPECIES/`.

5. **Core and Quasi-core Analysis**
   - Generates presence/absence tables and extracts `core.txt` and `qscore.txt` per species.

## 3. Output Directory & Files

After execution, a `results/` directory will be created inside `INPUT_DIR` with the following structure:

```bash
INPUT_DIR/results/
├── checkm/
│   └── Species_A/…  (CheckM outputs)
├── dRep/
│   └── Species_A/…  (dereplication outputs)
├── prokka/
│   └── Species_A/…  (Prokka outputs)
├── defense/
│   └── Species_A/
│       └── …  (DefenseFinder outputs)
└── core_qscore/
    └── Species_A/
        ├── Species_A_defense_systems_percentages_presence.tsv
        ├── Species_A_core.txt
        └── Species_A_qscore.txt
```

Note: The Species_A_core.txt and Species_A_qscore.txt files will only be created if at least one core or quasi‑core system is detected; otherwise, only the Species_A_defense_systems_percentages_presence.tsv file will be present.


**Key Output Files in `core_qscore/`:**
- `Species_A_defense_systems_percentages_presence.tsv`: A table showing the percentage presence of each defense system across the genomes of Species A.
- `Species_A_core.txt`: A list of defense systems present in all genomes of Species A (core defense systems).
- `Species_A_qscore.txt`: A list of defense systems present in a high percentage of genomes, but not necessarily all (quasi-core defense systems).
  
