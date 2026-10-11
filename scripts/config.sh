# ===== depth study: paths and versions (single source of truth) =====
export DS=/scratch/user/uqnzhai/depth_study          # working area (scratch)
export SIF=$DS/containers                             # apptainer images
export DB=$DS/db                                      # local copies of small DBs
export WORK=$DS/work                                  # per-batch temporary data
export RDM_OUT=/QRISdata/Q6636/depth_study_results    # final results (RDM)
export PLAN=$HOME/depth_study/sample_plan.tsv

# databases
export CARD_DIR=$DB                                   # contains localDB/ (CARD 4.0.1); rgi must run with cwd=$CARD_DIR or a symlinked localDB
export GENOMAD_DB=$DB/genomad_db                      # geNomad DB v1.9
export KRAKEN_DB=/scratch/opendata/genomics/Kraken2/Refseq/GTDB/v220
export BRACKEN_LEN=100
export CHECKM2_DIR=/scratch/opendata/genomics/CheckM2/version_3/CheckM2_database
export GTDBTK_DATA_PATH=/scratch/opendata/genomics/GTDB/releases/release220

# tool images (BioContainers tags recorded in pull_sif.sbatch)
export SIF_SEQTK=$SIF/seqtk_1.5.sif
export SIF_FASTP=$SIF/fastp_0.23.2.sif
export SIF_ARGSOAP=$SIF/args_oap_3.2.4_b0.sif
# ARGs-OAP cannot build its own DB inside the read-only image; use the pre-indexed copy
export ARGOAP_DB=$DB/argoap_db
export ARGOAP_DB_IN=/usr/local/lib/python3.12/site-packages/args_oap/db
export SIF_KRAKEN2=$SIF/kraken2_2.17.1.sif
export SIF_BRACKEN=$SIF/bracken_3.1p1.sif
export SIF_MEGAHIT=$SIF/megahit_1.2.9.sif
export SIF_RGI=$SIF/rgi_6.0.8.sif
export SIF_GENOMAD=$SIF/genomad_1.12.0.sif
export SIF_BOWTIE2=$SIF/bowtie2_2.5.5.sif
export SIF_SAMTOOLS=$SIF/samtools_1.24.sif
export SIF_METABAT2=$SIF/metabat2_2.18.sif
export SIF_CHECKM2=$SIF/checkm2_1.1.0.sif
export SIF_GTDBTK=$SIF/gtdbtk_2.4.1.sif

# run a tool inside its container with the storage areas visible
run() { local sif=$1; shift; apptainer exec -B /scratch,/QRISdata "$sif" "$@"; }
run_argoap() { apptainer exec -B /scratch,/QRISdata -B $ARGOAP_DB:$ARGOAP_DB_IN $SIF_ARGSOAP "$@"; }
