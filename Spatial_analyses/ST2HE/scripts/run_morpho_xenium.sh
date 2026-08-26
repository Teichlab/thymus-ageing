#!/bin/bash
#BSUB -J "morphoFM[8]"            # Array job: one task per .ndpi file
#BSUB -G s4152                      # Group name
#BSUB -q gpu-interactive                 # Queue name (batch, not interactive)
#BSUB -n 4                            # Number of CPU cores (per task)
#BSUB -M 32000                        # Memory in MB (per task)
#BSUB -R "span[hosts=1] select[mem>32000] rusage[mem=32000]"
#BSUB -gpu "num=1:gmem=10GB"                   # One GPU per task
#BSUB -o /nfs/team205/lm25/thymus_projects/thymus_ageing_atlas/Spatial_analyses/scripts/morpho/_logs/morphoFM_%J_%I.out   # %J=job ID, %I=array index
#BSUB -e /nfs/team205/lm25/thymus_projects/thymus_ageing_atlas/Spatial_analyses/scripts/morpho/_logs/morphoFM_%J_%I.err
#BSUB -W 01:30                      # Wallclock limit (per task)

# 1. Load necessary modules or environments
source ~/.bashrc
conda activate morpho-env

# Each array task processes a single slide end-to-end, selected by its 1-based
# array index: preprocess -> tiles, then generate UNI2-h embeddings for that same
# slide. The '&&' skips embedding if preprocessing fails for this task.
SCRIPT_DIR=/nfs/team205/lm25/thymus_projects/thymus_ageing_atlas/Spatial_analyses/scripts/morpho

python "${SCRIPT_DIR}/preprocess_morpho.py" --index "${LSB_JOBINDEX}" --file-list "${SCRIPT_DIR}/xenium_slides.txt" && \
python "${SCRIPT_DIR}/generate_embeddings_morpho.py" --index "${LSB_JOBINDEX}" --file-list "${SCRIPT_DIR}/xenium_slides.txt" --stain-norm --grayscale --batch-size 384