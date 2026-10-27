#!/bin/bash
#BSUB -J "morphoFM_nonXenium[1-38]"     # Array job: one task per line of nonXenium_slides.txt
#BSUB -q gpu-normal                 # GPU queue (site-specific; adjust)
#BSUB -n 4                          # CPU cores per task
#BSUB -M 80000                      # Memory in MB per task (whole-slide images are read into RAM)
#BSUB -R "span[hosts=1] select[mem>80000] rusage[mem=80000]"
#BSUB -gpu "num=1:gmem=10GB"        # One GPU per task
#BSUB -o _logs/morphoFM_%J_%I.out   # %J=job ID, %I=array index (relative to the submission dir)
#BSUB -e _logs/morphoFM_%J_%I.err
#BSUB -W 01:30                      # Wallclock limit per task

# Submit from this directory:  cd ST2HE/scripts && mkdir -p _logs && bsub < run_morpho_nonXenium.sh
# DATA_DIR = where <sample>_processed/ output folders are written (and where relative
# entries in the slide list are resolved). Override with: DATA_DIR=/my/dir bsub < ...
DATA_DIR=${DATA_DIR:-/lustre/scratch126/cellgen/teichmann/lm25/data/thymus_he/thyAgeing_HnE}

source ~/.bashrc
conda activate morpho-env

# Each array task processes a single slide end-to-end, selected by its 1-based
# array index: preprocess -> tiles, then generate UNI2-h embeddings for that same
# slide. The '&&' skips embedding if preprocessing fails for this task.
SCRIPT_DIR=${LSB_SUBCWD:-$(pwd)}

python "${SCRIPT_DIR}/preprocess_morpho.py" --index "${LSB_JOBINDEX}" --data-dir "${DATA_DIR}" --file-list "${SCRIPT_DIR}/nonXenium_slides.txt" && \
python "${SCRIPT_DIR}/generate_embeddings_morpho.py" --index "${LSB_JOBINDEX}" --data-dir "${DATA_DIR}" --file-list "${SCRIPT_DIR}/nonXenium_slides.txt" --grayscale --batch-size 384
