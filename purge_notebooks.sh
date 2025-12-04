#!/usr/bin/env bash
# Print cp statements to copy every .ipynb to a same-directory backup: name.bak.ipynb
set -euo pipefail

find . -type f -name '*.ipynb' -print0 |
while IFS= read -r -d '' file; do
    jupyter nbconvert --clear-output --inplace "$file"
done