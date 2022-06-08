#!/bin/sh
. /local/env/envconda.sh

conda activate wisp-env

# calling build
python wisp.py -b

conda deactivate