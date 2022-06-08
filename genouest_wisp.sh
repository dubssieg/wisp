#!/bin/sh
. /local/env/envconda.sh

conda activate /home/genouest/genscale/sdubois/wisp-env

# calling build
# python wisp.py -b -t 8

# calling dl
python utilities.py

conda deactivate