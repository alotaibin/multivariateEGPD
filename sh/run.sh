#!/bin/bash
unset R_HOME

set -e

if [[ -z "${quick_str+x}" ]]; then
  echo "Do you wish to use a very low number of parameter configurations and epochs to quickly establish that the code is working? (y/n)"
  read quick_str
fi

if [[ $quick_str == "y" ||  $quick_str == "Y" ]]; then
    quick=--quick
elif [[ $quick_str == "n" ||  $quick_str == "N" ]]; then
    quick=""
else
    echo "Please re-run and type y or n"
    exit 1
fi

Rscript R/ModelSimulationPlots.R
Rscript R/NeuralEstimation.R $quick
Rscript R/NaiveEstimation.R
Rscript R/ECDF.R
Rscript R/Results.R

# Delete extraneous files to reduce the total memory size
find . -type f -name "network_epoch*" -exec rm {} +  

# Convert pdf images to png
find img/ -type f -iname "*.pdf" -exec sh -c '
  for f do
    out="${f%.pdf}.png"
    gs -q -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r600 -o "$out" "$f" >/dev/null 2>&1
  done
' sh {} +
