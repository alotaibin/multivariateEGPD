#!/bin/bash
unset R_HOME
set -e

# Ask user about quick mode
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

###############################################################################
# Simulation + Estimation pipeline
###############################################################################
Rscript R/simulations/ModelSimulationPlots.R
Rscript R/simulations/NeuralEstimation.R $quick
Rscript R/simulations/NaiveEstimation.R
Rscript R/simulations/ECDF.R
Rscript R/simulations/Results.R

###############################################################################
# Application study pipeline
###############################################################################
Rscript R/data/Process_ecad_rainfall_data.R
Rscript R/application/ECA_data_diagnostic_plots.R
Rscript R/application/tale_compirson_plots.R

###############################################################################
# Clean-up and conversion
###############################################################################

# Delete extraneous network checkpoint files
find . -type f -name "network_epoch*" -exec rm {} +

# Convert PDFs to PNGs for simulation figures
find intermediates/Figures/Simulation/ -type f -iname "*.pdf" -exec sh -c '
  for f do
    out="${f%.pdf}.png"
    gs -q -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r600 -o "$out" "$f" >/dev/null 2>&1
  done
' sh {} +

# Convert PDFs in application figures
find intermediates/Figures/Application/ -type f -iname "*.pdf" -exec sh -c '
  for f do
    out="${f%.pdf}.png"
    gs -q -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r600 -o "$out" "$f" >/dev/null 2>&1
  done
' sh {} +
