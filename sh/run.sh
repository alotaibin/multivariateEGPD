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
Rscript R/ModelSimulationPlots.R
Rscript R/NeuralEstimation.R $quick
Rscript R/NaiveEstimation.R
Rscript R/ECDF.R
Rscript R/Results.R

###############################################################################
# Application study pipeline
###############################################################################
Rscript R/data/Process_ecad_rainfall_data.R
Rscript R/ECA_data_diagnostic_plots.R
Rscript R/tale_comparison_plots.R

###############################################################################
# Clean-up and conversion
###############################################################################

# Delete extraneous network checkpoint files
find . -type f -name "network_epoch*" -exec rm {} +

# Convert PDFs to PNGs for figures
find img/ -type f -iname "*.pdf" -exec sh -c '
  for f do
    out="${f%.pdf}.png"
    gs -q -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r600 -o "$out" "$f" >/dev/null 2>&1
  done
' sh {} +

# Convert PDFs in Application/Figures as well
find Figures/Application/ -type f -iname "*.pdf" -exec sh -c '
  for f do
    out="${f%.pdf}.png"
    gs -q -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r600 -o "$out" "$f" >/dev/null 2>&1
  done
' sh {} +
