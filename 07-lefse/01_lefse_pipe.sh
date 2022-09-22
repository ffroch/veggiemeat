#!/bin/bash
conda activate lefse

#data transformation
lefse_format_input.py lefse_data.txt lefse_data.in -f 'r' -c 2 -s -1 -o 1000000 -u 1

#analysis
lefse_run.py lefse_data.in lefse_data.res -l 4

#plot
lefse_plot_res.py lefse_data.res lefse_data.png

#cladogram
lefse_plot_cladogram.py lefse_data.res lefse_data.cladogram.png --format png

conda deactivate