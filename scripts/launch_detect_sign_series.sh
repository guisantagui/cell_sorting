#!/bin/bash

series_file="../test_data/sample_series.csv"
nboot=5000
out_dir="../results/series_stats/"

Rscript detect_sign_series_args.R $series_file --do_plots --smooth_loess --nboot $nboot --outDir $out_dir