#!/bin/sh
snakemake --profile slurm -j 500 --rerun-triggers mtime --rerun-incomplete
