#!/bin/bash
echo "Note: This script is designed to run with the amount of memory detected by BBMap."
echo "      If Samtools crashes, please ensure you are running on the same platform as BBMap,"
echo "      or reduce Samtools' memory setting (the -m flag)."
samtools sort -m 49G -@ 3 trimmed_Plate1-9H-9_S72.mapped.bam -o trimmed_Plate1-9H-9_S72.mapped_sorted.bam
samtools index trimmed_Plate1-9H-9_S72.mapped_sorted.bam
