#!/bin/sh
#
# Output generation script for lagan

# LAGAN=../../../../build/Debug/apps/lagan/lagan
LAGAN=../../../../../build/bin/lagan
# ============================================================
# First Section
# ============================================================


${LAGAN} -i 'input/prototyp_sequences.fa' -o /tmp/fileout_prototype.fa -la 4 2 2 -s 0 -1 -1 > output/editDistance.stdout
${LAGAN} -i 'input/prototyp_sequences.fa' -o /tmp/fileout_prototype.fa -la 4 2 2 -s 1 -1 -1 > output/noEditDistance.stdout