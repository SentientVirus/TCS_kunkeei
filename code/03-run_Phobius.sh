#!/bin/bash
mkdir -p results
mkdir -p logs
phobius ~/Akunkeei_files/faa/H3B1-04J_protein.faa > results/Phobius/H3B1-04J_phobius.txt 2> logs/phobius.log
