#!/bin/bash
#SBATCH --comment=773320000
#SBATCH --time=2-25:00:00
#SBATCH --mem=40000
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm.output_filter_merge%j.txt
#SBATCH --error=slurm.error_filter_merge%j.txt
#SBATCH --job-name=filter_merge
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matt.watts@wur.nl

cd /lustre/nobackup/WUR/ABGC/shared/Taurus/ARS-UCD1.2

ml legacy
ml python 

INPUT_FASTA="GCF_002263795.1_ARS-UCD1.2_genomic.fna"
OUTPUT_FASTA="renamed_ARS-UCD1.2.fna"
MAP_FILE="chr_name_map.txt"

python - <<EOF
name_map = dict(line.strip().split() for line in open("$MAP_FILE"))
with open("$INPUT_FASTA") as infile, open("$OUTPUT_FASTA", "w") as out:
    for line in infile:
        if line.startswith(">"):
            old = line[1:].split()[0]
            new = name_map.get(old, old)
            out.write(">" + new + "\\n")
        else:
            out.write(line)
EOF
