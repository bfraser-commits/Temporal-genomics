#!/bin/bash
#!/bin/bash
#SBATCH --time=00:45:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -A Research_Project-T109423
#SBATCH -p pq
#SBATCH --job-name=candidate_region
#SBATCH --error=candidate_region.err.txt
#SBATCH --output=candidate_region.out.txt
#SBATCH --export=All
#SBATCH -D .

module load SAMtools/1.17-GCC-12.2.0
SEL=candidate_window_sel2

awk 'NR > 1 {
  gsub(/"/, "", $1);                # remove quotes from first field
  split($1, parts, "_");            # split by underscore
  chrom = parts[1];
  start = parts[2];
  end = start + 49999;
  print chrom ":" start "-" end;
}' $SEL.txt > $SEL.regions

REF=/gpfs/ts0/projects/Research_Project-T109423/STAR/STAR.chromosomes.release.fasta

samtools faidx $REF -r $SEL.regions > $SEL.fa

~/minimap2-2.30_x64-linux/minimap2 $REF $SEL.fa > $SEL.paf

module load R/4.3.2-gfbf-2023a
Rscript Biomart.R $SEL.paf

