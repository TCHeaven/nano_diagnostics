```bash
screen -S blast
srun -p long -J blast  --mem-per-cpu 16G --cpus-per-task 4 --pty bash
conda activate blast+

mkdir /home/theaven/scratch/manuscript4/its-blast
cd /home/theaven/scratch/manuscript4/its-blast

cat ../genomes/*/*.fasta > genome.fasta
makeblastdb -in genome.fasta -dbtype nucl -out genome_db

blastn -query LC777859.1.fasta -db genome_db \
  -out hits.txt -outfmt 6 -evalue 1e-20

#LC777859.1      Pod_aph_SCOTT2020_contig_3456   99.765  1275    2       1       1       1275    4066    2793    0.0     2337
#LC777859.1      Pod_aph_DRT72021_contig_2814    99.608  1275    4       1       1       1275    2024    751     0.0     2326
#LC777859.1      Pod_aph_DRT72020_contig_4920    99.608  1275    4       1       1       1275    2024    751     0.0     2326

samtools faidx ../genomes/DRT72020/Pod_aph_DRT72020_clean_renamed.fasta Pod_aph_DRT72020_contig_4920:751-2024 > DRT72020_ITS_28S_region.fasta
samtools faidx ../genomes/DRT72021/Pod_aph_DRT72021_clean_renamed.fasta Pod_aph_DRT72021_contig_2814:751-2024 > DRT72021_ITS_28S_region.fasta
samtools faidx ../genomes/SCOTT2020/Pod_aph_SCOTT2020_clean_renamed.fasta Pod_aph_SCOTT2020_contig_3456:2793-4066 > SCOTT2020_ITS_28S_region.fasta

seqkit seq -r -p DRT72020_ITS_28S_region.fasta > DRT72020_ITS_28S_region.revcomp.fasta
seqkit seq -r -p DRT72021_ITS_28S_region.fasta > DRT72021_ITS_28S_region.revcomp.fasta
seqkit seq -r -p SCOTT2020_ITS_28S_region.fasta > SCOTT2020_ITS_28S_region.revcomp.fasta

cat DRT72020_ITS_28S_region.revcomp.fasta DRT72021_ITS_28S_region.revcomp.fasta SCOTT2020_ITS_28S_region.revcomp.fasta > newseqs-rev.fasta
singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/mafft_7.525--h031d066_1 mafft --thread 1 --6merpair --addfragments newseqs-rev.fasta Bradshaw-alignment.fasta > Bradshaw-plus3-rev.fasta

singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/iqtree_3.0.1--h503566f_0 iqtree3 -s Bradshaw-plus3-rev.fasta \
  -m MFP+MERGE \
  -B 1000 -alrt 1000 \
  -T 4


```
Enrichment of genes within expanded and contracted families across the mildew clade:
```bash
OutDir=~/scratch/nano_diagnostics/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/ortho/analysis/cafe/All_mildew-ara-pruned-species-tree2_results/29102025
OutFile=29102025
Orthogroups=~/scratch/nano_diagnostics/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/ortho/analysis/cafe/Orthogroups.GeneCount-7.tsv
Tree=~/scratch/nano_diagnostics/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/ortho/analysis/cafe/trimmed-6-r8s.tree
lamdaTree=NA
cpu=20
WorkDir=~/scratch/nano_diagnostics/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/ortho/analysis/cafe/All_mildew-ara-pruned-species-tree2_results/29102025

ln -s $Tree $WorkDir/tree.txt
ln -s $Orthogroups $WorkDir/counts.tsv
ln -s $lamdaTree $WorkDir/lamdatree.txt

singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/cafe_5.1.0--h43eeafb_0 cafe5 -i counts.tsv -t tree.txt -e
cp results/Base_error_model.txt error_model_047.txt

singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/cafe_5.1.0--h43eeafb_0 cafe5 -i counts.tsv -t tree.txt -o ${OutDir}/${OutFile}_results -eerror_model_047.txt --cores $cpu

singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven   /home/theaven/git_repos/Containers/python3.sif   bash -lc 'export PATH="$HOME/.local/bin:$PATH"; cafeplotter -i "'"${OutDir}/${OutFile}_results"'" -o "'"${OutDir}/${OutFile}_results"'" '

#FDR
#Perform fdr
WorkDir=/home/theaven/scratch/nano_diagnostics/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/ortho/analysis/cafe/All_mildew-ara-pruned-species-tree2_results/29102025/1
singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/python3.sif python3 /home/theaven/git_repos/Scripts/gruffalo/fdr_patch_branch_probs.py ${WorkDir}/Base_branch_probabilities.tab global
mv ${WorkDir}/Base_branch_probabilities_FDR.tab ${WorkDir}/Base_branch_probabilities.tab

#generate results_summary.tsv with fdr numbers
singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven   /home/theaven/git_repos/Containers/python3.sif   bash -lc 'export PATH="$HOME/.local/bin:$PATH"; cafeplotter -i "'"${WorkDir}"'" -o "'"${WorkDir}"'" '

#generate significant only plot
singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/python3.sif python3 /home/theaven/git_repos/Scripts/gruffalo/make_significant_clade_results_from_summary.py ${WorkDir}/result_summary.tsv 0.05 Base_clade_results_signif.txt
mv Base_clade_results_signif.txt Base_clade_results.txt
singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven   /home/theaven/git_repos/Containers/python3.sif   bash -lc 'export PATH="$HOME/.local/bin:$PATH"; cafeplotter -i "'"${OutDir}/1"'" -o "'"${OutDir}/1"'" --ignore_branch_length --expansion_color green --contraction_color red'

#View 4th line of Base_report.cafe in itol to see node assignments

cafe_results=/home/theaven/scratch/nano_diagnostics/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/ortho/analysis/cafe/All_mildew-ara-pruned-species-tree2_results/29102025/29102025_results/result_summary.tsv
orthofinder_results=/home/theaven/scratch/nano_diagnostics/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/ortho/analysis/orthology/orthofinder/All_mildew-ara/formatted/orthofinder52000-7/Results_all_mildew/Orthogroups/Orthogroups.tsv
OutDir=~/scratch/nano_diagnostics/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/ortho/analysis/cafe/All_mildew-ara-pruned-species-tree2_results/29102025

#Contracted:
#NOTE: a well annotated species from a sister clade should be used for contractions
awk -F'\t' '$2 ~ /(49|50|51|52|53|57|67|68|69|70|72|73|76|77|79|80|81|82|83|84|86|87|88|89|90|91)/ && $4 < 0 {print $1}' $cafe_results | sort | uniq -c | awk '{print $2"\t"$1}' > ${OutDir}/contracted_list.txt
awk -F'\t' '$2 ~ /(49|50|51|52|53|57|67|68|69|70|72|73|76|77|79|80|81|82|83|84|86|87|88|89|90|91)/ && $4 > 0 {print $1}' $cafe_results | sort | uniq -c | awk '{print $2"\t"$1}' > ${OutDir}/expanded_list.txt
grep -vFf <(cut -f1 ${OutDir}/expanded_list.txt) <(cut -f1 ${OutDir}/contracted_list.txt) > ${OutDir}/contracted_only_list.txt
grep -f ${OutDir}/contracted_only_list.txt $orthofinder_results | tr '\t' '\n' | grep '^Aa00' | tr ',' '\n'  | cut -d '|' -f2 > ${OutDir}/contracted.txt

#Expanded
#NOTE: use ingroup focal species for expansions
awk -F'\t' '$2 ~ /(57)/ && $4 < 0 {print $1}' $cafe_results | sort | uniq -c | awk '{print $2"\t"$1}' > ${OutDir}/graminis_contracted_list.txt
awk -F'\t' '$2 ~ /(57)/ && $4 > 0 {print $1}' $cafe_results | sort | uniq -c | awk '{print $2"\t"$1}' > ${OutDir}/graminis_expanded_list.txt
grep -vFf <(cut -f1 ${OutDir}/graminis_contracted_list.txt) <(cut -f1 ${OutDir}/graminis_expanded_list.txt) > ${OutDir}/graminis_expanded_only_list.txt
grep -f ${OutDir}/graminis_expanded_only_list.txt $orthofinder_results | tr '\t' '\n' | grep '^Bg06' | tr ',' '\n'  | cut -d '|' -f2 > ${OutDir}/graminis_expanded.txt

awk -F'\t' '$2 ~ /(53)/ && $4 < 0 {print $1}' $cafe_results | sort | uniq -c | awk '{print $2"\t"$1}' > ${OutDir}/hordei_contracted_list.txt
awk -F'\t' '$2 ~ /(53)/ && $4 > 0 {print $1}' $cafe_results | sort | uniq -c | awk '{print $2"\t"$1}' > ${OutDir}/hordei_expanded_list.txt
grep -vFf <(cut -f1 ${OutDir}/hordei_contracted_list.txt) <(cut -f1 ${OutDir}/hordei_expanded_list.txt) > ${OutDir}/hordei_expanded_only_list.txt
grep -f ${OutDir}/hordei_expanded_only_list.txt $orthofinder_results | tr '\t' '\n' | grep '^Bh15' | tr ',' '\n'  | cut -d '|' -f2 > ${OutDir}/hordei_expanded.txt

awk -F'\t' '$2 ~ /(73)/ && $4 < 0 {print $1}' $cafe_results | sort | uniq -c | awk '{print $2"\t"$1}' > ${OutDir}/leucotricha_contracted_list.txt
awk -F'\t' '$2 ~ /(73)/ && $4 > 0 {print $1}' $cafe_results | sort | uniq -c | awk '{print $2"\t"$1}' > ${OutDir}/leucotricha_expanded_list.txt
grep -vFf <(cut -f1 ${OutDir}/leucotricha_contracted_list.txt) <(cut -f1 ${OutDir}/leucotricha_expanded_list.txt) > ${OutDir}/leucotricha_expanded_only_list.txt
grep -f ${OutDir}/leucotricha_expanded_only_list.txt $orthofinder_results | tr '\t' '\n' | grep '^Pl44' | tr ',' '\n'  | cut -d '|' -f2 > ${OutDir}/leucotricha_expanded.txt

awk -F'\t' '$2 ~ /(77)/ && $4 < 0 {print $1}' $cafe_results | sort | uniq -c | awk '{print $2"\t"$1}' > ${OutDir}/aphanis_contracted_list.txt
awk -F'\t' '$2 ~ /(77)/ && $4 > 0 {print $1}' $cafe_results | sort | uniq -c | awk '{print $2"\t"$1}' > ${OutDir}/aphanis_expanded_list.txt
grep -vFf <(cut -f1 ${OutDir}/aphanis_contracted_list.txt) <(cut -f1 ${OutDir}/aphanis_expanded_list.txt) > ${OutDir}/aphanis_expanded_only_list.txt
grep -f ${OutDir}/aphanis_expanded_only_list.txt $orthofinder_results | tr '\t' '\n' | grep '^Pa36' | tr ',' '\n'  | cut -d '|' -f2 > ${OutDir}/aphanis_expanded.txt

awk -F'\t' '$2 ~ /(70)/ && $4 < 0 {print $1}' $cafe_results | sort | uniq -c | awk '{print $2"\t"$1}' > ${OutDir}/xanthii_contracted_list.txt
awk -F'\t' '$2 ~ /(70)/ && $4 > 0 {print $1}' $cafe_results | sort | uniq -c | awk '{print $2"\t"$1}' > ${OutDir}/xanthii_expanded_list.txt
grep -vFf <(cut -f1 ${OutDir}/xanthii_contracted_list.txt) <(cut -f1 ${OutDir}/xanthii_expanded_list.txt) > ${OutDir}/xanthii_expanded_only_list.txt
grep -f ${OutDir}/xanthii_expanded_only_list.txt $orthofinder_results | tr '\t' '\n' | grep '^Px47' | tr ',' '\n'  | cut -d '|' -f2 > ${OutDir}/xanthii_expanded.txt

awk -F'\t' '$2 ~ /(91)/ && $4 < 0 {print $1}' $cafe_results | sort | uniq -c | awk '{print $2"\t"$1}' > ${OutDir}/necator_contracted_list.txt
awk -F'\t' '$2 ~ /(91)/ && $4 > 0 {print $1}' $cafe_results | sort | uniq -c | awk '{print $2"\t"$1}' > ${OutDir}/necator_expanded_list.txt
grep -vFf <(cut -f1 ${OutDir}/necator_contracted_list.txt) <(cut -f1 ${OutDir}/necator_expanded_list.txt) > ${OutDir}/necator_expanded_only_list.txt
grep -f ${OutDir}/necator_expanded_only_list.txt $orthofinder_results | tr '\t' '\n' | grep '^En17' | tr ',' '\n'  | cut -d '|' -f2 > ${OutDir}/necator_expanded.txt

#Get the annotations used in the orthology to eggnoggmapper:
ls /home/theaven/scratch/nano_diagnostics/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/ortho/analysis/orthology/orthofinder/All_mildew-ara/formatted

```

```bash
conda activate trimmomatic

for ReadDir in /home/theaven/scratch/manuscript4/raw/DNA/*; do
    Freads=($ReadDir/*1.fq.gz)
    Rreads=($ReadDir/*2.fq.gz)
  ls "${Freads[@]}"
  ls "${Rreads[@]}"
    OutDir=$(echo "$ReadDir" | sed 's@raw@dna_qc@')
    Prefix=$(basename "$ReadDir")_
    echo "$Prefix"
  ProgDir=~/git_repos/Wrappers/gruffalo
  Adapters=TruSeq3-PE.fa
    sbatch "$ProgDir/srun_trimmomatic.sh" "${Freads[@]}" "${Rreads[@]}" "$Adapters" "$OutDir" "$Prefix"
done

conda deactivate
```
```bash
#!/bin/bash
#SBATCH -J trimmomatic
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=15G
#SBATCH --time=24:00:00

set -euo pipefail
WORK_DIR="/tmp/${SLURM_JOB_ID}"

# Parse all arguments
ALL_ARGS=("$@")
NUM_ARGS=${#ALL_ARGS[@]}

# Last 3 are fixed
ADAPTER_FILE="${ALL_ARGS[$NUM_ARGS-3]}"
OutDir="${ALL_ARGS[$NUM_ARGS-2]}"
Prefix="${ALL_ARGS[$NUM_ARGS-1]}"

# Remaining are input reads
NUM_READS=$((NUM_ARGS - 3))
HALF=$((NUM_READS / 2))
F_READS=("${ALL_ARGS[@]:0:$HALF}")
R_READS=("${ALL_ARGS[@]:$HALF:$HALF}")

LOGFILE="${Prefix}_trim_log.txt"
F_NO_ADAPT="${Prefix}_F_trim.fq.gz"
R_NO_ADAPT="${Prefix}_R_trim.fq.gz"
UNPAIRED_F_NO_ADAPT="${Prefix}_F_trim_unpaired.fq.gz"
UNPAIRED_R_NO_ADAPT="${Prefix}_R_trim_unpaired.fq.gz"

################################################################

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

zcat "${F_READS[@]}" | gzip > F.fq.gz
zcat "${R_READS[@]}" | gzip > R.fq.gz

trimmomatic PE -phred33 -threads 10 -trimlog "$LOGFILE" \
  F.fq.gz R.fq.gz \
  "$F_NO_ADAPT" "$UNPAIRED_F_NO_ADAPT" \
  "$R_NO_ADAPT" "$UNPAIRED_R_NO_ADAPT" \
  ILLUMINACLIP:"$ADAPTER_FILE":2:30:10 \
  LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36

OutDirF="$OutDir/F"
OutDirR="$OutDir/R"
mkdir -p "$OutDirF" "$OutDirR"

cp "$F_NO_ADAPT" "$UNPAIRED_F_NO_ADAPT" "$OutDirF/"
cp "$R_NO_ADAPT" "$UNPAIRED_R_NO_ADAPT" "$OutDirR/"
cp "$LOGFILE" "$OutDir/"

cd /
rm -r "$WORK_DIR"

echo "Done! Output saved to: $OutDir"
```
```bash
~/nextflow run ~/git_repos/Pipelines/gruffalo/kraken.nf \
 -c ~/git_repos/Pipelines/gruffalo/kraken.config   \
 --genome_dir ~/scratch/manuscript4/analysis   \
 --output_dir ~/scratch/manuscript4/analysis   \
 --kraken_db ~/scratch/manuscript4/db/krakendb   \
 --slurm_log ~/scratch/manuscript4/kraken_slurm_jobs.log   \
 -with-trace ~/scratch/manuscript4/kraken_trace.txt   \
 -with-report ~/scratch/manuscript4/kraken_report.html   \
 -with-timeline ~/scratch/manuscript4/kraken_timeline.html   \
 -with-dag ~/scratch/manuscript4/kraken_dag.png 2>&1 > screen_log.txt

srun -p long -J kraken --mem 80G -c 20 --pty bash
cd ~/scratch/manuscript4/db/krakendb-mildews
wget https://github.com/R-Wright-1/peptides/archive/refs/heads/master.zip
unzip master.zip
conda create -n biopython
conda activate biopython
conda update -n base conda
conda install -c conda-forge biopython
conda install pandas
mkdir tmp123
cd tmp123
DBNAME=/home/theaven/scratch/manuscript4/db/krakendb-mildews/nt-mildew
singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/kraken2.1.3.sif kraken2-build --download-taxonomy --use-ftp --threads 20 --db $DBNAME 
python ../peptides-master/download_domain_09082025.py --domain fungi --ext dna

mkdir ~/scratch/manuscript4/db/krakendb-mildews/limited-mildew
cd ~/scratch/manuscript4/db/krakendb-mildews/limited-mildew
#download GCA_000208805.1_ASM20880v1_genomic.fnanGCA_002918395.1_ASM291839v1_genomic.fna GCA_003610855.1_ASM361085v1_genomic.fna GCA_003611235.1_ASM361123v1_genomic.fna GCA_003957845.1_ASM395784v1_genomic.fna GCA_013170925.1_ASM1317092v1_genomic.fna GCA_019455505.1_ASM1945550v1_genomic.fna GCA_019455665.1_ASM1945566v1_genomic.fna GCA_028751805.1_YZU_CsPM_1.2_genomic.fna GCA_030378345.1_ASM3037834v1_genomic.fna GCA_044715735.1_Tps_LRv5b_mtDNAv350_genomic.fna GCA_050494535.1_ASM5049453v1_genomic.fna GCA_900239735.1_BGH_DH14_v4_genomic.fna GCA_902706625.1_Albugo_laibachii_Nc14_genomic.fna GCA_905067625.1_Bgtriticale_THUN12_genome_v1_2_genomic.fna

sed 's/ Timema/|kraken:taxid|170557 Timema/g' GCA_044715735.1_Tps_LRv5b_mtDNAv350_genomic.fna > GCA_044715735.1_Tps_LRv5b_mtDNAv350_genomic.mod.fna
sed 's/ Timema/|kraken:taxid|61476 Timema/g' GCA_050494535.1_ASM5049453v1_genomic.fna > GCA_050494535.1_ASM5049453v1_genomic.mod.fna
sed 's/ Albugo/|kraken:taxid|653948 Albugo/g' GCA_902706625.1_Albugo_laibachii_Nc14_genomic.fna > GCA_902706625.1_Albugo_laibachii_Nc14_genomic.mod.fna
sed 's/ Erysiphe/|kraken:taxid|36044 Erysiphe/g' GCA_000208805.1_ASM20880v1_genomic.fna > GCA_000208805.1_ASM20880v1_genomic.mod.fna
sed 's/ Erysiphe/|kraken:taxid|225359 Erysiphe/g' GCA_002918395.1_ASM291839v1_genomic.fna > GCA_002918395.1_ASM291839v1_genomic.mod.fna
sed 's/ Erysiphe/|kraken:taxid|212602 Erysiphe/g' GCA_003610855.1_ASM361085v1_genomic.fna > GCA_003610855.1_ASM361085v1_genomic.mod.fna
sed 's/ Golovinomyces/|kraken:taxid|62708 Golovinomyces/g' GCA_003611235.1_ASM361123v1_genomic.fna > GCA_003611235.1_ASM361123v1_genomic.mod.fna
sed 's/ Oidium/|kraken:taxid|299130 Oidium/g' GCA_003957845.1_ASM395784v1_genomic.fna > GCA_003957845.1_ASM395784v1_genomic.mod.fna
sed 's/ Podosphaera/|kraken:taxid|79249 Podosphaera/g' GCA_013170925.1_ASM1317092v1_genomic.fna > GCA_013170925.1_ASM1317092v1_genomic.mod.fna
sed 's/ Pleochaeta/|kraken:taxid|57462 Pleochaeta/g' GCA_019455505.1_ASM1945550v1_genomic.fna > GCA_019455505.1_ASM1945550v1_genomic.mod.fna
sed 's/ Phyllactinia/|kraken:taxid|57460 Phyllactinia/g' GCA_019455665.1_ASM1945566v1_genomic.fna > GCA_019455665.1_ASM1945566v1_genomic.mod.fna
sed 's/ Podosphaera/|kraken:taxid|135283 Podosphaera/g' GCA_028751805.1_YZU_CsPM_1.2_genomic.fna > GCA_028751805.1_YZU_CsPM_1.2_genomic.mod.fna
sed 's/ Podosphaera/|kraken:taxid|62727 Podosphaera/g' GCA_030378345.1_ASM3037834v1_genomic.fna > GCA_030378345.1_ASM3037834v1_genomic.mod.fna
sed 's/ Blumeria/|kraken:taxid|2867405 Blumeria/g' GCA_900239735.1_BGH_DH14_v4_genomic.fna > GCA_900239735.1_BGH_DH14_v4_genomic.mod.fna
sed 's/ Blumeria/|kraken:taxid|62688 Blumeria/g' GCA_905067625.1_Bgtriticale_THUN12_genome_v1_2_genomic.fna > GCA_905067625.1_Bgtriticale_THUN12_genome_v1_2_genomic.mod.fna

for assembly in $(ls *.mod.fna); do
singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/kraken2.1.3.sif kraken2-build --add-to-library $assembly --db $DBNAME 2>&1 | tee -a 2.log
done

singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/kraken2.1.3.sif kraken2-build --build --db $DBNAME 2>&1 | tee -a 3.log
singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/kraken2.1.3.sif kraken2-build --clean --threads 20 --db $DBNAME

srun -p long -J kraken --mem 16G -c 4 --pty bash
DBNAME=/home/theaven/scratch/manuscript4/db/krakendb-mildews/nt-mildew
for assembly in $(ls ~/scratch/manuscript4/analysis/*.fasta); do
  outprefix=$(dirname $assembly)/$(basename $assembly | cut -d '_' -f3)/kraken2/$(basename $assembly | cut -d '_' -f3)_$(basename $DBNAME)
    singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/kraken2.1.3.sif kraken2 \
      --db /home/theaven/scratch/manuscript4/db/krakendb-mildews/nt-mildew \
      --threads 4 \
      --output ${outprefix}_output.kraken.txt \
      --unclassified-out ${outprefix}_unclassified-out.kraken.txt \
      --classified-out ${outprefix}_classified-out.kraken.txt \
      --report ${outprefix}_report.kraken.txt \
      --use-names \
      $assembly 
done
```
Blobtools with the fcs filtered assemblies
```bash
#Collect the assemblies
ln -s ~/projects/niab/theaven/genomes/Pod_aph_DRT72020_clean_renamed.fasta genomes/DRT72020/.
ln -s ~/projects/niab/theaven/genomes/Pod_aph_DRT72021_clean_renamed.fasta genomes/DRT72021/.
ln -s ~/projects/niab/theaven/genomes/Pod_aph_SCOTT2020_clean_renamed.fasta genomes/SCOTT2020/.
ln -s ~/projects/niab/theaven/genomes/Pod_leu_OGB2019_clean_renamed.fasta genomes/OGB2019/.
ln -s ~/projects/niab/theaven/genomes/Pod_leu_OGB2021_clean_renamed.fasta genomes/OGB2021/.
ln -s ~/projects/niab/theaven/genomes/Pod_leu_OGBp112020_clean_renamed.fasta genomes/P112020/.

#Alignment of reads to the assemblies
srun -p long -J bowtie  --mem-per-cpu 8G --cpus-per-task 8 --pty bash

conda activate bowtie2
for genome in $(ls ~/scratch/manuscript4/genomes/DRT*/*.fasta); do
ID=$(echo $genome | cut -d '/' -f7)
OutDir=~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/bowtie2
mkdir -p "$OutDir"
cd "$OutDir" 
bowtie2-build "$genome" "${ID}_index"
FReads=$(ls ~/scratch/manuscript4/dna_qc/DNA/${ID}/F/*F_trim.fq.gz 2>/dev/null | paste -sd "," -)
RReads=$(ls ~/scratch/manuscript4/dna_qc/DNA/${ID}/R/*R_trim.fq.gz 2>/dev/null | paste -sd "," -)
Unpaired_FReads=$(ls ~/scratch/manuscript4/dna_qc/DNA/${ID}/F/*F_trim_unpaired.fq.gz 2>/dev/null | paste -sd "," -)
Unpaired_RReads=$(ls ~/scratch/manuscript4/dna_qc/DNA/${ID}/R/*R_trim_unpaired.fq.gz 2>/dev/null | paste -sd "," -)
bowtie2 \
-x "${ID}_index" -p 8\
-1 $FReads \
-2 $RReads \
-U $Unpaired_FReads,$Unpaired_RReads \
--un-gz "${ID}_s.fq.gz" \
--un-conc-gz "${ID}_fr.fq.gz" \
-S "${ID}_0.sam" | tee -a report_0.txt
cd ~/scratch/manuscript4
done
conda deactivate
exit
exit
echo finished

#blast
conda activate blast+
for Assembly in $(ls ~/scratch/manuscript4/genomes/*/*.fasta); do
ID=$(echo $Assembly | cut -d '/' -f7)
Outfile=${ID}.vs.nt.mts1.hsp1.1e25.megablast.out
db=~/scratch/manuscript4/db/ncbi-nt/nt
OutDir=~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/blastn2
mkdir -p $OutDir
ProgDir=~/git_repos/Wrappers/gruffalo
sbatch $ProgDir/blastn.sh $Assembly $OutDir/$Outfile $db
done
conda deactivate
#3253848-53,3256482-7

conda activate blast+
for Assembly in $(ls ~/scratch/manuscript4/genomes/*/*.fasta); do
ID=$(echo $Assembly | cut -d '/' -f7)
Outfile=${ID}.vs.nt.mts1.hsp1.1e25.megablast.out
db=/mnt/shared/datasets/databases/ncbi/nt
OutDir=~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/blastn3
mkdir -p $OutDir
ProgDir=~/git_repos/Wrappers/gruffalo
sbatch $ProgDir/blastn.sh $Assembly $OutDir/$Outfile $db
done
conda deactivate
#3883648-53

cd /home/theaven/scratch/manuscript4/db/blobtools
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P data/
cd data
tar zxf taxdump.tar.gz -C . nodes.dmp names.dmp
#./blobtools nodesdb --nodes nodes.dmp --names names.dmp

#diamond
conda activate diamond
for Assembly in $(ls ~/scratch/manuscript4/genomes/*/*.fasta); do
ID=$(echo $Assembly | cut -d '/' -f7)
Outfile=${ID}
Database=/mnt/shared/datasets/databases/diamond/nr.dmnd
Max_target=5
OutDir=~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/diamond
mkdir -p $OutDir
ProgDir=~/git_repos/Wrappers/gruffalo
sbatch $ProgDir/run_diamond_x.sh $Assembly $Database $OutDir $Outfile $Max_target
done
conda deactivate
#3883426-31

#busco
conda activate busco
for Assembly in $(ls ~/scratch/manuscript4/genomes/*/*.fasta); do
ID=$(echo $Assembly | cut -d '/' -f7)
Outfile=${ID}_busco
OutDir=~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/busco
mkdir -p $OutDir
Database=leotiomycetes_odb10
ProgDir=~/git_repos/Wrappers/gruffalo
sbatch $ProgDir/run_busco_keep.sh $Assembly $Database $OutDir $Outfile
done
conda deactivate
#3521412-7

#tiara
for Assembly in $(ls ~/scratch/manuscript4/genomes/*/*.fasta); do
ID=$(echo $Assembly | cut -d '/' -f7)
OutPrefix=${ID}
OutDir=~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/tiara
ProgDir=~/git_repos/Wrappers/gruffalo
mkdir $OutDir
sbatch $ProgDir/run_tiara.sh $Assembly $OutDir $OutPrefix
done
#3253859-64

for Assembly in $(ls /home/theaven/scratch/manuscript4/genomes/DRT72020/Pod_aph_DRT72020_clean_renamed.fasta /home/theaven/scratch/manuscript4/genomes/DRT72021/Pod_aph_DRT72021_clean_renamed.fasta /home/theaven/scratch/manuscript4/genomes/SCOTT2020/Pod_aph_SCOTT2020_clean_renamed.fasta); do
ID=$(echo $Assembly | cut -d '/' -f7)
Record_type=contig
MappingFile=$(ls ~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/bowtie2/*.sam)
BlastFile=$(ls ~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/blastn3/*.out2)
DiamondFile=$(ls ~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/diamond/*.out)
BUSCOFile=$(ls ~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/busco/run_leotiomycetes_odb10/full_table.tsv)
Tiara=$(ls ~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/tiara/${ID}.tiara)
OutDir=~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/02112025
OutPrefix=${ID}
Genus=Podosphaera
Species=aphanis
TaxID=79252
alias=$ID
ProgDir=~/git_repos/Wrappers/gruffalo
mkdir $OutDir
#sbatch $ProgDir/run_blobtools.sh $Assembly $MappingFile $BlastFile $OutDir/1.1.1 $OutPrefix
sbatch $ProgDir/run_blobtoolkit4.4.5.sh $Assembly $Record_type $MappingFile $BlastFile $DiamondFile $BUSCOFile $Tiara $OutDir $OutPrefix $Genus $Species $TaxID $alias
#sbatch $ProgDir/run_blobtoolkit4.4.5.sh $Assembly $Record_type $MappingFile $BlastFile $DiamondFile $BUSCOFile $Tiara $OutDir/no-tiara $OutPrefix $Genus $Species $TaxID $alias
#sbatch $ProgDir/run_blobtoolkit4.4.5.sh $Assembly $Record_type $MappingFile $BlastFile $DiamondFile $BUSCOFile $Tiara $OutDir/tiara-only $OutPrefix $Genus $Species $TaxID $alias
done
#3556117-9
#3879682-4
#4222140-2

for Assembly in $(ls /home/theaven/scratch/manuscript4/genomes/OGB2019/Pod_leu_OGB2019_clean_renamed.fasta /home/theaven/scratch/manuscript4/genomes/OGB2021/Pod_leu_OGB2021_clean_renamed.fasta /home/theaven/scratch/manuscript4/genomes/P112020/Pod_leu_OGBp112020_clean_renamed.fasta); do
ID=$(echo $Assembly | cut -d '/' -f7)
Record_type=contig
MappingFile=$(ls ~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/bowtie2/*.sam)
BlastFile=$(ls ~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/blastn3/*.out2)
DiamondFile=$(ls ~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/diamond/*.out)
BUSCOFile=$(ls ~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/busco/run_leotiomycetes_odb10/full_table.tsv)
Tiara=$(ls ~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/tiara/${ID}.tiara)
OutDir=~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/02112025
OutPrefix=${ID}
Genus=Podosphaera
Species=leucotricha
TaxID=79249
alias=$ID
ProgDir=~/git_repos/Wrappers/gruffalo
mkdir $OutDir
#sbatch $ProgDir/run_blobtools.sh $Assembly $MappingFile $BlastFile $OutDir/1.1.1 $OutPrefix
sbatch $ProgDir/run_blobtoolkit4.4.5.sh $Assembly $Record_type $MappingFile $BlastFile $DiamondFile $BUSCOFile $Tiara $OutDir $OutPrefix $Genus $Species $TaxID $alias
#sbatch $ProgDir/run_blobtoolkit4.4.5.sh $Assembly $Record_type $MappingFile $BlastFile $DiamondFile $BUSCOFile $Tiara $OutDir/no-tiara $OutPrefix $Genus $Species $TaxID $alias
#sbatch $ProgDir/run_blobtoolkit4.4.5.sh $Assembly $Record_type $MappingFile $BlastFile $DiamondFile $BUSCOFile $Tiara $OutDir/tiara-only $OutPrefix $Genus $Species $TaxID $alias
done
#3521771-3
#3522168-70
#3522189-91
#3879645-7
#3879785-90
#3880278-83
#3881699-704
#4222147-9
#4272943-8
#4494180-5

#Enter \\wsl$\Ubuntu\home\tcheaven into the address bar of file explorer

apptainer exec blobtoolkit.sif blobtools host /home/tcheaven/BlobDirs
apptainer exec blobtoolkit.sif blobtools host /home/tcheaven/no-tiara
apptainer exec blobtoolkit.sif blobtools host /home/tcheaven/tiara-only
apptainer exec blobtoolkit.sif blobtools host /home/tcheaven/down_12112025

singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/python3.sif python3 /home/theaven/git_repos/Scripts/NBI/seq_rm.py \
        --id_file ogb2021_rm.txt \
        --input /home/theaven/scratch/manuscript4/genomes/OGB2021/Pod_leu_OGB2021_clean_renamed.fasta \
        --output /home/theaven/scratch/manuscript4/genomes/OGB2021/Pod_leu_OGB2021-2_clean_renamed.fasta

singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/python3.sif python3 /home/theaven/git_repos/Scripts/NBI/seq_rm.py \
        --id_file ogb2019_rm.txt \
        --input /home/theaven/scratch/manuscript4/genomes/OGB2019/Pod_leu_OGB2019_clean_renamed.fasta \
        --output /home/theaven/scratch/manuscript4/genomes/OGB2019/Pod_leu_OGB2019-2_clean_renamed.fasta

singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/python3.sif python3 /home/theaven/git_repos/Scripts/NBI/seq_rm.py \
        --id_file scott2020_rm.txt \
        --input /home/theaven/scratch/manuscript4/genomes/SCOTT2020/Pod_aph_SCOTT2020_clean_renamed.fasta \
        --output /home/theaven/scratch/manuscript4/genomes/SCOTT2020/Pod_aph_SCOTT2020-2_clean_renamed.fasta

singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/python3.sif python3 /home/theaven/git_repos/Scripts/NBI/seq_rm.py \
        --id_file drt72021_rm.txt \
        --input /home/theaven/scratch/manuscript4/genomes/DRT72021/Pod_aph_DRT72021_clean_renamed.fasta \
        --output /home/theaven/scratch/manuscript4/genomes/DRT72021/Pod_aph_DRT72021-2_clean_renamed.fasta

#bowtie
srun -p long -J bowtie  --mem-per-cpu 8G --cpus-per-task 8 --pty bash
conda activate bowtie2
for genome in $(ls ~/scratch/manuscript4/genomes/O*/*-2_clean_renamed_fix.fasta ~/scratch/manuscript4/genomes/S*/*-2_clean_renamed_fix.fasta); do
ID=$(echo $genome | cut -d '/' -f7)
OutDir=~/scratch/manuscript4/analysis/${ID}-2/blobtools-fcs/bowtie2
mkdir -p "$OutDir"
cd "$OutDir" 
bowtie2-build "$genome" "${ID}_index"
FReads=$(ls ~/scratch/manuscript4/dna_qc/DNA/${ID}/F/*F_trim.fq.gz 2>/dev/null | paste -sd "," -)
RReads=$(ls ~/scratch/manuscript4/dna_qc/DNA/${ID}/R/*R_trim.fq.gz 2>/dev/null | paste -sd "," -)
Unpaired_FReads=$(ls ~/scratch/manuscript4/dna_qc/DNA/${ID}/F/*F_trim_unpaired.fq.gz 2>/dev/null | paste -sd "," -)
Unpaired_RReads=$(ls ~/scratch/manuscript4/dna_qc/DNA/${ID}/R/*R_trim_unpaired.fq.gz 2>/dev/null | paste -sd "," -)
ls $OutDir
ls $FReads $RReads $Unpaired_RReads $Unpaired_FReads
bowtie2 \
-x "${ID}_index" -p 8 \
-1 $FReads \
-2 $RReads \
-U $Unpaired_FReads,$Unpaired_RReads \
--un-gz "${ID}_s.fq.gz" \
--un-conc-gz "${ID}_fr.fq.gz" \
-S "${ID}_0.sam" | tee -a report_0.txt
cd ~/scratch/manuscript4
done

for genome in $(ls ~/scratch/manuscript4/genomes/D*/*-2_clean_renamed_fix.fasta); do
ID=$(echo $genome | cut -d '/' -f7)
OutDir=~/scratch/manuscript4/analysis/${ID}-2/blobtools-fcs/bowtie2
mkdir -p "$OutDir"
cd "$OutDir" 
bowtie2-build "$genome" "${ID}_index"
FReads=$(ls ~/scratch/manuscript4/dna_qc/DNA/DRCT72021/F/*F_trim.fq.gz 2>/dev/null | paste -sd "," -)
RReads=$(ls ~/scratch/manuscript4/dna_qc/DNA/DRCT72021/R/*R_trim.fq.gz 2>/dev/null | paste -sd "," -)
Unpaired_FReads=$(ls ~/scratch/manuscript4/dna_qc/DNA/DRCT72021/F/*F_trim_unpaired.fq.gz 2>/dev/null | paste -sd "," -)
Unpaired_RReads=$(ls ~/scratch/manuscript4/dna_qc/DNA/DRCT72021/R/*R_trim_unpaired.fq.gz 2>/dev/null | paste -sd "," -)
ls $OutDir
ls $FReads $RReads $Unpaired_RReads $Unpaired_FReads
bowtie2 \
-x "${ID}_index" -p 8 \
-1 $FReads \
-2 $RReads \
-U $Unpaired_FReads,$Unpaired_RReads \
--un-gz "${ID}_s.fq.gz" \
--un-conc-gz "${ID}_fr.fq.gz" \
-S "${ID}_0.sam" | tee -a report_0.txt
cd ~/scratch/manuscript4
done
conda deactivate
exit
exit
echo finished

#blast
conda activate blast+
for Assembly in $(ls ~/scratch/manuscript4/genomes/*/*-2_clean_renamed_fix.fasta); do
ID=$(echo $Assembly | cut -d '/' -f7)-2
Outfile=${ID}.vs.nt.mts1.hsp1.1e25.megablast.out
db=/mnt/shared/datasets/databases/ncbi/nt
OutDir=~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/blastn3
mkdir -p $OutDir
ProgDir=~/git_repos/Wrappers/gruffalo
sbatch $ProgDir/blastn.sh $Assembly $OutDir/$Outfile $db
done
conda deactivate
#4900126-9

#diamond
conda activate diamond
for Assembly in $(ls ~/scratch/manuscript4/genomes/*/*-2_clean_renamed_fix.fasta); do
ID=$(echo $Assembly | cut -d '/' -f7)-2
Outfile=${ID}
Database=/mnt/shared/datasets/databases/diamond/nr.dmnd
Max_target=5
OutDir=~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/diamond
mkdir -p $OutDir
ProgDir=~/git_repos/Wrappers/gruffalo
sbatch $ProgDir/run_diamond_x.sh $Assembly $Database $OutDir $Outfile $Max_target
done
conda deactivate
#4900121-4

conda activate busco
for Assembly in $(ls ~/scratch/manuscript4/genomes/*/*-2_clean_renamed_fix.fasta); do
ID=$(echo $Assembly | cut -d '/' -f7)-2
Outfile=${ID}_busco
OutDir=~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/busco
mkdir -p $OutDir
Database=leotiomycetes_odb10
ProgDir=~/git_repos/Wrappers/gruffalo
sbatch $ProgDir/run_busco_keep.sh $Assembly $Database $OutDir $Outfile
done
conda deactivate
#4900111-4

for Assembly in $(ls ~/scratch/manuscript4/genomes/*/*-2_clean_renamed_fix.fasta); do
ID=$(echo $Assembly | cut -d '/' -f7)-2
OutPrefix=${ID}
OutDir=~/scratch/manuscript4/analysis/${ID}/blobtools-fcs/tiara
ProgDir=~/git_repos/Wrappers/gruffalo
mkdir $OutDir
sbatch $ProgDir/run_tiara.sh $Assembly $OutDir $OutPrefix
done
#4900106-9

for file in $(ls ~/scratch/manuscript4/analysis/*/blobtools-fcs/busco/short_summary.*.txt); do
echo $file
cat $file
done
```






























```R
setwd("E:/R")

library(readxl)

df <- read_excel("P_aphanis_THeavenDRCT72020_annotations_master.xlsx")
df <- read_excel("P_aphanis_THeavenDRCT72021_annotations_master.xlsx")
df <- read_excel("P_aphanis_THeavenSCOTT2020_annotations_master.xlsx")
df <- read_excel("P_leucotricha_THeavenpOGB2019_annotations_master.xlsx")
df <- read_excel("P_leucotricha_THeavenp11_annotations_master.xlsx")
df <- read_excel("P_leucotricha_THeavenpOGB2021_annotations_master.xlsx")

subset_data <- df[df$te_group == 'Any TE', ]
subset_data$eka <- "Non EKA"
subset_data <- subset_data %>%
  mutate(eka = ifelse(grepl("BgtAVRk1|BgtAVRa10", effector_matches), "EKA", eka))
subset_data$eka[subset_data$is_secreted == 1] <- "Non EKA"
subset_data$ralph <- "Non RALPH"
subset_data <- subset_data %>%
  mutate(ralph = ifelse(grepl("BghBEC1011|BgtAVRa10|BgAVRA13|BgtAvrPm2|BgtSvrPm3a1f1", effector_matches), "RALPH", ralph))
sum(subset_data$eka == "EKA")
sum(subset_data$ralph == "RALPH")
```
DRCT72020 103 176
DRCT72021 91  161
SCOTT2020 72 133
OGB2019 155 429
p11 158 430
OGB2021 155 430

```bash
cd /home/theaven/scratch/manuscript4/mcag

for file in *.fsa; do
  [[ "$file" == GCA_* ]] && continue  # Skip files starting with GCA_
  
  awk '{ if ($0 ~ /^>/) { sub(/^>/, ">MCAG_"); print } else { print } }' "$file" > "${file}.tmp" && mv "${file}.tmp" "$file"
done

cat *.fsa >> proteins.faa

srun -p medium --mem 16G -c 4 --pty bash
conda activate blast+
makeblastdb -in proteins.faa -input_type fasta -dbtype prot  -title mcags_db  -parse_seqids -out mcags_db

/home/theaven/projects/niab/theaven/gene_pred/mildewreannotation/Oidiumheveae_GCA_003957845/braker/final_genes_renamed.pep.fasta





for file in $(ls /home/theaven/projects/niab/theaven/mcag/XP_752655_protein_results); do
  mcag_ID=$(basename "$file" | sed 's/_protein_results//')

  while read -r line; do
    # Extract genome ID from the second field
    genome_ID=$(echo "$line" | awk '{print $2}' | cut -d'_' -f3,4)
echo $genome_ID
    # Find the correct gene FASTA file using genome_ID
    gene_file=$(find /home/theaven/projects/niab/theaven/gene_pred/mildewreannotation/ -type f -name "final_genes_renamed.pep.fasta" | grep "$genome_ID")
echo $gene_file
    # Extract the 5th part of the second field (sequence ID) for lookup
    echo "$line" | awk '{print $2}' | awk -F'_' '{print $NF}' > temp_id.txt
cat temp_id.txt
    # Extract sequence from gene FASTA using Python script
    singularity exec \
      --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven \
      /home/theaven/git_repos/Containers/python3.sif \
      python3 /home/theaven/git_repos/Scripts/NBI/seq_get.py \
        --id_file temp_id.txt \
        --input "$gene_file" \
        --output temp.faa

    # Append results to output FASTA with renamed headers
    awk -v gid="$genome_ID" '{ if ($0 ~ /^>/) { sub(/^>/, ">" gid "_"); print } else { print } }' temp.faa >> "${mcag_ID}.faa"

  done < "$file"
done

input_dir="/home/theaven/projects/niab/theaven/mcag"
for file in "$input_dir"/*_protein_results; do
  while read -r line; do
    echo "$line" | awk '{
      split($2, a, "_")
      for (i = 1; i < length(a); i++) {
        printf i == 1 ? a[i] : "_" a[i]
      }
      print ""
    }'
  done < "$file"
done | grep -v "_file_1_file_1" | sort -u

Am_c_GCA_003019875
Am_c_GCA_018167515
Ar_a_GCA_003988855
A_s_GCA_000328965
B_g_GCA_000417025
B_g_GCA_000417865
B_g_GCA_000418435
B_g_GCA_000441875
B_g_GCA_900519115
B_g_GCA_905067625
B_h_GCA_000401675
B_h_GCA_900237765
B_h_GCA_900239735
B_h_GCA_900638725
C_a_GCA_002276475
D_b_GCA_000298775
E_alph
E_neca_GCA_000798715
E_neca_GCA_000798735
E_neca_GCA_000798755
E_neca_GCA_000798775
E_neca_GCA_000798795
E_neca_GCA_016906895
E_neol_GCA_003610855
E_pi_GCA_000208805
E_pi_GCA_000214055
E_pu_GCA_002918395
G_c_GCA_003611195
G_c_GCA_003611215
G_c_GCA_003611235
Gl_i_GCA_000409485
G_m_GCA_006912115
G_o
L_t_CADEPA01
M_s_GCA_001500285
N_a_GCA_003988965
O_h_GCA_003957845
Oi_m_GCA_000827325
P_a_frg2020
P_a_frg2021
P_a_rub2020
P_c_GCA_018398735
Phi_s-GCA_900073065
Phy_m_GCA_019455665
P_l2019
P_l2020
P_l2021
P_l_ganan
Pl_s_GCA_019455505
P_x_GCA_010015925
P_x_GCA_014884795
Scereviseae_GCF_000146045

   ```
```bash
cd /home/theaven/scratch/manuscript4/mcag
#theaven@gruffalo:~/scratch/manuscript4/mcag$ ls -lh GCA*fsa
#-rw-r--r-- 1 theaven theaven 6.4M Aug  4 07:04 GCA_000002495.2_MG8_genomic.fsa
#-rw-r--r-- 1 theaven theaven 7.5M Aug  4 07:03 GCA_000143535.4_ASM14353v4_genomic.fsa
#-rw-r--r-- 1 theaven theaven 3.3M Aug  4 07:04 GCA_000146045.2_R64_genomic.fsa
#-rw-r--r-- 1 theaven theaven 6.2M Aug  4 07:03 GCA_000146945.2_ASM14694v2_genomic.fsa
#-rw-r--r-- 1 theaven theaven 7.8M Aug  4 07:03 GCA_001672515.1_ASM167251v1_genomic.fsa

for file in *.fsa; do
  [[ "$file" == GCA_* ]] && continue  # Skip files starting with GCA_
  
  awk '{ if ($0 ~ /^>/) { sub(/^>/, ">MCAG_"); print } else { print } }' "$file" > "${file}.tmp" && mv "${file}.tmp" "$file"
done

cat *.fsa >> proteins.faa

srun -p medium --mem 16G -c 4 --pty bash
conda activate blast+
makeblastdb -in proteins.faa -input_type fasta -dbtype prot  -title mcags_db  -parse_seqids -out mcags_db

#Give species unique names to each predicted gene:
for file in $(find /home/theaven/projects/niab/theaven/gene_pred/mildewreannotation/ -type f -name "final_genes_renamed.pep.fasta"); do
ID=$(echo $file | cut -d '/' -f9 | cut -d '.' -f1 | rev | cut -d '_' -f1,2 | rev | sed 's@_AssemblyScaffolds@@g' | sed 's@phialocephalasubalpina-@@g')
echo $ID
cp $file ./${ID}_$(basename $file)
sed -i "s/>/>${ID}_/g" ${ID}_$(basename $file)
sed -i "s/file_1_file_1_//g" ${ID}_$(basename $file)
done

for file in $(ls *final_genes_renamed.pep.fasta); do
seqkit rmdup $file > temp.fasta && mv temp.fasta $file
done

#make a blast database from the proteomes and query each putative MCAG against it:
mkdir db
cat *final_genes_renamed.pep.fasta > db/db.faa
cd db
makeblastdb -in db.faa -input_type fasta -dbtype prot  -title mildew+db  -parse_seqids -out mildew+db

input_dir="/home/theaven/projects/niab/theaven/mcag/prot"
output_dir="/home/theaven/scratch/manuscript4/mcag/2/db"

for file in "$input_dir"/*.fsa; do
  filename=$(basename "$file")
  if [[ "$filename" == "out_genomes_protein_set.fsa" ]]; then
    continue
  fi
  name="${filename%.fsa}" 
  blastp -query "$file" \
         -db mildew+db \
         -out "${output_dir}/${name}_results" \
         -evalue 1e-5 \
         -outfmt 6 \
         -num_threads 4
done

for file in *_results; do
  while read -r line; do
Genome_ID=$(echo "$line" | awk '{print $2}' | awk -F'_' '{
  for (i = 1; i < NF; i++) {
    printf i == 1 ? $i : "_"$i
  }
  print ""
}')

Genome_ID2=$(ls /home/theaven/projects/niab/theaven/gene_pred/mildewreannotation/*/braker/final_genes_renamed.pep.fasta | grep "$Genome_ID" | cut -d '/' -f9)
echo $Genome_ID2
echo .
echo "$line" | sed "s@${Genome_ID}@${Genome_ID2}@g"  >> ${file}_2

done < "$file"
done

for file in $(ls /home/theaven/scratch/manuscript4/mcag/2/db/*_protein_results_2 | grep -v 'XP_752655'); do
  mcag_ID=$(basename "$file" | sed 's/_protein_results//')

  while read -r line; do
    # Extract genome ID from the second field
    Genome_ID=$(echo "$line" | awk '{print $2}' | awk -F'_' '{
  for (i = 1; i < NF; i++) {
    printf i == 1 ? $i : "_"$i
  }
  print ""
}')
    echo $Genome_ID
    # Find the correct gene FASTA file using genome_ID
    gene_file=$(ls /home/theaven/projects/niab/theaven/gene_pred/mildewreannotation/*/braker/final_genes_renamed.pep.fasta | grep "$Genome_ID")
echo $gene_file
    # Extract the 5th part of the second field (sequence ID) for lookup
    echo "$line" | awk '{print $2}' | awk -F'_' '{print $NF}' > temp_id.txt
cat temp_id.txt
    # Extract sequence from gene FASTA using Python script
    singularity exec \
      --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven \
      /home/theaven/git_repos/Containers/python3.sif \
      python3 /home/theaven/git_repos/Scripts/NBI/seq_get.py \
        --id_file temp_id.txt \
        --input "$gene_file" \
        --output temp.faa

    # Append results to output FASTA with renamed headers
    awk -v gid="$Genome_ID" '{ if ($0 ~ /^>/) { sub(/^>/, ">" gid "_"); print } else { print } }' temp.faa >> "${mcag_ID}.faa"

  done < "$file"
done

screen
srun -p long --mem 16G -c 4 --pty bash
conda activate blast+
cd /home/theaven/scratch/manuscript4/mcag
input_dir="/home/theaven/scratch/manuscript4/mcag/2/db"
output_dir="/home/theaven/scratch/manuscript4/mcag"
for file in "$input_dir"/*_2.faa; do
  filename=$(basename "$file")
  name="${filename%_2.faa}" 
  blastp -query "$file" \
         -db mcags_db \
         -out "${output_dir}/${name}_results_mcag" \
         -evalue 1e-5 \
         -outfmt 6 \
         -num_threads 4
done

for file in "$input_dir"/*_2.faa; do
  filename=$(basename "$file")
  name="${filename%_2.faa}" 
  blastp -query "$file" \
         -db mcags_db \
         -out "${output_dir}/${name}_results_mcag_top" \
         -evalue 1e-5 \
         -outfmt 6 \
         -max_target_seqs 1 \
         -num_threads 4
done

for file in *_results_mcag; do
  Out=${file}_mcag
  grep 'MCAG_' $file > $Out
done

for file in *_results_mcag_top; do
  Out=${file}_mcag
  grep 'MCAG_' $file > $Out
done


#Collect results into one file:
echo Gene: > out99_mcag.tsv

for file in $(find /home/theaven/projects/niab/theaven/gene_pred/mildewreannotation/ -type f -name "final_genes_renamed.pep.fasta"); do
ID=$(echo $file | cut -d '/' -f9 | cut -d '.' -f1 )
#ID=$(echo $file | cut -d '/' -f9 | cut -d '.' -f1 | rev | cut -d '_' -f1,2 | rev )
echo $ID
sed -i "\$s/\$/ $ID/" out99_mcag.tsv
done

for file in $(ls *_results_mcag_mcag); do
gene=$(echo $file | sed 's@_results_mcag_mcag@@g')
echo $gene >> out99_mcag.tsv
for file2 in $(find /home/theaven/projects/niab/theaven/gene_pred/mildewreannotation/ -type f -name "final_genes_renamed.pep.fasta"); do
ID=$(echo $file2 | cut -d '/' -f9 | cut -d '.' -f1 | rev | cut -d '_' -f1,2 | rev )
x=$(grep -m 1 $ID $file | awk '{print $11}')
if [ -z "$x" ]; then
    y=missing
else
    y=$x
fi
sed -i "\$s/\$/ $y/" out99_mcag.tsv
done
done

echo Gene: > out99_mcag2.tsv

for file in $(find /home/theaven/projects/niab/theaven/gene_pred/mildewreannotation/ -type f -name "final_genes_renamed.pep.fasta"); do
ID=$(echo $file | cut -d '/' -f9 | cut -d '.' -f1 )
#ID=$(echo $file | cut -d '/' -f9 | cut -d '.' -f1 | rev | cut -d '_' -f1,2 | rev )
echo $ID
sed -i "\$s/\$/ $ID/" out99_mcag2.tsv
done

for file in $(ls *_results_mcag_mcag); do
gene=$(echo $file | sed 's@_results_mcag_mcag@@g')
echo $gene >> out99_mcag2.tsv
for file2 in $(find /home/theaven/projects/niab/theaven/gene_pred/mildewreannotation/ -type f -name "final_genes_renamed.pep.fasta"); do
ID=$(echo $file2 | cut -d '/' -f9 | cut -d '.' -f1 | rev | cut -d '_' -f1,2 | rev )
x=$(grep -m 1 $ID $file | awk '{print $11}')
if [ -z "$x" ]; then
    y=0
else
    y=1
fi
sed -i "\$s/\$/ $y/" out99_mcag2.tsv
done
done

echo Gene: > out99_mcag_count.tsv

for file in $(find /home/theaven/projects/niab/theaven/gene_pred/mildewreannotation/ -type f -name "final_genes_renamed.pep.fasta"); do
ID=$(echo $file | cut -d '/' -f9 | cut -d '.' -f1 )
#ID=$(echo $file | cut -d '/' -f9 | cut -d '.' -f1 | rev | cut -d '_' -f1,2 | rev )
echo $ID
sed -i "\$s/\$/ $ID/" out99_mcag_count.tsv
done

for file in $(ls *_results_mcag_mcag); do
gene=$(echo $file | sed 's@_results_mcag_mcag@@g')
echo $gene >> out99_mcag_count.tsv
for file2 in $(find /home/theaven/projects/niab/theaven/gene_pred/mildewreannotation/ -type f -name "final_genes_renamed.pep.fasta"); do
ID=$(echo $file2 | cut -d '/' -f9 | cut -d '.' -f1 | rev | cut -d '_' -f1,2 | rev )
x=$(grep $ID $file | awk '{print $2}' | wc -l)
sed -i "\$s/\$/ $x/" out99_mcag_count.tsv
done
done

echo Gene: > out99_top_mcag.tsv

for file in $(find /home/theaven/projects/niab/theaven/gene_pred/mildewreannotation/ -type f -name "final_genes_renamed.pep.fasta"); do
ID=$(echo $file | cut -d '/' -f9 | cut -d '.' -f1 )
#ID=$(echo $file | cut -d '/' -f9 | cut -d '.' -f1 | rev | cut -d '_' -f1,2 | rev )
echo $ID
sed -i "\$s/\$/ $ID/" out99_top_mcag.tsv
done

for file in $(ls *_results_mcag_top_mcag); do
gene=$(echo $file | sed 's@_results_mcag_top_mcag@@g')
echo $gene >> out99_top_mcag.tsv
for file2 in $(find /home/theaven/projects/niab/theaven/gene_pred/mildewreannotation/ -type f -name "final_genes_renamed.pep.fasta"); do
ID=$(echo $file2 | cut -d '/' -f9 | cut -d '.' -f1 | rev | cut -d '_' -f1,2 | rev )
x=$(grep -m 1 $ID $file | awk '{print $11}')
if [ -z "$x" ]; then
    y=missing
else
    y=$x
fi
sed -i "\$s/\$/ $y/" out99_top_mcag.tsv
done
done
```
```bash
cd /home/theaven/scratch/manuscript4/mcag

#Append MCAG_ to the headers of MCAG proteins in order to make them easily identifiable in the blast results
for file in *.fsa; do
  [[ "$file" == GCA_* ]] && continue  
  
  awk '{ if ($0 ~ /^>/) { sub(/^>/, ">MCAG_"); print } else { print } }' "$file" > "${file}.tmp" && mv "${file}.tmp" "$file"
done

#Create a database with the downloaded MCAGs and yeast genes only for reciprocal blast
cat *protein.fsa >> proteins2.faa
cat GCA_000146045.2_R64_genomic.fsa >> proteins2.faa
srun -p medium --mem 16G -c 4 --pty bash
conda activate blast+
makeblastdb -in proteins2.faa -input_type fasta -dbtype prot  -title mcags_db2  -parse_seqids -out mcags_db2

#reciprocally BLAST the extracted MCAG hits from the mildew+ genomes against the MCAG and yeast proteome database
input_dir="/home/theaven/scratch/manuscript4/mcag/2/db"
output_dir="/home/theaven/scratch/manuscript4/mcag"
for file in "$input_dir"/*_2.faa; do
  filename=$(basename "$file")
  name="${filename%_2.faa}" 
  blastp -query "$file" \
         -db mcags_db2 \
         -out "${output_dir}/${name}_results_mcag2" \
         -evalue 1e-5 \
         -outfmt 6 \
         -num_threads 4
done

#Looking at the output files the cannonical MCAG sequences used must be slightly different from the proteins predicted de novo from the yeast genome as results returned for yeast commonly have 2 hits - the canonical MCAG and one other, presumably the predicted MCAG. 'Top' hits for MCAGs should therefore be considered to be either of these - either the canonical MCAG or the predicted MCAG.

#BLAST for the top hits only
for file in "$input_dir"/*_2.faa; do
  filename=$(basename "$file")
  name="${filename%_2.faa}" 
  blastp -query "$file" \
         -db mcags_db2 \
         -out "${output_dir}/${name}_results_mcag_top2" \
         -evalue 1e-5 \
         -outfmt 6 \
         -max_target_seqs 1 \
         -num_threads 4
done

#Extract reciprocal hits that are to the canonical MCAG or the predicted MCAG
grep 'MCAG_AAC03766' AAC03766_results_mcag_top2 > AAC03766_results_mcag_top2_mcag
grep 'MCAG_AAL85636' AAL85636_results_mcag_top2 > AAL85636_results_mcag_top2_mcag
grep 'MCAG_CAD28426' CAD28426_results_mcag_top2 > CAD28426_results_mcag_top2_mcag
grep 'MCAG_GDH3\|DAA06927.1' S288C_YAL062W_GDH3_results_mcag_top2 > S288C_YAL062W_GDH3_results_mcag_top2_mcag
grep 'MCAG_FUI1\|DAA07076.1' S288C_YBL042C_FUI1_results_mcag_top2 > S288C_YBL042C_FUI1_results_mcag_top2_mcag
grep 'MCAG_MNN2\|DAA07137.1' S288C_YBR015C_MNN2_results_mcag_top2 > S288C_YBR015C_MNN2_results_mcag_top2_mcag
grep 'MCAG_FUR4\|DAA07143.1' S288C_YBR021W_FUR4_results_mcag_top2 > S288C_YBR021W_FUR4_results_mcag_top2_mcag
grep 'MCAG_MET8\|DAA07329.1' S288C_YBR213W_MET8_results_mcag_top2 > S288C_YBR213W_MET8_results_mcag_top2_mcag
grep 'MCAG_SDS24\|DAA07330.1' S288C_YBR214W_SDS24_results_mcag_top2 > S288C_YBR214W_SDS24_results_mcag_top2_mcag
grep 'MCAG_MCX1\|DAA07342.1' S288C_YBR227C_MCX1_results_mcag_top2 > S288C_YBR227C_MCX1_results_mcag_top2_mcag
grep 'MCAG_PPS1\|DAA07392.1' S288C_YBR276C_PPS1_results_mcag_top2 > S288C_YBR276C_PPS1_results_mcag_top2_mcag
grep 'MCAG_APE3\|DAA07401.1' S288C_YBR286W_APE3_results_mcag_top2 > S288C_YBR286W_APE3_results_mcag_top2_mcag
grep 'MCAG_PHO89\|DAA07410.1' S288C_YBR296C_PHO89_results_mcag_top2 > S288C_YBR296C_PHO89_results_mcag_top2_mcag
grep 'MCAG_AAD3\|DAA07576.1' S288C_YCR107W_AAD3_results_mcag_top2 > S288C_YCR107W_AAD3_results_mcag_top2_mcag
grep 'MCAG_YDL144C\|DAA11714.1' S288C_YDL144C_YDL144C_results_mcag_top2 > S288C_YDL144C_YDL144C_results_mcag_top2_mcag

sort -k1,1 -k3,3nr S288C_YDL243C_AAD4_results_mcag2 | awk '!seen[$1]++' > S288C_YDL243C_AAD4_results_mcag_top2 #requires extra inspection as top hit by bitscore is not MCAG
grep 'MCAG_AAD4\|DAA11624.1' S288C_YDL243C_AAD4_results_mcag_top2 > S288C_YDL243C_AAD4_results_mcag_top2_mcag 

grep 'MCAG_RAD28\|DAA11874.1' S288C_YDR030C_RAD28_results_mcag_top2 > S288C_YDR030C_RAD28_results_mcag_top2_mcag
grep 'MCAG_MRX16\|DAA11978.1' S288C_YDR132C_MRX16_results_mcag_top2 > S288C_YDR132C_MRX16_results_mcag_top2_mcag
grep 'MCAG_AMD2\|DAA12082.1' S288C_YDR242W_AMD2_results_mcag_top2 > S288C_YDR242W_AMD2_results_mcag_top2_mcag
grep 'MCAG_RMT2\|DAA12300.1' S288C_YDR465C_RMT2_results_mcag_top2 > S288C_YDR465C_RMT2_results_mcag_top2_mcag
grep 'MCAG_YER152C\|DAA07813.1' S288C_YER152C_YER152C_results_mcag_top2 > S288C_YER152C_YER152C_results_mcag_top2_mcag
grep 'MCAG_PUG1\|DAA07848.1' S288C_YER185W_PUG1_results_mcag_top2 > S288C_YER185W_PUG1_results_mcag_top2_mcag

sort -k1,1 -k3,3nr S288C_YFL056C_AAD6_results_mcag2 | awk '!seen[$1]++' > S288C_YFL056C_AAD6_results_mcag_top2 #requires extra inspection as top hit by bitscore is not MCAG
grep 'MCAG_AAD6' S288C_YFL056C_AAD6_results_mcag_top2 > S288C_YFL056C_AAD6_results_mcag_top2_mcag 
#not the top hit for any yeast proteins

sort -k1,1 -k3,3nr S288C_YFL057C_AAD16_results_mcag2 | awk '!seen[$1]++' > S288C_YFL057C_AAD16_results_mcag_top2 #requires extra inspection as top hit by bitscore is not MCAG
grep 'MCAG_AAD16' S288C_YFL057C_AAD16_results_mcag_top2 > S288C_YFL057C_AAD16_results_mcag_top2_mcag 

grep 'MCAG_SDS23\|DAA08046.2' S288C_YGL056C_SDS23_results_mcag_top2 > S288C_YGL056C_SDS23_results_mcag_top2_mcag
grep 'MCAG_DSD1\|DAA07918.1' S288C_YGL196W_DSD1_results_mcag_top2 > S288C_YGL196W_DSD1_results_mcag_top2_mcag
grep 'MCAG_ARO8\|DAA07913.1' S288C_YGL202W_ARO8_results_mcag_top2 > S288C_YGL202W_ARO8_results_mcag_top2_mcag
grep 'MCAG_ADH4\|DAA07863.1' S288C_YGL256W_ADH4_results_mcag_top2 > S288C_YGL256W_ADH4_results_mcag_top2_mcag
grep 'MCAG_THI4\|DAA08235.1' S288C_YGR144W_THI4_results_mcag_top2 > S288C_YGR144W_THI4_results_mcag_top2_mcag
grep 'MCAG_GTO1\|DAA08245.1' S288C_YGR154C_GTO1_results_mcag_top2 > S288C_YGR154C_GTO1_results_mcag_top2_mcag
grep 'MCAG_RTA1\|DAA08307.1' S288C_YGR213C_RTA1_results_mcag_top2 > S288C_YGR213C_RTA1_results_mcag_top2_mcag
grep 'MCAG_YHB1\|DAA08325.1' S288C_YGR234W_YHB1_results_mcag_top2 > S288C_YGR234W_YHB1_results_mcag_top2_mcag
grep 'MCAG_DUR3\|DAA06669.1' S288C_YHL016C_DUR3_results_mcag_top2 > S288C_YHL016C_DUR3_results_mcag_top2_mcag
grep 'MCAG_DOG2\|DAA06735.1' S288C_YHR043C_DOG2_results_mcag_top2 > S288C_YHR043C_DOG2_results_mcag_top2_mcag
grep 'MCAG_DOG1\|DAA06736.1' S288C_YHR044C_DOG1_results_mcag_top2 > S288C_YHR044C_DOG1_results_mcag_top2_mcag
grep 'MCAG_ECM14\|DAA06824.1' S288C_YHR132C_ECM14_results_mcag_top2 > S288C_YHR132C_ECM14_results_mcag_top2_mcag
grep 'MCAG_ARO9\|DAA06830.1' S288C_YHR137W_ARO9_results_mcag_top2 > S288C_YHR137W_ARO9_results_mcag_top2_mcag
grep 'MCAG_FMO1\|DAA06869.1' S288C_YHR176W_FMO1_results_mcag_top2 > S288C_YHR176W_FMO1_results_mcag_top2_mcag
grep 'MCAG_YKE4\|DAA08522.1' S288C_YIL023C_YKE4_results_mcag_top2 > S288C_YIL023C_YKE4_results_mcag_top2_mcag
grep 'MCAG_GPP1\|DAA08494.1' S288C_YIL053W_GPP1_results_mcag_top2 > S288C_YIL053W_GPP1_results_mcag_top2_mcag
grep 'MCAG_YIL067C\|DAA08483.1' S288C_YIL067C_YIL067C_results_mcag_top2 > S288C_YIL067C_YIL067C_results_mcag_top2_mcag
grep 'MCAG_YIL108W\|DAA08445.1' S288C_YIL108W_YIL108W_results_mcag_top2 > S288C_YIL108W_YIL108W_results_mcag_top2_mcag
grep 'MCAG_AXL2\|DAA08412.1' S288C_YIL140W_AXL2_results_mcag_top2 > S288C_YIL140W_AXL2_results_mcag_top2_mcag
grep 'MCAG_SUC2\|DAA08390.1' S288C_YIL162W_SUC2_results_mcag_top2 > S288C_YIL162W_SUC2_results_mcag_top2_mcag
grep 'MCAG_DAL81\|DAA08570.1' S288C_YIR023W_DAL81_results_mcag_top2 > S288C_YIR023W_DAL81_results_mcag_top2_mcag
grep 'MCAG_DAL1\|DAA08574.1' S288C_YIR027C_DAL1_results_mcag_top2 > S288C_YIR027C_DAL1_results_mcag_top2_mcag
grep 'MCAG_DAL4\|DAA08575.1' S288C_YIR028W_DAL4_results_mcag_top2 > S288C_YIR028W_DAL4_results_mcag_top2_mcag
grep 'MCAG_DAL2\|DAA08576.1' S288C_YIR029W_DAL2_results_mcag_top2 > S288C_YIR029W_DAL2_results_mcag_top2_mcag
grep 'MCAG_TOK1\|DAA08707.1' S288C_YJL093C_TOK1_results_mcag_top2 > S288C_YJL093C_TOK1_results_mcag_top2_mcag
grep 'MCAG_SFH5\|DAA08656.1' S288C_YJL145W_SFH5_results_mcag_top2 > S288C_YJL145W_SFH5_results_mcag_top2_mcag
grep 'MCAG_MNN5\|DAA08620.2' S288C_YJL186W_MNN5_results_mcag_top2 > S288C_YJL186W_MNN5_results_mcag_top2_mcag
grep 'MCAG_MET3\|DAA08801.1' S288C_YJR010W_MET3_results_mcag_top2 > S288C_YJR010W_MET3_results_mcag_top2_mcag
grep 'MCAG_SOD1\|DAA08889.1' S288C_YJR104C_SOD1_results_mcag_top2 > S288C_YJR104C_SOD1_results_mcag_top2_mcag
grep 'MCAG_YJR124C\|DAA08909.1' S288C_YJR124C_YJR124C_results_mcag_top2 > S288C_YJR124C_YJR124C_results_mcag_top2_mcag
grep 'MCAG_AAD10\|DAA08940.1' S288C_YJR155W_AAD10_results_mcag_top2 > S288C_YJR155W_AAD10_results_mcag_top2_mcag
grep 'MCAG_MET14\|DAA09156.1' S288C_YKL001C_MET14_results_mcag_top2 > S288C_YKL001C_MET14_results_mcag_top2_mcag
grep 'MCAG_URA1\|DAA08953.1' S288C_YKL216W_URA1_results_mcag_top2 > S288C_YKL216W_URA1_results_mcag_top2_mcag
grep 'MCAG_MCH2\|DAA08948.2' S288C_YKL221W_MCH2_results_mcag_top2 > S288C_YKL221W_MCH2_results_mcag_top2_mcag
grep 'MCAG_MET1\|DAA09220.1' S288C_YKR069W_MET1_results_mcag_top2 > S288C_YKR069W_MET1_results_mcag_top2_mcag
grep 'MCAG_ECM4\|DAA09226.1' S288C_YKR076W_ECM4_results_mcag_top2 > S288C_YKR076W_ECM4_results_mcag_top2_mcag
grep 'MCAG_JLP1\|DAA09267.1' S288C_YLL057C_JLP1_results_mcag_top2 > S288C_YLL057C_JLP1_results_mcag_top2_mcag
grep 'MCAG_PDC1\|DAA09362.1' S288C_YLR044C_PDC1_results_mcag_top2 > S288C_YLR044C_PDC1_results_mcag_top2_mcag
grep 'MCAG_YLR046C\|DAA09364.1' S288C_YLR046C_YLR046C_results_mcag_top2 > S288C_YLR046C_YLR046C_results_mcag_top2_mcag
grep 'MCAG_FRE8\|DAA09365.1' S288C_YLR047C_FRE8_results_mcag_top2 > S288C_YLR047C_FRE8_results_mcag_top2_mcag
grep 'MCAG_YLR108C\|DAA09424.1' S288C_YLR108C_YLR108C_results_mcag_top2 > S288C_YLR108C_YLR108C_results_mcag_top2_mcag
grep 'MCAG_THI7\|DAA09551.1' S288C_YLR237W_THI7_results_mcag_top2 > S288C_YLR237W_THI7_results_mcag_top2_mcag
grep 'MCAG_YLR278C\|DAA09591.1' S288C_YLR278C_YLR278C_results_mcag_top2 > S288C_YLR278C_YLR278C_results_mcag_top2_mcag
grep 'MCAG_ECM38\|DAA09609.1' S288C_YLR299W_ECM38_results_mcag_top2 > S288C_YLR299W_ECM38_results_mcag_top2_mcag
grep 'MCAG_ALO1\|DAA09811.1' S288C_YML086C_ALO1_results_mcag_top2 > S288C_YML086C_ALO1_results_mcag_top2_mcag
grep 'MCAG_CCS1\|DAA09937.1' S288C_YMR038C_CCS1_results_mcag_top2 > S288C_YMR038C_CCS1_results_mcag_top2_mcag
grep 'MCAG_ARA2\|DAA09940.1' S288C_YMR041C_ARA2_results_mcag_top2 > S288C_YMR041C_ARA2_results_mcag_top2_mcag
grep 'MCAG_GTO3\|DAA10151.1' S288C_YMR251W_GTO3_results_mcag_top2 > S288C_YMR251W_GTO3_results_mcag_top2_mcag
grep 'MCAG_YME2\|DAA10203.1' S288C_YMR302C_YME2_results_mcag_top2 > S288C_YMR302C_YME2_results_mcag_top2_mcag
grep 'MCAG_URE2\|DAA10329.1' S288C_YNL229C_URE2_results_mcag_top2 > S288C_YNL229C_URE2_results_mcag_top2_mcag
grep 'MCAG_AAD14\|DAA10232.1' S288C_YNL331C_AAD14_results_mcag_top2 > S288C_YNL331C_AAD14_results_mcag_top2_mcag
grep 'MCAG_THI20\|DAA10728.1' S288C_YOL055C_THI20_results_mcag_top2 > S288C_YOL055C_THI20_results_mcag_top2_mcag
grep 'MCAG_ADH1\|DAA10699.1' S288C_YOL086C_ADH1_results_mcag_top2 > S288C_YOL086C_ADH1_results_mcag_top2_mcag
grep 'MCAG_BSC6\|DAA10648.1' S288C_YOL137W_BSC6_results_mcag_top2 > S288C_YOL137W_BSC6_results_mcag_top2_mcag

grep 'MCAG_YOL162W\|DAA10624.1' S288C_YOL162W_YOL162W_results_mcag2 > S288C_YOL162W_YOL162W_results_mcag_top2_mcag #requires extra inspection as top hit by  neither bitscore or pident is not MCAG
#there is no good hit represneting a mcag candidate from yeast

sort -k1,1 -k3,3nr S288C_YOL165C_AAD15_results_mcag2 | awk '!seen[$1]++' > S288C_YOL165C_AAD15_results_mcag_top2 #requires extra inspection as top hit by bitscore is not MCAG
grep 'MCAG_AAD15\|DAA10620.1' S288C_YOL165C_AAD15_results_mcag_top2 > S288C_YOL165C_AAD15_results_mcag_top2_mcag 

grep 'MCAG_NRT1\|DAA10850.1' S288C_YOR071C_NRT1_results_mcag_top2 > S288C_YOR071C_NRT1_results_mcag_top2_mcag
grep 'MCAG_THI72\|DAA10964.1' S288C_YOR192C_THI72_results_mcag_top2 > S288C_YOR192C_THI72_results_mcag_top2_mcag
grep 'MCAG_HEM4\|DAA11043.1' S288C_YOR278W_HEM4_results_mcag_top2 > S288C_YOR278W_HEM4_results_mcag_top2_mcag
grep 'MCAG_GDH1\|DAA11134.1' S288C_YOR375C_GDH1_results_mcag_top2 > S288C_YOR375C_GDH1_results_mcag_top2_mcag
grep 'MCAG_FDH1\|DAA11147.1' S288C_YOR388C_FDH1_results_mcag_top2 > S288C_YOR388C_FDH1_results_mcag_top2_mcag
grep 'MCAG_YPL088W\|DAA11345.1' S288C_YPL088W_YPL088W_results_mcag_top2 > S288C_YPL088W_YPL088W_results_mcag_top2_mcag
grep 'MCAG_PNG1\|DAA11337.1' S288C_YPL096W_PNG1_results_mcag_top2 > S288C_YPL096W_PNG1_results_mcag_top2_mcag
grep 'MCAG_FMP30\|DAA11330.1' S288C_YPL103C_FMP30_results_mcag_top2 > S288C_YPL103C_FMP30_results_mcag_top2_mcag
grep 'MCAG_THI6\|DAA11222.1' S288C_YPL214C_THI6_results_mcag_top2 > S288C_YPL214C_THI6_results_mcag_top2_mcag
grep 'MCAG_THI21\|DAA11177.1' S288C_YPL258C_THI21_results_mcag_top2 > S288C_YPL258C_THI21_results_mcag_top2_mcag

sort -k1,1 -k3,3nr S288C_YPL277C_YPL277C_results_mcag2 | awk '!seen[$1]++' > S288C_YPL277C_YPL277C_results_mcag_top2 #requires extra inspection as top hit by bitscore is not MCAG
grep 'MCAG_YPL277C\|DAA11161.1' S288C_YPL277C_YPL277C_results_mcag_top2 > S288C_YPL277C_YPL277C_results_mcag_top2_mcag 
#The % identities outside yeast are very low.

grep 'MCAG_SDD4\|DAA11448.1' S288C_YPR022C_SDD4_results_mcag_top2 > S288C_YPR022C_SDD4_results_mcag_top2_mcag
grep 'MCAG_THI22\|DAA11536.2' S288C_YPR121W_THI22_results_mcag_top2 > S288C_YPR121W_THI22_results_mcag_top2_mcag
grep 'MCAG_YPR127W\|DAA11540.1' S288C_YPR127W_YPR127W_results_mcag_top2 > S288C_YPR127W_YPR127W_results_mcag_top2_mcag
grep 'MCAG_MET16\|DAA11584.1' S288C_YPR167C_MET16_results_mcag_top2 > S288C_YPR167C_MET16_results_mcag_top2_mcag
grep 'MCAG_ARR3\|DAA11615.1' S288C_YPR201W_ARR3_results_mcag_top2 > S288C_YPR201W_ARR3_results_mcag_top2_mcag
grep 'MCAG_XP_1548821' XP_001548821_results_mcag_top2 > XP_001548821_results_mcag_top2_mcag
grep 'MCAG_XP_001560697' XP_001560697_results_mcag_top2 > XP_001560697_results_mcag_top2_mcag
grep 'MCAG_XP_752655' XP_752655_results_mcag_top2 > XP_752655_results_mcag_top2_mcag

#Collect results into one file:
echo Gene: > out99_mcag2-2.tsv

for file in $(find /home/theaven/projects/niab/theaven/gene_pred/mildewreannotation/ -type f -name "final_genes_renamed.pep.fasta"); do
ID=$(echo $file | cut -d '/' -f9 | cut -d '.' -f1 )
#ID=$(echo $file | cut -d '/' -f9 | cut -d '.' -f1 | rev | cut -d '_' -f1,2 | rev )
echo $ID
sed -i "\$s/\$/ $ID/" out99_mcag2-2.tsv
done

for file in $(ls *_results_mcag_top2_mcag); do
gene=$(echo $file | sed 's@_results_mcag_top2_mcag@@g')
echo $gene >> out99_mcag2-2.tsv
for file2 in $(find /home/theaven/projects/niab/theaven/gene_pred/mildewreannotation/ -type f -name "final_genes_renamed.pep.fasta"); do
ID=$(echo $file2 | cut -d '/' -f9 | cut -d '.' -f1 | rev | cut -d '_' -f1,2 | rev )
x=$(grep -m 1 $ID $file | awk '{print $11}')
if [ -z "$x" ]; then
    y=0
else
    y=1
fi
sed -i "\$s/\$/ $y/" out99_mcag2-2.tsv
done
done

```



```bash
singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/python3.sif python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file temp_id.txt --input  /home/theaven/projects/niab/theaven/genomes/Pod_leu_OGBp112020_clean_renamed.fasta --output  ~/scratch/manuscript4/test_genome.fasta

srun -p medium --mem 16G -c 4 --pty bash
conda activate earlgreyte

# first, change directory to the famdb library location
cd /mnt/apps/users/theaven/conda/envs/earlgreyte/share/RepeatMasker/Libraries/famdb/

# download the partitions you require from Dfam 3.9. In the below, change the numbers or range inside the square brackets to choose your subsets.
# e.g. to download partitions 0 to 10: [0-10]; or to download partitions 3,5, and 7: [3,5,7]; [0-16] is ALL PARTITIONS
curl -O https://dfam.org/releases/current/families/FamDB/dfam39_full.16.h5.gz
curl -O https://dfam.org/releases/current/families/FamDB/dfam39_full.0.h5.gz

# decompress Dfam 3.9 paritions
gunzip *.gz

# move up to RepeatMasker main directory
cd /mnt/apps/users/theaven/conda/envs/earlgreyte/share/RepeatMasker/

# save the min_init partition as a backup, just in case!
mv /mnt/apps/users/theaven/conda/envs/earlgreyte/share/RepeatMasker/Libraries/famdb/min_init.0.h5 /mnt/apps/users/theaven/conda/envs/earlgreyte/share/RepeatMasker/Libraries/famdb/min_init.0.h5.bak

# Rerun RepeatMasker configuration
perl ./configure -libdir /mnt/apps/users/theaven/conda/envs/earlgreyte/share/RepeatMasker/Libraries/ -trf_prgm /mnt/apps/users/theaven/conda/envs/earlgreyte/bin/trf -rmblast_dir /mnt/apps/users/theaven/conda/envs/earlgreyte/bin -hmmer_dir /mnt/apps/users/theaven/conda/envs/earlgreyte/bin -abblast_dir /mnt/apps/users/theaven/conda/envs/earlgreyte/bin  


earlGreyLibConstruct -g /home/theaven/projects/niab/theaven/genomes/Ust_may_GCA_000328475.2_clean_renamed.fasta -s test -o . -t 4 -r leotiomycetes

earlGreyLibConstruct -g /home/theaven/projects/niab/theaven/genomes/Ery_pis_GCA_000208805.1_clean_renamed.fasta -s test -o test -t 4

earlGrey -g ~/scratch/manuscript4/test_genome.fasta -s test -o . -t 4 -d yes -r leotiomycetes


conda activate transposonpsi
/mnt/apps/users/theaven/conda/envs/transposonpsi/share/transposonPSI/transposonPSI.pl ~/scratch/manuscript4/test_genome.fasta nuc

formatdb.log
test_genome.fasta.TPSI.allHits.chains.gff3
test_genome.fasta.TPSI.allHits.chains.bestPerLocus.gff3
test_genome.fasta.TPSI.allHits.chains.bestPerLocus
test_genome.fasta.TPSI.allHits.chains
test_genome.fasta.TPSI.allHits

bedtools getfasta \
  -fi test_genome.fasta \
  -bed test_genome.fasta.TPSI.allHits.chains.bestPerLocus.gff3 \
  -s -name+ \
  -fo TPSI_bestPerLocus.fa

  vsearch --cluster_fast TPSI_bestPerLocus.fa \
        --id 0.8 \
        --centroids TPSI_clusters.fa \
        --consout TPSI_consensus.fasta










screen -S reannotate
cd /home/theaven/scratch/manuscript4
conda activate earlgreyte
for Assembly in $(ls /home/theaven/projects/niab/theaven/genomes/*.fasta ~/scratch/manuscript4/genomes/*/*_fix.fasta); do
ID=$(basename $Assembly | sed 's@_clean_renamed.fasta@@g' | sed 's@_clean_renamed_fix.fasta@@g')
Jobs=$(squeue -u theaven| grep 'Earlgrey'  | wc -l)
echo $ID
while [ $Jobs -gt 2 ]; do
    sleep 3600s
    printf "."
    Jobs=$(squeue -u theaven| grep 'Earlgrey'  | wc -l)
done

OutDir=/home/theaven/scratch/manuscript4/reannotation/$ID/earlgrey_Lib
  jobid=$(sbatch --parsable ~/git_repos/Wrappers/gruffalo/run_earlgrey_Lib.sh "$Assembly" "$ID" "$OutDir")
  printf "%s\t%s\tearlgrey_Lib\t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> slurm_log.tsv
done 
conda deactivate

conda activate transposonpsi
for Assembly in $(ls /home/theaven/projects/niab/theaven/genomes/*.fasta); do
ID=$(basename $Assembly | sed 's@_clean_renamed.fasta@@g')
Jobs=$(squeue -u theaven| grep 'TPSI'  | wc -l)
echo $ID
while [ $Jobs -gt 2 ]; do
    sleep 3600s
    printf "."
    Jobs=$(squeue -u theaven| grep 'TPSI'  | wc -l)
done

OutDir=/home/theaven/scratch/manuscript4/reannotation/$ID/TPSI
  jobid=$(sbatch --parsable ~/git_repos/Wrappers/gruffalo/run_transposonpsi.sh "$Assembly" "$ID" "$OutDir")
  printf "%s\t%s\tTPSI\t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> slurm_log.tsv
done
conda deactivate



screen -S reannotate
cd /home/theaven/scratch/manuscript4
# ----- Earlgrey -----
conda activate earlgreyte
for Assembly in $(ls /home/theaven/projects/niab/theaven/genomes/*.fasta ~/scratch/manuscript4/genomes/*/*_fix.fasta); do
  [[ -e "$Assembly" ]] || { echo "No FASTA files found"; break; }
  ID="$(basename "$Assembly" | sed 's@_clean_renamed\.fasta@@' | sed 's@_clean_renamed_fix\.fasta@@')"
  echo "$ID"

  Jobs=$(squeue -h -u theaven -n EarlgreyTE | wc -l)
  while [ "$Jobs" -gt 5 ]; do
    sleep 3600s        
    printf "."
    Jobs=$(squeue -h -u theaven -n EarlgreyTE | wc -l)
  done
  echo

  OutDir="/home/theaven/scratch/manuscript4/reannotation/$ID/earlgrey_Lib"
  if [ ! -f "$OutDir/EG_lib/${ID}_EarlGrey/${ID}_summaryFiles/${ID}-families.fa.strained" ]; then
  jobid=$(sbatch --parsable ~/git_repos/Wrappers/gruffalo/run_earlgrey_Lib.sh "$Assembly" "$ID" "$OutDir")
  printf "%s\t%s\tearlgrey_Lib\t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> slurm_log.tsv
  else
  echo "Already run for ${ID}"
  fi
done
conda deactivate

for Assembly in /home/theaven/projects/niab/theaven/genomes/*.fasta; do
  ID="$(basename "$Assembly" | sed 's@_clean_renamed\.fasta@@')"
  EG="/home/theaven/scratch/manuscript4/reannotation/${ID}/earlgrey_Lib/EG_lib/${ID}_EarlGrey/${ID}_summaryFiles/${ID}-families.fa.strained"
  sed -i 's@^>@\>EG_@' $EG
done

# ----- TPSI -----
conda activate transposonpsi
for Assembly in $(ls /home/theaven/projects/niab/theaven/genomes/*.fasta ~/scratch/manuscript4/genomes/*/*_fix.fasta); do
  [[ -e "$Assembly" ]] || { echo "No FASTA files found"; break; }
  ID="$(basename "$Assembly" | sed 's@_clean_renamed\.fasta@@' | sed 's@_clean_renamed_fix\.fasta@@')"
  echo "$ID"

  Jobs=$(squeue -h -u theaven -n TPSI | wc -l)
  while [ "$Jobs" -gt 5 ]; do
    sleep 3600s
    printf "."
    Jobs=$(squeue -h -u theaven -n TPSI | wc -l)
  done
  echo

  OutDir="/home/theaven/scratch/manuscript4/reannotation/$ID/TPSI"
  if [ ! -f "$OutDir/${ID}_clean_renamed.fasta.TPSI.allHits.chains.bestPerLocus.gff3" ]; then
  jobid=$(sbatch --parsable ~/git_repos/Wrappers/gruffalo/run_transposonpsi.sh "$Assembly" "$ID" "$OutDir")
  printf "%s\t%s\tTPSI\t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> slurm_log.tsv
  else
  echo "Already run for ${ID}"
  fi
done
conda deactivate
echo Done

#Get consensus TE sequences from TransposonPSI output
#Earl Grey follows the “80–80–80 rule” for TE family clustering:≥ 80% sequence identity; across ≥ 80% of the aligned length; for sequences that are ≥ 80 bp long

srun -p long -J bedtools --mem-per-cpu 4G --cpus-per-task 8 --pty bash
conda activate bedtools
for Assembly in $(ls /home/theaven/projects/niab/theaven/genomes/Phy_mor_GCA_019455665.1_clean_renamed.fasta); do
  [[ -e "$Assembly" ]] || { echo "No FASTA files found"; break; }
  ID="$(basename "$Assembly" | sed 's@_clean_renamed\.fasta@@')"
  echo "$ID"
  TPSI_DIR="/home/theaven/scratch/manuscript4/reannotation/${ID}/TPSI"
bedtools getfasta \
  -fi "$Assembly" \
  -bed "$TPSI_DIR"/"$ID"_clean_renamed.fasta.TPSI.allHits.chains.bestPerLocus.gff3 \
  -fo "$TPSI_DIR"/"$ID"_TPSI_hits.fa
cd-hit-est -i "$TPSI_DIR"/"$ID"_TPSI_hits.fa -o "$TPSI_DIR"/"$ID"_TPSI_clusters.fa -c 0.8 -n 5 -d 0 -T 8 -M 30000

mkdir "$TPSI_DIR"/cluster_members
awk -v out="$TPSI_DIR/cluster_members" '
  BEGIN { c = -1 }
  /^>Cluster[[:space:]][0-9]+/ { c++; next }
  {
    s = $0
    sub(/^.*>/, "", s)          # drop everything up to and including the first ">"
    sub(/\.\.\..*/, "", s)      # drop from the first "..." to end
    gsub(/[ \t\r]+$/, "", s)    # trim trailing whitespace just in case
    if (length(s) > 0) {
      printf("%s\n", s) >> sprintf("%s/cluster_%05d.ids", out, c)
    }
  }
' "$TPSI_DIR/${ID}_TPSI_clusters.fa.clstr"

  for ids in "$TPSI_DIR"/cluster_members/*.ids; do
    cbase="$(basename "$ids" .ids)"
    seqkit grep -f "$ids" "$TPSI_DIR/${ID}_TPSI_hits.fa" \
      > "$TPSI_DIR/cluster_members/${cbase}.fa" || true
  done

mkdir "$TPSI_DIR/cluster_consensus"
find "$TPSI_DIR/cluster_members" -maxdepth 1 -name '*.fa' -print0 |
while IFS= read -r -d '' fa; do
  [[ -s "$fa" ]] || continue
  echo "Running for $fa"
  nseq=$(grep -c '^>' "$fa" || true)
  cbase=$(basename "$fa" .fa)
  aln="$TPSI_DIR/cluster_members/${cbase}.aln.fa"
  cons="$TPSI_DIR/cluster_consensus/${ID}_${cbase}_cons.fa"

  if (( nseq >= 2 )); then
    echo "multiple sequences"
    mafft --auto --thread 8 "$fa" > "$aln"
    consambig -sequence "$aln" -outseq "$cons"
    echo "output = $cons"
  else
    echo "single sequence"
    cp "$fa" "$cons"
    sed -i "1s/^>.*/>${ID}_${cbase}_cons/" "$cons"
    echo "output = $cons"
  fi
done

for cons in "$TPSI_DIR"/cluster_consensus/*.fa; do
  newhdr=">TPSI_"$(basename "$cons" .fa)
  awk -v h="$newhdr" 'NR==1{print h; next} {print}' "$cons" > "$cons.tmp" && mv "$cons.tmp" "$cons"
done

cat "$TPSI_DIR"/cluster_consensus/*.fa > "$TPSI_DIR/TPSI_consensus.fa"
echo "Wrote: $TPSI_DIR/TPSI_consensus.fa"

done

# Combine
for Assembly in /home/theaven/projects/niab/theaven/genomes/*.fasta ~/scratch/manuscript4/genomes/*/*_fix.fasta; do
  ID="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@' | sed 's@_clean_renamed\.fasta@@')"
TPSI=/home/theaven/scratch/manuscript4/reannotation/${ID}/TPSI/TPSI_consensus.fa
EG="/home/theaven/scratch/manuscript4/reannotation/${ID}/earlgrey_Lib/EG_lib/${ID}_EarlGrey/${ID}_summaryFiles/${ID}-families.fa.strained"
ls -lh $TPSI $EG
cat $TPSI $EG > /home/theaven/scratch/manuscript4/reannotation/${ID}/combined_TE_library.fasta
done

# Annotate TEs
conda activate earlgreyte
for Assembly in /home/theaven/projects/niab/theaven/genomes/Phy_mor_GCA_019455665.1_clean_renamed.fasta ~/scratch/manuscript4/genomes/*/*_fix.fasta; do
  ID="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@' | sed 's@_clean_renamed\.fasta@@')"
  Lib="/home/theaven/scratch/manuscript4/reannotation/${ID}/combined_TE_library.fasta"
  OutDir="/home/theaven/scratch/manuscript4/reannotation/$ID/earlgrey_annotation"
  jobid=$(sbatch --parsable ~/git_repos/Wrappers/gruffalo/run_earlgrey_annotation.sh "$Assembly" "$ID" "$OutDir" "$Lib")
  printf "%s\t%s\tearlgrey_annotation\t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> slurm_log.tsv
done

# Mask the genomes
srun -p long -J bedtools  --mem-per-cpu 16G --cpus-per-task 1 --pty bash
for Assembly in /home/theaven/projects/niab/theaven/genomes/Phy_mor_GCA_019455665.1_clean_renamed.fasta ~/scratch/manuscript4/genomes/*/*_fix.fasta; do
  ID="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@' | sed 's@_clean_renamed\.fasta@@')"
  Repeats="/home/theaven/scratch/manuscript4/reannotation/$ID/earlgrey_annotation/${ID}_EarlGrey/${ID}_summaryFiles/${ID}.filteredRepeats.gff"
  OutDir="/home/theaven/scratch/manuscript4/reannotation/$ID/earlgrey_annotation/${ID}_EarlGrey"
bedtools maskfasta \
  -fi "$Assembly" \
  -bed "$Repeats" \
  -fo ${OutDir}/${ID}.softmasked.fa \
  -soft

bedtools maskfasta \
  -fi "$Assembly" \
  -bed "$Repeats" \
  -fo ${OutDir}/${ID}.hardmasked.fa 
done

for Assembly in $(ls /home/theaven/projects/niab/theaven/genomes/*_renamed.fasta ~/scratch/manuscript4/genomes/*/*_fix.fasta); do
  ID="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@' | sed 's@_clean_renamed\.fasta@@')"
  contigs="$(grep '>' "$Assembly" | wc -l)"
  echo "$ID $contigs"
done

# Gene prediction - abinitio
for Assembly in $(ls /home/theaven/projects/niab/theaven/genomes/*.fasta); do
  ID="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@' | sed 's@_clean_renamed\.fasta@@')"

    Jobs=$(squeue -h -u theaven -n braker | wc -l)
  while [ "$Jobs" -gt 10 ]; do
    sleep 900s
    printf "."
    Jobs=$(squeue -h -u theaven -n braker | wc -l)
  done
  echo

  Masked_Genome="/home/theaven/scratch/manuscript4/reannotation/$ID/earlgrey_annotation/${ID}_EarlGrey/${ID}.softmasked.fa"
  RNA_alignment=NA
  Protein_database=NA
  GeneMark=NA
  Min_contig=50000
  Braker_run=$(echo $ID)_ab_initio
  OutDir="/home/theaven/scratch/manuscript4/reannotation/$ID/braker_abinitio"
  mkdir "$OutDir"
  if [ ! -f "$OutDir/braker.aa" ]; then
  jobid=$(sbatch --parsable ~/git_repos/Wrappers/gruffalo/run_braker3.sh "$Masked_Genome" "$Braker_run" "$OutDir" "$RNA_alignment" "$Protein_database" "$Min_contig" "$GeneMark")
  printf "%s\t%s\tbraker_abinitio_annotation\t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> slurm_log.tsv
  else
  echo "Found $OutDir/braker.aa"
  fi
done

#Gene prediction - abinitio - fragmented genomes
#GeneMark-ES does self-training: it extracts long, “clean” fragments from the input and uses them to build parameters. Genemark sees masked regions and long stretches of NNNs as gaps. Fragments of greater than the Min_contig are required for genemark training. When nothing passes its preprocessing filters (most commonly: too much is masked/lowercase or N) there’s nothing to make a training set from. Can either use an unmasked genome or smaller min_contig length during the genemark step to fix this - 50000 is the default.
#NOTE: with very many tiny scaffolds, the parallel “chunking” AUGUSTUS does can thrash—too many tiny jobs, file handles, and scheduler overhead. The safe fix is to run in linear mode (single thread).
for Assembly in $(ls /home/theaven/projects/niab/theaven/genomes/*.fasta); do
  ID="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@' | sed 's@_clean_renamed\.fasta@@')"
  Masked_Genome="/home/theaven/scratch/manuscript4/reannotation/$ID/earlgrey_annotation/${ID}_EarlGrey/${ID}.softmasked.fa"
  OutDir="/home/theaven/scratch/manuscript4/reannotation/$ID/GMES-soft"
  mkdir "$OutDir"
  if [ ! -f "$OutDir/genemark.gtf" ]; then
  jobid=$(sbatch --parsable ~/git_repos/Wrappers/gruffalo/run_genemarkES-fragmented.sh "$Masked_Genome" "$OutDir")
  printf "%s\t%s\tGMES_soft_frag\t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> slurm_log.tsv
  else
  echo "Found $OutDir/genemark.gtf"
  fi
done

cp /home/theaven/projects/niab/theaven/genomes/Pod_aph_SCOTT2020-2_clean_renamed_fix.fasta /home/theaven/projects/niab/theaven/genomes/Pod_aph_SCOTT2020-2_clean_renamed_fix.fasta_backup
awk 'BEGIN{FS="[ \t]"}
     /^>/{print $1; next}
     {print}' /home/theaven/projects/niab/theaven/genomes/Pod_aph_SCOTT2020-2_clean_renamed_fix.fasta > masked.simple.fa && mv masked.simple.fa /home/theaven/projects/niab/theaven/genomes/Pod_aph_SCOTT2020-2_clean_renamed_fix.fasta

cp /home/theaven/scratch/manuscript4/reannotation/Pod_aph_SCOTT2020-2/earlgrey_annotation/Pod_aph_SCOTT2020-2_EarlGrey/Pod_aph_SCOTT2020-2.softmasked.fa /home/theaven/scratch/manuscript4/reannotation/Pod_aph_SCOTT2020-2/earlgrey_annotation/Pod_aph_SCOTT2020-2_EarlGrey/Pod_aph_SCOTT2020-2.softmasked.fa_backup
awk 'BEGIN{FS="[ \t]"} 
     /^>/{h=$1; sub(/^>/,"",h); sub(/^lcl\|/,"",h); print ">"h; next}
     {print}' /home/theaven/scratch/manuscript4/reannotation/Pod_aph_SCOTT2020-2/earlgrey_annotation/Pod_aph_SCOTT2020-2_EarlGrey/Pod_aph_SCOTT2020-2.softmasked.fa > masked.simple.fa && mv masked.simple.fa /home/theaven/scratch/manuscript4/reannotation/Pod_aph_SCOTT2020-2/earlgrey_annotation/Pod_aph_SCOTT2020-2_EarlGrey/Pod_aph_SCOTT2020-2.softmasked.fa

cp /home/theaven/scratch/manuscript4/reannotation/Pod_aph_SCOTT2020-2/GMES-soft/genemark.gtf /home/theaven/scratch/manuscript4/reannotation/Pod_aph_SCOTT2020-2/GMES-soft/genemark.gtf_backup
awk 'BEGIN{OFS="\t"} {sub(/^lcl\|/,"",$1); print}' /home/theaven/scratch/manuscript4/reannotation/Pod_aph_SCOTT2020-2/GMES-soft/genemark.gtf > genemark.nolcl.gtf && mv genemark.nolcl.gtf /home/theaven/scratch/manuscript4/reannotation/Pod_aph_SCOTT2020-2/GMES-soft/genemark.gtf

sed -E 's/^lcl\|//' /home/theaven/scratch/manuscript4/reannotation/Pod_aph_SCOTT2020-2/GMES-soft/genemark.gtf_backup > genemark.nolcl.gtf && mv genemark.nolcl.gtf /home/theaven/scratch/manuscript4/reannotation/Pod_aph_SCOTT2020-2/GMES-soft/genemark.gtf

for Assembly in $(ls /home/theaven/projects/niab/theaven/genomes/Pod_aph_SCOTT2020-2*.fasta); do
  ID="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@' | sed 's@_clean_renamed\.fasta@@')"

    Jobs=$(squeue -h -u theaven -n braker | wc -l)
  while [ "$Jobs" -gt 10 ]; do
    sleep 900s
    printf "."
    Jobs=$(squeue -h -u theaven -n braker | wc -l)
  done
  echo

  Masked_Genome="/home/theaven/scratch/manuscript4/reannotation/$ID/earlgrey_annotation/${ID}_EarlGrey/${ID}.softmasked.fa"
  RNA_alignment=NA
  Protein_database=NA
  Braker_run=$(echo $ID)_ab_initio_soft_10000
  GeneMark="/home/theaven/scratch/manuscript4/reannotation/$ID/GMES-soft/genemark.gtf"
  OutDir="/home/theaven/scratch/manuscript4/reannotation/$ID/braker_abinitio_10000"
  #mkdir "$OutDir/rsync"
  if [ ! -f "$OutDir/braker.aa" ]; then
  echo "Not Found................. $OutDir/braker.aa"    
  jobid=$(sbatch --parsable ~/git_repos/Wrappers/gruffalo/run_braker3-aug.sh "$Masked_Genome" "$RNA_alignment" "$Protein_database" "$Braker_run" "$OutDir" "$GeneMark")
  printf "%s\t%s\tbraker_abinitio_annotation\t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> slurm_log.tsv
  else
  echo "Found $OutDir/braker.aa"
  fi
done

#Gene prediction - abinitio - hardmask for GMES step
for Assembly in $(ls /home/theaven/projects/niab/theaven/genomes/*.fasta); do
  ID="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@' | sed 's@_clean_renamed\.fasta@@')"
  Masked_Genome="/home/theaven/scratch/manuscript4/reannotation/$ID/earlgrey_annotation/${ID}_EarlGrey/${ID}.hardmasked.fa"
  OutDir="/home/theaven/scratch/manuscript4/reannotation/$ID/GMES-hard-frag"
  mkdir "$OutDir"
  if [ ! -f "$OutDir/genemark.gtf" ]; then
  jobid=$(sbatch --parsable ~/git_repos/Wrappers/gruffalo/run_genemarkES_hard.sh "$Masked_Genome" "$OutDir")
  printf "%s\t%s\tGMES_hard\t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> slurm_log.tsv
  echo "missing for $ID"
  else
  echo "Found $OutDir/genemark.gtf"
  fi
done

cp /home/theaven/scratch/manuscript4/reannotation/Pod_aph_SCOTT2020-2/GMES-hard-frag/genemark.gtf /home/theaven/scratch/manuscript4/reannotation/Pod_aph_SCOTT2020-2/GMES-hard-frag/genemark.gtf_backup
sed -E 's/^lcl\|//' /home/theaven/scratch/manuscript4/reannotation/Pod_aph_SCOTT2020-2/GMES-hard-frag/genemark.gtf > genemark.nolcl.gtf && mv genemark.nolcl.gtf /home/theaven/scratch/manuscript4/reannotation/Pod_aph_SCOTT2020-2/GMES-hard-frag/genemark.gtf

for Assembly in $(ls /home/theaven/projects/niab/theaven/genomes/*.fasta); do
  ID="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@' | sed 's@_clean_renamed\.fasta@@')"
  Masked_Genome="/home/theaven/scratch/manuscript4/reannotation/$ID/earlgrey_annotation/${ID}_EarlGrey/${ID}.softmasked.fa"
  RNA_alignment=NA
  Protein_database=NA
  Braker_run=$(echo $ID)_ab_initio_hard_10000
  GeneMark="/home/theaven/scratch/manuscript4/reannotation/$ID/GMES-hard-frag/genemark.gtf"
  OutDir="/home/theaven/scratch/manuscript4/reannotation/$ID/braker_abinitio_hard_10000"
  mkdir -p "$OutDir"
  if [ ! -f "$OutDir/braker.aa" ]; then
  rm -rf $OutDir/*
  rm -rf /mnt/shared/projects/niab/theaven/augustus_config/species/${Braker_run}
  jobid=$(sbatch --parsable ~/git_repos/Wrappers/gruffalo/run_braker3-aug.sh "$Masked_Genome" "$RNA_alignment" "$Protein_database" "$Braker_run" "$OutDir" "$GeneMark")
  printf "%s\t%s\tbraker_abinitio_annotation_hard\t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> slurm_log.tsv
  else
  echo "Found $OutDir/braker.aa"
  fi
done

cp /home/theaven/scratch/manuscript4/reannotation/Pod_aph_SCOTT2020-2/earlgrey_annotation/Pod_aph_SCOTT2020-2_EarlGrey/Pod_aph_SCOTT2020-2.hardmasked.fa /home/theaven/scratch/manuscript4/reannotation/Pod_aph_SCOTT2020-2/earlgrey_annotation/Pod_aph_SCOTT2020-2_EarlGrey/Pod_aph_SCOTT2020-2.hardmasked.fa_backup
awk 'BEGIN{FS="[ \t]"} 
     /^>/{h=$1; sub(/^>/,"",h); sub(/^lcl\|/,"",h); print ">"h; next}
     {print}' /home/theaven/scratch/manuscript4/reannotation/Pod_aph_SCOTT2020-2/earlgrey_annotation/Pod_aph_SCOTT2020-2_EarlGrey/Pod_aph_SCOTT2020-2.hardmasked.fa > masked.simple.fa && mv masked.simple.fa /home/theaven/scratch/manuscript4/reannotation/Pod_aph_SCOTT2020-2/earlgrey_annotation/Pod_aph_SCOTT2020-2_EarlGrey/Pod_aph_SCOTT2020-2.hardmasked.fa

#Gene prediction - abinitio - hardmask for both steps
for Assembly in $(ls /home/theaven/projects/niab/theaven/genomes/*.fasta); do
  ID="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@' | sed 's@_clean_renamed\.fasta@@')"

      Jobs=$(squeue -h -u theaven -n braker | wc -l)
  while [ "$Jobs" -gt 30 ]; do
    sleep 900s
    printf "."
    Jobs=$(squeue -h -u theaven -n braker | wc -l)
  done
  echo

  Masked_Genome="/home/theaven/scratch/manuscript4/reannotation/$ID/earlgrey_annotation/${ID}_EarlGrey/${ID}.hardmasked.fa"
  RNA_alignment=NA
  Protein_database=NA
  Braker_run=$(echo $ID)_ab_initio_hard2_10000
  GeneMark="/home/theaven/scratch/manuscript4/reannotation/$ID/GMES-hard-frag/genemark.gtf"
  OutDir="/home/theaven/scratch/manuscript4/reannotation/$ID/braker_abinitio_hard2_10000"
  mkdir -p "$OutDir"
  if [ ! -f "$OutDir/braker.aa" ]; then
  rm -rf $OutDir/*
  rm -rf /mnt/shared/projects/niab/theaven/augustus_config/species/${Braker_run}
  jobid=$(sbatch --parsable ~/git_repos/Wrappers/gruffalo/run_braker3-aug.sh "$Masked_Genome" "$RNA_alignment" "$Protein_database" "$Braker_run" "$OutDir" "$GeneMark")
  printf "%s\t%s\tbraker_abinitio_annotation_hard2\t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> slurm_log.tsv
  else
  echo "Found $OutDir/braker.aa"
  fi
done

ls /home/theaven/scratch/manuscript4/reannotation/*/braker_abinitio/braker.aa | wc -l
ls /home/theaven/scratch/manuscript4/reannotation/*/braker_abinitio_10000/braker.aa | wc -l
ls /home/theaven/scratch/manuscript4/reannotation/*/braker_abinitio_hard_10000/braker.aa | wc -l
ls /home/theaven/scratch/manuscript4/reannotation/*/braker_abinitio_hard2_10000/braker.aa | wc -l
ls /home/theaven/scratch/manuscript4/reannotation/*/braker_prot/braker.aa | wc -l

Blu_gra_GCA_000417865.1 3814925 3877442 3879389
Blu_gra_GCA_000441875.1 
Blu_gra_GCA_900519115.1 

Oid_hev_GCA_003957845.1 3879174
Fus_oxy_GCF_013085055.1 3869890
Gla_loz_GCA_000409485.1 
Gol_cic_GCA_003611195.1 

theaven@gruffalo:~/scratch/manuscript4$ squeue -u theaven
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
           4272754      long   braker  theaven  R      23:20      1 n24-64-384-giles
           4272751      long   braker  theaven  R      23:23      1 n24-64-384-giles
           4272752      long   braker  theaven  R      23:23      1 n24-64-384-giles
           4272753      long   braker  theaven  R      23:23      1 n24-64-384-giles
           4272748      long   braker  theaven  R      23:36      1 n24-64-384-gunn
           4272749      long   braker  theaven  R      23:36      1 n24-64-384-giles
           4272750      long   braker  theaven  R      23:36      1 n24-64-384-giles


# Gene prediction - orthodb11
for Assembly in $(ls /home/theaven/projects/niab/theaven/genomes/*.fasta); do
  ID="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@' | sed 's@_clean_renamed\.fasta@@')"

    Jobs=$(squeue -h -u theaven -n braker | wc -l)
  while [ "$Jobs" -gt 9 ]; do
    sleep 900s
    printf "."
    Jobs=$(squeue -h -u theaven -n braker | wc -l)
  done
  echo

  Masked_Genome="/home/theaven/scratch/manuscript4/reannotation/$ID/earlgrey_annotation/${ID}_EarlGrey/${ID}.softmasked.fa"
  RNA_alignment=NA
  Protein_database=~/scratch/manuscript4/db/orthodb11/Fungi.fa
  GeneMark=NA
  Min_contig=10000
  Braker_run=$(echo $ID)_uniprot
  OutDir="/home/theaven/scratch/manuscript4/reannotation/$ID/braker_prot"
  mkdir -p "$OutDir"
  if [ ! -f "$OutDir/braker.aa" ]; then
  rm -rf $OutDir/*
  rm -rf /mnt/shared/projects/niab/theaven/augustus_config/species/${Braker_run}
  jobid=$(sbatch --parsable ~/git_repos/Wrappers/gruffalo/run_braker3.sh "$Masked_Genome" "$Braker_run" "$OutDir" "$RNA_alignment" "$Protein_database" "$Min_contig" "$GeneMark")
  printf "%s\t%s\tbraker_prot_annotation\t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> slurm_log.tsv
  else
  echo "Found $OutDir/braker.aa"
  fi
done

sacct -j 4272789 --format=JobID,JobName,ReqMem,MaxRSS,TotalCPU,AllocCPUS,Elapsed,State,ExitCode,Ntasks,NCPUS,User

for job in $(squeue -u theaven | awk '{print $1}'); do
    cat slurm.$job.out
done

for job in $(squeue -u theaven | awk 'NR>1 && $6 ~ /^[0-9]+:[0-9][0-9]$/ {
    # Split the time field (mm:ss) into minutes and seconds
    split($6, t, ":")
    total_minutes = t[1] + t[2]/60
    if (total_minutes < 20) print $1
}'); do
    cat slurm.$job.out | grep -A 1 'Species:' | tail -n 1 | sed 's@_uniprot:@@g'
done


grep 'braker_abinitio_annotation' slurm_log.tsv | awk '{print $4}' | while read -r jobid; do
    sacct -j "$jobid" \
      --format=JobID,JobName,ReqMem,MaxRSS,TotalCPU,AllocCPUS,Elapsed,State,ExitCode,NTasks,NCPUS,User
done










conda activate busco
for predictions in $(ls /home/theaven/scratch/manuscript4/reannotation/*/*/braker.aa); do
  ID=$(echo "$predictions" | cut -d '/' -f7)
  Database="leotiomycetes_odb10"
  OutDir=$(dirname "$predictions")

  Jobs=$(squeue -h -u theaven -n busco | wc -l)
  while [ "$Jobs" -gt 9 ]; do
    sleep 900s
    printf "."
    Jobs=$(squeue -h -u theaven -n busco | wc -l)
  done
  echo

  if [ ! -f "$OutDir/${ID}_short_summary.txt" ]; then
  echo "Running for............. $OutDir/${ID}_short_summary.txt"
  jobid=$(sbatch --parsable ~/git_repos/Wrappers/gruffalo/run_busco_p.sh "$predictions" "$Database" "$OutDir" "$ID")
  printf "%s\t%s\tbusco\t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> slurm_log.tsv
  else
  echo "Found $OutDir/${ID}_short_summary.txt"
  fi
done
conda deactivate

for braker in $(ls /home/theaven/scratch/manuscript4/reannotation/*/*/braker.aa); do
  Dir=$(dirname $braker)
  ID=$(echo "$braker" | cut -d '/' -f7)
if [ ! -f "$Dir/${ID}_short_summary.txt" ]; then
  echo "$braker missing"
fi
done

conda activate omark
for predictions in $(ls /home/theaven/scratch/manuscript4/reannotation/*/*/braker.aa); do
  ID=$(echo "$predictions" | cut -d '/' -f7)
  Database="/home/theaven/scratch/manuscript4/db/oma/LUCA.h5"

  Jobs=$(squeue -h -u theaven -n omark | wc -l)
  while [ "$Jobs" -gt 9 ]; do
    sleep 60s
    printf "."
    Jobs=$(squeue -h -u theaven -n omark | wc -l)
  done
  echo

  OutDir=$(dirname "$predictions")/omark
  if [ ! -f "$OutDir/${ID}_detailed_summary.txt" ]; then
    mkdir "$OutDir"
  jobid=$(sbatch --parsable ~/git_repos/Wrappers/gruffalo/run_omark.sh "$predictions" "$Database" "$OutDir" "$ID")
  printf "%s\t%s\tomark\t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> slurm_log.tsv
  else
  echo "Found $OutDir/${ID}_detailed_summary.txt"
  fi
done
conda deactivate

python3 ~/git_repos/Scripts/gruffalo/parse_reannotation_results_2.py /home/theaven/scratch/manuscript4/reannotation --pivot -o results2.tsv

for predictions in $(ls /home/theaven/scratch/manuscript4/reannotation/*/*/braker.aa);
echo $predictions
grep '>' $predictions | wc -l
done
  ```

/home/theaven/scratch/manuscript4/genomes/DRT72021/Pod_aph_DRT72021-2_clean_renamed_fix.fasta
  ```bash
for Assembly in ~/scratch/manuscript4/genomes/*/*_fix.fasta; do
  [[ -e "$Assembly" ]] || { echo "No FASTA files found"; break; }
  ID="$(basename "$Assembly" | sed 's@_clean_renamed\.fasta@@')"
  ID2="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@')"
  echo "$ID"
  TPSI_DIR="/home/theaven/scratch/manuscript4/reannotation/${ID2}/TPSI"
bedtools getfasta \
  -fi "$Assembly" \
  -bed "$TPSI_DIR"/"$ID".TPSI.allHits.chains.bestPerLocus.gff3 \
  -fo "$TPSI_DIR"/"$ID2"_TPSI_hits.fa
cd-hit-est -i "$TPSI_DIR"/"$ID2"_TPSI_hits.fa -o "$TPSI_DIR"/"$ID2"_TPSI_clusters.fa -c 0.8 -n 5 -d 0 -T 8 -M 30000

mkdir "$TPSI_DIR"/cluster_members
awk -v out="$TPSI_DIR/cluster_members" '
  BEGIN { c = -1 }
  /^>Cluster[[:space:]][0-9]+/ { c++; next }
  {
    s = $0
    sub(/^.*>/, "", s)          # drop everything up to and including the first ">"
    sub(/\.\.\..*/, "", s)      # drop from the first "..." to end
    gsub(/[ \t\r]+$/, "", s)    # trim trailing whitespace just in case
    if (length(s) > 0) {
      printf("%s\n", s) >> sprintf("%s/cluster_%05d.ids", out, c)
    }
  }
' "$TPSI_DIR/${ID2}_TPSI_clusters.fa.clstr"

  for ids in "$TPSI_DIR"/cluster_members/*.ids; do
    cbase="$(basename "$ids" .ids)"
    seqkit grep -f "$ids" "$TPSI_DIR/${ID2}_TPSI_hits.fa" \
      > "$TPSI_DIR/cluster_members/${cbase}.fa" || true
  done

mkdir "$TPSI_DIR/cluster_consensus"
for fa in "$TPSI_DIR"/cluster_members/*.fa; do
  [[ -s "$fa" ]] || continue
  nseq=$(grep -c '^>' "$fa" || true)
  cbase=$(basename "$fa" .fa)
  aln="$TPSI_DIR/cluster_members/${cbase}.aln.fa"
  cons="$TPSI_DIR/cluster_consensus/${ID2}_${cbase}_cons.fa"

  if (( nseq >= 2 )); then
    mafft --auto --thread 8 "$fa" > "$aln"
    consambig -sequence "$aln" -outseq "$cons"
  else
    cp "$fa" "$cons"
    sed -i "1s/^>.*/>${ID2}_${cbase}_cons/" "$cons"
  fi
done

for cons in "$TPSI_DIR"/cluster_consensus/*.fa; do
  newhdr=">TPSI_"$(basename "$cons" .fa)
  awk -v h="$newhdr" 'NR==1{print h; next} {print}' "$cons" > "$cons.tmp" && mv "$cons.tmp" "$cons"
done

cat "$TPSI_DIR"/cluster_consensus/*.fa > "$TPSI_DIR/TPSI_consensus.fa"
echo "Wrote: $TPSI_DIR/TPSI_consensus.fa"

done


  ```

  ```bash
sbatch ~/git_repos/Wrappers/gruffalo/srun_trimmomatic.sh \
  /home/theaven/projects/niab/theaven/raw_data/RNA/s30008508_1.fq.gz \
  /home/theaven/projects/niab/theaven/raw_data/RNA/s30008508_2.fq.gz \
  ~/git_repos/Scripts/gruffalo/ncbi_adapters.fa \
  ~/scratch/manuscript4/rna_qc/RNA/strawberry \
  strawberry

sbatch ~/git_repos/Wrappers/gruffalo/srun_trimmomatic.sh \
  /home/theaven/projects/niab/theaven/raw_data/RNA/THMLAP201_1.fq.gz /home/theaven/projects/niab/theaven/raw_data/RNA/THMLAP202_1.fq.gz \
  /home/theaven/projects/niab/theaven/raw_data/RNA/THMLAP201_2.fq.gz /home/theaven/projects/niab/theaven/raw_data/RNA/THMLAP202_2.fq.gz \
  ~/git_repos/Scripts/gruffalo/ncbi_adapters.fa \
  ~/scratch/manuscript4/rna_qc/RNA/apple \
  apple



/home/theaven/scratch/manuscript4/genomes/OGB2019/Pod_leu_OGB2019-2_clean_renamed_fix.fasta
/home/theaven/scratch/manuscript4/genomes/OGB2021/Pod_leu_OGB2021-2_clean_renamed_fix.fasta
/home/theaven/projects/niab/theaven/genomes/Pod_leu_OGBp112020_clean_renamed.fasta

/home/theaven/scratch/manuscript4/genomes/SCOTT2020/Pod_aph_SCOTT2020-2_clean_renamed_fix.fasta
/home/theaven/scratch/manuscript4/genomes/DRT72021/Pod_aph_DRT72021-2_clean_renamed_fix.fasta
/home/theaven/projects/niab/theaven/genomes/Pod_aph_DRT72020_clean_renamed.fasta

for Assembly in $(ls /home/theaven/scratch/manuscript4/genomes/OGB2019/Pod_leu_OGB2019-2_clean_renamed_fix.fasta /home/theaven/scratch/manuscript4/genomes/OGB2021/Pod_leu_OGB2021-2_clean_renamed_fix.fasta /home/theaven/projects/niab/theaven/genomes/Pod_leu_OGBp112020_clean_renamed.fasta); do
ID="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@' | sed 's@_clean_renamed\.fasta@@')"
OutDir=/home/theaven/scratch/manuscript4/reannotation/"$ID" 
mkdir $OutDir
sbatch ~/git_repos/Wrappers/gruffalo/star_alignment.sh \
  $Assembly \
  /home/theaven/scratch/manuscript4/rna_qc/RNA/apple/F/apple_F_trim.fq.gz \
  /home/theaven/scratch/manuscript4/rna_qc/RNA/apple/R/apple_R_trim.fq.gz \
  $OutDir
done

for Assembly in $(ls /home/theaven/scratch/manuscript4/genomes/SCOTT2020/Pod_aph_SCOTT2020-2_clean_renamed_fix.fasta /home/theaven/scratch/manuscript4/genomes/DRT72021/Pod_aph_DRT72021-2_clean_renamed_fix.fasta /home/theaven/projects/niab/theaven/genomes/Pod_aph_DRT72020_clean_renamed.fasta); do
ID="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@' | sed 's@_clean_renamed\.fasta@@')"
OutDir=/home/theaven/scratch/manuscript4/reannotation/"$ID" 
mkdir $OutDir
sbatch ~/git_repos/Wrappers/gruffalo/star_alignment.sh \
  $Assembly \
  /home/theaven/scratch/manuscript4/rna_qc/RNA/strawberry/F/strawberry_F_trim.fq.gz \
  /home/theaven/scratch/manuscript4/rna_qc/RNA/strawberry/R/strawberry_R_trim.fq.gz \
  $OutDir
done

for Assembly in $(ls /home/theaven/scratch/manuscript4/genomes/SCOTT2020/Pod_aph_SCOTT2020-2_clean_renamed_fix.fasta /home/theaven/scratch/manuscript4/genomes/DRT72021/Pod_aph_DRT72021-2_clean_renamed_fix.fasta /home/theaven/projects/niab/theaven/genomes/Pod_aph_DRT72020_clean_renamed.fasta /home/theaven/scratch/manuscript4/genomes/OGB2019/Pod_leu_OGB2019-2_clean_renamed_fix.fasta /home/theaven/scratch/manuscript4/genomes/OGB2021/Pod_leu_OGB2021-2_clean_renamed_fix.fasta /home/theaven/projects/niab/theaven/genomes/Pod_leu_OGBp112020_clean_renamed.fasta); do
  ID="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@' | sed 's@_clean_renamed\.fasta@@')"

    Jobs=$(squeue -h -u theaven -n braker | wc -l)
  while [ "$Jobs" -gt 90 ]; do
    sleep 900s
    printf "."
    Jobs=$(squeue -h -u theaven -n braker | wc -l)
  done
  echo

  Masked_Genome="/home/theaven/scratch/manuscript4/reannotation/$ID/earlgrey_annotation/${ID}_EarlGrey/${ID}.softmasked.fa_backup"
  RNA_alignment=/home/theaven/scratch/manuscript4/reannotation/"$ID"/star_aligmentAligned.sortedByCoord.out.bam 
  Protein_database=NA
  GeneMark=NA
  Min_contig=10000
  Braker_run=$(echo $ID)_RNA
  OutDir="/home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA"
  mkdir -p "$OutDir"
  if [ ! -f "$OutDir/braker.aa" ]; then
  rm -rf $OutDir/*
  rm -rf /mnt/shared/projects/niab/theaven/augustus_config/species/${Braker_run}
  jobid=$(sbatch --parsable ~/git_repos/Wrappers/gruffalo/run_braker3.sh "$Masked_Genome" "$Braker_run" "$OutDir" "$RNA_alignment" "$Protein_database" "$Min_contig" "$GeneMark")
  printf "%s\t%s\tbraker_RNA_annotation\t%s\n" "$(date -Iseconds)" "$ID" "$jobid" >> slurm_log.tsv
  else
  echo "Found $OutDir/braker.aa"
  fi
done

#Check which genes have RNA support:
for Assembly in $(ls /home/theaven/scratch/manuscript4/genomes/SCOTT2020/Pod_aph_SCOTT2020-2_clean_renamed_fix.fasta /home/theaven/scratch/manuscript4/genomes/DRT72021/Pod_aph_DRT72021-2_clean_renamed_fix.fasta /home/theaven/projects/niab/theaven/genomes/Pod_aph_DRT72020_clean_renamed.fasta /home/theaven/scratch/manuscript4/genomes/OGB2019/Pod_leu_OGB2019-2_clean_renamed_fix.fasta /home/theaven/scratch/manuscript4/genomes/OGB2021/Pod_leu_OGB2021-2_clean_renamed_fix.fasta /home/theaven/projects/niab/theaven/genomes/Pod_leu_OGBp112020_clean_renamed.fasta); do
  ID="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@' | sed 's@_clean_renamed\.fasta@@')"

#get introns from the braker output:
grep -P "\tintron\t" /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/braker.gff3 > /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/braker.introns.gff
awk -F'\t' '$3=="intron"{print $1"\t"$4-1"\t"$5"\t.\t"$6"\t"$7}' /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/braker.introns.gff > /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/braker.introns.bed

#get introns from the hintsfile:
grep -P "\tintron\t|\tintron_hint\t" /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/hintsfile.gff > /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/rna_intron_hints.gff
awk -F'\t' '$3=="intron" || $3=="intron_hint"{print $1"\t"$4-1"\t"$5"\t.\t"$6"\t"$7}' /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/rna_intron_hints.gff > /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/rna_intron_hints.bed

#intersect the braker and hints introns:
bedtools intersect -u -s -f 1.0 -r -a /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/braker.introns.bed -b /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/rna_intron_hints.bed > /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/introns.supported.bed

#reverse intersect the braker and hints introns:
bedtools intersect -v -s -f 1.0 -r -a /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/braker.introns.bed -b /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/rna_intron_hints.bed > /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/introns.unsupported.bed

#get genes from the braker output:
awk -F'\t' '$3=="gene" {split($9,a,";"); for(i in a){if(a[i]~/^ID=/){gsub(/^ID=/,"",a[i]); print $1"\t"$4-1"\t"$5"\t"a[i]"\t0\t"$7;}}}' /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/braker.gff3 > /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/braker.genes.bed

#intersect the braker genes and hints introns:
bedtools intersect -u -s -a /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/braker.genes.bed -b /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/introns.supported.bed > /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/genes.with_rna_junction_support.bed

#reverse intersect the braker genes and hints introns:
bedtools intersect -v -s -a /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/braker.genes.bed -b /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/introns.supported.bed > /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/genes.no_rna_junction_support.bed

echo "$ID":
echo "$(wc -l < /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/introns.supported.bed) supported / $(wc -l < /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/braker.introns.bed) total introns"
echo "$(wc -l < /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/introns.unsupported.bed) unsupported / $(wc -l < /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/braker.introns.bed) total introns"
echo "$(wc -l < /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/genes.with_rna_junction_support.bed) supported / $(wc -l < /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/braker.genes.bed) total genes"
echo "$(wc -l < /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/genes.no_rna_junction_support.bed) unsupported / $(wc -l < /home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/braker.genes.bed) total genes"
done

#Pod_aph_DRT72020:
#13942 supported / 17743 total introns
#3801 unsupported / 17743 total introns
#5277 supported / 7328 total genes
#2051 unsupported / 7328 total genes
#Pod_leu_OGBp112020:
#18222 supported / 22866 total introns
#4644 unsupported / 22866 total introns
#6321 supported / 8074 total genes
#1753 unsupported / 8074 total genes
#Pod_aph_DRT72021-2:
#13746 supported / 17394 total introns
#3648 unsupported / 17394 total introns
#5210 supported / 7181 total genes
#1971 unsupported / 7181 total genes
#Pod_leu_OGB2019-2:
#17978 supported / 22756 total introns
#4778 unsupported / 22756 total introns
#6259 supported / 8038 total genes
#1779 unsupported / 8038 total genes
#Pod_leu_OGB2021-2:
#18371 supported / 23128 total introns
#4757 unsupported / 23128 total introns
#6360 supported / 8078 total genes
#1718 unsupported / 8078 total genes
#Pod_aph_SCOTT2020-2:
#13125 supported / 16689 total introns
#3564 unsupported / 16689 total introns
#4983 supported / 6943 total genes
#1960 unsupported / 6943 total genes

  for Assembly in $(ls /home/theaven/scratch/manuscript4/genomes/SCOTT2020/Pod_aph_SCOTT2020-2_clean_renamed_fix.fasta /home/theaven/scratch/manuscript4/genomes/DRT72021/Pod_aph_DRT72021-2_clean_renamed_fix.fasta /home/theaven/projects/niab/theaven/genomes/Pod_aph_DRT72020_clean_renamed.fasta /home/theaven/scratch/manuscript4/genomes/OGB2019/Pod_leu_OGB2019-2_clean_renamed_fix.fasta /home/theaven/scratch/manuscript4/genomes/OGB2021/Pod_leu_OGB2021-2_clean_renamed_fix.fasta /home/theaven/projects/niab/theaven/genomes/Pod_leu_OGBp112020_clean_renamed.fasta); do
    ID="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@' | sed 's@_clean_renamed\.fasta@@')"
    OutDir="/home/theaven/scratch/manuscript4/reannotation/$ID/stringtie"
    AcceptedHits=/home/theaven/scratch/manuscript4/reannotation/"$ID"/star_aligmentAligned.sortedByCoord.out.bam 
    ProgDir=~/git_repos/Wrappers/gruffalo
    sbatch $ProgDir/stringtie.sh $AcceptedHits $OutDir
   done

  for Assembly in $(ls /data/users/theaven/down_02012025/Pod_aph_SCOTT2020-2/Pod_aph_SCOTT2020-2_clean_renamed_fix.fasta /data/users/theaven/down_02012025/Pod_aph_DRT72021-2/Pod_aph_DRT72021-2_clean_renamed_fix.fasta /data/users/theaven/down_02012025/Pod_aph_DRT72020/Pod_aph_DRT72020_clean_renamed.fasta /data/users/theaven/down_02012025/Pod_leu_OGB2019-2/Pod_leu_OGB2019-2_clean_renamed_fix.fasta /data/users/theaven/down_02012025/Pod_leu_OGB2021-2/Pod_leu_OGB2021-2_clean_renamed_fix.fasta /data/users/theaven/down_02012025/Pod_leu_OGBp112020/Pod_leu_OGBp112020_clean_renamed.fasta); do
    ID="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@' | sed 's@_clean_renamed\.fasta@@')"
    OutDir="/data/users/theaven/down_02012025/$ID/codingquarry"
    GTF=/data/users/theaven/down_02012025/$ID/out.gtf
    mkdir $OutDir
    cd $OutDir
    CufflinksGTF_to_CodingQuarryGFF3.py $GTF > transcripts.gff3
    echo "Running the standard CodingQuarry predictions..."
    CodingQuarry -f $Assembly -t transcripts.gff3 -p 1
    echo "Translating and quality filtering the coding sequence..."
    python /mnt/apps/users/theaven/conda/envs/py2/opt/codingquarry-2.0/QuarryFiles/scripts/fastaTranslate.py out/Predicted_CDS.fa | sed 's/*$//g' > CQ_Proteins.fa
    python /mnt/apps/users/theaven/conda/envs/py2/opt/codingquarry-2.0/QuarryFiles/scripts/gene_errors_Xs.py CQ_Proteins.fa CQPMtmp.fa
    mv CQPMtmp.fa CQ_Proteins.fa
    echo "Running signalP..."
    python /mnt/apps/users/theaven/conda/envs/py2/opt/codingquarry-2.0/QuarryFiles/scripts/split_fasta.py CQ_Proteins.fa 200
i=0

SIGNALP="~/signalp/signalp-4.1"

    for FILE in CQ_Proteins.fa-*
    do
       $SIGNALP/signalp $FILE > CQ_Proteins_out_$i
       rm $FILE
        i=$(($i+1))
    done
    cat CQ_Proteins_out_* | grep -v "#" | awk '($10 == "Y"){print $1" "$5}' > Secretome.txt
    rm CQ_Proteins_out_*
    rm CQ_Proteins.fa
    echo "Running CodingQuarry-PM..."
    CodingQuarry -f $Assembly -t transcripts.gff3 -2 out/PredictedPass.gff3 -p 1 -g Secretome.txt -h
    cd /home/theaven/scratch/manuscript4
  done

srun -p long -J CQ  --mem-per-cpu 4G --cpus-per-task 1 --pty bash
conda activate py2
for Assembly in $(ls /home/theaven/scratch/manuscript4/genomes/SCOTT2020/Pod_aph_SCOTT2020-2_clean_renamed_fix.fasta /home/theaven/projects/niab/theaven/genomes/Pod_aph_DRT72020_clean_renamed.fasta /home/theaven/scratch/manuscript4/genomes/OGB2019/Pod_leu_OGB2019-2_clean_renamed_fix.fasta /home/theaven/scratch/manuscript4/genomes/OGB2021/Pod_leu_OGB2021-2_clean_renamed_fix.fasta /home/theaven/projects/niab/theaven/genomes/Pod_leu_OGBp112020_clean_renamed.fasta); do
  ID="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@' | sed 's@_clean_renamed\.fasta@@')"
  OutDir="/home/theaven/scratch/manuscript4/reannotation/$ID/codingquarry"
  GTF=/home/theaven/scratch/manuscript4/reannotation/$ID/stringtie/out.gtf
  mkdir $OutDir
  cd $OutDir
  CufflinksGTF_to_CodingQuarryGFF3.py $GTF > transcripts.gff3
  echo "Running the standard CodingQuarry predictions..."
  CodingQuarry -f $Assembly -t transcripts.gff3 -p 1
  echo "Translating and quality filtering the coding sequence..."
  python /mnt/apps/users/theaven/conda/envs/py2/opt/codingquarry-2.0/QuarryFiles/scripts/fastaTranslate.py out/Predicted_CDS.fa | sed 's/*$//g' > CQ_Proteins.fa
  python /mnt/apps/users/theaven/conda/envs/py2/opt/codingquarry-2.0/QuarryFiles/scripts/gene_errors_Xs.py CQ_Proteins.fa CQPMtmp.fa
  mv CQPMtmp.fa CQ_Proteins.fa
  echo "Running signalP..."
done

conda activate signalp6
for In in $(ls /home/theaven/scratch/manuscript4/reannotation/*/codingquarry/CQ_Proteins.fa); do
  Out=$(dirname $In)/signalp6
sbatch ~/git_repos/Wrappers/gruffalo/signalp6.sh $In $Out
done

for In in $(ls /home/theaven/scratch/manuscript4/reannotation/*/braker_RNA/braker.aa); do
  Out=$(dirname $In)/signalp6
sbatch ~/git_repos/Wrappers/gruffalo/signalp6-fast.sh $In $Out
done #8498550-4,8508717
conda deactivate

for In in $(ls /home/theaven/scratch/manuscript4/reannotation/*/braker_RNA/braker.aa); do
  cd $(dirname "$In")
  P6=$(dirname "$In")/signalp6/prediction_results.txt
gawk -F'\t' '$1 !~ /^#/ && $2=="SP" {if (match($0, /CS pos:[[:space:]]*([0-9]+)-[0-9]+/, m)) print $1" "m[1]; else print $1" 0";}' "$P6" > Secretome2.txt

done

for Assembly in $(ls /home/theaven/scratch/manuscript4/genomes/SCOTT2020/Pod_aph_SCOTT2020-2_clean_renamed_fix.fasta /home/theaven/projects/niab/theaven/genomes/Pod_aph_DRT72020_clean_renamed.fasta /home/theaven/scratch/manuscript4/genomes/OGB2019/Pod_leu_OGB2019-2_clean_renamed_fix.fasta /home/theaven/scratch/manuscript4/genomes/OGB2021/Pod_leu_OGB2021-2_clean_renamed_fix.fasta /home/theaven/projects/niab/theaven/genomes/Pod_leu_OGBp112020_clean_renamed.fasta); do
  ID="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@' | sed 's@_clean_renamed\.fasta@@')"
  echo $ID
  OutDir="/home/theaven/scratch/manuscript4/reannotation/$ID/codingquarry"
  cd "$OutDir"
  echo "Running CodingQuarry-PM..."
  CodingQuarry -f $Assembly -t transcripts.gff3 -2 /home/theaven/scratch/manuscript4/reannotation/"$ID"/braker_RNA/braker.gff3 -p 1 -g /home/theaven/scratch/manuscript4/reannotation/"$ID"/braker_RNA/Secretome2.txt -h
done

/home/theaven/scratch/manuscript4/reannotation/Pod_aph_DRT72020/codingquarry/CQ_Proteins.fa
16163
/home/theaven/scratch/manuscript4/reannotation/Pod_aph_DRT72021-2/codingquarry/CQ_Proteins.fa
13967
/home/theaven/scratch/manuscript4/reannotation/Pod_aph_SCOTT2020-2/codingquarry/CQ_Proteins.fa
14403
/home/theaven/scratch/manuscript4/reannotation/Pod_leu_OGB2019-2/codingquarry/CQ_Proteins.fa
18923
/home/theaven/scratch/manuscript4/reannotation/Pod_leu_OGB2021-2/codingquarry/CQ_Proteins.fa
17917
/home/theaven/scratch/manuscript4/reannotation/Pod_leu_OGBp112020/codingquarry/CQ_Proteins.fa
18943



cd apha
orthofinder -f proteomes -t 16 -a 16

ulimit -n 52000
orthofinder -S blast -f ./ -t 32 -a 8 -n $prefix -o $OutDir

for file in $(ls apha/signalp/DRT72020/prediction_results.txt apha/signalp/DRT72021/prediction_results.txt apha/signalp/SCOTT2020/prediction_results.txt leuco/signalp/OGB2019/prediction_results.txt leuco/signalp/OGB2021/prediction_results.txt leuco/signalp/P112020/prediction_results.txt
); do
awk -F'\t' '$1!~/^#/ && $2=="SP"{ if(match($(NF),/CS pos: *([0-9]+)-[0-9]+/,m)) print $1"\t"m[1]; else print $1"\t0" }' \
  "$file" > "$(dirname $file)"/sp_cs.tsv
done

OG_TSV=$(ls -d apha/proteomes/OrthoFinder/Results_Jan09/Orthogroups/Orthogroups.tsv | head -n1)
echo "$OG_TSV"
head -n 1 "$OG_TSV" | tr '\t' '\n'
python make_cq_secretome_from_orthofinder.py "$OG_TSV" apha/signalp DRT72020 apha/DRT72020.Secretome.txt
python make_cq_secretome_from_orthofinder.py "$OG_TSV" apha/signalp DRT72021 apha/DRT72021.Secretome.txt
python make_cq_secretome_from_orthofinder.py "$OG_TSV" apha/signalp SCOTT2020 apha/SCOTT2020.Secretome.txt

OG_TSV=$(ls -d leuco/proteomes/OrthoFinder/Results_Jan09/Orthogroups/Orthogroups.tsv | head -n1)
echo "$OG_TSV"
head -n 1 "$OG_TSV" | tr '\t' '\n'
python make_cq_secretome_from_orthofinder.py "$OG_TSV" leuco/signalp OGB2019 leuco/OGB2019.Secretome.txt
python make_cq_secretome_from_orthofinder.py "$OG_TSV" leuco/signalp OGB2021 leuco/OGB2021.Secretome.txt
python make_cq_secretome_from_orthofinder.py "$OG_TSV" leuco/signalp P112020 leuco/P112020.Secretome.txt

cp apha/DRT72020.Secretome.txt /home/theaven/scratch/manuscript4/reannotation/Pod_aph_DRT72020/codingquarry/Secretome2.txt
cp apha/DRT72021.Secretome.txt /home/theaven/scratch/manuscript4/reannotation/Pod_aph_DRT72021-2/codingquarry/Secretome2.txt
cp apha/SCOTT2020.Secretome.txt /home/theaven/scratch/manuscript4/reannotation/Pod_aph_SCOTT2020-2/codingquarry/Secretome2.txt
cp leuco/OGB2019.Secretome.txt /home/theaven/scratch/manuscript4/reannotation/Pod_leu_OGB2019-2/codingquarry/Secretome2.txt
cp leuco/OGB2021.Secretome.txt /home/theaven/scratch/manuscript4/reannotation/Pod_leu_OGB2021-2/codingquarry/Secretome2.txt
cp leuco/P112020.Secretome.txt /home/theaven/scratch/manuscript4/reannotation/Pod_leu_OGBp112020/codingquarry/Secretome2.txt

for Assembly in $(ls /home/theaven/scratch/manuscript4/genomes/SCOTT2020/Pod_aph_SCOTT2020-2_clean_renamed_fix.fasta /home/theaven/projects/niab/theaven/genomes/Pod_aph_DRT72020_clean_renamed.fasta /home/theaven/scratch/manuscript4/genomes/OGB2019/Pod_leu_OGB2019-2_clean_renamed_fix.fasta /home/theaven/scratch/manuscript4/genomes/OGB2021/Pod_leu_OGB2021-2_clean_renamed_fix.fasta /home/theaven/projects/niab/theaven/genomes/Pod_leu_OGBp112020_clean_renamed.fasta); do
  ID="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@' | sed 's@_clean_renamed\.fasta@@')"
  echo $ID
  OutDir="/home/theaven/scratch/manuscript4/reannotation/$ID/codingquarry"
  cd "$OutDir"
  echo "Running CodingQuarry-PM..."
  CodingQuarry -f $Assembly -t transcripts.gff3 -2 out/PredictedPass.gff3 -p 1 -g Secretome.txt -h
done


for In in $(ls /home/theaven/scratch/manuscript4/reannotation/*/codingquarry/CQ_Proteins.fa); do
  Out=$(dirname $In)/signalp6-fast
sbatch ~/git_repos/Wrappers/gruffalo/signalp6-fast.sh $In $Out
done #8455918-23

signalp6 --fastafile  --organism eukarya --mode fast --output_dir .

signalp6 --fastafile /home/theaven/scratch/manuscript4/reannotation/Pod_aph_DRT72020/codingquarry/CQ_Proteins.fa --organism eukarya --mode fast --model_dir /home/theaven/signalp6_fast/signalp-6-package/models --torch_num_threads 1 --bsize 4 --write_procs 1 --output_dir .



    cat CQ_Proteins_out_* | grep -v "#" | awk '($10 == "Y"){print $1" "$5}' > Secretome.txt
    rm CQ_Proteins_out_*
    rm CQ_Proteins.fa
    echo "Running CodingQuarry-PM..."
    CodingQuarry -f $Assembly -t transcripts.gff3 -2 out/PredictedPass.gff3 -p 1 -g Secretome.txt -h
    cd /home/theaven/scratch/manuscript4
  done

singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/codingquarry_2.0--py312hf731ba3_11 CufflinksGTF_to_CodingQuarryGFF3.py $GTF > transcripts.gff3
cp $Assembly assembly.fa

mkdir out
singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/codingquarry_2.0--py312hf731ba3_11 run_CQ-PM_stranded.sh transcripts.gff3 $Threads 2>&1 | tee codingquaryPM_log.txt


#!/bin/bash

#This is the prediction from transcript run script
#for stranded RNA-seq. This script automates
#CodingQuarry's pathogen mode (see manual).

#Provide a genome and gff of transcripts to run:
#run_CQ-PM_stranded.sh myGenome.fa myTranscripts.gff
#See the manual for input format requirements.

echo "Running the standard CodingQuarry predictions..."
CodingQuarry -f $1 -t $2 -p 1

echo "Translating and quality filtering the coding sequence..."
python /mnt/apps/users/theaven/conda/envs/py2/opt/codingquarry-2.0/QuarryFiles/scripts/fastaTranslate.py out/Predicted_CDS.fa | sed 's/*$//g' > CQ_Proteins.fa
python /mnt/apps/users/theaven/conda/envs/py2/opt/codingquarry-2.0/QuarryFiles/scripts/gene_errors_Xs.py CQ_Proteins.fa CQPMtmp.fa
mv CQPMtmp.fa CQ_Proteins.fa

echo "Running signalP..."
#Split the protein file up into smaller chunks
python /mnt/apps/users/theaven/conda/envs/py2/opt/codingquarry-2.0/QuarryFiles/scripts/split_fasta.py CQ_Proteins.fa 200
i=0
for FILE in CQ_Proteins.fa-*
do
    #If signalP is not in your path, change the line below to specify its location
    signalp $FILE > CQ_Proteins_out_$i
    rm $FILE
    i=$(($i+1))
done
cat CQ_Proteins_out_* | grep -v "#" | awk '($10 == "Y"){print $1" "$5}' > Secretome.txt
rm CQ_Proteins_out_*
rm CQ_Proteins.fa

echo "Running CodingQuarry-PM..."
CodingQuarry -f $1 -t $2 -2 out/PredictedPass.gff3 -p 1 -g Secretome.txt -h


  ```
  ```bash
srun -p short  --mem 4G -c 1 --pty bash
cpanm Bio::Perl

conda activate bedtools
for Assembly in $(ls /home/theaven/scratch/manuscript4/genomes/SCOTT2020/Pod_aph_SCOTT2020-2_clean_renamed_fix.fasta /home/theaven/scratch/manuscript4/genomes/DRT72021/Pod_aph_DRT72021-2_clean_renamed_fix.fasta /home/theaven/projects/niab/theaven/genomes/Pod_aph_DRT72020_clean_renamed.fasta /home/theaven/scratch/manuscript4/genomes/OGB2019/Pod_leu_OGB2019-2_clean_renamed_fix.fasta /home/theaven/scratch/manuscript4/genomes/OGB2021/Pod_leu_OGB2021-2_clean_renamed_fix.fasta /home/theaven/projects/niab/theaven/genomes/Pod_leu_OGBp112020_clean_renamed.fasta); do
  ID="$(basename "$Assembly" | sed 's@_clean_renamed_fix\.fasta@@' | sed 's@_clean_renamed\.fasta@@')"
  BrakerGff=/home/theaven/scratch/manuscript4/reannotation/$ID/braker_RNA/braker.gff3
    if [ ! -f "/home/theaven/scratch/manuscript4/reannotation/$ID/earlgrey_annotation/${ID}_EarlGrey/${ID}.softmasked.fa_backup" ]; then
  Masked_Genome="/home/theaven/scratch/manuscript4/reannotation/$ID/earlgrey_annotation/${ID}_EarlGrey/${ID}.softmasked.fa_backup"
else
  Masked_Genome="/home/theaven/scratch/manuscript4/reannotation/$ID/earlgrey_annotation/${ID}_EarlGrey/${ID}.softmasked.fa"
fi
  CodingQuaryGff=/home/theaven/scratch/manuscript4/reannotation/${ID}/codingquarry/out/PredictedPass.gff3
  PGNGff=/home/theaven/scratch/manuscript4/reannotation/${ID}/codingquarry/out/PGN_predictedPass.gff3
  AddDir=$(dirname $CodingQuaryGff)/additional
  FinalDir=$(dirname $CodingQuaryGff)/final
  AddGenesList=$AddDir/additional_genes.txt
  AddGenesGff=$AddDir/additional_genes.gff
  FinalGff=$AddDir/combined_genes.gff
  mkdir -p $AddDir
  mkdir -p $FinalDir

#Create a list with the additional transcripts in CodingQuarry gff vs Braker gene models
bedtools intersect -v -s -a $CodingQuaryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
bedtools intersect -v -s -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList

#Create Gff file with the additional transcripts
ProgDir=/home/theaven/git_repos/Scripts/gruffalo
perl $ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff

#Create a final Gff file with gene features
perl $ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuary.gff3

#Create fasta files from each gene feature in the CodingQuarry gff3
perl $ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary

#Create fasta files from each gene feature in the Braker gff3
cp $BrakerGff $FinalDir/final_genes_Braker.gff3
perl $ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker

#Combine both fasta files
cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta

#Combine both gff3 files
GffBraker=$FinalDir/final_genes_Braker.gff3
GffQuary=$FinalDir/final_genes_CodingQuary.gff3
GffAppended=$FinalDir/final_genes_appended.gff3
cat $GffBraker $GffQuary > $GffAppended
done

#Check the final number of genes 
for DirPath in $(ls -d $FinalDir); do echo $DirPath; cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l; cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l; cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l; echo ""; 
done
#gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final
#Braker genes: 15642
#Coding quarry: 2987
#Combined: 18629
  ```

```bash
nextflow run \
  -with-singularity "~/git_repos/Containers/predector.sif" \
  -profile c16,r32 \
  -resume \
  -r 1.2.6 \
  -with-conda /mnt/shared/scratch/theaven/apps/conda/envs/predector \
  ccdmb/predector \
  --proteome $proteome 2>&1 | tee -a predector_report.txt

/home/theaven/nextflow run \
  ccdmb/predector \
  -with-singularity "~/git_repos/Containers/predector.sif" \
  -profile c16,r32 \
  -resume \
  -r 1.2.6 \
  --proteome /home/theaven/scratch/manuscript4/reannotation/Pod_leu_OGBp112020/braker_RNA/braker.aa 2>&1 | tee -a predector_report.txt

  /home/theaven/nextflow run -profile test -with-singularity ~/git_repos/Containers/predector.sif -resume -r 1.2.7 ccdmb/predector

# or if you've build the container using docker and it's in your local docker registry.
nextflow run -profile test,singularity -resume -r 1.2.7 ccdmb/predector
```