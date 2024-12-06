# Powdery mildew comparative genomics pipeline

Documentation of analysis and commands used for comparative genomics work as part of Nano Diagnostics PhD. - This file documents effector prediciton and comparisons across the powdery mildew clade.

Genomes assembled and annotated in nano_diagnostics/P_aphanis, nano_diagnostics/P_leucotricha and nano_diagnostics/B_graminis were used alongside genomes downloaded from NCBI for the following analysis.

All the following commands were executed on the Crop Diversity HPC from the folder /home/theaven/scratch unless stated otherwise.

The original notes for this work have become unbelievably convoluted and hard to follow as pieces of work ahave been restarted from mid way through and new genomes have been included etc. This is an edited version that will hopefully have better readability whilst ommiting now redundant work.

Many output files are now located within /home/theaven/projects/niab/theaven/

## Collecting data

Genomes were initially downloaded to the NIAB HPC, these were subsequently moved to the Crop Diversity HPC and all data in the /data/scratch/heavet directory of the NIAB HPC has since been deleted to save space.

Our powdery mildew genomes:
```bash
ls /home/theaven/scratch/uncompressed/genomes/*/*/*/*a
/home/theaven/scratch/uncompressed/genomes/Amorphotheca/resinae/GCA_003019875.1/amorphothecaresinae_GCA_003019875.1_Amore1_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Amorphotheca/resinae/GCA_018167515.1/Amorphothecaresinae_GCA_018167515.1_ASM1816751v1_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Arachnopeziza/araneosa/GCA_003988855.1/arachnopezizaaraneosa_GCA_003988855.1_ASM398885v1_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Ascocoryne/sarcoides/GCA_000328965.1/ascocorynesarcoides_GCA_000328965.1_ASM32896v1_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Aspergillus/nidulans/GCA_000149205.2/Aspergillusnidulans_GCA_000149205.2_ASM14920v2_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000151065.3/GCA_000151065.3_ASM15106v3_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000417025.1/Blumeriagraminis_GCA_000417025.1_ASM41702v1_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000417865.1/Blumeriagraminis_GCA_000417865.1_Bgt_94202_1_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000418435.1/Blumeriagraminis_GCA_000418435.1_Bgt_454_newbler_assembly_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000441875.1/Blumeriagraminis_GCA_000441875.1_Bgt_70_1_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_900519115.1/Blumeriagraminis_GCA_900519115.1_Bgt_genome_v3.16_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_905067625.1/Blumeriagraminis_GCA_905067625.1_Bgtriticale_THUN12_genome_v1_2_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153116/S_1459_contigs_min_500bp.c.fna
/home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153117/S_1201_contigs_min_500bp.c.fna
/home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153118/S_1203_contigs_min_500bp.c.fna
/home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153119/S_1391_contigs_min_500bp.c.fna
/home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153120/S_1400_contigs_min_500bp.c.fna
/home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_000401675.1/Blumeriahordei_GCA_000401675.1_A6_trimmed_clc_assembly_flt.fa_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_900237765.1/Blumeriahordei_GCA_900237765.1_BghRACE1_v1_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_900239735.1/Blumeriahordei_GCA_900239735.1_BGH_DH14_v4_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_900638725.1/Blumeriahordei_GCA_900638725.1_BGH_K1_v2_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Botrytis/cineria/GCA_000143535.4/GCA_000143535.4_ASM14353v4_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Chlorociboria/aeruginascens/GCA_002276475.2/chlorociboriaaeruginascens_GCA_002276475.2_Chloro_aeru_v2_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Colletotrichum/higginsianum/broad_KN1394/KN-1394_23MarchAllContigs.fna
/home/theaven/scratch/uncompressed/genomes/Colletotrichum/higginsianum/GCA_001672515.1/Colletotrichumhigginsianum_GCA_001672515.1_ASM167251v1_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Coprinus/cinereus/GCA_000182895.1/GCA_000182895.1_CC3_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Cryptococcus/neoformans/GCF_000091045.1/GCF_000091045.1_ASM9104v1_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Drepanopeziza/brunnea/GCA_000298775.1/drepanopezizabrunnea_GCA_000298775.1_ASM29877v1_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Erysiphe/alphitoides/CLCBIO/Erysiphealphitoides_CLCBIO_assembly_cdhitest_0.95_morethan500bp.fna
/home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798715.1/Erysiphenecator_GCA_000798715.1_ASM79871v1_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798735.1/Erysiphenecator_GCA_000798735.1_ASM79873v1_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798755.1/Erysiphenecator_GCA_000798755.1_ASM79875v1_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798775.1/Erysiphenecator_GCA_000798775.1_ASM79877v1_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798795.1/Erysiphenecator_GCA_000798795.1_ASM79879v1_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_016906895.1/Erysiphenecator_GCA_016906895.1_ASM1690689v1_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_024703715.1/GCA_024703715.1_EnE101_chr_assembly_v3_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Erysiphe/neolycopersici/GCA_003610855.1/Erysipheneolycopersici_GCA_003610855.1_ASM361085v1_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Erysiphe/pisi/GCA_000208805.1/Erysiphepisi_GCA_000208805.1_ASM20880v1_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Erysiphe/pisi/GCA_000214055.1/Erysiphepisi_GCA_000214055.1_ASM21405v1_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Erysiphe/pulchra/GCA_002918395.1/Erysiphepulchra_GCA_002918395.1_ASM291839v1_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Fusarium/graminearum/GCA_000240135.3/GCA_000240135.3_ASM24013v3_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Fusarium/oxysporum/GCF_013085055.1/GCF_013085055.1_ASM1308505v1_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Glarea/lozoyensis/GCA_000409485.1/glareaiozoyensis_GCA_000409485.1_GLAREA_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Golovinomyces/cichoracearum/GCA_003611195.1/Golovinomycescichoracearum_GCA_003611195.1_ASM361119v1_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Golovinomyces/cichoracearum/GCA_003611215.1/Golovinomycescichoracearum_GCA_003611215.1_ASM361121v1_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Golovinomyces/cichoracearum/GCA_003611235.1/Golovinomycescichoracearum_GCA_003611235.1_ASM361123v11_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Golovinomyces/magnicellulatus/GCA_006912115.1/Golovinomycesmagnicellulatus_GCA_006912115.1_ASM691211v1_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Golovinomyces/orontii/MGH1/Golovinomycesorontii_AssemblyScaffolds.fasta
/home/theaven/scratch/uncompressed/genomes/Leveillula/taurrica/CADEPA01/leveillulataurrica_CADEPA01.fasta
/home/theaven/scratch/uncompressed/genomes/Magnaporthe/oryzae/GCF_000002495.2/GCF_000002495.2_MG8_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Melampsora/laricis-populina/GCA_000204055.1/Melampsoralaricis-populina_GCA_000204055.1_v1.0_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Mollisia/scopiformis/GCA_001500285.1/mollisiascopiformis_GCA_001500285.1_Phisc1_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Neobulgaria/alba/GCA_003988965.1/neobulgariaablba_GCA_003988965.1_ASM398896v1_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Neurospora/crassa/GCA_000182925.2/Neurosporacrassa_GCA_000182925.2_NC12_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Oidiodendron/maius/GCA_000827325.1/Oidiodendron_maius_GCA_000827325.1_Zn_v1.0_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Oidium/heveae/GCA_003957845.1/Oidiumheveae_GCA_003957845.1_ASM395784v1_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Parauncinula/polyspora/Parp01/Parpo1_AssemblyScaffolds.fna
/home/theaven/scratch/uncompressed/genomes/Phialocephala/subalpina/GCA_900073065.1/phialocephalasubalpina-GCA_900073065.1_PAC_version_1_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Phyllactinia/moricola/GCA_019455665.1/phyllactiniamoricola_GCA_019455665.1_ASM1945566v1_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Pleochaeta/shiraiana/GCA_019455505.1/pleochaetashiraiana_GCA_019455505.1_ASM1945550v1_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Pleurotus/ostreatus/GCA_014466165.1/GCA_014466165.1_ASM1446616v1_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/DRT72020/HEAVEN_strawberry2020.fna
/home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/DRT72021/HEAVEN_strawberry2021.fna
/home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/SCOTT2020/HEAVEN_raspberry2020.fna
/home/theaven/scratch/uncompressed/genomes/Podosphaera/cerasii/GCA_018398735.1/Podosphaeracerasii_GCA_018398735.1_WSU_PS_Poc_1.0_genomic.1c.fna
/home/theaven/scratch/uncompressed/genomes/Podosphaera/fusca/GCA_030378345.1/GCA_030378345.1_ASM3037834v1_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/GCA_013170925.1/GANAN_apple.fna
/home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGB2019/HEAVEN_apple2019.fna
/home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGB2021/HEAVEN_apple2021.fna
/home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGBp112020/HEAVEN_apple2020.fna
/home/theaven/scratch/uncompressed/genomes/Podosphaera/xanthii/GCA_010015925.1/Podosphaeraxanthii_GCA_010015925.1_POXAN_niab_assembly_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Podosphaera/xanthii/GCA_014884795.1/Podosphaeraxanthii_GCA_014884795.1_ASM1488479v1_genomic.c.fna
/home/theaven/scratch/uncompressed/genomes/Podosphaera/xanthii/GCA_028751805.1/GCA_028751805.1_YZU_CsPM_1.2_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Psilocybe/cubensis/GCA_017499595.2/GCA_017499595.2_MGC_Penvy_1_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Puccinia/graminis/GCA_000149925.1/GCA_000149925.1_ASM14992v1_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Puccinia/striiformis/GCA_021901695.1/GCA_021901695.1_Pst134E36_v1_pri_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Puccinia/triticina/GCA_000151525.2/GCA_000151525.2_P_triticina_1_1_V2_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Pyricularia/oryzae/GCA_000002495.2/Pyriculariaoryzae_GCA_000002495.2_MG8_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Saccharomyces/cerevisiae/GCF_000146045.2/Scereviseae_GCF_000146045.2_R64_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Schizophyllum/commune/GCA_000143185.2/GCA_000143185.2_Schco3_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Sclerotinia/sclerotiorum/GCA_000146945.2/GCA_000146945.2_ASM14694v2_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Sclerotinia/sclerotiorum/GCF_000146945.1/GCF_000146945.1_ASM14694v1_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Tuber/melanosporum/GCF_000151645.1/GCF_000151645.1_ASM15164v1_genomic.fna
/home/theaven/scratch/uncompressed/genomes/Ustilago/maydis/GCA_000328475.2/GCA_000328475.2_Umaydis521_2.0_genomic.fna
```

## Phylogeny

#### BUSCO
```bash
for assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/Genomes/*/*/*/*a); do
    Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/busco/leotiomycetes_odb10
    ProgDir=~/git_repos/Wrappers/NBI
    OutDir=$(dirname $assembly)/BUSCO
    OutFile=$(echo $assembly | cut -d '/' -f9,10,11 | sed 's@/@_@g')
    Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    echo x
    while [ $Jobs -gt 3 ]; do
      sleep 900s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'busco'  | wc -l)
    done
    mkdir $OutDir 
    sleep 30s
    sbatch $ProgDir/run_busco_keep.sh $assembly $Database $OutDir $OutFile
done

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/Genomes/*/*/*/BUSCO/leotiomycetes_odb10/run_leotiomycetes_odb10/full_table.tsv); do
grep -v "^#" $file | awk '$2=="Complete" {print $1}' >> /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leotiomycetes_complete_busco_ids.txt;
done

sort /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leotiomycetes_complete_busco_ids.txt |uniq -c > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leotiomycetes_complete_busco_ids_with_counts.txt
grep -v " 2 " /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leotiomycetes_complete_busco_ids_with_counts.txt | grep -v " 1 " > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leotiomycetes_complete_busco_ids.txt
awk '{print $2}' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leotiomycetes_complete_busco_ids.txt > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leotiomycetes_complete_busco_ids_3.txt


mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leo_busco_nt
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/Genomes/*/*/*/BUSCO/leotiomycetes_odb10/run_leotiomycetes_odb10/busco_sequences/single_copy_busco_sequences.tar.gz); do
cd $(dirname $file)
tar -xzvf single_copy_busco_sequences.tar.gz
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/Genomes
done

#Give unique names to the complete busco genes from each assembly:
for dir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/Genomes/*/*/*/BUSCO/leotiomycetes_odb10/run_leotiomycetes_odb10/busco_sequences/single_copy_busco_sequences); do
  sppname=$(echo $dir |cut -f9,10,11 -d "/" | sed 's@/@_@g');
  abbrv=$(echo $dir | cut -d '/' -f9 | cut -c 1-3)_$(echo $dir | cut -d '/' -f10 | cut -c 1-3)_$(echo $dir | cut -d '/' -f11)
  echo $sppname
  echo $abbrv
  for file in ${dir}/*.fna; do
    out=$(echo $file |rev |cut -f 1 -d "/"|rev)
    cp $file /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leo_busco_nt/${sppname}_${out}
    sed -i 's/^>/>'${abbrv}'|/g' /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leo_busco_nt/${sppname}_${out}
  cut -f 1 -d ":" /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leo_busco_nt/${sppname}_${out} | tr '[:lower:]' '[:upper:]' > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leo_busco_nt/${sppname}_${out}.1 && mv /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leo_busco_nt/${sppname}_${out}.1 /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leo_busco_nt/${sppname}_${out}  
  done
done

#Combine genes from each assembly into a single file per gene:
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leo_busco_nt
buscos=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leotiomycetes_complete_busco_ids_3.txt
lines=$(cat $buscos)
for line in $lines; do
  for fna in $(ls *_$line.fna); do
  output=$(echo $line)_nt.fasta
  cat $fna >> $output
  done
done
rm *.fna
```
#### MAFFT
```bash
#Align the gene sequences;
AlignDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leo_busco_nt
for file in $(ls ${AlignDir}/*_nt.fasta); do
OutFile=$(basename $file | sed 's@_nt.fasta@_nt_aligned.fasta@g')
Jobs=$(squeue -u did23faz| grep 'mafft'  | wc -l)
while [ $Jobs -gt 150 ]; do
    sleep 300s
    printf "."
    Jobs=$(squeue -u did23faz| grep 'mafft'| wc -l)
done
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leo_busco_nt
ProgDir=~/git_repos/Wrappers/NBI
echo "$file" >> mafft_log.txt
sbatch $ProgDir/sub_mafft_alignment.sh $file $OutDir $OutFile 2>&1 >> mafft_log.txt
done

for gene in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leo_busco_nt/*_aligned.fasta); do
  ID=$(basename $gene |sed 's@_nt_aligned.fasta@@g')
  echo $ID
  mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leo_busco_nt/$ID
  cp $gene /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leo_busco_nt/$ID
done
```
#### TRIMAL

```bash
#Trim the alignments:
for Alignment in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leo_busco_nt/*/*_aligned.fasta); do
  OutDir=$(dirname $Alignment)
  TrimmedName=$(basename $Alignment .fasta)"_trimmed.fasta"
  echo $Alignment
  echo $OutDir
  echo $TrimmedName
  sed -i '/^>/! s/[yrmwn]/N/g' "$Alignment" #trimal does not like all IUPAC symbols, just atgcn
  singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/trimal1.4.1.sif trimal -in $Alignment -out $OutDir/$TrimmedName -keepheader -automated1
done

#Trim header names as RAxML need <60 characters in length:
for Alignment in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leo_busco_nt/*/*_aligned_trimmed.fasta); do
New=$(dirname $Alignment)/$(basename $Alignment .fasta)_edit.fasta
cat $Alignment  | cut -f1 -d '|'  > $New
done

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/leo_busco_nt/*/*_edit.fasta); do
while IFS= read -r line; do
    if [[ "$line" =~ ^\>.+ ]]; then
        echo "$line" | wc -c
    fi
done < $file
done
```
#### IQTREE
```bash
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/busco_nt/AlignDir-100
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/busco_nt/AlignDir/ | shuf -n 100 | xargs -I {} ln -s /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/busco_nt/AlignDir/{} /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/busco_nt/AlignDir-100/{}
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/busco_nt/iqtree2
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/busco_nt/iqtree2
AlignDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/busco_nt/AlignDir
cpu=32
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/iqtree_2.3.0.sif iqtree2 -s $AlignDir -m MFA -T AUTO --threads-max $cpu
#59524339, 59531606 (all)
#Best-fit model: SYM+I+R9 chosen according to BIC, with 100 only SYM+I+R7 was chosen

mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/busco_nt/AlignDir-6
for gene in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/leo_busco_nt/*_nt_aligned_trimmed_edit.fasta); do
ln -s $gene /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/busco_nt/AlignDir-6/.
done
mkdir /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/busco_nt/iqtree2-6
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/busco_nt/iqtree2-6
AlignDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/busco_nt/AlignDir-6
cpu=32
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/iqtree_2.3.0.sif iqtree2 -s $AlignDir -m SYM+I+R9 -B 1000 -T AUTO --threads-max $cpu
#2042808
```
#### RAXML
```bash
#Run RAxML
#In practice, one should use specialized tools such as ModelTest-NG, IQTree/ModelFinder or PartitionFinder for model selection
conda activate raxml
for Alignment in $(ls /home/theaven/scratch/uncompressed/analysis/phylogeny/busco/leo_busco_nt/*/*_edit.fasta); do
Prefix=$(basename $Alignment | cut -f1 -d '_')
OutDir=/home/theaven/scratch/uncompressed/analysis/phylogeny/busco/mildew_leo_busco_nt/RAxML/$Prefix
ProgDir=/home/theaven/scratch/apps/phylogeny
mkdir -p $OutDir
if [ ! -f ${OutDir}/${Prefix}.raxml.bestTree ]; then
Jobs=$(squeue -u theaven| grep 'RAxML'  | wc -l)
while [ $Jobs -gt 16 ]; do
    sleep 300s
    printf "."
    Jobs=$(squeue -u theaven| grep 'RAxML'| wc -l)
done
model=SYM+I+R9
sbatch $ProgDir/run_RAxML_msa.sh $Alignment $OutDir $Prefix $model 2>&1 >> raxmllog.txt
echo $Prefix >> raxmllog.txt
tail -n 2 raxmllog.txt
else
echo $Prefix already done
fi
done
conda deactivate
```
#### ASTRAL
```bash
#Combine individual gene trees into a consensus tree:
srun -p medium --mem 64G --cpus-per-task 32 --pty bash
conda activate astral
ProgDir=../apps/astral/Astral/
cat /home/theaven/scratch/uncompressed/analysis/phylogeny/busco/mildew_leo_busco_nt/RAxML/*/*.raxml.bestTree > /home/theaven/scratch/uncompressed/analysis/phylogeny/tree-files.txt
ls /home/theaven/scratch/uncompressed/analysis/phylogeny/busco/mildew_leo_busco_nt/RAxML/*/*.raxml.bootstraps > /home/theaven/scratch/uncompressed/analysis/phylogeny/bs-files.txt
java -Xmx60000M -jar $ProgDir/astral.5.7.8.jar -i /home/theaven/scratch/uncompressed/analysis/phylogeny/tree-files.txt -o /home/theaven/scratch/uncompressed/analysis/phylogeny/leotiomycetes_phylogeny+++.astral.tre
java -Xmx60000M -jar $ProgDir/astral.5.7.8.jar -t 4 -i /home/theaven/scratch/uncompressed/analysis/phylogeny/tree-files.txt -b /home/theaven/scratch/uncompressed/analysis/phylogeny/bs-files.txt -r 100 -o /home/theaven/scratch/uncompressed/analysis/phylogeny/leotiomycetes_phylogeny+++.bootstrapped.astral.tre-4

tail -n 1 /home/theaven/scratch/uncompressed/analysis/phylogeny/leotiomycetes_phylogeny+++.bootstrapped.astral.tre-2 > /home/theaven/scratch/uncompressed/analysis/phylogeny/leotiomycetes_phylogeny+++.bootstrapped.consensus2.astral.tre
head -n 100 /home/theaven/scratch/uncompressed/analysis/phylogeny/leotiomycetes_phylogeny+++.bootstrapped.astral.tre-2 > /home/theaven/scratch/uncompressed/analysis/phylogeny/leotiomycetes_phylogeny+++.bootstraps.astral.tre
java -Xmx60000M -jar $ProgDir/astral.5.7.8.jar -t 32 -q /home/theaven/scratch/uncompressed/analysis/phylogeny/leotiomycetes_phylogeny+++.bootstrapped.consensus2.astral.tre -i /home/theaven/scratch/uncompressed/analysis/phylogeny/leotiomycetes_phylogeny+++.bootstraps.astral.tre -o /home/theaven/scratch/uncompressed/analysis/phylogeny/leotiomycetes_phylogeny+++.bootstrapped.scored2.astral.tre 
java -Xmx60000M -jar $ProgDir/astral.5.7.8.jar -t 32 -i /home/theaven/scratch/uncompressed/analysis/phylogeny/tree-files.txt -q /home/theaven/scratch/uncompressed/analysis/phylogeny/leotiomycetes_phylogeny+++.bootstrapped.consensus2.astral.tre -o /home/theaven/scratch/uncompressed/analysis/phylogeny/leotiomycetes_phylogeny+++.bootstrapped.localposteriorprobabilityscored2.astral.tre 2> astral-scored.log
```
#### R8s
The ASTRAL tree could not be made ultrametric as all terminal branch length are zero, the iqtree was used instead.

Making the species tree ultrametric CAFE requires a tree that is ultramatric. There are many ways to obtain ultrametric trees (also known as timetrees, these are phylogenetic trees scaled to time, where all paths from root to tips have the same length). Here, we use a fast program called r8s. You will need to know the number of sites in the alignment used to estimate the species tree (the one you want to make ultrametric), and then you can specify one or more calibration points (ideally, the age or age window of a documented fossil) to scale branch lengths into time units. 

```bash
singularity exec /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/containers/python3.sif python3 ~/git_repos/Scripts/NBI/collpse_newick_single_node_leaves.py /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/trimmed-6.txt /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/ortho/analysis/cafe/trimmed-6.tree
nano /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/analysis/phylogeny/busco_nt/iqtree2-6/AlignDir-6.iqtree #3926774 sites
source package /nbi/software/testing/bin/r8s-1.80
echo "#NEXUS" > r8s_ctl_file-m.txt
echo "begin trees;" >> r8s_ctl_file-m.txt
echo "tree erysiphales_tree = [&R] $(cat /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/ortho/analysis/cafe/trimmed-6.tree)" >> r8s_ctl_file-m.txt
echo "end;" >> r8s_ctl_file-m.txt
echo "begin rates;" >> r8s_ctl_file-m.txt
echo "blformat nsites=3926774 lengths=persite ultrametric=no;" >> r8s_ctl_file-m.txt
echo "collapse;" >> r8s_ctl_file-m.txt
echo "mrca callibration_node_1 BLU_HOR_GCA_900237765.1 POD_APH_DRCT72020;" >> r8s_ctl_file-m.txt
echo "mrca callibration_node_2 BLU_HOR_GCA_900237765.1 PAR_POL_PARPO1;" >> r8s_ctl_file-m.txt
echo "mrca callibration_node_3 BLU_HOR_GCA_900237765.1 BLU_HOR_GCA_900519115.1;" >> r8s_ctl_file-m.txt
echo "mrca callibration_node_4 ERY_NEC_GCA_000798795.1 GOL_CIC_GCA_003611215.1;" >> r8s_ctl_file-m.txt
echo "mrca callibration_node_5 ERY_NEC_GCA_000798795.1 BLU_HOR_GCA_900237765.1;" >> r8s_ctl_file-m.txt
echo "constrain taxon=callibration_node_1 min_age=75.0 max_age=75.0;" >> r8s_ctl_file-m.txt
echo "constrain taxon=callibration_node_2 min_age=79.0 max_age=79.0;" >> r8s_ctl_file-m.txt
echo "constrain taxon=callibration_node_3 min_age=0.01 max_age=11.0;" >> r8s_ctl_file-m.txt
echo "constrain taxon=callibration_node_4 min_age=57.0 max_age=57.0;" >> r8s_ctl_file-m.txt
echo "constrain taxon=callibration_node_5 min_age=53.2 max_age=75.0;" >> r8s_ctl_file-m.txt
echo "divtime method=lf algorithm=tn cvStart=0 cvInc=0.5 cvNum=8 crossv=yes;" >> r8s_ctl_file-m.txt
echo "describe plot=chronogram;" >> r8s_ctl_file-m.txt
echo "describe plot=tree_description;" >> r8s_ctl_file-m.txt
echo "end;" >> r8s_ctl_file-m.txt

r8s -b -f r8s_ctl_file-m.txt > temp_r8s-m.txt
tail -n 1 temp_r8s-m.txt > /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/ortho/analysis/cafe/trimmed-6-r8s.tree
```

















































## CSEP Prediction from in house genome assemblies
Duplicate gene predictions have been spotted in our gene files, these genes will need to be removed from future analysis
```bash
for Logfile in $(ls gene_pred/*/*/codingquarry/rep_modeling/final/final_genes_appended_renamed.log); do
Output=$(dirname $Logfile)/$(echo $Logfile |cut -f3 -d '/')_duplicated_genes
echo $Output
rm $Output
for prediction in $(cat $Logfile |cut -f1|sort|uniq -c|sort -nr | awk '$1>1' | awk '{print $2}'); do
    printf ">\n" >> $Output
    grep -w $prediction $Logfile | awk '{print $3}' >> $Output
done
sed -e '/>/,+1d' ${Output} > ${Output}.txt #removed one of each duplicate from the removal list
done

scp gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/THeavenDRCT72020_1_duplicated_genes.txt theaven@gruffalo.cropdiversity.ac.uk:/home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/THeavenDRCT72020_1_duplicated_genes.txt
scp gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/THeavenDRCT72021_1_duplicated_genes.txt theaven@gruffalo.cropdiversity.ac.uk:/home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/THeavenDRCT72021_1_duplicated_genes.txt
scp gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/THeavenSCOTT2020_1_duplicated_genes.txt theaven@gruffalo.cropdiversity.ac.uk:/home/theaven/scratch/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/THeavenSCOTT2020_1_duplicated_genes.txt
scp gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/THeavenp11_1_duplicated_genes.txt theaven@gruffalo.cropdiversity.ac.uk:/home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/THeavenp11_1_duplicated_genes.txt
scp gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/THeavenpOGB2019_1_duplicated_genes.txt theaven@gruffalo.cropdiversity.ac.uk:/home/theaven/scratch/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/THeavenpOGB2019_1_duplicated_genes.txt
scp gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/THeavenpOGB2021_1_duplicated_genes.txt theaven@gruffalo.cropdiversity.ac.uk:/home/theaven/scratch/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/THeavenpOGB2021_1_duplicated_genes.txt
```
#### Predector
```bash
scp /projects/nano_diagnostics/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta theaven@gruffalo.cropdiversity.ac.uk:/home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/
scp /projects/nano_diagnostics/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta theaven@gruffalo.cropdiversity.ac.uk:/home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/
scp /projects/nano_diagnostics/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta theaven@gruffalo.cropdiversity.ac.uk:/home/theaven/scratch/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/
scp /projects/nano_diagnostics/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta theaven@gruffalo.cropdiversity.ac.uk:/home/theaven/scratch/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/
scp /projects/nano_diagnostics/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta theaven@gruffalo.cropdiversity.ac.uk:/home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/
scp /projects/nano_diagnostics/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta theaven@gruffalo.cropdiversity.ac.uk:/home/theaven/scratch/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/

conda activate predector
for proteome in $(ls /home/theaven/scratch/gene_pred/*/*/codingquarry/rep_modeling/final/*.fasta); do
ProgDir=/home/theaven/scratch/apps/predector
OutDir=$(echo $proteome|cut -f1,2,3,4,5,6,7 -d '/')/predector_singularity2
sbatch $ProgDir/predector_singularity.sh $proteome 1.2.6
sleep 30
done
conda deactivate
#2860718,21,25,26,27,28, 2864314-19
#3038774-9
#3038780-5,6,7

for file in $(ls gene_pred/*/*/predector_singularity2/results/final_genes_appended_renamed.pep/final_genes_appended_renamed.pep-ranked.tsv); do
Organism=$(echo $file|cut -f2 -d '/')
Strain=$(echo $file|cut -f3 -d '/')
OutDir=$(dirname $file)
echo ${Organism}-${Strain}
cp $file ${OutDir}/${Organism}-${Strain}-ranked.tsv
done

#Remove discovered duplicates:
for file in $(ls /home/theaven/scratch/gene_pred/*/*/codingquarry/rep_modeling/final/*_duplicated_genes.txt); do
  Organism=$(echo $file|cut -f6 -d '/')
  Strain=$(echo $file|cut -f7 -d '/')
  tsv=$(echo $file | cut -f1,2,3,4,5,6,7 -d '/')/predector_singularity2/results/final_genes_appended_renamed.pep/${Organism}-${Strain}-ranked.tsv
  echo $tsv
  for var1 in $(cat $file); do
  sed -i "/$var1/d" "$tsv"
done
done
```

#### Looking for transcription factors
```bash
conda activate predector2.7
  for Interpro in $(ls gene_pred/*/*/interproscan/NRI/*_interproscan.tsv); do
    Organism=$(echo $Interpro | rev | cut -f5 -d '/' | rev)
    Strain=$(echo $Interpro | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/$Organism/$Strain/transcription_factors
    mkdir -p $OutDir
    ProgDir=/home/theaven/scratch/apps/tools
    python2.7 $ProgDir/interpro2TFs.py --InterPro $Interpro > $OutDir/"$Strain"_TF_domains.tsv
        for var1 in $(cat /home/theaven/scratch/gene_pred/${Organism}/${Strain}/codingquarry/rep_modeling/final/${Strain}_duplicated_genes.txt); do
        sed -i "/$var1/d" "$OutDir/"$Strain"_TF_domains.tsv"
        done #remove discovered duplicates
    echo "$Organism - $Strain " >> TFreport.txt
    printf "total number of transcription factors: \t" >> TFreport.txt
    cat $OutDir/"$Strain"_TF_domains.tsv | cut -f1 | sort | uniq > $OutDir/"$Strain"_TF_gene_headers.txt
    cat $OutDir/"$Strain"_TF_gene_headers.txt | wc -l >> TFreport.txt 
    # Gene ID rather than transcript ID
    cat $OutDir/"$Strain"_TF_gene_headers.txt | sed -e "s/.t.*//g" > $OutDir/"$Strain"_TF_geneid_headers.txt
  done
conda deactivate
mv *report.txt reports/.
```
P_aphanis - THeavenDRCT72020_1
total number of transcription factors:  1162
P_aphanis - THeavenDRCT72021_1
total number of transcription factors:  1034
P_aphanis - THeavenSCOTT2020_1
total number of transcription factors:  1035
P_leucotricha - THeavenp11_1
total number of transcription factors:  994
P_leucotricha - THeavenpOGB2019_1
total number of transcription factors:  1015
P_leucotricha - THeavenpOGB2021_1
total number of transcription factors:  991

#### Sectretome
Proteins with a predicted signal peptide and no predicted transmembrane domains were extracted.
```bash
conda activate predector2.7
for file in $(ls gene_pred/*/*/predector_singularity2/results/final_genes_appended_renamed.pep/final_genes_appended_renamed.pep-ranked.tsv); do
  Organism=$(echo $file|cut -f2 -d '/')
  Strain=$(echo $file|cut -f3 -d '/')
  echo "$Organism - $Strain"
  Proteome=$(ls /home/theaven/scratch/gene_pred/${Organism}/${Strain}/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta)
  Sequential_Proteome=$(dirname $Proteome)/${Strain}.pep.fasta
  Outfile=$(dirname $file)/${Strain}_secretome.tsv
  Outfile2=$(dirname $file)/${Strain}_secretome.pep.fasta
  rm $Outfile
  rm $Outfile2
  head -n 1 $file >> $Outfile
  awk '$25==1' $file |awk '$27==0' |awk '$28==0' >> $Outfile
#  perl -MBio::SeqIO -e 'my $seqin = Bio::SeqIO->new(-fh => \*STDIN, -format => 'fasta'); while (my $seq = $seqin->next_seq) { print ">",$seq->id,"\n",$seq->seq,"\n"; }' < $Proteome > $Sequential_Proteome
  cat $Outfile | cut -f1 | grep -A 1 -f - $Sequential_Proteome >> $Outfile2
done
conda  deactivate
```

#### SSCP
Small cysteine rich proteins were identified within proteomes. These proteins may be identified by EffectorP, but this approach allows direct control over what constitutes a SSCP.

```bash
conda activate predector2.7
for Proteome in $(ls /home/theaven/scratch/gene_pred/*/*/predector_singularity2/results/final_genes_appended_renamed.pep/*_secretome.pep.fasta); do
Organism=$(echo $Proteome|cut -f6 -d '/')
Strain=$(echo $Proteome|cut -f7 -d '/')
echo "$Organism - $Strain"
OutDir=gene_pred/$Organism/$Strain/sscp
mkdir $OutDir
ProgDir=/home/theaven/scratch/apps/tools
python2.7 $ProgDir/sscp_filter.py --inp_fasta $Proteome --max_length 300 --threshold 3 --out_fasta $OutDir/"$Strain"_sscp_all_results.fa
cat $OutDir/"$Strain"_sscp_all_results.fa | grep 'Yes' > $OutDir/"$Strain"_sscp.fa
echo "$Organism - $Strain " >> sscpreport.txt
printf "number of SSC-rich genes:\t" >> sscpreport.txt
cat $OutDir/"$Strain"_sscp.fa | grep '>' | tr -d '>' | cut -f1 -d '.' | sort | uniq | wc -l >> sscpreport.txt
done
conda deactivate
mv *report.txt reports/.
```
P_aphanis - THeavenDRCT72020_1
number of SSC-rich genes:       52
P_aphanis - THeavenDRCT72021_1
number of SSC-rich genes:       59
P_aphanis - THeavenSCOTT2020_1
number of SSC-rich genes:       44
P_leucotricha - THeavenp11_1
number of SSC-rich genes:       85
P_leucotricha - THeavenpOGB2019_1
number of SSC-rich genes:       89
P_leucotricha - THeavenpOGB2021_1
number of SSC-rich genes:       90

#### CAZY proteins
Carbohydrte active enzymes were idnetified using CAZYfollowing recomendations at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt :
```bash
conda activate cazy
for Proteome in $(ls gene_pred/*/*/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta); do
Organism=$(echo $Proteome|cut -f2 -d '/')
Strain=$(echo $Proteome|cut -f3 -d '/')
echo "$Organism - $Strain"
OutDir=gene_pred/$Organism/$Strain/cazy
mkdir $OutDir
Prefix="$Strain"_CAZY
CazyHmm=apps/tools/dbCAN/dbCAN/dbCAN-fam-HMMs.txt
ProgDir=/home/theaven/scratch/apps/tools
sbatch $ProgDir/hmm_scan.sh $CazyHmm $Proteome $Prefix $OutDir
done
#3054426-31
conda deactivate
mv *report.txt reports/.


```
The Hmm parser was used to filter hits by an E-value of E1x10-5 or E1x10-e3 if they had a hit over a length of X %.

Those proteins with a signal peptide were extracted from the list and gff files representing these proteins made.
```bash
for File in $(ls gene_pred/*/*/cazy/*CAZY.out.dm); do
Organism=$(echo $File|cut -f2 -d '/')
Strain=$(echo $File|cut -f3 -d '/')
echo "$Organism - $Strain"
    for var1 in $(cat /home/theaven/scratch/gene_pred/${Organism}/${Strain}/codingquarry/rep_modeling/final/${Strain}_duplicated_genes.txt); do
    sed -i "/$var1/d" "gene_pred/"${Organism}"/"${Strain}"/cazy/"$Strain"_CAZY.out.dm"  
    done #remove discovered duplicates
OutDir=$(dirname $File)
ProgDir=/home/theaven/scratch/apps/tools/dbCAN/dbCAN
$ProgDir/hmmscan-parser.sh $OutDir/"$Strain"_CAZY.out.dm > $OutDir/"$Strain"_CAZY.out.dm.ps
CazyHeaders=$(echo $File | sed 's/.out.dm/_headers.txt/g')
rm $OutDir/"$Strain"_CAZY.out.dm.ps.f
awk '$5<0.00000000000000001' $OutDir/"$Strain"_CAZY.out.dm.ps |awk '$10>0.45' >> $OutDir/"$Strain"_CAZY.out.dm.ps.f #For fungi, use E-value < 1e-17 and coverage > 0.45
cat $OutDir/"$Strain"_CAZY.out.dm.ps.f | cut -f3 | sort | uniq > $CazyHeaders
echo $CazyHeaders
head $CazyHeaders
Gff=$(ls gene_pred/${Organism}/${Strain}/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3)
CazyGff=${OutDir}/${Strain}_CAZY.gff
ProgDir=/home/theaven/scratch/apps/tools
$ProgDir/extract_gff_for_sigP_hits.pl $CazyHeaders $Gff CAZyme ID > $CazyGff
echo $CazyGff
head $CazyGff
SecretedProts=gene_pred/${Organism}/${Strain}/predector_singularity2/results/final_genes_appended_renamed.pep/*_secretome.pep.fasta
SecretedHeaders=$(echo $SecretedProts | sed 's/.pep.fasta/_headers.txt/g')
cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
CazyGffSecreted=$OutDir/${Strain}_CAZY_secreted.gff
$ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $CazyGff Secreted_CAZyme ID > $CazyGffSecreted
echo $CazyGffSecreted
head $CazyGffSecreted
cat $CazyGffSecreted | grep -w 'mRNA' | cut -f9 | tr -d 'ID=' | cut -f1 -d ';' > $OutDir/"$Strain"_CAZY_secreted_headers.txt

echo "$Organism - $Strain" >> cazyreport.txt
echo "number of CAZY proteins identified:" >> cazyreport.txt
cat $CazyHeaders | wc -l >> cazyreport.txt
echo "number of CAZY genes identified:" >> cazyreport.txt
cat $CazyGff | grep -w 'gene' | wc -l >> cazyreport.txt
echo "number of Secreted CAZY proteins identified:">> cazyreport.txt
cat $OutDir/"$Strain"_CAZY_secreted_headers.txt | wc -l >> cazyreport.txt
echo "number of Secreted CAZY genes identified:">> cazyreport.txt
cat $CazyGffSecreted | grep -w 'gene' | wc -l >> cazyreport.txt
done
```

P_aphanis - THeavenDRCT72020_1
number of CAZY proteins identified:
138
number of CAZY genes identified:
138
number of Secreted CAZY proteins identified:
58
number of Secreted CAZY genes identified:
58
P_aphanis - THeavenDRCT72021_1
number of CAZY proteins identified:
190
number of CAZY genes identified:
190
number of Secreted CAZY proteins identified:
89
number of Secreted CAZY genes identified:
89
P_aphanis - THeavenSCOTT2020_1
number of CAZY proteins identified:
139
number of CAZY genes identified:
139
number of Secreted CAZY proteins identified:
58
number of Secreted CAZY genes identified:
58
P_leucotricha - THeavenp11_1
number of CAZY proteins identified:
144
number of CAZY genes identified:
144
number of Secreted CAZY proteins identified:
57
number of Secreted CAZY genes identified:
57
P_leucotricha - THeavenpOGB2019_1
number of CAZY proteins identified:
170
number of CAZY genes identified:
170
number of Secreted CAZY proteins identified:
66
number of Secreted CAZY genes identified:
66
P_leucotricha - THeavenpOGB2021_1
number of CAZY proteins identified:
145
number of CAZY genes identified:
145
number of Secreted CAZY proteins identified:
58
number of Secreted CAZY genes identified:
58

```bash
conda activate predector2.7
for File in $(ls gene_pred/*/*/cazy/*_CAZY.out.dm.ps.f); do
Organism=$(echo $File|cut -f2 -d '/')
Strain=$(echo $File|cut -f3 -d '/')
echo "$Organism - $Strain"
OutDir=$(dirname $File)
Secreted=$(ls gene_pred/${Organism}/${Strain}/predector_singularity2/results/final_genes_appended_renamed.pep/*_secretome_headers.txt)
Gff=$(ls gene_pred/${Organism}/${Strain}/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3)
ProgDir=/home/theaven/scratch/apps/tools/CAZY
python2.7 $ProgDir/summarise_CAZY.py --cazy $File --inp_secreted $Secreted --inp_gff $Gff --summarise_family --trim_gene_id 2 --kubicek_2014
done
conda deactivate
```
P_aphanis - THeavenDRCT72020_1
A-Galactosidases - 1
Cellulases - 2
other - 54
Polygalacturonase - 1

P_aphanis - THeavenDRCT72021_1
B-Galactosidases - 1
Polygalacturonase - 1
A-Arabinosidases - 2
Xylanases - 2
Polygalacturonate lyases - 2
A-Galactosidases - 1
B-Glycosidases - 1
Cellulases - 4
other - 75

P_aphanis - THeavenSCOTT2020_1
A-Galactosidases - 1
Cellulases - 2
other - 54
Polygalacturonase - 1

P_leucotricha - THeavenp11_1
A-Galactosidases - 1
Cellulases - 2
other - 53
Polygalacturonase - 1

P_leucotricha - THeavenpOGB2019_1
A-Galactosidases - 1
Cellulases - 2
other - 62
Polygalacturonase - 1

P_leucotricha - THeavenpOGB2021_1
A-Galactosidases - 1
Cellulases - 2
other - 54
Polygalacturonase - 1

#### BLAST - CSEPs and Mildew orthologues
```bash
conda activate blast+
for proteome in $(ls /home/theaven/scratch/gene_pred/P*/T*/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta); do 
OutDir=$(dirname $proteome)
Database=db/CSEPdb
ProgDir=/home/theaven/scratch/apps/blast
sbatch $ProgDir/blastp.sh $proteome $Database $OutDir/outblast.tsv 2>&1 >> blastoutlog.txt
Database=db2/blumeriadb
sbatch $ProgDir/blastp.sh $proteome $Database $OutDir/inblast.tsv 2>&1 >> blastinlog.txt
done
```
#### SecretomeP
Abinitio naming schema is wrong, these are from RNA evidences protein predictions
```bash
conda activate jellyfish
#linearize and edit 
for proteome in $(ls gene_pred/P*/T*/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta); do
  Outfile=$(echo $proteome | cut -d '/' -f3 | cut -d '.' -f1)_
  OutDir=$(dirname $proteome)
  ProgDir=/home/theaven/scratch/apps/secretomep
  echo ${OutDir}/${Outfile}abinitio_genes_40.fsa
  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $proteome | sed "s@file_1_@@g" | sed -e "s@-@N@g" | sed -e "s@X@@g" | awk '/^>/{val=$0;next} length($0)>40{print val ORS $0}' > ${OutDir}/${Outfile}abinitio_genes_40.fsa
  sbatch $ProgDir/secretomep.sh ${OutDir}/${Outfile}abinitio_genes_40.fsa ${OutDir} ${Outfile}abinitio_genes_40_secretomep.out 2>&1 >> abinitiosecretomeplog.txt
done

for ome in $(ls /home/theaven/scratch/gene_pred/P*/T*/codingquarry/rep_modeling/final/*_genes_40_secretomep.out); do 
sed -i -n '/SecretomeP 1.0g predictions/,$p' $ome 
sed -i '$d' $ome 
sed -i '1d' $ome 
sed -i '1d' $ome 
sed -i '1d' $ome 
sed -i '1d' $ome 

for proteome in $(ls gene_pred/P*/T*/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta); do 
 ID=$( echo $proteome | cut -d '/' -f3)
 Annotations=$(ls gene_pred/P*/${ID}/predector_singularity3/results/final_genes_appended_renamed.pep/P*_1-ranked.tsv)
 echo $Annotations
 blast=$(dirname $proteome)
 OutDir=$(echo gene_pred/P*/${ID}/predector_singularity3)
 echo $OutDir
 cat ${blast}/inblast.tsv | cut -f1 |sort | uniq > ${blast}/abinitio_ingenes.tsv
 cat ${blast}/outblast.tsv | cut -f1 |sort | uniq > ${blast}/abinitio_outgenes.tsv
 for f in ${blast}/abinitio_ingenes.tsv; do sed -i "s/$/\t1/" $f; done
 for f in ${blast}/abinitio_outgenes.tsv; do sed -i "s/$/\t1/" $f; done
 printf 'name\tmildew orthology match' >> ${blast}/abinitio_ingenes.tsv
 printf 'name\toutgroup orthology match' >> ${blast}/abinitio_outgenes.tsv
 join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84 <(sort $Annotations) <(sort ${blast}/abinitio_outgenes.tsv) > ${OutDir}/${ID}_abinitio_annotations.tsv
 join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84,1.85 <(sort ${OutDir}/${ID}_abinitio_annotations.tsv) <(sort ${blast}/abinitio_ingenes.tsv) > ${OutDir}/${ID}_abinitio_annotations2.tsv
 Secretomep=$(dirname $proteome)/$(echo $proteome | cut -d '/' -f3 | cut -d '.' -f1)_abinitio_genes_40_secretomep
 ls ${Secretomep}.out
 printf 'name\tNN-score\tOdds\tWeighted\tWarning' >> ${Secretomep}.out
 awk -v OFS='\t' '{$1=$1}1' ${Secretomep}.out > ${Secretomep}.tsv
 join -t $'\t' -a 1 -e 0 -o 0,2.2,2.3,2.4,2.5,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84,1.85,1.86 <(sort ${OutDir}/${ID}_abinitio_annotations2.tsv) <(sort ${Secretomep}.tsv) > ${OutDir}/${ID}_abinitio_annotations3.tsv
done
```
Unfortunately, whilst abinitio_ingenes.tsv etc. are in /home/theaven/projects/niab/theaven/gene_pred///codingquarry/rep_modeling/final/ and _1-ranked.tsv are in /home/theaven/projects/niab/theaven/gene_pred//predector_singularity3/results/final_genes_appended_renamed.pep/ the combined files seem to have been deleted.

#### Report Results
```bash
 echo Species+Assembly Predicted_genes pfamID pfam_virulence PHIbase_effector effector_ortholog canonical_SP noncanonical_SP SSCP CSEP Secreted_effectorP3 Secreted_CAZY newSSCP newCSEP newSecreted_effectorP3 > all6EffectorPredictionReport_22092024.txt
for proteome in $(ls /home/theaven/projects/niab/theaven/gene_pred/P*/THeaven*/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta); do 
 ID=$( echo $proteome | cut -d '/' -f9)
 OutDir=$(dirname $proteome)
 Annotations=${OutDir}/${ID}_abinitio_annotations3.tsv
 echo $Annotations
 Predictedgenes=$(cut -f1 $Annotations| grep -v 'name' | wc -l)
 pfamID=$(awk '$18 != "." {print $1}' $Annotations| grep -v 'name' | wc -l)
 pfamvirulence=$(awk '$20 == "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 PHIbaseeffector=$(awk '$15 == "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 effectorortholog=$(awk '$11 != "." {print $1}' $Annotations| grep -v 'name' | wc -l)
 canonicalsp=$(awk '$31 == "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 noncanonical=$(awk '$31 != "1" && $5 != 'signal' && $2 > 0.6 {print $1}' $Annotations| grep -v 'name' | wc -l)
 sscp=$(awk '$31 == "1" && $36 <= 300 && $39 > 3 {print $1}' $Annotations| grep -v 'name' | wc -l)
 csep=$(awk '$31 == "1" && $7 == "0" {print $1}' $Annotations| grep -v 'name' | wc -l)
 secretedeffector=$(awk '$31 == "1" && $27 == "." && $7 == "0" {print $1}' $Annotations| grep -v 'name' | wc -l)
 secretedcazy=$(awk '$31 == "1" && $21 != "." {print $1}' $Annotations| grep -v 'name' | wc -l)
 newsscp=$(awk '$31 == "1" && $36 <= 300 && $39 > 3 && $11 == "." && $15 != "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 newcsep=$(awk '$31 == "1" && $7 == "0" && $11 == "." && $15 != "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 newsecretedeffector=$(awk '$31 == "1" && $27 == "." && $11 == "." && $15 != "1" && $7 == "0" {print $1}'  $Annotations| grep -v 'name' | wc -l)
 echo $ID $Predictedgenes $pfamID $pfamvirulence $PHIbaseeffector $effectorortholog $canonicalsp $noncanonical $sscp $csep $secretedeffector $secretedcazy $newsscp $newcsep $newsecretedeffector >> all6EffectorPredictionReport_22092024.txt
  ID=$( echo $proteome | cut -d '/' -f9)
 OutDir=$(dirname $proteome)
 Annotations=${OutDir}/${ID}_abinitio_annotations3.tsv
 echo $Annotations
 Predictedgenes=$(awk '$1 ~ /\.t1$/ {print $1}' $Annotations| grep -v 'name' | wc -l)
 pfamID=$(awk '$18 != "." && $1 ~ /\.t1$/ {print $1}' $Annotations| grep -v 'name' | wc -l)
 pfamvirulence=$(awk '$20 == "1" && $1 ~ /\.t1$/ {print $1}' $Annotations| grep -v 'name' | wc -l)
 PHIbaseeffector=$(awk '$15 == "1" && $1 ~ /\.t1$/ {print $1}' $Annotations| grep -v 'name' | wc -l)
 effectorortholog=$(awk '$11 != "." && $1 ~ /\.t1$/ {print $1}' $Annotations| grep -v 'name' | wc -l)
 canonicalsp=$(awk '$31 == "1" && $1 ~ /\.t1$/ {print $1}' $Annotations| grep -v 'name' | wc -l)
 noncanonical=$(awk '$31 != "1" && $5 != 'signal' && $2 > 0.6 && $1 ~ /\.t1$/ {print $1}' $Annotations| grep -v 'name' | wc -l)
 sscp=$(awk '$31 == "1" && $36 <= 300 && $39 > 3 && $1 ~ /\.t1$/ {print $1}' $Annotations| grep -v 'name' | wc -l)
 csep=$(awk '$31 == "1" && $7 == "0" && $1 ~ /\.t1$/ {print $1}' $Annotations| grep -v 'name' | wc -l)
 secretedeffector=$(awk '$31 == "1" && $27 == "." && $7 == "0" && $1 ~ /\.t1$/ {print $1}' $Annotations| grep -v 'name' | wc -l)
 secretedcazy=$(awk '$31 == "1" && $21 != "." && $1 ~ /\.t1$/ {print $1}' $Annotations| grep -v 'name' | wc -l)
 newsscp=$(awk '$31 == "1" && $36 <= 300 && $39 > 3 && $11 == "." && $15 != "1" && $1 ~ /\.t1$/ {print $1}' $Annotations| grep -v 'name' | wc -l)
 newcsep=$(awk '$31 == "1" && $7 == "0" && $11 == "." && $15 != "1" && $1 ~ /\.t1$/ {print $1}' $Annotations| grep -v 'name' | wc -l)
 newsecretedeffector=$(awk '$31 == "1" && $27 == "." && $11 == "." && $15 != "1" && $7 == "0" && $1 ~ /\.t1$/ {print $1}'  $Annotations| grep -v 'name' | wc -l)
 echo $ID $Predictedgenes $pfamID $pfamvirulence $PHIbaseeffector $effectorortholog $canonicalsp $noncanonical $sscp $csep $secretedeffector $secretedcazy $newsscp $newcsep $newsecretedeffector >> all6EffectorPredictionReport_22092024.txt
done
```
Report files are in /home/theaven/projects/niab/theaven/reports/


Identifying mating types
```bash
#P.leucotricha MAT primers from ganan et al.
for assembly in $(ls mildews/HEA*.fna); do
grep -B 1 'CGTCCCAACCAAGGTTTAAG\|TGCTGCTGAAATTATTGACCA' $assembly > MAT1-1-1forward.txt
grep -B 1 'TGGTCAATAATTTCAGCAGCA\|CTTAAACCTTGGTTGGGACG' $assembly > MAT1-1-1reverse.txt
grep '>' MAT1-1-1forward.txt > temp.txt
echo $assembly >> MAT1-1-1.txt
grep -f temp.txt MAT1-1-1reverse.txt >> MAT1-1-1.txt
done

for assembly in $(ls mildews/HEA*.fna); do
grep -B 1 'CCCTCCAAACTCGTGGATTT\|AAATCCACGAGTTTGGAGGG' $assembly > MAT1-2-1forward.txt
grep -B 1 'CCCGCATACGGTTTCTCTTT\|AAAGAGAAACCGTATGCGGG' $assembly > MAT1-2-1reverse.txt
grep '>' MAT1-2-1forward.txt > temp.txt
echo $assembly >> MAT1-2-1.txt
grep -f temp.txt MAT1-2-1reverse.txt >> MAT1-2-1.txt
done

#DRCT72021:
#DRCT72020:
#SCOTT2020:
#OGB2021:contig_1528 MAT1-2 CCCGCATACGGTTTCTCTTTTGACGCCGCTTCACTTCCGATGATTTTCGTGGTTTATAACGGTATCCAGGATTTTGAACAGCATGTTGACGCTTTGCACGATCTGCGAGATCCTTGTACGTATTTCGAACTTCGTTTGGCTCTCTAGCCCACATTGTAGAGATTAAGGTTGCTAGAAAATACATATGAGAATAAGAATCATACTAGATAGACAAGCTTACATATTTGGTTATTACTTAGACCAGGGTTCTGAGCCACAGTTTCTTCATGTTTAGCCTTTCGATATAAAATCCACGAGTTTGGAGGG
#P112020:contig_1488 MAT1-2 CCCTCCAAACTCGTGGATTTTATATCGAAAGGCTAAACATGAAGAAACTGTGGCTCAGAACCCTGGTCTAAGTAATAACCAAATATGTAAGCTTGTCTATCTAGTATGATTCTTATTCTCATATGTATTTTCTAGCAACCTTAATCTCTACAATGTGGGCTAGAGAGCCAAACGAAGTTCGAAATACGTACAAGGATCTCGCAGATCGTGCAAAGCGTCAACATGCTGTTCAAAATCCTGGATACCGTTATAAACCACGAAAATCATCGGAAGTGAAGCGGCGTCAAAAGAGAAACCGTATGCGGG
#OGB2019:contig_1598 MAT1-2 CCCTCCAAACTCGTGGATTTTATATCGAAAGGCTAAACATGAAGAAACTGTGGCTCAGAACCCTGGTCTAAGTAATAACCAAATATGTAAGCTTGTCTATCTAGTATGATTCTTATTCTCATATGTATTTTCTAGCAACCTTAATCTCTACAATGTGGGCTAGAGAGCCAAACGAAGTTCGAAATACGTACAAGGATCTCGCAGATCGTGCAAAGCGTCAACATGCTGTTCAAAATCCTGGATACCGTTATAAACCACGAAAATCATCGGAAGTGAAGCGGCGTCAAAAGAGAAACCGTATGCGGG

##############################################################################
#P.aphanis MAT primers from asalf et al.
for assembly in $(ls mildews/HEA*.fna); do
grep -B 1 'GGCTTCAGAAGTATGTCTCAATCAA\|TTGATTGAGACATACTTCTGAAGCC' $assembly > aMAT1-1-1forward.txt
grep -B 1 'TTCCACAATTTCGTCAGAAGAG\|CTCTTCTGACGAAATTGTGGAA' $assembly > aMAT1-1-1reverse.txt
grep '>' aMAT1-1-1forward.txt > atemp.txt
echo $assembly >> aMAT1-1-1.txt
grep -f atemp.txt aMAT1-1-1reverse.txt >> aMAT1-1-1.txt
echo ${assembly}_1-1
done

for assembly in $(ls mildews/HEA*.fna); do
grep -B 1 'AAAGGCTAAACACGAAGAAACTGT\|ACAGTTTCTTCGTGTTTAGCCTTT' $assembly > aMAT1-2-1forward.txt
grep -B 1 'TTGTAACGGTAACCAGGATTTTC\|GAAAATCCTGGTTACCGTTACAA' $assembly > aMAT1-2-1reverse.txt
grep '>' aMAT1-2-1forward.txt > atemp.txt
echo $assembly >> aMAT1-2-1.txt
grep -f atemp.txt aMAT1-2-1reverse.txt >> aMAT1-2-1.txt
echo ${assembly}_1-2
done

#DRCT72021:contig_3816 MAT1-2 AAAGGCTAAACACGAAGAAACTGTAGCCCGCAACCCTGGTCTTAGCAACAATCAAATATGTAAGCTTGTTTTTTTAGTATCAATTTTATTCTAATCAAATTCTTCTTTAGCAACCTTGATCTCTTCAATGTGGGCTACAGAGTCGAACGAAGTTCGGAATGTCTACAAGGATCTCGCCGATCGTGCAAAGCGCCAACATGCCGTTGAAAATCCTGGTTACCGTTACAA
#DRCT72020:contig_3763 MAT1-2 TTGTAACGGTAACCAGGATTTTCAACGGCATGTTGGCGCTTTGCACGATCGGCGAGATCCTTGTAGACATTCCGAACTTCGTTCGACTCTGTAGCCCACATTGAAGAGATCAAGGTTGCTAAAGAAGAATTTGATTAGAATAAAATTGATACTAAAAAAACAAGCTTACATATTTGATTGTTGCTAAGACCAGGGTTGCGGGCTACAGTTTCTTCGTGTTTAGCCTTT
#SCOTT2020:contig_4023 MAT1-2 TTGTAACGGTAACCAGGATTTTCAACGGCATGTTGGCGCTTTGCACGATCGGCGAGATCCTTGTAGGCATGTCGAACTTCGTTCGACTCTGTAGCCCACATTGAAGAGATCAAGGTTGCTAAAGAAGAATTTGATTAGAATAAAATTGATACTAAAAAAACAAGCTTACATATTTGATTGTTGCTAAGACCAGGGTTGCGGGCTACAGTTTCTTCGTGTTTAGCCTTT
#OGB2021:
#P112020:
#OGB2019:

#BLASTP with published 1-1 and 1-2 P.leucotricha MAT loci
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' mat.fna > MAT.fna
makeblastdb -in MAT.fna -dbtype nucl -title matdb -out db/matdb
for assembly in $(ls mildews/HEA*.fna); do
  output=$(echo $assembly | cut -d '/' -f2 | cut -d '.' -f1)
  echo $output
blastn -query $assembly -db db/matdb -out ${output}_mat.tsv -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq'
done

#DRCT72021:
#DRCT72020:
#SCOTT2020:
#OGB2021:contig_1528 MAT1-2
#P112020:contig_1488 MAT1-2
#OGB2019:contig_1598 MAT1-2

grep -A 1 '>contig_3763' mildews/HEAVEN_strawberry2020.fna > contig_3763.txt
grep -A 1 '>contig_3816' mildews/HEAVEN_strawberry2021.fna > contig_3816.txt
grep -A 1 '>contig_4023' mildews/HEAVEN_raspberry2020.fna > contig_4023.txt
grep -A 1 '>contig_1598' mildews/HEAVEN_apple2019.fna > contig_1598.txt
grep -A 1 '>contig_1488' mildews/HEAVEN_apple2020.fna > contig_1488.txt
grep -A 1 '>contig_1528' mildews/HEAVEN_apple2021.fna > contig_1528.txt
```
## Comparmentalised genome?

See Frantzeskakis et al. (2019) (Parauncinula) and Armitage et al. (2018).

Gene density plot for Podosphaera - mark BUSCOs and CSEPs + Violin plot of CSEP vs non-CSEP gene to TE distance for Podosphaera.

### Gene-to-TE distances

Müller et al. (2018) found that in the B. graminis genome the upstream regions of short CSEPs were enriched for LTRs of the Copia and Gypsy superfamilies, whereas upstream regions of longer CSEPs were enriched for LINEs and SINEs.

#### Calculte gene-to-TE distances
Put gene and TE annotations into one file and keep only annotations for 'cleaned' contigs:
```bash
awk '$3 == "gene"' /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3 > /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/te.gff
cat /home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/DRT72020/fcs/earlgreyte/Pod_aph_DRT72020/Pod_aph_DRT72020_EarlGrey/Pod___summaryFiles/Pod__.filteredRepeats.gff >> /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/te.gff

awk '$3 == "gene"' /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3 > /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/te.gff  
cat /home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/DRT72021/fcs/earlgreyte/Pod_aph_DRT72021/Pod_aph_DRT72021_EarlGrey/Pod___summaryFiles/Pod__.filteredRepeats.gff >> /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/te.gff 

awk '$3 == "gene"' /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3 > /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/te.gff 
cat /home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/SCOTT2020/fcs/earlgreyte/Pod_aph_SCOTT2020/Pod_aph_SCOTT2020_EarlGrey/Pod___summaryFiles/Pod__.filteredRepeats.gff >> /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/te.gff

awk '$3 == "gene"' /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3 > /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/te.gff 
cat /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGBp112020/fcs/earlgreyte/Pod_leu_OGBp112020/Pod_leu_OGBp112020_EarlGrey/Pod___summaryFiles/Pod__.filteredRepeats.gff >> /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/te.gff

awk '$3 == "gene"' /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3 > /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/te.gff 
cat /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGB2019/fcs/earlgreyte/Pod_leu_OGB2019/Pod_leu_OGB2019_EarlGrey/Pod___summaryFiles/Pod__.filteredRepeats.gff >> /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/te.gff 

awk '$3 == "gene"' /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3 > /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/te.gff
cat /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGB2021/fcs/earlgreyte/Pod_leu_OGB2021/Pod_leu_OGB2021_EarlGrey/Pod___summaryFiles/Pod__.filteredRepeats.gff >> /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/te.gff

sort /home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/DRT72020/fcs/earlgreyte/Pod_aph_DRT72020/Pod_aph_DRT72020_EarlGrey/Pod___summaryFiles/Pod__.filteredRepeats.gff | uniq -c
awk '{count[$1]++} END {for (word in count) print word, count[word]}' /home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/DRT72020/fcs/earlgreyte/Pod_aph_DRT72020/Pod_aph_DRT72020_EarlGrey/Pod___summaryFiles/Pod__.filteredRepeats.gff | wc -l
awk -F'\t' '{split($NF, arr, "|"); for (i in arr) {count[arr[i]]++}} END {for (word in count) print word, count[word]}' OGB2021_sub200-5-GO.tsv
```
Get the gene-to-TE distances for different TE classes
```bash
#Find the 3' and 5' te distances + match names to the predector output files
#script is written for earlgrey output only, see below for runnign of earlgrey on all assemblies
conda activate predector2.7
for Gff in $(ls /home/theaven/projects/niab/theaven/gene_pred/P_*/*/codingquarry/rep_modeling/final/te.gff); do
  ID=$(echo $Gff | cut -d '/' -f9)
OutFile=$(dirname $Gff)/${ID}_flanking_tes.txt
echo -e "ID\tfive_prime_lgth\tfive_prime_family\tfive_prime_ID\tfive_prime_strand\tthree_prime_lgth\tthree_prime_family\tthree_prime_ID\tthree_prime_strand\tstrand" > $OutFile
python2.7 ../apps/tools/find_intergenic_te_regions.py --Gff $Gff >> $OutFile
head -n 1 $OutFile > temp.temp
tail -n +2 $OutFile | awk -F'\t' '{print $1 ".t1\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' >> temp.temp && mv temp.temp $OutFile #NOTE: this line adds .t1 to the gene name, the script itself does not considers different variants, only genes
x=$(grep 'gene' $Gff | wc -l)
y=$(cat $OutFile | wc -l)
echo "Genes in $ID GFF: $x, TEs found for: $y"
done

for Gff in $(ls /home/theaven/projects/niab/theaven/gene_pred/P_*/*/codingquarry/rep_modeling/final/te.gff); do
  ID=$(echo $Gff | cut -d '/' -f9)
OutFile=$(dirname $Gff)/${ID}_flanking_LINEs.txt
echo -e "ID\tfive_prime_lgth\tfive_prime_family\tfive_prime_ID\tfive_prime_strand\tthree_prime_lgth\tthree_prime_family\tthree_prime_ID\tthree_prime_strand\tstrand"  > $OutFile
python2.7 ../apps/tools/find_intergenic_te_distances.py --Gff $Gff --Family LINE >> $OutFile
head -n 1 $OutFile > temp.temp
tail -n +2 $OutFile | awk -F'\t' '{print $1 ".t1\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' >> temp.temp && mv temp.temp $OutFile
x=$(grep 'gene' $Gff | wc -l)
y=$(cat $OutFile | wc -l)
echo "Genes in $ID GFF: $x, distances found for: $y"
done

for Gff in $(ls /home/theaven/projects/niab/theaven/gene_pred/P_*/*/codingquarry/rep_modeling/final/te.gff); do
  ID=$(echo $Gff | cut -d '/' -f9)
OutFile=$(dirname $Gff)/${ID}_flanking_LTRs.txt
echo -e "ID\tfive_prime_lgth\tfive_prime_family\tfive_prime_ID\tfive_prime_strand\tthree_prime_lgth\tthree_prime_family\tthree_prime_ID\tthree_prime_strand\tstrand"  > $OutFile
python2.7 ../apps/tools/find_intergenic_te_distances.py --Gff $Gff --Family LTR >> $OutFile
head -n 1 $OutFile > temp.temp
tail -n +2 $OutFile | awk -F'\t' '{print $1 ".t1\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' >> temp.temp && mv temp.temp $OutFile
x=$(grep 'gene' $Gff | wc -l)
y=$(cat $OutFile | wc -l)
echo "Genes in $ID GFF: $x, distances found for: $y"
done

for Gff in $(ls /home/theaven/projects/niab/theaven/gene_pred/P_*/*/codingquarry/rep_modeling/final/te.gff); do
  ID=$(echo $Gff | cut -d '/' -f9)
OutFile=$(dirname $Gff)/${ID}_flanking_DNAs.txt
echo -e "ID\tfive_prime_lgth\tfive_prime_family\tfive_prime_ID\tfive_prime_strand\tthree_prime_lgth\tthree_prime_family\tthree_prime_ID\tthree_prime_strand\tstrand"  > $OutFile
python2.7 ../apps/tools/find_intergenic_te_distances.py --Gff $Gff --Family DNA >> $OutFile
head -n 1 $OutFile > temp.temp
tail -n +2 $OutFile | awk -F'\t' '{print $1 ".t1\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' >> temp.temp && mv temp.temp $OutFile
x=$(grep 'gene' $Gff | wc -l)
y=$(cat $OutFile | wc -l)
echo "Genes in $ID GFF: $x, distances found for: $y"
done

for Gff in $(ls /home/theaven/projects/niab/theaven/gene_pred/P_*/*/codingquarry/rep_modeling/final/te.gff); do
  ID=$(echo $Gff | cut -d '/' -f9)
OutFile=$(dirname $Gff)/${ID}_flanking_Pogos.txt
echo -e "ID\tfive_prime_lgth\tfive_prime_family\tfive_prime_ID\tfive_prime_strand\tthree_prime_lgth\tthree_prime_family\tthree_prime_ID\tthree_prime_strand\tstrand"  > $OutFile
python2.7 ../apps/tools/find_intergenic_te_distances.py --Gff $Gff --ID DNA/TcMar-Pogo >> $OutFile
head -n 1 $OutFile > temp.temp
tail -n +2 $OutFile | awk -F'\t' '{print $1 ".t1\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' >> temp.temp && mv temp.temp $OutFile
x=$(grep 'gene' $Gff | wc -l)
y=$(cat $OutFile | wc -l)
echo "Genes in $ID GFF: $x, distances found for: $y"
done

for Gff in $(ls /home/theaven/projects/niab/theaven/gene_pred/P_*/*/codingquarry/rep_modeling/final/te.gff); do
  ID=$(echo $Gff | cut -d '/' -f9)
OutFile=$(dirname $Gff)/${ID}_flanking_Mariners.txt
echo -e "ID\tfive_prime_lgth\tfive_prime_family\tfive_prime_ID\tfive_prime_strand\tthree_prime_lgth\tthree_prime_family\tthree_prime_ID\tthree_prime_strand\tstrand"  > $OutFile
python2.7 ../apps/tools/find_intergenic_te_distances.py --Gff $Gff --ID DNA/TcMar >> $OutFile
head -n 1 $OutFile > temp.temp
tail -n +2 $OutFile | awk -F'\t' '{print $1 ".t1\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' >> temp.temp && mv temp.temp $OutFile
x=$(grep 'gene' $Gff | wc -l)
y=$(cat $OutFile | wc -l)
echo "Genes in $ID GFF: $x, distances found for: $y"
done

for Gff in $(ls /home/theaven/projects/niab/theaven/gene_pred/P_*/*/codingquarry/rep_modeling/final/te.gff); do
  ID=$(echo $Gff | cut -d '/' -f9)
OutFile=$(dirname $Gff)/${ID}_flanking_Tads.txt
echo -e "ID\tfive_prime_lgth\tfive_prime_family\tfive_prime_ID\tfive_prime_strand\tthree_prime_lgth\tthree_prime_family\tthree_prime_ID\tthree_prime_strand\tstrand"  > $OutFile
python2.7 ../apps/tools/find_intergenic_te_distances.py --Gff $Gff --ID LINE/Tad1 >> $OutFile
head -n 1 $OutFile > temp.temp
tail -n +2 $OutFile | awk -F'\t' '{print $1 ".t1\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' >> temp.temp && mv temp.temp $OutFile
x=$(grep 'gene' $Gff | wc -l)
y=$(cat $OutFile | wc -l)
echo "Genes in $ID GFF: $x, distances found for: $y"
done

for Gff in $(ls /home/theaven/projects/niab/theaven/gene_pred/P_*/*/codingquarry/rep_modeling/final/te.gff); do
  ID=$(echo $Gff | cut -d '/' -f9)
OutFile=$(dirname $Gff)/${ID}_flanking_Gypsys.txt
echo -e "ID\tfive_prime_lgth\tfive_prime_family\tfive_prime_ID\tfive_prime_strand\tthree_prime_lgth\tthree_prime_family\tthree_prime_ID\tthree_prime_strand\tstrand"  > $OutFile
python2.7 ../apps/tools/find_intergenic_te_distances.py --Gff $Gff --ID LTR/Gypsy >> $OutFile
head -n 1 $OutFile > temp.temp
tail -n +2 $OutFile | awk -F'\t' '{print $1 ".t1\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' >> temp.temp && mv temp.temp $OutFile
x=$(grep 'gene' $Gff | wc -l)
y=$(cat $OutFile | wc -l)
echo "Genes in $ID GFF: $x, distances found for: $y"
done

for Gff in $(ls /home/theaven/projects/niab/theaven/gene_pred/P_*/*/codingquarry/rep_modeling/final/te.gff); do
  ID=$(echo $Gff | cut -d '/' -f9)
OutFile=$(dirname $Gff)/${ID}_flanking_Copias.txt
echo -e "ID\tfive_prime_lgth\tfive_prime_family\tfive_prime_ID\tfive_prime_strand\tthree_prime_lgth\tthree_prime_family\tthree_prime_ID\tthree_prime_strand\tstrand"  > $OutFile
python2.7 ../apps/tools/find_intergenic_te_distances.py --Gff $Gff --ID LTR/Copia >> $OutFile
head -n 1 $OutFile > temp.temp
tail -n +2 $OutFile | awk -F'\t' '{print $1 ".t1\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' >> temp.temp && mv temp.temp $OutFile
x=$(grep 'gene' $Gff | wc -l)
y=$(cat $OutFile | wc -l)
echo "Genes in $ID GFF: $x, distances found for: $y"
done
conda deactivate

for file in $(ls /home/theaven/projects/niab/theaven/gene_pred/P_*/*/codingquarry/rep_modeling/final/*.txt | grep -v 'intergenic'); do
cp $file tmp08092024/.
done

for file in $(ls /home/theaven/projects/niab/theaven/gene_pred/*/THeaven*/codingquarry/rep_modeling/final/THeaven*_abinitio_annotations3.tsv); do
tail -n 1 $file > temp.tsv
head -n -1 $file >> temp.tsv
cp temp.tsv tmp08092024/$(basename $file)
done
```
g15196
```python
from collections import Counter

# Open the file
with open('/home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/THeavenp11_1_flanking_tes.txt', 'r') as file:
    next(file)
    nested_strings_counter = Counter()
    for line in file:
        columns = line.split()
        if "Nested" in columns[2]:
            nested_strings_counter[columns[4]] += 1

# Print the occurrences of each unique string in the fourth column
for string, count in nested_strings_counter.items():
    print(f"{string}: {count} occurrences")

```
#### Plot gene-to-TE distances
Note-copia and gypsy are LTRs and tad1 is an LINE, pogo is a DNA and mariner.
```R
setwd("C:/Users/did23faz/OneDrive - Norwich Bioscience Institutes/Desktop/R")

#Use merged files created in permutation_vs_busco.R as base
#prefix <- "THeavenp11"
#prefix <- "THeavenpOGB2021"
#prefix <- "THeavenpOGB2019"
#species <- "P_leucotricha"
#prefix <- "THeavenSCOTT2020"
prefix <- "THeavenDRCT72020"
#prefix <- "THeavenDRCT72021"
species <- "P_aphanis"

df3 <- merged_data4 #DRCT72020
#df3 <- merged_data5 #DRCT72021
#df3 <- merged_data6 #SCOTT2020
#df3 <- merged_data7 #OGB2019
#df3 <- merged_data8 #p11
#df3 <- merged_data9 #OGB2021

#To allow facet plotting add lines for each TE category measurement

any <- read.table(paste0("download/tmp08092024/", prefix, "_1_flanking_tes.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(any)[1] <- "name"
any <- subset(any, select = -strand)

copia_te <- read.table(paste0("download/tmp08092024/", prefix, "_1_flanking_Copias.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(copia_te)[1] <- "name"
copia_te <- subset(copia_te, select = -strand)

DNA_te <- read.table(paste0("download/tmp08092024/", prefix, "_1_flanking_DNAs.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(DNA_te)[1] <- "name"
DNA_te <- subset(DNA_te, select = -strand)

Gypsy_te <- read.table(paste0("download/tmp08092024/", prefix, "_1_flanking_Gypsys.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(Gypsy_te)[1] <- "name"
Gypsy_te <- subset(Gypsy_te, select = -strand)

LINE_te <- read.table(paste0("download/tmp08092024/", prefix, "_1_flanking_LINEs.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(LINE_te)[1] <- "name"
LINE_te <- subset(LINE_te, select = -strand)

LTR_te <- read.table(paste0("download/tmp08092024/", prefix, "_1_flanking_LTRs.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(LTR_te)[1] <- "name"
LTR_te <- subset(LTR_te, select = -strand)

Mariner_te <- read.table(paste0("download/tmp08092024/", prefix, "_1_flanking_Mariners.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(Mariner_te)[1] <- "name"
Mariner_te <- subset(Mariner_te, select = -strand)

Pogo_te <- read.table(paste0("download/tmp08092024/", prefix, "_1_flanking_Pogos.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(Pogo_te)[1] <- "name"
Pogo_te <- subset(Pogo_te, select = -strand)

Tad_te <- read.table(paste0("download/tmp08092024/", prefix, "_1_flanking_Tads.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(Tad_te)[1] <- "name"
Tad_te <- subset(Tad_te, select = -strand)

df3$sec_type <- 'Non CSEP'
df3$sec_type[df3$outgroup_orthology_match == 0 & df3$is_secreted == 1 & df3$effectorp3_noneffector == "." & df3$residue_number < 201] <- "Short CSEP"
df3$sec_type[df3$outgroup_orthology_match == 0 & df3$is_secreted == 1 & df3$effectorp3_noneffector == "." & df3$residue_number > 200] <- "Long CSEP"
df3$sec_type <- factor(df3$sec_type, levels = c("Non CSEP", "Short CSEP", "Long CSEP"))
count_short_csep <- sum(df3$sec_type == "Short CSEP")
count_long_csep <- sum(df3$sec_type == "Long CSEP")
count_non_csep <- sum(df3$sec_type == "Non CSEP")

any_te <- merge(df3, any, by = "name", all = FALSE)
any_te$te_group <- 'Any TE'
copia_te <- merge(df3, copia_te, by = "name", all = FALSE)
copia_te$te_group <- 'Copia'
DNA_te <- merge(df3, DNA_te, by = "name", all = FALSE)
DNA_te$te_group <- 'DNA'
Gypsy_te <- merge(df3, Gypsy_te, by = "name", all = FALSE)
Gypsy_te$te_group <- 'Gypsy'
LINE_te <- merge(df3, LINE_te, by = "name", all = FALSE)
LINE_te$te_group <- 'LINE'
LTR_te <- merge(df3, LTR_te, by = "name", all = FALSE)
LTR_te$te_group <- 'LTR'
Mariner_te <- merge(df3, Mariner_te, by = "name", all = FALSE)
Mariner_te$te_group <- 'Mariner'
Pogo_te <- merge(df3, Pogo_te, by = "name", all = FALSE)
Pogo_te$te_group <- 'Pogo'
Tad_te <- merge(df3, Tad_te, by = "name", all = FALSE)
Tad_te$te_group <- 'Tad'

gene_te <- merge(df3, any, by = "name", all = FALSE)
gene_te$five_prime_lgth <- gene_te$IG_five_prime
gene_te$five_prime_family <- 'gene'
gene_te$five_prime_ID <- 'gene'
gene_te$three_prime_lgth <- gene_te$IG_three_prime
gene_te$three_prime_family <- 'gene'
gene_te$three_prime_ID <- 'gene'
gene_te$te_group <- 'Gene'
gene_te$three_prime_lgth[gene_te$three_prime_lgth == '99999'] <- NA
gene_te$five_prime_lgth[gene_te$five_prime_lgth == '99999'] <- NA

df <- rbind(any_te,copia_te,DNA_te,Gypsy_te,LINE_te,LTR_te,Mariner_te,Pogo_te,Tad_te,gene_te)
allte <- df

allte$three_prime_lgth[allte$three_prime_lgth == 'Nested TE'] <- 0
allte$five_prime_lgth[allte$five_prime_lgth == 'Nested TE'] <- 0

allte$three_prime_lgth[allte$three_prime_lgth == 'gene'] <- NA
allte$five_prime_lgth[allte$five_prime_lgth == 'gene'] <- NA

allte$three_prime_lgth[allte$three_prime_lgth == 'end'] <- NA
allte$five_prime_lgth[allte$five_prime_lgth == 'end'] <- NA

allte$three_prime_lgth <- as.numeric(allte$three_prime_lgth)
allte$five_prime_lgth <- as.numeric(allte$five_prime_lgth)

te_group_order <- c("Gene", "Any TE", "LINE", "Tad", "LTR", "Copia", "Gypsy", "DNA", "Mariner", "Pogo")
allte$te_group <- factor(allte$te_group, levels = te_group_order)

mean_3_gene_value <- mean(allte$three_prime_lgth[grepl("Gene", allte$te_group) & grepl("Non CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_gene_value <- mean(allte$five_prime_lgth[grepl("Gene", allte$te_group) & grepl("Non CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_AnyTE_value <- mean(allte$three_prime_lgth[grepl("Any TE", allte$te_group) & grepl("Non CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_AnyTE_value <- mean(allte$five_prime_lgth[grepl("Any TE", allte$te_group) & grepl("Non CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_Copia_value <- mean(allte$three_prime_lgth[grepl("Copia", allte$te_group) & grepl("Non CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_Copia_value <- mean(allte$five_prime_lgth[grepl("Copia", allte$te_group) & grepl("Non CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_DNA_value <- mean(allte$three_prime_lgth[grepl("DNA", allte$te_group) & grepl("Non CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_DNA_value <- mean(allte$five_prime_lgth[grepl("DNA", allte$te_group) & grepl("Non CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_Gypsy_value <- mean(allte$three_prime_lgth[grepl("Gypsy", allte$te_group) & grepl("Non CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_Gypsy_value <- mean(allte$five_prime_lgth[grepl("Gypsy", allte$te_group) & grepl("Non CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_LINE_value <- mean(allte$three_prime_lgth[grepl("LINE", allte$te_group) & grepl("Non CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_LINE_value <- mean(allte$five_prime_lgth[grepl("LINE", allte$te_group) & grepl("Non CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_LTR_value <- mean(allte$three_prime_lgth[grepl("LTR", allte$te_group) & grepl("Non CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_LTR_value <- mean(allte$five_prime_lgth[grepl("LTR", allte$te_group) & grepl("Non CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_Mariner_value <- mean(allte$three_prime_lgth[grepl("Mariner", allte$te_group) & grepl("Non CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_Mariner_value <- mean(allte$five_prime_lgth[grepl("Mariner", allte$te_group) & grepl("Non CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_Pogo_value <- mean(allte$three_prime_lgth[grepl("Pogo", allte$te_group) & grepl("Non CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_Pogo_value <- mean(allte$five_prime_lgth[grepl("Pogo", allte$te_group) & grepl("Non CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_Tad_value <- mean(allte$three_prime_lgth[grepl("Tad", allte$te_group) & grepl("Non CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_Tad_value <- mean(allte$five_prime_lgth[grepl("Tad", allte$te_group) & grepl("Non CSEP", allte$sec_type)], na.rm = TRUE)

text <- print(paste(prefix, "of", count_non_csep, "Non CSEP", "Gene 3", mean_3_gene_value, "Gene 5", mean_5_gene_value, "Any TE 3", mean_3_AnyTE_value, "Any TE 5", mean_5_AnyTE_value, "Copia 3", mean_3_Copia_value, "Copia 5", mean_5_Copia_value, "DNA 3", mean_3_DNA_value, "DNA 5", mean_5_DNA_value, "Gypsy 3", mean_3_Gypsy_value, "Gypsy 5", mean_5_Gypsy_value, "LINE 3", mean_3_LINE_value, "LINE 5", mean_5_LINE_value, "LTR 3", mean_3_LTR_value, "LTR 5", mean_5_LTR_value, "Mariner 3", mean_3_Mariner_value, "Mariner 5", mean_5_Mariner_value, "Pogo 3", mean_3_Pogo_value, "Pogo 5", mean_5_Pogo_value, "Tad 3", mean_3_Tad_value, "Tad 5", mean_5_Tad_value))

mean_3_gene_value <- mean(allte$three_prime_lgth[grepl("Gene", allte$te_group) & grepl("Short CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_gene_value <- mean(allte$five_prime_lgth[grepl("Gene", allte$te_group) & grepl("Short CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_AnyTE_value <- mean(allte$three_prime_lgth[grepl("Any TE", allte$te_group) & grepl("Short CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_AnyTE_value <- mean(allte$five_prime_lgth[grepl("Any TE", allte$te_group) & grepl("Short CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_Copia_value <- mean(allte$three_prime_lgth[grepl("Copia", allte$te_group) & grepl("Short CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_Copia_value <- mean(allte$five_prime_lgth[grepl("Copia", allte$te_group) & grepl("Short CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_DNA_value <- mean(allte$three_prime_lgth[grepl("DNA", allte$te_group) & grepl("Short CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_DNA_value <- mean(allte$five_prime_lgth[grepl("DNA", allte$te_group) & grepl("Short CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_Gypsy_value <- mean(allte$three_prime_lgth[grepl("Gypsy", allte$te_group) & grepl("Short CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_Gypsy_value <- mean(allte$five_prime_lgth[grepl("Gypsy", allte$te_group) & grepl("Short CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_LINE_value <- mean(allte$three_prime_lgth[grepl("LINE", allte$te_group) & grepl("Short CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_LINE_value <- mean(allte$five_prime_lgth[grepl("LINE", allte$te_group) & grepl("Short CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_LTR_value <- mean(allte$three_prime_lgth[grepl("LTR", allte$te_group) & grepl("Short CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_LTR_value <- mean(allte$five_prime_lgth[grepl("LTR", allte$te_group) & grepl("Short CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_Mariner_value <- mean(allte$three_prime_lgth[grepl("Mariner", allte$te_group) & grepl("Short CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_Mariner_value <- mean(allte$five_prime_lgth[grepl("Mariner", allte$te_group) & grepl("Short CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_Pogo_value <- mean(allte$three_prime_lgth[grepl("Pogo", allte$te_group) & grepl("Short CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_Pogo_value <- mean(allte$five_prime_lgth[grepl("Pogo", allte$te_group) & grepl("Short CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_Tad_value <- mean(allte$three_prime_lgth[grepl("Tad", allte$te_group) & grepl("Short CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_Tad_value <- mean(allte$five_prime_lgth[grepl("Tad", allte$te_group) & grepl("Short CSEP", allte$sec_type)], na.rm = TRUE)

text2 <- print(paste(prefix, "of", count_short_csep, "Short CSEP", "Gene 3", mean_3_gene_value, "Gene 5", mean_5_gene_value, "Any TE 3", mean_3_AnyTE_value, "Any TE 5", mean_5_AnyTE_value, "Copia 3", mean_3_Copia_value, "Copia 5", mean_5_Copia_value, "DNA 3", mean_3_DNA_value, "DNA 5", mean_5_DNA_value, "Gypsy 3", mean_3_Gypsy_value, "Gypsy 5", mean_5_Gypsy_value, "LINE 3", mean_3_LINE_value, "LINE 5", mean_5_LINE_value, "LTR 3", mean_3_LTR_value, "LTR 5", mean_5_LTR_value, "Mariner 3", mean_3_Mariner_value, "Mariner 5", mean_5_Mariner_value, "Pogo 3", mean_3_Pogo_value, "Pogo 5", mean_5_Pogo_value, "Tad 3", mean_3_Tad_value, "Tad 5", mean_5_Tad_value))

mean_3_gene_value <- mean(allte$three_prime_lgth[grepl("Gene", allte$te_group) & grepl("Long CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_gene_value <- mean(allte$five_prime_lgth[grepl("Gene", allte$te_group) & grepl("Long CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_AnyTE_value <- mean(allte$three_prime_lgth[grepl("Any TE", allte$te_group) & grepl("Long CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_AnyTE_value <- mean(allte$five_prime_lgth[grepl("Any TE", allte$te_group) & grepl("Long CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_Copia_value <- mean(allte$three_prime_lgth[grepl("Copia", allte$te_group) & grepl("Long CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_Copia_value <- mean(allte$five_prime_lgth[grepl("Copia", allte$te_group) & grepl("Long CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_DNA_value <- mean(allte$three_prime_lgth[grepl("DNA", allte$te_group) & grepl("Long CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_DNA_value <- mean(allte$five_prime_lgth[grepl("DNA", allte$te_group) & grepl("Long CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_Gypsy_value <- mean(allte$three_prime_lgth[grepl("Gypsy", allte$te_group) & grepl("Long CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_Gypsy_value <- mean(allte$five_prime_lgth[grepl("Gypsy", allte$te_group) & grepl("Long CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_LINE_value <- mean(allte$three_prime_lgth[grepl("LINE", allte$te_group) & grepl("Long CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_LINE_value <- mean(allte$five_prime_lgth[grepl("LINE", allte$te_group) & grepl("Long CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_LTR_value <- mean(allte$three_prime_lgth[grepl("LTR", allte$te_group) & grepl("Long CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_LTR_value <- mean(allte$five_prime_lgth[grepl("LTR", allte$te_group) & grepl("Long CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_Mariner_value <- mean(allte$three_prime_lgth[grepl("Mariner", allte$te_group) & grepl("Long CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_Mariner_value <- mean(allte$five_prime_lgth[grepl("Mariner", allte$te_group) & grepl("Long CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_Pogo_value <- mean(allte$three_prime_lgth[grepl("Pogo", allte$te_group) & grepl("Long CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_Pogo_value <- mean(allte$five_prime_lgth[grepl("Pogo", allte$te_group) & grepl("Long CSEP", allte$sec_type)], na.rm = TRUE)

mean_3_Tad_value <- mean(allte$three_prime_lgth[grepl("Tad", allte$te_group) & grepl("Long CSEP", allte$sec_type)], na.rm = TRUE)
mean_5_Tad_value <- mean(allte$five_prime_lgth[grepl("Tad", allte$te_group) & grepl("Long CSEP", allte$sec_type)], na.rm = TRUE)

text3 <- print(paste(prefix, "of", count_long_csep, "Long CSEP", "Gene 3", mean_3_gene_value, "Gene 5", mean_5_gene_value, "Any TE 3", mean_3_AnyTE_value, "Any TE 5", mean_5_AnyTE_value, "Copia 3", mean_3_Copia_value, "Copia 5", mean_5_Copia_value, "DNA 3", mean_3_DNA_value, "DNA 5", mean_5_DNA_value, "Gypsy 3", mean_3_Gypsy_value, "Gypsy 5", mean_5_Gypsy_value, "LINE 3", mean_3_LINE_value, "LINE 5", mean_5_LINE_value, "LTR 3", mean_3_LTR_value, "LTR 5", mean_5_LTR_value, "Mariner 3", mean_3_Mariner_value, "Mariner 5", mean_5_Mariner_value, "Pogo 3", mean_3_Pogo_value, "Pogo 5", mean_5_Pogo_value, "Tad 3", mean_3_Tad_value, "Tad 5", mean_5_Tad_value))


cat(text, "\n", file = "TE_distance.txt", append = TRUE)
cat(text2, "\n", file = "TE_distance.txt", append = TRUE)
cat(text3, "\n", file = "TE_distance.txt", append = TRUE)
wb <- createWorkbook()
addWorksheet(wb, "Sheet1")
writeData(wb, "Sheet1", df)
saveWorkbook(wb, paste0(species, "_",  prefix, "_annotations_master2.xlsx"), overwrite = TRUE)

# Define custom colors
te_group_colors <- c("Gene" = "#eeeeee", "Any TE" = "#b6d7a8" , "LINE" = "#9fc5e8", "Tad" = "#9fc5e8", "LTR" = "#e06666", "Copia" = "#e06666", "Gypsy" = "#e06666", "DNA" = "#ffd966", "Mariner" = "#ffd966", "Pogo" = "#ffd966")

library(ggplot2)
library(hrbrthemes)
library(scales)
library(dplyr)

counts <- allte %>%
  group_by(te_group, sec_type) %>%
  summarise(n = sum(!is.na(three_prime_lgth)), 
            median_y = median(as.numeric(three_prime_lgth[as.numeric(three_prime_lgth) != 0]), na.rm = TRUE))

ggplot(allte, aes(x = te_group, y = three_prime_lgth)) +
    geom_point(alpha = 0.1, position = 'jitter') +
    geom_violin(aes(fill = te_group)) +
    geom_text(data = counts, aes(x = te_group, y = median_y, label = n), position = position_dodge(width = 0.9), angle = 90, vjust = 0.3, size = 3) +
    facet_wrap(~sec_type) +
    scale_fill_manual(values = te_group_colors, name = "TE Group") +
    scale_x_discrete(limits = te_group_order) +
    scale_y_log10(labels = scales::label_number(accuracy = 1)) +  
    theme_ipsum() +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),  
        axis.title.x = element_text(hjust = 0.5),  
        axis.title.y = element_text(hjust = 0.5),  
        plot.title = element_text(hjust = 0.5),  
        strip.text = element_text(hjust = 0.5),  
        plot.margin = margin(10, 10, 10, 10)  
    ) +
    ylab("3' intergenic distance (bp)") +
    xlab("Repetitive element type (bp)") 

counts <- allte %>%
  group_by(te_group, sec_type) %>%
  summarise(n = sum(!is.na(five_prime_lgth)), 
            median_y = median(as.numeric(five_prime_lgth[as.numeric(five_prime_lgth) != 0]), na.rm = TRUE))

ggplot(allte, aes(x = te_group, y = five_prime_lgth)) +
    geom_point(alpha = 0.1, position = 'jitter') +
    geom_violin(aes(fill = te_group)) +
    geom_text(data = counts, aes(x = te_group, y = median_y, label = n), position = position_dodge(width = 0.9), angle = 90, vjust = 0.3, size = 3) +
    facet_wrap(~sec_type) +
    scale_fill_manual(values = te_group_colors, name = "TE Group") +
    scale_x_discrete(limits = te_group_order) +
    scale_y_log10(labels = scales::label_number(accuracy = 1)) +  
    theme_ipsum() +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),  
        axis.title.x = element_text(hjust = 0.5),  
        axis.title.y = element_text(hjust = 0.5),  
        plot.title = element_text(hjust = 0.5),  
        strip.text = element_text(hjust = 0.5), 
        plot.margin = margin(10, 10, 10, 10)  
    ) +
    ylab("5' intergenic distance (bp)") +
    xlab("Repetitive element type (bp)") 
```
```R
library(readxl)

df <- read_excel("P_aphanis_THeavenDRCT72020_annotations_master.xlsx")
df <- read_excel("P_aphanis_THeavenDRCT72021_annotations_master.xlsx")
df <- read_excel("P_aphanis_THeavenSCOTT2020_annotations_master.xlsx")
df <- read_excel("P_leucotricha_THeavenpOGB2019_annotations_master.xlsx")
df <- read_excel("P_leucotricha_THeavenp11_annotations_master.xlsx")
df <- read_excel("P_leucotricha_THeavenpOGB2021_annotations_master.xlsx")

df$BUSCO <- "Non"
df$BUSCO[df$BUSCO_blast_match == 1] <- "BUSCO"

subset_data <- df[df$BUSCO == 'BUSCO' & df$te_group == 'Any TE' & df$IG_five_prime != 99999, ]
shapiro.test(subset_data$IG_five_prime)
subset_data <- df[df$BUSCO == 'BUSCO' & df$te_group == 'Any TE' & df$IG_three_prime != 99999, ]
shapiro.test(subset_data$IG_five_prime)
```

### Intergenic distances
#### Calculate integenic distances for the six genomes

```bash
ls ~/projects/niab/theaven/gene_pred/*/*/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3
ls ~/projects/niab/theaven/gene_pred/*/*/predector_singularity3/results/final_genes_appended_renamed.pep/P*-ranked.tsv

../apps/tools/find_intergenic_regions.py --Gff $GeneGff > $OutDir/10300_intergenic_regions.txt

#Find the 3' and 5' intergenic distances + match names to the predector output files
conda activate predector2.7
for Gff in $(ls ~/projects/niab/theaven/gene_pred/*/*/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3); do
  ID=$(echo $Gff | cut -d '/' -f9)
OutFile=$(dirname $Gff)/${ID}_intergenic_regions_nonest+.txt
echo -e "name\tfive_prime\tthree_prime\tstrand\tfive_prime_neighbour_strand\tthree_prime_neighbour_strand" > $OutFile
python2.7 ../apps/tools/find_intergenic_regions_nonest.py --Gff $Gff >> $OutFile
head -n 1 $OutFile > temp.temp
tail -n +2 $OutFile | awk -F'\t' '{print $1 ".t1\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' >> temp.temp && mv temp.temp $OutFile
x=$(grep 'gene' $Gff | wc -l)
y=$(cat $OutFile | wc -l)
echo "Genes in $ID GFF: $x, intergenic regions found for: $y"
done
conda deactivate

cp ~/projects/niab/theaven/gene_pred/*/*/codingquarry/rep_modeling/final/*intergenic_regions_nonest.txt tmp08092024/.

#Investigate missing genes:
for Gff in $(ls ~/projects/niab/theaven/gene_pred/*/*/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3); do
file=$(dirname $Gff)/intergenic_regions.txt
while IFS= read -r value; do
    # Extract the value from the first column
    search_value=$(echo "$value" | grep 'mRNA' | awk '{print $9}' | cut -d ';' -f1 | cut -d '=' -f2)
    # Search for the value in the first column of the GFF file
    if ! grep -q "$search_value" "$file"; then
        # If no match found, print the value
        echo "$value" >> $(dirname $file)/missing.txt
    fi
done < "$Gff"
done

for file in $(ls ~/projects/niab/theaven/gene_pred/*/*/codingquarry/rep_modeling/final/intergenic_regions_nonest.txt); do 
echo $Gff | cut -d '/' -f9
max_column_2=$(awk 'NR>1 { if ($2 > max) max = $2 } END { print max }' "$file")
max_column_3=$(awk 'NR>1 { if ($3 > max) max = $3 } END { print max }' "$file")
echo "Maximum 5' value: $max_column_2"
echo "Maximum 3' value: $max_column_3"
done

for file in $(ls ~/projects/niab/theaven/gene_pred/*/*/codingquarry/rep_modeling/final/*intergenic_regions_nonest+.txt); do
echo $(echo $file | cut -d'/' -f9)
python2.7 ../apps/tools/find_intergenic_regions_orient.py $file
done
```
THeavenDRCT72020_1
('Average different five_prime value:', 1644.81777333999)
('Average same five_prime value:', 1380.3653256051823)
('Average different three_prime value:', 1701.2392550143265)
('Average same three_prime value:', 1380.3653256051823)
THeavenDRCT72021_1
('Average different five_prime value:', 1486.19154676259)
('Average same five_prime value:', 1364.3916256157636)
('Average different three_prime value:', 1468.642551369863)
('Average same three_prime value:', 1364.3916256157636)
THeavenSCOTT2020_1
('Average different five_prime value:', 1290.6614173228347)
('Average same five_prime value:', 1193.3056271477662)
('Average different three_prime value:', 1529.4132543103449)
('Average same three_prime value:', 1193.3056271477662)
THeavenp11_1
('Average different five_prime value:', 1690.9040952047603)
('Average same five_prime value:', 1358.8400804289545)
('Average different three_prime value:', 1659.2817294281729)
('Average same three_prime value:', 1358.8400804289545)
THeavenpOGB2019_1
('Average different five_prime value:', 1599.4325410110478)
('Average same five_prime value:', 1322.509015971149)
('Average different three_prime value:', 1585.5731667214732)
('Average same three_prime value:', 1322.509015971149)
THeavenpOGB2021_1
('Average different five_prime value:', 1698.1465360899529)
('Average same five_prime value:', 1405.0219162686149)
('Average different three_prime value:', 1711.6103019538189)
('Average same three_prime value:', 1405.0219162686149)
```R
library(readxl)

DRCT72020 <- read_excel("P_aphanis_THeavenDRCT72020_annotations_master.xlsx")
DRCT72021 <- read_excel("P_aphanis_THeavenDRCT72021_annotations_master.xlsx")
SCOTT2020 <- read_excel("P_aphanis_THeavenSCOTT2020_annotations_master.xlsx")
OGB2019 <- read_excel("P_leucotricha_THeavenpOGB2019_annotations_master.xlsx")
P112020 <- read_excel("P_leucotricha_THeavenp11_annotations_master.xlsx")
OGB2021 <- read_excel("P_leucotricha_THeavenpOGB2021_annotations_master.xlsx")

DRCT72020_nonest <- read.table("download/download11092024/THeavenDRCT72020_1_intergenic_regions_nonest+.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
DRCT72021_nonest <- read.table("download/download11092024/THeavenDRCT72021_1_intergenic_regions_nonest+.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
SCOTT2020_nonest <- read.table("download/download11092024/THeavenSCOTT2020_1_intergenic_regions_nonest+.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
OGB2019_nonest <- read.table("download/download11092024/THeavenpOGB2019_1_intergenic_regions_nonest+.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
P112020_nonest <- read.table("download/download11092024/THeavenp11_1_intergenic_regions_nonest+.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
OGB2021_nonest <- read.table("download/download11092024/THeavenpOGB2021_1_intergenic_regions_nonest+.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

DRCT72020_1 <- merge(DRCT72020, DRCT72020_nonest,  by = "name", all = TRUE)
DRCT72021_1 <- merge(DRCT72021, DRCT72021_nonest,  by = "name", all = TRUE)
SCOTT2020_1 <- merge(SCOTT2020, SCOTT2020_nonest,  by = "name", all = TRUE)
OGB2019_1 <- merge(OGB2019, OGB2019_nonest,  by = "name", all = TRUE)
P112020_1 <- merge(P112020, P112020_nonest,  by = "name", all = TRUE)
OGB2021_1 <- merge(OGB2021, OGB2021_nonest,  by = "name", all = TRUE)

genomes_list2 <- list(DRCT72020_1, DRCT72021_1, SCOTT2020_1, OGB2019_1, P112020_1, OGB2021_1)
genome_names2 <- c("DRCT72020", "DRCT72021", "SCOTT2020", "OGB2019", "P112020", "OGB2021")

for (k in seq_along(genomes_list2)) {
df1 <- genomes_list2[[k]]
df_name <- genome_names2[k]

df <- df1 %>% filter(te_group == "Any TE")

df$IG_three_prime[df$IG_three_prime == 99999] <- NA
df$IG_five_prime[df$IG_five_prime == 99999] <- NA

head_to_head_csep <- mean(df$IG_five_prime[
                            df$sec_type != "Non CSEP" & 
                            df$strand.x != df$five_prime_neighbour_strand & !is.na(df$five_prime_neighbour_strand) &
                            !is.na(df$IG_five_prime)], na.rm = TRUE)

parallel_5_csep <- mean(df$IG_five_prime[
                            df$sec_type != "Non CSEP" & 
                            df$strand.x == df$five_prime_neighbour_strand & !is.na(df$five_prime_neighbour_strand) &
                            !is.na(df$IG_five_prime)], na.rm = TRUE)

tail_to_tail_csep <- mean(df$IG_three_prime[
                            df$sec_type != "Non CSEP" & 
                            df$strand.x != df$three_prime_neighbour_strand & !is.na(df$three_prime_neighbour_strand) &
                            !is.na(df$IG_three_prime)], na.rm = TRUE)

parallel_3_csep <- mean(df$IG_three_prime[
                            df$sec_type != "Non CSEP" & 
                            df$strand.x == df$three_prime_neighbour_strand & !is.na(df$three_prime_neighbour_strand) &
                            !is.na(df$IG_three_prime)], na.rm = TRUE)

head_to_head_noncsep <- mean(df$IG_five_prime[
                            df$sec_type == "Non CSEP" & 
                            df$strand.x != df$five_prime_neighbour_strand & !is.na(df$five_prime_neighbour_strand) &
                            !is.na(df$IG_five_prime)], na.rm = TRUE)

parallel_5_noncsep <- mean(df$IG_five_prime[
                            df$sec_type == "Non CSEP" & 
                            df$strand.x == df$five_prime_neighbour_strand & !is.na(df$five_prime_neighbour_strand) &
                            !is.na(df$IG_five_prime)], na.rm = TRUE)

tail_to_tail_noncsep <- mean(df$IG_three_prime[
                            df$sec_type == "Non CSEP" & 
                            df$strand.x != df$three_prime_neighbour_strand & !is.na(df$three_prime_neighbour_strand) &
                            !is.na(df$IG_three_prime)], na.rm = TRUE)

parallel_3_noncsep <- mean(df$IG_three_prime[
                            df$sec_type == "Non CSEP" & 
                            df$strand.x == df$three_prime_neighbour_strand & !is.na(df$three_prime_neighbour_strand) &
                            !is.na(df$IG_three_prime)], na.rm = TRUE)

text1 <- print(paste("For ", df_name, " CSEP head-to-head distance:", head_to_head_csep))
text2 <- print(paste("For ", df_name, " CSEP parallel-5' distance:", parallel_5_csep))
text3 <- print(paste("For ", df_name, " CSEP tail-to-tail distance:", tail_to_tail_csep))
text4 <- print(paste("For ", df_name, " CSEP parallel-3' distance:", parallel_3_csep))
text5 <- print(paste("For ", df_name, " non-CSEP head-to-head distance:", head_to_head_noncsep))
text6 <- print(paste("For ", df_name, " non-CSEP parallel-5' distance:", parallel_5_noncsep))
text7 <- print(paste("For ", df_name, " non-CSEP tail-to-tail distance:", tail_to_tail_noncsep))
text8 <- print(paste("For ", df_name, " non-CSEP parallel-3' distance:", parallel_3_noncsep))
cat(text1, "\n", file = "IGorientationreport.txt", append = TRUE)
cat(text2, "\n", file = "IGorientationreport.txt", append = TRUE)
cat(text3, "\n", file = "IGorientationreport.txt", append = TRUE)
cat(text4, "\n", file = "IGorientationreport.txt", append = TRUE)
cat(text5, "\n", file = "IGorientationreport.txt", append = TRUE)
cat(text6, "\n", file = "IGorientationreport.txt", append = TRUE)
cat(text7, "\n", file = "IGorientationreport.txt", append = TRUE)
cat(text8, "\n", file = "IGorientationreport.txt", append = TRUE)

}



mean(df$IG_five_prime[df$sec_type != "Non CSEP" & df$strand != df$five_prime_neighbour_strand], na.rm = TRUE)


[
                            df$sec_type != "Non CSEP" & 
                            df$strand != df$five_prime_neighbour_strand & !is.na(df$five_prime_neighbour_strand) &
                            !is.na(df$IG_five_prime)], na.rm = TRUE)
```

```python
import sys
filename = sys.argv[1]  
with open(filename, 'r') as file:
    lines = file.readlines()

five_prime_values = []

for line in lines[1:]:
    columns = line.strip().split()
    name = columns[0]
    five_prime = int(columns[1])
    strand = columns[3]
    five_prime_neighbour_strand = columns[4]
    if five_prime != 99999 and strand != five_prime_neighbour_strand:
        five_prime_values.append(five_prime)

if five_prime_values:
    average_five_prime = sum(five_prime_values) / float(len(five_prime_values))
    print("Average different five_prime value:", average_five_prime)
else:
    print("No valid rows found for averaging.")

with open(filename, 'r') as file:
    lines = file.readlines()

five_prime_values = []

for line in lines[1:]:
    columns = line.strip().split()
    name = columns[0]
    five_prime = int(columns[1])
    strand = columns[3]
    five_prime_neighbour_strand = columns[4]
    if five_prime != 99999 and strand == five_prime_neighbour_strand:
        five_prime_values.append(five_prime)

# Calculate the average of the filtered values
if five_prime_values:
    average_five_prime = sum(five_prime_values) / float(len(five_prime_values))
    print("Average same five_prime value:", average_five_prime)
else:
    print("No valid rows found for averaging.")

with open(filename, 'r') as file:
    lines = file.readlines()

three_prime_values = []

for line in lines[1:]:
    columns = line.strip().split()
    name = columns[0]
    three_prime = int(columns[2])
    strand = columns[3]
    three_prime_neighbour_strand = columns[5]
    if three_prime != 99999 and strand != three_prime_neighbour_strand:
        three_prime_values.append(three_prime)

if three_prime_values:
    average_three_prime = sum(three_prime_values) / float(len(three_prime_values))
    print("Average different three_prime value:", average_three_prime)
else:
    print("No valid rows found for averaging.")

with open(filename, 'r') as file:
    lines = file.readlines()

three_prime_values = []

for line in lines[1:]:
    columns = line.strip().split()
    name = columns[0]
    three_prime = int(columns[2])
    strand = columns[3]
    three_prime_neighbour_strand = columns[5]
    if three_prime != 99999 and strand == three_prime_neighbour_strand:
        three_prime_values.append(three_prime)

if three_prime_values:
    average_three_prime = sum(three_prime_values) / float(len(three_prime_values))
    print("Average same three_prime value:", average_three_prime)
else:
    print("No valid rows found for averaging.")
```
#### Identify gene groups

Identify BUSCO genes:
```bash
screen -S blast
srun -p long --mem 10G --pty bash
conda activate blast+
mkdir busco_db
find gene_pred/Podosphaera/leucotricha/THeavenpOGB2019_1.c/BUSCO/leotiomycetes_odb10/THeavenpOGB2019_1.c/run_leotiomycetes_odb10/busco_sequences/single_copy_busco_sequences/ \
     gene_pred/Podosphaera/leucotricha/THeavenpOGB2021_1.c/BUSCO/leotiomycetes_odb10/THeavenpOGB2021_1.c/run_leotiomycetes_odb10/busco_sequences/single_copy_busco_sequences/ \
     gene_pred/Podosphaera/leucotricha/THeavenp11_1.c/BUSCO/leotiomycetes_odb10/THeavenp11_1.c/run_leotiomycetes_odb10/busco_sequences/single_copy_busco_sequences/ \
     gene_pred/Podosphaera/aphanis/THeavenDRCT72020_1.c/BUSCO/leotiomycetes_odb10/THeavenDRCT72020_1.c/run_leotiomycetes_odb10/busco_sequences/single_copy_busco_sequences/ \
     gene_pred/Podosphaera/aphanis/THeavenDRCT72021_1.c/BUSCO/leotiomycetes_odb10/THeavenDRCT72021_1.c/run_leotiomycetes_odb10/busco_sequences/single_copy_busco_sequences/ \
     gene_pred/Podosphaera/aphanis/THeavenSCOTT2020_1.c/BUSCO/leotiomycetes_odb10/THeavenSCOTT2020_1.c/run_leotiomycetes_odb10/busco_sequences/single_copy_busco_sequences/ \
     -name '*.faa' -exec cat {} >> busco_db/my_busco_db.faa \;

cd busco_db
makeblastdb -in my_busco_db.faa -dbtype prot -title buscodb -out buscodb

for file in $(ls ~/projects/niab/theaven/gene_pred/*/*/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta); do
  ID=$(echo $file | cut -d '/' -f9)
  blastp -query $file -db buscodb -out ${ID}_busco_blast.tsv -max_target_seqs 1 -evalue 1e-15 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq'
  echo -e "name\tmatch\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqseq\tsseq" > temp.tsv
  cat ${ID}_busco_blast.tsv >> temp.tsv && mv temp.tsv ${ID}_busco_blast2.tsv
done
```
As combined files have been lost (see above), CSEPs will need to be recalculated.

Identify genes with non-homology (requirement for CSEPs):
```bash
#Make blast database as in Praz et al., 2017
cat podosporaanserina_GCF_000226545.1_ASM22654v1_protein.faa neurosporacrassa_GCF_000182925.2_NC12_protein.faa > mydb.faa
mkdir /home/theaven/scratch/uncompressed/csep_db/db
makeblastdb -in mydb.faa -dbtype prot -title CSEPdb -out /home/theaven/scratch/uncompressed/csep_db/db/CSEPdb

#blast for matches
for file in $(ls ~/projects/niab/theaven/gene_pred/*/*/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta); do
  ID=$(echo $file | cut -d '/' -f9)
blastp -query $file -db /home/theaven/scratch/uncompressed/csep_db/db/CSEPdb -out ${ID}_csep_blast.tsv -max_target_seqs 1 -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq'
  echo -e "name\tmatch_1\tpident_1\tlength_1\tmismatch_1\tgapopen_1\tqstart_1\tqend_1\tsstart_1\tsend_1\tevalue_1\tbitscore_1\tqseq_1\tsseq_1" > temp.tsv
  cat ${ID}_csep_blast.tsv >> temp.tsv && mv temp.tsv ${ID}_csep_blast.tsv
done

conda deactivate
```
#### Permutation test
See permutation_vs_BUSCO.R, genes not in 'clean' contigs are excluded therein.

#### Plot intergenic distance distribution

```R
library(ggplot2)
library(patchwork)
library(dplyr)

#See permutation_vs_busco.R for generation of merged_data4 - merged_data9

df_all <- merged_data4 #DRCT72020
#df_all <- merged_data5 #DRCT72021
#df_all <- merged_data6 #SCOTT2020
#df_all <- merged_data7 #OGB2019
#df_all <- merged_data8 #p112020
#df_all <- merged_data9 #OGB2021

df_all$BUSCO <- "Non"
df_all$BUSCO[df_all$BUSCO_blast_match == 1] <- "BUSCO"

df_all$CSEP <- "Non SP"
df_all$CSEP[df_all$outgroup_orthology_match == 0 & df_all$is_secreted == 1 & df_all$effectorp3_noneffector == "."] <- "CSEP"

df_all$cazy <- "Non CAZY"
df_all$cazy[df_all$dbcan_matches != "."] <- "CAZY"

df_all$ralph <- "Non RALPH"
df_all <- df_all %>%
  mutate(ralph = ifelse(grepl("BghBEC1011|BgtAVRa10|BgAVRA13|BgtAvrPm2|BgtSvrPm3a1f1", effector_matches), "RALPH", ralph))

df_all$eka <- "Non EKA"
df_all <- df_all %>%
  mutate(eka = ifelse(grepl("BgtAVRk1|BgtAVRa10", effector_matches), "EKA", eka))
df_all$eka[df_all$is_secreted == 1] <- "Non EKA"

CSEP_group <- df_all[df_all$CSEP == 'CSEP' & df_all$IG_five_prime != 99999, 'IG_five_prime']
non_CSEP_group <- df_all[df_all$CSEP == 'Non SP' & df_all$BUSCO == 'BUSCO' & df_all$IG_five_prime != 99999, 'IG_five_prime']
t_test_result <- t.test(CSEP_group, non_CSEP_group)
print(t_test_result)

CSEP_group <- df_all[df_all$CSEP == 'CSEP' & df_all$IG_three_prime != 99999, 'IG_three_prime']
non_CSEP_group <- df_all[df_all$CSEP == 'Non SP' & df_all$BUSCO == 'BUSCO' & df_all$IG_three_prime != 99999, 'IG_three_prime']
t_test_result <- t.test(CSEP_group, non_CSEP_group)
print(t_test_result)

cazy_group <- df_all[df_all$cazy == 'CAZY' & df_all$IG_five_prime != 99999, 'IG_five_prime']
non_cazy_group <- df_all[df_all$cazy == 'Non CAZY' & df_all$BUSCO == 'BUSCO' & df_all$IG_five_prime != 99999, 'IG_five_prime']
t_test_result <- t.test(cazy_group, non_cazy_group)
print(t_test_result)

cazy_group <- df_all[df_all$cazy == 'CAZY' & df_all$IG_three_prime != 99999, 'IG_three_prime']
non_cazy_group <- df_all[df_all$cazy == 'Non CAZY' & df_all$BUSCO == 'BUSCO' & df_all$IG_three_prime != 99999, 'IG_three_prime']
t_test_result <- t.test(cazy_group, non_cazy_group)
print(t_test_result)

ralph_group <- df_all[df_all$ralph == 'RALPH' & df_all$IG_five_prime != 99999, 'IG_five_prime']
non_ralph_group <- df_all[df_all$ralph == 'Non RALPH' & df_all$BUSCO == 'BUSCO' & df_all$IG_five_prime != 99999, 'IG_five_prime']
t_test_result <- t.test(ralph_group, non_ralph_group)
print(t_test_result)

ralph_group <- df_all[df_all$ralph == 'RALPH' & df_all$IG_three_prime != 99999, 'IG_three_prime']
non_ralph_group <- df_all[df_all$ralph == 'Non RALPH' & df_all$BUSCO == 'BUSCO' & df_all$IG_three_prime != 99999, 'IG_three_prime']
t_test_result <- t.test(ralph_group, non_ralph_group)
print(t_test_result)

eka_group <- df_all[df_all$eka == 'EKA' & df_all$IG_five_prime != 99999, 'IG_five_prime']
non_eka_group <- df_all[df_all$eka == 'Non EKA' & df_all$BUSCO == 'BUSCO' & df_all$IG_five_prime != 99999, 'IG_five_prime']
t_test_result <- t.test(eka_group, non_eka_group)
print(t_test_result)

eka_group <- df_all[df_all$eka == 'EKA' & df_all$IG_three_prime != 99999, 'IG_three_prime']
non_eka_group <- df_all[df_all$eka == 'Non EKA' & df_all$BUSCO == 'BUSCO' & df_all$IG_three_prime != 99999, 'IG_three_prime']
t_test_result <- t.test(eka_group, non_eka_group)
print(t_test_result)

subset_data <- df_all[df_all$BUSCO == 'BUSCO' & df_all$IG_five_prime != 99999, ]
shapiro.test(subset_data$IG_five_prime)
hist(subset_data$IG_five_prime, 
     breaks = 1000, 
     main = "Histogram of IG_five_prime (BUSCO == 'BUSCO')", 
     xlab = "IG_five_prime", 
     probability = TRUE, 
     col = "lightblue")

subset_data <- df_all[df_all$BUSCO == 'BUSCO' & df_all$IG_three_prime != 99999, ]
shapiro.test(subset_data$IG_three_prime)
hist(subset_data$IG_three_prime, 
     breaks = 1000, 
     main = "Histogram of IG_five_prime (BUSCO == 'BUSCO')", 
     xlab = "IG_five_prime", 
     probability = TRUE, 
     col = "lightblue")

# Plot 1
df_all$tempvar <- "Total"
plot101 <- ggplot(df_all, aes(x= IG_five_prime, y=IG_three_prime)) +
  geom_hex(binwidth=0.05) +
  scale_fill_distiller(palette = "Spectral", name="Gene\ncount", trans = "log10", breaks = c(0, 10, 100, 1000, 4000)) +
  scale_x_log10(labels = scales::number_format(), breaks = NULL) +
  scale_y_log10(labels = scales::number_format(), breaks = c(1, 10, 100, 1000, 10000, 100000)) +
  ylab("") +
  xlab("") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "grey93"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot101 <- plot101 + facet_grid(. ~ tempvar) +
  theme(strip.background = element_rect(fill="grey35"),
        strip.text = element_text(size=15, colour="white"))

# Plot 2
df_all$tempvar <- "BUSCO"
plot102 <- ggplot(df_all, aes(x = IG_five_prime, y = IG_three_prime)) +
  geom_hex(binwidth = 0.05, fill = "grey") +
  scale_x_log10(labels = scales::number_format(), breaks = NULL) +
  scale_y_log10(labels = scales::number_format(), breaks = NULL) +
  ylab("") +
  xlab("") +
  geom_point(data=subset(df_all, BUSCO == "BUSCO"), color = 'black', fill = 'firebrick2',shape = 21,alpha = 0.7, size = 1) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "grey93"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot102 <- plot102 + facet_grid(. ~ tempvar) +
  theme(strip.background = element_rect(fill="grey35"),
        strip.text = element_text(size=15, colour="white"))

# Plot 3
df_all$tempvar <- "CAZY"
plot103 <- ggplot(df_all, aes(x = IG_five_prime, y = IG_three_prime)) +
  geom_hex(binwidth = 0.05, fill = "grey") +
  scale_x_log10(labels = scales::number_format(), breaks = NULL) +
  scale_y_log10(labels = scales::number_format(), breaks = c(1, 10, 100, 1000, 10000, 100000)) +
  ylab("5' prime intergenic length (bp)") +
  xlab("") +
  geom_point(data=subset(df_all, cazy == "CAZY"), color = 'black', fill = 'firebrick2',shape = 21,alpha = 0.7, size = 1.5) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "grey93"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

plot103 <- plot103 + facet_grid(. ~ tempvar) +
  theme(strip.background = element_rect(fill="grey35"),
        strip.text = element_text(size=15, colour="white"))

# Plot 4
df_all$tempvar <- "CSEP"
plot104 <- ggplot(df_all, aes(x = IG_five_prime, y = IG_three_prime)) +
  geom_hex(binwidth = 0.05, fill = "grey") +
  scale_x_log10(labels = scales::number_format(), breaks = NULL) +
  scale_y_log10(labels = scales::number_format(), breaks = NULL) +
  ylab("") +
  xlab("") +
  geom_point(data=subset(df_all, CSEP == "CSEP"), color = 'black', fill = 'firebrick2',shape = 21,alpha = 0.7, size = 1.5) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "grey93"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

plot104 <- plot104 + facet_grid(. ~ tempvar) +
  theme(strip.background = element_rect(fill="grey35"),
        strip.text = element_text(size=15, colour="white"))

# Plot 5
df_all$tempvar <- "EKA"
plot105 <- ggplot(df_all, aes(x = IG_five_prime, y = IG_three_prime)) +
  geom_hex(binwidth = 0.05, fill = "grey") +
  scale_x_log10(labels = scales::number_format(), breaks = c(1, 10, 100, 1000, 10000, 100000)) +
  scale_y_log10(labels = scales::number_format(), breaks = c(1, 10, 100, 1000, 10000, 100000)) +
  ylab("") +
  xlab("") +
  geom_point(data=subset(df_all, eka == "EKA"), color = 'black', fill = 'firebrick2',shape = 21,alpha = 0.7, size = 1.5) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "grey93"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

plot105 <- plot105 + facet_grid(. ~ tempvar) +
  theme(strip.background = element_rect(fill="grey35"),
        strip.text = element_text(size=15, colour="white"))

# Plot 6
df_all$tempvar <- "RALPH"
plot106 <- ggplot(df_all, aes(x = IG_five_prime, y = IG_three_prime)) +
  geom_hex(binwidth = 0.05, fill = "grey") +
  scale_x_log10(labels = scales::number_format(), breaks = c(1, 10, 100, 1000, 10000, 100000)) +
  scale_y_log10(labels = scales::number_format(), breaks = NULL) +
  ylab("") +
  xlab("3' prime intergenic length (bp)") +
  geom_point(data=subset(df_all, ralph == "RALPH"), color = 'black', fill = 'firebrick2',shape = 21,alpha = 0.7, size = 1.5) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "grey93"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(hjust = -0.6, vjust = -1))  

plot106 <- plot106 + facet_grid(. ~ tempvar) +
  theme(strip.background = element_rect(fill="grey35"),
        strip.text = element_text(size=15, colour="white"))


# Combine plots with axis titles
combined_plot <- plot101 + plot102 + plot103 + plot104 + plot105 + plot106 +
  plot_layout(ncol = 2, nrow = 3, guides = 'collect') 

# Display combined plot
combined_plot

#Export as pdf to retain hexagonal bin shape
```
#### Investigate short IG genes:
```R
df <- merged_data4 #DRCT72020
search <- read.table("DRCT72020_good.txt", header = FALSE)
df_all <- df %>%
  filter(name %in% search$V1)
filtered_rows <- df_all %>%
  filter(five_prime_2 >= 0 & five_prime_2 <= 200)
write.table(filtered_rows$name, "DRCT72020_low5.txt", row.names = FALSE, col.names = FALSE)
filtered_rows <- df_all %>%
  filter(five_prime_2 >= 0 & three_prime_2 <= 200)
write.table(filtered_rows$name, "DRCT72020_low3.txt", row.names = FALSE, col.names = FALSE)

df <- merged_data5 #DRCT72021
search <- read.table("DRCT72021_good.txt", header = FALSE)
df_all <- df %>%
  filter(name %in% search$V1)
filtered_rows <- df_all %>%
  filter(five_prime_2 >= 0 & five_prime_2 <= 200)
write.table(filtered_rows$name, "DRCT72021_low5.txt", row.names = FALSE, col.names = FALSE)
filtered_rows <- df_all %>%
  filter(five_prime_2 >= 0 & three_prime_2 <= 200)
write.table(filtered_rows$name, "DRCT72021_low3.txt", row.names = FALSE, col.names = FALSE)

df <- merged_data6 #SCOTT2020
search <- read.table("SCOTT2020_good.txt", header = FALSE)
df_all <- df %>%
  filter(name %in% search$V1)
filtered_rows <- df_all %>%
  filter(five_prime_2 >= 0 & five_prime_2 <= 200)
write.table(filtered_rows$name, "SCOTT2020_low5.txt", row.names = FALSE, col.names = FALSE)
filtered_rows <- df_all %>%
  filter(five_prime_2 >= 0 & three_prime_2 <= 200)
write.table(filtered_rows$name, "SCOTT2020_low3.txt", row.names = FALSE, col.names = FALSE)

df <- merged_data7 #OGB2019
search <- read.table("OGB2019_good.txt", header = FALSE)
df_all <- df %>%
  filter(name %in% search$V1)
filtered_rows <- df_all %>%
  filter(five_prime_2 >= 0 & five_prime_2 <= 200)
write.table(filtered_rows$name, "OGB2019_low5.txt", row.names = FALSE, col.names = FALSE)
filtered_rows <- df_all %>%
  filter(five_prime_2 >= 0 & three_prime_2 <= 200)
write.table(filtered_rows$name, "OGB2019_low3.txt", row.names = FALSE, col.names = FALSE)

df <- merged_data8 #p112020
search <- read.table("p11_1_good.txt", header = FALSE)
df_all <- df %>%
  filter(name %in% search$V1)
filtered_rows <- df_all %>%
  filter(five_prime_2 >= 0 & five_prime_2 <= 200)
write.table(filtered_rows$name, "p112020_low5.txt", row.names = FALSE, col.names = FALSE)
filtered_rows <- df_all %>%
  filter(five_prime_2 >= 0 & three_prime_2 <= 200)
write.table(filtered_rows$name, "p112020_low3.txt", row.names = FALSE, col.names = FALSE)

df <- merged_data9 #OGB2021
search <- read.table("OGB2021_good.txt", header = FALSE)
df_all <- df %>%
  filter(name %in% search$V1)
filtered_rows <- df_all %>%
  filter(five_prime_2 >= 0 & five_prime_2 <= 200)
write.table(filtered_rows$name, "OGB2021_low5.txt", row.names = FALSE, col.names = FALSE)
filtered_rows <- df_all %>%
  filter(three_prime_2 >= 0 & three_prime_2 <= 200)
write.table(filtered_rows$name, "OGB2021_low3.txt", row.names = FALSE, col.names = FALSE)
```
```bash
sed -i 's/"//g' *low3.txt
sed -i 's/"//g' *low5.txt
sed -i 's/[^[:print:]\t]//g' *low3.txt
sed -i 's/[^[:print:]\t]//g' *low5.txt

grep -f DRCT72020_low5.txt gene_pred/P_aphanis/THeavenDRCT72020_1/interproscan/NRI/P.aphanis_strawberry2020_interproscan.tsv | grep 'GO:' | awk -F'\t' '{split($NF, arr, "|"); for (i in arr) {print $1, arr[i]}}' | sort -u | awk '{print $2}' | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | awk '{print $2, $1}' | sort -nr | awk '{print $2, $1}' > DRCT72020_sub200-5.tsv
grep -f DRCT72021_low5.txt gene_pred/P_aphanis/THeavenDRCT72021_1/interproscan/NRI/P.aphanis_strawberry2021_interproscan.tsv | grep 'GO:' | awk -F'\t' '{split($NF, arr, "|"); for (i in arr) {print $1, arr[i]}}' | sort -u | awk '{print $2}' | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | awk '{print $2, $1}' | sort -nr | awk '{print $2, $1}' > DRCT72021_sub200-5.tsv
grep -f SCOTT2020_low5.txt gene_pred/P_aphanis/THeavenSCOTT2020_1/interproscan/NRI/P.aphanis_raspberry2020_interproscan.tsv | grep 'GO:' | awk -F'\t' '{split($NF, arr, "|"); for (i in arr) {print $1, arr[i]}}' | sort -u | awk '{print $2}' | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | awk '{print $2, $1}' | sort -nr | awk '{print $2, $1}' > SCOTT2020_sub200-5.tsv
grep -f OGB2019_low5.txt gene_pred/P_leucotricha/THeavenpOGB2019_1/interproscan/NRI/P.leucotricha_apple2019_interproscan.tsv | grep 'GO:' | awk -F'\t' '{split($NF, arr, "|"); for (i in arr) {print $1, arr[i]}}' | sort -u | awk '{print $2}' | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | awk '{print $2, $1}' | sort -nr | awk '{print $2, $1}' > OGB2019_sub200-5.tsv
grep -f p112020_low5.txt gene_pred/P_leucotricha/THeavenp11_1/interproscan/NRI/P.leucotricha_apple2020_interproscan.tsv | grep 'GO:' | awk -F'\t' '{split($NF, arr, "|"); for (i in arr) {print $1, arr[i]}}' | sort -u | awk '{print $2}' | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | awk '{print $2, $1}' | sort -nr | awk '{print $2, $1}' > p112020_sub200-5.tsv
grep -f OGB2021_low5.txt gene_pred/P_leucotricha/THeavenpOGB2021_1/interproscan/NRI/P.leucotricha_apple2021_interproscan.tsv | grep 'GO:' | awk -F'\t' '{split($NF, arr, "|"); for (i in arr) {print $1, arr[i]}}' | sort -u | awk '{print $2}' | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | awk '{print $2, $1}' | sort -nr | awk '{print $2, $1}' > OGB2021_sub200-5.tsv

grep -f DRCT72020_low3.txt gene_pred/P_aphanis/THeavenDRCT72020_1/interproscan/NRI/P.aphanis_strawberry2020_interproscan.tsv | grep 'GO:' | awk -F'\t' '{split($NF, arr, "|"); for (i in arr) {print $1, arr[i]}}' | sort -u | awk '{print $2}' | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | awk '{print $2, $1}' | sort -nr | awk '{print $2, $1}' > DRCT72020_sub200-3.tsv
grep -f DRCT72021_low3.txt gene_pred/P_aphanis/THeavenDRCT72021_1/interproscan/NRI/P.aphanis_strawberry2021_interproscan.tsv | grep 'GO:' | awk -F'\t' '{split($NF, arr, "|"); for (i in arr) {print $1, arr[i]}}' | sort -u | awk '{print $2}' | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | awk '{print $2, $1}' | sort -nr | awk '{print $2, $1}' > DRCT72021_sub200-3.tsv
grep -f SCOTT2020_low3.txt gene_pred/P_aphanis/THeavenSCOTT2020_1/interproscan/NRI/P.aphanis_raspberry2020_interproscan.tsv | grep 'GO:' | awk -F'\t' '{split($NF, arr, "|"); for (i in arr) {print $1, arr[i]}}' | sort -u | awk '{print $2}' | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | awk '{print $2, $1}' | sort -nr | awk '{print $2, $1}' > SCOTT2020_sub200-3.tsv
grep -f OGB2019_low3.txt gene_pred/P_leucotricha/THeavenpOGB2019_1/interproscan/NRI/P.leucotricha_apple2019_interproscan.tsv | grep 'GO:' | awk -F'\t' '{split($NF, arr, "|"); for (i in arr) {print $1, arr[i]}}' | sort -u | awk '{print $2}' | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | awk '{print $2, $1}' | sort -nr | awk '{print $2, $1}' > OGB2019_sub200-3.tsv
grep -f p112020_low3.txt gene_pred/P_leucotricha/THeavenp11_1/interproscan/NRI/P.leucotricha_apple2020_interproscan.tsv | grep 'GO:' | awk -F'\t' '{split($NF, arr, "|"); for (i in arr) {print $1, arr[i]}}' | sort -u | awk '{print $2}' | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | awk '{print $2, $1}' | sort -nr | awk '{print $2, $1}' > p112020_sub200-3.tsv
grep -f OGB2021_low3.txt gene_pred/P_leucotricha/THeavenpOGB2021_1/interproscan/NRI/P.leucotricha_apple2021_interproscan.tsv | grep 'GO:' | awk -F'\t' '{split($NF, arr, "|"); for (i in arr) {print $1, arr[i]}}' | sort -u | awk '{print $2}' | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | awk '{print $2, $1}' | sort -nr | awk '{print $2, $1}' > OGB2021_sub200-3.tsv

files=(DRCT72020_sub200-5.tsv DRCT72020_sub200-3.tsv DRCT72021_sub200-5.tsv DRCT72021_sub200-3.tsv SCOTT2020_sub200-5.tsv SCOTT2020_sub200-3.tsv OGB2019_sub200-5.tsv OGB2019_sub200-3.tsv p112020_sub200-5.tsv p112020_sub200-3.tsv OGB2021_sub200-5.tsv OGB2021_sub200-3.tsv)
for file in "${files[@]:1}"; do
    paste -d'\t' combined_lowIG.tsv <(cut -f1 "$file") > temp_merged_file.txt
    mv temp_merged_file.txt combined_lowIG.tsv
done

```















## In house genome assembly Orthology analysis

```bash
gene_file=P_aphanis-THeavenSCOTT2020_1-ranked.tsv
taxon_code=RA20
echo -e "name\torthogroup" > temp.tsv
while IFS=$'\t' read -r line; do
    gene=$(echo "$line" | awk '{print $1}')
    gene_search="$taxon_code|$gene"
    orthogroup=$(grep "$gene_search" analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | cut -d ':' -f1)
    if [ $gene != "name" ]; then
    echo -e "$gene\t$orthogroup" >> temp.tsv
    fi
done < $gene_file
join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84,1.85 <(sort $gene_file) <(sort temp.tsv) | sort | uniq > temp2.tsv && mv temp2.tsv $gene_file

gene_file=P_aphanis-THeavenDRCT72020_1-ranked.tsv
taxon_code=ST20
echo -e "name\torthogroup" > temp1.tsv
while IFS=$'\t' read -r line; do
    gene=$(echo "$line" | awk '{print $1}')
    gene_search="$taxon_code|$gene"
    orthogroup=$(grep "$gene_search" analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | cut -d ':' -f1)
    if [ $gene != "name" ]; then
    echo -e "$gene\t$orthogroup" >> temp1.tsv
    fi
done < $gene_file
join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84,1.85 <(sort $gene_file) <(sort temp1.tsv) | sort | uniq > temp12.tsv && mv temp12.tsv $gene_file

gene_file=P_aphanis-THeavenDRCT72021_1-ranked.tsv
taxon_code=ST21
echo -e "name\torthogroup" > temp1.tsv
while IFS=$'\t' read -r line; do
    gene=$(echo "$line" | awk '{print $1}')
    gene_search="$taxon_code|$gene"
    orthogroup=$(grep "$gene_search" analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | cut -d ':' -f1)
    if [ $gene != "name" ]; then
    echo -e "$gene\t$orthogroup" >> temp1.tsv
    fi
done < $gene_file
join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84,1.85 <(sort $gene_file) <(sort temp1.tsv) | sort | uniq > temp12.tsv && mv temp12.tsv $gene_file

gene_file=P_leucotricha-THeavenpOGB2019_1-ranked.tsv
taxon_code=AP19
echo -e "name\torthogroup" > temp.tsv
while IFS=$'\t' read -r line; do
    gene=$(echo "$line" | awk '{print $1}')
    gene_search="$taxon_code|$gene"
    orthogroup=$(grep "$gene_search" analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | cut -d ':' -f1)
    if [ $gene != "name" ]; then
    echo -e "$gene\t$orthogroup" >> temp.tsv
    fi
done < $gene_file
join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84,1.85 <(sort $gene_file) <(sort temp.tsv) | sort | uniq > temp2.tsv && mv temp2.tsv $gene_file

gene_file=P_leucotricha-THeavenp11_1-ranked.tsv
taxon_code=AP20
echo -e "name\torthogroup" > temp1.tsv
while IFS=$'\t' read -r line; do
    gene=$(echo "$line" | awk '{print $1}')
    gene_search="$taxon_code|$gene"
    orthogroup=$(grep "$gene_search" analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | cut -d ':' -f1)
    if [ $gene != "name" ]; then
    echo -e "$gene\t$orthogroup" >> temp1.tsv
    fi
done < $gene_file
join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84,1.85 <(sort $gene_file) <(sort temp1.tsv) | sort | uniq > temp12.tsv && mv temp12.tsv $gene_file

gene_file=P_leucotricha-THeavenpOGB2021_1-ranked.tsv
taxon_code=AP21
echo -e "name\torthogroup" > temp.tsv
while IFS=$'\t' read -r line; do
    gene=$(echo "$line" | awk '{print $1}')
    gene_search="$taxon_code|$gene"
    orthogroup=$(grep "$gene_search" analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | cut -d ':' -f1)
    if [ $gene != "name" ]; then
    echo -e "$gene\t$orthogroup" >> temp.tsv
    fi
done < $gene_file
join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84,1.85 <(sort $gene_file) <(sort temp.tsv) | sort | uniq > temp2.tsv && mv temp2.tsv $gene_file

while IFS=$'\t' read -r line; do
    awk -F'\t' '$27 == "1" && $23 == "." && $3 == "0" {print $2}' <<< "$line" | grep -v 'orthogroup' >> temp.txt
done < <(cat P_leucotricha-THeavenpOGB2021_1-ranked.tsv P_leucotricha-THeavenp11_1-ranked.tsv P_leucotricha-THeavenpOGB2019_1-ranked.tsv P_aphanis-THeavenDRCT72021_1-ranked.tsv P_aphanis-THeavenDRCT72020_1-ranked.tsv P_aphanis-THeavenSCOTT2020_1-ranked.tsv)


gene_file=P_leucotricha-THeavenpOGB2021_1-ranked.tsv
while IFS=$'\t' read -r line; do
    awk -F'\t' '$27 == "1" && $23 == "." && $3 == "0" {print $2}' <<< "$line" | grep -v 'orthogroup' 
done < $gene_file | awk '{count[$1]++} END {for (term in count) print term, count[term]}' > temp1.tsv

gene_file=P_leucotricha-THeavenp11_1-ranked.tsv
while IFS=$'\t' read -r line; do
    awk -F'\t' '$27 == "1" && $23 == "." && $3 == "0" {print $2}' <<< "$line" | grep -v 'orthogroup' 
done < $gene_file | awk '{count[$1]++} END {for (term in count) print term, count[term]}' > temp2.tsv

gene_file=P_leucotricha-THeavenpOGB2019_1-ranked.tsv
while IFS=$'\t' read -r line; do
    awk -F'\t' '$27 == "1" && $23 == "." && $3 == "0" {print $2}' <<< "$line" | grep -v 'orthogroup' 
done < $gene_file | awk '{count[$1]++} END {for (term in count) print term, count[term]}' > temp3.tsv

gene_file=P_aphanis-THeavenDRCT72021_1-ranked.tsv
while IFS=$'\t' read -r line; do
    awk -F'\t' '$27 == "1" && $23 == "." && $3 == "0" {print $2}' <<< "$line" | grep -v 'orthogroup' 
done < $gene_file | awk '{count[$1]++} END {for (term in count) print term, count[term]}' > temp4.tsv

gene_file=P_aphanis-THeavenDRCT72020_1-ranked.tsv
while IFS=$'\t' read -r line; do
    awk -F'\t' '$27 == "1" && $23 == "." && $3 == "0" {print $2}' <<< "$line" | grep -v 'orthogroup' 
done < $gene_file | awk '{count[$1]++} END {for (term in count) print term, count[term]}' > temp5.tsv

gene_file=P_aphanis-THeavenSCOTT2020_1-ranked.tsv
while IFS=$'\t' read -r line; do
    awk -F'\t' '$27 == "1" && $23 == "." && $3 == "0" {print $2}' <<< "$line" | grep -v 'orthogroup' 
done < $gene_file | awk '{count[$1]++} END {for (term in count) print term, count[term]}' > temp6.tsv

tr ' ' '\t' < temp1.tsv > temp_file && mv temp_file temp1.tsv
tr ' ' '\t' < temp2.tsv > temp_file && mv temp_file temp2.tsv
tr ' ' '\t' < temp3.tsv > temp_file && mv temp_file temp3.tsv
tr ' ' '\t' < temp4.tsv > temp_file && mv temp_file temp4.tsv
tr ' ' '\t' < temp5.tsv > temp_file && mv temp_file temp5.tsv
tr ' ' '\t' < temp6.tsv > temp_file && mv temp_file temp6.tsv

cat temp1.tsv temp2.tsv temp3.tsv temp4.tsv temp5.tsv temp6.tsv | awk '{print $1}' | sort | uniq > temp.tsv

join -t $'\t' -a 1 -e 0 -o 0,2.2 <(sort temp.tsv) <(sort temp1.tsv) > tempx.txt && mv tempx.txt temp.tsv
join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2 <(sort temp.tsv) <(sort temp2.tsv) > tempx.txt && mv tempx.txt temp.tsv
join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3 <(sort temp.tsv) <(sort temp3.tsv) > tempx.txt && mv tempx.txt temp.tsv
join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3,1.4 <(sort temp.tsv) <(sort temp4.tsv) > tempx.txt && mv tempx.txt temp.tsv
join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3,1.4,1.5 <(sort temp.tsv) <(sort temp5.tsv) > tempx.txt && mv tempx.txt temp.tsv
join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3,1.4,1.5,1.6 <(sort temp.tsv) <(sort temp6.tsv) > tempx.txt && mv tempx.txt 6_OGs.tsv
```

### Format fasta files
```bash
  ProjDir=/projects/nano_diagnostics
  cd $ProjDir
  IsolateAbrv=All6_isolates
  WorkDir=analysis/orthology/orthofinder/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
  ProgDir=/home/heavet/git_repos/tools/gene_prediction

  ProjDir=/projects/nano_diagnostics
  cd $ProjDir
  IsolateAbrv=Aphanis_isolates
  WorkDir=analysis/orthology/orthofinder/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
  ProgDir=/home/heavet/git_repos/tools/gene_prediction

  ProjDir=/projects/nano_diagnostics
  cd $ProjDir
  IsolateAbrv=Leucotricha_isolates
  WorkDir=analysis/orthology/orthofinder/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
  ProgDir=/home/heavet/git_repos/tools/gene_prediction

  ProjDir=/projects/nano_diagnostics
  cd $ProjDir
  IsolateAbrv=Besthost_isolates
  WorkDir=analysis/orthology/orthofinder/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
  ProgDir=/home/heavet/git_repos/tools/gene_prediction

  ProjDir=/projects/nano_diagnostics
  cd $ProjDir
  IsolateAbrv=Bestspecies_isolates
  WorkDir=analysis/orthology/orthofinder/$IsolateAbrv
  mkdir -p $WorkDir
  mkdir -p $WorkDir/formatted
  mkdir -p $WorkDir/goodProteins
  mkdir -p $WorkDir/badProteins  
  ProgDir=/home/heavet/git_repos/tools/gene_prediction

Taxon_code=ST20
Fasta_file=$(ls /projects/nano_diagnostics/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta)
Id_field=1
$ProgDir/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=ST21
Fasta_file=$(ls /projects/nano_diagnostics/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta)
Id_field=1
$ProgDir/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=RA20
Fasta_file=$(ls /projects/nano_diagnostics/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta)
Id_field=1
$ProgDir/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=AP19
Fasta_file=$(ls /projects/nano_diagnostics/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta)
Id_field=1
$ProgDir/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=AP20
Fasta_file=$(ls /projects/nano_diagnostics/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta)
Id_field=1
$ProgDir/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta

Taxon_code=AP21
Fasta_file=$(ls /projects/nano_diagnostics/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta)
Id_field=1
$ProgDir/orthomclAdjustFasta.pl  $Taxon_code $Fasta_file $Id_field
mv "$Taxon_code".fasta $WorkDir/formatted/"$Taxon_code".fasta
```
For all 6:
OrthoFinder assigned 102145 genes (96.4% of total) to 18296 orthogroups. Fifty percent of all genes were in orthogroups with 6 or more genes (G50 was 6) and were contained in the largest 5190 orthogroups (O50 was 5190). There were 7207 orthogroups with all species present and 4800 of these consisted entirely of single-copy genes.

For the 3 P.aphanis:
OrthoFinder assigned 46655 genes (95.0% of total) to 12465 orthogroups. Fifty percent of all genes were in orthogroups with 3 or more genes (G50 was 3) and were contained in the largest 4271 orthogroups (O50 was 4271). There were 8930 orthogroups with all species present and 6370 of these consisted entirely of single-copy genes.

For the 3 P.leucotricha:
OrthoFinder assigned 55528 genes (97.7% of total) to 16217 orthogroups. Fifty percent of all genes were in orthogroups with 3 or more genes (G50 was 3) and were contained in the largest 6625 orthogroups (O50 was 6625). There were 14129 orthogroups with all species present and 11463 of these consisted entirely of single-copy genes.

For the best aphanis and leucotricha assembly
OrthoFinder assigned 31680 genes (87.4% of total) to 9601 orthogroups. Fifty percent of all genes were in orthogroups with 3 or more genes (G50 was 3) and were contained in the largest 2854 orthogroups (O50 was 2854). There were 8142 orthogroups with all species present and 6069 of these consisted entirely of single-copy genes.

For the best genome from apple, strawberry and raspberry
OrthoFinder assigned 46194 genes (89.8% of total) to 11482 orthogroups. Fifty percent of all genes were in orthogroups with 3 or more genes (G50 was 3) and were contained in the largest 3951 orthogroups (O50 was 3951). There were 7504 orthogroups with all species present and 5537 of these consisted entirely of single-copy genes.

### Filter proteins into good and poor sets
```bash
for WorkDir in $(ls -d analysis/orthology/orthofinder/*); do
ProgDir=/home/heavet/git_repos/tools/gene_prediction
Input_dir=$WorkDir/formatted
Min_length=10
Max_percent_stops=20
Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
$ProgDir/orthomclFilterFasta.pl $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file
done
```
### Orthofinder
```bash
scp -r /projects/nano_diagnostics/analysis/orthology theaven@gruffalo.cropdiversity.ac.uk:/home/theaven/scratch/analysis/.
```
```bash
  ProjDir=/projects/nano_diagnostics
  cd $ProjDir

screen -S orthofinder
srun -p himem -J orthofinder --cpus-per-task=144 --pty bash
conda activate orthofinder2.7
for analysis in $(ls -d analysis/orthology/orthofinder/*); do
cd $analysis
orthofinder -f formatted -t 144 -a 10 2>&1 | tee -a orthofinder_report.txt
cd $ProjDir
done
conda deactivate


```
### Extracting fasta files for orthogroups
```bash
conda activate orthofinder2.7
for WorkDir in $(ls -d analysis/orthology/orthofinder/All6_isolates); do
ProgDir=/home/heavet/git_repos/tools/gene_prediction
OrthogroupsTxt=$(ls $WorkDir/formatted/OrthoFinder/*/Orthogroups/Orthogroups.txt)
GoodProts=$WorkDir/goodProteins/goodProteins.fasta
OutDir=$WorkDir/orthogroups_fasta
mkdir -p $OutDir
python $ProgDir/orthoMCLgroups2fasta.py --orthogroups $OrthogroupsTxt --fasta $GoodProts --out_dir $OutDir
done
conda deactivate
```
```bash
#Discovered duplicate gene should not affect the number of orthogroups? As there is still one legitimate gene of each. However will affect the number of genes in a given orthogroup if not removed:

WorkDir=analysis/orthology/orthofinder/All6_isolates

for Proteome in $(ls /home/theaven/scratch/gene_pred/*/*/predector_singularity2/results/final_genes_appended_renamed.pep/*_secretome.pep.fasta); do
Organism=$(echo $Proteome|cut -f6 -d '/')
Strain=$(echo $Proteome|cut -f7 -d '/')
echo $Orangism - $Strain
    for var1 in $(cat /home/theaven/scratch/gene_pred/${Organism}/${Strain}/codingquarry/rep_modeling/final/${Strain}_duplicated_genes.txt); do
    sed -i "s/ AP19|${var1}.t1//g" analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt
    sed -i "s/ AP20|${var1}.t1//g" analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt
    sed -i "s/ AP21|${var1}.t1//g" analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt
    sed -i "s/ ST20|${var1}.t1//g" analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt
    sed -i "s/ ST21|${var1}.t1//g" analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt
    sed -i "s/ RA20|${var1}.t1//g" analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt  
    done #remove discovered duplicates
done
#From this ^ discovered that genes can be in multiple orthogroups?

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < $WorkDir/goodProteins/goodProteins.fasta > linearized.fasta
for Proteome in $(ls /home/theaven/scratch/gene_pred/*/*/predector_singularity2/results/final_genes_appended_renamed.pep/*_secretome.pep.fasta); do
Organism=$(echo $Proteome|cut -f6 -d '/')
Strain=$(echo $Proteome|cut -f7 -d '/')
echo $Orangism - $Strain
    for var1 in $(cat /home/theaven/scratch/gene_pred/${Organism}/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/${Strain}_duplicated_genes.txt); do
    sed -i "/>RA20|${var1}.t1/d" linearized.fasta
    done #remove discovered duplicates
        for var1 in $(cat /home/theaven/scratch/gene_pred/${Organism}/THeavenDRCT72020_1/codingquarry/rep_modeling/final/${Strain}_duplicated_genes.txt); do
    sed -i "/>ST20|${var1}.t1/d" linearized.fasta
    done #remove discovered duplicates
        for var1 in $(cat /home/theaven/scratch/gene_pred/${Organism}/THeavenDRCT72021_1/codingquarry/rep_modeling/final/${Strain}_duplicated_genes.txt); do
    sed -i "/>ST21|${var1}.t1/d" linearized.fasta
    done #remove discovered duplicates
        for var1 in $(cat /home/theaven/scratch/gene_pred/${Organism}/THeavenp11_1/codingquarry/rep_modeling/final/${Strain}_duplicated_genes.txt); do
    sed -i "/>AP20|${var1}.t1/d" linearized.fasta
    done #remove discovered duplicates
        for var1 in $(cat /home/theaven/scratch/gene_pred/${Organism}/THeavenpOGB2019_1/codingquarry/rep_modeling/final/${Strain}_duplicated_genes.txt); do
    sed -i "/>AP19|${var1}.t1/d" linearized.fasta
    done #remove discovered duplicates
        for var1 in $(cat /home/theaven/scratch/gene_pred/${Organism}/THeavenpOGB2021_1/codingquarry/rep_modeling/final/${Strain}_duplicated_genes.txt); do
    sed -i "/>AP21|${var1}.t1/d" linearized.fasta
    done #remove discovered duplicates
done
tr "\t" "\n" < linearized.fasta | fold -w 60 > goodProteins.fasta 
cp goodProteins.fasta $WorkDir/goodProteins/goodProteins.fasta

WorkDir=analysis/orthology/orthofinder/Aphanis_isolates

for Proteome in $(ls /home/theaven/scratch/gene_pred/*/*/predector_singularity2/results/final_genes_appended_renamed.pep/*_secretome.pep.fasta); do
Organism=$(echo $Proteome|cut -f6 -d '/')
Strain=$(echo $Proteome|cut -f7 -d '/')
echo $Orangism - $Strain
    for var1 in $(cat /home/theaven/scratch/gene_pred/${Organism}/${Strain}/codingquarry/rep_modeling/final/${Strain}_duplicated_genes.txt); do
    sed -i "s/ AP19|${var1}.t1//g" analysis/orthology/orthofinder/Aphanis_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt
    sed -i "s/ AP20|${var1}.t1//g" analysis/orthology/orthofinder/Aphanis_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt
    sed -i "s/ AP21|${var1}.t1//g" analysis/orthology/orthofinder/Aphanis_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt
    sed -i "s/ ST20|${var1}.t1//g" analysis/orthology/orthofinder/Aphanis_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt
    sed -i "s/ ST21|${var1}.t1//g" analysis/orthology/orthofinder/Aphanis_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt
    sed -i "s/ RA20|${var1}.t1//g" analysis/orthology/orthofinder/Aphanis_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt  
    done #remove discovered duplicates
done
```
#### Orthogroups only containing P. leucotricha genes were extracted:
```bash
WorkDir=analysis/orthology/orthofinder/All6_isolates
Orthogroups=$(ls $WorkDir/formatted/OrthoFinder/*/Orthogroups/Orthogroups.txt)
GoodProts=$WorkDir/goodProteins/goodProteins.fasta

for num in 1; do
UniqDir=$WorkDir/leucotricha_unique
Uniq_leucotricha_groups=$UniqDir/leucotricha_unique.txt
mkdir -p $UniqDir
cat $Orthogroups | grep -v -e 'RA20' | grep -v -e 'ST20' | grep -v -e 'ST21' | grep -e 'AP19' | grep -e 'AP20' | grep -e 'AP21'> $Uniq_leucotricha_groups
echo "The number of orthogroups unique to podosphaera leucotricha are:"
cat $Uniq_leucotricha_groups | wc -l
echo "The following number genes from apple 2019 are contained in these orthogorups:"
cat $Uniq_leucotricha_groups | grep -o -e 'AP19' | wc -l
echo "The following number genes from apple 2020 2020 are contained in these orthogorups:"
cat $Uniq_leucotricha_groups | grep -o -e 'AP20' | wc -l
echo "The following number genes from apple 2021 2021 are contained in these orthogorups:"
cat $Uniq_leucotricha_groups | grep -o -e 'AP21' | wc -l
done
```
The number of orthogroups unique to podosphaera leucotricha are:
4405
The following number genes from apple 2019 are contained in these orthogorups:
5046
The following number genes from apple 2020 2020 are contained in these orthogorups:
5066
The following number genes from apple 2021 2021 are contained in these orthogorups:
5056

#### Orthogroups only containing P. aphanis genes were extracted:
```bash
WorkDir=analysis/orthology/orthofinder/All6_isolates
Orthogroups=$(ls $WorkDir/formatted/OrthoFinder/*/Orthogroups/Orthogroups.txt)
GoodProts=$WorkDir/goodProteins/goodProteins.fasta

for num in 1; do
UniqDir=$WorkDir/aphanis_unique
Uniq_aphanis_groups=$UniqDir/aphanis_unique.txt
mkdir -p $UniqDir
cat $Orthogroups | grep -e 'RA20' | grep -e 'ST20' | grep -e 'ST21' | grep -v -e 'AP19' | grep -v -e 'AP20' | grep -v -e 'AP21'> $Uniq_aphanis_groups
echo "The number of orthogroups unique to podosphaera aphanis are:"
cat $Uniq_aphanis_groups | wc -l
echo "The following number genes from raspberry are contained in these orthogorups:"
cat $Uniq_aphanis_groups | grep -o -e 'RA20' | wc -l
echo "The following number genes from Strawberry 2020 are contained in these orthogorups:"
cat $Uniq_aphanis_groups | grep -o -e 'ST20' | wc -l
echo "The following number genes from strawberry 2021 are contained in these orthogorups:"
cat $Uniq_aphanis_groups | grep -o -e 'ST21' | wc -l
done
```
The number of orthogroups unique to podosphaera aphanis are:
1332
The following number genes from raspberry are contained in these orthogorups:
2028
The following number genes from Strawberry 2020 are contained in these orthogorups:
2018
The following number genes from strawberry 2021 are contained in these orthogorups:
1819

#### Orthogroups only containing P. aphanis strawberry genes were extracted:
```bash
WorkDir=analysis/orthology/orthofinder/All6_isolates
Orthogroups=$(ls $WorkDir/formatted/OrthoFinder/*/Orthogroups/Orthogroups.txt)
GoodProts=$WorkDir/goodProteins/goodProteins.fasta

for num in 1; do
UniqDir=$WorkDir/strawberryaphanis_unique
Uniq_strawberryaphanis_groups=$UniqDir/strawberryaphanis_unique.txt
mkdir -p $UniqDir
cat $Orthogroups | grep -v -e 'RA20' | grep -e 'ST20' | grep -e 'ST21' | grep -v -e 'AP19' | grep -v -e 'AP20' | grep -v -e 'AP21'> $Uniq_strawberryaphanis_groups
echo "The number of orthogroups unique to strawberry podosphaera aphanis are:"
cat $Uniq_strawberryaphanis_groups | wc -l
echo "The following number genes from Strawberry 2020 are contained in these orthogorups:"
cat $Uniq_strawberryaphanis_groups | grep -o -e 'ST20' | wc -l
echo "The following number genes from strawberry 2021 are contained in these orthogorups:"
cat $Uniq_strawberryaphanis_groups | grep -o -e 'ST21' | wc -l
done
```
The number of orthogroups unique to strawberry podosphaera aphanis are:
1390
The following number genes from Strawberry 2020 are contained in these orthogorups:
1569
The following number genes from strawberry 2021 are contained in these orthogorups:
1561

#### Orthogroups only containing P. aphanis raspberry genes were extracted:
```bash
WorkDir=analysis/orthology/orthofinder/All6_isolates
Orthogroups=$(ls $WorkDir/formatted/OrthoFinder/*/Orthogroups/Orthogroups.txt)
GoodProts=$WorkDir/goodProteins/goodProteins.fasta

for num in 1; do
UniqDir=$WorkDir/raspberryaphanis_unique
Uniq_raspberryaphanis_groups=$UniqDir/raspberryaphanis_unique.txt
mkdir -p $UniqDir
cat $Orthogroups | grep -e 'RA20' | grep -v -e 'ST20' | grep -v -e 'ST21' | grep -v -e 'AP19' | grep -v -e 'AP20' | grep -v -e 'AP21'> $Uniq_raspberryaphanis_groups
echo "The number of orthogroups unique to raspberry podosphaera aphanis are:"
cat $Uniq_raspberryaphanis_groups | wc -l
echo "The following number genes from raspberry are contained in these orthogorups:"
cat $Uniq_raspberryaphanis_groups | grep -o -e 'RA20' | wc -l
done
```
The number of orthogroups unique to raspberry podosphaera aphanis are:
894
The following number genes from raspberry are contained in these orthogorups:
991

#### Venn diagrams
```bash
conda activate orthofinder2.7
#ProgDir=/home/heavet/git_repos/tools/gene_prediction
ProgDir=/home/theaven/scratch/apps/tools

WorkDir=analysis/orthology/orthofinder/All6_isolates
GoodProts=$WorkDir/goodProteins/goodProteins.fasta
$ProgDir/orthoMCLgroups2tab.py $GoodProts $WorkDir/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt > $WorkDir/formatted/All6_isolates_orthogroups.tab

IsolateAbrv=Aphanis_isolates
WorkDir=analysis/orthology/orthofinder/$IsolateAbrv
GoodProts=$WorkDir/goodProteins/goodProteins.fasta
$ProgDir/orthoMCLgroups2tab.py $GoodProts $WorkDir/formatted/OrthoFinder/Results_Jul20/Orthogroups/Orthogroups.txt > $WorkDir/formatted/${IsolateAbrv}_orthogroups.tab
```
```R
setwd("C:/Users/tchea/Desktop")
install.packages("optparse")
install.packages("colorspace")
install.packages("VennDiagram")
install.packages("grid")
install.packages("ggVennDiagram")
install.packages("ggvenn")
install.packages("devtools")
library("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library("ggvenn")
library("ggVennDiagram")
library(optparse)
library(colorspace)
library(VennDiagram)
require(dendextend)
require(installr)
require(colorspace)

orthotabs <- read.table("Leucotricha_isolates_orthogroups.tab")
df1 <- t(orthotabs)
df <- data.frame(df1)
ap19<-which(df$AP19 == 1)
ap20<-which(df$AP20 == 1)
ap21<-which(df$AP21 == 1)
x <- list(ap19,ap20,ap21)
ggVennDiagram(x)
names(x) <- c("Apple 2019","Apple 2020","Apple 2021")
plot1 <- ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, text_size = 6, set_name_size = 5
)

orthotabs <- read.table("Aphanis_isolates_orthogroups.tab")
df1 <- t(orthotabs)
df2 <- data.frame(df1)
ra20<-which(df2$RA20 == 1)
st20<-which(df2$ST20 == 1)
st21<-which(df2$ST21 == 1)
y <- list(ra20,st20,st21)
ggVennDiagram(y)
names(y) <- c("Raspberry 2020","Strawberry 2020","Strawberry 2021")
plot2 <- ggvenn(
  y, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, text_size = 6, set_name_size = 5 
)

orthotabs <- read.table("All6_isolates_orthogroups.tab")
df1 <- t(orthotabs)
df3 <- data.frame(df1)
APH<-which(df3$RA20 == 1 & df3$ST20 == 1 & df3$ST21 == 1)
LEU<-which(df3$AP20 == 1 & df3$AP21 == 1 & df3$AP19 == 1)
z <- list(APH, LEU)
length(APH)
length(LEU)
ggVennDiagram(z)
names(z) <- c("P. aphanis","P. leucotricha")
plot3 <- ggvenn(
  z, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, text_size = 4, set_name_size = 5, auto_scale = TRUE
)

STR<-which(df3$ST20 == 1 & df3$ST21 == 1)
RAS<-which(df3$RA20 == 1)
a <- list(STR, RAS)
ggVennDiagram(a)
names(a) <- c("Stawberry","Raspberry")
plot4 <- ggvenn(
  a, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 5, auto_scale = TRUE
)

STR20<-which(df3$ST20 == 1)
STR21<-which(df3$ST21 == 1)
b <- list(STR20, STR21)
ggVennDiagram(b)
names(b) <- c("Strawberry 2020","Strawberry 2021")
plot5 <- ggvenn(
  b, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 5, auto_scale = TRUE
)

APP19<-which(df3$AP19 == 1)
APP20<-which(df3$AP20 == 1)
APP21<-which(df3$AP21 == 1)
c <- list(APP20, APP21)
ggVennDiagram(c)
names(c) <- c("Apple 2020","Apple 2021")
plot6 <- ggvenn(
  c, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 5, auto_scale = TRUE
)

d <- list(APP19, APP20)
ggVennDiagram(c)
names(d) <- c("Apple 2019","Apple 2020")
plot7 <- ggvenn(
  d, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 5, auto_scale = TRUE
)

e <- list(APP19, APP21)
ggVennDiagram(e)
names(e) <- c("Apple 2019","Apple 2021")
plot8 <- ggvenn(
  e, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, text_size = 4, set_name_size = 5, auto_scale = TRUE
)

jpeg(file="venn1.jpeg", width=700, height=700)
plot(plot1)
dev.off()
jpeg(file="venn2.jpeg", width=700, height=700)
plot(plot2)
dev.off()
jpeg(file="venn3.jpeg", width=500, height=500)
plot(plot3)
dev.off()
jpeg(file="venn4.jpeg", width=500, height=500)
plot(plot4)
dev.off()
jpeg(file="venn5.jpeg", width=500, height=500)
plot(plot5)
dev.off()
jpeg(file="venn6.jpeg", width=500, height=500)
plot(plot6)
dev.off()
jpeg(file="venn7.jpeg", width=500, height=500)
plot(plot7)
dev.off()
jpeg(file="venn8.jpeg", width=500, height=500)
plot(plot8)
dev.off()
```
### Compare best examples of leucotricha and aphanis:
Extract P. lecuotricha and P. aphanis only gene names for the apple 2020 and strawberry 2020 samples.
```bash
WorkDir=analysis/orthology/orthofinder/All6_isolates
Orthogroups=$(ls $WorkDir/formatted/OrthoFinder/*/Orthogroups/Orthogroups.txt)
GoodProts=$WorkDir/goodProteins/goodProteins.fasta

for num in 1; do
UniqDir=$WorkDir/speciescomp
mkdir -p $UniqDir
cat $Orthogroups | grep -e 'ST20' | grep -e 'ST21' | grep -e 'RA20' |  grep -e 'AP20' |  grep -e 'AP21' |  grep -e 'AP19' > $UniqDir/both.txt
cat $UniqDir/both.txt | wc -l #7,207
cat $Orthogroups | grep -v -e 'ST20' | grep -v -e 'ST21' | grep -v -e 'RA20' |  grep -e 'AP20' |  grep -e 'AP21' |  grep -e 'AP19' > $UniqDir/leucotricha.txt
cat $UniqDir/leucotricha.txt | wc -l #4,478
cat $Orthogroups |  grep -v -e 'AP20' |  grep -v -e 'AP21' |  grep -v -e 'AP19' | grep -e 'ST20' | grep -e 'ST21' | grep -e 'RA20' > $UniqDir/aphanis.txt
cat $UniqDir/aphanis.txt | wc -l #1,361

for gene in $(cat $UniqDir/both.txt | grep -o -P '(?<=ST20).*?(?=.t1)'| sed -e "s@|@@g"); do
echo ${gene}.t1 >>$UniqDir/both_ST20_genes.txt
done
for gene in $(cat $UniqDir/both.txt | grep -o -P '(?<=AP20).*?(?=.t1)'| sed -e "s@|@@g"); do
echo ${gene}.t1 >>$UniqDir/both_AP20_genes.txt
done
for gene in $(cat $UniqDir/aphanis.txt | grep -o -P '(?<=ST20).*?(?=.t1)'| sed -e "s@|@@g"); do
echo ${gene}.t1 >>$UniqDir/aphanis_ST20_genes.txt
done
for gene in $(cat $UniqDir/leucotricha.txt | grep -o -P '(?<=AP20).*?(?=.t1)'| sed -e "s@|@@g"); do
echo ${gene}.t1 >>$UniqDir/leucotricha_AP20_genes.txt
done
done



```
SecretomeP 2.0 for non-canonically secreted proteins:
```bash
mkdir split
cd split
for line in $(cat P_leucotricha-THeavenp11_1-leucotricha-ranked.tsv | grep -P '(?<=g).*?(?=.t1)' | cut -f1); do
echo $line>>temp.txt
grep -A 1 -f temp.txt /home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/THeavenp11_1.pep.fasta | sed -e "s@-@N@g">> P_leucotricha-unique.fasta #10126
rm temp.txt
done

for line in $(cat P_aphanis-THeavenDRCT72020_1-aphanis-ranked.tsv | grep -P '(?<=g).*?(?=.t1)' | cut -f1); do
echo $line>>temp.txt
grep -A 1 -f temp.txt /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/THeavenDRCT72020_1.pep.fasta | sed -e "s@-@N@g">> P_aphanis-unique.fasta #4104
rm temp.txt
done

for line in $(cat P_leucotricha-THeavenp11_1-common-ranked.tsv | grep -P '(?<=g).*?(?=.t1)' | cut -f1); do
echo $line>>temp.txt
grep -A 1 -f temp.txt /home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/THeavenp11_1.pep.fasta | sed -e "s@-@N@g">> P_leucotricha-common.fasta #19114
rm temp.txt
done

for line in $(cat P_aphanis-THeavenDRCT72020_1-common-ranked.tsv | grep -P '(?<=g).*?(?=.t1)' | cut -f1); do
echo $line>>temp.txt
grep -A 1 -f temp.txt /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/THeavenDRCT72020_1.pep.fasta | sed -e "s@-@N@g">> P_aphanis-common.fasta #21298
rm temp.txt
done

awk '/^>/{val=$0;next} length($0)>40{print val ORS $0}' P_leucotricha-unique.fasta > P_leucotricha-unique_40.fasta #9884
awk '/^>/{val=$0;next} length($0)>40{print val ORS $0}' P_aphanis-unique.fasta > P_aphanis-unique_40.fasta #4090
awk '/^>/{val=$0;next} length($0)>40{print val ORS $0}' P_leucotricha-common.fasta > P_leucotricha-common_40.fasta #19102
awk '/^>/{val=$0;next} length($0)>40{print val ORS $0}' P_aphanis-common.fasta > P_aphanis-common_40.fasta #21284

split -l 200 P_leucotricha-unique_40.fasta P_leucotricha-unique_split_
split -l 200 P_leucotricha-common_40.fasta P_leucotricha-common_split_
split -l 200 P_aphanis-unique_40.fasta P_aphanis-unique_split_
split -l 200 P_aphanis-common_40.fasta P_aphanis-common_split_
mkdir x
cp *_split_* x/.
#split files downloaded for submission to SecretomeP 2.0 server

#Secretome1:
screen -S secretomep
srun -J SecretomeP -p long --mem-per-cpu 8G --cpus-per-task 8 --pty bash
conda activate jellyfish
cat /home/theaven/scratch/split/P_leucotricha-unique_40.fasta | sed -e "s@X@@g" > /home/theaven/scratch/split/P_leucotricha-unique_40.fsa
cat /home/theaven/scratch/split/P_leucotricha-common_40.fasta | sed -e "s@X@@g" > /home/theaven/scratch/split/P_leucotricha-common_40.fsa
cat /home/theaven/scratch/split/P_aphanis-unique_40.fasta | sed -e "s@X@@g" > /home/theaven/scratch/split/P_aphanis-unique_40.fsa
cat /home/theaven/scratch/split/P_aphanis-common_40.fasta | sed -e "s@X@@g" > /home/theaven/scratch/split/P_aphanis-common_40.fsa

/home/theaven/scratch/apps/secretomep/secretomep-1.0/secretomep -v -s /home/theaven/scratch/split/P_leucotricha-unique_40.fsa >& /home/theaven/scratch/P_leucotricha-unique_40_secretomep.out
/home/theaven/scratch/apps/secretomep/secretomep-1.0/secretomep -v -s /home/theaven/scratch/split/P_leucotricha-common_40.fsa >& /home/theaven/scratch/P_leucotricha-common_40_secretomep.out
/home/theaven/scratch/apps/secretomep/secretomep-1.0/secretomep -v -s /home/theaven/scratch/split/P_aphanis-unique_40.fsa >& /home/theaven/scratch/P_aphanis-unique_40_secretomep.out
/home/theaven/scratch/apps/secretomep/secretomep-1.0/secretomep -v -s /home/theaven/scratch/split/P_aphanis-common_40.fsa >& /home/theaven/scratch/P_aphanis-common_40_secretomep.out
echo finished
exit
exit

awk '/^>/{val=$0;next} length($0)<4000{print val ORS $0}' P_aphanis-common_split_ah

P_leucotricha-unique.fasta > P_leucotricha-unique_40.fasta
```
BLAST proteins:
```bash
#Extract name of gene from predector file to temp.txt, then get the protein fasta entry for that gene
conda activate blast+
cat P_leucotricha-THeavenp11_1-leucotricha-ranked.tsv | grep -v 'name'| cut -f1>temp.txt
seqtk subseq /home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/THeavenp11_1.pep.fasta temp.txt > P_leucotricha-unique.fasta 

cat P_aphanis-THeavenDRCT72020_1-aphanis-ranked.tsv | grep -v 'name'| cut -f1>temp.txt
seqtk subseq /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/THeavenDRCT72020_1.pep.fasta temp.txt > P_aphanis-unique.fasta 

cat P_leucotricha-THeavenp11_1-common-ranked.tsv | grep -v 'name'| cut -f1>temp.txt
seqtk subseq /home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/THeavenp11_1.pep.fasta temp.txt > P_leucotricha-common.fasta 

cat P_aphanis-THeavenDRCT72020_1-common-ranked.tsv | grep -v 'name'| cut -f1>temp.txt
seqtk subseq /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/THeavenDRCT72020_1.pep.fasta temp.txt > P_aphanis-common.fasta 
rm temp.txt


#Make blast database as in Praz et al., 2017
cat podosporaanserina_GCF_000226545.1_ASM22654v1_protein.faa neurosporacrassa_GCF_000182925.2_NC12_protein.faa > mydb.faa
mkdir db
makeblastdb -in mydb.faa -dbtype prot -title CSEPdb -out db/CSEPdb

#blast for matches
blastp -query P_leucotricha-unique.fasta -db db/CSEPdb -out P_leucotricha-unique_blast.tsv -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq'
blastp -query P_leucotricha-common.fasta -db db/CSEPdb -out P_leucotricha-common_blast.tsv -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq'
blastp -query P_aphanis-unique.fasta -db db/CSEPdb -out P_aphanis-unique_blast.tsv -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq'
blastp -query P_aphanis-common.fasta -db db/CSEPdb -out P_aphanis-common_blast.tsv -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq'

#Extract the gene name from blast results
cat P_leucotricha-unique_blast.tsv | cut -f1 |sort | uniq > P_leucotricha-unique_out_blast.txt
cat P_leucotricha-common_blast.tsv | cut -f1 |sort | uniq  > P_leucotricha-common_out_blast.txt
cat P_aphanis-unique_blast.tsv | cut -f1 |sort | uniq  > P_aphanis-unique_out_blast.txt
cat P_aphanis-common_blast.tsv | cut -f1 |sort | uniq > P_aphanis-common_out_blast.txt
#Add information to predector file
for f in P_leucotricha-unique_out_blast.txt; do sed -i "s/$/\t1/" $f; done
for f in P_leucotricha-common_out_blast.txt; do sed -i "s/$/\t1/" $f; done
for f in P_aphanis-unique_out_blast.txt; do sed -i "s/$/\t1/" $f; done
for f in P_aphanis-common_out_blast.txt; do sed -i "s/$/\t1/" $f; done
printf 'name\torthology match' >> P_leucotricha-unique_out_blast.txt
printf 'name\torthology match' >> P_leucotricha-common_out_blast.txt
printf 'name\torthology match' >> P_aphanis-unique_out_blast.txt
printf 'name\torthology match' >> P_aphanis-common_out_blast.txt
join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84 <(sort P_leucotricha-THeavenp11_1-leucotricha-ranked.tsv) <(sort P_leucotricha-unique_out_blast.txt) > P_leucotricha-THeavenp11_1-leucotricha-ranked2.tsv
join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84 <(sort P_aphanis-THeavenDRCT72020_1-aphanis-ranked.tsv) <(sort P_aphanis-unique_out_blast.txt) > P_aphanis-THeavenDRCT72020_1-aphanis-ranked2.tsv
join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84 <(sort P_leucotricha-THeavenp11_1-common-ranked.tsv) <(sort P_leucotricha-common_out_blast.txt) > P_leucotricha-THeavenp11_1-common-ranked2.tsv
join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84 <(sort P_aphanis-THeavenDRCT72020_1-common-ranked.tsv) <(sort P_aphanis-common_out_blast.txt) > P_aphanis-THeavenDRCT72020_1-common-ranked2.tsv
####################################################################################################################################################################################################

#Make blast database for mildews: B.hordei, B. graminis, B. graminis f.sp. triticale, G. cichoracearum, E. necator - downloaded proteins predictions from NCBI genome assembly projects

cat GCA_000151065.3_ASM15106v3_protein.faa GCA_900519115.1_Bgt_genome_v3.16_protein.faa GCA_905067625.1_Bgtriticale_THUN12_genome_v1_2_protein.faa GCA_024703715.1_EnE101_chr_assembly_v3_protein.faa GCA_003611235.1_ASM361123v1_protein.faa > mydb.faa
mkdir db2
makeblastdb -in mydb.faa -dbtype prot -title blumeriadb -out db2/blumeriadb

#blast for matches
blastp -query P_leucotricha-unique.fasta -db db2/blumeriadb -out temp1.tsv -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq' 
blastp -query P_leucotricha-common.fasta -db db2/blumeriadb -out temp2.tsv -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq' 
blastp -query P_aphanis-unique.fasta -db db2/blumeriadb -out temp3.tsv -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq' 
blastp -query P_aphanis-common.fasta -db db2/blumeriadb -out temp4.tsv -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq' 

#Extract the gene name from blast results
cat temp1.tsv | cut -f1 |sort | uniq  > P_leucotricha-unique_in_blast.txt
cat temp2.tsv | cut -f1 |sort | uniq  > P_leucotricha-common_in_blast.txt
cat temp3.tsv | cut -f1 |sort | uniq  > P_aphanis-unique_in_blast.txt
cat temp4.tsv | cut -f1 |sort | uniq  > P_aphanis-common_in_blast.txt
#Add information to predector file
for f in P_leucotricha-unique_in_blast.txt; do sed -i "s/$/\t1/" $f; done
for f in P_leucotricha-common_in_blast.txt; do sed -i "s/$/\t1/" $f; done
for f in P_aphanis-unique_in_blast.txt; do sed -i "s/$/\t1/" $f; done
for f in P_aphanis-common_in_blast.txt; do sed -i "s/$/\t1/" $f; done
printf 'name\tmildew match' >> P_leucotricha-unique_in_blast.txt
printf 'name\tmildew match' >> P_leucotricha-common_in_blast.txt
printf 'name\tmildew match' >> P_aphanis-unique_in_blast.txt
printf 'name\tmildew match' >> P_aphanis-common_in_blast.txt
join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84,1.85 <(sort P_leucotricha-THeavenp11_1-leucotricha-ranked2.tsv) <(sort P_leucotricha-unique_in_blast.txt) > P_leucotricha-THeavenp11_1-leucotricha-ranked3.tsv
join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84,1.85 <(sort P_aphanis-THeavenDRCT72020_1-aphanis-ranked2.tsv) <(sort P_aphanis-unique_in_blast.txt) > P_aphanis-THeavenDRCT72020_1-aphanis-ranked3.tsv
join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84,1.85 <(sort P_leucotricha-THeavenp11_1-common-ranked2.tsv) <(sort P_leucotricha-common_in_blast.txt) > P_leucotricha-THeavenp11_1-common-ranked3.tsv
join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84,1.85 <(sort P_aphanis-THeavenDRCT72020_1-common-ranked2.tsv) <(sort P_aphanis-common_in_blast.txt) > P_aphanis-THeavenDRCT72020_1-common-ranked3.tsv

for ome in $(ls /home/theaven/scratch/P_leucotricha-unique_40_secretomep.out); do 
sed -i -n '/SecretomeP 1.0g predictions/,$p' $ome 
sed -i '$d' $ome 
sed -i '1d' $ome 
sed -i '1d' $ome 
sed -i '1d' $ome 
sed -i '1d' $ome 
done
for ome in $(ls /home/theaven/scratch/P_leucotricha-common_40_secretomep.out); do 
sed -i -n '/SecretomeP 1.0g predictions/,$p' $ome 
sed -i '$d' $ome 
sed -i '1d' $ome 
sed -i '1d' $ome 
sed -i '1d' $ome 
sed -i '1d' $ome 
done
for ome in $(ls /home/theaven/scratch/P_aphanis-unique_40_secretomep.out); do 
sed -i -n '/SecretomeP 1.0g predictions/,$p' $ome 
sed -i '$d' $ome 
sed -i '1d' $ome 
sed -i '1d' $ome 
sed -i '1d' $ome 
sed -i '1d' $ome 
done
for ome in $(ls /home/theaven/scratch/P_aphanis-common_40_secretomep.out); do 
sed -i -n '/SecretomeP 1.0g predictions/,$p' $ome 
sed -i '$d' $ome 
sed -i '1d' $ome 
sed -i '1d' $ome 
sed -i '1d' $ome 
sed -i '1d' $ome 
done
printf 'name\tNN-score\tOdds\tWeighted\tWarning' >> /home/theaven/scratch/P_leucotricha-unique_40_secretomep.out
printf 'name\tNN-score\tOdds\tWeighted\tWarning' >> /home/theaven/scratch/P_leucotricha-common_40_secretomep.out
printf 'name\tNN-score\tOdds\tWeighted\tWarning' >> /home/theaven/scratch/P_aphanis-unique_40_secretomep.out
printf 'name\tNN-score\tOdds\tWeighted\tWarning' >> /home/theaven/scratch/P_aphanis-common_40_secretomep.out
awk -v OFS='\t' '{$1=$1}1' /home/theaven/scratch/P_leucotricha-unique_40_secretomep.out > /home/theaven/scratch/P_leucotricha-unique_40_secretomep.tsv
awk -v OFS='\t' '{$1=$1}1' /home/theaven/scratch/P_leucotricha-common_40_secretomep.out > /home/theaven/scratch/P_leucotricha-common_40_secretomep.tsv
awk -v OFS='\t' '{$1=$1}1' /home/theaven/scratch/P_aphanis-unique_40_secretomep.out > /home/theaven/scratch/P_aphanis-unique_40_secretomep.tsv
awk -v OFS='\t' '{$1=$1}1' /home/theaven/scratch/P_aphanis-common_40_secretomep.out > /home/theaven/scratch/P_aphanis-common_40_secretomep.tsv
join -t $'\t' -a 1 -e 0 -o 0,2.2,2.3,2.4,2.5,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84,1.85,1.86 <(sort P_leucotricha-THeavenp11_1-leucotricha-ranked3.tsv) <(sort /home/theaven/scratch/P_leucotricha-unique_40_secretomep.tsv) > P_leucotricha-THeavenp11_1-leucotricha-ranked4.tsv
join -t $'\t' -a 1 -e 0 -o 0,2.2,2.3,2.4,2.5,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84,1.85,1.86 <(sort P_aphanis-THeavenDRCT72020_1-aphanis-ranked3.tsv) <(sort /home/theaven/scratch/P_aphanis-common_40_secretomep.tsv) > P_aphanis-THeavenDRCT72020_1-aphanis-ranked4.tsv
join -t $'\t' -a 1 -e 0 -o 0,2.2,2.3,2.4,2.5,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84,1.85,1.86 <(sort P_leucotricha-THeavenp11_1-common-ranked3.tsv) <(sort /home/theaven/scratch/P_leucotricha-unique_40_secretomep.tsv) > P_leucotricha-THeavenp11_1-common-ranked4.tsv
join -t $'\t' -a 1 -e 0 -o 0,2.2,2.3,2.4,2.5,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84,1.85,1.86 <(sort P_aphanis-THeavenDRCT72020_1-common-ranked3.tsv) <(sort /home/theaven/scratch/P_aphanis-common_40_secretomep.tsv) > P_aphanis-THeavenDRCT72020_1-common-ranked4.tsv
```
Y/F/WxC
```bash
for proteome in $(echo GCA_000151065.3_ASM15106v3_protein.fsa); do
grep -B1 'YGC\|YAC\|YLC\|YMC\|YFC\|YWC\|YKC\|YQC\|YEC\|YSC\|YPC\|YVC\|YIC\|YYC\|YHC\|YRC\|YNC\|YDC\|YTC\|FGC\|FAC\|FLC\|FMC\|FFC\|FWC\|FKC\|FQC\|FEC\|FSC\|FPC\|FVC\|FIC\|FYC\|FHC\|FRC\|FNC\|FDC\|FTC\|WGC\|WAC\|WLC\|WMC\|WFC\|WWC\|WKC\|WQC\|WEC\|WSC\|WPC\|WVC\|WIC\|WYC\|WHC\|WRC\|WNC\|WDC\|WTC' $proteome | grep -v 'YGC\|YAC\|YLC\|YMC\|YFC\|YWC\|YKC\|YQC\|YEC\|YSC\|YPC\|YVC\|YIC\|YYC\|YHC\|YRC\|YNC\|YDC\|YTC\|FGC\|FAC\|FLC\|FMC\|FFC\|FWC\|FKC\|FQC\|FEC\|FSC\|FPC\|FVC\|FIC\|FYC\|FHC\|FRC\|FNC\|FDC\|FTC\|WGC\|WAC\|WLC\|WMC\|WFC\|WWC\|WKC\|WQC\|WEC\|WSC\|WPC\|WVC\|WIC\|WYC\|WHC\|WRC\|WNC\|WDC\|WTC\|--' > temp.txt
done
#returns about 1/3 regardless of input species
```
```bash
#Effectors:
#P. aphanis strawberry 2020
for file in $(ls /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-ranked.tsv); do
  head -n 1 /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-ranked.tsv > /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-common-ranked.tsv
  head -n 1 /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-ranked.tsv > /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-aphanis-ranked.tsv
  grep -w -f analysis/orthology/orthofinder/All6_isolates/speciescomp/both_ST20_genes.txt $file >> /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-common-ranked.tsv
  grep -w -f analysis/orthology/orthofinder/All6_isolates/speciescomp/aphanis_ST20_genes.txt $file >> /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-aphanis-ranked.tsv
done

#P. leucotricha apple 2020
for file in $(ls /home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_leucotricha-THeavenp11_1-ranked.tsv); do
  head -n 1 /home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_leucotricha-THeavenp11_1-ranked.tsv > /home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_leucotricha-THeavenp11_1-common-ranked.tsv
  head -n 1 /home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_leucotricha-THeavenp11_1-ranked.tsv > /home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_leucotricha-THeavenp11_1-leucotricha-ranked.tsv
  grep -w -f analysis/orthology/orthofinder/All6_isolates/speciescomp/both_AP20_genes.txt $file >> /home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_leucotricha-THeavenp11_1-common-ranked.tsv
  grep -w -f analysis/orthology/orthofinder/All6_isolates/speciescomp/leucotricha_AP20_genes.txt $file >> /home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_leucotricha-THeavenp11_1-leucotricha-ranked.tsv
done

#Transcription factors:
#P.aphanis
for file in $(ls gene_pred/P_aphanis/THeavenDRCT72020_1/transcription_factors/THeavenDRCT72020_1_TF_gene_headers.txt); do
  aphanis_ST20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/speciescomp/aphanis_ST20_genes.txt)
  both_ST20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/speciescomp/both_ST20_genes.txt)
  aphanis_ST20_genes_out=$(dirname $file)/$(basename $file .txt)_common.txt
  both_ST20_genes_out=$(dirname $file)/$(basename $file .txt)_aphanis.txt
  echo $aphanis_ST20_genes_out
  echo $both_ST20_genes_out
  grep -w -f $aphanis_ST20_genes $file > $aphanis_ST20_genes_out
  grep -w -f $both_ST20_genes $file > $both_ST20_genes_out
done
echo number of transcription factors in THeavenDRCT72020_1 >> TForthologyreport.txt 
echo unique to P.aphanis: >> TForthologyreport.txt  >> TForthologyreport.txt 
cat $aphanis_ST20_genes_out| wc -l >> TForthologyreport.txt 
echo shared with P.leucotricha: >> TForthologyreport.txt 
cat $both_ST20_genes_out| wc -l >> TForthologyreport.txt 
printf '\n'  >> TForthologyreport.txt
#P.leucotricha
for file in $(ls gene_pred/P_leucotricha/THeavenp11_1/transcription_factors/THeavenp11_1_TF_gene_headers.txt); do
  leucotricha_AP20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/speciescomp/leucotricha_AP20_genes.txt)
  both_AP20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/speciescomp/both_AP20_genes.txt)
  leucotricha_AP20_genes_out=$(dirname $file)/$(basename $file .txt)_common.txt
  both_AP20_genes_out=$(dirname $file)/$(basename $file .txt)_leucotricha.txt
  echo $leucotricha_AP20_genes_out
  echo $both_AP20_genes_out
  grep -w -f $leucotricha_AP20_genes $file > $leucotricha_AP20_genes_out
  grep -w -f $both_AP20_genes $file > $both_AP20_genes_out
done
echo number of transcription factors THeavenp11_1 >> TForthologyreport.txt 
echo unique to P.leucotricha: >> TForthologyreport.txt 
cat $leucotricha_AP20_genes_out| wc -l >> TForthologyreport.txt 
echo shared with P.aphanis: >> TForthologyreport.txt 
cat $both_AP20_genes_out| wc -l >> TForthologyreport.txt 
printf '\n'  >> TForthologyreport.txt

#SSCPs
#P.aphanis
for file in $(ls gene_pred/P_aphanis/THeavenDRCT72020_1/sscp/THeavenDRCT72020_1_sscp.fa); do
  aphanis_ST20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/speciescomp/aphanis_ST20_genes.txt)
  both_ST20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/speciescomp/both_ST20_genes.txt)
  aphanis_ST20_genes_out=$(dirname $file)/$(basename $file .txt)_common.txt
  both_ST20_genes_out=$(dirname $file)/$(basename $file .txt)_aphanis.txt
  echo $aphanis_ST20_genes_out
  echo $both_ST20_genes_out
  grep -w -f $aphanis_ST20_genes $file > $aphanis_ST20_genes_out
  grep -w -f $both_ST20_genes $file > $both_ST20_genes_out
done
echo number of SSCPs in THeavenDRCT72020_1 >> SSCPorthologyreport.txt 
echo unique to P.aphanis: >> TForthologyreport.txt  >> SSCPorthologyreport.txt 
cat $aphanis_ST20_genes_out| wc -l >> SSCPorthologyreport.txt 
echo shared with P.leucotricha: >> SSCPorthologyreport.txt 
cat $both_ST20_genes_out| wc -l >> SSCPorthologyreport.txt 
printf '\n'  >> SSCPorthologyreport.txt
#P.leucotricha
for file in $(ls gene_pred/P_leucotricha/THeavenp11_1/sscp/THeavenp11_1_sscp.fa); do
  leucotricha_AP20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/speciescomp/leucotricha_AP20_genes.txt)
  both_AP20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/speciescomp/both_AP20_genes.txt)
  leucotricha_AP20_genes_out=$(dirname $file)/$(basename $file .txt)_common.txt
  both_AP20_genes_out=$(dirname $file)/$(basename $file .txt)_leucotricha.txt
  echo $leucotricha_AP20_genes_out
  echo $both_AP20_genes_out
  grep -w -f $leucotricha_AP20_genes $file > $leucotricha_AP20_genes_out
  grep -w -f $both_AP20_genes $file > $both_AP20_genes_out
done
echo number of SSCPs THeavenp11_1 >> SSCPorthologyreport.txt 
echo unique to P.leucotricha: >> SSCPorthologyreport.txt 
cat $leucotricha_AP20_genes_out| wc -l >> SSCPorthologyreport.txt 
echo shared with P.aphanis: >> SSCPorthologyreport.txt 
cat $both_AP20_genes_out| wc -l >> SSCPorthologyreport.txt 
printf '\n'  >> SSCPorthologyreport.txt

#CAZYs
#P.aphanis
for file in $(ls gene_pred/P_aphanis/THeavenDRCT72020_1/cazy/THeavenDRCT72020_1_CAZY_secreted_headers.txt); do
  aphanis_ST20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/speciescomp/aphanis_ST20_genes.txt)
  both_ST20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/speciescomp/both_ST20_genes.txt)
  aphanis_ST20_genes_out=$(dirname $file)/$(basename $file .txt)_common.txt
  both_ST20_genes_out=$(dirname $file)/$(basename $file .txt)_aphanis.txt
  echo $aphanis_ST20_genes_out
  echo $both_ST20_genes_out
  grep -w -f $aphanis_ST20_genes $file > $aphanis_ST20_genes_out
  grep -w -f $both_ST20_genes $file > $both_ST20_genes_out
done
echo number of secreted CAZYs in THeavenDRCT72020_1 >> CAZYorthologyreport.txt 
echo unique to P.aphanis: >> TForthologyreport.txt  >> CAZYorthologyreport.txt 
cat $aphanis_ST20_genes_out| wc -l >> CAZYorthologyreport.txt 
echo shared with P.leucotricha: >> CAZYorthologyreport.txt 
cat $both_ST20_genes_out| wc -l >> CAZYorthologyreport.txt 
printf '\n'  >> CAZYorthologyreport.txt
#P.leucotricha
for file in $(ls gene_pred/P_leucotricha/THeavenp11_1/cazy/THeavenp11_1_CAZY_secreted_headers.txt); do
  leucotricha_AP20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/speciescomp/leucotricha_AP20_genes.txt)
  both_AP20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/speciescomp/both_AP20_genes.txt)
  leucotricha_AP20_genes_out=$(dirname $file)/$(basename $file .txt)_common.txt
  both_AP20_genes_out=$(dirname $file)/$(basename $file .txt)_leucotricha.txt
  echo $leucotricha_AP20_genes_out
  echo $both_AP20_genes_out
  grep -w -f $leucotricha_AP20_genes $file > $leucotricha_AP20_genes_out
  grep -w -f $both_AP20_genes $file > $both_AP20_genes_out
done
echo number of secreted CAZYs in THeavenp11_1 >> CAZYorthologyreport.txt 
echo unique to P.leucotricha: >> CAZYorthologyreport.txt 
cat $leucotricha_AP20_genes_out| wc -l >> CAZYorthologyreport.txt 
echo shared with P.aphanis: >> CAZYorthologyreport.txt 
cat $both_AP20_genes_out| wc -l >> CAZYorthologyreport.txt 
printf '\n'  >> CAZYorthologyreport.txt

#Interpro
#P.aphanis
for file in $(ls gene_pred/P_aphanis/THeavenDRCT72020_1/interproscan/NRI/P.aphanis_strawberry2020_interproscan.tsv); do
  aphanis_ST20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/speciescomp/aphanis_ST20_genes.txt)
  both_ST20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/speciescomp/both_ST20_genes.txt)
  aphanis_ST20_genes_out=$(dirname $file)/$(basename $file .tsv)_common.tsv
  both_ST20_genes_out=$(dirname $file)/$(basename $file .tsv)_aphanis.tsv
  echo $aphanis_ST20_genes_out
  echo $both_ST20_genes_out
  grep -w -f $aphanis_ST20_genes $file > $aphanis_ST20_genes_out
  grep -w -f $both_ST20_genes $file > $both_ST20_genes_out
done

#P.leucotricha
for file in $(ls gene_pred/P_leucotricha/THeavenp11_1/interproscan/NRI/P.leucotricha_apple2020_interproscan.tsv); do
  leucotricha_AP20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/speciescomp/leucotricha_AP20_genes.txt)
  both_AP20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/speciescomp/both_AP20_genes.txt)
  leucotricha_AP20_genes_out=$(dirname $file)/$(basename $file .tsv)_common.tsv
  both_AP20_genes_out=$(dirname $file)/$(basename $file .tsv)_leucotricha.tsv
  echo $leucotricha_AP20_genes_out
  echo $both_AP20_genes_out
  grep -w -f $leucotricha_AP20_genes $file > $leucotricha_AP20_genes_out
  grep -w -f $both_AP20_genes $file > $both_AP20_genes_out
done

#Effector interproscan results:
#P.aphanis
for file in $(ls gene_pred/P_aphanis/THeavenDRCT72020_1/interproscan/NRI/P.aphanis_strawberry2020_interproscan.tsv); do
  echo $file
  cat /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-aphanis-ranked.tsv | cut -f1 > x.txt
grep -w -f x.txt $file > $(dirname $file)/$(basename $file .tsv)_aphanis_effectors.tsv
done

#P.leucotricha
for file in $(ls gene_pred/P_leucotricha/THeavenp11_1/interproscan/NRI/P.leucotricha_apple2020_interproscan.tsv); do
  echo $file
  cat /home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_leucotricha-THeavenp11_1-leucotricha-ranked.tsv | cut -f1 > x.txt
grep -w -f x.txt $file > $(dirname $file)/$(basename $file .tsv)_leucotricha_effectors.tsv
done
```
```
gene_pred/P_aphanis/THeavenDRCT72020_1/interproscan/NRI/P.aphanis_strawberry2020_interproscan.tsv
gene_pred/P_aphanis/THeavenDRCT72021_1/interproscan/NRI/P.aphanis_strawberry2021_interproscan.tsv
gene_pred/P_aphanis/THeavenSCOTT2020_1/interproscan/NRI/P.aphanis_raspberry2020_interproscan.tsv
gene_pred/P_leucotricha/THeavenp11_1/interproscan/NRI/P.leucotricha_apple2020_interproscan.tsv
gene_pred/P_leucotricha/THeavenpOGB2019_1/interproscan/NRI/P.leucotricha_apple2019_interproscan.tsv
gene_pred/P_leucotricha/THeavenpOGB2021_1/interproscan/NRI/P.leucotricha_apple2021_interproscan.tsv
```
### Compare best of apple to best of strawberry only:
Extract Strawberry and Apple only gene names for the 2020 samples.
```bash
WorkDir=analysis/orthology/orthofinder/All6_isolates
Orthogroups=$(ls $WorkDir/formatted/OrthoFinder/*/Orthogroups/Orthogroups.txt)
GoodProts=$WorkDir/goodProteins/goodProteins.fasta

for num in 1; do
UniqDir=$WorkDir/minusrasphostcomp
mkdir -p $UniqDir
cat $Orthogroups | grep -e 'ST20' | grep -e 'ST21' |  grep -e 'AP20' |  grep -e 'AP21' |  grep -e 'AP19' > $UniqDir/minusraspboth.txt
cat $UniqDir/minusraspboth.txt | wc -l #7,416
cat $Orthogroups | grep -v -e 'ST20' | grep -v -e 'ST21' | grep -e 'AP20' |  grep -e 'AP21' |  grep -e 'AP19' > $UniqDir/minusraspapple.txt
cat $UniqDir/minusraspapple.txt | wc -l #4,666
cat $Orthogroups |  grep -e 'ST20' |  grep -e 'ST21' |  grep -v -e 'AP20' |  grep -v -e 'AP21' |  grep -v -e 'AP19' > $UniqDir/minusraspstrawb.txt
cat $UniqDir/minusraspstrawb.txt | wc -l #2,722

for gene in $(cat $UniqDir/minusraspboth.txt | grep -o -P '(?<=ST20).*?(?=.t1)'| sed -e "s@|@@g"); do
echo ${gene}.t1 >>$UniqDir/minusraspboth_ST20_genes.txt
done
for gene in $(cat $UniqDir/minusraspboth.txt | grep -o -P '(?<=AP20).*?(?=.t1)'| sed -e "s@|@@g"); do
echo ${gene}.t1 >>$UniqDir/minusraspboth_AP20_genes.txt
done
for gene in $(cat $UniqDir/minusraspstrawb.txt | grep -o -P '(?<=ST20).*?(?=.t1)'| sed -e "s@|@@g"); do
echo ${gene}.t1 >>$UniqDir/minusraspstrawb_ST20_genes.txt
done
for gene in $(cat $UniqDir/minusraspapple.txt | grep -o -P '(?<=AP20).*?(?=.t1)'| sed -e "s@|@@g"); do
echo ${gene}.t1 >>$UniqDir/minusraspapple_AP20_genes.txt
done
done
```
```bash
#Effectors:
#P. aphanis strawberry 2020
for file in $(ls /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-ranked.tsv); do
  head -n 1 /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-ranked.tsv > /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-minusraspcommon-ranked.tsv
  head -n 1 /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-ranked.tsv > /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-minusraspaphanis-ranked.tsv
  grep -w -f analysis/orthology/orthofinder/All6_isolates/minusrasphostcomp/minusraspboth_ST20_genes.txt $file >> /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-minusraspcommon-ranked.tsv
  grep -w -f analysis/orthology/orthofinder/All6_isolates/minusrasphostcomp/minusraspstrawb_ST20_genes.txt $file >> /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-minusraspaphanis-ranked.tsv
done

#P. leucotricha apple 2020
for file in $(ls /home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_leucotricha-THeavenp11_1-ranked.tsv); do
  head -n 1 /home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_leucotricha-THeavenp11_1-ranked.tsv > /home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_leucotricha-THeavenp11_1-minusraspcommon-ranked.tsv
  head -n 1 /home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_leucotricha-THeavenp11_1-ranked.tsv > /home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_leucotricha-THeavenp11_1-minusraspleucotricha-ranked.tsv
  grep -w -f analysis/orthology/orthofinder/All6_isolates/minusrasphostcomp/minusraspboth_AP20_genes.txt $file >> /home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_leucotricha-THeavenp11_1-minusraspcommon-ranked.tsv
  grep -w -f analysis/orthology/orthofinder/All6_isolates/minusrasphostcomp/minusraspapple_AP20_genes.txt $file >> /home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_leucotricha-THeavenp11_1-minusraspleucotricha-ranked.tsv
done
```

### Compare best examples of raspberry and strawberry:
Extract Strawberry and Raspberry only gene names for the 2020 samples.
```bash
WorkDir=analysis/orthology/orthofinder/All6_isolates
Orthogroups=$(ls $WorkDir/formatted/OrthoFinder/*/Orthogroups/Orthogroups.txt)
GoodProts=$WorkDir/goodProteins/goodProteins.fasta

for num in 1; do
UniqDir=$WorkDir/aphanishostcomp
mkdir -p $UniqDir
cat $Orthogroups | grep -e 'ST20' | grep -e 'ST21' |  grep -e 'RA20' |  grep -v -e 'AP20' |  grep -v -e 'AP21' |  grep -v -e 'AP19' > $UniqDir/both.txt
cat $UniqDir/both.txt | wc -l #
cat $Orthogroups | grep -v -e 'ST20' | grep -v -e 'ST21' | grep -e 'RA20' |  grep -v -e 'AP20' |  grep -v -e 'AP21' |  grep -v -e 'AP19' > $UniqDir/rasp.txt
cat $UniqDir/rasp.txt | wc -l #
cat $Orthogroups |  grep -v -e 'RA20' |  grep -e 'ST20' |  grep -e 'ST21' |  grep -v -e 'AP20' |  grep -v -e 'AP21' |  grep -v -e 'AP19' > $UniqDir/strawb.txt
cat $UniqDir/strawb.txt | wc -l #

for gene in $(cat $UniqDir/both.txt | grep -o -P '(?<=ST20).*?(?=.t1)'| sed -e "s@|@@g"); do
echo ${gene}.t1 >>$UniqDir/both_ST20_genes.txt
done
for gene in $(cat $UniqDir/both.txt | grep -o -P '(?<=RA20).*?(?=.t1)'| sed -e "s@|@@g"); do
echo ${gene}.t1 >>$UniqDir/both_RA20_genes.txt
done
for gene in $(cat $UniqDir/strawb.txt | grep -o -P '(?<=ST20).*?(?=.t1)'| sed -e "s@|@@g"); do
echo ${gene}.t1 >>$UniqDir/strawb_ST20_genes.txt
done
for gene in $(cat $UniqDir/rasp.txt | grep -o -P '(?<=RA20).*?(?=.t1)'| sed -e "s@|@@g"); do
echo ${gene}.t1 >>$UniqDir/rasp_RA20_genes.txt
done
done
```
```bash
#Effectors:
#P. aphanis strawberry 2020
for file in $(ls /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-ranked.tsv); do
  head -n 1 /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-ranked.tsv > /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-aphaniscommon-ranked.tsv
  head -n 1 /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-ranked.tsv > /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-strawberryaphanis-ranked.tsv
  grep -w -f analysis/orthology/orthofinder/All6_isolates/aphanishostcomp/both_ST20_genes.txt $file >> /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-aphaniscommon-ranked.tsv
  grep -w -f analysis/orthology/orthofinder/All6_isolates/aphanishostcomp/strawb_ST20_genes.txt $file >> /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-strawberryaphanis-ranked.tsv
done

#P. aphanis raspberry 2020
for file in $(ls /home/theaven/scratch/gene_pred/P_aphanis/THeavenSCOTT2020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenSCOTT2020_1-ranked.tsv); do
  head -n 1 /home/theaven/scratch/gene_pred/P_aphanis/THeavenSCOTT2020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenSCOTT2020_1-ranked.tsv > /home/theaven/scratch/gene_pred/P_aphanis/THeavenSCOTT2020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenSCOTT2020_1-aphaniscommon-ranked.tsv
  head -n 1 /home/theaven/scratch/gene_pred/P_aphanis/THeavenSCOTT2020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenSCOTT2020_1-ranked.tsv > /home/theaven/scratch/gene_pred/P_aphanis/THeavenSCOTT2020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenSCOTT2020_1-raspberryaphanis-ranked.tsv
  grep -w -f analysis/orthology/orthofinder/All6_isolates/aphanishostcomp/both_RA20_genes.txt $file >> /home/theaven/scratch/gene_pred/P_aphanis/THeavenSCOTT2020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenSCOTT2020_1-aphaniscommon-ranked.tsv
  grep -w -f analysis/orthology/orthofinder/All6_isolates/aphanishostcomp/rasp_RA20_genes.txt $file >> /home/theaven/scratch/gene_pred/P_aphanis/THeavenSCOTT2020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenSCOTT2020_1-raspberryaphanis-ranked.tsv
done

#Transcription factors:
#Strawberry
for file in $(ls gene_pred/P_aphanis/THeavenDRCT72020_1/transcription_factors/THeavenDRCT72020_1_TF_gene_headers.txt); do
  strawb_ST20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/aphanishostcomp/strawb_ST20_genes.txt)
  both_ST20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/aphanishostcomp/both_ST20_genes.txt)
  strawb_ST20_genes_out=$(dirname $file)/$(basename $file .txt)_common.txt
  both_ST20_genes_out=$(dirname $file)/$(basename $file .txt)_aphanis.txt
  echo $strawb_ST20_genes_out
  echo $both_ST20_genes_out
  grep -w -f $strawb_ST20_genes $file > $strawb_ST20_genes_out
  grep -w -f $both_ST20_genes $file > $both_ST20_genes_out
done
echo number of transcription factors in THeavenDRCT72020_1 >> TForthologyreport.txt 
echo unique to strawberry: >> TForthologyreport.txt 
cat $aphanis_ST20_genes_out| wc -l >> TForthologyreport.txt 
echo shared with raspberry: >> TForthologyreport.txt 
cat $both_ST20_genes_out| wc -l >> TForthologyreport.txt
printf '\n'  >> TForthologyreport.txt
#Raspberry
for file in $(ls gene_pred/P_aphanis/THeavenSCOTT2020_1/transcription_factors/THeavenSCOTT2020_1_TF_gene_headers.txt); do
  rasp_RA20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/aphanishostcomp/rasp_RA20_genes.txt)
  both_RA20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/aphanishostcomp/both_RA20_genes.txt)
  rasp_RA20_genes_out=$(dirname $file)/$(basename $file .txt)_common.txt
  both_RA20_genes_out=$(dirname $file)/$(basename $file .txt)_aphanis.txt
  echo $rasp_RA20_genes_out
  echo $both_RA20_genes_out
  grep -w -f $rasp_RA20_genes $file > $rasp_RA20_genes_out
  grep -w -f $both_RA20_genes $file > $both_RA20_genes_out
done
echo number of transcription factors in THeavenSCOTT2020_1 >> TForthologyreport.txt 
echo unique to raspberry: >> TForthologyreport.txt 
cat $rasp_RA20_genes_out| wc -l >> TForthologyreport.txt 
echo shared with strawberry: >> TForthologyreport.txt 
cat $both_RA20_genes_out| wc -l >> TForthologyreport.txt
printf '\n'  >> TForthologyreport.txt

#SSCPs:
#Strawberry
for file in $(ls gene_pred/P_aphanis/THeavenDRCT72020_1/sscp/THeavenDRCT72020_1_sscp.fa); do
  strawb_ST20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/aphanishostcomp/strawb_ST20_genes.txt)
  both_ST20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/aphanishostcomp/both_ST20_genes.txt)
  strawb_ST20_genes_out=$(dirname $file)/$(basename $file .txt)_common.txt
  both_ST20_genes_out=$(dirname $file)/$(basename $file .txt)_aphanis.txt
  echo $strawb_ST20_genes_out
  echo $both_ST20_genes_out
  grep -w -f $strawb_ST20_genes $file > $strawb_ST20_genes_out
  grep -w -f $both_ST20_genes $file > $both_ST20_genes_out
done
echo number of SSCPs in THeavenDRCT72020_1 >> SSCPorthologyreport.txt 
echo unique to strawberry: >> TForthologyreport.txt  >> SSCPorthologyreport.txt 
cat $aphanis_ST20_genes_out| wc -l >> SSCPorthologyreport.txt 
echo shared with raspberry: >> SSCPorthologyreport.txt 
cat $both_ST20_genes_out| wc -l >> SSCPorthologyreport.txt
printf '\n'  >> SSCPorthologyreport.txt
#Raspberry
for file in $(ls gene_pred/P_aphanis/THeavenSCOTT2020_1/sscp/THeavenSCOTT2020_1_sscp.fa); do
  rasp_RA20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/aphanishostcomp/rasp_RA20_genes.txt)
  both_RA20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/aphanishostcomp/both_RA20_genes.txt)
  rasp_RA20_genes_out=$(dirname $file)/$(basename $file .txt)_common.txt
  both_RA20_genes_out=$(dirname $file)/$(basename $file .txt)_aphanis.txt
  echo $rasp_RA20_genes_out
  echo $both_RA20_genes_out
  grep -w -f $rasp_RA20_genes $file > $rasp_RA20_genes_out
  grep -w -f $both_RA20_genes $file > $both_RA20_genes_out
done
echo number of SSCPs in THeavenSCOTT2020_1 >> SSCPorthologyreport.txt 
echo unique to raspberry: >> SSCPorthologyreport.txt 
cat $rasp_RA20_genes_out| wc -l >> SSCPorthologyreport.txt 
echo shared with strawberry: >> SSCPorthologyreport.txt 
cat $both_RA20_genes_out| wc -l >> SSCPorthologyreport.txt
printf '\n'  >> SSCPorthologyreport.txt

#CAZYs:
#Strawberry
for file in $(ls gene_pred/P_aphanis/THeavenDRCT72020_1/cazy/THeavenDRCT72020_1_CAZY_secreted_headers.txt); do
  strawb_ST20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/aphanishostcomp/strawb_ST20_genes.txt)
  both_ST20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/aphanishostcomp/both_ST20_genes.txt)
  strawb_ST20_genes_out=$(dirname $file)/$(basename $file .txt)_common.txt
  both_ST20_genes_out=$(dirname $file)/$(basename $file .txt)_aphanis.txt
  echo $strawb_ST20_genes_out
  echo $both_ST20_genes_out
  grep -w -f $strawb_ST20_genes $file > $strawb_ST20_genes_out
  grep -w -f $both_ST20_genes $file > $both_ST20_genes_out
done
echo number of secreted CAZYs in THeavenDRCT72020_1 >> CAZYorthologyreport.txt 
echo unique to strawberry: >> TForthologyreport.txt  >> CAZYorthologyreport.txt 
cat $aphanis_ST20_genes_out| wc -l >> CAZYorthologyreport.txt 
echo shared with raspberry: >> CAZYorthologyreport.txt 
cat $both_ST20_genes_out| wc -l >> CAZYorthologyreport.txt
printf '\n'  >> CAZYorthologyreport.txt
#Raspberry
for file in $(ls gene_pred/P_aphanis/THeavenSCOTT2020_1/cazy/THeavenSCOTT2020_1_CAZY_secreted_headers.txt); do
  rasp_RA20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/aphanishostcomp/rasp_RA20_genes.txt)
  both_RA20_genes=$(echo analysis/orthology/orthofinder/All6_isolates/aphanishostcomp/both_RA20_genes.txt)
  rasp_RA20_genes_out=$(dirname $file)/$(basename $file .txt)_common.txt
  both_RA20_genes_out=$(dirname $file)/$(basename $file .txt)_aphanis.txt
  echo $rasp_RA20_genes_out
  echo $both_RA20_genes_out
  grep -w -f $rasp_RA20_genes $file > $rasp_RA20_genes_out
  grep -w -f $both_RA20_genes $file > $both_RA20_genes_out
done
echo number of secreted CAZYs in THeavenSCOTT2020_1 >> CAZYorthologyreport.txt 
echo unique to raspberry: >> CAZYorthologyreport.txt 
cat $rasp_RA20_genes_out| wc -l >> CAZYorthologyreport.txt 
echo shared with strawberry: >> CAZYorthologyreport.txt 
cat $both_RA20_genes_out| wc -l >> CAZYorthologyreport.txt
printf '\n'  >> CAZYorthologyreport.txt
```
#### Find the no. of orthogroups for gene features
```bash
ls /home/theaven/projects/niab/theaven/results/

Predictedgenes=$(cut -f1 P_aphanis-THeavenDRCT72020_1-aphanis-ranked4.tsv| grep -v 'name' | wc -l)
pfamID=$(awk '$18 != "." {print $1}' P_aphanis-THeavenDRCT72020_1-aphanis-ranked4.tsv| grep -v 'name' | wc -l)
pfamvirulence=$(awk '$20 == "1" {print $1}' P_aphanis-THeavenDRCT72020_1-aphanis-ranked4.tsv| grep -v 'name' | wc -l)
PHIbaseeffector=$(awk '$15 == "1" {print $1}' P_aphanis-THeavenDRCT72020_1-aphanis-ranked4.tsv| grep -v 'name' | wc -l)
effectorortholog=$(awk '$11 != "." {print $1}' P_aphanis-THeavenDRCT72020_1-aphanis-ranked4.tsv| grep -v 'name' | wc -l)
canonicalsp=$(awk '$31 == "1" {print $1}' P_aphanis-THeavenDRCT72020_1-aphanis-ranked4.tsv| grep -v 'name' | wc -l)
non-canonical=$(awk '$31 != "1" && $5 != 'signal' && $2 > 0.6 {print $1}' P_aphanis-THeavenDRCT72020_1-aphanis-ranked4.tsv| grep -v 'name' | wc -l)
sscp=$(awk '$31 == "1" && $36 <= 300 && $39 > 3 {print $1}' P_aphanis-THeavenDRCT72020_1-aphanis-ranked4.tsv| grep -v 'name' | wc -l)
 csep=$(awk '$31 == "1" && $7 == "0" {print $1}' $Annotations| grep -v 'name' | wc -l)
 secretedeffector=$(awk '$31 == "1" && $27 != "." {print $1}' $Annotations| grep -v 'name' | wc -l)
 secretedcazy=$(awk '$31 == "1" && $21 != "." {print $1}' $Annotations| grep -v 'name' | wc -l)
 
 echo Species+Assembly Predicted_genes pfamID pfam_virulence PHIbase_effector effector_ortholog canonical_SP noncanonical_SP SSCP CSEP Secreted_effectorP3 Secreted_CAZY newSSCP newCSEP newSecreted_effectorP3 > OwnEffectorPredictionReport.txt
 for file in $(ls *-ranked4.tsv); do
 ID=$(echo $file | cut -d '-' -f2,3)
 echo $ID
 Annotations=$file
  echo $Annotations
 Predictedgenes=$(cut -f1 $Annotations| grep -v 'name' | wc -l)
 pfamID=$(awk '$18 != "." {print $1}' $Annotations| grep -v 'name' | wc -l)
 pfamvirulence=$(awk '$20 == "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 PHIbaseeffector=$(awk '$15 == "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 effectorortholog=$(awk '$11 != "." {print $1}' $Annotations| grep -v 'name' | wc -l)
 canonicalsp=$(awk '$31 == "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 noncanonical=$(awk '$31 != "1" && $5 != 'signal' && $2 > 0.6 {print $1}' $Annotations| grep -v 'name' | wc -l)
 sscp=$(awk '$31 == "1" && $36 <= 300 && $39 > 3 {print $1}' $Annotations| grep -v 'name' | wc -l)
 csep=$(awk '$31 == "1" && $7 == "0" {print $1}' $Annotations| grep -v 'name' | wc -l)
 secretedeffector=$(awk '$31 == "1" && $27 == "." && $7 == "0" {print $1}' $Annotations| grep -v 'name' | wc -l)
 secretedcazy=$(awk '$31 == "1" && $21 != "." {print $1}' $Annotations| grep -v 'name' | wc -l)
 newsscp=$(awk '$31 == "1" && $36 <= 300 && $39 > 3 && $11 == "." && $15 != "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 newcsep=$(awk '$31 == "1" && $7 == "0" && $11 == "." && $15 != "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 newsecretedeffector=$(awk '$31 == "1" && $27 == "." && $11 == "." && $15 != "1" && $7 == "0" {print $1}'  $Annotations| grep -v 'name' | wc -l)
 echo $ID $Predictedgenes $pfamID $pfamvirulence $PHIbaseeffector $effectorortholog $canonicalsp $noncanonical $sscp $csep $secretedeffector $secretedcazy $newsscp $newcsep $newsecretedeffector >> OwnEffectorPredictionReport.txt
 done






echo Orthogroups >> OwnEffectorPredictionReport.txt
#Predicted genes=
for Annotations in $(ls P_aphanis-THeavenDRCT72020_1-*-ranked4.tsv); do
ID=$(echo $Annotations | cut -d '-' -f2,3)
cut -f1 $Annotations| grep -v 'name' > ob1.txt
sed -i -e 's/^/ST20|/' ob1.txt
x1=$(grep -f ob1.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#TFs
#pfamID=
awk '$18 != "." {print $1}' $Annotations| grep -v 'name' > ob1x.txt
sed -i -e 's/^/ST20|/' ob1x.txt
x2=$(grep -f ob1x.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#pfamvirulence=
awk '$20 == "1" {print $1}' $Annotations| grep -v 'name' > ob1xx.txt
sed -i -e 's/^/ST20|/' ob1xx.txt
x3=$(grep -f ob1xx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#PHIbaseeffector=
awk '$15 == "1" {print $1}' $Annotations| grep -v 'name' > ob1xxx.txt
sed -i -e 's/^/ST20|/' ob1xxx.txt
x4=$(grep -f ob1xxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#effectorortholog=
awk '$11 != "." {print $1}' $Annotations| grep -v 'name' > ob1xxxx.txt
sed -i -e 's/^/ST20|/' ob1xxxx.txt
x5=$(grep -f ob1xxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#canonicalsp=
awk '$31 == "1" {print $1}' $Annotations| grep -v 'name' > ob1xxxxx.txt
sed -i -e 's/^/ST20|/' ob1xxxxx.txt
x6=$(grep -f ob1xxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#non-canonical=
awk '$31 != "1" && $5 != 'signal' && $2 > 0.6 {print $1}' $Annotations| grep -v 'name' > ob1xxxxxx.txt
sed -i -e 's/^/ST20|/' ob1xxxxxx.txt
x7=$(grep -f ob1xxxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#sscp=
awk '$31 == "1" && $36 <= 300 && $39 > 3 {print $1}' $Annotations| grep -v 'name' > ob1xxxxxxx.txt
sed -i -e 's/^/ST20|/' ob1xxxxxxx.txt
x8=$(grep -f ob1xxxxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#csep=
awk '$31 == "1" && $7 == "0" {print $1}' $Annotations| grep -v 'name' > ob1xxxxxxxx.txt
sed -i -e 's/^/ST20|/' ob1xxxxxxxx.txt
x9=$(grep -f ob1xxxxxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#secretedeffector=
awk '$31 == "1" && $27 != "." {print $1}' $Annotations| grep -v 'name' > ob1xxxxxxxxx.txt
sed -i -e 's/^/ST20|/' ob1xxxxxxxxx.txt
x10=$(grep -f ob1xxxxxxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#secretedcazy=
awk '$31 == "1" && $21 != "." {print $1}' $Annotations| grep -v 'name' > ob1xxxxxxxxxx.txt
sed -i -e 's/^/ST20|/' ob1xxxxxxxxxx.txt
x11=$(grep -f ob1xxxxxxxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
echo $ID $x1 $x2 $x3 $x4 $x5 $x6 $x7 $x8 $x9 $x10 $x11 >> OwnEffectorPredictionReport.txt
done

for Annotations in $(ls P_leucotricha-THeavenp11_1*-ranked4.tsv); do
ID=$(echo $Annotations | cut -d '-' -f2,3)
cut -f1 $Annotations| grep -v 'name' > ob1.txt
sed -i -e 's/^/AP20|/' ob1.txt
x1=$(grep -f ob1.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#TFs
#pfamID=
awk '$18 != "." {print $1}' $Annotations| grep -v 'name' > ob1x.txt
sed -i -e 's/^/AP20|/' ob1x.txt
x2=$(grep -f ob1x.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#pfamvirulence=
awk '$20 == "1" {print $1}' $Annotations| grep -v 'name' > ob1xx.txt
sed -i -e 's/^/AP20|/' ob1xx.txt
x3=$(grep -f ob1xx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#PHIbaseeffector=
awk '$15 == "1" {print $1}' $Annotations| grep -v 'name' > ob1xxx.txt
sed -i -e 's/^/AP20|/' ob1xxx.txt
x4=$(grep -f ob1xxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#effectorortholog=
awk '$11 != "." {print $1}' $Annotations| grep -v 'name' > ob1xxxx.txt
sed -i -e 's/^/AP20|/' ob1xxxx.txt
x5=$(grep -f ob1xxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#canonicalsp=
awk '$31 == "1" {print $1}' $Annotations| grep -v 'name' > ob1xxxxx.txt
sed -i -e 's/^/AP20|/' ob1xxxxx.txt
x6=$(grep -f ob1xxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#non-canonical=
awk '$31 != "1" && $5 != 'signal' && $2 > 0.6 {print $1}' $Annotations| grep -v 'name' > ob1xxxxxx.txt
sed -i -e 's/^/AP20|/' ob1xxxxxx.txt
x7=$(grep -f ob1xxxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#sscp=
awk '$31 == "1" && $36 <= 300 && $39 > 3 {print $1}' $Annotations| grep -v 'name' > ob1xxxxxxx.txt
sed -i -e 's/^/AP20|/' ob1xxxxxxx.txt
x8=$(grep -f ob1xxxxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#csep=
awk '$31 == "1" && $7 == "0" {print $1}' $Annotations| grep -v 'name' > ob1xxxxxxxx.txt
sed -i -e 's/^/AP20|/' ob1xxxxxxxx.txt
x9=$(grep -f ob1xxxxxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#secretedeffector=
awk '$31 == "1" && $27 != "." {print $1}' $Annotations| grep -v 'name' > ob1xxxxxxxxx.txt
sed -i -e 's/^/AP20|/' ob1xxxxxxxxx.txt
x10=$(grep -f ob1xxxxxxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#secretedcazy=
awk '$31 == "1" && $21 != "." {print $1}' $Annotations| grep -v 'name' > ob1xxxxxxxxxx.txt
sed -i -e 's/^/AP20|/' ob1xxxxxxxxxx.txt
x11=$(grep -f ob1xxxxxxxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
echo $ID $x1 $x2 $x3 $x4 $x5 $x6 $x7 $x8 $x9 $x10 $x11 >> OwnEffectorPredictionReport.txt
done


####################################################################################################################
#MINUs RASPBERRY gene orthologies

cut -f1 /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-minusraspcommon-ranked.tsv > aphanis-rasp-common.txt
cut -f1 /home/theaven/scratch/gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_aphanis-THeavenDRCT72020_1-minusraspaphanis-ranked.tsv > aphanis-rasp-unique.txt
cut -f1 /home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_leucotricha-THeavenp11_1-minusraspcommon-ranked.tsv > leucotricha-rasp-common.txt
cut -f1 /home/theaven/scratch/gene_pred/P_leucotricha/THeavenp11_1/predector_singularity2/results/final_genes_appended_renamed.pep/P_leucotricha-THeavenp11_1-minusraspleucotricha-ranked.tsv > leucotricha-rasp-unique.txt

ls gene_pred/P*/*/predector_singularity3/*_abinitio_annotations3.tsv

grep -f aphanis-rasp-common.txt gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity3/THeavenDRCT72020_1_abinitio_annotations3.tsv > aphanis-rasp-common.tsv
grep -f aphanis-rasp-unique.txt gene_pred/P_aphanis/THeavenDRCT72020_1/predector_singularity3/THeavenDRCT72020_1_abinitio_annotations3.tsv > aphanis-rasp-unique.tsv
grep -f leucotricha-rasp-common.txt gene_pred/P_leucotricha/THeavenp11_1/predector_singularity3/THeavenp11_1_abinitio_annotations3.tsv > leucotricha-rasp-common.tsv
grep -f leucotricha-rasp-unique.txt gene_pred/P_leucotricha/THeavenp11_1/predector_singularity3/THeavenp11_1_abinitio_annotations3.tsv > leucotricha-rasp-unique.tsv

 echo Species+Assembly Predicted_genes pfamID pfam_virulence PHIbase_effector effector_ortholog canonical_SP noncanonical_SP SSCP CSEP Secreted_effectorP3 Secreted_CAZY newSSCP newCSEP newSecreted_effectorP3 > raspPredictionReport.txt
 for file in $(ls *-rasp-*.tsv); do
 ID=$(echo $file | cut -d '.' -f1)
 echo $ID
 Annotations=$file
 Predictedgenes=$(cut -f1 $Annotations| grep -v 'name' | wc -l)
 pfamID=$(awk '$18 != "." {print $1}' $Annotations| grep -v 'name' | wc -l)
 pfamvirulence=$(awk '$20 == "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 PHIbaseeffector=$(awk '$15 == "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 effectorortholog=$(awk '$11 != "." {print $1}' $Annotations| grep -v 'name' | wc -l)
 canonicalsp=$(awk '$31 == "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 noncanonical=$(awk '$31 != "1" && $5 != 'signal' && $2 > 0.6 {print $1}' $Annotations| grep -v 'name' | wc -l)
 sscp=$(awk '$31 == "1" && $36 <= 300 && $39 > 3 {print $1}' $Annotations| grep -v 'name' | wc -l)
 csep=$(awk '$31 == "1" && $7 == "0" {print $1}' $Annotations| grep -v 'name' | wc -l)
 secretedeffector=$(awk '$31 == "1" && $27 == "." && $7 == "0" {print $1}' $Annotations| grep -v 'name' | wc -l)
 secretedcazy=$(awk '$31 == "1" && $21 != "." {print $1}' $Annotations| grep -v 'name' | wc -l)
 newsscp=$(awk '$31 == "1" && $36 <= 300 && $39 > 3 && $11 == "." && $15 != "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 newcsep=$(awk '$31 == "1" && $7 == "0" && $11 == "." && $15 != "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 newsecretedeffector=$(awk '$31 == "1" && $27 == "." && $11 == "." && $15 != "1" && $7 == "0" {print $1}'  $Annotations| grep -v 'name' | wc -l)
 echo $ID $Predictedgenes $pfamID $pfamvirulence $PHIbaseeffector $effectorortholog $canonicalsp $noncanonical $sscp $csep $secretedeffector $secretedcazy $newsscp $newcsep $newsecretedeffector >> raspPredictionReport.txt
 done

for Annotations in $(ls leucotricha-rasp*.tsv); do
ID=$(echo $Annotations | cut -d '.' -f1)
cut -f1 $Annotations| grep -v 'name' > ob1.txt
sed -i -e 's/^/AP20|/' ob1.txt
x1=$(grep -f ob1.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#TFs
#pfamID=
awk '$18 != "." {print $1}' $Annotations| grep -v 'name' > ob1x.txt
sed -i -e 's/^/AP20|/' ob1x.txt
x2=$(grep -f ob1x.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#pfamvirulence=
awk '$20 == "1" {print $1}' $Annotations| grep -v 'name' > ob1xx.txt
sed -i -e 's/^/AP20|/' ob1xx.txt
x3=$(grep -f ob1xx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#PHIbaseeffector=
awk '$15 == "1" {print $1}' $Annotations| grep -v 'name' > ob1xxx.txt
sed -i -e 's/^/AP20|/' ob1xxx.txt
x4=$(grep -f ob1xxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#effectorortholog=
awk '$11 != "." {print $1}' $Annotations| grep -v 'name' > ob1xxxx.txt
sed -i -e 's/^/AP20|/' ob1xxxx.txt
x5=$(grep -f ob1xxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#canonicalsp=
awk '$31 == "1" {print $1}' $Annotations| grep -v 'name' > ob1xxxxx.txt
sed -i -e 's/^/AP20|/' ob1xxxxx.txt
x6=$(grep -f ob1xxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#non-canonical=
awk '$31 != "1" && $5 != 'signal' && $2 > 0.6 {print $1}' $Annotations| grep -v 'name' > ob1xxxxxx.txt
sed -i -e 's/^/AP20|/' ob1xxxxxx.txt
x7=$(grep -f ob1xxxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#sscp=
awk '$31 == "1" && $36 <= 300 && $39 > 3 {print $1}' $Annotations| grep -v 'name' > ob1xxxxxxx.txt
sed -i -e 's/^/AP20|/' ob1xxxxxxx.txt
x8=$(grep -f ob1xxxxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#csep=
awk '$31 == "1" && $7 == "0" {print $1}' $Annotations| grep -v 'name' > ob1xxxxxxxx.txt
sed -i -e 's/^/AP20|/' ob1xxxxxxxx.txt
x9=$(grep -f ob1xxxxxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#secretedeffector=
awk '$31 == "1" && $27 != "." {print $1}' $Annotations| grep -v 'name' > ob1xxxxxxxxx.txt
sed -i -e 's/^/AP20|/' ob1xxxxxxxxx.txt
x10=$(grep -f ob1xxxxxxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#secretedcazy=
awk '$31 == "1" && $21 != "." {print $1}' $Annotations| grep -v 'name' > ob1xxxxxxxxxx.txt
sed -i -e 's/^/AP20|/' ob1xxxxxxxxxx.txt
x11=$(grep -f ob1xxxxxxxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
echo $ID $x1 $x2 $x3 $x4 $x5 $x6 $x7 $x8 $x9 $x10 $x11 >> raspPredictionReport.txt
done

for Annotations in $(ls aphanis-rasp*.tsv); do
ID=$(echo $Annotations | cut -d '.' -f1)
cut -f1 $Annotations| grep -v 'name' > ob1.txt
sed -i -e 's/^/ST20|/' ob1.txt
x1=$(grep -f ob1.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#TFs
#pfamID=
awk '$18 != "." {print $1}' $Annotations| grep -v 'name' > ob1x.txt
sed -i -e 's/^/ST20|/' ob1x.txt
x2=$(grep -f ob1x.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#pfamvirulence=
awk '$20 == "1" {print $1}' $Annotations| grep -v 'name' > ob1xx.txt
sed -i -e 's/^/ST20|/' ob1xx.txt
x3=$(grep -f ob1xx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#PHIbaseeffector=
awk '$15 == "1" {print $1}' $Annotations| grep -v 'name' > ob1xxx.txt
sed -i -e 's/^/ST20|/' ob1xxx.txt
x4=$(grep -f ob1xxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#effectorortholog=
awk '$11 != "." {print $1}' $Annotations| grep -v 'name' > ob1xxxx.txt
sed -i -e 's/^/ST20|/' ob1xxxx.txt
x5=$(grep -f ob1xxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#canonicalsp=
awk '$31 == "1" {print $1}' $Annotations| grep -v 'name' > ob1xxxxx.txt
sed -i -e 's/^/ST20|/' ob1xxxxx.txt
x6=$(grep -f ob1xxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#non-canonical=
awk '$31 != "1" && $5 != 'signal' && $2 > 0.6 {print $1}' $Annotations| grep -v 'name' > ob1xxxxxx.txt
sed -i -e 's/^/ST20|/' ob1xxxxxx.txt
x7=$(grep -f ob1xxxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#sscp=
awk '$31 == "1" && $36 <= 300 && $39 > 3 {print $1}' $Annotations| grep -v 'name' > ob1xxxxxxx.txt
sed -i -e 's/^/ST20|/' ob1xxxxxxx.txt
x8=$(grep -f ob1xxxxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#csep=
awk '$31 == "1" && $7 == "0" {print $1}' $Annotations| grep -v 'name' > ob1xxxxxxxx.txt
sed -i -e 's/^/ST20|/' ob1xxxxxxxx.txt
x9=$(grep -f ob1xxxxxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#secretedeffector=
awk '$31 == "1" && $27 != "." {print $1}' $Annotations| grep -v 'name' > ob1xxxxxxxxx.txt
sed -i -e 's/^/ST20|/' ob1xxxxxxxxx.txt
x10=$(grep -f ob1xxxxxxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
#secretedcazy=
awk '$31 == "1" && $21 != "." {print $1}' $Annotations| grep -v 'name' > ob1xxxxxxxxxx.txt
sed -i -e 's/^/ST20|/' ob1xxxxxxxxxx.txt
x11=$(grep -f ob1xxxxxxxxxx.txt analysis/orthology/orthofinder/All6_isolates/formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt | uniq | wc -l)
echo $ID $x1 $x2 $x3 $x4 $x5 $x6 $x7 $x8 $x9 $x10 $x11 >> raspPredictionReport.txt
done
```










































## Comparing to wider mildew clade:
#### FCS
```bash
for Genome in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/*a); do
Jobs=$(squeue -u theaven| grep 'fcs' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 120s
printf "."
Jobs=$(squeue -u theaven| grep 'fcs' | wc -l)
done
Species=$(echo $Genome | cut -d '/' -f7,8)
if [ "$Species" = "Amorphotheca/resinae" ]; then
TAXID=5101
elif [ "$Species" = "Arachnopeziza/araneosa" ]; then
  TAXID=2259763
elif [ "$Species" = "Ascocoryne/sarcoides" ]; then
  TAXID=139061
elif [ "$Species" = "Aspergillus/nidulans" ]; then
  TAXID=162425
elif [ "$Species" = "Blumeria/graminis" ]; then
  TAXID=34373
elif [ "$Species" = "Blumeria/graminis-secale" ]; then
  TAXID=34373
elif [ "$Species" = "Blumeria/hordei" ]; then
  TAXID=2867405
elif [ "$Species" = "Botrytis/cineria" ]; then
  TAXID=40559
elif [ "$Species" = "Chlorociboria/aeruginascens" ]; then
  TAXID=296797
elif [ "$Species" = "Colletotrichum/higginsianum" ]; then
  TAXID=80884
elif [ "$Species" = "Coprinus/cinereus" ]; then
  TAXID=5346
elif [ "$Species" = "Cryptococcus/neoformans" ]; then
  TAXID=5207
elif [ "$Species" = "Drepanopeziza/brunnea" ]; then
  TAXID=698440
elif [ "$Species" = "Erysiphe/alphitoides" ]; then
  TAXID=157594
elif [ "$Species" = "Erysiphe/necator" ]; then
  TAXID=52586
elif [ "$Species" = "Erysiphe/neolycopersici" ]; then
  TAXID=212602
elif [ "$Species" = "Erysiphe/pisi" ]; then
  TAXID=36044
elif [ "$Species" = "Erysiphe/pulchra" ]; then
  TAXID=225359
elif [ "$Species" = "Fusarium/graminearum" ]; then
  TAXID=5518
elif [ "$Species" = "Fusarium/oxysporum" ]; then
  TAXID=5507
elif [ "$Species" = "Glarea/lozoyensis" ]; then
  TAXID=101852
elif [ "$Species" = "Golovinomyces/cichoracearum" ]; then
  TAXID=62708
elif [ "$Species" = "Golovinomyces/magnicellulatus" ]; then
  TAXID=62714
elif [ "$Species" = "Golovinomyces/orontii" ]; then
  TAXID=62715
elif [ "$Species" = "Leveillula/taurrica" ]; then
  TAXID=62718
elif [ "$Species" = "Magnaporthe/oryzae" ]; then
  TAXID=318829
elif [ "$Species" = "Melampsora/laricis-populina" ]; then
  TAXID=203908
elif [ "$Species" = "Mollisia/scopiformis" ]; then
  TAXID=149040
elif [ "$Species" = "Neobulgaria/alba" ]; then
  TAXID=1461596
elif [ "$Species" = "Neurospora/crassa" ]; then
  TAXID=5141
elif [ "$Species" = "Oidiodendron/maius" ]; then
  TAXID=78148
elif [ "$Species" = "Oidium/heveae" ]; then
  TAXID=299130
elif [ "$Species" = "Parauncinula/polyspora" ]; then
  TAXID=2016243
elif [ "$Species" = "Phialocephala/subalpina" ]; then
  TAXID=576137
elif [ "$Species" = "Phyllactinia/moricola" ]; then
  TAXID=57460
elif [ "$Species" = "Pleochaeta/shiraiana" ]; then
  TAXID=57462
elif [ "$Species" = "Pleurotus/ostreatus" ]; then
  TAXID=5322
elif [ "$Species" = "Podosphaera/aphanis" ]; then
  TAXID=79252
elif [ "$Species" = "Podosphaera/cerasii" ]; then
  TAXID=2086344
elif [ "$Species" = "Podosphaera/fusca" ]; then
  TAXID=62727
elif [ "$Species" = "Podosphaera/leucotricha" ]; then
  TAXID=79249
elif [ "$Species" = "Podosphaera/xanthii" ]; then
  TAXID=135283
elif [ "$Species" = "Psilocybe/cubensis" ]; then
  TAXID=181762
elif [ "$Species" = "Puccinia/graminis" ]; then
  TAXID=5297
elif [ "$Species" = "Puccinia/striiformis" ]; then
  TAXID=27350
elif [ "$Species" = "Puccinia/triticina" ]; then
  TAXID=208348
elif [ "$Species" = "Pyricularia/oryzae" ]; then
  TAXID=318829
elif [ "$Species" = "Saccharomyces/cerevisiae" ]; then
  TAXID=4932
elif [ "$Species" = "Schizophyllum/commune" ]; then
  TAXID=5334
elif [ "$Species" = "Sclerotinia/sclerotiorum" ]; then
  TAXID=5180
elif [ "$Species" = "Tuber/melanosporum" ]; then
  TAXID=39416
elif [ "$Species" = "Ustilago/maydis" ]; then
  TAXID=5270
else
  TAXID=NA
fi
OutDir=$(dirname $Genome)/fcs
OutFile=$(basename $Genome | sed 's@.c.fna@@g'| sed 's@.fna@@g'| sed 's@.fasta@@g')
mkdir $OutDir
ProgDir=/home/theaven/scratch/apps
sbatch $ProgDir/fcs/run_fcs.sh $Genome $TAXID $OutDir $OutFile
#
done
```

#### EarlGreyTE
```bash
head -n 2 /home/theaven/scratch/uncompressed/genomes/HEAVEN_apple2020.fna > temp_earlgreytest.fna

conda activate earlgrey
for Genome in $(ls /home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_900237765.1/fcs/Blumeriahordei_GCA_900237765.1_BghRACE1_v1_genomic_clean.fasta); do
  OutFile=$(echo $Genome | cut -d '/' -f7 | cut -c 1-3)_$(echo $Genome | cut -d '/' -f8 | cut -c 1-3)_$(echo $Genome | cut -d '/' -f9)
  OutDir=$(dirname $Genome)/earlgreyte/${OutFile}
  ProgDir=~/scratch/apps
  RMsearch=erysiphales
  mkdir -p $OutDir
  sbatch ${ProgDir}/earlgrey/run_earlgrey.sh $Genome $OutFile $OutDir $RMsearch
done
conda deactivate
```
### Predict genes - no RNA evidence
#### Repeatmasker
```bash
for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/temp_genomes/Parauncinula/*/*/*a); do
OutFile=$(basename $Genome | sed 's@.fa@@g')
OutDir=$(dirname $Genome)/repeatmodeler
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_repeatmodeler.sh $Genome $OutFile $OutDir
done
#58170701, 58940849

for Genome in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/temp_genomes/Parauncinula/*/*/*a); do
Species=$(echo $Genome | cut -d '/' -f9,10 | sed 's@/@_@g')
Repeat_library=$(ls $(dirname $Genome)/repeatmodeler/*/consensi.fa.classified)
OutFile=$(basename $Genome | sed 's@.fa@@g')
OutDir=$(dirname $Genome)/repeatmasker
ProgDir=~/git_repos/Wrappers/NBI
mkdir $OutDir
sbatch $ProgDir/run_repeatmasker4.1.5.sh $Genome $OutFile $OutDir $Species $Repeat_library
done

#move to gruffalo
```
#### TransposonPSI
```bash
conda activate transposonpsi
for assembly in $(ls temp_genomes/*/*/*/*na); do
    ID=$(echo $assembly | cut -d '/' -f2 | cut -c 1-3)_$(echo $assembly | cut -d '/' -f3 | cut -c 1-3)_$(echo $assembly | cut -d '/' -f4)
    echo $ID
    ProgDir=~/scratch/apps/repeat_masking
    #ProgDir=~/git_repos/tools/seq_tools/repeat_masking
    OutDir=$(dirname $assembly)/transposonpsi
    mkdir -p $OutDir
    if [[ ! -f ${OutDir}/formatdb.log ]]; then
    echo missing for $ID
    sbatch $ProgDir/gomez_transposonPSI.sh $assembly $OutDir $ID
    else 
    echo done for $ID
    fi
done
conda deactivate

conda activate bedtools
#Merge the masking:
for Assembly in $(ls temp_genomes/*/*/*/*.f*); do
    ID=$(dirname $Assembly | cut -d '/' -f2,3,4)
File=$(ls temp_download/$ID/*.fna.masked)
OutDir=$(dirname $Assembly)/repeatmasking/combined
TPSI=$(ls temp_genomes/${ID}/transposonpsi/*allHits.chains.gff3)
OutFile=${OutDir}/$(basename $File | sed 's/.fna.masked/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
if [[ ! -f ${OutFile} ]]; then
if [[ -f ${TPSI} && -f ${File} ]]; then
mkdir -p $OutDir
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo ${OutFile}
No=$(cat $OutFile | grep -v '>' | tr -d '\n' | wc -c)
No2=$(cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' ')
echo $ID >> mildews/combinedrepeatreport.txt
echo "Number of bases:" >> mildews/combinedrepeatreport.txt
echo $No >> mildews/combinedrepeatreport.txt
echo "Number of masked bases:" >> mildews/combinedrepeatreport.txt
echo $No2 >> mildews/combinedrepeatreport.txt
echo "% of masked bases:" >> mildews/combinedrepeatreport.txt
echo "100/$No*$No2" | bc -l >> mildews/combinedrepeatreport.txt
printf "\n" >> mildews/combinedrepeatreport.txt
else 
echo "inputs missing for $ID"
fi
else
echo "already run for $ID"
fi 
done
conda deactivate
```
#### BRAKER
```bash
conda activate braker
for assembly in $(ls temp_genomes/*/*/*/repeatmasking/combined/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
  ID=$(echo $assembly | cut -d '/' -f2 | cut -c 1-3)_$(echo $assembly | cut -d '/' -f3 | cut -c 1-3)_$(echo $assembly | cut -d '/' -f4)
  OutDir=mildews/${ID}/gene_pred/prothint
  ProgDir=/home/theaven/scratch/apps/braker
  if [[ ! -f ${OutDir}/prothint.gff ]]; then
  if [[ -f ${assembly} ]]; then
    echo "running for $ID"
  sbatch $ProgDir/prothint_fungi.sh $assembly $OutDir
  else 
  echo "inputs missing for $ID"
  fi
  else
  echo "already run for $ID"
  fi 
done
conda deactivate

conda activate braker
for Assembly in $(ls /home/theaven/scratch/uncompressed/genomes/Sclerotinia/sclerotiorum/GCF_000146945.1/fcs/repeatmasking/combined/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
  ID=$(echo $Assembly | cut -d '/' -f7 | cut -c 1-3)_$(echo $Assembly | cut -d '/' -f8 | cut -c 1-3)_$(echo $Assembly | cut -d '/' -f9)
  echo $ID
  OutDir=/home/theaven/scratch/uncompressed/genomes/Sclerotinia/sclerotiorum/GCF_000146945.1/fcs/gene_pred/braker2
  GeneModelName=${ID}_genemodels4
  Hintfile=/home/theaven/scratch/uncompressed/genomes/Sclerotinia/sclerotiorum/GCF_000146945.1/fcs/gene_pred/prothint/prothint_augustus.gff
  ProgDir=/home/theaven/scratch/apps/braker
  if [[ ! -f ${OutDir}/braker.gff3 ]]; then
  if [[ -f ${Hintfile} ]]; then
    echo "running for $ID"
    mkdir -p $OutDir
  sbatch $ProgDir/braker_fungi_noRNA.sh $Assembly $OutDir $GeneModelName $Hintfile
  else 
  echo "inputs missing for $ID"
  fi
  else
  echo "already run for $ID"
  fi 
done 
conda deactivate
```
```bash
conda activate predector2.7
for Assembly in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/repeatmasking/combined/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
  #remove duplicates and rename:
  ID=$(basename $Assembly|cut -f1,2,3 -d '_'| sed 's@_contigs_softmasked_repeatmasker_TPSI_appended@@g')
  echo $ID
  Gff=$(echo $Assembly | cut -d '/' -f1,2,3,4,5,6,7,8,9,10 )/gene_pred/braker/braker.gff3
  ab_gff=$(echo $Assembly | cut -d '/' -f1,2,3,4,5,6,7,8,9,10 )/gene_pred/braker/augustus.ab_initio.gff3
  FinalDir=$(echo $Assembly | cut -d '/' -f1,2,3,4,5,6,7,8,9,10 )/gene_pred/braker/final-02042924
  echo $FinalDir
  if [ ! -e "${FinalDir}/final_genes_renamed.pep.fasta" ]; then
  mkdir $FinalDir
  GffFiltered=$FinalDir/filtered_duplicates.gff
  ab_GffFiltered=$FinalDir/ab_initio_filtered_duplicates.gff
  ProgDir=/home/theaven/scratch/apps/tools
  #ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Gene_prediction 
  $ProgDir/remove_dup_features.py --inp_gff $Gff --out_gff $GffFiltered
  $ProgDir/remove_dup_features.py --inp_gff $ab_gff --out_gff $ab_GffFiltered
  awk '{ gsub("\\.","_",$1); print}' $GffFiltered > ${GffFiltered}.2
  awk '{ gsub("\\.","_",$1); print}' $ab_GffFiltered > ${ab_GffFiltered}.2
  GffRenamed=$FinalDir/final_genes_renamed.gff3 
  LogFile=$FinalDir/final_genes_renamed.log 
  ab_GffRenamed=$FinalDir/ab_initio_final_genes_renamed.gff3 
  ab_LogFile=$FinalDir/ab_initio_final_genes_renamed.log
  $ProgDir/gff_rename_genes.py --inp_gff ${GffFiltered}.2 --conversion_log $LogFile > $GffRenamed 
  $ProgDir/gff_rename_genes.py --inp_gff ${ab_GffFiltered}.2 --conversion_log $ab_LogFile > $ab_GffRenamed  
  awk '{ gsub("\\.","_"); print}' $Assembly > ${Assembly}.2
  $ProgDir/gff2fasta.pl ${Assembly}.2 $GffRenamed $FinalDir/final_genes_renamed
  $ProgDir/gff2fasta.pl ${Assembly}.2 $ab_GffRenamed $FinalDir/ab_initio_final_genes_renamed
  sed -i 's/*/X/g' $FinalDir/final_genes_renamed.pep.fasta
  sed -i 's/*/X/g' $FinalDir/ab_initio_final_genes_renamed.pep.fasta
  #remove any discovered duplicates
  Output=$(dirname $LogFile)/$(echo $LogFile |cut -f2 -d '/')_duplicated_genes
  ab_Output=$(dirname $ab_LogFile)/$(echo $ab_LogFile |cut -f2 -d '/')ab_initio_duplicated_genes
  for prediction in $(cat $LogFile |cut -f1|sort|uniq -c|sort -nr | awk '$1>1' | awk '{print $2}'); do
    printf ">\n" >> $Output
    grep -w $prediction $LogFile | awk '{print $3}' >> $Output
  done
  sed -e '/>/,+1d' ${Output} > ${Output}.txt #removed one of each duplicate from the removal list
  awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < $FinalDir/final_genes_renamed.cdna.fasta > linearized.fasta
  for var1 in $(cat ${Output}.txt); do
    sed -i "/>${var1}.t1/d" linearized.fasta
  done
  tr "\t" "\n" < linearized.fasta | fold -w 60 > $FinalDir/final_genes_renamed.cdna.fasta
  awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < $FinalDir/final_genes_renamed.pep.fasta > linearized.fasta
  for var1 in $(cat ${Output}.txt); do
    sed -i "/>${var1}.t1/d" linearized.fasta
  done
  tr "\t" "\n" < linearized.fasta | fold -w 60 > $FinalDir/final_genes_renamed.pep.fasta
  for prediction in $(cat $ab_LogFile |cut -f1|sort|uniq -c|sort -nr | awk '$1>1' | awk '{print $2}'); do
    printf ">\n" >> $ab_Output
    grep -w $prediction $ab_LogFile | awk '{print $3}' >> $ab_Output
  done
  sed -e '/>/,+1d' ${ab_Output} > ${ab_Output}.txt #removed one of each duplicate from the removal list
  awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < $FinalDir/ab_initio_final_genes_renamed.cdna.fasta > linearized.fasta
  for var1 in $(cat ${ab_Output}.txt); do
    sed -i "/>${var1}.t1/d" linearized.fasta
  done
  tr "\t" "\n" < linearized.fasta | fold -w 60 > $FinalDir/ab_initio_final_genes_renamed.cdna.fasta
  awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < $FinalDir/ab_initio_final_genes_renamed.pep.fasta > linearized.fasta
  for var1 in $(cat ${ab_Output}.txt); do
    sed -i "/>${var1}.t1/d" linearized.fasta
  done
  tr "\t" "\n" < linearized.fasta | fold -w 60 > $FinalDir/ab_initio_final_genes_renamed.pep.fasta
  rm $Output
  rm $ab_Output
  #count predicated genes
  cat $FinalDir/final_genes_renamed.cdna.fasta | grep '>' > $FinalDir/genenames.txt
  echo ""  >> other/brakerreport.txt
  echo $ID >> other/brakerreport.txt
  echo No. of prothint guided genes: >> other/brakerreport.txt
  grep -c -i '>' $FinalDir/final_genes_renamed.cdna.fasta >> other/brakerreport.txt
  cat $FinalDir/ab_initio_final_genes_renamed.cdna.fasta | grep '>' > $FinalDir/ab_initio_genenames.txt
  echo No. of ab initio predicted genes: >> other/brakerreport.txt
  grep -c -i '>' $FinalDir/ab_initio_final_genes_renamed.cdna.fasta >> other/brakerreport.txt
  rm ${Assembly}.2
  rm ${GffFiltered}.2
  rm ${ab_GffFiltered}.2
  fi
done
```
Ensure that protein sequences are only from genes in 'clean' contigs
```bash
for file in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/*_clean.fasta); do
grep '>' $file | cut -d ' ' -f1 | sed 's@>@@g'| sed 's@\.@_@g' > temp_contig_ids.txt
Dir=$(dirname $file)/gene_pred/braker
if [ -d "${Dir}/final-02042924" ]; then
gff=${Dir}/final-02042924/final_genes_renamed.gff3
fasta=${Dir}/final-02042924/final_genes_renamed.pep.fasta
out=${Dir}/final-02042924/clean_final_genes_renamed.pep.fasta
grep -f temp_contig_ids.txt $gff | grep 'mRNA'| cut -d ' ' -f9  | cut -d ';' -f1 | cut -d '=' -f2 > temp_gene_ids.txt
python3 /home/theaven/scratch/apps/tools/seq_get.py --id_file temp_gene_ids.txt --input $fasta --output $out
ls -lh $out
else
echo error
fi
#rm temp_contig_ids.txt
#rm temp_gene_ids.txt
done

for file in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/*_clean.fasta); do
grep '>' $file | cut -d ' ' -f1 | sed 's@>@@g'| sed 's@\.@_@g' > temp_contig_ids.txt
Dir=$(dirname $file)/gene_pred/braker
if [ -d "${Dir}/final-02042924" ]; then
gff=${Dir}/final-02042924/ab_initio_final_genes_renamed.gff3
fasta=${Dir}/final-02042924/ab_initio_final_genes_renamed.pep.fasta
out=${Dir}/final-02042924/clean_ab_initio_final_genes_renamed.pep.fasta
grep -f temp_contig_ids.txt $gff | grep 'mRNA'| cut -d ' ' -f9  | cut -d ';' -f1 | cut -d '=' -f2 > temp_gene_ids.txt
python3 /home/theaven/scratch/apps/tools/seq_get.py --id_file temp_gene_ids.txt --input $fasta --output $out
ls -lh $out
else
echo error
fi
#rm temp_contig_ids.txt
#rm temp_gene_ids.txt
done

ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta
ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/clean_ab_initio_final_genes_renamed.pep.fasta
```
#### Interproscan
```bash
#tar and move between clusters
list="/home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000151065.3/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153116/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153117/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153118/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153119/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153120/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Parauncinula/polyspora/Parp01/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/fusca/GCA_030378345.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/xanthii/GCA_028751805.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000417025.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000417865.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000418435.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000441875.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_900519115.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_905067625.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_000401675.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_900237765.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_900239735.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_900638725.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/alphitoides/CLCBIO/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798715.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798735.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798755.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798775.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798795.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_016906895.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/neolycopersici/GCA_003610855.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/pisi/GCA_000208805.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/pisi/GCA_000214055.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/pulchra/GCA_002918395.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/cichoracearum/GCA_003611195.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/cichoracearum/GCA_003611215.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/cichoracearum/GCA_003611235.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/magnicellulatus/GCA_006912115.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/orontii/MGH1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Leveillula/taurrica/CADEPA01/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Phyllactinia/moricola/GCA_019455665.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Pleochaeta/shiraiana/GCA_019455505.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/DRT72020/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/DRT72021/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/SCOTT2020/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/cerasii/GCA_018398735.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/GCA_013170925.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGB2019/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGB2021/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGBp112020/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/xanthii/GCA_010015925.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/xanthii/GCA_014884795.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Oidium/heveae/GCA_003957845.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_024703715.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta"
for File in $list; do 
InFile=$(echo $File | sed 's@/home@/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/download2/home@g')
ls $InFile
OutDir=$(dirname $InFile)/interproscan/Prothint
ProgDir=~/git_repos/Wrappers/NBI
mkdir -p $OutDir
sbatch $ProgDir/run_interproscan.sh $InFile $OutDir
done

ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/interproscan/Prothint/
```
#### Predector
```bash
screen -S predector
conda activate predector
for proteome in $(ls /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000151065.3/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153116/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153117/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153118/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153119/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153120/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Parauncinula/polyspora/Parp01/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/fusca/GCA_030378345.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/xanthii/GCA_028751805.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000417025.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000417865.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000418435.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000441875.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_900519115.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_905067625.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_000401675.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_900237765.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_900239735.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_900638725.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/alphitoides/CLCBIO/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798715.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798735.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798755.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798775.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798795.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_016906895.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/neolycopersici/GCA_003610855.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/pisi/GCA_000208805.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/pisi/GCA_000214055.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/pulchra/GCA_002918395.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/cichoracearum/GCA_003611195.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/cichoracearum/GCA_003611215.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/cichoracearum/GCA_003611235.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/magnicellulatus/GCA_006912115.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/orontii/MGH1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Leveillula/taurrica/CADEPA01/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Phyllactinia/moricola/GCA_019455665.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Pleochaeta/shiraiana/GCA_019455505.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/DRT72020/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/DRT72021/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/SCOTT2020/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/cerasii/GCA_018398735.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/GCA_013170925.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGB2019/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGB2021/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGBp112020/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/xanthii/GCA_010015925.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/xanthii/GCA_014884795.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Oidium/heveae/GCA_003957845.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_024703715.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta); do
#  ls $proteome  
    OutDir=$(dirname $proteome)/predector/Prothint
    ProgDir=/home/theaven/scratch/apps/predector
#    ls $OutDir/results/clean_final_genes_renamed.pep/clean_final_genes_renamed.pep-ranked.tsv
if [ ! -f "$OutDir/results/clean_final_genes_renamed.pep/clean_final_genes_renamed.pep-ranked.tsv" ]; then
echo $proteome >> predectorpathlog.txt
sbatch $ProgDir/predector_singularity.sh $proteome 1.2.6 $OutDir  2>&1 >> predectorpathlog.txt
echo $proteome
else 
echo $OutDir
fi
#sleep 30
done 

for proteome2 in $(ls /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000151065.3/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153116/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153117/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153118/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153119/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153120/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Parauncinula/polyspora/Parp01/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/fusca/GCA_030378345.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/xanthii/GCA_028751805.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000417025.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000417865.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000418435.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000441875.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_900519115.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_905067625.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_000401675.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_900237765.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_900239735.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_900638725.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/alphitoides/CLCBIO/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798715.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798735.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798755.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798775.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798795.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_016906895.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/neolycopersici/GCA_003610855.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/pisi/GCA_000208805.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/pisi/GCA_000214055.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/pulchra/GCA_002918395.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/cichoracearum/GCA_003611195.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/cichoracearum/GCA_003611215.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/cichoracearum/GCA_003611235.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/magnicellulatus/GCA_006912115.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/orontii/MGH1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Leveillula/taurrica/CADEPA01/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Phyllactinia/moricola/GCA_019455665.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Pleochaeta/shiraiana/GCA_019455505.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/DRT72020/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/DRT72021/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/SCOTT2020/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/cerasii/GCA_018398735.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/GCA_013170925.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGB2019/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGB2021/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGBp112020/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/xanthii/GCA_010015925.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/xanthii/GCA_014884795.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Oidium/heveae/GCA_003957845.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_024703715.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta); do
  proteome=$(echo $proteome2 | sed 's@clean_final_genes_renamed.pep.fasta@clean_ab_initio_final_genes_renamed.pep.fasta@g') 
    OutDir=$(dirname $proteome)/predector/abinitio
    ProgDir=/home/theaven/scratch/apps/predector
    ls $OutDir/results/clean_ab_initio_final_genes_renamed.pep/clean_ab_initio_final_genes_renamed.pep-ranked.tsv
if [ ! -f "$OutDir/results/clean_ab_initio_final_genes_renamed.pep/clean_ab_initio_final_genes_renamed.pep-ranked.tsv" ]; then
echo $proteome >> predectorpathlog.txt
sbatch $ProgDir/predector_singularity.sh $proteome 1.2.6 $OutDir  2>&1 >> predectorpathlog.txt
echo $proteome
else 
echo $OutDir
fi
sleep 30
done 
```
SecretomeP 2.0 for non-canonically secreted proteins:
```bash
#Secretome1:
screen -S secretomep
conda activate jellyfish
#linearize and edit 
for proteome in $(ls genomes/*/*/*/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta ); do
  Outfile=$(basename $proteome | cut -d '.' -f1)
  OutDir=$(dirname $proteome)
  ProgDir=/home/theaven/scratch/apps/secretomep
  echo ${OutDir}/${Outfile}_genes_40.fsa
  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $proteome | sed "s@file_1_@@g" | sed -e "s@-@N@g" | sed -e "s@X@@g" | awk '/^>/{val=$0;next} length($0)>40{print val ORS $0}' > ${OutDir}/${Outfile}_genes_40.fsa
  sbatch $ProgDir/secretomep.sh ${OutDir}/${Outfile}_genes_40.fsa ${OutDir} ${Outfile}_genes_40_secretomep.out 2>&1 >> secretomeplog.txt
done
#4152856

for ome in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*_genes_40_secretomep.out); do 
sed -i -n '/SecretomeP 1.0g predictions/,$p' $ome 
sed -i '$d' $ome 
sed -i '1d' $ome 
sed -i '1d' $ome 
sed -i '1d' $ome 
sed -i '1d' $ome 
done

for x in $(ls */*/gene_pred/braker/final/final_genes_renamed.pep.fasta); do
ID=$(echo $x | cut -d '/' -f2 | cut -d '.' -f1)
cp $x downloadprothint/$ID.pep.fasta
done
for x in $(ls */*/gene_pred/braker/final/ab_initio_final_genes_renamed.pep.fasta); do
ID=$(echo $x | cut -d '/' -f2 | cut -d '.' -f1)
cp $x downloadabinitio/$ID.pep.fasta
done
```
```bash
#Secretome1:
screen -S secretomep2
conda activate jellyfish
#linearize and edit 
for proteome in $(ls genomes/*/*/*/fcs/gene_pred/braker/final-02042924/clean_ab_initio_final_genes_renamed.pep.fasta ); do
  Outfile=$(basename $proteome | cut -d '.' -f1)
  OutDir=$(dirname $proteome)
  ProgDir=/home/theaven/scratch/apps/secretomep
  echo ${OutDir}/${Outfile}_genes_40.fsa
  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $proteome | sed "s@file_1_@@g" | sed -e "s@-@N@g" | sed -e "s@X@@g" | awk '/^>/{val=$0;next} length($0)>40{print val ORS $0}' > ${OutDir}/${Outfile}_genes_40.fsa
  sbatch $ProgDir/secretomep.sh ${OutDir}/${Outfile}_genes_40.fsa ${OutDir} ${Outfile}_genes_40_secretomep.out 2>&1 >> secretomeplog.txt
done
#4152856

for ome in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*_genes_40_secretomep.out); do 
sed -i -n '/SecretomeP 1.0g predictions/,$p' $ome 
sed -i '$d' $ome 
sed -i '1d' $ome 
sed -i '1d' $ome 
sed -i '1d' $ome 
sed -i '1d' $ome 
done

for x in $(ls */*/gene_pred/braker/final/final_genes_renamed.pep.fasta); do
ID=$(echo $x | cut -d '/' -f2 | cut -d '.' -f1)
cp $x downloadprothint/$ID.pep.fasta
done
for x in $(ls */*/gene_pred/braker/final/ab_initio_final_genes_renamed.pep.fasta); do
ID=$(echo $x | cut -d '/' -f2 | cut -d '.' -f1)
cp $x downloadabinitio/$ID.pep.fasta
done
```
Blast for homology to non-mildew species
```bash
screen -S blast
conda activate blast+
for proteome in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta ); do 
    Jobs=$(squeue -u theaven |grep 'blast'| wc -l)
    echo x
    while [ $Jobs -gt 10 ]; do
      sleep 120s
      printf "."
      Jobs=$(squeue -u theaven |grep 'blast'| wc -l)
    done
OutDir=$(dirname $proteome)
Database=/home/theaven/scratch/uncompressed/csep_db/db/CSEPdb
ProgDir=/home/theaven/scratch/apps/blast
sbatch $ProgDir/blastp.sh $proteome $Database $OutDir/outblast.tsv 2>&1 >> blastoutlog.txt
done
```
```bash
screen -S blast2
conda activate blast+
for proteome in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/clean_ab_initio_final_genes_renamed.pep.fasta ); do 
    Jobs=$(squeue -u theaven |grep 'blast'| wc -l)
    echo x
    while [ $Jobs -gt 10 ]; do
      sleep 120s
      printf "."
      Jobs=$(squeue -u theaven |grep 'blast'| wc -l)
    done
OutDir=$(dirname $proteome)
Database=/home/theaven/scratch/uncompressed/csep_db/db/CSEPdb
ProgDir=/home/theaven/scratch/apps/blast
sbatch $ProgDir/blastp.sh $proteome $Database $OutDir/ab_initio_outblast.tsv 2>&1 >> blastoutlog.txt
done
```
ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/predector/Prothint/results/clean_final_genes_renamed.pep/clean_final_genes_renamed.pep-ranked.tsv
ls genomes/*/*/*/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta
ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*_genes_40_secretomep.out | grep -v 'ab_initio'

ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/predector/abinitio/results/clean_ab_initio_final_genes_renamed.pep/clean_ab_initio_final_genes_renamed.pep-ranked.tsv
ls genomes/*/*/*/fcs/gene_pred/braker/final-02042924/clean_ab_initio_final_genes_renamed.pep.fasta
ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/clean_ab_initio*_genes_40_secretomep.out 


Joining data into predector files - remember not to include mildew ortholog 'ingenes'
```bash
for Annotations in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/predector/Prothint/results/clean_final_genes_renamed.pep/clean_final_genes_renamed.pep-ranked.tsv); do
ID=$(echo $Annotations | cut -d '/' -f7 | cut -c 1-3)_$(echo $Annotations | cut -d '/' -f8 | cut -c 1-3)_$(echo $Annotations | cut -d '/' -f9)
blast=$(echo $Annotations | cut -d '/' -f1,2,3,4,5,6,7,8,9,10,11,12,13)
Secretomep=$(echo $Annotations | cut -d '/' -f1,2,3,4,5,6,7,8,9,10,11,12,13)/clean_final_genes_renamed_genes_40_secretomep
OutDir=$(echo $Annotations | cut -d '/' -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
cat ${blast}/outblast.tsv | cut -f1 |sort | uniq > ${blast}/outgenes.tsv
for f in ${blast}/outgenes.tsv; do sed -i "s/$/\t1/" $f; done
printf 'name\toutgroup_orthology_match' >> ${blast}/outgenes.tsv
 join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84 <(sort $Annotations) <(sort ${blast}/outgenes.tsv) > ${OutDir}/${ID}_annotations.tsv
printf 'name\tNN-score\tOdds\tWeighted\tWarning' >> ${Secretomep}.out
awk -v OFS='\t' '{$1=$1}1' ${Secretomep}.out > ${Secretomep}.tsv
join -t $'\t' -a 1 -e 0 -o 0,2.2,2.3,2.4,2.5,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84,1.85 <(sort ${OutDir}/${ID}_annotations.tsv) <(sort ${Secretomep}.tsv) > ${OutDir}/${ID}_annotations3.tsv
done

for Annotations in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/predector/abinitio/results/clean_ab_initio_final_genes_renamed.pep/clean_ab_initio_final_genes_renamed.pep-ranked.tsv); do
ID=$(echo $Annotations | cut -d '/' -f7 | cut -c 1-3)_$(echo $Annotations | cut -d '/' -f8 | cut -c 1-3)_$(echo $Annotations | cut -d '/' -f9)
blast=$(echo $Annotations | cut -d '/' -f1,2,3,4,5,6,7,8,9,10,11,12,13)
Secretomep=$(echo $Annotations | cut -d '/' -f1,2,3,4,5,6,7,8,9,10,11,12,13)/clean_ab_initio_final_genes_renamed_genes_40_secretomep
OutDir=$(echo $Annotations | cut -d '/' -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
cat ${blast}/ab_initio_outblast.tsv | cut -f1 |sort | uniq > ${blast}/ab_initio_outgenes.tsv
for f in ${blast}/ab_initio_outgenes.tsv; do sed -i "s/$/\t1/" $f; done
printf 'name\toutgroup_orthology_match' >> ${blast}/ab_initio_outgenes.tsv
 join -t $'\t' -a 1 -e 0 -o 0,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84 <(sort $Annotations) <(sort ${blast}/ab_initio_outgenes.tsv) > ${OutDir}/${ID}_ab_initio_annotations.tsv
printf 'name\tNN-score\tOdds\tWeighted\tWarning' >> ${Secretomep}.out
awk -v OFS='\t' '{$1=$1}1' ${Secretomep}.out > ${Secretomep}.tsv
join -t $'\t' -a 1 -e 0 -o 0,2.2,2.3,2.4,2.5,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84,1.85 <(sort ${OutDir}/${ID}_ab_initio_annotations.tsv) <(sort ${Secretomep}.tsv) > ${OutDir}/${ID}_ab_initio_annotations3.tsv
done

#Making report files
 echo Species+Assembly Predicted_genes pfamID pfam_virulence PHIbase_effector effector_ortholog canonical_SP noncanonical_SP SSCP CSEP Secreted_effectorP3 Secreted_CAZY newSSCP newCSEP newSecreted_effectorP3 all csep_val phiboref > EffectorPredictionReport.txt
for Annotations in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/predector/Prothint/*annotations3.tsv); do
 ID=$(echo $Annotations | cut -d '/' -f7 | cut -c 1-3)_$(echo $Annotations | cut -d '/' -f8 | cut -c 1-3)_$(echo $Annotations | cut -d '/' -f9)
 Predictedgenes=$(cut -f1 $Annotations| grep -v 'name' | wc -l)
 pfamID=$(awk '$17 != "." {print $1}' $Annotations| grep -v 'name' | wc -l)
 pfamvirulence=$(awk '$19 == "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 PHIbaseeffector=$(awk '$14 == "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 effectorortholog=$(awk '$10 != "." {print $1}' $Annotations| grep -v 'name' | wc -l)
 canonicalsp=$(awk '$30 == "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 noncanonical=$(awk '$30 != "1" && $5 != 'signal' && $2 > 0.6 {print $1}' $Annotations| grep -v 'name' | wc -l)
 sscp=$(awk '$30 == "1" && $35 <= 300 && $38 > 3 {print $1}' $Annotations| grep -v 'name' | wc -l)
 csep=$(awk '$30 == "1" && $6 == "0" {print $1}' $Annotations| grep -v 'name' | wc -l)
 secretedeffector=$(awk '$30 == "1" && $26 == "." && $6 == "0" {print $1}' $Annotations| grep -v 'name' | wc -l)
 secretedcazy=$(awk '$30 == "1" && $20 != "." {print $1}' $Annotations| grep -v 'name' | wc -l)
 newsscp=$(awk '$30 == "1" && $35 <= 300 && $38 > 3 && $10 == "." && $14 != "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 newcsep=$(awk '$30 == "1" && $6 == "0" && $10 == "." && $14 != "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 newsecretedeffector=$(awk '$30 == "1" && $26 == "." && $10 == "." && $14 != "1" && $6 == "0" {print $1}'  $Annotations| grep -v 'name' | wc -l)
 all=$(awk '$14 == "1" || $10 != "." || ($30 == "1" && $26 == "." && $6 == "0") {print $1}' $Annotations| grep -v 'name' | wc -l)
 csep_val=$(awk '$10 != "." || ($30 == "1" && $26 == "." && $6 == "0") {print $1}' $Annotations| grep -v 'name' | wc -l)
 PHIbaseororthologeffector=$(awk '$14 == "1" || $10 != "." || ($30 == "1" && $26 == "." && $6 == "0") {print $1}' $Annotations| grep -v 'name' | wc -l)
 echo $ID $Predictedgenes $pfamID $pfamvirulence $PHIbaseeffector $effectorortholog $canonicalsp $noncanonical $sscp $csep $secretedeffector $secretedcazy $newsscp $newcsep $newsecretedeffector $PHIbaseororthologeffector >> EffectorPredictionReport.txt
done

echo Species+Assembly Predicted_genes pfamID pfam_virulence PHIbase_effector effector_ortholog canonical_SP noncanonical_SP SSCP CSEP Secreted_effectorP3 Secreted_CAZY newSSCP newCSEP newSecreted_effectorP3 csep_val all csep_val phiboref > AbinitioEffectorPredictionReport.txt
for Annotations in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/predector/abinitio/*annotations3.tsv); do
 ID=$(echo $Annotations | cut -d '/' -f7 | cut -c 1-3)_$(echo $Annotations | cut -d '/' -f8 | cut -c 1-3)_$(echo $Annotations | cut -d '/' -f9)
 Predictedgenes=$(cut -f1 $Annotations| grep -v 'name' | wc -l)
 pfamID=$(awk '$17 != "." {print $1}' $Annotations| grep -v 'name' | wc -l)
 pfamvirulence=$(awk '$19 == "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 PHIbaseeffector=$(awk '$14 == "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 effectorortholog=$(awk '$10 != "." {print $1}' $Annotations| grep -v 'name' | wc -l)
 canonicalsp=$(awk '$30 == "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 noncanonical=$(awk '$30 != "1" && $5 != 'signal' && $2 > 0.6 {print $1}' $Annotations| grep -v 'name' | wc -l)
 sscp=$(awk '$30 == "1" && $35 <= 300 && $38 > 3 {print $1}' $Annotations| grep -v 'name' | wc -l)
 csep=$(awk '$30 == "1" && $6 == "0" {print $1}' $Annotations| grep -v 'name' | wc -l)
 secretedeffector=$(awk '$30 == "1" && $26 == "." && $6 == "0" {print $1}' $Annotations| grep -v 'name' | wc -l)
 secretedcazy=$(awk '$30 == "1" && $20 != "." {print $1}' $Annotations| grep -v 'name' | wc -l)
 newsscp=$(awk '$30 == "1" && $35 <= 300 && $38 > 3 && $10 == "." && $14 != "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 newcsep=$(awk '$30 == "1" && $6 == "0" && $10 == "." && $14 != "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 newsecretedeffector=$(awk '$30 == "1" && $26 == "." && $10 == "." && $14 != "1" && $6 == "0" {print $1}'  $Annotations| grep -v 'name' | wc -l)
 all=$(awk '$14 == "1" || $10 != "." || ($30 == "1" && $26 == "." && $6 == "0") {print $1}' $Annotations| grep -v 'name' | wc -l)
 csep_val=$(awk '$10 != "." || ($30 == "1" && $26 == "." && $6 == "0") {print $1}' $Annotations| grep -v 'name' | wc -l)
 PHIbaseororthologeffector=$(awk '$14 == "1" || $10 != "." || ($30 == "1" && $26 == "." && $6 == "0") {print $1}' $Annotations| grep -v 'name' | wc -l)
 echo $ID $Predictedgenes $pfamID $pfamvirulence $PHIbaseeffector $effectorortholog $canonicalsp $noncanonical $sscp $csep $secretedeffector $secretedcazy $newsscp $newcsep $newsecretedeffector $csep_val $all >> AbinitioEffectorPredictionReport.txt
done

screen -S busco
conda activate BUSCO
for assembly in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/*_clean.fasta); do
      Jobs=$(squeue -u theaven |grep 'busco'| wc -l)
    echo x
    while [ $Jobs -gt 5 ]; do
      sleep 120s
      printf "."
      Jobs=$(squeue -u theaven |grep 'busco'| wc -l)
    done
  ID=$(echo $assembly | cut -d '/' -f7 | cut -c 1-3)_$(echo $assembly | cut -d '/' -f8 | cut -c 1-3)_$(echo $assembly | cut -d '/' -f9)
  OutDir=$(dirname $assembly)/gene_pred/braker/final-02042924/predector/Prothint/results/BUSCO
  echo $OutDir
  ProgDir=~/scratch/apps/busco
  BuscoDB=ascomycota_odb10
  echo $BuscoDB
  if [ ! -d $OutDir ]
  then
    echo running for $ID
    echo $assembly 2>&1 >> buscolog.txt
    sbatch $ProgDir/sub_busco.sh $assembly $BuscoDB $OutDir 2>&1 >> buscolog.txt
  else
    echo already run
  fi
done


for file in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/predector/Prothint/results/BUSCO/ascomycota_odb10/*/short_summary*.txt); do
ID=$(echo $file | cut -d '/' -f7 | cut -c 1-3)_$(echo $file | cut -d '/' -f8 | cut -c 1-3)_$(echo $file | cut -d '/' -f9)
percent=$(grep 'C:' $file)
no1=$(grep 'Complete and single-copy BUSCOs (S)' $file | sed 's@Complete and single-copy BUSCOs (S)@@g')
cent1=$(echo $percent|cut -d ':' -f2 | sed 's@\[S@@g')
no2=$(grep 'Complete and duplicated BUSCOs (D)' $file | sed 's@Complete and duplicated BUSCOs (D)@@g')
cent2=$(echo $percent|cut -d ':' -f4 | sed 's@\],F@@g')
no3=$(grep 'Missing BUSCOs (M)' $file| sed 's@Missing BUSCOs (M)@@g')
cent3=$(echo $percent|cut -d ':' -f6 | sed 's@,n@@g')
echo $ID $no1 $cent1 $no2 $cent2 $no3 $cent3 >> BUSCOreportallmildews.txt
done
conda deactivate
exit
echo finished
```
#### Orthofinder
```bash
for Annotations in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/predector/Prothint/*annotations3.tsv); do
fasta=$(echo $Annotations | cut -d '/' -f1,2,3,4,5,6,7,8,9,10,11,12,13)/clean_final_genes_renamed.pep.fasta
out=$(dirname $Annotations)/cseps-xtra.fasta
awk '($30 == "1" && $6 == "0") || $10 ~ /BgtAVRa10|BgtAVRk1/ {print $1}' $Annotations > ids.txt
python3 /home/theaven/scratch/apps/tools/seq_get.py --id_file ids.txt --input $fasta --output $out
done

for Annotations in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/predector/abinitio/*annotations3.tsv); do
fasta=$(echo $Annotations | cut -d '/' -f1,2,3,4,5,6,7,8,9,10,11,12,13)/clean_ab_initio_final_genes_renamed.pep.fasta
out=$(dirname $Annotations)/ab_initio_cseps-xtra.fasta
awk '($30 == "1" && $6 == "0") || $10 ~ /BgtAVRa10|BgtAVRk1/ {print $1}' $Annotations > ids2.txt
python3 /home/theaven/scratch/apps/tools/seq_get.py --id_file ids2.txt --input $fasta --output $out
done

ProjDir=/home/theaven/scratch/uncompressed/analysis/orthology/orthofinder
cd $ProjDir
IsolateAbrv=All_mildew_cseps
WorkDir=/home/theaven/scratch/uncompressed/analysis/orthology/orthofinder/$IsolateAbrv
mkdir -p $WorkDir
mkdir -p $WorkDir/formatted
mkdir -p $WorkDir/goodProteins
mkdir -p $WorkDir/badProteins  
cd $WorkDir/formatted

for Fasta_file in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/predector/Prothint/cseps-xtra.fasta); do
Id_field=1
ID=$(echo $Fasta_file | cut -d '/' -f7 | cut -c 1-3)_$(echo $Fasta_file | cut -d '/' -f8 | cut -c 1-3)_$(echo $Fasta_file | cut -d '/' -f9)
~/git_repos/Scripts/gruffalo/orthomclAdjustFasta.pl $ID $Fasta_file $Id_field
done

cd $ProjDir

Input_dir=$WorkDir/formatted
Min_length=10
Max_percent_stops=20
Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
~/git_repos/Scripts/gruffalo/orthomclFilterFasta.pl $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file


conda activate orthofinder
Prefix=$IsolateAbrv
OutDir=$IsolateAbrv
sbatch ~/git_repos/Wrappers/gruffalo/run_orthofinder.sh $Input_dir $Prefix $OutDir
#23731163

ProjDir=/home/theaven/scratch/uncompressed/analysis/orthology/orthofinder
cd $ProjDir
IsolateAbrv=All_mildew_cseps_abinitio
WorkDir=/home/theaven/scratch/uncompressed/analysis/orthology/orthofinder/$IsolateAbrv
mkdir -p $WorkDir
mkdir -p $WorkDir/formatted
mkdir -p $WorkDir/goodProteins
mkdir -p $WorkDir/badProteins  
cd $WorkDir/formatted

for Fasta_file in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/predector/abinitio/ab_initio_cseps-xtra.fasta); do
Id_field=1
ID=$(echo $Fasta_file | cut -d '/' -f7 | cut -c 1-3)_$(echo $Fasta_file | cut -d '/' -f8 | cut -c 1-3)_$(echo $Fasta_file | cut -d '/' -f9)
~/git_repos/Scripts/gruffalo/orthomclAdjustFasta.pl $ID $Fasta_file $Id_field
done

cd $ProjDir

Input_dir=$WorkDir/formatted
Min_length=10
Max_percent_stops=20
Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
~/git_repos/Scripts/gruffalo/orthomclFilterFasta.pl $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file


conda activate orthofinder
Prefix=$IsolateAbrv
OutDir=$IsolateAbrv
sbatch ~/git_repos/Wrappers/gruffalo/run_orthofinder.sh $Input_dir $Prefix $OutDir
#23730145

ProjDir=/home/theaven/scratch/uncompressed/analysis/orthology/orthofinder
cd $ProjDir
IsolateAbrv=All_mildew_prothint
WorkDir=/home/theaven/scratch/uncompressed/analysis/orthology/orthofinder/$IsolateAbrv
mkdir -p $WorkDir
mkdir -p $WorkDir/formatted
mkdir -p $WorkDir/goodProteins
mkdir -p $WorkDir/badProteins  
cd $WorkDir/formatted

for Fasta_file in $(ls /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000151065.3/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153116/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153117/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153118/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153119/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153120/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Parauncinula/polyspora/Parp01/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/fusca/GCA_030378345.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/xanthii/GCA_028751805.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000417025.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000417865.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000418435.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000441875.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_900519115.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_905067625.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_000401675.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_900237765.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_900239735.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_900638725.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/alphitoides/CLCBIO/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798715.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798735.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798755.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798775.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798795.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_016906895.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/neolycopersici/GCA_003610855.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/pisi/GCA_000208805.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/pisi/GCA_000214055.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/pulchra/GCA_002918395.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/cichoracearum/GCA_003611195.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/cichoracearum/GCA_003611215.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/cichoracearum/GCA_003611235.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/magnicellulatus/GCA_006912115.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/orontii/MGH1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Leveillula/taurrica/CADEPA01/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Phyllactinia/moricola/GCA_019455665.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Pleochaeta/shiraiana/GCA_019455505.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/DRT72020/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/DRT72021/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/SCOTT2020/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/cerasii/GCA_018398735.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/GCA_013170925.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGB2019/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGB2021/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGBp112020/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/xanthii/GCA_010015925.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/xanthii/GCA_014884795.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Oidium/heveae/GCA_003957845.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_024703715.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta); do
Id_field=1
ID=$(echo $Fasta_file | cut -d '/' -f7 | cut -c 1-3)_$(echo $Fasta_file | cut -d '/' -f8 | cut -c 1-3)_$(echo $Fasta_file | cut -d '/' -f9)
~/git_repos/Scripts/gruffalo/orthomclAdjustFasta.pl $ID $Fasta_file $Id_field
done

cd $ProjDir

Input_dir=$WorkDir/formatted
Min_length=10
Max_percent_stops=20
Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
~/git_repos/Scripts/gruffalo/orthomclFilterFasta.pl $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file


conda activate orthofinder
Prefix=$IsolateAbrv
OutDir=$IsolateAbrv
sbatch ~/git_repos/Wrappers/gruffalo/run_orthofinder.sh $Input_dir $Prefix $OutDir
#23726846

ProjDir=/home/theaven/scratch/uncompressed/analysis/orthology/orthofinder
cd $ProjDir
IsolateAbrv=All_mildew_ab_initio
WorkDir=/home/theaven/scratch/uncompressed/analysis/orthology/orthofinder/$IsolateAbrv
mkdir -p $WorkDir
mkdir -p $WorkDir/formatted
mkdir -p $WorkDir/goodProteins
mkdir -p $WorkDir/badProteins  
cd $WorkDir/formatted

for Fasta_file in $(ls /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000151065.3/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153116/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153117/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153118/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153119/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis-secale/SRR2153120/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Parauncinula/polyspora/Parp01/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/fusca/GCA_030378345.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/xanthii/GCA_028751805.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000417025.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000417865.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000418435.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_000441875.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_900519115.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/graminis/GCA_905067625.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_000401675.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_900237765.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_900239735.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Blumeria/hordei/GCA_900638725.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/alphitoides/CLCBIO/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798715.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798735.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798755.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798775.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_000798795.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_016906895.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/neolycopersici/GCA_003610855.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/pisi/GCA_000208805.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/pisi/GCA_000214055.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/pulchra/GCA_002918395.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/cichoracearum/GCA_003611195.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/cichoracearum/GCA_003611215.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/cichoracearum/GCA_003611235.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/magnicellulatus/GCA_006912115.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Golovinomyces/orontii/MGH1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Leveillula/taurrica/CADEPA01/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Phyllactinia/moricola/GCA_019455665.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Pleochaeta/shiraiana/GCA_019455505.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/DRT72020/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/DRT72021/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/SCOTT2020/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/cerasii/GCA_018398735.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/GCA_013170925.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGB2019/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGB2021/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGBp112020/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/xanthii/GCA_010015925.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Podosphaera/xanthii/GCA_014884795.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Oidium/heveae/GCA_003957845.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta /home/theaven/scratch/uncompressed/genomes/Erysiphe/necator/GCA_024703715.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta | sed 's@clean_final_genes_renamed.pep.fasta@clean_ab_initio_final_genes_renamed.pep.fasta@g'); do
#echo $Fasta_file
#grep '>;' $Fasta_file
Id_field=1
ID=$(echo $Fasta_file | cut -d '/' -f7 | cut -c 1-3)_$(echo $Fasta_file | cut -d '/' -f8 | cut -c 1-3)_$(echo $Fasta_file | cut -d '/' -f9)
~/git_repos/Scripts/gruffalo/orthomclAdjustFasta.pl $ID $Fasta_file $Id_field
done

cd $ProjDir

Input_dir=$WorkDir/formatted
Min_length=10
Max_percent_stops=20
Good_proteins_file=$WorkDir/goodProteins/goodProteins.fasta
Poor_proteins_file=$WorkDir/badProteins/poorProteins.fasta
~/git_repos/Scripts/gruffalo/orthomclFilterFasta.pl $Input_dir $Min_length $Max_percent_stops $Good_proteins_file $Poor_proteins_file


conda activate orthofinder
Prefix=$IsolateAbrv
OutDir=$IsolateAbrv
sbatch ~/git_repos/Wrappers/gruffalo/run_orthofinder.sh $Input_dir $Prefix $OutDir
#23726964
```




















## MCAGs

#### Kraken
Run the 'clean' genomes through kraken to further screen for contaminants
```bash
#Kraken 
for Assembly in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/Genomes/*/*/*/*.fna /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/Genomes/*/*/*/*.fasta); do
sleep 30s
      while  squeue -u did23faz | grep -q '(AssocMaxJobsLimit)\|(Resources)\|(Priority)'; do
        sleep 300s
        printf "."
    done
Jobs=$(squeue -u did23faz| grep 'kraken'  | wc -l)
      while [ $Jobs -gt 4 ]; do
          sleep 300s
            printf "."
          Jobs=$(squeue -u did23faz| grep 'kraken'| wc -l)
      done
OutPrefix=$(echo $Assembly | cut -d '/' -f9,10,11 | sed 's@/@_@g')_kraken2nt
OutDir=$(dirname $Assembly)/kraken2.1.3
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/kraken/nt_14092023
if [ ! -e "${OutDir}/${OutPrefix}_report.txt" ] || [ ! -s "${OutDir}/${OutPrefix}_report.txt" ]; then
echo $OutPrefix
mkdir $OutDir
sbatch $ProgDir/run_kraken2.sh $Assembly $Database $OutDir $OutPrefix
else
echo already done for $OutPrefix
fi
done

for report in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/Genomes/*/*/*/kraken2.1.3/*_kraken2nt_report.txt); do
echo $(basename $report | sed 's@_kraken2nt_report.txt@@g')
grep 'Fungi' $report
grep 'Erysiphales' $report
done
```
### Screen for MCAGs
#### BLAST - search for MCAGs
```bash
#Give species unique names to each predicted gene:
mkdir /home/theaven/scratch/uncompressed/mcag/db/99
cd /home/theaven/scratch/uncompressed/mcag/db/99
for file in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta); do
ID=$(echo $file | cut -d '/' -f7 | cut -c 1-3)_$(echo $file | cut -d '/' -f8 | cut -c 1-3)_$(echo $file | cut -d '/' -f9)
echo $ID
cp $file ./${ID}_$(basename $file)
sed -i "s/>/>${ID}_/g" ${ID}_$(basename $file)
done

for file in $(ls *final_genes_renamed.pep.fasta); do
seqkit rmdup $file > temp.fasta && mv temp.fasta $file
done

#make a blast database from the proteomes and query each putative MCAG against it:
mkdir db
cat *final_genes_renamed.pep.fasta > db/db.faa
cd db
makeblastdb -in db.faa -input_type fasta -dbtype prot  -title mildew+db  -parse_seqids -out mildew+db

for file in $(ls /home/theaven/scratch/uncompressed/mcag/prot/*.fsa); do
name=$(basename $file | cut -d '.' -f1)
blastp -query $file -db mildew+db -out ../${name}_results -evalue 1e-5 -outfmt 6 -num_threads 1
done

ls ../*_results | wc -l

#Collect results into one file:
echo Gene: > out99.tsv

for file2 in $(ls ../*_final_genes_renamed.pep.fasta); do
  name=$(basename $file2 | sed 's@_clean_final_genes_renamed.pep.fasta@@g')
  sed -i "\$s/\$/ $name/" out99.tsv
done

for file in $(ls ../*_results); do
gene=$(awk '{print $1; exit}' $file)
echo $gene >> out99.tsv
for file2 in $(ls ../*_final_genes_renamed.pep.fasta); do
  name=$(basename $file2 | sed 's@_clean_final_genes_renamed.pep.fasta@@g')
x=$(grep -m 1 $name $file | awk '{print $11}')
if [ -z "$x" ]; then
    y=missing
else
    y=$x
fi
sed -i "\$s/\$/ $y/" out99.tsv
done
done
```

### Screen for new MCAGs

#### BLAST - Find CAGs Shared by B.cinerea, M. oryzae, C. higginsianum, S. sclerotiorum and S. cerevisiae
```bash
cd /home/theaven/scratch/uncompressed/mcag/db
makeblastdb -in /home/theaven/scratch/uncompressed/genomes/Botrytis/cineria/GCA_000143535.4/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta -input_type fasta -dbtype prot  -title Bot_cin  -parse_seqids -out Bot_cin
makeblastdb -in /home/theaven/scratch/uncompressed/genomes/Magnaporthe/oryzae/GCF_000002495.2/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta -input_type fasta -dbtype prot  -title Mag_ory  -parse_seqids -out Mag_ory
makeblastdb -in /home/theaven/scratch/uncompressed/genomes/Colletotrichum/higginsianum/GCA_001672515.1/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta -input_type fasta -dbtype prot  -title Col_hig  -parse_seqids -out Col_hig
makeblastdb -in /home/theaven/scratch/uncompressed/genomes/Sclerotinia/sclerotiorum/GCA_000146945.2/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta -input_type fasta -dbtype prot  -title Scl_scl  -parse_seqids -out Scl_scl
makeblastdb -in /home/theaven/scratch/uncompressed/genomes/Saccharomyces/cerevisiae/GCF_000146045.2/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta -input_type fasta -dbtype prot  -title Sac_cer  -parse_seqids -out Sac_cer

blastp -query /home/theaven/scratch/uncompressed/genomes/Saccharomyces/cerevisiae/GCF_000146045.2/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta -db Scl_scl -out ../Scl_scl_results -evalue 1e-5 -outfmt 6 -num_threads 1 -max_target_seqs 1
awk '{print $1}' ../Scl_scl_results | sort | uniq > hits.txt
python3 /home/theaven/scratch/apps/tools/seq_get.py --id_file hits.txt --input /home/theaven/scratch/uncompressed/genomes/Saccharomyces/cerevisiae/GCF_000146045.2/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta --output Sac_cer_Scl_scl.faa
blastp -query Sac_cer_Scl_scl.faa -db Col_hig -out ../Col_hig_results -evalue 1e-5 -outfmt 6 -num_threads 1 -max_target_seqs 1
awk '{print $1}' ../Col_hig_results | sort | uniq > hits.txt
python3 /home/theaven/scratch/apps/tools/seq_get.py --id_file hits.txt --input Sac_cer_Scl_scl.faa --output Sac_cer_Scl_scl_Col_hig.faa
blastp -query Sac_cer_Scl_scl_Col_hig.faa -db Mag_ory -out ../Mag_ory_results -evalue 1e-5 -outfmt 6 -num_threads 1 -max_target_seqs 1
awk '{print $1}' ../Mag_ory_results | sort | uniq > hits.txt
python3 /home/theaven/scratch/apps/tools/seq_get.py --id_file hits.txt --input Sac_cer_Scl_scl_Col_hig.faa --output Sac_cer_Scl_scl_Col_hig_Mag_ory.faa
blastp -query Sac_cer_Scl_scl_Col_hig_Mag_ory.faa -db Bot_cin -out ../Bot_cin_results -evalue 1e-5 -outfmt 6 -num_threads 1 -max_target_seqs 1
awk '{print $1}' ../Bot_cin_results | sort | uniq > hits.txt
python3 /home/theaven/scratch/apps/tools/seq_get.py --id_file hits.txt --input Sac_cer_Scl_scl_Col_hig_Mag_ory.faa --output Sac_cer_Scl_scl_Col_hig_Mag_ory_Bot_cin.faa
grep '>' Sac_cer_Scl_scl_Col_hig_Mag_ory_Bot_cin.faa | wc -l #4,048

blastp -query /home/theaven/scratch/uncompressed/genomes/Saccharomyces/cerevisiae/GCF_000146045.2/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta -db Scl_scl -out ../Scl_scl_results -evalue 1e-10 -outfmt 6 -num_threads 1 -max_target_seqs 1
awk '{print $1}' ../Scl_scl_results | sort | uniq > hits.txt
python3 /home/theaven/scratch/apps/tools/seq_get.py --id_file hits.txt --input /home/theaven/scratch/uncompressed/genomes/Saccharomyces/cerevisiae/GCF_000146045.2/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta --output Sac_cer_Scl_scl.faa
blastp -query Sac_cer_Scl_scl.faa -db Col_hig -out ../Col_hig_results -evalue 1e-10 -outfmt 6 -num_threads 1 -max_target_seqs 1
awk '{print $1}' ../Col_hig_results | sort | uniq > hits.txt
python3 /home/theaven/scratch/apps/tools/seq_get.py --id_file hits.txt --input Sac_cer_Scl_scl.faa --output Sac_cer_Scl_scl_Col_hig.faa
blastp -query Sac_cer_Scl_scl_Col_hig.faa -db Mag_ory -out ../Mag_ory_results -evalue 1e-10 -outfmt 6 -num_threads 1 -max_target_seqs 1
awk '{print $1}' ../Mag_ory_results | sort | uniq > hits.txt
python3 /home/theaven/scratch/apps/tools/seq_get.py --id_file hits.txt --input Sac_cer_Scl_scl_Col_hig.faa --output Sac_cer_Scl_scl_Col_hig_Mag_ory.faa
blastp -query Sac_cer_Scl_scl_Col_hig_Mag_ory.faa -db Bot_cin -out ../Bot_cin_results -evalue 1e-10 -outfmt 6 -num_threads 1 -max_target_seqs 1
awk '{print $1}' ../Bot_cin_results | sort | uniq > hits.txt
python3 /home/theaven/scratch/apps/tools/seq_get.py --id_file hits.txt --input Sac_cer_Scl_scl_Col_hig_Mag_ory.faa --output Sac_cer_Scl_scl_Col_hig_Mag_ory_Bot_cin2.faa
grep '>' Sac_cer_Scl_scl_Col_hig_Mag_ory_Bot_cin2.faa | wc -l #3,660
```
#### BLAST - search mildew proteomes for this new CAG set
```bash
screen -S mcag
srun -p short  -c 1 --mem 50G --pty bash
mkdir -p test/out
cd /home/theaven/scratch/uncompressed/mcag/db/test
for proteome in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta); do
ID=$(echo $proteome | cut -d '/' -f7 | cut -c 1-3)_$(echo $proteome | cut -d '/' -f8 | cut -c 1-3)_$(echo $proteome | cut -d '/' -f9)
echo $ID
makeblastdb -in $proteome -input_type fasta -dbtype prot  -title $ID -parse_seqids -out $ID
blastp -query ../Sac_cer_Scl_scl_Col_hig_Mag_ory_Bot_cin2.faa -db $ID -out out/${ID}_results -evalue 1e-6 -outfmt 6 -num_threads 32 -max_target_seqs 1
done

cd out
grep '>' ../../Sac_cer_Scl_scl_Col_hig_Mag_ory_Bot_cin2.faa | sed 's@>@@g' > CAGs.txt
echo CAGs > mcags.txt
cat CAGs.txt >> mcags.txt
for proteome in *_results; do
    echo "$proteome" > temp.txt
    while IFS= read -r line; do
      echo $line
        # Check if the line exists in the current proteome file
        if grep -qF "$line" "$proteome"; then
            echo "1" >> temp.txt
        else
            echo "0" >> temp.txt
        fi
    done < "CAGs.txt"
    paste mcags.txt temp.txt > temp2.txt && mv temp2.txt mcags.txt
done

head -n 1 mcags.txt > MCAGs.tsv
awk -F'\t' 'BEGIN {OFS="\t"} {count=0; for (i=7; i<=22; i++) if ($i == "0") count++; for (i=30; i<=50; i++) if ($i == "0") count++; for (i=56; i<=74; i++) if ($i == "0") count++; if (count >= 3) print}' mcags.txt >> MCAGs.tsv
#572 are missing from at least 3 mildews
#927 are missing from at least 2 mildews
#1,622 are missing from at least 1 mildew
```
#### Look for overlap between existing and new CAG sets
```bash
cd /home/theaven/scratch/uncompressed/mcag/db/test
makeblastdb -in /home/theaven/scratch/uncompressed/mcag/db/Sac_cer_Scl_scl_Col_hig_Mag_ory_Bot_cin2.faa -input_type fasta -dbtype prot  -title spanu -parse_seqids -out spanu
for file in $(ls /home/theaven/scratch/uncompressed/mcag/prot/*); do
name=$(basename $file | cut -d '.' -f1)
blastp -query $file -db spanu -out ../${name}_results -evalue 1e-6 -outfmt 6 -num_threads 1 
done

cat ../*results > temp.txt
awk -F'\t' '$3 > 90' temp.txt > temp2.txt
```
#### Check to see if MCAGs that are present blast to distantly related species and are therefore likely contaminants rather than truely present MCAGs:
```bash
#Collect gne IDs for select de novo mcags:
mkdir /home/theaven/scratch/uncompressed/mcag/db/blast
extramcags=/home/theaven/scratch/uncompressed/mcags+.txt
for blastresult in /home/theaven/scratch/uncompressed/mcag/db/test/out/*_results; do
  assembly_id=$(basename "$blastresult" | sed 's@_results@@g')
  awk 'NR==FNR { patterns[$1]; next } $1 in patterns && !($1 in printed) { print; printed[$1] }' "$extramcags" "$blastresult" | awk '{print $2}' > temp_genenames.txt
  while read -r gene; do
    echo "${assembly_id}_${gene}" >> "/home/theaven/scratch/uncompressed/mcag/db/blast/${assembly_id}_proteins.txt"
  done < temp_genenames.txt
  rm temp_genenames.txt
done

#Collect gene IDs for 99 predefined MCAGS:
ls /home/theaven/scratch/uncompressed/mcag/db/99/*_final_genes_renamed.pep.fasta | rev | cut -d '/' -f1 |  rev | sed 's@_clean_final_genes_renamed.pep.fasta@@g' > assemblies.txt
for assembly_id in $(cat assemblies.txt); do
for protein in $(ls /home/theaven/scratch/uncompressed/mcag/db/99/*_protein_results); do
grep -m 1 "$assembly_id" $protein | awk '{print $2}' >> /home/theaven/scratch/uncompressed/mcag/db/blast/${assembly_id}_proteins.txt
done
done

#Get the corresponding protein sequences for the mcags
for file in $(ls /home/theaven/scratch/uncompressed/mcag/db/blast/*_proteins.txt); do
sort $file > temp-temp.txt && mv temp-temp.txt $file
done

for file in $(ls /home/theaven/scratch/uncompressed/mcag/db/blast/*_proteins.txt); do
  ID=$(basename $file | sed 's@_proteins.txt@@g')
  ls /home/theaven/scratch/uncompressed/mcag/db/99/${ID}_clean_final_genes_renamed.pep.fasta
  python3 /home/theaven/scratch/apps/tools/seq_get.py --id_file $file --input /home/theaven/scratch/uncompressed/mcag/db/99/${ID}_clean_final_genes_renamed.pep.fasta --output /home/theaven/scratch/uncompressed/mcag/db/blast/${ID}_final_genes_renamed.pep.fasta
done

#Get search the mcag proteins against the uniprot database
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/MCAG/blast2/*_final_genes_renamed.pep.fasta); do
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/databases/blast/uniprot_01102023/Uniprot_01102023_reference_proteomes.dmnd
OutDir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/MCAG/blast2
OutFile=$(basename $file | sed 's@_final_genes_renamed.pep.fasta@@g')
hits=5
ProgDir=~/git_repos/Wrappers/NBI
sbatch $ProgDir/run_diamond_blastp.sh $file $Database $OutDir $OutFile $hits
done
#1272157-240

#Collect the outputs into a single file for each assembly:
for assembly_id in $(cat assemblies.txt); do
echo -e "MCAG\tTop_hit\tblast_1\tblast_2\tblast_3\tblast_4\tblast_5\tID_1\tID_2\tID_3\tID_4\tID_5" > /home/theaven/scratch/uncompressed/mcag/db/blast/${assembly_id}_proteins_notes.txt
for protein in $(ls /home/theaven/scratch/uncompressed/mcag/db/99/*_protein_results); do
x=$(basename $protein | sed 's@_protein_results@@g'| sed 's@results/@@g') 
y=$(grep -m 1 "$assembly_id" $protein | awk '{print $2}')
echo -e "$x\t$y" >> /home/theaven/scratch/uncompressed/mcag/db/blast/${assembly_id}_proteins_notes.txt
done
done

extramcags=/home/theaven/scratch/uncompressed/mcags+.txt
mapfile -t patterns < "$extramcags"
for blastresult in /home/theaven/scratch/uncompressed/mcag/db/test/out/*_results; do
  assembly_id=$(basename "$blastresult" | sed 's@_results@@g')
  
  for pattern in "${patterns[@]}"; do
    y=$(awk -v pat="$pattern" '$1 == pat { print $2; exit }' "$blastresult")
    if [ -n "$y" ]; then
      echo -e "${pattern}\t${assembly_id}_${y}" >> "/home/theaven/scratch/uncompressed/mcag/db/blast/${assembly_id}_proteins_notes.txt"
    fi
  done
done

for file in $(ls /home/theaven/scratch/uncompressed/mcag/db/blast/*_proteins_notes.txt); do
awk 'BEGIN {FS=OFS="\t"} $2 == "" {$2="###"}1' $file > temp.txt && mv temp.txt $file
done

for file in /home/theaven/scratch/uncompressed/mcag/db/blast/*_proteins_notes.txt; do
    ID=$(basename "$file" | sed 's@_proteins_notes.txt@@g')
    while IFS=$'\t' read -r col1 col2 rest; do
        matches=$(grep "^$col2" "/home/theaven/scratch/uncompressed/mcag/db/blast/${ID}.vs.Uniprot_01102023_reference_proteomes.dmnd.diamond_blastp.out" | awk -v col2="$col2" '$1 == col2 {printf "%s\t", $5}')
        if [ -n "$matches" ]; then
            matches="${matches%?}"  # Remove the last character (trailing tab)
        fi
        echo -e "$col1\t$col2\t$rest\t$matches"
    done < "$file" > "${ID}merged_notes.txt"
    awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $4, $5, $6, $7, $8, $9}' ${ID}merged_notes.txt > $file
    rm ${ID}merged_notes.txt
done

#Get the uniprot codes for lookup:
for file in /home/theaven/scratch/uncompressed/mcag/db/blast/*_proteins_notes.txt; do
awk -F'\t' '{split($3, a, "|"); split($4, b, "|"); split($5, c, "|"); split($6, d, "|"); split($7, e, "|"); print a[2] "\n" b[2] "\n" c[2] "\n" d[2] "\n" e[2]}' $file | sort | uniq > $(dirname $file)/$(basename $file | sed 's@_proteins_notes.txt@@g')_protcodes.txt
done

#Edit notes file to contain species info:
for notes in $(ls /home/theaven/scratch/uncompressed/mcag/db/blast/*_proteins_notes.txt); do
tsv=/home/theaven/scratch/uncompressed/mcag/db/blast/uniprot/$(basename "$notes" | sed 's@_proteins_notes.txt@_protcodes.tsv@g')
output_file="edited_file.txt"
rm "$output_file"
touch "$output_file"
while IFS= read -r line; do
    id=$(echo "$line" | awk -F'\t' '{split($3, arr, "|"); print arr[2]}')
    fifth_column=$(awk -v var="$id" -F'\t' '$1 == var {print $5}' $tsv | cut -d ' ' -f1,2)
    edited_line="$line\t$fifth_column"
    echo -e "$edited_line" >> "$output_file"
done < $notes
mv "$output_file" $notes
rm "$output_file"
touch "$output_file"
while IFS= read -r line; do
    id=$(echo "$line" | awk -F'\t' '{split($4, arr, "|"); print arr[2]}')
    fifth_column=$(awk -v var="$id" -F'\t' '$1 == var {print $5}' $tsv | cut -d ' ' -f1,2)
    edited_line="$line\t$fifth_column"
    echo -e "$edited_line" >> "$output_file"
done < $notes
mv "$output_file" $notes
rm "$output_file"
touch "$output_file"
while IFS= read -r line; do
    id=$(echo "$line" | awk -F'\t' '{split($5, arr, "|"); print arr[2]}')
    fifth_column=$(awk -v var="$id" -F'\t' '$1 == var {print $5}' $tsv | cut -d ' ' -f1,2)
    edited_line="$line\t$fifth_column"
    echo -e "$edited_line" >> "$output_file"
done < $notes
mv "$output_file" $notes
rm "$output_file"
touch "$output_file"
while IFS= read -r line; do
    id=$(echo "$line" | awk -F'\t' '{split($6, arr, "|"); print arr[2]}')
    fifth_column=$(awk -v var="$id" -F'\t' '$1 == var {print $5}' $tsv | cut -d ' ' -f1,2)
    edited_line="$line\t$fifth_column"
    echo -e "$edited_line" >> "$output_file"
done < $notes
mv "$output_file" $notes
rm "$output_file"
touch "$output_file"
while IFS= read -r line; do
    id=$(echo "$line" | awk -F'\t' '{split($7, arr, "|"); print arr[2]}')
    fifth_column=$(awk -v var="$id" -F'\t' '$1 == var {print $5}' $tsv | cut -d ' ' -f1,2)
    edited_line="$line\t$fifth_column"
    echo -e "$edited_line" >> "$output_file"
done < $notes
mv "$output_file" $notes
done
```














## Unique Effectors


Blast for aphanis only effectors across wider mildews:
```bash
for file in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/clean_final_genes_renamed.pep.fasta); do
output_dir=$(dirname $file)
ID=$(echo $file | cut -d '/' -f7 | cut -c 1-3)_$(echo $file | cut -d '/' -f8 | cut -c 1-3)_$(echo $file | cut -d '/' -f9)
echo $ID
sed -E "/^>/ s/>[[:space:]]*/&$ID /" $file | sed 's/ \([^\t]\)/_\1/' > ${output_dir}/${ID}_clean_final_genes.pep.fasta
done

mkdir db
for file in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  cat $file >> db/clean_mildew.faa
done

conda activate seqkit
seqkit rmdup db/clean_mildew.faa > temp.fasta && mv temp.fasta db/clean_mildew.faa
conda deactivate

conda activate blast+
cd /home/theaven/scratch/uncompressed/db
makeblastdb -in clean_mildew.faa -input_type fasta -dbtype prot  -title mildew+db  -parse_seqids -out mildew+db

while IFS=$'\t' read -r line; do
    awk -F'\t' '$31 == "1" && $27 == "." && $7 == "0" {print $1}' <<< "$line" | grep -v 'name' >> aphanis_only_secreted_effectors.txt
done < ~/projects/niab/theaven/results/P_aphanis-THeavenDRCT72020_1-aphanis-ranked4.tsv

while IFS=$'\t' read -r line; do
echo $line > hits.txt
python3 /home/theaven/scratch/apps/tools/seq_get.py --id_file hits.txt --input /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta --output effectors/${line}.faa
done < aphanis_only_secreted_effectors.txt

for file in $(ls effectors/*.faa); do
name=$(basename $file | cut -d '.' -f1)
blastp -query $file -db db/mildew+db -out effectors/${name}_results -evalue 1e-5 -outfmt 6 -num_threads 16 -max_target_seqs 99999
done

echo Gene: > effectors/out.tsv

for file in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file | cut -d '/' -f7 | cut -c 1-3)_$(echo $file | cut -d '/' -f8 | cut -c 1-3)_$(echo $file | cut -d '/' -f9)
sed -i "\$s/\$/ $ID/" effectors/out.tsv
done


for file in $(ls effectors/*_results); do
gene=$(basename $file | sed 's@_results@@g')
echo $gene
echo $gene >> effectors/out.tsv
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep -m 1 $ID $file | awk '{print $11}')
if [ -z "$x" ]; then
    y=missing
else
    y=$x
fi
sed -i "\$s/\$/ $y/" effectors/out.tsv
done
done

ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/*clean.fasta

for file in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/*clean.fasta); do
output_dir=$(dirname $file)
ID=$(echo $file | cut -d '/' -f7 | cut -c 1-3)_$(echo $file | cut -d '/' -f8 | cut -c 1-3)_$(echo $file | cut -d '/' -f9)
echo $ID
sed -E "/^>/ s/>[[:space:]]*/&$ID /" $file | sed 's/ \([^\t]\)/_\1/' | awk '{print $1}' | sed 's@_Oidium_NoIndex_L002_R1_R2_paired_trimmed@@g' | awk -F'_length' '{print $1}' > ${output_dir}/${ID}_clean_renamed.fasta
done

for file in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/*clean_renamed.fasta); do
  cat $file >> db/clean_mildew.fa
done

awk 'length($0) > 50 && /^>/ {print}' db/clean_mildew.fa

conda activate seqkit
seqkit rmdup db/clean_mildew.fa > temp.fasta && mv temp.fasta db/clean_mildew.fa
conda deactivate

srun -p short --mem-per-cpu 4G --cpus-per-task 8 --pty bash
conda activate blast+
cd /home/theaven/scratch/uncompressed/db
makeblastdb -in clean_mildew.fa -input_type fasta -dbtype nucl  -title mildew+db2  -parse_seqids -out mildew+db2

for file in $(ls effectors/*.faa); do
name=$(basename $file | cut -d '.' -f1)
tblastn -query $file -db db/mildew+db2 -out effectors/${name}_resultsx -evalue 1e-5 -outfmt 6 -num_threads 16 -max_target_seqs 99999
done

g6562.t1        Pod_aph_SCOTT2020_contig_55     65.000  80      10      2       1       63      10354   10118
g12955.t1       Pod_aph_DRT72020_contig_3643    100.000 248     0       0       48      295     3388    4131    2.33e-170       514


echo Gene: >> effectors/out7.tsv

for file in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file | cut -d '/' -f7 | cut -c 1-3)_$(echo $file | cut -d '/' -f8 | cut -c 1-3)_$(echo $file | cut -d '/' -f9)
sed -i "\$s/\$/ $ID/" effectors/out7.tsv
done

#E value of top hit:
for file in $(ls effectors/*_resultsx); do
gene=$(basename $file | sed 's@_resultsx@@g')
echo $gene
echo $gene >> effectors/out2.tsv
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep -m 1 $ID $file | awk '{print $11}')
if [ -z "$x" ]; then
    y=missing
else
    y=$x
fi
sed -i "\$s/\$/ $y/" effectors/out2.tsv
done
done

#E value of top hit, filter for % query length and % identity:
for file in $(ls effectors/*_resultsx); do
gene=$(basename $file | sed 's@_resultsx@@g')
echo $gene
echo $gene >> effectors/out6.tsv
gene_file=$(ls effectors/${gene}.t1.faa)
gene_length=$(grep -v '>' $gene_file | tr -d '\n' | wc -c)
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep -m 1 $ID $file | awk '$3 > 70' | awk -v gene_length="$gene_length" '$4 / gene_length > 0.4' | awk '{print $11}')
if [ -z "$x" ]; then
    y=missing
else
    y=$x
fi
sed -i "\$s/\$/ $y/" effectors/out6.tsv
done
done

#No. of hits
for file in $(ls effectors/*_resultsx); do
gene=$(basename $file | sed 's@_resultsx@@g')
echo $gene
echo $gene >> effectors/out7.tsv
gene_file=$(ls effectors/${gene}.t1.faa)
gene_length=$(grep -v '>' $gene_file | tr -d '\n' | wc -c)
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 > 60' | awk -v gene_length="$gene_length" '$4 / gene_length > 0.4' | wc -l)
sed -i "\$s/\$/ $x/" effectors/out7.tsv
done
done


for file in $(ls effectors/*_resultsx); do
  grep $ID $file | awk '$3 > 60' | awk -v gene_length="$gene_length" '$4 / gene_length > 0.4' >> temp.out.txt
done



#Check for stop codons???
for file in $(ls effectors/*_resultsx); do
gene=$(basename $file | sed 's@_resultsx@@g')
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file >> effectors/species/${ID}2.tsv
#sed -i "s@${ID}_@@g" effectors/species/${ID}.tsv
#sed -i 's@_contig@contig@g' effectors/species/${ID}.tsv
done
done

for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
  out=$(dirname $file2)/aphanis_only2.bed
awk -F'\t' '{print $2"\t"$9"\t"$10"\t"$1}' effectors/species/${ID}2.tsv > $out
done

for file in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/aphanis_only2.bed); do
  genome=$(ls $(echo $file | cut -d '/' -f1,2,3,4,5,6,7,8,9,10)/*clean_renamed.fasta)
singularity exec ~/scratch/apps/containers/bedtools_2.31.1--hf5e1c6e_1 bedtools getfasta -fi $genome -bed $file > $(echo $file | sed 's@.bed@.fa@g')
done










for file in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/final_genes_renamed.gff3); do
  out=$(echo $file | sed 's@.gff3@.bed@g')
singularity exec ~/scratch/apps/containers/bedops_2.4.41--h4ac6f70_2 gff2bed < $file > $out
done



for file in $(ls effectors/*_resultsx); do
gene=$(basename $file | sed 's@_resultsx@@g')
echo $gene
echo $gene >> effectors/out7.tsv
gene_file=$(ls effectors/${gene}.t1.faa)
gene_length=$(grep -v '>' $gene_file | tr -d '\n' | wc -c)
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 > 60' | awk -v gene_length="$gene_length" '$4 / gene_length > 0.4' | grep 'Pod_aph_DRT72020' >> temp.out.txt
grep $ID $file | awk '$3 > 60' | awk -v gene_length="$gene_length" '$4 / gene_length > 0.4' | grep 'Pod_aph_DRT72021' >> temp2.out.txt
grep $ID $file | awk '$3 > 60' | awk -v gene_length="$gene_length" '$4 / gene_length > 0.4' | grep 'Pod_aph_SCOTT2020' >> temp3.out.txt
sed -i "\$s/\$/ $x/" effectors/out7.tsv
done
done


makeblastdb -in data/assembly/genome/fragaria/ananassa/F_ana_Camarosa_6-28-17.c.fna -input_type fasta -dbtype nucl  -title fananasadb2  -parse_seqids -out fananasadb2
for file in $(ls effectors/*.t1.faa); do
name=$(basename $file | cut -d '.' -f1)
tblastn -query $file -db fananasadb2 -out ${name}_resultsx_fragan -evalue 1e-5 -outfmt 6 -num_threads 12 -max_target_seqs 99999
done


for file in $(ls *_resultsx_fragan); do
gene=$(basename $file | sed 's@_resultsx_fragan@@g')
echo $gene
gene_file=$(ls effectors/${gene}.t1.faa)
gene_length=$(grep -v '>' $gene_file | tr -d '\n' | wc -c)
echo $file | awk '$3 > 60' | awk -v gene_length="$gene_length" '$4 / gene_length > 0.4'
done

makeblastdb -in data/assembly/genome/rubus/idaeus/AnitraCuratedED_HiC.c.fna -input_type fasta -dbtype nucl  -title ridaeusdb2  -parse_seqids -out ridaeusdb2
for file in $(ls effectors/*.t1.faa); do
name=$(basename $file | cut -d '.' -f1)
tblastn -query $file -db ridaeusdb2 -out ${name}_resultsx_ridaeus -evalue 1e-5 -outfmt 6 -num_threads 12 -max_target_seqs 99999
done

for file in $(ls *_resultsx_ridaeus); do
gene=$(basename $file | sed 's@_resultsx_ridaeus@@g')
echo $gene
gene_file=$(ls effectors/${gene}.t1.faa)
gene_length=$(grep -v '>' $gene_file | tr -d '\n' | wc -c)
echo $file | awk '$3 > 60' | awk -v gene_length="$gene_length" '$4 / gene_length > 0.4'
done
```
```bash
grep 'g10232;\|g11054;\|g11278;\|g1209;\|g12302;\|g12450;\|g12453;\|g12628;\|g12667;\|g12955;\|g13367;\|g13490;\|g1355;\|g13885;\|g14348;\|g14731;\|g14832;\|g15923;\|g2113;\|g2919;\|g2928;\|g3129;\|g324;\|g422;\|g548;\|g5605;\|g5645;\|g5721;\|g6055;\|g6450;\|g6520;\|g6562;\|g7278;\|g757;\|g8019;\|g9359;\|g10232.t1\|g11054.t1\|g11278.t1\|g1209.t1\|g12302.t1\|g12450.t1\|g12453.t1\|g12628.t1\|g12667.t1\|g12955.t1\|g13367.t1\|g13490.t1\|g1355.t1\|g13885.t1\|g14348.t1\|g14731.t1\|g14832.t1\|g15923.t1\|g2113.t1\|g2919.t1\|g2928.t1\|g3129.t1\|g324.t1\|g422.t1\|g548.t1\|g5605.t1\|g5645.t1\|g5721.t1\|g6055.t1\|g6450.t1\|g6520.t1\|g6562.t1\|g7278.t1\|g757.t1\|g8019.t1\|g9359.t1' /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3 > effectors/aphanis_only.gff
```
Aphanis_Gene: Exon_no.: Lengths:
g10232  1 1289
g11054  4 321 506 437 157
g11278  3 160 407 873
g1209 2 253 228
g12302  2 115 378
g12450  1 764
g12453  2 483 37
g12628  1 257
g12667* 2 142 141
g12955  2 141 745
g13367  4 157 407 509 321
g13490  1 500
g1355 1 491
g13885  1 647
g14348  1 899
g14731  2 439 123
g14832  1 786
g15923  1 611
g2113 2 479 44
g2919 1 230
g2928 1 605
g3129 1 446
g324  1 539
g422  3 138 226 668
g548  1 332
g5605 2 274 357
g5645 1 665
g5721 2 109 360
g6055 3 68 847 363
g6450 1 1232
g6520*  2 463 66
g6562 2 80 183
g7278 3 13 71 369
g757  1 386
g8019 2 152 386
g9359 1 389

```bash
#No. of hits
Output=effectors/out13.tsv
for file in $(ls effectors/*_resultsxyz); do
gene=$(basename $file | sed 's@_resultsxyz@@g')
echo $gene
echo $gene >> $Output

if [[ $gene == "g10232" || $gene == "g12450" || $gene == "g12628" || $gene == "g13490" || $gene == "g1355" || $gene == "g13885" || $gene == "g14348" || $gene == "g14832" || $gene == "g15923" || $gene == "g2919" || $gene == "g2928" || $gene == "g3129" || $gene == "g324" || $gene == "g548" || $gene == "g5645" || $gene == "g6450" || $gene == "g757" || $gene == "g9359" ]]; then
gene_file=$(ls effectors/${gene}.t1.faa)
gene_length=$(grep -v '>' $gene_file | tr -d '\n' | wc -c)
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk -v gene_length="$gene_length" '$4 / gene_length > 0.7' | awk '$11 < 1e-23' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g11054" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 118' | awk '$11 < 1e-23' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g11278" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 203' | awk '$11 < 1e-23' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g1209" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 59' | awk '$11 < 1e-23' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g12302" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 88' | awk '$11 < 1e-23' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g12453" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 112' | awk '$11 < 1e-23' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g12667" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 32' | awk '$11 < 1e-23' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g12955" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 173' | awk '$11 < 1e-23' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g13367" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 118' | awk '$11 < 1e-23' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g14731" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 102' | awk '$11 < 1e-23' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g2113" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 111' | awk '$11 < 1e-23' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g422" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 155' | awk '$11 < 1e-23' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g5605" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 83' | awk '$11 < 1e-23' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g5721" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 84' | awk '$11 < 1e-23' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g6055" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 197' | awk '$11 < 1e-23' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g6520" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 108' | awk '$11 < 1e-23' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g6562" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 42' | awk '$11 < 1e-23' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g7278" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 86' | awk '$11 < 1e-23' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g8019" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 90' | awk '$11 < 1e-23' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

fi

done

for file in $(ls effectors/*.faa); do
name=$(basename $file | cut -d '.' -f1)
tblastn -query $file -db db/mildew+db2 -out effectors/${name}_resultsxyz -evalue 1e-5 -outfmt "6 std qseq sseq" -num_threads 12 -max_target_seqs 99999
done

for file in $(ls effectors/*_out13.tsv); do
x=$(basename $file | sed 's@.tsv@.fasta@g')
mv $file effectors/$x
done

for file in $(ls effectors/*_resultsxyz); do
gene=$(basename $file | sed 's@_resultsxyz@@g')
Output=effectors/${gene}_out13.fasta

  if [[ $gene == "g10232" || $gene == "g12450" || $gene == "g12628" || $gene == "g13490" || $gene == "g1355" || $gene == "g13885" || $gene == "g14348" || $gene == "g14832" || $gene == "g15923" || $gene == "g2919" || $gene == "g2928" || $gene == "g3129" || $gene == "g324" || $gene == "g548" || $gene == "g5645" || $gene == "g6450" || $gene == "g757" || $gene == "g9359" ]]; then
    gene_file=$(ls effectors/${gene}.t1.faa)
    gene_length=$(grep -v '>' "$gene_file" | tr -d '\n' | wc -c)
    for file2 in /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta; do
      ID=$(echo "$file2" | cut -d '/' -f7 | cut -c 1-3)_$(echo "$file2" | cut -d '/' -f8 | cut -c 1-3)_$(echo "$file2" | cut -d '/' -f9)
      grep "$ID" "$file" | awk '$3 >= 60' | awk -v gene_length="$gene_length" '$4 / gene_length > 0.7' | awk '$11 < 1e-23' |
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

  elif [[ $gene == "g11054" ]]; then
    for file2 in /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta; do
      ID=$(echo "$file2" | cut -d '/' -f7 | cut -c 1-3)_$(echo "$file2" | cut -d '/' -f8 | cut -c 1-3)_$(echo "$file2" | cut -d '/' -f9)
      grep "$ID" "$file" | awk '$3 >= 60' | awk '$4 > 118' | awk '$11 < 1e-23' |
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g11278" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 203' | awk '$11 < 1e-23' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g1209" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 59' | awk '$11 < 1e-23' |
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g12302" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 88' | awk '$11 < 1e-23' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g12453" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 112' | awk '$11 < 1e-23' |
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g12667" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 32' | awk '$11 < 1e-23' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g12955" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 173' | awk '$11 < 1e-23' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g13367" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 118' | awk '$11 < 1e-23' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g14731" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 102' | awk '$11 < 1e-23' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g2113" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 111' | awk '$11 < 1e-23' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g422" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 155' | awk '$11 < 1e-23' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g5605" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 83' | awk '$11 < 1e-23' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g5721" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 84' | awk '$11 < 1e-23' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g6055" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 197' | awk '$11 < 1e-23' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g6520" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 108' | awk '$11 < 1e-23' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g6562" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 42' | awk '$11 < 1e-23' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g7278" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 86' | awk '$11 < 1e-23' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g8019" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 90' | awk '$11 < 1e-23' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

fi

done



```
```bash
for file in $(ls effectors/*.faa); do
name=$(basename $file | cut -d '.' -f1)
tblastn -query $file -db db/mildew+db2 -out effectors/${name}_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length -num_threads 12 -max_target_seqs 99999
done

for file in $(ls slurm-19869203.out); do
tblastn -query effectors/g10232.t1.faa -db db/mildew+db2 -out effectors/g10232_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g11054.t1.faa -db db/mildew+db2 -out effectors/g11054_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 94 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g11278.t1.faa -db db/mildew+db2 -out effectors/g11278_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 93 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g1209.t1.faa -db db/mildew+db2 -out effectors/g1209_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 304 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g12302.t1.faa -db db/mildew+db2 -out effectors/g12302_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 73 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g12450.t1.faa -db db/mildew+db2 -out effectors/g12450_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g12453.t1.faa -db db/mildew+db2 -out effectors/g12453_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 59 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g12628.t1.faa -db db/mildew+db2 -out effectors/g12628_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g12667.t1.faa -db db/mildew+db2 -out effectors/g12667_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 55 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g12955.t1.faa -db db/mildew+db2 -out effectors/g12955_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 248 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g13367.t1.faa -db db/mildew+db2 -out effectors/g13367_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 100 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g13490.t1.faa -db db/mildew+db2 -out effectors/g13490_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g1355.t1.faa -db db/mildew+db2 -out effectors/g1355_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g13885.t1.faa -db db/mildew+db2 -out effectors/g13885_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g14348.t1.faa -db db/mildew+db2 -out effectors/g14348_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g14731.t1.faa -db db/mildew+db2 -out effectors/g14731_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 59 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g14832.t1.faa -db db/mildew+db2 -out effectors/g14832_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g15923.t1.faa -db db/mildew+db2 -out effectors/g15923_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g2113.t1.faa -db db/mildew+db2 -out effectors/g2113_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 48 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g2919.t1.faa -db db/mildew+db2 -out effectors/g2919_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g2928.t1.faa -db db/mildew+db2 -out effectors/g2928_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g3129.t1.faa -db db/mildew+db2 -out effectors/g3129_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g324.t1.faa -db db/mildew+db2 -out effectors/g324_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g422.t1.faa -db db/mildew+db2 -out effectors/g422_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 349 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g548.t1.faa -db db/mildew+db2 -out effectors/g548_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g5605.t1.faa -db db/mildew+db2 -out effectors/g5605_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 72 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g5645.t1.faa -db db/mildew+db2 -out effectors/g5645_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g5721.t1.faa -db db/mildew+db2 -out effectors/g5721_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 64 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g6055.t1.faa -db db/mildew+db2 -out effectors/g6055_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 67 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g6450.t1.faa -db db/mildew+db2 -out effectors/g6450_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g6520.t1.faa -db db/mildew+db2 -out effectors/g6520_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 73 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g6562.t1.faa -db db/mildew+db2 -out effectors/g6562_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 56 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g7278.t1.faa -db db/mildew+db2 -out effectors/g7278_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 56 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g757.t1.faa -db db/mildew+db2 -out effectors/g757_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g8019.t1.faa -db db/mildew+db2 -out effectors/g8019_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 129 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors/g9359.t1.faa -db db/mildew+db2 -out effectors/g9359_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
done


for file in $(ls effectors/*_resultsxyza); do
gene=$(basename $file | sed 's@_resultsxyza@@g')
Output=effectors/${gene}_out15.fasta

  if [[ $gene == "g10232" || $gene == "g12450" || $gene == "g12628" || $gene == "g13490" || $gene == "g1355" || $gene == "g13885" || $gene == "g14348" || $gene == "g14832" || $gene == "g15923" || $gene == "g2919" || $gene == "g2928" || $gene == "g3129" || $gene == "g324" || $gene == "g548" || $gene == "g5645" || $gene == "g6450" || $gene == "g757" || $gene == "g9359" ]]; then
    gene_file=$(ls effectors/${gene}.t1.faa)
    gene_length=$(grep -v '>' "$gene_file" | tr -d '\n' | wc -c)
    for file2 in /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta; do
      ID=$(echo "$file2" | cut -d '/' -f7 | cut -c 1-3)_$(echo "$file2" | cut -d '/' -f8 | cut -c 1-3)_$(echo "$file2" | cut -d '/' -f9)
      grep "$ID" "$file" | awk '$3 >= 60' | awk -v gene_length="$gene_length" '$4 / gene_length > 0.7' | awk '$11 < 1e-20' |
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

  elif [[ $gene == "g11054" ]]; then
    for file2 in /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta; do
      ID=$(echo "$file2" | cut -d '/' -f7 | cut -c 1-3)_$(echo "$file2" | cut -d '/' -f8 | cut -c 1-3)_$(echo "$file2" | cut -d '/' -f9)
      grep "$ID" "$file" | awk '$3 >= 60' | awk '$4 > 118' | awk '$11 < 1e-20' |
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g11278" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 203' | awk '$11 < 1e-20' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g1209" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 59' | awk '$11 < 1e-20' |
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g12302" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 88' | awk '$11 < 1e-20' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g12453" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 112' | awk '$11 < 1e-20' |
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g12667" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 32' | awk '$11 < 1e-20' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g12955" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 173' | awk '$11 < 1e-20' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g13367" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 118' | awk '$11 < 1e-20' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g14731" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 102' | awk '$11 < 1e-20' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g2113" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 111' | awk '$11 < 1e-20' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g422" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 155' | awk '$11 < 1e-20' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g5605" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 83' | awk '$11 < 1e-20' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g5721" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 84' | awk '$11 < 1e-20' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g6055" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 197' | awk '$11 < 1e-20' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g6520" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 108' | awk '$11 < 1e-20' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g6562" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 42' | awk '$11 < 1e-20' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g7278" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 86' | awk '$11 < 1e-20' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

elif [[ $gene == "g8019" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
grep $ID $file | awk '$3 >= 60' | awk '$4 > 90' | awk '$11 < 1e-20' | 
      while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
      done
    done

fi

done


#No. of hits
Output=effectors/out15.tsv
for file in $(ls effectors/*_resultsxyza); do
gene=$(basename $file | sed 's@_resultsxyza@@g')
echo $gene
echo $gene >> $Output

if [[ $gene == "g10232" || $gene == "g12450" || $gene == "g12628" || $gene == "g13490" || $gene == "g1355" || $gene == "g13885" || $gene == "g14348" || $gene == "g14832" || $gene == "g15923" || $gene == "g2919" || $gene == "g2928" || $gene == "g3129" || $gene == "g324" || $gene == "g548" || $gene == "g5645" || $gene == "g6450" || $gene == "g757" || $gene == "g9359" ]]; then
gene_file=$(ls effectors/${gene}.t1.faa)
gene_length=$(grep -v '>' $gene_file | tr -d '\n' | wc -c)
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk -v gene_length="$gene_length" '$4 / gene_length > 0.7' | awk '$11 < 1e-20' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g11054" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 118' | awk '$11 < 1e-20' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g11278" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 203' | awk '$11 < 1e-20' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g1209" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 59' | awk '$11 < 1e-20' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g12302" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 88' | awk '$11 < 1e-20' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g12453" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 112' | awk '$11 < 1e-20' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g12667" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 32' | awk '$11 < 1e-20' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g12955" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 173' | awk '$11 < 1e-20' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g13367" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 118' | awk '$11 < 1e-20' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g14731" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 102' | awk '$11 < 1e-20' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g2113" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 111' | awk '$11 < 1e-20' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g422" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 155' | awk '$11 < 1e-20' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g5605" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 83' | awk '$11 < 1e-20' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g5721" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 84' | awk '$11 < 1e-20' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g6055" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 197' | awk '$11 < 1e-20' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g6520" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 108' | awk '$11 < 1e-20' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g6562" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 42' | awk '$11 < 1e-20' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g7278" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 86' | awk '$11 < 1e-20' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

elif [[ $gene == "g8019" ]]; then
for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
  ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
x=$(grep $ID $file | awk '$3 >= 60' | awk '$4 > 90' | awk '$11 < 1e-20' | wc -l)
sed -i "\$s/\$/ $x/" $Output
done

fi

done

Output="effectors/out25.tsv"
> $Output 
for file in effectors/*_resultsxyza; do
    gene=$(basename "$file" | sed 's@_resultsxyza@@g')
    echo $gene
    echo -n "$gene" >> $Output  
    while read -r line; do
        x=$(echo "$line" | awk '{print $2}')
        sed -i "\$s/\$/ $x/" $Output
    done < "$file"
    echo "" >> $Output 
done
```
blastp against good gene models
```bash
cat /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta | sed 's@>@>Pod_aph_DRT72020@g' >> good_preds.faa
cat /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta | sed 's@>@>Pod_aph_DRT72021@g' >> good_preds.faa
cat /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta | sed 's@>@>Pod_aph_SCOTT2020@g' >> good_preds.faa
cat /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta | sed 's@>@>Pod_leu_OGBp112020@g' >> good_preds.faa
cat /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta | sed 's@>@>Pod_leu_OGB2019@g' >> good_preds.faa
cat /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta | sed 's@>@>Pod_leu_OGB2021@g' >> good_preds.faa

makeblastdb -in good_preds.faa -dbtype prot -title gooddb2 -out db/gooddb2

for file in $(ls effectors/*.faa); do
name=$(basename $file | cut -d '.' -f1)
blastp -query $file -db db/gooddb2 -out effectors/${name}_resultsxyzab -evalue 1e-5 -outfmt "6 std qseq sseq" -num_threads 12 -max_target_seqs 99999
done

#No. of hits
Output=effectors/out17.tsv
for file in $(ls effectors/*_resultsxyzab); do
gene=$(basename $file | sed 's@_resultsxyzab@@g')
echo $gene
echo $gene >> $Output

  gene_file=$(ls effectors/${gene}.t1.faa)
  gene_length=$(grep -v '>' $gene_file | tr -d '\n' | wc -c)
  for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
    ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
    x=$(grep $ID $file | awk '$3 >= 60' | awk -v gene_length="$gene_length" '$4 / gene_length > 0.7' | awk '$11 < 1e-20' | wc -l)
    sed -i "\$s/\$/ $x/" $Output
  done

done

/home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/DRT72020/fcs/gene_pred/braker/final-02042924/Pod_aph_DRT72020_clean_final_genes.pep.fasta
/home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/DRT72021/fcs/gene_pred/braker/final-02042924/Pod_aph_DRT72021_clean_final_genes.pep.fasta
/home/theaven/scratch/uncompressed/genomes/Podosphaera/aphanis/SCOTT2020/fcs/gene_pred/braker/final-02042924/Pod_aph_SCOTT2020_clean_final_genes.pep.fasta
/home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/GCA_013170925.1/fcs/gene_pred/braker/final-02042924/Pod_leu_GCA_013170925.1_clean_final_genes.pep.fasta
/home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGB2019/fcs/gene_pred/braker/final-02042924/Pod_leu_OGB2019_clean_final_genes.pep.fasta
/home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGB2021/fcs/gene_pred/braker/final-02042924/Pod_leu_OGB2021_clean_final_genes.pep.fasta
/home/theaven/scratch/uncompressed/genomes/Podosphaera/leucotricha/OGBp112020/fcs/gene_pred/braker/final-02042924/Pod_leu_OGBp112020_clean_final_genes.pep.fasta

```
```bash
ls ~/projects/niab/theaven/results/P_aphanis-THeavenDRCT72020_1-aphanis-ranked4.tsv
ls /home/theaven/scratch/uncompressed/gene_pred/P_aphanis/THeavenDRCT72020_1/interproscan/NRI/P.aphanis_strawberry2020_interproscan.tsv

grep 'g10232;\|g11054;\|g11278;\|g1209;\|g12302;\|g12450;\|g12453;\|g12628;\|g12667;\|g12955;\|g13367;\|g13490;\|g1355;\|g13885;\|g14348;\|g14731;\|g14832;\|g15923;\|g2113;\|g2919;\|g2928;\|g3129;\|g324;\|g422;\|g548;\|g5605;\|g5645;\|g5721;\|g6055;\|g6450;\|g6520;\|g6562;\|g7278;\|g757;\|g8019;\|g9359;\|g10232.t1\|g11054.t1\|g11278.t1\|g1209.t1\|g12302.t1\|g12450.t1\|g12453.t1\|g12628.t1\|g12667.t1\|g12955.t1\|g13367.t1\|g13490.t1\|g1355.t1\|g13885.t1\|g14348.t1\|g14731.t1\|g14832.t1\|g15923.t1\|g2113.t1\|g2919.t1\|g2928.t1\|g3129.t1\|g324.t1\|g422.t1\|g548.t1\|g5605.t1\|g5645.t1\|g5721.t1\|g6055.t1\|g6450.t1\|g6520.t1\|g6562.t1\|g7278.t1\|g757.t1\|g8019.t1\|g9359.t1' /home/theaven/scratch/uncompressed/gene_pred/P_aphanis/THeavenDRCT72020_1/interproscan/NRI/P.aphanis_strawberry2020_interproscan.tsv > /home/theaven/scratch/uncompressed/gene_pred/P_aphanis/THeavenDRCT72020_1/interproscan/NRI/aphanis_cseps.tsv
```
leucotricha only effectors
```bash
while IFS=$'\t' read -r line; do
    awk -F'\t' '$31 == "1" && $27 == "." && $7 == "0" {print $1}' <<< "$line" | grep -v 'name' >> leuco_only_secreted_effectors.txt
done < /home/theaven/projects/niab/theaven/results/P_leucotricha-THeavenp11_1-leucotricha-ranked4.tsv

while IFS=$'\t' read -r line; do
echo $line > hits.txt
python3 /home/theaven/scratch/apps/tools/seq_get.py --id_file hits.txt --input /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta --output effectors2/${line}.faa
done < leuco_only_secreted_effectors.txt

grep 'g4176;\|g4600;\|g4897;\|g5141;\|g5142;\|g5372;\|g5674;\|g5692;\|g5745;\|g5831;\|g6026;\|g6027;\|g6105;\|g6162;\|g6243;\|g6281;\|g6315;\|g6354;\|g6355;\|g6406;\|g6458;\|g6590;\|g6999;\|g7051;\|g7120;\|g7121;\|g7132;\|g7198;\|g7200;\|g7361;\|g7603;\|g7646;\|g7647;\|g7696;\|g7802;\|g7927;\|g7928;\|g7935;\|g7937;\|g8009;\|g8035;\|g8081;\|g8248;\|g8278;\|g8734;\|g8760;\|g8897;\|g8974;\|g9008;\|g9065;\|g911;\|g9132;\|g9487;\|g978;\|g981;\|g9883;\|g9979;\|g10095;\|g1010;\|g10145;\|g10372;\|g10489;\|g10497;\|g10500;\|g10675;\|g10692;\|g10783;\|g10784;\|g10794;\|g1106;\|g11286;\|g11369;\|g11409;\|g11416;\|g11565;\|g11652;\|g11660;\|g11719;\|g11727;\|g11778;\|g11806;\|g11807;\|g11808;\|g11809;\|g1186;\|g11906;\|g11998;\|g12135;\|g12284;\|g12296;\|g12452;\|g12489;\|g12500;\|g12653;\|g12716;\|g12776;\|g13051;\|g13093;\|g13170;\|g13217;\|g13219;\|g13226;\|g13566;\|g13811;\|g13892;\|g13921;\|g14015;\|g14026;\|g14153;\|g14154;\|g14161;\|g14208;\|g14239;\|g14266;\|g14430;\|g14456;\|g14463;\|g14471;\|g14546;\|g14551;\|g14690;\|g14900;\|g14908;\|g15093;\|g150;\|g1518;\|g15306;\|g15465;\|g15827;\|g16211;\|g16340;\|g16440;\|g16601;\|g16680;\|g172;\|g1765;\|g17771;\|g18158;\|g1899;\|g1901;\|g2067;\|g2071;\|g2085;\|g2280;\|g2297;\|g254;\|g255;\|g258;\|g259;\|g2761;\|g2902;\|g3072;\|g3170;\|g3172;\|g3204;\|g3237;\|g3319;\|g3921;\|g4016;\|g4176.t1\|g4600.t1\|g4897.t1\|g5141.t1\|g5142.t1\|g5372.t1\|g5674.t1\|g5692.t1\|g5745.t1\|g5831.t1\|g6026.t1\|g6027.t1\|g6105.t1\|g6162.t1\|g6243.t1\|g6281.t1\|g6315.t1\|g6354.t1\|g6355.t1\|g6406.t1\|g6458.t1\|g6590.t1\|g6999.t1\|g7051.t1\|g7120.t1\|g7121.t1\|g7132.t1\|g7198.t1\|g7200.t1\|g7361.t1\|g7603.t1\|g7646.t1\|g7647.t1\|g7696.t1\|g7802.t1\|g7927.t1\|g7928.t1\|g7935.t1\|g7937.t1\|g8009.t1\|g8035.t1\|g8081.t1\|g8248.t1\|g8278.t1\|g8734.t1\|g8760.t1\|g8897.t1\|g8974.t1\|g9008.t1\|g9065.t1\|g911.t1\|g9132.t1\|g9487.t1\|g978.t1\|g981.t1\|g9883.t1\|g9979.t1\|g10095.t1\|g1010.t1\|g10145.t1\|g10372.t1\|g10489.t1\|g10497.t1\|g10500.t1\|g10675.t1\|g10692.t1\|g10783.t1\|g10784.t1\|g10794.t1\|g1106.t1\|g11286.t1\|g11369.t1\|g11409.t1\|g11416.t1\|g11565.t1\|g11652.t1\|g11660.t1\|g11719.t1\|g11727.t1\|g11778.t1\|g11806.t1\|g11807.t1\|g11808.t1\|g11809.t1\|g1186.t1\|g11906.t1\|g11998.t1\|g12135.t1\|g12284.t1\|g12296.t1\|g12452.t1\|g12489.t1\|g12500.t1\|g12653.t1\|g12716.t1\|g12776.t1\|g13051.t1\|g13093.t1\|g13170.t1\|g13217.t1\|g13219.t1\|g13226.t1\|g13566.t1\|g13811.t1\|g13892.t1\|g13921.t1\|g14015.t1\|g14026.t1\|g14153.t1\|g14154.t1\|g14161.t1\|g14208.t1\|g14239.t1\|g14266.t1\|g14430.t1\|g14456.t1\|g14463.t1\|g14471.t1\|g14546.t1\|g14551.t1\|g14690.t1\|g14900.t1\|g14908.t1\|g15093.t1\|g150.t1\|g1518.t1\|g15306.t1\|g15465.t1\|g15827.t1\|g16211.t1\|g16340.t1\|g16440.t1\|g16601.t1\|g16680.t1\|g172.t1\|g1765.t1\|g17771.t1\|g18158.t1\|g1899.t1\|g1901.t1\|g2067.t1\|g2071.t1\|g2085.t1\|g2280.t1\|g2297.t1\|g254.t1\|g255.t1\|g258.t1\|g259.t1\|g2761.t1\|g2902.t1\|g3072.t1\|g3170.t1\|g3172.t1\|g3204.t1\|g3237.t1\|g3319.t1\|g3921.t1\|g4016.t1' /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3 > effectors2/leuco_only.gff

#max intron length 1.1x the largest intron in the gene prediction or 40, whichever is larger
for file in $(ls slurm-19869203.out); do
tblastn -query effectors2/g10095.t1.faa -db db/mildew+db2 -out effectors2/g10095_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g1010.t1.faa -db db/mildew+db2 -out effectors2/g1010_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 49 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g10145.t1.faa -db db/mildew+db2 -out effectors2/g10145_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g10372.t1.faa -db db/mildew+db2 -out effectors2/g10372_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g10489.t1.faa -db db/mildew+db2 -out effectors2/g10489_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 81 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g10497.t1.faa -db db/mildew+db2 -out effectors2/g10497_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 73 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g10500.t1.faa -db db/mildew+db2 -out effectors2/g10500_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 132 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g10675.t1.faa -db db/mildew+db2 -out effectors2/g10675_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g10692.t1.faa -db db/mildew+db2 -out effectors2/g10692_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g10783.t1.faa -db db/mildew+db2 -out effectors2/g10783_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 104 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g10784.t1.faa -db db/mildew+db2 -out effectors2/g10784_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 95 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g10794.t1.faa -db db/mildew+db2 -out effectors2/g10794_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 65 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g1106.t1.faa -db db/mildew+db2 -out effectors2/g1106_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 76 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g11286.t1.faa -db db/mildew+db2 -out effectors2/g11286_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 72 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g11369.t1.faa -db db/mildew+db2 -out effectors2/g11369_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 72 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g11409.t1.faa -db db/mildew+db2 -out effectors2/g11409_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g11416.t1.faa -db db/mildew+db2 -out effectors2/g11416_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 91 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g11565.t1.faa -db db/mildew+db2 -out effectors2/g11565_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g11652.t1.faa -db db/mildew+db2 -out effectors2/g11652_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g11660.t1.faa -db db/mildew+db2 -out effectors2/g11660_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 80 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g11719.t1.faa -db db/mildew+db2 -out effectors2/g11719_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 85 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g11727.t1.faa -db db/mildew+db2 -out effectors2/g11727_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 88 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g11778.t1.faa -db db/mildew+db2 -out effectors2/g11778_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g11806.t1.faa -db db/mildew+db2 -out effectors2/g11806_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 59 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g11807.t1.faa -db db/mildew+db2 -out effectors2/g11807_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 59 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g11808.t1.faa -db db/mildew+db2 -out effectors2/g11808_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 66 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g11809.t1.faa -db db/mildew+db2 -out effectors2/g11809_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 65 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g1186.t1.faa -db db/mildew+db2 -out effectors2/g1186_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 54 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g11906.t1.faa -db db/mildew+db2 -out effectors2/g11906_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 82 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g11998.t1.faa -db db/mildew+db2 -out effectors2/g11998_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g12135.t1.faa -db db/mildew+db2 -out effectors2/g12135_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 189 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g12284.t1.faa -db db/mildew+db2 -out effectors2/g12284_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 63 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g12296.t1.faa -db db/mildew+db2 -out effectors2/g12296_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 83 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g12452.t1.faa -db db/mildew+db2 -out effectors2/g12452_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g12489.t1.faa -db db/mildew+db2 -out effectors2/g12489_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g12500.t1.faa -db db/mildew+db2 -out effectors2/g12500_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 55 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g12653.t1.faa -db db/mildew+db2 -out effectors2/g12653_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 198 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g12716.t1.faa -db db/mildew+db2 -out effectors2/g12716_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 62 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g12776.t1.faa -db db/mildew+db2 -out effectors2/g12776_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 110 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g13051.t1.faa -db db/mildew+db2 -out effectors2/g13051_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 65 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g13093.t1.faa -db db/mildew+db2 -out effectors2/g13093_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g13170.t1.faa -db db/mildew+db2 -out effectors2/g13170_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 61 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g13217.t1.faa -db db/mildew+db2 -out effectors2/g13217_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 62 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g13219.t1.faa -db db/mildew+db2 -out effectors2/g13219_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 61 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g13226.t1.faa -db db/mildew+db2 -out effectors2/g13226_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g13566.t1.faa -db db/mildew+db2 -out effectors2/g13566_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 66 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g13811.t1.faa -db db/mildew+db2 -out effectors2/g13811_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 105 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g13892.t1.faa -db db/mildew+db2 -out effectors2/g13892_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g13921.t1.faa -db db/mildew+db2 -out effectors2/g13921_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 169 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g14015.t1.faa -db db/mildew+db2 -out effectors2/g14015_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 106 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g14026.t1.faa -db db/mildew+db2 -out effectors2/g14026_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 71 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g14153.t1.faa -db db/mildew+db2 -out effectors2/g14153_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 104 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g14154.t1.faa -db db/mildew+db2 -out effectors2/g14154_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 82 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g14161.t1.faa -db db/mildew+db2 -out effectors2/g14161_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g14208.t1.faa -db db/mildew+db2 -out effectors2/g14208_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 65 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g14239.t1.faa -db db/mildew+db2 -out effectors2/g14239_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g14266.t1.faa -db db/mildew+db2 -out effectors2/g14266_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 70 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g14430.t1.faa -db db/mildew+db2 -out effectors2/g14430_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 66 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g14456.t1.faa -db db/mildew+db2 -out effectors2/g14456_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 59 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g14463.t1.faa -db db/mildew+db2 -out effectors2/g14463_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 70 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g14471.t1.faa -db db/mildew+db2 -out effectors2/g14471_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 70 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g14546.t1.faa -db db/mildew+db2 -out effectors2/g14546_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 74 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g14551.t1.faa -db db/mildew+db2 -out effectors2/g14551_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 58 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g14690.t1.faa -db db/mildew+db2 -out effectors2/g14690_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 70 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g14900.t1.faa -db db/mildew+db2 -out effectors2/g14900_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 69 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g14908.t1.faa -db db/mildew+db2 -out effectors2/g14908_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 62 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g150.t1.faa -db db/mildew+db2 -out effectors2/g150_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g15093.t1.faa -db db/mildew+db2 -out effectors2/g15093_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g1518.t1.faa -db db/mildew+db2 -out effectors2/g1518_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g15306.t1.faa -db db/mildew+db2 -out effectors2/g15306_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 104 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g15465.t1.faa -db db/mildew+db2 -out effectors2/g15465_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 75 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g15827.t1.faa -db db/mildew+db2 -out effectors2/g15827_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g16211.t1.faa -db db/mildew+db2 -out effectors2/g16211_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 68 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g16340.t1.faa -db db/mildew+db2 -out effectors2/g16340_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g16440.t1.faa -db db/mildew+db2 -out effectors2/g16440_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 102 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g16601.t1.faa -db db/mildew+db2 -out effectors2/g16601_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 69 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g16680.t1.faa -db db/mildew+db2 -out effectors2/g16680_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g172.t1.faa -db db/mildew+db2 -out effectors2/g172_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g1765.t1.faa -db db/mildew+db2 -out effectors2/g1765_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 101 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g17771.t1.faa -db db/mildew+db2 -out effectors2/g17771_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g18158.t1.faa -db db/mildew+db2 -out effectors2/g18158_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g1899.t1.faa -db db/mildew+db2 -out effectors2/g1899_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g1901.t1.faa -db db/mildew+db2 -out effectors2/g1901_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 54 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g2067.t1.faa -db db/mildew+db2 -out effectors2/g2067_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 98 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g2071.t1.faa -db db/mildew+db2 -out effectors2/g2071_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g2085.t1.faa -db db/mildew+db2 -out effectors2/g2085_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 97 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g2280.t1.faa -db db/mildew+db2 -out effectors2/g2280_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 70 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g2297.t1.faa -db db/mildew+db2 -out effectors2/g2297_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 61 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g254.t1.faa -db db/mildew+db2 -out effectors2/g254_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 52 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g255.t1.faa -db db/mildew+db2 -out effectors2/g255_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 109 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g258.t1.faa -db db/mildew+db2 -out effectors2/g258_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 86 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g259.t1.faa -db db/mildew+db2 -out effectors2/g259_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 95 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g2761.t1.faa -db db/mildew+db2 -out effectors2/g2761_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g2902.t1.faa -db db/mildew+db2 -out effectors2/g2902_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g3072.t1.faa -db db/mildew+db2 -out effectors2/g3072_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 66 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g3170.t1.faa -db db/mildew+db2 -out effectors2/g3170_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 51 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g3172.t1.faa -db db/mildew+db2 -out effectors2/g3172_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 51 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g3204.t1.faa -db db/mildew+db2 -out effectors2/g3204_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 104 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g3237.t1.faa -db db/mildew+db2 -out effectors2/g3237_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 55 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g3319.t1.faa -db db/mildew+db2 -out effectors2/g3319_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 91 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g3921.t1.faa -db db/mildew+db2 -out effectors2/g3921_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 57 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g4016.t1.faa -db db/mildew+db2 -out effectors2/g4016_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g4176.t1.faa -db db/mildew+db2 -out effectors2/g4176_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 94 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g4600.t1.faa -db db/mildew+db2 -out effectors2/g4600_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 168 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g4897.t1.faa -db db/mildew+db2 -out effectors2/g4897_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 73 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g5141.t1.faa -db db/mildew+db2 -out effectors2/g5141_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 68 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g5142.t1.faa -db db/mildew+db2 -out effectors2/g5142_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 51 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g5372.t1.faa -db db/mildew+db2 -out effectors2/g5372_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 87 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g5674.t1.faa -db db/mildew+db2 -out effectors2/g5674_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 75 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g5692.t1.faa -db db/mildew+db2 -out effectors2/g5692_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g5745.t1.faa -db db/mildew+db2 -out effectors2/g5745_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 49 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g5831.t1.faa -db db/mildew+db2 -out effectors2/g5831_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 75 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g6026.t1.faa -db db/mildew+db2 -out effectors2/g6026_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 65 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g6027.t1.faa -db db/mildew+db2 -out effectors2/g6027_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 59 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g6105.t1.faa -db db/mildew+db2 -out effectors2/g6105_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 51 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g6162.t1.faa -db db/mildew+db2 -out effectors2/g6162_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 107 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g6243.t1.faa -db db/mildew+db2 -out effectors2/g6243_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 91 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g6281.t1.faa -db db/mildew+db2 -out effectors2/g6281_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 57 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g6315.t1.faa -db db/mildew+db2 -out effectors2/g6315_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 68 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g6354.t1.faa -db db/mildew+db2 -out effectors2/g6354_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 106 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g6355.t1.faa -db db/mildew+db2 -out effectors2/g6355_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 96 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g6406.t1.faa -db db/mildew+db2 -out effectors2/g6406_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g6458.t1.faa -db db/mildew+db2 -out effectors2/g6458_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g6590.t1.faa -db db/mildew+db2 -out effectors2/g6590_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 58 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g6999.t1.faa -db db/mildew+db2 -out effectors2/g6999_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 58 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g7051.t1.faa -db db/mildew+db2 -out effectors2/g7051_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 87 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g7120.t1.faa -db db/mildew+db2 -out effectors2/g7120_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 71 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g7121.t1.faa -db db/mildew+db2 -out effectors2/g7121_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 57 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g7132.t1.faa -db db/mildew+db2 -out effectors2/g7132_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g7198.t1.faa -db db/mildew+db2 -out effectors2/g7198_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g7200.t1.faa -db db/mildew+db2 -out effectors2/g7200_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 60 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g7361.t1.faa -db db/mildew+db2 -out effectors2/g7361_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 103 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g7603.t1.faa -db db/mildew+db2 -out effectors2/g7603_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 97 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g7646.t1.faa -db db/mildew+db2 -out effectors2/g7646_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 63 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g7647.t1.faa -db db/mildew+db2 -out effectors2/g7647_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 60 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g7696.t1.faa -db db/mildew+db2 -out effectors2/g7696_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g7802.t1.faa -db db/mildew+db2 -out effectors2/g7802_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 90 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g7927.t1.faa -db db/mildew+db2 -out effectors2/g7927_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 116 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g7928.t1.faa -db db/mildew+db2 -out effectors2/g7928_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 95 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g7935.t1.faa -db db/mildew+db2 -out effectors2/g7935_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 70 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g7937.t1.faa -db db/mildew+db2 -out effectors2/g7937_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 70 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g8009.t1.faa -db db/mildew+db2 -out effectors2/g8009_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 76 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g8035.t1.faa -db db/mildew+db2 -out effectors2/g8035_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g8081.t1.faa -db db/mildew+db2 -out effectors2/g8081_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 57 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g8248.t1.faa -db db/mildew+db2 -out effectors2/g8248_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g8278.t1.faa -db db/mildew+db2 -out effectors2/g8278_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g8734.t1.faa -db db/mildew+db2 -out effectors2/g8734_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 91 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g8760.t1.faa -db db/mildew+db2 -out effectors2/g8760_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 91 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g8897.t1.faa -db db/mildew+db2 -out effectors2/g8897_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 62 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g8974.t1.faa -db db/mildew+db2 -out effectors2/g8974_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g9008.t1.faa -db db/mildew+db2 -out effectors2/g9008_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 69 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g9065.t1.faa -db db/mildew+db2 -out effectors2/g9065_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 73 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g911.t1.faa -db db/mildew+db2 -out effectors2/g911_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 75 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g9132.t1.faa -db db/mildew+db2 -out effectors2/g9132_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 82 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g9487.t1.faa -db db/mildew+db2 -out effectors2/g9487_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 53 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g978.t1.faa -db db/mildew+db2 -out effectors2/g978_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g981.t1.faa -db db/mildew+db2 -out effectors2/g981_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 95 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g9883.t1.faa -db db/mildew+db2 -out effectors2/g9883_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 98 -num_threads 12 -max_target_seqs 99999
tblastn -query effectors2/g9979.t1.faa -db db/mildew+db2 -out effectors2/g9979_resultsxyza -evalue 1e-5 -outfmt "6 std qseq sseq" -max_intron_length 40 -num_threads 12 -max_target_seqs 99999
done

tblastn -query effectors2/g8974.t1.faa -db db/mildew+db2 -out effectors2/g8974_resultsxyza_2 -evalue 1e-5 -outfmt "6 std qseq sseq" -num_threads 12 -max_target_seqs 99999
#this returns no hits, possibly because the gene is so short?

#length of hit 70% of lingest exon:
declare -A exon_lengths
exon_lengths["g1010"]=371
exon_lengths["g10497"]=221
exon_lengths["g10500"]=208
exon_lengths["g10783"]=242
exon_lengths["g10784"]=244
exon_lengths["g10794"]=233
exon_lengths["g11286"]=454
exon_lengths["g11369"]=485
exon_lengths["g11719"]=310
exon_lengths["g11807"]=223
exon_lengths["g11808"]=242
exon_lengths["g11809"]=217
exon_lengths["g1186"]=165
exon_lengths["g12284"]=214
exon_lengths["g12500"]=212
exon_lengths["g12653"]=445
exon_lengths["g12716"]=896
exon_lengths["g12776"]=202
exon_lengths["g13051"]=240
exon_lengths["g13170"]=206
exon_lengths["g13217"]=99
exon_lengths["g13219"]=212
exon_lengths["g13566"]=240
exon_lengths["g14026"]=210
exon_lengths["g14208"]=204
exon_lengths["g14266"]=292
exon_lengths["g14430"]=223
exon_lengths["g14463"]=196
exon_lengths["g14471"]=235
exon_lengths["g14690"]=558
exon_lengths["g14900"]=225
exon_lengths["g14908"]=227
exon_lengths["g15465"]=242
exon_lengths["g16211"]=661
exon_lengths["g16440"]=882
exon_lengths["g16601"]=135
exon_lengths["g1901"]=214
exon_lengths["g2280"]=79
exon_lengths["g2297"]=263
exon_lengths["g254"]=225
exon_lengths["g3072"]=248
exon_lengths["g3170"]=233
exon_lengths["g3172"]=269
exon_lengths["g3204"]=206
exon_lengths["g3237"]=208
exon_lengths["g3921"]=169
exon_lengths["g5142"]=217
exon_lengths["g5674"]=246
exon_lengths["g5745"]=231
exon_lengths["g5831"]=214
exon_lengths["g6026"]=227
exon_lengths["g6027"]=227
exon_lengths["g6105"]=408
exon_lengths["g6281"]=206
exon_lengths["g6590"]=212
exon_lengths["g6999"]=263
exon_lengths["g7051"]=214
exon_lengths["g7120"]=246
exon_lengths["g7121"]=294
exon_lengths["g7200"]=248
exon_lengths["g7646"]=248
exon_lengths["g7647"]=231
exon_lengths["g7935"]=244
exon_lengths["g7937"]=238
exon_lengths["g8081"]=217
exon_lengths["g8760"]=240
exon_lengths["g8897"]=128
exon_lengths["g9008"]=210
exon_lengths["g9065"]=395
exon_lengths["g911"]=235
exon_lengths["g9487"]=156
exon_lengths["g10489"]=231
exon_lengths["g1106"]=316
exon_lengths["g11660"]=245
exon_lengths["g11806"]=252
exon_lengths["g13921"]=210
exon_lengths["g14154"]=475
exon_lengths["g14456"]=114
exon_lengths["g14551"]=520
exon_lengths["g4897"]=235
exon_lengths["g6315"]=293
exon_lengths["g8009"]=156
exon_lengths["g8734"]=80
exon_lengths["g11416"]=270
exon_lengths["g11727"]=283
exon_lengths["g11906"]=266
exon_lengths["g12135"]=270
exon_lengths["g12296"]=347
exon_lengths["g13811"]=264
exon_lengths["g14153"]=258
exon_lengths["g14546"]=275
exon_lengths["g15306"]=266
exon_lengths["g1765"]=270
exon_lengths["g2067"]=264
exon_lengths["g2085"]=258
exon_lengths["g255"]=283
exon_lengths["g258"]=266
exon_lengths["g259"]=275
exon_lengths["g3319"]=258
exon_lengths["g4176"]=264
exon_lengths["g5141"]=144
exon_lengths["g5372"]=266
exon_lengths["g6162"]=277
exon_lengths["g6243"]=170
exon_lengths["g6354"]=291
exon_lengths["g6355"]=289
exon_lengths["g7361"]=270
exon_lengths["g7603"]=270
exon_lengths["g7802"]=273
exon_lengths["g7927"]=277
exon_lengths["g7928"]=277
exon_lengths["g9132"]=266
exon_lengths["g981"]=258
exon_lengths["g9883"]=270
exon_lengths["g14015"]=270
exon_lengths["g4600"]=124

#get fastas
for file in $(ls effectors2/*_resultsxyza); do
gene=$(basename $file | sed 's@_resultsxyza@@g')
Output=effectors2/${gene}_out5.fasta
echo $gene

if [[ $gene == "g10675" || $gene == "g10692"  || $gene == "g10095" || $gene == "g10145" || $gene == "g10372" || $gene == "g11409" || $gene == "g11565" || $gene == "g11652" || $gene == "g11778" || $gene == "g11998" || $gene == "g12452" || $gene == "g12489" || $gene == "g13093" || $gene == "g13226" || $gene == "g13892" || $gene == "g14161" || $gene == "g14239" || $gene == "g150" || $gene == "g15093" || $gene == "g1518" || $gene == "g15827" || $gene == "g16340" || $gene == "g16680" || $gene == "g172" || $gene == "g17771" || $gene == "g18158" || $gene == "g1899" || $gene == "g2071" || $gene == "g2761" || $gene == "g2902" || $gene == "g4016" || $gene == "g5692" || $gene == "g6406" || $gene == "g6458" || $gene == "g7132" || $gene == "g7198" || $gene == "g7696" || $gene == "g8035" || $gene == "g8248" || $gene == "g8278" || $gene == "g8974" || $gene == "g978" || $gene == "g9979" ]]; then
  gene_file=$(ls effectors2/${gene}.t1.faa)
  gene_length=$(grep -v '>' "$gene_file" | tr -d '\n' | wc -c)
  for file2 in /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta; do
    ID=$(echo "$file2" | cut -d '/' -f7 | cut -c 1-3)_$(echo "$file2" | cut -d '/' -f8 | cut -c 1-3)_$(echo "$file2" | cut -d '/' -f9)
    grep "$ID" "$file" | awk '$3 >= 60' | awk -v gene_length="$gene_length" '$4 / gene_length > 0.7' | awk '$11 < 1e-20' |
    while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
    done
  done
else 
  exon_length=${exon_lengths[$gene]}
  for file2 in /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta; do
    ID=$(echo "$file2" | cut -d '/' -f7 | cut -c 1-3)_$(echo "$file2" | cut -d '/' -f8 | cut -c 1-3)_$(echo "$file2" | cut -d '/' -f9)
    grep "$ID" "$file" | awk '$3 >= 60' | awk -v exon_length="$exon_length" '$4 > exon_length / 3' | awk '$11 < 1e-20' |
    while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12 col13 col14; do
        echo ">${col2}_${col9}" >> "$Output"
        echo "$col14" >> "$Output"
    done
  done
fi
done

#No. of hits
Output=effectors2/out5.tsv
for file in $(ls effectors2/*_resultsxyza); do
gene=$(basename $file | sed 's@_resultsxyza@@g')
echo $gene
echo $gene >> $Output

if [[ $gene == "g10675" || $gene == "g10692"  || $gene == "g10095" || $gene == "g10145" || $gene == "g10372" || $gene == "g11409" || $gene == "g11565" || $gene == "g11652" || $gene == "g11778" || $gene == "g11998" || $gene == "g12452" || $gene == "g12489" || $gene == "g13093" || $gene == "g13226" || $gene == "g13892" || $gene == "g14161" || $gene == "g14239" || $gene == "g150" || $gene == "g15093" || $gene == "g1518" || $gene == "g15827" || $gene == "g16340" || $gene == "g16680" || $gene == "g172" || $gene == "g17771" || $gene == "g18158" || $gene == "g1899" || $gene == "g2071" || $gene == "g2761" || $gene == "g2902" || $gene == "g4016" || $gene == "g5692" || $gene == "g6406" || $gene == "g6458" || $gene == "g7132" || $gene == "g7198" || $gene == "g7696" || $gene == "g8035" || $gene == "g8248" || $gene == "g8278" || $gene == "g8974" || $gene == "g978" || $gene == "g9979" ]]; then
  gene_file=$(ls effectors2/${gene}.t1.faa)
  gene_length=$(grep -v '>' $gene_file | tr -d '\n' | wc -c)
  for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
    ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
    x=$(grep $ID $file | awk '$3 >= 60' | awk -v gene_length="$gene_length" '$4 / gene_length > 0.7' | awk '$11 < 1e-20' | wc -l)
    sed -i "\$s/\$/ $x/" $Output
  done
else 
  exon_length=${exon_lengths[$gene]}
  for file2 in /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta; do
    ID=$(echo "$file2" | cut -d '/' -f7 | cut -c 1-3)_$(echo "$file2" | cut -d '/' -f8 | cut -c 1-3)_$(echo "$file2" | cut -d '/' -f9)
    x=$(grep "$ID" "$file" | awk '$3 >= 60' | awk -v exon_length="$exon_length" '$4 > exon_length / 3' | awk '$11 < 1e-20' | wc -l)
    sed -i "\$s/\$/ $x/" $Output
  done
fi
done


```
blastp against good gene models
```bash
for file in $(ls effectors2/*.faa); do
name=$(basename $file | cut -d '.' -f1)
blastp -query $file -db db/gooddb2 -out effectors2/${name}_resultsxyzab -evalue 1e-5 -outfmt "6 std qseq sseq" -num_threads 12 -max_target_seqs 99999
done

#No. of hits
Output=effectors2/out6.tsv
for file in $(ls effectors2/*_resultsxyzab); do
gene=$(basename $file | sed 's@_resultsxyzab@@g')
echo $gene
echo $gene >> $Output
  gene_file=$(ls effectors2/${gene}.t1.faa)
  gene_length=$(grep -v '>' $gene_file | tr -d '\n' | wc -c)
  for file2 in $(ls /home/theaven/scratch/uncompressed/genomes/*/*/*/fcs/gene_pred/braker/final-02042924/*clean_final_genes.pep.fasta); do
    ID=$(echo $file2 | cut -d '/' -f7 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f8 | cut -c 1-3)_$(echo $file2 | cut -d '/' -f9)
    x=$(grep $ID $file | awk '$3 >= 60' | awk -v gene_length="$gene_length" '$4 / gene_length > 0.7' | awk '$11 < 1e-20' | wc -l)
    sed -i "\$s/\$/ $x/" $Output
  done
done
```

```bash
ls /home/theaven/scratch/uncompressed/gene_pred/P_leucotricha/THeavenp11_1/interproscan/NRI/P.leucotricha_apple2020_interproscan.tsv

grep 'g4176;\|g4600;\|g4897;\|g5141;\|g5142;\|g5372;\|g5674;\|g5692;\|g5745;\|g5831;\|g6026;\|g6027;\|g6105;\|g6162;\|g6243;\|g6281;\|g6315;\|g6354;\|g6355;\|g6406;\|g6458;\|g6590;\|g6999;\|g7051;\|g7120;\|g7121;\|g7132;\|g7198;\|g7200;\|g7361;\|g7603;\|g7646;\|g7647;\|g7696;\|g7802;\|g7927;\|g7928;\|g7935;\|g7937;\|g8009;\|g8035;\|g8081;\|g8248;\|g8278;\|g8734;\|g8760;\|g8897;\|g8974;\|g9008;\|g9065;\|g911;\|g9132;\|g9487;\|g978;\|g981;\|g9883;\|g9979;\|g10095;\|g1010;\|g10145;\|g10372;\|g10489;\|g10497;\|g10500;\|g10675;\|g10692;\|g10783;\|g10784;\|g10794;\|g1106;\|g11286;\|g11369;\|g11409;\|g11416;\|g11565;\|g11652;\|g11660;\|g11719;\|g11727;\|g11778;\|g11806;\|g11807;\|g11808;\|g11809;\|g1186;\|g11906;\|g11998;\|g12135;\|g12284;\|g12296;\|g12452;\|g12489;\|g12500;\|g12653;\|g12716;\|g12776;\|g13051;\|g13093;\|g13170;\|g13217;\|g13219;\|g13226;\|g13566;\|g13811;\|g13892;\|g13921;\|g14015;\|g14026;\|g14153;\|g14154;\|g14161;\|g14208;\|g14239;\|g14266;\|g14430;\|g14456;\|g14463;\|g14471;\|g14546;\|g14551;\|g14690;\|g14900;\|g14908;\|g15093;\|g150;\|g1518;\|g15306;\|g15465;\|g15827;\|g16211;\|g16340;\|g16440;\|g16601;\|g16680;\|g172;\|g1765;\|g17771;\|g18158;\|g1899;\|g1901;\|g2067;\|g2071;\|g2085;\|g2280;\|g2297;\|g254;\|g255;\|g258;\|g259;\|g2761;\|g2902;\|g3072;\|g3170;\|g3172;\|g3204;\|g3237;\|g3319;\|g3921;\|g4016;\|g4176.t1\|g4600.t1\|g4897.t1\|g5141.t1\|g5142.t1\|g5372.t1\|g5674.t1\|g5692.t1\|g5745.t1\|g5831.t1\|g6026.t1\|g6027.t1\|g6105.t1\|g6162.t1\|g6243.t1\|g6281.t1\|g6315.t1\|g6354.t1\|g6355.t1\|g6406.t1\|g6458.t1\|g6590.t1\|g6999.t1\|g7051.t1\|g7120.t1\|g7121.t1\|g7132.t1\|g7198.t1\|g7200.t1\|g7361.t1\|g7603.t1\|g7646.t1\|g7647.t1\|g7696.t1\|g7802.t1\|g7927.t1\|g7928.t1\|g7935.t1\|g7937.t1\|g8009.t1\|g8035.t1\|g8081.t1\|g8248.t1\|g8278.t1\|g8734.t1\|g8760.t1\|g8897.t1\|g8974.t1\|g9008.t1\|g9065.t1\|g911.t1\|g9132.t1\|g9487.t1\|g978.t1\|g981.t1\|g9883.t1\|g9979.t1\|g10095.t1\|g1010.t1\|g10145.t1\|g10372.t1\|g10489.t1\|g10497.t1\|g10500.t1\|g10675.t1\|g10692.t1\|g10783.t1\|g10784.t1\|g10794.t1\|g1106.t1\|g11286.t1\|g11369.t1\|g11409.t1\|g11416.t1\|g11565.t1\|g11652.t1\|g11660.t1\|g11719.t1\|g11727.t1\|g11778.t1\|g11806.t1\|g11807.t1\|g11808.t1\|g11809.t1\|g1186.t1\|g11906.t1\|g11998.t1\|g12135.t1\|g12284.t1\|g12296.t1\|g12452.t1\|g12489.t1\|g12500.t1\|g12653.t1\|g12716.t1\|g12776.t1\|g13051.t1\|g13093.t1\|g13170.t1\|g13217.t1\|g13219.t1\|g13226.t1\|g13566.t1\|g13811.t1\|g13892.t1\|g13921.t1\|g14015.t1\|g14026.t1\|g14153.t1\|g14154.t1\|g14161.t1\|g14208.t1\|g14239.t1\|g14266.t1\|g14430.t1\|g14456.t1\|g14463.t1\|g14471.t1\|g14546.t1\|g14551.t1\|g14690.t1\|g14900.t1\|g14908.t1\|g15093.t1\|g150.t1\|g1518.t1\|g15306.t1\|g15465.t1\|g15827.t1\|g16211.t1\|g16340.t1\|g16440.t1\|g16601.t1\|g16680.t1\|g172.t1\|g1765.t1\|g17771.t1\|g18158.t1\|g1899.t1\|g1901.t1\|g2067.t1\|g2071.t1\|g2085.t1\|g2280.t1\|g2297.t1\|g254.t1\|g255.t1\|g258.t1\|g259.t1\|g2761.t1\|g2902.t1\|g3072.t1\|g3170.t1\|g3172.t1\|g3204.t1\|g3237.t1\|g3319.t1\|g3921.t1\|g4016.t1' /home/theaven/scratch/uncompressed/gene_pred/P_leucotricha/THeavenp11_1/interproscan/NRI/P.leucotricha_apple2020_interproscan.tsv > /home/theaven/scratch/uncompressed/gene_pred/P_leucotricha/THeavenp11_1/interproscan/NRI/leucotricha_cseps.tsv
```

#### Alphafold
```bash
cat effectors/*.faa > aphanis-unique-cseps.faa
cat effectors2/*.faa > leucotricha-unique-cseps.faa

/home/theaven/scratch/apps/tools/signalp-4.1/signalp -t euk -f short aphanis-unique-cseps.faa > aphanis-unique-cseps-signalp4.1.out
/home/theaven/scratch/apps/tools/signalp-4.1/signalp -t euk -f short leucotricha-unique-cseps.faa > leucotricha-unique-cseps-signalp4.1.out

#mot all proteins have signal peptide as predictd by signalp4.1, some are only predicted by signalp3, phobius or deepsig. Where there was dissagreement between different prediction methods signalp4.1 was considered most reliable, followed by singalp3NN, then signalp3HMM, then phobius, then deepsig - however inmost cases the different tools were in agreement with one another.
```
Aphanis_Gene: Signal predicted:
g10232  SignalP3NN:23
g11054  SignalP3NN:22,SignalP3HMM:22,SignalP4:22,SignalP5:22,DeepSig:22,Phobius:22,SignalP6:22
g11278  SignalP3NN:22,SignalP3HMM:22,SignalP4:22,SignalP5:22,DeepSig:22,Phobius:22,SignalP6:22
g1209 SignalP3NN:21
g12302  SignalP3NN:20,SignalP3HMM:20,SignalP4:20,SignalP5:20,DeepSig:20,Phobius:20,SignalP6:20
g12450  SignalP3NN:19,Phobius:19
g12453  Phobius:20
g12628  SignalP3NN:17,SignalP4:17,DeepSig:17
g12667* SignalP3HMM:19
g12955  Phobius:17
g13367  SignalP3NN:22,SignalP3HMM:22,SignalP4:22,SignalP5:22,DeepSig:22,Phobius:22,SignalP6:22
g13490  SignalP3NN:22
g1355 SignalP3NN:17,SignalP3HMM:17,SignalP4:17,Phobius:17
g13885  Phobius:23
g14348  Phobius:19,SignalP3NN:24,SignalP3HMM:24,SignalP4:24,SignalP5:24,DeepSig:24
g14731  SignalP3NN:14
g14832  SignalP3NN:19,Phobius:19
g15923  SignalP3NN:27,SignalP3HMM:27,Phobius:27
g2113 Phobius:18
g2919 Phobius:29
g2928 Phobius:16
g3129 SignalP3NN:21,DeepSig:21
g324  SignalP3NN:27
g422  SignalP3NN:30
g548  SignalP3NN:20
g5605 SignalP3NN:20,SignalP3HMM:20,SignalP4:20,SignalP5:20,DeepSig:20,Phobius:20,SignalP6:20
g5645 SignalP3HMM:16,Phobius:16
g5721 SignalP3NN:28,SignalP3HMM:28,SignalP4:28,SignalP5:28,Phobius:28,SignalP6:28
g6055 SignalP3NN:24,DeepSig:24,Phobius:24
g6450 SignalP3NN:28
g6520*  SignalP3NN:19
g6562 SignalP3HMM:18,SignalP5:18
g7278 SignalP3HMM:19,SignalP3NN:23
g757  SignalP4:10
g8019 SignalP3NN:20,SignalP3HMM:20,SignalP4:20,SignalP5:20,DeepSig:20,Phobius:20,SignalP6:30
g9359 SignalP3HMM:31,Phobius:31
  
Leuco_Gene: 
g10095  SignalP3HMM:24
g1010 SignalP3HMM:18,Phobius:18
g10145  SignalP3NN:23,SignalP4:23,Phobius:23
g10372  SignalP3NN:20,SignalP4:20,SignalP5:20,DeepSig:20,Phobius:20,SignalP6:20
g10489  SignalP3NN:15,SignalP3HMM:15,Phobius:15,SignalP4:19,SignalP6:19,SignalP5:23,DeepSig:23
g10497  SignalP3NN:19,SignalP3HMM:19,SignalP4:19,SignalP5:19,DeepSig:19,Phobius:19,SignalP6:19
g10500  SignalP3NN:21,SignalP3HMM:21,SignalP4:21,SignalP5:21,SignalP6:21
g10675  SignalP3NN:23
g10692  SignalP3HMM:18,Phobius:18
g10783  SignalP3NN:18,SignalP3HMM:18,SignalP4:18,SignalP5:18,DeepSig:18,Phobius:18,SignalP6:18
g10784  SignalP3NN:18,SignalP3HMM:18,SignalP4:18,SignalP5:18,DeepSig:18,Phobius:18,SignalP6:18
g10794  SignalP3NN:33
g1106 SignalP3HMM:16,Phobius:19,SignalP3NN:26
g11286  SignalP3NN:19,SignalP3HMM:20,SignalP4:20,SignalP5:20,DeepSig:20,Phobius:20,SignalP6:20
g11369  SignalP3NN:20,SignalP3HMM:20,SignalP5:20,Phobius:36
g11409  SignalP3NN:17,SignalP4:17,DeepSig:17,SignalP5:19,SignalP6:19,SignalP3HMM:21,Phobius:21
g11416  SignalP3NN:25,SignalP3HMM:25,SignalP4:25,SignalP5:25,DeepSig:25,Phobius:25,SignalP6:25
g11565  SignalP3NN:19,SignalP3HMM:19,SignalP4:19,DeepSig:19
g11652  SignalP3HMM:30
g11660  SignalP3NN:24,SignalP5:24,Phobius:24,SignalP3HMM:28,SignalP6:28
g11719  Phobius:34
g11727  DeepSig:22,SignalP3NN:23,SignalP3HMM:23,SignalP4:23,SignalP5:23,Phobius:23,SignalP6:23
g11778  SignalP3NN:16,SignalP3HMM:16,SignalP4:16,DeepSig:16,Phobius:16
g11806  SignalP4:19,DeepSig:19,SignalP3NN:23,SignalP5:23
g11807  SignalP3NN:17,Phobius:17,SignalP3HMM:21,SignalP5:21,DeepSig:21,SignalP6:24
g11808  SignalP3NN:19,SignalP3HMM:19,SignalP4:19,SignalP5:19,DeepSig:19,Phobius:19,SignalP6:19
g11809  SignalP3NN:19,SignalP3HMM:19,SignalP4:19,SignalP5:19,DeepSig:19,Phobius:19,SignalP6:19
g1186 SignalP3NN:19,SignalP3HMM:19,SignalP4:19,SignalP5:19,DeepSig:19,Phobius:19,SignalP6:19
g11906  SignalP3NN:22,SignalP3HMM:22,SignalP4:22,SignalP5:22,DeepSig:22,Phobius:22,SignalP6:22
g11998  SignalP3NN:16
g12135  SignalP3NN:24,SignalP3HMM:24,SignalP4:24,SignalP5:24,SignalP6:24,DeepSig:30,Phobius:30
g12284  SignalP3NN:24,SignalP3HMM:24,SignalP4:24,SignalP5:24,DeepSig:24,Phobius:24,SignalP6:24
g12296  SignalP3NN:21,SignalP3HMM:21,SignalP4:21,SignalP5:21,DeepSig:21,Phobius:21,SignalP6:21
g12452  SignalP3NN:20,SignalP3HMM:20,SignalP4:20,DeepSig:20,Phobius:20,SignalP5:22
g12489  Phobius:29
g12500  SignalP3NN:19,SignalP3HMM:19,SignalP4:19,SignalP5:19,DeepSig:19,Phobius:19,SignalP6:19
g12653  SignalP4:15,Phobius:15
g12716  SignalP3NN:19,SignalP3HMM:19,SignalP4:19,SignalP5:19,DeepSig:19,Phobius:19,SignalP6:19
g12776  SignalP3NN:21,SignalP3HMM:21,SignalP4:21,SignalP5:21,DeepSig:21,Phobius:21,SignalP6:21
g13051  SignalP3NN:20,SignalP3HMM:20,SignalP4:20,SignalP5:20,DeepSig:20,Phobius:20,SignalP6:25
g13093  SignalP3HMM:16
g13170  SignalP3NN:16,SignalP3HMM:16,SignalP4:16,SignalP5:16,DeepSig:16,Phobius:16,SignalP6:16
g13217  SignalP3HMM:22
g13219  SignalP3NN:23,SignalP3HMM:23,SignalP4:23,SignalP5:23,DeepSig:23,Phobius:23,SignalP6:23
g13226  Phobius:21,SignalP6:22,SignalP3NN:24,SignalP4:24
g13566  SignalP3NN:21,SignalP3HMM:21,SignalP4:21,SignalP5:21,DeepSig:21,Phobius:21,SignalP6:21
g13811  SignalP3NN:21,SignalP3HMM:21,SignalP4:21,SignalP5:21,DeepSig:21,Phobius:21,SignalP6:21
g13892  Phobius:18
g13921  SignalP3NN:18,SignalP3HMM:18,SignalP5:18,Phobius:18,SignalP6:18
g14015  SignalP4:20,SignalP6:22,SignalP3NN:24,SignalP3HMM:24,SignalP5:24,DeepSig:24,Phobius:24
g14026  SignalP3NN:20,SignalP3HMM:20,SignalP4:20,SignalP5:20,DeepSig:20,Phobius:20,SignalP6:20
g14153  SignalP3NN:19,SignalP3HMM:19,SignalP4:19,SignalP5:19,DeepSig:19,Phobius:19,SignalP6:19
g14154  SignalP3NN:19,SignalP3HMM:19,SignalP4:19,SignalP5:19,DeepSig:19,Phobius:19,SignalP6:19
g14161  SignalP3NN:20,SignalP3HMM:20,SignalP4:20,SignalP5:20,DeepSig:20,Phobius:20,SignalP6:20
g14208  SignalP3NN:16,SignalP3HMM:16,SignalP4:16,SignalP5:16,DeepSig:16,Phobius:16,SignalP6:16
g14239  SignalP3NN:22
g14266  SignalP4:19,SignalP3NN:20,SignalP3HMM:20,SignalP5:20,DeepSig:20,Phobius:20,SignalP6:20
g14430  SignalP3NN:19,SignalP3HMM:19,SignalP4:19,SignalP5:19,DeepSig:19,Phobius:19,SignalP6:19
g14456  Phobius:21
g14463  SignalP3NN:21,SignalP3HMM:21,SignalP4:21,SignalP5:21,DeepSig:21,Phobius:21,SignalP6:21
g14471  SignalP3NN:21,SignalP3HMM:21,SignalP4:21,SignalP5:21,DeepSig:21,Phobius:21,SignalP6:21
g14546  SignalP3NN:24,SignalP3HMM:24,SignalP4:24,SignalP5:24,DeepSig:24,Phobius:24,SignalP6:24
g14551  SignalP3HMM:21,SignalP4:21,SignalP6:21,DeepSig:23,SignalP3NN:24,SignalP5:24
g14690  SignalP3NN:20,SignalP4:20,SignalP5:20,Phobius:20,SignalP6:20
g14900  SignalP3NN:19,SignalP3HMM:19,SignalP4:19,SignalP5:19,DeepSig:19,Phobius:19,SignalP6:29
g14908  SignalP3NN:21,SignalP3HMM:21,SignalP4:21,SignalP5:21,DeepSig:21,SignalP6:21,Phobius:27
g150  Phobius:40
g15093  SignalP3NN:21,SignalP3HMM:21
g1518 SignalP4:11,Phobius:14,DeepSig:20,SignalP3NN:28
g15306  SignalP3NN:19,SignalP3HMM:19,SignalP4:19,SignalP5:19,DeepSig:19,Phobius:19,SignalP6:19
g15465  SignalP3NN:19,SignalP3HMM:19,SignalP4:19,SignalP5:19,DeepSig:19,Phobius:24,SignalP6:24
g15827  Phobius:17
g16211  SignalP3HMM:21,Phobius:21,SignalP6:22,SignalP3NN:24,SignalP4:24
g16340  SignalP3NN:19,SignalP4:19,DeepSig:19,SignalP3HMM:21,SignalP5:21,Phobius:21
g16440  DeepSig:21
g16601  SignalP3HMM:21,Phobius:21
g16680  Phobius:45
g172  SignalP3NN:18,SignalP4:18
g1765 SignalP6:18,SignalP3NN:25,SignalP3HMM:25,SignalP4:25,SignalP5:25,DeepSig:25,Phobius:25
g17771  Phobius:34
g18158  SignalP3HMM:22,Phobius:22
g1899 Phobius:28
g1901 SignalP3NN:19,SignalP3HMM:19,SignalP5:19,DeepSig:19,SignalP6:19
g2067 Phobius:18,SignalP3NN:19,SignalP3HMM:19,SignalP4:19,SignalP6:19,SignalP5:21,DeepSig:21
g2071 Phobius:24
g2085 SignalP3NN:21,SignalP3HMM:21,SignalP4:21,SignalP5:21,DeepSig:21,Phobius:21,SignalP6:21
g2280 SignalP3HMM:20,Phobius:23
g2297 SignalP3NN:20,SignalP3HMM:20,SignalP4:20,SignalP5:20,DeepSig:20,Phobius:20,SignalP6:20
g254  SignalP3NN:22,SignalP4:22,SignalP5:22,DeepSig:22,SignalP6:22,SignalP3HMM:25,Phobius:25
g255  SignalP3NN:20,SignalP3HMM:20,SignalP4:20,SignalP5:20,DeepSig:20,Phobius:20,SignalP6:20
g258  SignalP3NN:20,SignalP3HMM:20,SignalP4:20,SignalP5:20,DeepSig:20,Phobius:20,SignalP6:20
g259  SignalP3NN:20,SignalP3HMM:20,SignalP4:20,SignalP5:20,DeepSig:20,Phobius:20,SignalP6:20
g2761 SignalP3NN:29,SignalP4:29,SignalP5:29,SignalP6:29
g2902 SignalP3NN:24,Phobius:27
g3072 SignalP3NN:21,SignalP3HMM:21,SignalP4:21,SignalP5:21,DeepSig:21,SignalP6:21,Phobius:26
g3170 SignalP3NN:22,SignalP3HMM:22,SignalP4:22,SignalP5:22,DeepSig:22,Phobius:22,SignalP6:22
g3172 SignalP3NN:22,SignalP3HMM:22,SignalP4:22,SignalP5:22,DeepSig:22,Phobius:22,SignalP6:22
g3204 SignalP3NN:21,SignalP3HMM:21,SignalP4:21,SignalP5:21,DeepSig:21,Phobius:21,SignalP6:21
g3237 SignalP3NN:19,SignalP3HMM:19,SignalP4:19,SignalP5:19,DeepSig:19,Phobius:19,SignalP6:19
g3319 SignalP3NN:19,SignalP3HMM:19,SignalP4:19,SignalP5:19,DeepSig:19,Phobius:19,SignalP6:19
g3921 SignalP3NN:19,SignalP3HMM:19,SignalP4:19,SignalP5:19,DeepSig:19,Phobius:19,SignalP6:19
g4016 SignalP3NN:18
g4176 SignalP3NN:20,SignalP3HMM:20,SignalP4:20,SignalP5:20,DeepSig:20,Phobius:20,SignalP6:20
g4600 SignalP3NN:15,SignalP3HMM:15,DeepSig:18,Phobius:18,SignalP5:24
g4897 SignalP3NN:23,SignalP3HMM:23,SignalP4:23,SignalP5:23,DeepSig:23,Phobius:23,SignalP6:23
g5141 SignalP3NN:24,SignalP3HMM:24,SignalP4:24,SignalP5:24,DeepSig:24,Phobius:24,SignalP6:24
g5142 SignalP3NN:23,SignalP3HMM:23,SignalP4:23,SignalP5:23,DeepSig:23,Phobius:23,SignalP6:23
g5372 SignalP3NN:21,SignalP3HMM:21,SignalP4:21,SignalP5:21,DeepSig:21,Phobius:21,SignalP6:21
g5674 SignalP3NN:23,SignalP3HMM:23,SignalP4:23,SignalP5:23,Phobius:23,SignalP6:23
g5692 SignalP3NN:22,SignalP3HMM:22,SignalP4:22
g5745 SignalP3NN:17,SignalP3HMM:17,SignalP4:17,SignalP5:17,DeepSig:17,SignalP6:17,Phobius:19
g5831 SignalP3NN:23,SignalP3HMM:23,SignalP4:23,SignalP5:23,DeepSig:23,Phobius:23,SignalP6:23
g6026 SignalP3NN:20,SignalP3HMM:20,SignalP4:20,SignalP5:20,DeepSig:20,Phobius:20,SignalP6:20
g6027 SignalP3NN:22,SignalP3HMM:22,SignalP4:22,SignalP5:22,DeepSig:22,Phobius:22,SignalP6:22
g6105 SignalP3NN:18,SignalP3HMM:18,SignalP4:18,SignalP5:18,DeepSig:18,Phobius:18,SignalP6:18
g6162 SignalP6:18,SignalP3NN:21,SignalP3HMM:21,SignalP4:21,SignalP5:21,DeepSig:21,Phobius:21
g6243 SignalP3NN:26,DeepSig:26
g6281 Phobius:13
g6315 SignalP3HMM:17,Phobius:17,DeepSig:18
g6354 SignalP3NN:20,SignalP4:20,DeepSig:20,Phobius:20,SignalP5:25,SignalP3HMM:27,SignalP6:27
g6355 SignalP3NN:20,SignalP3HMM:20,SignalP4:20,SignalP5:20,DeepSig:20,Phobius:20,SignalP6:20
g6406 SignalP3HMM:17
g6458 SignalP3NN:28
g6590 SignalP3NN:19,SignalP3HMM:19,SignalP4:19,SignalP5:19,DeepSig:19,Phobius:19,SignalP6:19
g6999 SignalP3NN:19,SignalP3HMM:19,SignalP4:19,SignalP5:19,DeepSig:19,Phobius:19,SignalP6:19
g7051 SignalP3NN:20,SignalP3HMM:20,SignalP4:20,SignalP5:20,DeepSig:20,SignalP6:20
g7120 SignalP3NN:20,SignalP3HMM:20,SignalP4:20,SignalP5:20,DeepSig:20,Phobius:20,SignalP6:20
g7121 SignalP3NN:20,SignalP3HMM:20,SignalP4:20,SignalP5:20,DeepSig:20,Phobius:20,SignalP6:20
g7132 DeepSig:25
g7198 SignalP3NN:18
g7200 SignalP3NN:19,SignalP3HMM:19,SignalP4:19,SignalP5:19,DeepSig:19,Phobius:19,SignalP6:19
g7361 SignalP3NN:23,SignalP3HMM:24,SignalP4:24,SignalP5:24,DeepSig:24,Phobius:24,SignalP6:24
g7603 SignalP3NN:21,SignalP3HMM:21,SignalP4:21,SignalP5:21,DeepSig:21,Phobius:21,SignalP6:21
g7646 SignalP3NN:21,SignalP3HMM:21,SignalP4:21,SignalP5:21,DeepSig:21,Phobius:21,SignalP6:21
g7647 Phobius:18,SignalP3NN:20,SignalP3HMM:20,SignalP4:20,SignalP5:20,DeepSig:20,SignalP6:20
g7696 SignalP3NN:23,SignalP4:23,DeepSig:23
g7802 SignalP3NN:20,SignalP3HMM:20,SignalP4:20,SignalP5:20,DeepSig:20,Phobius:20,SignalP6:20
g7927 SignalP3NN:21,SignalP4:21,SignalP5:21,DeepSig:21,Phobius:21,SignalP6:21,SignalP3HMM:23
g7928 SignalP3NN:21,SignalP3HMM:21,SignalP4:21,SignalP5:21,DeepSig:21,Phobius:21,SignalP6:21
g7935 SignalP3NN:20,SignalP3HMM:21,SignalP4:21,SignalP5:21,DeepSig:21,Phobius:21,SignalP6:21
g7937 SignalP3NN:21,SignalP3HMM:21,SignalP4:21,SignalP5:21,DeepSig:21,Phobius:21,SignalP6:21
g8009 SignalP3NN:20,SignalP3HMM:20,SignalP4:20,DeepSig:20,SignalP5:22
g8035 SignalP3NN:17,SignalP3HMM:17,SignalP4:17,SignalP5:17,DeepSig:17,SignalP6:17,Phobius:19
g8081 SignalP3HMM:19,DeepSig:19,Phobius:19,SignalP3NN:21,SignalP4:21,SignalP5:21,SignalP6:21
g8248 SignalP4:21,SignalP5:21,DeepSig:21,SignalP6:21,SignalP3HMM:22,SignalP3NN:25,Phobius:25
g8278 SignalP3NN:17,SignalP4:17,SignalP5:17,DeepSig:17,SignalP6:17,SignalP3HMM:21,Phobius:21
g8734 Phobius:20
g8760 SignalP3NN:19,SignalP3HMM:19,SignalP4:19,SignalP5:19,DeepSig:19,Phobius:19,SignalP6:19
g8897 Phobius:27
g8974 SignalP3NN:14,SignalP3HMM:14,SignalP5:14,Phobius:14
g9008 SignalP3NN:18,SignalP3HMM:18,SignalP5:18,Phobius:18
g9065 SignalP3NN:27
g911  SignalP3NN:19,SignalP3HMM:19,SignalP4:19,SignalP5:19,DeepSig:19,SignalP6:19,Phobius:22
g9132 SignalP3NN:22,SignalP3HMM:22,SignalP4:22,SignalP5:22,Phobius:22
g9487 SignalP3HMM:19,Phobius:19
g978  SignalP3HMM:17,Phobius:17,SignalP3NN:18
g981  SignalP3NN:20,SignalP3HMM:20,SignalP4:20,SignalP5:20,DeepSig:20,Phobius:20,SignalP6:20
g9883 Phobius:21,SignalP3NN:23,SignalP3HMM:23,SignalP4:23,SignalP5:23,DeepSig:23
g9979 SignalP3NN:13,SignalP3HMM:13,SignalP4:13,SignalP5:13,Phobius:13
```txt
Aphanis_Gene: Signal predicted:
g10232  23
g11054  22
g11278  22
g1209 21
g12302  20
g12450  19
g12453  20
g12628  17
g12667 19
g12955  17
g13367  22
g13490  22
g1355 17
g13885  23
g14348  24
g14731  14
g14832  19
g15923  27
g2113 18
g2919 29
g2928 16
g3129 21
g324  27
g422  30
g548  20
g5605 20
g5645 16
g5721 28
g6055 24
g6450 28
g6520  19
g6562 18
g7278 23
g757  10
g8019 20
g9359 31
  
Leuco_Gene: 
g10095  24
g1010 18
g10145  23
g10372  20
g10489  19
g10497  19
g10500  21
g10675  23
g10692  18
g10783  18
g10784  18
g10794  33
g1106 26
g11286  20
g11369  20
g11409  17
g11416  25
g11565  19
g11652  30
g11660  24
g11719  34
g11727  23
g11778  16
g11806  23
g11807  17
g11808  19
g11809  19
g1186 19
g11906  22
g11998  16
g12135  24
g12284  24
g12296  21
g12452  20
g12489  29
g12500  19
g12653  15
g12716  19
g12776  21
g13051  20
g13093  16
g13170  16
g13217  22
g13219  23
g13226  24
g13566  21
g13811  21
g13892  18
g13921  18
g14015  20
g14026  20
g14153  19
g14154  19
g14161  20
g14208  16
g14239  22
g14266  19
g14430  19
g14456  21
g14463  21
g14471  21
g14546  24
g14551  21
g14690  20
g14900  19
g14908  21
g150  40
g15093  21
g1518 11
g15306  19
g15465  19
g15827  17
g16211  24
g16340  19
g16440  21
g16601  21
g16680  45
g172  18
g1765 25
g17771  34
g18158  22
g1899 28
g1901 19
g2067 19
g2071 24
g2085 21
g2280 20
g2297 20
g254  22
g255  20
g258  20
g259  20
g2761 29
g2902 24
g3072 21
g3170 22
g3172 22
g3204 21
g3237 19
g3319 19
g3921 19
g4016 18
g4176 20
g4600 15
g4897 23
g5141 24
g5142 23
g5372 21
g5674 23
g5692 22
g5745 17
g5831 23
g6026 20
g6027 22
g6105 18
g6162 21
g6243 26
g6281 13
g6315 17
g6354 20
g6355 20
g6406 17
g6458 28
g6590 19
g6999 19
g7051 20
g7120 20
g7121 20
g7132 25
g7198 18
g7200 19
g7361 24
g7603 21
g7646 21
g7647 20
g7696 23
g7802 20
g7927 21
g7928 21
g7935 21
g7937 21
g8009 20
g8035 17
g8081 21
g8248 21
g8278 17
g8734 20
g8760 19
g8897 27
g8974 14
g9008 18
g9065 27
g911  199
g9132 22
g9487 19
g978  18
g981  20
g9883 23
g9979 13

```
```bash
nano aphanis-unique-cseps-manual.out
nano leucotricha-unique-cseps-manual.out

cat effectors/*.faa > aphanis-unique-cseps.faa
cat effectors2/*.faa > leucotricha-unique-cseps.faa

mkdir aphanis-unique-cseps-cleaved
mkdir leucotricha-unique-cseps-cleaved
```
```python
from Bio import SeqIO

#signalp_file = "aphanis-unique-cseps-manual.out"
#input_fasta = "aphanis-unique-cseps.faa"
#output_fasta = "aphanis-unique-cseps-cleaved.faa"
signalp_file = "leucotricha-unique-cseps-manual.out"
input_fasta = "leucotricha-unique-cseps.faa"
output_fasta = "leucotricha-unique-cseps-cleaved.faa"

def load_signal_peptide_info(signalp_file):
    cleavage_sites = {}
    with open(signalp_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            seq_id = parts[0]
            cleavage_site = int(parts[1])
            cleavage_sites[seq_id] = cleavage_site
    return cleavage_sites

def remove_signal_peptide_fasta(input_fasta, output_fasta, cleavage_sites):
    with open(output_fasta, 'w') as output_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq_id = record.id.split('.')[0]  # To match with IDs in cleavage_sites
            if seq_id in cleavage_sites:
                cleavage_site = cleavage_sites[seq_id]
                mature_sequence = record.seq[cleavage_site:]
                record.seq = mature_sequence
            SeqIO.write(record, output_handle, "fasta")

cleavage_sites = load_signal_peptide_info(signalp_file)
remove_signal_peptide_fasta(input_fasta, output_fasta, cleavage_sites)



from Bio import SeqIO
import os
def split_fasta(input_fasta, output_directory):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    for record in SeqIO.parse(input_fasta, "fasta"):
        output_file = os.path.join(output_directory, f"{record.id}.fasta")
        with open(output_file, 'w') as output_handle:
            SeqIO.write(record, output_handle, "fasta")
        print(f"Written {record.id} to {output_file}")

#input_fasta = "aphanis-unique-cseps-cleaved.faa"
#output_directory = "aphanis-unique-cseps-cleaved"
input_fasta = "leucotricha-unique-cseps-cleaved.faa"
output_directory = "leucotricha-unique-cseps-cleaved"

split_fasta(input_fasta, output_directory)
```
```bash
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/aphanis-unique-cseps-cleaved/*fasta); do
echo $(basename $file)
grep -v '>' $file | wc -c
done

for protein in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/aphanis-unique-cseps-cleaved/*.fasta); do
sed -i 's/X*$//' $protein
OutDir=$(dirname $protein)/$(basename $protein | sed 's@.t1.fasta@@g')
mkdir $OutDir
sbatch ~/git_repos/Wrappers/NBI/alphafold.sh $protein $OutDir
done
#3555477-3555523 - gpu node with 8cpu, multiple jobs are running on each node

for protein in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/leucotricha-unique-cseps-cleaved/g911.t1.fasta); do
sed -i 's/X*$//' $protein
OutDir=$(dirname $protein)/$(basename $protein | sed 's@.t1.fasta@@g')
mkdir $OutDir
sbatch ~/git_repos/Wrappers/NBI/alphafold.sh $protein $OutDir
done
#3556341-3556499,6431854 - gpu node with 8cpu, multiple jobs are running on each node

mkdir temp_AF2  
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/*-unique-cseps-cleaved/*/*/ranked_0.pdb); do
ID=$(echo $file | cut -d '/' -f8 | cut -d '-' -f1 )_$(echo $file | cut -d '/' -f9)_$(basename $file)
cp $file temp_AF2/$ID
done

#Download .pdb files from seong and krasileva paper:
cd /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/Seong+Krasileva
for file in $(ls *.tar.gz); do
tar -xvzf $file
done
mkdir all
cp *PDBs/*.pdb all/.
cd ..


for file in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/Seong+Krasileva/all/*.pdb
do
    ln -s "$file" temp_AF2/.
done

# -c and --cov-mode refer to structural alignments not amino acid alignments
# --min-seq-id is for amio acid similiarity
#--lddt-threshold sets the minimum lDDT (Local Distance Difference Test) score required for an alignment to be considered for clustering, --tmscore-threshold sets the minimum TM-score (Template Modeling score) required for an alignment to be considered for clustering. TM-score focuses on the overall structural resemblance between two proteins, while lDDT evaluates the similarity at a more localized level, considering the arrangement of atoms within specific regions of the structures.
source package cda29b6a-320e-4d73-83c6-240ed7a6201e
pdb_dir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/temp_AF2
cpu=1  
Prefix=cluster-0.9
foldseek easy-cluster --threads $cpu -s 7.5 --max-seqs 2000 --min-seq-id 0.6 -c 0.9 $pdb_dir $Prefix tmp
Prefix=cluster-0.8
foldseek easy-cluster --threads $cpu  -s 7.5 --max-seqs 2000 --min-seq-id 0.6 -c 0.8 $pdb_dir $Prefix tmp
Prefix=cluster-0.7
foldseek easy-cluster --threads $cpu  -s 7.5 --max-seqs 2000 --min-seq-id 0.6 -c 0.7 $pdb_dir $Prefix tmp
Prefix=cluster-0.6
foldseek easy-cluster --threads $cpu  -s 7.5 --max-seqs 2000 --min-seq-id 0.6 -c 0.6 $pdb_dir $Prefix tmp

Prefix=cluster-0.9-5
foldseek easy-cluster --threads $cpu -s 7.5 --max-seqs 2000 --min-seq-id 0.5 -c 0.9 $pdb_dir $Prefix tmp
Prefix=cluster-0.8-5
foldseek easy-cluster --threads $cpu  -s 7.5 --max-seqs 2000 --min-seq-id 0.5 -c 0.8 $pdb_dir $Prefix tmp
Prefix=cluster-0.7-5
foldseek easy-cluster --threads $cpu  -s 7.5 --max-seqs 2000 --min-seq-id 0.5 -c 0.7 $pdb_dir $Prefix tmp
Prefix=cluster-0.6-5
foldseek easy-cluster --threads $cpu  -s 7.5 --max-seqs 2000 --min-seq-id 0.5 -c 0.6 $pdb_dir $Prefix tmp
Prefix=cluster-0.5-5
foldseek easy-cluster --threads $cpu  -s 7.5 --max-seqs 2000 --min-seq-id 0.5 -c 0.5 $pdb_dir $Prefix tmp

Prefix=cluster-0.9-4
foldseek easy-cluster --threads $cpu -s 7.5 --max-seqs 2000 --min-seq-id 0.4 -c 0.9 $pdb_dir $Prefix tmp
Prefix=cluster-0.8-4
foldseek easy-cluster --threads $cpu  -s 7.5 --max-seqs 2000 --min-seq-id 0.4 -c 0.8 $pdb_dir $Prefix tmp
Prefix=cluster-0.7-4
foldseek easy-cluster --threads $cpu  -s 7.5 --max-seqs 2000 --min-seq-id 0.4 -c 0.7 $pdb_dir $Prefix tmp
Prefix=cluster-0.6-4
foldseek easy-cluster --threads $cpu  -s 7.5 --max-seqs 2000 --min-seq-id 0.4 -c 0.6 $pdb_dir $Prefix tmp

Prefix=cluster
foldseek easy-cluster --threads $cpu --alignment-type 2 -s 7.5 --max-seqs 99999 --min-seq-id 0.000 --cov-mode 0 -c 0.5 --lddt-threshold 0 --tmscore-threshold 0.5 $pdb_dir $Prefix tmp
Prefix=clusterx-4
foldseek easy-cluster --threads $cpu --alignment-type 2 -s 7.5 --max-seqs 99999 --min-seq-id 0.000 --cov-mode 0 -c 0.4 --lddt-threshold 0 --tmscore-threshold 0.5 $pdb_dir $Prefix tmp
Prefix=clusterx-4-4
foldseek easy-cluster --threads $cpu --alignment-type 2 -s 7.5 --max-seqs 99999 --min-seq-id 0.000 --cov-mode 0 -c 0.4 --lddt-threshold 0 --tmscore-threshold 0.4 $pdb_dir $Prefix tmp
Prefix=clusterx-0-5
foldseek easy-cluster --threads $cpu --alignment-type 2 -s 7.5 --max-seqs 99999 --min-seq-id 0.000 --cov-mode 0 -c 0.0 --lddt-threshold 0 --tmscore-threshold 0.5 $pdb_dir $Prefix tmp

#make a foldseek database from the Seong and Krasileva pdbs
source package cda29b6a-320e-4d73-83c6-240ed7a6201e
foldseek createdb /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/Seong+Krasileva/all /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/Seong+Krasileva/S+Kdb

#The cost to run with the afdb database is ~£2 per protein, therefore run only for the custom database and submit through the web server - time consuming :(
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/*-unique-cseps-cleaved/*/*/ranked_0.pdb | grep -v 'leucotricha-unique-cseps-cleaved/g9979/g9979.t1'); do
InFile=$file
OutDir=$(dirname $file)/foldseek
mkdir $OutDir
#Database=/nbi/Reference-Data/Foldseek/afdb
#OutFile=$(echo $file | cut -d '/' -f8 | cut -d '-' -f1 )_$(echo $file | cut -d '/' -f9)_$(basename $file | sed 's@.pdb@@g')-$(basename $Database)-foldseek-output.txt
#sbatch ~/git_repos/Wrappers/NBI/foldseek.sh $InFile $Database $OutDir $OutFile
#Database=/nbi/Reference-Data/Foldseek/afdb_proteome
#OutFile=$(echo $file | cut -d '/' -f8 | cut -d '-' -f1 )_$(echo $file | cut -d '/' -f9)_$(basename $file | sed 's@.pdb@@g')-$(basename $Database)-foldseek-output.txt
#sbatch ~/git_repos/Wrappers/NBI/foldseek.sh $InFile $Database $OutDir $OutFile
#Database=/nbi/Reference-Data/Foldseek/afdb_swissprot
#OutFile=$(echo $file | cut -d '/' -f8 | cut -d '-' -f1 )_$(echo $file | cut -d '/' -f9)_$(basename $file | sed 's@.pdb@@g')-$(basename $Database)-foldseek-output.txt
#sbatch ~/git_repos/Wrappers/NBI/foldseek.sh $InFile $Database $OutDir $OutFile
#Database=/nbi/Reference-Data/Foldseek/afdb_uniprot50
#OutFile=$(echo $file | cut -d '/' -f8 | cut -d '-' -f1 )_$(echo $file | cut -d '/' -f9)_$(basename $file | sed 's@.pdb@@g')-$(basename $Database)-foldseek-output.txt
#sbatch ~/git_repos/Wrappers/NBI/foldseek.sh $InFile $Database $OutDir $OutFile
#Database=/nbi/Reference-Data/Foldseek/cath50
#OutFile=$(echo $file | cut -d '/' -f8 | cut -d '-' -f1 )_$(echo $file | cut -d '/' -f9)_$(basename $file | sed 's@.pdb@@g')-$(basename $Database)-foldseek-output.txt
#sbatch ~/git_repos/Wrappers/NBI/foldseek.sh $InFile $Database $OutDir $OutFile
#Database=/nbi/Reference-Data/Foldseek/esmatlas30
#OutFile=$(echo $file | cut -d '/' -f8 | cut -d '-' -f1 )_$(echo $file | cut -d '/' -f9)_$(basename $file | sed 's@.pdb@@g')-$(basename $Database)-foldseek-output.txt
#sbatch ~/git_repos/Wrappers/NBI/foldseek.sh $InFile $Database $OutDir $OutFile
#Database=/nbi/Reference-Data/Foldseek/pdb
#OutFile=$(echo $file | cut -d '/' -f8 | cut -d '-' -f1 )_$(echo $file | cut -d '/' -f9)_$(basename $file | sed 's@.pdb@@g')-$(basename $Database)-foldseek-output.txt
#sbatch ~/git_repos/Wrappers/NBI/foldseek.sh $InFile $Database $OutDir $OutFile
Database=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/Seong+Krasileva/S+Kdb
OutFile=$(echo $file | cut -d '/' -f8 | cut -d '-' -f1 )_$(echo $file | cut -d '/' -f9)_$(basename $file | sed 's@.pdb@@g')-$(basename $Database)-foldseek-output.txt
sbatch ~/git_repos/Wrappers/NBI/foldseek.sh $InFile $Database $OutDir $OutFile
done
#3751068-3751075
#3780842-3782679
#3782680-3782872

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/*-unique-cseps-cleaved/*/*/ranked_0.pdb); do
result=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/temp_AF2/foldseek_web_server/mmseqs_results_$(echo $file | cut -d '/' -f8 | cut -d '-' -f1 )_$(echo $file | cut -d '/' -f9)_ranked_0.gz
mkdir $(dirname $file)/foldseek_web_server
mv $result $(dirname $file)/foldseek_web_server/$(basename $result | sed 's@.gz@.tar.gz@g')
ls $(dirname $file)/foldseek_web_server/$(basename $result)
cd $(dirname $file)/foldseek_web_server/
tar -xvzf $(basename $result | sed 's@.gz@.tar.gz@g')
done

for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/*-unique-cseps-cleaved/*/*/foldseek_web_server/*.m8); do
db=$(basename $file | cut -d '_' -f2 | cut -d '.' -f1)
awk -v db="$db" 'BEGIN{FS=OFS="\t"} {print $0, db}' $file > temp.tsv && mv temp.tsv $file
done

for dir in $(ls -d /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/*-unique-cseps-cleaved/*/*/foldseek_web_server); do
cat $dir/*.m8 > $dir/all_db.tsv
sort -k11,11nr -k12,12g $dir/all_db.tsv > temp.tsv && mv temp.tsv $dir/all_db.tsv
done

ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/*-unique-cseps-cleaved/*/*/foldseek_web_server/all_db.tsv
```
```bash
#mot all proteins have signal peptide as predictd by signalp4.1, some are only predicted by signalp3, phobius or deepsig. Where there was dissagreement between different prediction methods signalp4.1 was considered most reliable, followed by singalp3NN, then signalp3HMM, then phobius, then deepsig - however inmost cases the different tools were in agreement with one another.

mkdir effectors3
mkdir effectors4
mkdir effectors5
mkdir effectors6
mkdir effectors7
mkdir effectors8

while IFS=$'\t' read -r line; do
    awk -F'\t' '$31 == "1" && $27 == "." && $7 == "0" {print $1}' <<< "$line" | grep -v 'name' >> aphan_common_secreted_effectors.txt
done < /home/theaven/projects/niab/theaven/results/P_aphanis-THeavenDRCT72020_1-common-ranked4.tsv

while IFS=$'\t' read -r line; do
echo $line > hits.txt
python3 /home/theaven/scratch/apps/tools/seq_get.py --id_file hits.txt --input /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta --output effectors3/${line}.faa
done < aphan_common_secreted_effectors.txt

while IFS=$'\t' read -r line; do
    awk -F'\t' '$31 == "1" && $7 == "0" {print $1}' <<< "$line" | grep -v 'name' | grep -v -F -f aphan_common_secreted_effectors.txt >> aphan_common_secreted_effectors2.txt
done < /home/theaven/projects/niab/theaven/results/P_aphanis-THeavenDRCT72020_1-common-ranked4.tsv


while IFS=$'\t' read -r line; do
    awk -F'\t' '$31 == "1" && $7 == "0" {print $1}' <<< "$line" | grep -v 'name' | grep -v -F -f aphanis_only_secreted_effectors.txt >> aphanis_only_secreted_effectors2.txt 
done < /home/theaven/projects/niab/theaven/results/P_aphanis-THeavenDRCT72020_1-aphanis-ranked4.tsv

while IFS=$'\t' read -r line; do
echo $line > hits.txt
python3 /home/theaven/scratch/apps/tools/seq_get.py --id_file hits.txt --input /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta --output effectors5/${line}.faa
done < aphan_common_secreted_effectors2.txt

while IFS=$'\t' read -r line; do
echo $line > hits.txt
python3 /home/theaven/scratch/apps/tools/seq_get.py --id_file hits.txt --input /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta --output effectors7/${line}.faa
done < aphanis_only_secreted_effectors2.txt 

awk -F'\t' '{value = ""; if ($54 ~ /SignalP4:/) {match($54, /SignalP4:([0-9]+)/, arr); value = arr[1];} else if ($54 ~ /SignalP3NN:/) {match($54, /SignalP3NN:([0-9]+)/, arr); value = arr[1];} else if ($54 ~ /SignalP3HMM:/) {match($54, /SignalP3HMM:([0-9]+)/, arr);value = arr[1];} else if ($54 ~ /Phobius:/) {match($54, /Phobius:([0-9]+)/, arr);value = arr[1];} else if ($54 ~ /DeepSig:/) {match($54, /DeepSig:([0-9]+)/, arr);value = arr[1];}print $1 "\t" value;}' /home/theaven/projects/niab/theaven/results/P_aphanis-THeavenDRCT72020_1-aphanis-ranked4.tsv | grep -F -f aphanis_only_secreted_effectors2.txt > aphanis_only_secreted_effectors2.out
awk -F'\t' '{value = ""; if ($54 ~ /SignalP4:/) {match($54, /SignalP4:([0-9]+)/, arr); value = arr[1];} else if ($54 ~ /SignalP3NN:/) {match($54, /SignalP3NN:([0-9]+)/, arr); value = arr[1];} else if ($54 ~ /SignalP3HMM:/) {match($54, /SignalP3HMM:([0-9]+)/, arr);value = arr[1];} else if ($54 ~ /Phobius:/) {match($54, /Phobius:([0-9]+)/, arr);value = arr[1];} else if ($54 ~ /DeepSig:/) {match($54, /DeepSig:([0-9]+)/, arr);value = arr[1];}print $1 "\t" value;}' /home/theaven/projects/niab/theaven/results/P_aphanis-THeavenDRCT72020_1-common-ranked4.tsv | grep -F -f aphan_common_secreted_effectors.txt > aphan_common_secreted_effectors.out
awk -F'\t' '{value = ""; if ($54 ~ /SignalP4:/) {match($54, /SignalP4:([0-9]+)/, arr); value = arr[1];} else if ($54 ~ /SignalP3NN:/) {match($54, /SignalP3NN:([0-9]+)/, arr); value = arr[1];} else if ($54 ~ /SignalP3HMM:/) {match($54, /SignalP3HMM:([0-9]+)/, arr);value = arr[1];} else if ($54 ~ /Phobius:/) {match($54, /Phobius:([0-9]+)/, arr);value = arr[1];} else if ($54 ~ /DeepSig:/) {match($54, /DeepSig:([0-9]+)/, arr);value = arr[1];} else if ($54 ~ /SignalP6:/) {match($54, /SignalP6:([0-9]+)/, arr);value = arr[1];} else if ($54 ~ /SignalP5:/) {match($54, /SignalP5:([0-9]+)/, arr);value = arr[1];}print $1 "\t" value;}' /home/theaven/projects/niab/theaven/results/P_aphanis-THeavenDRCT72020_1-common-ranked4.tsv | grep -F -f aphan_common_secreted_effectors2.txt > aphan_common_secreted_effectors2.out

cat effectors3/*.faa > aphan_common_secreted_effectors.faa
cat effectors5/*.faa > aphan_common_secreted_effectors2.faa
cat effectors7/*.faa > aphanis_only_secreted_effectors2.faa

python3 ../apps/tools/remove_signalp.py aphan_common_secreted_effectors.out aphan_common_secreted_effectors.faa
python3 ../apps/tools/remove_signalp.py aphan_common_secreted_effectors2.out aphan_common_secreted_effectors2.faa 
python3 ../apps/tools/remove_signalp.py aphanis_only_secreted_effectors2.out aphanis_only_secreted_effectors2.faa

for protein in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/aphan_common_secreted_effectors-cleaved/*.fasta); do
      Jobs=$(squeue -u did23faz| grep 'AF2'  | wc -l)
    echo x
    while [ $Jobs -gt 9 ]; do
      sleep 900s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'AF2'  | wc -l)
    done
OutDir=$(dirname $protein)/$(basename $protein | sed 's@.t1.fasta@@g')
mkdir $OutDir
if [ ! ${OutDir}/$(basename $protein | sed 's@.fasta@@g')/ranked_0.pdb ]; then
sed -i 's/X*$//' $protein
sbatch ~/git_repos/Wrappers/NBI/alphafold.sh $protein $OutDir
else
echo done for $(basename $protein | sed 's@.fasta@@g')
fi
done
############################################################################################################################
############################################################################################################################

while IFS=$'\t' read -r line; do
    awk -F'\t' '$31 == "1" && $27 == "." && $7 == "0" {print $1}' <<< "$line" | grep -v 'name' >> leuco_common_secreted_effectors.txt
done < /home/theaven/projects/niab/theaven/results/P_leucotricha-THeavenp11_1-common-ranked4.tsv

while IFS=$'\t' read -r line; do
echo $line > hits.txt
python3 /home/theaven/scratch/apps/tools/seq_get.py --id_file hits.txt --input /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta --output effectors4/${line}.faa
done < leuco_common_secreted_effectors.txt

while IFS=$'\t' read -r line; do
    awk -F'\t' '$31 == "1" && $7 == "0" {print $1}' <<< "$line" | grep -v 'name' | grep -v -F -f leuco_common_secreted_effectors.txt >> leuco_common_secreted_effectors2.txt
done < /home/theaven/projects/niab/theaven/results/P_leucotricha-THeavenp11_1-common-ranked4.tsv

while IFS=$'\t' read -r line; do
    awk -F'\t' '$31 == "1" && $7 == "0" {print $1}' <<< "$line" | grep -v 'name' | grep -v -F -f leuco_only_secreted_effectors.txt >> leuco_only_secreted_effectors2.txt
done < /home/theaven/projects/niab/theaven/results/P_leucotricha-THeavenp11_1-leucotricha-ranked4.tsv

while IFS=$'\t' read -r line; do
echo $line > hits.txt
python3 /home/theaven/scratch/apps/tools/seq_get.py --id_file hits.txt --input /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta --output effectors6/${line}.faa
done < leuco_common_secreted_effectors2.txt

while IFS=$'\t' read -r line; do
echo $line > hits.txt
python3 /home/theaven/scratch/apps/tools/seq_get.py --id_file hits.txt --input /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.pep.fasta --output effectors8/${line}.faa
done < leuco_only_secreted_effectors2.txt

awk -F'\t' '{value = ""; if ($54 ~ /SignalP4:/) {match($54, /SignalP4:([0-9]+)/, arr); value = arr[1];} else if ($54 ~ /SignalP3NN:/) {match($54, /SignalP3NN:([0-9]+)/, arr); value = arr[1];} else if ($54 ~ /SignalP3HMM:/) {match($54, /SignalP3HMM:([0-9]+)/, arr);value = arr[1];} else if ($54 ~ /Phobius:/) {match($54, /Phobius:([0-9]+)/, arr);value = arr[1];} else if ($54 ~ /DeepSig:/) {match($54, /DeepSig:([0-9]+)/, arr);value = arr[1];}print $1 "\t" value;}' /home/theaven/projects/niab/theaven/results/P_leucotricha-THeavenp11_1-leucotricha-ranked4.tsv | grep -F -f leuco_only_secreted_effectors2.txt > leuco_only_secreted_effectors2.out
awk -F'\t' '{value = ""; if ($54 ~ /SignalP4:/) {match($54, /SignalP4:([0-9]+)/, arr); value = arr[1];} else if ($54 ~ /SignalP3NN:/) {match($54, /SignalP3NN:([0-9]+)/, arr); value = arr[1];} else if ($54 ~ /SignalP3HMM:/) {match($54, /SignalP3HMM:([0-9]+)/, arr);value = arr[1];} else if ($54 ~ /Phobius:/) {match($54, /Phobius:([0-9]+)/, arr);value = arr[1];} else if ($54 ~ /DeepSig:/) {match($54, /DeepSig:([0-9]+)/, arr);value = arr[1];}print $1 "\t" value;}' /home/theaven/projects/niab/theaven/results/P_leucotricha-THeavenp11_1-common-ranked4.tsv | grep -F -f leuco_common_secreted_effectors.txt > leuco_common_secreted_effectors.out
awk -F'\t' '{value = ""; if ($54 ~ /SignalP4:/) {match($54, /SignalP4:([0-9]+)/, arr); value = arr[1];} else if ($54 ~ /SignalP3NN:/) {match($54, /SignalP3NN:([0-9]+)/, arr); value = arr[1];} else if ($54 ~ /SignalP3HMM:/) {match($54, /SignalP3HMM:([0-9]+)/, arr);value = arr[1];} else if ($54 ~ /Phobius:/) {match($54, /Phobius:([0-9]+)/, arr);value = arr[1];} else if ($54 ~ /DeepSig:/) {match($54, /DeepSig:([0-9]+)/, arr);value = arr[1];} else if ($54 ~ /SignalP6:/) {match($54, /SignalP6:([0-9]+)/, arr);value = arr[1];} else if ($54 ~ /SignalP5:/) {match($54, /SignalP5:([0-9]+)/, arr);value = arr[1];}print $1 "\t" value;}' /home/theaven/projects/niab/theaven/results/P_leucotricha-THeavenp11_1-common-ranked4.tsv | grep -F -f leuco_common_secreted_effectors2.txt > leuco_common_secreted_effectors2.out

cat effectors4/*.faa > leuco_common_secreted_effectors.faa
cat effectors6/*.faa > leuco_common_secreted_effectors2.faa
cat effectors8/*.faa > leuco_only_secreted_effectors2.faa

python3 ../apps/tools/remove_signalp.py leuco_only_secreted_effectors2.out leuco_only_secreted_effectors2.faa
python3 ../apps/tools/remove_signalp.py leuco_common_secreted_effectors2.out leuco_common_secreted_effectors2.faa
python3 ../apps/tools/remove_signalp.py leuco_common_secreted_effectors.out leuco_common_secreted_effectors.faa

for protein in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/leuco_common_secreted_effectors-cleaved/*.fasta); do
      Jobs=$(squeue -u did23faz| grep 'AF2'  | wc -l)
    echo x
    while [ $Jobs -gt 9 ]; do
      sleep 900s
      printf "."
      Jobs=$(squeue -u did23faz| grep 'AF2'  | wc -l)
    done
OutDir=$(dirname $protein)/$(basename $protein | sed 's@.t1.fasta@@g')
mkdir $OutDir
if [ ! ${OutDir}/$(basename $protein | sed 's@.fasta@@g')/ranked_0.pdb ]; then
sed -i 's/X*$//' $protein
sbatch ~/git_repos/Wrappers/NBI/alphafold.sh $protein $OutDir
else
echo done for $(basename $protein | sed 's@.fasta@@g')
fi
done
```
```bash
mkdir temp_AF2-2 
for file in $(ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/*-unique-cseps-cleaved/*/*/ranked_0.pdb /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/*_common_secreted_effectors-cleaved/*/*/ranked_0.pdb); do
ID=$(echo $file | cut -d '/' -f8 | cut -d '-' -f1 )_$(echo $file | cut -d '/' -f9)_$(basename $file)
cp $file temp_AF2-2/$ID
done

for file in /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/Seong+Krasileva/all/*.pdb
do
    ln -s "$file" temp_AF2-2/.
done

source package cda29b6a-320e-4d73-83c6-240ed7a6201e
pdb_dir=/jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/temp_AF2-2
cpu=1  

Prefix=cluster
foldseek easy-cluster --threads $cpu --alignment-type 2 -s 7.5 --max-seqs 99999 --min-seq-id 0.000 --cov-mode 0 -c 0.5 --lddt-threshold 0 --tmscore-threshold 0.5 $pdb_dir $Prefix tmp
Prefix=clusterxy-4
foldseek easy-cluster --threads $cpu --alignment-type 2 -s 7.5 --max-seqs 99999 --min-seq-id 0.000 --cov-mode 0 -c 0.4 --lddt-threshold 0 --tmscore-threshold 0.5 $pdb_dir $Prefix tmp
Prefix=clusterxy-4-4
foldseek easy-cluster --threads $cpu --alignment-type 2 -s 7.5 --max-seqs 99999 --min-seq-id 0.000 --cov-mode 0 -c 0.4 --lddt-threshold 0 --tmscore-threshold 0.4 $pdb_dir $Prefix tmp
Prefix=clusterxy-0-5
foldseek easy-cluster --threads $cpu --alignment-type 2 -s 7.5 --max-seqs 99999 --min-seq-id 0.000 --cov-mode 0 -c 0.0 --lddt-threshold 0 --tmscore-threshold 0.5 $pdb_dir $Prefix tmp

ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/aphan_common_secreted_effectors-cleaved/*/*/ranked_0.pdb | wc -l
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/aphanis-unique-cseps-cleaved/*/*/ranked_0.pdb | wc -l
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/leuco_common_secreted_effectors-cleaved/*/*/ranked_0.pdb | wc -l
ls /jic/scratch/groups/Saskia-Hogenhout/tom_heaven/nano_diagnostics/leucotricha-unique-cseps-cleaved/*/*/ranked_0.pdb | wc -l
```
########################################################

```R
DRCT72020 <- read_excel("P_aphanis_THeavenDRCT72020_annotations_master.xlsx")
DRCT72021 <- read_excel("P_aphanis_THeavenDRCT72021_annotations_master.xlsx")
scott2020 <- read_excel("P_aphanis_THeavenSCOTT2020_annotations_master.xlsx")
OGB2019 <- read_excel("P_leucotricha_THeavenpOGB2019_annotations_master.xlsx")
OGB2021 <- read_excel("P_leucotricha_THeavenpOGB2021_annotations_master.xlsx")
p112020 <- read_excel("P_leucotricha_THeavenp11_annotations_master.xlsx")
OGB2019 <- read.xlsx("P_leucotricha_THeavenpOGB2019_annotations_master.xlsx", sheet = 1)

data <- DRCT72020
data <- DRCT72021
data <- scott2020
data <- OGB2019
data <- p112020
data <- OGB2021

shapiro_test_cleaned <- data %>%
  filter(te_group == "Any TE") %>%
  filter(BUSCO_blast_match == "1") %>%
  filter(IG_five_prime != 99999) %>%  # Remove end of contig values
  filter(!is.na(IG_five_prime)) %>%  
  pull(IG_five_prime)

if (length(shapiro_test_cleaned) >= 3) {
  shapiro_result <- shapiro.test(shapiro_test_cleaned)
  print(shapiro_result)
} else {
  print("Sample size is too small for the Shapiro-Wilk test.")
}

shapiro_test_cleaned <- data %>%
  filter(te_group == "Any TE") %>%
  filter(BUSCO_blast_match == "1") %>%
  filter(IG_three_prime != 99999) %>%  # Remove end of contig values
  filter(!is.na(IG_three_prime)) %>%  
  pull(IG_three_prime)

if (length(shapiro_test_cleaned) >= 3) {
  shapiro_result <- shapiro.test(shapiro_test_cleaned)
  print(shapiro_result)
} else {
  print("Sample size is too small for the Shapiro-Wilk test.")
}
```
Wilcoxon rank sum test - BUSCOs, is 3' larger than 5'
```R
datasets <- list(
  DRCT72020 = DRCT72020,
  DRCT72021 = DRCT72021,
  scott2020 = scott2020,
  OGB2019 = OGB2019,
  p112020 = p112020,
  OGB2021 = OGB2021
)

for (dataset_name in names(datasets)) {
  data <- datasets[[dataset_name]]
  
  three_prime_data <- data %>%
    filter(te_group == "Gene") %>%
    filter(BUSCO_blast_match == "1") %>%
    filter(!is.na(IG_three_prime), IG_three_prime != 99999) %>%
    select(IG_three_prime) %>%
    mutate(group = "three_prime")
  
  five_prime_data <- data %>%
    filter(te_group == "Gene") %>%
    filter(BUSCO_blast_match == "1") %>%
    filter(!is.na(IG_five_prime), IG_five_prime != 99999) %>%
    select(IG_five_prime) %>%
    mutate(group = "five_prime")
  
  combined_data <- bind_rows(
    three_prime_data %>% rename(value = IG_three_prime),
    five_prime_data %>% rename(value = IG_five_prime)
  )

  combined_data$group <- factor(combined_data$group, levels = c("three_prime", "five_prime"))
  
  mann_whitney_result <- wilcox.test(value ~ group, 
                                     data = combined_data, 
                                     alternative = "greater")
  
  cat("\nDataset:", dataset_name, "\n")
  print(mann_whitney_result)
}
```
Dataset: DRCT72020 
W = 2769357, p-value = 2.813e-11
alternative hypothesis: true location shift is less than 0

Dataset: DRCT72021 
W = 3335718, p-value = 1.122e-10
alternative hypothesis: true location shift is less than 0

Dataset: scott2020 
W = 2080894, p-value = 8.793e-16
alternative hypothesis: true location shift is less than 0

Dataset: OGB2019 
W = 4545421, p-value = 8.949e-06
alternative hypothesis: true location shift is less than 0

Dataset: p112020 
W = 4093420, p-value = 2.799e-06
alternative hypothesis: true location shift is less than 0

Dataset: OGB2021 
W = 4015754, p-value = 1.191e-07
alternative hypothesis: true location shift is less than 0

Wilcoxon rank sum test - CSEPs, is 5' larger than 3'
```R
datasets <- list(
  DRCT72020 = DRCT72020,
  DRCT72021 = DRCT72021,
  scott2020 = scott2020,
  OGB2019 = OGB2019,
  p112020 = p112020,
  OGB2021 = OGB2021
)

for (dataset_name in names(datasets)) {
  data <- datasets[[dataset_name]]
  
  three_prime_data <- data %>%
    filter(te_group == "Gene") %>%
    filter(outgroup_orthology_match == 0) %>%
    filter(is_secreted == 1) %>%
    filter(effectorp3_noneffector == ".") %>%
    filter(!is.na(IG_three_prime), IG_three_prime != 99999) %>%
    select(IG_three_prime) %>%
    mutate(group = "three_prime")
  
  five_prime_data <- data %>%
    filter(te_group == "Gene") %>%
    filter(outgroup_orthology_match == 0) %>%
    filter(is_secreted == 1) %>%
    filter(effectorp3_noneffector == ".") %>%
    filter(!is.na(IG_five_prime), IG_five_prime != 99999) %>%
    select(IG_five_prime) %>%
    mutate(group = "five_prime")
  
  combined_data <- bind_rows(
    three_prime_data %>% rename(value = IG_three_prime),
    five_prime_data %>% rename(value = IG_five_prime)
  )

  combined_data$group <- factor(combined_data$group, levels = c("five_prime", "three_prime"))
  
  mann_whitney_result <- wilcox.test(value ~ group, 
                                     data = combined_data, 
                                     alternative = "greater")
  
  cat("\nDataset:", dataset_name, "\n")
  print(mann_whitney_result)
}
```
Dataset: DRCT72020 
W = 7014.5, p-value = 0.04858
alternative hypothesis: true location shift is greater than 0

Dataset: DRCT72021 
W = 6624, p-value = 0.5341
alternative hypothesis: true location shift is greater than 0

Dataset: scott2020 
W = 2663.5, p-value = 0.6066
alternative hypothesis: true location shift is greater than 0

Dataset: OGB2019 
W = 28646, p-value = 0.1317
alternative hypothesis: true location shift is greater than 0

Dataset: p112020 
W = 25213, p-value = 0.431
alternative hypothesis: true location shift is greater than 0

Dataset: OGB2021 
W = 26229, p-value = 0.4327
alternative hypothesis: true location shift is greater than 0

Wilcoxon rank sum test - CAZYs, is 5' larger than 3'
```R
datasets <- list(
  DRCT72020 = DRCT72020,
  DRCT72021 = DRCT72021,
  scott2020 = scott2020,
  OGB2019 = OGB2019,
  p112020 = p112020,
  OGB2021 = OGB2021
)

for (dataset_name in names(datasets)) {
  data <- datasets[[dataset_name]]
  
  three_prime_data <- data %>%
    filter(te_group == "Gene") %>%
    filter(dbcan_matches != ".") %>%
    filter(!is.na(IG_three_prime), IG_three_prime != 99999) %>%
    select(IG_three_prime) %>%
    mutate(group = "three_prime")
  
  five_prime_data <- data %>%
    filter(te_group == "Gene") %>%
    filter(dbcan_matches != ".") %>%
    filter(!is.na(IG_five_prime), IG_five_prime != 99999) %>%
    select(IG_five_prime) %>%
    mutate(group = "five_prime")
  
  combined_data <- bind_rows(
    three_prime_data %>% rename(value = IG_three_prime),
    five_prime_data %>% rename(value = IG_five_prime)
  )

  combined_data$group <- factor(combined_data$group, levels = c("five_prime", "three_prime"))
  
  mann_whitney_result <- wilcox.test(value ~ group, 
                                     data = combined_data, 
                                     alternative = "greater")
  
  cat("\nDataset:", dataset_name, "\n")
  print(mann_whitney_result)
}
```
Dataset: DRCT72020 
W = 6666, p-value = 0.3778
alternative hypothesis: true location shift is greater than 0

Dataset: DRCT72021 
W = 18316, p-value = 0.01473
alternative hypothesis: true location shift is greater than 0

Dataset: scott2020 
W = 4974, p-value = 0.4667
alternative hypothesis: true location shift is greater than 0

Dataset: OGB2019 
W = 20066, p-value = 0.002017
alternative hypothesis: true location shift is greater than 0

Dataset: p112020 
W = 13004, p-value = 0.02906
alternative hypothesis: true location shift is greater than 0

Dataset: OGB2021 
W = 12370, p-value = 0.04242
alternative hypothesis: true location shift is greater than 0

Wilcoxon rank sum test - RALPHs, is 5' larger than 3'
```R
datasets <- list(
  DRCT72020 = DRCT72020,
  DRCT72021 = DRCT72021,
  scott2020 = scott2020,
  OGB2019 = OGB2019,
  p112020 = p112020,
  OGB2021 = OGB2021
)

for (dataset_name in names(datasets)) {
  data <- datasets[[dataset_name]]
  
  three_prime_data <- data %>%
    filter(te_group == "Gene") %>%
    filter(str_detect(effector_matches, "BghBEC1011|BgtAVRa10|BgAVRA13|BgtAvrPm2|BgtSvrPm3a1f1")) %>%
    filter(!is.na(IG_three_prime), IG_three_prime != 99999) %>%
    select(IG_three_prime) %>%
    mutate(group = "three_prime")
  
  five_prime_data <- data %>%
    filter(te_group == "Gene") %>%
    filter(str_detect(effector_matches, "BghBEC1011|BgtAVRa10|BgAVRA13|BgtAvrPm2|BgtSvrPm3a1f1")) %>%
    filter(!is.na(IG_five_prime), IG_five_prime != 99999) %>%
    select(IG_five_prime) %>%
    mutate(group = "five_prime")
  
  combined_data <- bind_rows(
    three_prime_data %>% rename(value = IG_three_prime),
    five_prime_data %>% rename(value = IG_five_prime)
  )

  combined_data$group <- factor(combined_data$group, levels = c("five_prime", "three_prime"))
  
  mann_whitney_result <- wilcox.test(value ~ group, 
                                     data = combined_data, 
                                     alternative = "greater")
  
  cat("\nDataset:", dataset_name, "\n")
  print(mann_whitney_result)
}
```
Dataset: DRCT72020 
W = 5007.5, p-value = 0.006883
alternative hypothesis: true location shift is greater than 0

Dataset: DRCT72021 
W = 3506.5, p-value = 0.006431
alternative hypothesis: true location shift is greater than 0

Dataset: scott2020 
W = 2222.5, p-value = 0.04384
alternative hypothesis: true location shift is greater than 0

Dataset: OGB2019 
W = 52128, p-value = 0.01667
alternative hypothesis: true location shift is greater than 0

Dataset: p112020 
W = 53784, p-value = 0.01695
alternative hypothesis: true location shift is greater than 0

Dataset: OGB2021 
W = 49915, p-value = 0.01254
alternative hypothesis: true location shift is greater than 0

Wilcoxon rank sum test - EKAs, is 5' larger than 3'
```R
datasets <- list(
  DRCT72020 = DRCT72020,
  DRCT72021 = DRCT72021,
  scott2020 = scott2020,
  OGB2019 = OGB2019,
  p112020 = p112020,
  OGB2021 = OGB2021
)

for (dataset_name in names(datasets)) {
  data <- datasets[[dataset_name]]
  
  three_prime_data <- data %>%
    filter(te_group == "Gene") %>%
    filter(str_detect(effector_matches, "BgtAVRk1|BgtAVRa10")) %>%
    filter(is_secreted != 1) %>%
    filter(!is.na(IG_three_prime), IG_three_prime != 99999) %>%
    select(IG_three_prime) %>%
    mutate(group = "three_prime")
  
  five_prime_data <- data %>%
    filter(te_group == "Gene") %>%
    filter(str_detect(effector_matches, "BgtAVRk1|BgtAVRa10")) %>%
    filter(is_secreted != 1) %>%
    filter(!is.na(IG_five_prime), IG_five_prime != 99999) %>%
    select(IG_five_prime) %>%
    mutate(group = "five_prime")
  
  combined_data <- bind_rows(
    three_prime_data %>% rename(value = IG_three_prime),
    five_prime_data %>% rename(value = IG_five_prime)
  )

  combined_data$group <- factor(combined_data$group, levels = c("five_prime", "three_prime"))
  
  mann_whitney_result <- wilcox.test(value ~ group, 
                                     data = combined_data, 
                                     alternative = "greater")
  
  cat("\nDataset:", dataset_name, "\n")
  print(mann_whitney_result)
}
```
Dataset: DRCT72020 
W = 2102.5, p-value = 0.005381
alternative hypothesis: true location shift is greater than 0

Dataset: DRCT72021 
W = 1499.5, p-value = 0.02297
alternative hypothesis: true location shift is greater than 0

Dataset: scott2020 
W = 1226.5, p-value = 0.008693
alternative hypothesis: true location shift is greater than 0

Dataset: OGB2019 
W = 6506.5, p-value = 0.00451
alternative hypothesis: true location shift is greater than 0

Dataset: p112020 
W = 7486, p-value = 0.0009135
alternative hypothesis: true location shift is greater than 0

Dataset: OGB2021 
W = 5897.5, p-value = 0.001609
alternative hypothesis: true location shift is greater than 0