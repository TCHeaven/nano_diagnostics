```bash
Orthogroups=formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.tsv

cat $Orthogroups | grep -e 'AP19' | grep -v -e 'AP20' | grep -v -e 'AP21'| wc -l
cat $Orthogroups | grep -e 'AP19' | grep -e 'AP20' | grep -v -e 'AP21'| wc -l
cat $Orthogroups | grep -v -e 'AP19' | grep -e 'AP20' | grep -v -e 'AP21'| wc -l
cat $Orthogroups | grep -e 'AP19' | grep -v -e 'AP20' | grep -e 'AP21'| wc -l
cat $Orthogroups | grep -e 'AP19' | grep -e 'AP20' | grep -e 'AP21'| wc -l
cat $Orthogroups | grep -v -e 'AP19' | grep -e 'AP20' | grep -e 'AP21'| wc -l
cat $Orthogroups | grep -v -e 'AP19' | grep -v -e 'AP20' | grep -e 'AP21'| wc -l

cat $Orthogroups | grep -e 'RA20' | grep -v -e 'ST20' | grep -v -e 'ST21' | wc -l
cat $Orthogroups | grep -e 'RA20' | grep -e 'ST20' | grep -v -e 'ST21' | wc -l
cat $Orthogroups | grep -v -e 'RA20' | grep -e 'ST20' | grep -v -e 'ST21' | wc -l
cat $Orthogroups | grep -e 'RA20' | grep -v -e 'ST20' | grep -e 'ST21' | wc -l
cat $Orthogroups | grep -e 'RA20' | grep -e 'ST20' | grep -e 'ST21' | wc -l
cat $Orthogroups | grep -v -e 'RA20' | grep -e 'ST20' | grep -e 'ST21' | wc -l
cat $Orthogroups | grep -v -e 'RA20' | grep -v -e 'ST20' | grep -e 'ST21' | wc -l

Orthogroups=formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.txt
cat $Orthogroups | grep -e 'RA20' | grep -e 'ST20' | grep -e 'ST21' | grep -v -e 'AP19' | grep -v -e 'AP20' | grep -v -e 'AP21'| awk '{print $1}' | sed 's@:@@g' > txt-api.txt
cat $Orthogroups | grep -e 'RA20' | grep -e 'ST20' | grep -e 'ST21' | grep -e 'AP19' | grep -e 'AP20' | grep -e 'AP21'| awk '{print $1}' | sed 's@:@@g' > txt-both.txt
cat $Orthogroups | grep -v -e 'RA20' | grep -v -e 'ST20' | grep -v -e 'ST21' | grep -e 'AP19' | grep -e 'AP20' | grep -e 'AP21'| awk '{print $1}' | sed 's@:@@g' > txt-leu.txt

Orthogroups=formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.tsv
cat $Orthogroups | grep -e 'RA20' | grep -e 'ST20' | grep -e 'ST21' | grep -v -e 'AP19' | grep -v -e 'AP20' | grep -v -e 'AP21'| awk '{print $1}' | sed 's@:@@g' > tsv-api.txt
cat $Orthogroups | grep -e 'RA20' | grep -e 'ST20' | grep -e 'ST21' | grep -e 'AP19' | grep -e 'AP20' | grep -e 'AP21'|  awk '{print $1}' | sed 's@:@@g' > tsv-both.txt
cat $Orthogroups | grep -v -e 'RA20' | grep -v -e 'ST20' | grep -v -e 'ST21' | grep -e 'AP19' | grep -e 'AP20' | grep -e 'AP21'| awk '{print $1}' | sed 's@:@@g' > tsv-leu.txt

grep -F -v -f txt-api.txt tsv-api.txt > apidif.txt
grep -F -v -f txt-leu.txt tsv-leu.txt > leudif.txt
grep -F -v -f txt-both.txt tsv-both.txt > bothdif.txt

grep -F -v -f tsv-api.txt txt-api.txt > apidif.txt2
grep -F -v -f tsv-leu.txt txt-leu.txt > leudif.txt2
grep -F -v -f tsv-both.txt txt-both.txt > bothdif.txt2

grep -f apidif.txt formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.tsv | grep -o -P '(?<=ST20).*?(?=.t1)'| sed -e "s@|@@g" | sed 's/$/.t1/'> aphanis_ST20_genes-missing.txt
grep -f leudif.txt formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.tsv | grep -o -P '(?<=AP20).*?(?=.t1)'| sed -e "s@|@@g" | sed 's/$/.t1/'> leuco_AP20_genes-missing.txt
grep -f bothdif.txt formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.tsv | grep -o -P '(?<=ST20).*?(?=.t1)'| sed -e "s@|@@g" | sed 's/$/.t1/'> both_ST20_genes-missing.txt
grep -f bothdif.txt formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.tsv | grep -o -P '(?<=AP20).*?(?=.t1)'| sed -e "s@|@@g" | sed 's/$/.t1/'> both_AP20_genes-missing.txt

grep -f apidif.txt2 formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.tsv | grep -o -P '(?<=ST20).*?(?=.t1)'| sed -e "s@|@@g" | sed 's/$/.t1/'> aphanis_ST20_genes-mistaken.txt
grep -f leudif.txt2 formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.tsv | grep -o -P '(?<=AP20).*?(?=.t1)'| sed -e "s@|@@g" | sed 's/$/.t1/'> leuco_AP20_genes-mistaken.txt
grep -f bothdif.txt2 formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.tsv | grep -o -P '(?<=ST20).*?(?=.t1)'| sed -e "s@|@@g" | sed 's/$/.t1/'> both_ST20_genes-mistaken.txt
grep -f bothdif.txt2 formatted/OrthoFinder/Results_Jul20_2/Orthogroups/Orthogroups.tsv | grep -o -P '(?<=AP20).*?(?=.t1)'| sed -e "s@|@@g" | sed 's/$/.t1/'> both_AP20_genes-mistaken.txt

grep -f aphanis_ST20_genes-mistaken.txt P_aphanis-THeavenDRCT72020_1-aphanis-ranked4.tsv
grep -f leuco_AP20_genes-mistaken.txt P_leucotricha-THeavenp11_1-leucotricha-ranked4.tsv
grep -f both_ST20_genes-mistaken.txt P_aphanis-THeavenDRCT72020_1-common-ranked4.tsv
grep -f both_AP20_genes-mistaken.txt P_leucotricha-THeavenp11_1-common-ranked4.tsv

cut -f1 P_aphanis-THeavenDRCT72020_1-aphanis-ranked4.tsv | grep -F -v -f - aphanis_ST20_genes-missing.txt
cut -f1 P_leucotricha-THeavenp11_1-leucotricha-ranked4.tsv | grep -F -v -f - leuco_AP20_genes-missing.txt
#g13819.t1 - on the dups list
#g13821.t1 - on the dups list
cut -f1 P_aphanis-THeavenDRCT72020_1-common-ranked4.tsv | grep -F -v -f - both_ST20_genes-missing.txt
cut -f1 P_leucotricha-THeavenp11_1-common-ranked4.tsv | grep -F -v -f - both_AP20_genes-missing.txt
#g12419.t1

cat $Orthogroups | grep -v -e 'RA20' | grep -v -e 'ST20' | grep -v -e 'ST21' | grep -v -e 'AP19' | grep -v -e 'AP20' | grep -e 'AP21'| wc -l
cat $Orthogroups | grep -v -e 'RA20' | grep -v -e 'ST20' | grep -v -e 'ST21' | grep -v -e 'AP19' | grep -e 'AP20' | grep -v -e 'AP21'| wc -l
cat $Orthogroups | grep -v -e 'RA20' | grep -v -e 'ST20' | grep -v -e 'ST21' | grep -e 'AP19' | grep -v -e 'AP20' | grep -v -e 'AP21'| wc -l
cat $Orthogroups | grep -v -e 'RA20' | grep -v -e 'ST20' | grep -e 'ST21' | grep -v -e 'AP19' | grep -v -e 'AP20' | grep -v -e 'AP21'| wc -l
cat $Orthogroups | grep -v -e 'RA20' | grep -e 'ST20' | grep -v -e 'ST21' | grep -v -e 'AP19' | grep -v -e 'AP20' | grep -v -e 'AP21'| wc -l
cat $Orthogroups | grep -e 'RA20' | grep -v -e 'ST20' | grep -v -e 'ST21' | grep -v -e 'AP19' | grep -v -e 'AP20' | grep -v -e 'AP21'| wc -l

```

```R

setwd("E:/R")
install.packages("devtools")
install.packages("optparse")
install.packages("colorspace")
install.packages("VennDiagram")
install.packages("grid")
install.packages("ggVennDiagram")
install.packages("ggvenn")
install.packages("installr")
install.packages("dendextend")
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

orthotabs <- read.table("All6_isolates_orthogroups-09052025.tab")
df1 <- t(orthotabs)
df2 <- data.frame(df1)

df3 <- head(df2, 22071)
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

svg("aphanis_V_leucotricha.svg")
plot3
dev.off()

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
```bash
for file in $(ls /home/theaven/projects/niab/theaven/gene_pred/P_*/*/predector_singularity3/results/final_genes_appended_renamed.pep/P*-ranked.tsv); do
file2=$(dirname $file)/final_genes_appended_renamed.pep-ranked.tsv
ID=missing_$(echo $file | rev | cut -d '/' -f1 | rev)
head -n 1 $file > $ID
awk 'FNR==NR {exclude[$1]; next} !($1 in exclude)' $file $file2 >> $ID
done

echo Species+Assembly Predicted_genes pfamID pfam_virulence PHIbase_effector effector_ortholog canonical_SP noncanonical_SP SSCP CSEP Secreted_effectorP3 Secreted_CAZY newSSCP newCSEP newSecreted_effectorP3 csep_val all csep_val phiboref > missingEffectorPredictionReport.txt
for Annotations in $(ls ~/projects/niab/theaven/gene_pred/P_*/*/codingquarry/rep_modeling/final/*annotations3.tsv); do
 ID=$(basename $Annotations | sed 's@_1_abinitio_annotations3.tsv@@g')
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
 all=$(awk '$15 == "1" || $11 != "." || ($31 == "1" && $27 == "." && $7 == "0") {print $1}' $Annotations| grep -v 'name' | wc -l)
 csep_val=$(awk '$11 != "." || ($31 == "1" && $27 == "." && $7 == "0") {print $1}' $Annotations| grep -v 'name' | wc -l)
 PHIbaseororthologeffector=$(awk '$15 == "1" || $11 != "." || ($31 == "1" && $27 == "." && $7 == "0") {print $1}' $Annotations| grep -v 'name' | wc -l)
 echo $ID $Predictedgenes $pfamID $pfamvirulence $PHIbaseeffector $effectorortholog $canonicalsp $noncanonical $sscp $csep $secretedeffector $secretedcazy $newsscp $newcsep $newsecretedeffector $csep_val $all >> missingEffectorPredictionReport.txt
done

for file in $(ls ~/projects/niab/theaven/gene_pred/P_*/*/codingquarry/rep_modeling/final/*annotations3.tsv); do
head -n 1 $file > duplicated_$(basename $file)
file2=$(ls $(dirname $file)/*_duplicated_genes.txt)
sed 's/$/.t1/' $file2 | grep -Ff - $file >> duplicated_$(basename $file)
done
~/projects/niab/theaven/gene_pred/P*/*/codingquarry/rep_modeling/final/*_duplicated_genes.txt
```
