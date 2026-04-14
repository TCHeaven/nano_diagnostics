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
```R
install.packages("readxl")
install.packages("magrittr")
install.packages("dplyr")
library(dplyr)
library(magrittr)
library(readxl)
library(ggplot2)
library(tidyr)

mydata <- read_excel("EDbySNP_zz.xlsx", sheet = 1, col_names = TRUE)
mydata <- mydata %>%
  rename_with(~ paste0("marker_", .), -root_EstimateED10, -isolate_name)

model <- lm(root_EstimateED10 ~ marker_009.1_102387, data = mydata)
anova_table <- anova(model)
SS_marker <- anova_table["marker_009.1_102387", "Sum Sq"]
SS_resid <- anova_table["Residuals", "Sum Sq"]
TSS <- SS_marker + SS_resid
PVE <- SS_marker / TSS
PVE

model <- lm(root_EstimateED10 ~ marker_005.1_648012, data = mydata)
anova_table <- anova(model)
SS_marker <- anova_table["marker_005.1_648012", "Sum Sq"]
SS_resid <- anova_table["Residuals", "Sum Sq"]
TSS <- SS_marker + SS_resid
PVE <- SS_marker / TSS
PVE

glass_delta_func <- function(data, outcome, group_var, resistant_label, susceptible_label) {
  resistant <- data[[outcome]][data[[group_var]] == resistant_label]
  susceptible <- data[[outcome]][data[[group_var]] == susceptible_label]
  
  (mean(resistant, na.rm = TRUE) - mean(susceptible, na.rm = TRUE)) /
    sd(susceptible, na.rm = TRUE)
}

glass_delta_func(mydata, "root_EstimateED10", "marker_009.1_102387", "0", "1")
glass_delta_func(mydata, "root_EstimateED10", "marker_005.1_648012", "0", "1")

# Cohen's d function
cohens_d <- function(data, outcome, group_var, group1_label, group2_label) {
  # Extract the two groups
  x1 <- data[[outcome]][data[[group_var]] == group1_label]
  x2 <- data[[outcome]][data[[group_var]] == group2_label]
  
  # Sample sizes
  n1 <- length(x1)
  n2 <- length(x2)
  
  # Means
  m1 <- mean(x1, na.rm = TRUE)
  m2 <- mean(x2, na.rm = TRUE)
  
  # Standard deviations
  sd1 <- sd(x1, na.rm = TRUE)
  sd2 <- sd(x2, na.rm = TRUE)
  
  # Pooled SD
  pooled_sd <- sqrt(((n1 - 1)*sd1^2 + (n2 - 1)*sd2^2) / (n1 + n2 - 2))
  
  # Cohen's d
  d <- (m1 - m2) / pooled_sd
  return(d)
}

cohens_d(mydata, "root_EstimateED10", "marker_009.1_102387", "0", "1")
cohens_d(mydata, "root_EstimateED10", "marker_005.1_648012", "0", "1")

mydata$marker_009.1_102387 <- as.factor(mydata$marker_009.1_102387)
mydata$marker_005.1_648012 <- as.factor(mydata$marker_005.1_648012)

model <- aov(root_EstimateED10 ~ marker_009.1_102387 * marker_005.1_648012, data = mydata)
shapiro.test(residuals(model))
summary(model)

var_G <- (0.5177 + 0.3617 + 0.0461) / 3
var_E <- 0.0259
H2 <- var_G / (var_G + var_E)
H2

anova_table <- anova(model)
TSS <- sum(anova_table$`Sum Sq`)
SS_marker1 <- anova_table["marker_009.1_102387", "Sum Sq"]
SS_marker2 <- anova_table["marker_005.1_648012", "Sum Sq"]
SS_interaction <- anova_table["marker_009.1_102387:marker_005.1_648012", "Sum Sq"]
PVE_marker1 <- SS_marker1 / TSS
PVE_marker2 <- SS_marker2 / TSS
PVE_interaction <- SS_interaction / TSS
PVE_total_model <- (SS_marker1 + SS_marker2 + SS_interaction) / TSS
cat("PVE Marker 1:", round(PVE_marker1, 4), "\n")
cat("PVE Marker 2:", round(PVE_marker2, 4), "\n")
cat("PVE Interaction:", round(PVE_interaction, 4), "\n")
cat("PVE Combined:", round(PVE_total_model, 4), "\n")

mydata <- read_excel("EDbySNP_zz.xlsx", sheet = 2, col_names = TRUE)
mydata <- mydata %>%
  rename_with(~ paste0("marker_", .), -root_EstimateED50, -isolate_name)

model <- lm(root_EstimateED50 ~ marker_009.1_102387, data = mydata)
anova_table <- anova(model)
SS_marker <- anova_table["marker_009.1_102387", "Sum Sq"]
SS_resid <- anova_table["Residuals", "Sum Sq"]
TSS <- SS_marker + SS_resid
PVE <- SS_marker / TSS
PVE

model <- lm(root_EstimateED50 ~ marker_005.1_648012, data = mydata)
anova_table <- anova(model)
SS_marker <- anova_table["marker_005.1_648012", "Sum Sq"]
SS_resid <- anova_table["Residuals", "Sum Sq"]
TSS <- SS_marker + SS_resid
PVE <- SS_marker / TSS
PVE

model <- lm(root_EstimateED50 ~ marker_016.1_750503, data = mydata)
anova_table <- anova(model)
SS_marker <- anova_table["marker_016.1_750503", "Sum Sq"]
SS_resid <- anova_table["Residuals", "Sum Sq"]
TSS <- SS_marker + SS_resid
PVE <- SS_marker / TSS
PVE

glass_delta_func(mydata, "root_EstimateED50", "marker_009.1_102387", "0", "1")
glass_delta_func(mydata, "root_EstimateED50", "marker_005.1_648012", "0", "1")
glass_delta_func(mydata, "root_EstimateED50", "marker_016.1_750503", "0", "1")

cohens_d(mydata, "root_EstimateED50", "marker_009.1_102387", "0", "1")
cohens_d(mydata, "root_EstimateED50", "marker_005.1_648012", "0", "1")
cohens_d(mydata, "root_EstimateED50", "marker_016.1_750503", "0", "1")

mydata$marker_009.1_102387 <- as.factor(mydata$marker_009.1_102387)
mydata$marker_005.1_648012 <- as.factor(mydata$marker_005.1_648012)

model <- aov(root_EstimateED50 ~ marker_009.1_102387 * marker_005.1_648012, data = mydata)
shapiro.test(residuals(model))
summary(model)

var_G <- (0.6539 + 1.8998 + 0.0121) / 3
var_E <- 0.0750
H2 <- var_G / (var_G + var_E)
H2

anova_table <- anova(model)
TSS <- sum(anova_table$`Sum Sq`)
SS_marker1 <- anova_table["marker_009.1_102387", "Sum Sq"]
SS_marker2 <- anova_table["marker_005.1_648012", "Sum Sq"]
SS_interaction <- anova_table["marker_009.1_102387:marker_005.1_648012", "Sum Sq"]
PVE_marker1 <- SS_marker1 / TSS
PVE_marker2 <- SS_marker2 / TSS
PVE_interaction <- SS_interaction / TSS
PVE_total_model <- (SS_marker1 + SS_marker2 + SS_interaction) / TSS
cat("PVE Marker 1:", round(PVE_marker1, 4), "\n")
cat("PVE Marker 2:", round(PVE_marker2, 4), "\n")
cat("PVE Interaction:", round(PVE_interaction, 4), "\n")
cat("PVE Combined:", round(PVE_total_model, 4), "\n")

mydata <- read_excel("EDbySNP_zz.xlsx", sheet = 3, col_names = TRUE)
mydata <- mydata %>%
  rename_with(~ paste0("marker_", .), -root_EstimateED90, -isolate_name)

mydata$marker_009.1_102387 <- as.factor(mydata$marker_009.1_102387)
mydata$marker_005.1_648012 <- as.factor(mydata$marker_005.1_648012)

model <- aov(root_EstimateED90 ~ marker_009.1_102387 * marker_005.1_648012, data = mydata)
shapiro.test(residuals(model))
summary(model)

anova_table <- anova(model)
TSS <- sum(anova_table$`Sum Sq`)
SS_marker1 <- anova_table["marker_009.1_102387", "Sum Sq"]
SS_marker2 <- anova_table["marker_005.1_648012", "Sum Sq"]
SS_interaction <- anova_table["marker_009.1_102387:marker_005.1_648012", "Sum Sq"]
PVE_marker1 <- SS_marker1 / TSS
PVE_marker2 <- SS_marker2 / TSS
PVE_interaction <- SS_interaction / TSS
PVE_total_model <- (SS_marker1 + SS_marker2 + SS_interaction) / TSS
cat("PVE Marker 1:", round(PVE_marker1, 4), "\n")
cat("PVE Marker 2:", round(PVE_marker2, 4), "\n")
cat("PVE Interaction:", round(PVE_interaction, 4), "\n")
cat("PVE Combined:", round(PVE_total_model, 4), "\n")

ed10 <- read_excel("EDbySNP_zz.xlsx", sheet = 1)
ed10 <- ed10 %>%
  rename_with(~ paste0("marker_", .), -root_EstimateED10, -isolate_name)
ed10 <- ed10 %>%
  rename(ED = root_EstimateED10)

ed50 <- read_excel("EDbySNP_zz.xlsx", sheet = 2)
ed50 <- ed50 %>%
  rename_with(~ paste0("marker_", .), -root_EstimateED50, -isolate_name)
ed50 <- ed50 %>%
  rename(ED = root_EstimateED50)

ed90 <- read_excel("EDbySNP_zz.xlsx", sheet = 3)
ed90 <- ed90 %>%
  rename_with(~ paste0("marker_", .), -root_EstimateED90, -isolate_name)
ed90 <- ed90 %>%
  rename(ED = root_EstimateED90)

ed10$measurement <- "ED10"
ed50$measurement <- "ED50"
ed90$measurement <- "ED90"

combined_data <- bind_rows(ed10, ed50, ed90)

combined_data <- combined_data %>%
  mutate(marker_group = case_when(
    marker_009.1_102387 == 0 & marker_005.1_648012 == 1 ~ "Neither",
    marker_009.1_102387 == 1 & marker_005.1_648012 == 1 ~ "Marker 009 only",
    marker_009.1_102387 == 0 & marker_005.1_648012 == 0 ~ "Marker 005 only",
    marker_009.1_102387 == 1 & marker_005.1_648012 == 0 ~ "Both markers"
  ))

combined_data$ED <- combined_data$ED^2

ggplot(combined_data, aes(x = "", y = ED, color = measurement)) +
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  facet_grid(rows = vars(measurement), cols = vars(marker_group), scales = "free_y") +
  labs(title = "Effective Dose by Marker Group",
       x = NULL,
       y = "ED value") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)  # <- adds plot boxes
  )




```
```R
setwd("E:/R")

library(readxl)

df1 <- read_excel("P_aphanis_THeavenDRCT72020_annotations_master.xlsx")
df2 <- read_excel("P_aphanis_THeavenDRCT72021_annotations_master.xlsx")
df3 <- read_excel("P_aphanis_THeavenSCOTT2020_annotations_master.xlsx")
df4 <- read_excel("P_leucotricha_THeavenpOGB2019_annotations_master.xlsx")
df5 <- read_excel("P_leucotricha_THeavenp11_annotations_master.xlsx")
df6 <- read_excel("P_leucotricha_THeavenpOGB2021_annotations_master.xlsx")

df <- df1

subset_data <- df[df$te_group == 'Any TE' & df$IG_five_prime != 99999, ]


model <- aov(IG_five_prime ~ BUSCO_blast_match, data = subset_data)
plot(model)                # Residual plots
shapiro.test(residuals(model))  # Test normality of residuals

t.test(IG_five_prime ~ BUSCO_blast_match, data = subset_data)
model <- lm(IG_five_prime ~ BUSCO_blast_match, data = subset_data)
qqnorm(residuals(model))
qqline(residuals(model), col = "red")

shapiro.test(residuals(model))

res <- residuals(model)
shapiro.test(sample(res, 5000))

subset_data <- df[df$te_group == "Any TE" & grepl("[0-9]", df$five_prime_lgth), ]
model <- aov(five_prime_lgth ~ sec_type, data = subset_data)
plot(model)     
shapiro.test(residuals(model))

subset_data <- df[df$te_group == "Any TE" & grepl("[0-9]", df$three_prime_lgth), ]
model <- aov(three_prime_lgth ~ sec_type, data = subset_data)
plot(model) 

subset_data <- df[df$te_group == "LINE" & grepl("[0-9]", df$five_prime_lgth), ]
model <- aov(five_prime_lgth ~ sec_type, data = subset_data)
plot(model)     
shapiro.test(residuals(model))

subset_data <- df[df$te_group == "LINE" & grepl("[0-9]", df$three_prime_lgth), ]
model <- aov(three_prime_lgth ~ sec_type, data = subset_data)
plot(model) 

subset_data <- df[df$te_group == "Tad" & grepl("[0-9]", df$five_prime_lgth), ]
model <- aov(five_prime_lgth ~ sec_type, data = subset_data)
plot(model)     
shapiro.test(residuals(model))

subset_data <- df[df$te_group == "Tad" & grepl("[0-9]", df$three_prime_lgth), ]
model <- aov(three_prime_lgth ~ sec_type, data = subset_data)
plot(model) 

subset_data <- df[df$te_group == "LTR" & grepl("[0-9]", df$five_prime_lgth), ]
model <- aov(five_prime_lgth ~ sec_type, data = subset_data)
plot(model)     
shapiro.test(residuals(model))

subset_data <- df[df$te_group == "LTR" & grepl("[0-9]", df$three_prime_lgth), ]
model <- aov(three_prime_lgth ~ sec_type, data = subset_data)
plot(model) 

subset_data <- df[df$te_group == "Copia" & grepl("[0-9]", df$five_prime_lgth), ]
model <- aov(five_prime_lgth ~ sec_type, data = subset_data)
plot(model)     
shapiro.test(residuals(model))

subset_data <- df[df$te_group == "Copia" & grepl("[0-9]", df$three_prime_lgth), ]
model <- aov(three_prime_lgth ~ sec_type, data = subset_data)
plot(model) 

subset_data <- df[df$te_group == "Gypsy" & grepl("[0-9]", df$five_prime_lgth), ]
model <- aov(five_prime_lgth ~ sec_type, data = subset_data)
plot(model)     
shapiro.test(residuals(model))

subset_data <- df[df$te_group == "Gypsy" & grepl("[0-9]", df$three_prime_lgth), ]
model <- aov(three_prime_lgth ~ sec_type, data = subset_data)
plot(model) 

subset_data <- df[df$te_group == "DNA" & grepl("[0-9]", df$five_prime_lgth), ]
model <- aov(five_prime_lgth ~ sec_type, data = subset_data)
plot(model)     
shapiro.test(residuals(model))

subset_data <- df[df$te_group == "DNA" & grepl("[0-9]", df$three_prime_lgth), ]
model <- aov(three_prime_lgth ~ sec_type, data = subset_data)
plot(model) 

subset_data <- df[df$te_group == "Mariner" & grepl("[0-9]", df$five_prime_lgth), ]
model <- aov(five_prime_lgth ~ sec_type, data = subset_data)
plot(model)     
shapiro.test(residuals(model))

subset_data <- df[df$te_group == "Mariner" & grepl("[0-9]", df$three_prime_lgth), ]
model <- aov(three_prime_lgth ~ sec_type, data = subset_data)
plot(model) 

subset_data <- df[df$te_group == "Pogo" & grepl("[0-9]", df$five_prime_lgth), ]
model <- aov(five_prime_lgth ~ sec_type, data = subset_data)
plot(model)     
shapiro.test(residuals(model))

subset_data <- df[df$te_group == "Pogo" & grepl("[0-9]", df$three_prime_lgth), ]
model <- aov(three_prime_lgth ~ sec_type, data = subset_data)
plot(model) 


```





























































Frozen genomes:
```bash
/home/theaven/projects/niab/theaven/genomes/Pod_leu_OGBp112020_clean_renamed.fasta
/home/theaven/projects/niab/theaven/genomes/Pod_leu_OGB2021-2_clean_renamed_fix.fasta
/home/theaven/projects/niab/theaven/genomes/Pod_leu_OGB2019-2_clean_renamed_fix.fasta
/home/theaven/projects/niab/theaven/genomes/Pod_aph_DRT72021-2_clean_renamed_fix.fasta
/home/theaven/projects/niab/theaven/genomes/Pod_aph_SCOTT2020-2_clean_renamed_fix.fasta
/home/theaven/projects/niab/theaven/genomes/Pod_aph_DRT72020_clean_renamed.fasta
```
Combined gene and TE predictions:
```bash
/home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/te.gff
/home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/te.gff
/home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/te.gff
/home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/te.gff
/home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/te.gff
/home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/te.gff
```
Genes:
```bash
/home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3
/home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3
/home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3
/home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3
/home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3
/home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3
```
Contigs in the frozen genomes:
```bash
grep '>' /home/theaven/projects/niab/theaven/genomes/Pod_leu_OGBp112020_clean_renamed.fasta | cut -d '_' -f4,5 | sed 's/^/^/; s/$/\\b/' > Pod_leu_OGBp112020_clean_renamed.txt
grep '>' /home/theaven/projects/niab/theaven/genomes/Pod_leu_OGB2021-2_clean_renamed_fix.fasta | cut -d '_' -f4,5 | sed 's/^/^/; s/$/\\b/' > Pod_leu_OGB2021-2_clean_renamed_fix.txt
grep '>' /home/theaven/projects/niab/theaven/genomes/Pod_leu_OGB2019-2_clean_renamed_fix.fasta | cut -d '_' -f4,5 | sed 's/^/^/; s/$/\\b/' > Pod_leu_OGB2019-2_clean_renamed_fix.txt
grep '>' /home/theaven/projects/niab/theaven/genomes/Pod_aph_DRT72021-2_clean_renamed_fix.fasta | cut -d '_' -f4,5 | sed 's/^/^/; s/$/\\b/' > Pod_aph_DRT72021-2_clean_renamed_fix.txt
grep '>' /home/theaven/projects/niab/theaven/genomes/Pod_aph_SCOTT2020-2_clean_renamed_fix.fasta | cut -d '_' -f4,5 | sed 's/^/^/; s/$/\\b/' > Pod_aph_SCOTT2020-2_clean_renamed_fix.txt
grep '>' /home/theaven/projects/niab/theaven/genomes/Pod_aph_DRT72020_clean_renamed.fasta | cut -d '_' -f4,5 | sed 's/^/^/; s/$/\\b/' > Pod_aph_DRT72020_clean_renamed.txt
```
Frozen genes gff:
```bash
grep -f Pod_aph_DRT72020_clean_renamed.txt /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3 > /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3
grep -f Pod_aph_DRT72021-2_clean_renamed_fix.txt /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3 > /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3
grep -f Pod_aph_SCOTT2020-2_clean_renamed_fix.txt /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3 > /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3
grep -f Pod_leu_OGBp112020_clean_renamed.txt /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3 > /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3
grep -f Pod_leu_OGB2019-2_clean_renamed_fix.txt /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3 > /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3
grep -f Pod_leu_OGB2021-2_clean_renamed_fix.txt /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/final_genes_appended_renamed.gff3 > /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3

for file in $(ls /home/theaven/projects/niab/theaven/gene_pred/*/*/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3); do
  ID=$(echo $file | cut -d '/' -f9)
  cp $file /home/theaven/scratch/manuscript4/06012025/${ID}_$(basename $file)

done
```
Frozen Combined gene and TE predictions:
```bash
grep -f Pod_aph_DRT72020_clean_renamed.txt /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/te.gff > /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/frozen_te.gff
grep -f Pod_aph_DRT72021-2_clean_renamed_fix.txt /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/te.gff > /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/frozen_te.gff
grep -f Pod_aph_SCOTT2020-2_clean_renamed_fix.txt /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/te.gff > /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/frozen_te.gff
grep -f Pod_leu_OGBp112020_clean_renamed.txt /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/te.gff > /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/frozen_te.gff
grep -f Pod_leu_OGB2019-2_clean_renamed_fix.txt /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/te.gff > /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/frozen_te.gff
grep -f Pod_leu_OGB2021-2_clean_renamed_fix.txt /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/te.gff > /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/frozen_te.gff

for file in $(ls /home/theaven/projects/niab/theaven/gene_pred/*/*/codingquarry/rep_modeling/final/frozen_te.gff); do
  ID=$(echo $file | cut -d '/' -f9)
  cp $file /home/theaven/scratch/manuscript4/06012025/${ID}_$(basename $file)

done
```
Frozen gene names:
```bash
awk -F'\t' '$3 == "gene" {split($9, id_parts, /;|ID=/); for (i in id_parts) { if (id_parts[i] ~ /^[^;]+$/) { print id_parts[i]; break; } } }' /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3 | sed 's/$/.t1/' > Pod_aph_DRT72020_clean_renamed_genes.txt
awk -F'\t' '$3 == "gene" {split($9, id_parts, /;|ID=/); for (i in id_parts) { if (id_parts[i] ~ /^[^;]+$/) { print id_parts[i]; break; } } }' /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3 | sed 's/$/.t1/' > Pod_aph_DRT72021-2_clean_renamed_fix_genes.txt
awk -F'\t' '$3 == "gene" {split($9, id_parts, /;|ID=/); for (i in id_parts) { if (id_parts[i] ~ /^[^;]+$/) { print id_parts[i]; break; } } }' /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3 | sed 's/$/.t1/' > Pod_aph_SCOTT2020-2_clean_renamed_fix_genes.txt
awk -F'\t' '$3 == "gene" {split($9, id_parts, /;|ID=/); for (i in id_parts) { if (id_parts[i] ~ /^[^;]+$/) { print id_parts[i]; break; } } }' /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3 | sed 's/$/.t1/' > Pod_leu_OGBp112020_clean_renamed_genes.txt
awk -F'\t' '$3 == "gene" {split($9, id_parts, /;|ID=/); for (i in id_parts) { if (id_parts[i] ~ /^[^;]+$/) { print id_parts[i]; break; } } }' /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3 | sed 's/$/.t1/' > Pod_leu_OGB2019-2_clean_renamed_fix_genes.txt
awk -F'\t' '$3 == "gene" {split($9, id_parts, /;|ID=/); for (i in id_parts) { if (id_parts[i] ~ /^[^;]+$/) { print id_parts[i]; break; } } }' /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3 | sed 's/$/.t1/' > Pod_leu_OGB2021-2_clean_renamed_fix_genes.txt
```
Frozen gene fastas and annotations
```bash
python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file Pod_leu_OGBp112020_clean_renamed_genes.txt --input /home/theaven/projects/niab/theaven/THeavenp11_1.pep.fasta --output /home/theaven/projects/niab/theaven/Pod_leu_OGBp112020_clean_renamed.pep.fasta
python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file Pod_leu_OGB2019-2_clean_renamed_fix_genes.txt --input /home/theaven/projects/niab/theaven/THeavenpOGB2019_1.pep.fasta --output /home/theaven/projects/niab/theaven/Pod_leu_OGB2019-2_clean_renamed.pep.fasta
python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file Pod_leu_OGB2021-2_clean_renamed_fix_genes.txt --input /home/theaven/projects/niab/theaven/THeavenpOGB2021_1.pep.fasta --output /home/theaven/projects/niab/theaven/Pod_leu_OGB2021-2_clean_renamed.pep.fasta
python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file Pod_aph_SCOTT2020-2_clean_renamed_fix_genes.txt --input /home/theaven/projects/niab/theaven/THeavenSCOTT2020_1.pep.fasta --output /home/theaven/projects/niab/theaven/Pod_aph_SCOTT2020-2_clean_renamed.pep.fasta
python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file Pod_aph_DRT72021-2_clean_renamed_fix_genes.txt --input /home/theaven/projects/niab/theaven/THeavenDRCT72021_1.pep.fasta --output /home/theaven/projects/niab/theaven/Pod_aph_DRT72021-2_clean_renamed.pep.fasta
python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file Pod_aph_DRT72020_clean_renamed_genes.txt --input /home/theaven/projects/niab/theaven/THeavenDRCT72020_1.pep.fasta --output /home/theaven/projects/niab/theaven/Pod_aph_DRT72020_clean_renamed.pep.fasta

grep -f Pod_leu_OGBp112020_clean_renamed_genes.txt /home/theaven/projects/niab/theaven/THeavenp11_1_annotations3.tsv > /home/theaven/projects/niab/theaven/frozen_THeavenp11_1_annotations3.tsv 
grep -f Pod_leu_OGB2019-2_clean_renamed_fix_genes.txt /home/theaven/projects/niab/theaven/THeavenpOGB2019_1_annotations3.tsv > /home/theaven/projects/niab/theaven/frozen_THeavenpOGB2019_1_annotations3.tsv 
grep -f Pod_leu_OGB2021-2_clean_renamed_fix_genes.txt /home/theaven/projects/niab/theaven/THeavenpOGB2021_1_annotations3.tsv > /home/theaven/projects/niab/theaven/frozen_THeavenpOGB2021_1_annotations3.tsv 
grep -f Pod_aph_SCOTT2020-2_clean_renamed_fix_genes.txt /home/theaven/projects/niab/theaven/THeavenSCOTT2020_1_annotations3.tsv > /home/theaven/projects/niab/theaven/frozen_THeavenSCOTT2020_1_annotations3.tsv 
grep -f Pod_aph_DRT72021-2_clean_renamed_fix_genes.txt /home/theaven/projects/niab/theaven/THeavenDRCT72021_1_annotations3.tsv > /home/theaven/projects/niab/theaven/frozen_THeavenDRCT72021_1_annotations3.tsv
grep -f Pod_aph_DRT72020_clean_renamed_genes.txt /home/theaven/projects/niab/theaven/THeavenDRCT72020_1_annotations3.tsv > /home/theaven/projects/niab/theaven/frozen_THeavenDRCT72020_1_annotations3.tsv
```
```bash
echo Species+Assembly Predicted_genes pfamID SP CSEP newCSEP PHIbase_effector effector_ortholog EKA RALPH All > SixEffectorPredictionReport.txt
for Annotations in $(ls /home/theaven/projects/niab/theaven/frozen_*_annotations3.tsv); do
 ID=$(basename $Annotations | cut -d '_' -f1,2,3)
 Predictedgenes=$(cut -f1 $Annotations | grep -v 'name' | wc -l)
 pfamID=$(awk -F'\t' '$27 != "." {print $1}' $Annotations | grep -v 'name' | wc -l)
 SPs=$(awk -F'\t' '$40 == "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 csep=$(awk -F'\t' '$40 == "1" && $16 == "0" && $36 == "." {print $1}' $Annotations | grep -v 'name' | wc -l)
 newcsep=$(awk -F'\t' '$40 == "1" && $16 == "0" && $36 == "." && $20 == "." && $24 != "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 PHIbaseeffector=$(awk -F'\t' '$24 == "1" {print $1}' $Annotations| grep -v 'name' | wc -l)
 effectorortholog=$(awk -F'\t' '$20 != "." {print $1}' $Annotations| grep -v 'name' | wc -l)
 EKA=$(awk -F'\t' '$20 ~ /BgtAVRa10|BgtAVRk1/ {print $1}' $Annotations| grep -v 'name' | wc -l)
 RALPH=$(awk -F'\t' '$20 ~ /BghBEC1011|BgtAVRa10|BgAVRA13|BgtAvrPm2|BgtSvrPm3a1f1/ {print $1}' $Annotations| grep -v 'name' | wc -l)
 all=$(awk -F'\t' '$24 == "1" || $20 != "." || ($40 == "1" && $36 == "." && $16 == "0") {print $1}' $Annotations| grep -v 'name' | wc -l)
 echo $ID $Predictedgenes $pfamID $SPs $csep $newcsep $PHIbaseeffector $effectorortholog $EKA $RALPH $all >> SixEffectorPredictionReport.txt
done
```
Frozen all mildew annotations:
```bash
for file in $(ls Pod_leu_OGBp112020_clean_renamed_genes.txt Pod_leu_OGB2019-2_clean_renamed_fix_genes.txt Pod_leu_OGB2021-2_clean_renamed_fix_genes.txt Pod_aph_SCOTT2020-2_clean_renamed_fix_genes.txt Pod_aph_DRT72021-2_clean_renamed_fix_genes.txt Pod_aph_DRT72020_clean_renamed_genes.txt); do
  ID=$(echo $file | cut -d '_' -f1,2,3 | sed 's@-2@@g')
  search2=$(ls ${ID}*_clean_renamed*.txt | grep -v 'genes')

grep -f $file /home/theaven/scratch/manuscript4/06012025/${ID}_ab_initio_annotations3.tsv > temp.tsv && mv temp.tsv /home/theaven/scratch/manuscript4/06012025/${ID}_ab_initio_annotations3.tsv

python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file $file --input /home/theaven/scratch/manuscript4/06012025/${ID}_ab_initio_cseps-xtra.fasta --output temp.fasta && mv temp.fasta /home/theaven/scratch/manuscript4/06012025/${ID}_ab_initio_cseps-xtra.fasta

grep -f $search2 /home/theaven/scratch/manuscript4/06012025/${ID}_ab_initio_final_genes_renamed.gff3 > temp.gff3 && mv temp.gff3 /home/theaven/scratch/manuscript4/06012025/${ID}_ab_initio_final_genes_renamed.gff3

grep -f $file /home/theaven/scratch/manuscript4/06012025/${ID}_annotations3.tsv > temp.tsv && mv temp.tsv /home/theaven/scratch/manuscript4/06012025/${ID}_annotations3.tsv

python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file $file --input /home/theaven/scratch/manuscript4/06012025/${ID}_clean_ab_initio_final_genes_renamed.pep.fasta --output temp.fasta && mv temp.fasta /home/theaven/scratch/manuscript4/06012025/${ID}_clean_ab_initio_final_genes_renamed.pep.fasta

python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file $file --input /home/theaven/scratch/manuscript4/06012025/${ID}_clean_final_genes_renamed.pep.fasta --output temp.fasta && mv  temp.fasta /home/theaven/scratch/manuscript4/06012025/${ID}_clean_final_genes_renamed.pep.fasta

python3 ~/git_repos/Scripts/NBI/seq_get.py --id_file $file --input /home/theaven/scratch/manuscript4/06012025/${ID}_cseps-xtra.fasta --output temp.fasta && mv  temp.fasta /home/theaven/scratch/manuscript4/06012025/${ID}_cseps-xtra.fasta

grep -f $search2 /home/theaven/scratch/manuscript4/06012025/${ID}_final_genes_renamed.gff3 > temp.gff3 && mv temp.gff3 /home/theaven/scratch/manuscript4/06012025/${ID}_final_genes_renamed.gff3

done
```
```bash
echo Species+Assembly Predicted_genes pfamID pfam_virulence PHIbase_effector effector_ortholog canonical_SP noncanonical_SP SSCP CSEP Secreted_effectorP3 Secreted_CAZY newSSCP newCSEP newSecreted_effectorP3 csep_val all csep_val phiboref > AbinitioEffectorPredictionReport.txt
for Annotations in $(ls /mnt/shared/scratch/theaven/manuscript4/06012025/*annotations3.tsv); do
 ID=$(basename $Annotations | cut -d '_' -f1,2,3)
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
```
```bash
conda activate funannotate

funannotate remote -i /home/theaven/projects/niab/theaven/Pod_aph_DRT72020_clean_renamed.pep.fasta -m interproscan -e tcheaven@googlemail.com
funannotate remote -i /home/theaven/projects/niab/theaven/Pod_aph_DRT72021-2_clean_renamed.pep.fasta -m interproscan -e tcheaven@googlemail.com
funannotate remote -i /home/theaven/projects/niab/theaven/Pod_aph_SCOTT2020-2_clean_renamed.pep.fasta -m interproscan -e tcheaven@googlemail.com
funannotate remote -i /home/theaven/projects/niab/theaven/Pod_leu_OGB2019-2_clean_renamed.pep.fasta -m interproscan -e tcheaven@googlemail.com
funannotate remote -i /home/theaven/projects/niab/theaven/Pod_leu_OGB2021-2_clean_renamed.pep.fasta -m interproscan -e tcheaven@googlemail.com
funannotate remote -i /home/theaven/projects/niab/theaven/Pod_leu_OGBp112020_clean_renamed.pep.fasta -m interproscan -e tcheaven@googlemail.com

cd ~/scratch/manuscript4/my_interproscan/interproscan-5.76-107.0
srun  -p long  -c 16 --mem=64G --pty bash
for file in $(ls /home/theaven/projects/niab/theaven/*_clean_renamed.pep.fasta); do
sbatch ~/git_repos/Wrappers/gruffalo/interproscan.sh "$file"
done #8432718-23

./interproscan.sh \
  -cpu 16 \
  -i "$file" \
  -f xml \
  -goterms \
  -iprlookup \
  -dp \
  -o "$file".interproscan.xml
done

funannotate species -d /home/theaven/scratch/manuscript4/db/funannotate
funannotate setup -b fungi -d /home/theaven/scratch/manuscript4/db/funannotate

cd ~/scratch/manuscript4
srun  -p long  -c 16 --mem=64G --pty bash

for FA in $(ls /home/theaven/projects/niab/theaven/genomes/Pod_leu_OGB2021-2_clean_renamed_fix.fasta /home/theaven/projects/niab/theaven/genomes/Pod_leu_OGB2019-2_clean_renamed_fix.fasta /home/theaven/projects/niab/theaven/genomes/Pod_leu_OGBp112020_clean_renamed.fasta /home/theaven/projects/niab/theaven/genomes/Pod_aph_SCOTT2020-2_clean_renamed_fix.fasta /home/theaven/projects/niab/theaven/genomes/Pod_aph_DRT72021-2_clean_renamed_fix.fasta /home/theaven/projects/niab/theaven/genomes/Pod_aph_DRT72020_clean_renamed.fasta); do
FA_OUT=$(echo $FA | sed 's@.fasta@_renamed2.fasta@g')

awk '
  /^>/ {
    match($0, /(contig_[^[:space:]]+)/, a)
    if (a[1] != "") {
      print ">" a[1]
    } else {
      print $0
    }
    next
  }
  { print }
' "$FA" > "$FA_OUT"
echo "$FA_OUT"
done

conda activate signalp6
for In in $(ls /home/theaven/projects/niab/theaven/*_clean_renamed.pep.fasta); do
  ID=$(basename $In | cut -d '_' -f1,2,3)
  Out=$(dirname $In)/signalp6/$ID
sbatch ~/git_repos/Wrappers/gruffalo/signalp6-fast.sh $In $Out
done #8612838-43
conda deactivate

for In in $(ls /home/theaven/projects/niab/theaven/*_clean_renamed.pep.fasta); do
  ID=$(basename $In | cut -d '_' -f1,2,3)
    Out=$(dirname $In)/signalp6/$ID
awk 'BEGIN{OFS="\t"} /^#/ {next} {yn = ($2 ~ /^SP/) ? "Y" : "N"; cs = "."; if (yn=="Y" && match($0,/CS pos: ([0-9]+-[0-9]+)/,a)) cs=a[1]; print $1,0,0,0,0,0,0,0,0,yn,cs}' "$Out"/prediction_results.txt > "$Out"/signalp6_shortform.tsv
done

conda activate diamond
for In in $(ls /home/theaven/projects/niab/theaven/*_clean_renamed.pep.fasta); do
  ID=$(basename $In | cut -d '_' -f1,2,3)
  Out=~/scratch/manuscript4/swissprot/$ID
  mkdir -p "$Out"
diamond blastp \
  -d /home/theaven/scratch/manuscript4/db/diamond/uniprot_sprot \
  -q "$In" \
  -o "$Out"/"$ID"_vs_swissprot.diamond.tsv \
  --max-target-seqs 1 \
  --evalue 1e-5 \
  --threads 16 \
  --outfmt 6
done

for In in $(ls /home/theaven/projects/niab/theaven/*_clean_renamed.pep.fasta); do
  ID=$(basename $In | cut -d '_' -f1,2,3)
  Out=~/scratch/manuscript4/eggnog/$ID
  mkdir -p "$Out"
  /home/theaven/eggnog-mapper/emapper.py \
  -i "$In" \
  --data_dir /home/theaven/eggnog-mapper \
  --itype proteins \
  --cpu 16 \
  -o "$ID" \
  --output_dir "$Out" --override
done

emapper.py -i proteins.faa --data_dir /home/theaven/eggnog-mapper -o annotations


export FUNANNOTATE_DB=/home/theaven/scratch/manuscript4/db/funannotate

funannotate annotate \
  --gff /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3 \
  --fasta /home/theaven/projects/niab/theaven/genomes/Pod_leu_OGB2021-2_clean_renamed_fix_renamed2.fasta \
  --species "Podosphaera leucotricha" --strain OGB2021 \
  --iprscan /home/theaven/projects/niab/theaven/Pod_leu_OGB2021-2_clean_renamed.pep.fasta.interproscan.xml \
  --eggnog /home/theaven/scratch/manuscript4/eggnog/Pod_leu_OGB2021-2/Pod_leu_OGB2021-2.emapper.annotations \
  --busco_db ascomycota \
  --cpus 16 \
  -d /home/theaven/scratch/manuscript4/db/funannotate \
  -o /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/funannotate3 \
  --sbt /home/theaven/scratch/manuscript4/OGB2021_template.sbt \
  --rename MP843_ --force

funannotate annotate \
  --gff /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3 \
  --fasta /home/theaven/projects/niab/theaven/genomes/Pod_leu_OGB2019-2_clean_renamed_fix_renamed2.fasta \
  --species "Podosphaera leucotricha" --strain OGB2019 \
  --iprscan /home/theaven/projects/niab/theaven/Pod_leu_OGB2019-2_clean_renamed.pep.fasta.interproscan.xml \
  --eggnog /home/theaven/scratch/manuscript4/eggnog/Pod_leu_OGB2019-2/Pod_leu_OGB2019-2.emapper.annotations \
  --busco_db ascomycota \
  --cpus 16 \
  -d /home/theaven/scratch/manuscript4/db/funannotate \
  -o /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/funannotate3 \
  --sbt /home/theaven/scratch/manuscript4/OGB2019_template.sbt \
  --rename MP842_ --force

funannotate annotate \
  --gff /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3 \
  --fasta /home/theaven/projects/niab/theaven/genomes/Pod_leu_OGBp112020_clean_renamed_renamed2.fasta \
  --species "Podosphaera leucotricha" --strain P112020 \
  --iprscan /home/theaven/projects/niab/theaven/Pod_leu_OGBp112020_clean_renamed.pep.fasta.interproscan.xml \
  --eggnog /home/theaven/scratch/manuscript4/eggnog/Pod_leu_OGBp112020/Pod_leu_OGBp112020.emapper.annotations \
  --busco_db ascomycota \
  --cpus 16 \
  -d /home/theaven/scratch/manuscript4/db/funannotate \
  -o /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/funannotate3 \
  --sbt /home/theaven/scratch/manuscript4/P112020_template.sbt \
  --rename K3496_ --force

funannotate annotate \
  --gff /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3 \
  --fasta /home/theaven/projects/niab/theaven/genomes/Pod_aph_SCOTT2020-2_clean_renamed_fix_renamed2.fasta \
  --species "Podosphaera aphanis" --strain SCOTT2021 \
  --iprscan /home/theaven/projects/niab/theaven/Pod_aph_SCOTT2020-2_clean_renamed.pep.fasta.interproscan.xml \
  --eggnog /home/theaven/scratch/manuscript4/eggnog/Pod_aph_SCOTT2020-2/Pod_aph_SCOTT2020-2.emapper.annotations \
  --busco_db ascomycota \
  --cpus 16 \
  -d /home/theaven/scratch/manuscript4/db/funannotate \
  -o /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/funannotate3 \
  --sbt /home/theaven/scratch/manuscript4/SCOTT2020_template.sbt \
  --rename LBW10_ --force

funannotate annotate \
  --gff /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3 \
  --fasta /home/theaven/projects/niab/theaven/genomes/Pod_aph_DRT72021-2_clean_renamed_fix_renamed2.fasta \
  --species "Podosphaera aphanis" --strain DRCT72021 \
  --iprscan /home/theaven/projects/niab/theaven/Pod_aph_DRT72021-2_clean_renamed.pep.fasta.interproscan.xml \
  --eggnog /home/theaven/scratch/manuscript4/eggnog/Pod_aph_DRT72021-2/Pod_aph_DRT72021-2.emapper.annotations \
  --busco_db ascomycota \
  --cpus 16 \
  -d /home/theaven/scratch/manuscript4/db/funannotate \
  -o /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/funannotate3 \
  --sbt /home/theaven/scratch/manuscript4/DRCT2021_template.sbt \
  --rename MP841_ --force

funannotate annotate \
  --gff /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3 \
  --fasta /home/theaven/projects/niab/theaven/genomes/Pod_aph_DRT72020_clean_renamed_renamed2.fasta \
  --species "Podosphaera aphanis" --strain DRCT72020 \
  --iprscan /home/theaven/projects/niab/theaven/Pod_aph_DRT72020_clean_renamed.pep.fasta.interproscan.xml \
  --eggnog /home/theaven/scratch/manuscript4/eggnog/Pod_aph_DRT72020/Pod_aph_DRT72020.emapper.annotations \
  --busco_db ascomycota \
  --cpus 16 \
  -d /home/theaven/scratch/manuscript4/db/funannotate \
  -o /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/funannotate3 \
  --sbt /home/theaven/scratch/manuscript4/DRCT72020_template.sbt \
  --rename K3495_ --force

TBL=/home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_DRCT72020.tbl
egg=/home/theaven/scratch/manuscript4/eggnog/Pod_aph_DRT72020/Pod_aph_DRT72020.emapper.annotations
singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/python3.sif python3 ~/git_repos/Scripts/gruffalo/ec_unknown_from_tbl2.py "$TBL" --emapper "$egg" --join-on protein_id -o "$(dirname $TBL)"
TBL=/home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_DRCT72021.tbl
egg=/home/theaven/scratch/manuscript4/eggnog/Pod_aph_DRT72021-2/Pod_aph_DRT72021-2.emapper.annotations
singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/python3.sif python3 ~/git_repos/Scripts/gruffalo/ec_unknown_from_tbl2.py "$TBL" --emapper "$egg" --join-on protein_id -o "$(dirname $TBL)"
TBL=/home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_SCOTT2021.tbl
egg=/home/theaven/scratch/manuscript4/eggnog/Pod_aph_SCOTT2020-2/Pod_aph_SCOTT2020-2.emapper.annotations
singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/python3.sif python3 ~/git_repos/Scripts/gruffalo/ec_unknown_from_tbl2.py "$TBL" --emapper "$egg" --join-on protein_id -o "$(dirname $TBL)"
TBL=/home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_leucotricha_OGB2019.tbl
egg=/home/theaven/scratch/manuscript4/eggnog/Pod_leu_OGB2019-2/Pod_leu_OGB2019-2.emapper.annotations
singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/python3.sif python3 ~/git_repos/Scripts/gruffalo/ec_unknown_from_tbl2.py "$TBL" --emapper "$egg" --join-on protein_id -o "$(dirname $TBL)"
TBL=/home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_leucotricha_OGB2021.tbl
egg=/home/theaven/scratch/manuscript4/eggnog/Pod_leu_OGB2021-2/Pod_leu_OGB2021-2.emapper.annotations
singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/python3.sif python3 ~/git_repos/Scripts/gruffalo/ec_unknown_from_tbl2.py "$TBL" --emapper "$egg" --join-on protein_id -o "$(dirname $TBL)"
TBL=/home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_leucotricha_P112020.tbl
egg=/home/theaven/scratch/manuscript4/eggnog/Pod_leu_OGBp112020/Pod_leu_OGBp112020.emapper.annotations
singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/python3.sif python3 ~/git_repos/Scripts/gruffalo/ec_unknown_from_tbl2.py "$TBL" --emapper "$egg" --join-on protein_id -o "$(dirname $TBL)"

for TBL in $(ls /home/theaven/projects/niab/theaven/gene_pred/*/*/codingquarry/rep_modeling/final/funannotate3/annotate_results/*.tbl | grep -v 'part' | grep -v 'fix'); do
  Fun_Dir=$(dirname "$TBL")
singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/python3.sif python3 ~/git_repos/Scripts/gruffalo/remove_ec_from_problem_tbl.py \
  --tbl "$TBL" \
  --hits "$(ls $Fun_Dir/*_EC_UNKNOWN_with_eggnog.tsv)" \
  --out "$Fun_Dir"/"$(basename $TBL | sed 's@.tbl@_fix.tbl@g')"
done


tbl2asn \
  -t /home/theaven/scratch/manuscript4/DRCT2021_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_DRCT72021_fix.tbl \
  -a s -V vb \
  -Z /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/discrep_report_after.txt \
  -o /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_DRCT72021_tbl2asn.sqn

tbl2asn \
  -t /home/theaven/scratch/manuscript4/DRCT72020_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_DRCT72020_fix.tbl \
  -a s -V vb \
  -Z /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/discrep_report_after.txt \
  -o /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_DRCT72020_tbl2asn.sqn

tbl2asn \
  -t /home/theaven/scratch/manuscript4/OGB2019_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_leucotricha_OGB2019_fix.tbl \
  -a s -V vb \
  -Z /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/discrep_report_after.txt \
  -o /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_OGB2019_tbl2asn.sqn

tbl2asn \
  -t /home/theaven/scratch/manuscript4/OGB2021_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_leucotricha_OGB2021_fix.tbl \
  -a s -V vb \
  -Z /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/discrep_report_after.txt \
  -o /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_OGB2021_tbl2asn.sqn

tbl2asn \
  -t /home/theaven/scratch/manuscript4/P112020_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_leucotricha_P112020_fix.tbl \
  -a s -V vb \
  -Z /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/discrep_report_after.txt \
  -o /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_P112020_tbl2asn.sqn

tbl2asn \
  -t /home/theaven/scratch/manuscript4/SCOTT2020_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_SCOTT2021_fix.tbl \
  -a s -V vb \
  -Z /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/discrep_report_after.txt \
  -o /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_SCOTT2020_tbl2asn.sqn

#######################################################################################

conda activate table2asn
cd /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results

table2asn \
  -t /home/theaven/scratch/manuscript4/DRCT72020_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f Podosphaera_aphanis_DRCT72020_fix.tbl \
  -a s -V vb -M n -euk -c w -j "[organism=Podosphaera aphanis] [isolate=DRCT72020] [moltype=genomic DNA] [gcode=1]" -Z \
  -o Podosphaera_aphanis_DRCT72020_table2asn.sqn > table2asn.log 2>&1

grep "SEQ_FEAT.InternalStop" Podosphaera_aphanis_DRCT72020_table2asn.val \
| sed -n 's/.*-> \[\(gnl|ncbi|[^]]*\)\].*/\1/p' \
| sort -u > internal_stop_protein_ids.txt

awk -F'\t' '
$4=="locus_tag"  { lt=$5 }
$4=="protein_id" { if (lt!="") print $5 "\t" lt }
' Podosphaera_aphanis_DRCT72020_fix.tbl \
| sort -t $'\t' -k1,1 -u > protein_to_locus.tsv

join -t $'\t' internal_stop_protein_ids.txt protein_to_locus.tsv \
| cut -f2 \
| sort -u > internal_stop_locus_tags.txt


awk -f /home/theaven/git_repos/Scripts/gruffalo/mark_pseudogenes.awk Podosphaera_aphanis_DRCT72020_fix.tbl \
> Podosphaera_aphanis_DRCT72020_fix_pseudo.tbl

awk -F'\t' 'BEGIN{OFS="\t"}
{
  if ($(NF-1)=="EC_number") {
    v=$NF
    gsub(/[[:space:]]+$/,"",v)
    n=split(v,a,".")
    if (n==1)      $NF=v".-.-.-"
    else if (n==2) $NF=v".-.-"
    else if (n==3) $NF=v".-"
  }
  print
}' Podosphaera_aphanis_DRCT72020_fix_pseudo.tbl > Podosphaera_aphanis_DRCT72020_fix_pseudo_ECfixed.tbl

table2asn \
  -t /home/theaven/scratch/manuscript4/DRCT72020_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f Podosphaera_aphanis_DRCT72020_fix_pseudo_ECfixed.tbl \
  -a s -V vb -M n -euk -c w -j "[organism=Podosphaera aphanis] [isolate=DRCT72020] [moltype=genomic DNA] [gcode=1]" -Z \
  -o Podosphaera_aphanis_DRCT72020_table2asn2.sqn > table2asn.log 2>&1

grep "SEQ_FEAT.InternalStop" Podosphaera_aphanis_DRCT72020_table2asn2.val \
| sed -n 's/.*-> \[\(gnl|ncbi|[^]]*\)\].*/\1/p' \
| sort -u 

gnl|ncbi|g13673.t1
gnl|ncbi|g14591.t1

grep "SEQ_INST.StopInProtein" Podosphaera_aphanis_DRCT72020_table2asn2.val \
| sed -n 's/.*-> \[\(gnl|ncbi|[^]]*\)\].*/\1/p' \
| sort -u 

gnl|ncbi|g13673.t1
gnl|ncbi|g14591.t1

grep -n "MissingGeneXref" -n Podosphaera_aphanis_DRCT72020_table2asn2.val 

cp Podosphaera_aphanis_DRCT72020_fix_pseudo_ECfixed.tbl Podosphaera_aphanis_DRCT72020_fix_pseudo_ECfixed_fix.tbl
nano Podosphaera_aphanis_DRCT72020_fix_pseudo_ECfixed_fix.tbl

table2asn \
  -t /home/theaven/scratch/manuscript4/DRCT72020_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f Podosphaera_aphanis_DRCT72020_fix_pseudo_ECfixed_fix.tbl \
  -a s -V vb -M n -euk -c w -j "[organism=Podosphaera aphanis] [isolate=DRCT72020] [moltype=genomic DNA] [gcode=1]" -Z \
  -o Podosphaera_aphanis_DRCT72020_table2asn3.sqn > table2asn.log 2>&1

#######################################################################################

conda activate table2asn
cd /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/funannotate3/annotate_results

table2asn \
  -t /home/theaven/scratch/manuscript4/DRCT2021_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f Podosphaera_aphanis_DRCT72021_fix.tbl \
  -a s -V vb -M n -euk -c w -j "[organism=Podosphaera aphanis] [isolate=DRCT72021] [moltype=genomic DNA] [gcode=1]" -Z \
  -o Podosphaera_aphanis_DRCT72021_table2asn.sqn > table2asn.log 2>&1

grep "SEQ_FEAT.InternalStop" Podosphaera_aphanis_DRCT72021_table2asn.val \
| sed -n 's/.*-> \[\(gnl|ncbi|[^]]*\)\].*/\1/p' \
| sort -u > internal_stop_protein_ids.txt

awk -F'\t' '
$4=="locus_tag"  { lt=$5 }
$4=="protein_id" { if (lt!="") print $5 "\t" lt }
' Podosphaera_aphanis_DRCT72021_fix.tbl \
| sort -t $'\t' -k1,1 -u > protein_to_locus.tsv

join -t $'\t' internal_stop_protein_ids.txt protein_to_locus.tsv \
| cut -f2 \
| sort -u > internal_stop_locus_tags.txt


awk -f /home/theaven/git_repos/Scripts/gruffalo/mark_pseudogenes.awk Podosphaera_aphanis_DRCT72021_fix.tbl \
> Podosphaera_aphanis_DRCT72021_fix_pseudo.tbl

awk -F'\t' 'BEGIN{OFS="\t"}
{
  if ($(NF-1)=="EC_number") {
    v=$NF
    gsub(/[[:space:]]+$/,"",v)
    n=split(v,a,".")
    if (n==1)      $NF=v".-.-.-"
    else if (n==2) $NF=v".-.-"
    else if (n==3) $NF=v".-"
  }
  print
}' Podosphaera_aphanis_DRCT72021_fix_pseudo.tbl > Podosphaera_aphanis_DRCT72021_fix_pseudo_ECfixed.tbl

table2asn \
  -t /home/theaven/scratch/manuscript4/DRCT2021_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f Podosphaera_aphanis_DRCT72021_fix_pseudo_ECfixed.tbl \
  -a s -V vb -M n -euk -c w -j "[organism=Podosphaera aphanis] [isolate=DRCT72021] [moltype=genomic DNA] [gcode=1]" -Z \
  -o Podosphaera_aphanis_DRCT72021_table2asn2.sqn > table2asn.log 2>&1

grep "SEQ_FEAT.InternalStop" Podosphaera_aphanis_DRCT72021_table2asn2.val \
| sed -n 's/.*-> \[\(gnl|ncbi|[^]]*\)\].*/\1/p' \
| sort -u 

gnl|ncbi|g14345.t1

grep "SEQ_INST.StopInProtein" Podosphaera_aphanis_DRCT72021_table2asn2.val \
| sed -n 's/.*-> \[\(gnl|ncbi|[^]]*\)\].*/\1/p' \
| sort -u 

gnl|ncbi|g14345.t1

grep -n "MissingGeneXref" -n Podosphaera_aphanis_DRCT72021_table2asn2.val 

cp Podosphaera_aphanis_DRCT72021_fix_pseudo_ECfixed.tbl Podosphaera_aphanis_DRCT72021_fix_pseudo_ECfixed_fix.tbl
nano Podosphaera_aphanis_DRCT72021_fix_pseudo_ECfixed_fix.tbl

table2asn \
  -t /home/theaven/scratch/manuscript4/DRCT2021_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f Podosphaera_aphanis_DRCT72021_fix_pseudo_ECfixed_fix.tbl \
  -a s -V vb -M n -euk -c w -j "[organism=Podosphaera aphanis] [isolate=DRCT72021] [moltype=genomic DNA] [gcode=1]" -Z \
  -o Podosphaera_aphanis_DRCT72021_table2asn3.sqn > table2asn.log 2>&1

#######################################################################################

conda activate table2asn
cd /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/funannotate3/annotate_results

table2asn \
  -t /home/theaven/scratch/manuscript4/OGB2019_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f Podosphaera_leucotricha_OGB2019_fix.tbl \
  -a s -V vb -M n -euk -c w -j "[organism=Podosphaera leucotricha] [isolate=OGB2019] [moltype=genomic DNA] [gcode=1]" -Z \
  -o Podosphaera_leucotricha_OGB2019_table2asn.sqn > table2asn.log 2>&1

grep "SEQ_FEAT.InternalStop" Podosphaera_leucotricha_OGB2019_table2asn.val \
| sed -n 's/.*-> \[\(gnl|ncbi|[^]]*\)\].*/\1/p' \
| sort -u > internal_stop_protein_ids.txt

awk -F'\t' '
$4=="locus_tag"  { lt=$5 }
$4=="protein_id" { if (lt!="") print $5 "\t" lt }
' Podosphaera_leucotricha_OGB2019_fix.tbl \
| sort -t $'\t' -k1,1 -u > protein_to_locus.tsv

join -t $'\t' internal_stop_protein_ids.txt protein_to_locus.tsv \
| cut -f2 \
| sort -u > internal_stop_locus_tags.txt


awk -f /home/theaven/git_repos/Scripts/gruffalo/mark_pseudogenes.awk Podosphaera_leucotricha_OGB2019_fix.tbl \
> Podosphaera_leucotricha_OGB2019_fix_pseudo.tbl

awk -F'\t' 'BEGIN{OFS="\t"}
{
  if ($(NF-1)=="EC_number") {
    v=$NF
    gsub(/[[:space:]]+$/,"",v)
    n=split(v,a,".")
    if (n==1)      $NF=v".-.-.-"
    else if (n==2) $NF=v".-.-"
    else if (n==3) $NF=v".-"
  }
  print
}' Podosphaera_leucotricha_OGB2019_fix_pseudo.tbl > Podosphaera_leucotricha_OGB2019_fix_pseudo_ECfixed.tbl

table2asn \
  -t /home/theaven/scratch/manuscript4/OGB2019_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f Podosphaera_leucotricha_OGB2019_fix_pseudo_ECfixed.tbl \
  -a s -V vb -M n -euk -c w -j "[organism=Podosphaera leucotricha] [isolate=OGB2019] [moltype=genomic DNA] [gcode=1]" -Z \
  -o Podosphaera_leucotricha_OGB2019_table2asn2.sqn > table2asn.log 2>&1

grep "SEQ_FEAT.InternalStop" Podosphaera_leucotricha_OGB2019_table2asn2.val \
| sed -n 's/.*-> \[\(gnl|ncbi|[^]]*\)\].*/\1/p' \
| sort -u 

gnl|ncbi|g568.t1

grep "SEQ_INST.StopInProtein" Podosphaera_leucotricha_OGB2019_table2asn2.val \
| sed -n 's/.*-> \[\(gnl|ncbi|[^]]*\)\].*/\1/p' \
| sort -u 

gnl|ncbi|g568.t1

grep -n "MissingGeneXref" -n Podosphaera_leucotricha_OGB2019_table2asn2.val 

cp Podosphaera_leucotricha_OGB2019_fix_pseudo_ECfixed.tbl Podosphaera_leucotricha_OGB2019_fix_pseudo_ECfixed_fix.tbl
nano Podosphaera_leucotricha_OGB2019_fix_pseudo_ECfixed_fix.tbl

table2asn \
  -t /home/theaven/scratch/manuscript4/OGB2019_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f Podosphaera_leucotricha_OGB2019_fix_pseudo_ECfixed_fix.tbl \
  -a s -V vb -M n -euk -c w -j "[organism=Podosphaera leucotricha] [isolate=OGB2019] [moltype=genomic DNA] [gcode=1]" -Z \
  -o Podosphaera_leucotricha_OGB2019_table2asn3.sqn > table2asn.log 2>&1

#######################################################################################

conda activate table2asn
cd /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/funannotate3/annotate_results

table2asn \
  -t /home/theaven/scratch/manuscript4/OGB2021_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f Podosphaera_leucotricha_OGB2021_fix.tbl \
  -a s -V vb -M n -euk -c w -j "[organism=Podosphaera leucotricha] [isolate=OGB2021] [moltype=genomic DNA] [gcode=1]" -Z \
  -o Podosphaera_leucotricha_OGB2021_table2asn.sqn > table2asn.log 2>&1

grep "SEQ_FEAT.InternalStop" Podosphaera_leucotricha_OGB2021_table2asn.val \
| sed -n 's/.*-> \[\(gnl|ncbi|[^]]*\)\].*/\1/p' \
| sort -u > internal_stop_protein_ids.txt

awk -F'\t' '
$4=="locus_tag"  { lt=$5 }
$4=="protein_id" { if (lt!="") print $5 "\t" lt }
' Podosphaera_leucotricha_OGB2021_fix.tbl \
| sort -t $'\t' -k1,1 -u > protein_to_locus.tsv

join -t $'\t' internal_stop_protein_ids.txt protein_to_locus.tsv \
| cut -f2 \
| sort -u > internal_stop_locus_tags.txt


awk -f /home/theaven/git_repos/Scripts/gruffalo/mark_pseudogenes.awk Podosphaera_leucotricha_OGB2021_fix.tbl \
> Podosphaera_leucotricha_OGB2021_fix_pseudo.tbl

awk -F'\t' 'BEGIN{OFS="\t"}
{
  if ($(NF-1)=="EC_number") {
    v=$NF
    gsub(/[[:space:]]+$/,"",v)
    n=split(v,a,".")
    if (n==1)      $NF=v".-.-.-"
    else if (n==2) $NF=v".-.-"
    else if (n==3) $NF=v".-"
  }
  print
}' Podosphaera_leucotricha_OGB2021_fix_pseudo.tbl > Podosphaera_leucotricha_OGB2021_fix_pseudo_ECfixed.tbl

table2asn \
  -t /home/theaven/scratch/manuscript4/OGB2021_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f Podosphaera_leucotricha_OGB2021_fix_pseudo_ECfixed.tbl \
  -a s -V vb -M n -euk -c w -j "[organism=Podosphaera leucotricha] [isolate=OGB2021] [moltype=genomic DNA] [gcode=1]" -Z \
  -o Podosphaera_leucotricha_OGB2021_table2asn2.sqn > table2asn.log 2>&1

grep "SEQ_FEAT.InternalStop" Podosphaera_leucotricha_OGB2021_table2asn2.val \
| sed -n 's/.*-> \[\(gnl|ncbi|[^]]*\)\].*/\1/p' \
| sort -u 

gnl|ncbi|g1231.t1
gnl|ncbi|g15782.t1
gnl|ncbi|g16148.t1
gnl|ncbi|g294.t1
gnl|ncbi|g5256.t1
gnl|ncbi|g9246.t3

grep "SEQ_INST.StopInProtein" Podosphaera_leucotricha_OGB2021_table2asn2.val \
| sed -n 's/.*-> \[\(gnl|ncbi|[^]]*\)\].*/\1/p' \
| sort -u 

gnl|ncbi|g1231.t1
gnl|ncbi|g15782.t1
gnl|ncbi|g16148.t1
gnl|ncbi|g294.t1
gnl|ncbi|g5256.t1
gnl|ncbi|g9246.t3

grep -n "MissingGeneXref" -n Podosphaera_leucotricha_OGB2021_table2asn2.val 

cp Podosphaera_leucotricha_OGB2021_fix_pseudo_ECfixed.tbl Podosphaera_leucotricha_OGB2021_fix_pseudo_ECfixed_fix.tbl
nano Podosphaera_leucotricha_OGB2021_fix_pseudo_ECfixed_fix.tbl

table2asn \
  -t /home/theaven/scratch/manuscript4/OGB2021_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f Podosphaera_leucotricha_OGB2021_fix_pseudo_ECfixed_fix.tbl \
  -a s -V vb -M n -euk -c w -j "[organism=Podosphaera leucotricha] [isolate=OGB2021] [moltype=genomic DNA] [gcode=1]" -Z \
  -o Podosphaera_leucotricha_OGB2021_table2asn3.sqn > table2asn.log 2>&1

#######################################################################################

conda activate table2asn
cd /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/funannotate3/annotate_results

table2asn \
  -t /home/theaven/scratch/manuscript4/P112020_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f Podosphaera_leucotricha_P112020_fix.tbl \
  -a s -V vb -M n -euk -c w -j "[organism=Podosphaera leucotricha] [isolate=P112020] [moltype=genomic DNA] [gcode=1]" -Z \
  -o Podosphaera_leucotricha_P112020_table2asn.sqn > table2asn.log 2>&1

grep "SEQ_FEAT.InternalStop" Podosphaera_leucotricha_P112020_table2asn.val \
| sed -n 's/.*-> \[\(gnl|ncbi|[^]]*\)\].*/\1/p' \
| sort -u > internal_stop_protein_ids.txt

awk -F'\t' '
$4=="locus_tag"  { lt=$5 }
$4=="protein_id" { if (lt!="") print $5 "\t" lt }
' Podosphaera_leucotricha_P112020_fix.tbl \
| sort -t $'\t' -k1,1 -u > protein_to_locus.tsv

join -t $'\t' internal_stop_protein_ids.txt protein_to_locus.tsv \
| cut -f2 \
| sort -u > internal_stop_locus_tags.txt


awk -f /home/theaven/git_repos/Scripts/gruffalo/mark_pseudogenes.awk Podosphaera_leucotricha_P112020_fix.tbl \
> Podosphaera_leucotricha_P112020_fix_pseudo.tbl

awk -F'\t' 'BEGIN{OFS="\t"}
{
  if ($(NF-1)=="EC_number") {
    v=$NF
    gsub(/[[:space:]]+$/,"",v)
    n=split(v,a,".")
    if (n==1)      $NF=v".-.-.-"
    else if (n==2) $NF=v".-.-"
    else if (n==3) $NF=v".-"
  }
  print
}' Podosphaera_leucotricha_P112020_fix_pseudo.tbl > Podosphaera_leucotricha_P112020_fix_pseudo_ECfixed.tbl

table2asn \
  -t /home/theaven/scratch/manuscript4/P112020_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f Podosphaera_leucotricha_P112020_fix_pseudo_ECfixed.tbl \
  -a s -V vb -M n -euk -c w -j "[organism=Podosphaera leucotricha] [isolate=P112020] [moltype=genomic DNA] [gcode=1]" -Z \
  -o Podosphaera_leucotricha_P112020_table2asn2.sqn > table2asn.log 2>&1

grep "SEQ_FEAT.InternalStop" Podosphaera_leucotricha_P112020_table2asn2.val \
| sed -n 's/.*-> \[\(gnl|ncbi|[^]]*\)\].*/\1/p' \
| sort -u 

gnl|ncbi|g12354.t1
gnl|ncbi|g14776.t1

grep "SEQ_INST.StopInProtein" Podosphaera_leucotricha_P112020_table2asn2.val \
| sed -n 's/.*-> \[\(gnl|ncbi|[^]]*\)\].*/\1/p' \
| sort -u 

gnl|ncbi|g12354.t1
gnl|ncbi|g14776.t1

grep -n "MissingGeneXref" -n Podosphaera_leucotricha_P112020_table2asn2.val 

cp Podosphaera_leucotricha_P112020_fix_pseudo_ECfixed.tbl Podosphaera_leucotricha_P112020_fix_pseudo_ECfixed_fix.tbl
nano Podosphaera_leucotricha_P112020_fix_pseudo_ECfixed_fix.tbl

table2asn \
  -t /home/theaven/scratch/manuscript4/P112020_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f Podosphaera_leucotricha_P112020_fix_pseudo_ECfixed_fix.tbl \
  -a s -V vb -M n -euk -c w -j "[organism=Podosphaera leucotricha] [isolate=P112020 [moltype=genomic DNA] [gcode=1]" -Z \
  -o Podosphaera_leucotricha_P112020_table2asn3.sqn > table2asn.log 2>&1

note    nonfunctional due to internal stop codon

#NCBI curator requests more removals:
cp Podosphaera_leucotricha_P112020_fix_pseudo_ECfixed_fix.tbl Podosphaera_leucotricha_P112020_fix_pseudo_ECfixed_fix2.tbl
nano Podosphaera_leucotricha_P112020_fix_pseudo_ECfixed_fix2.tbl

table2asn \
  -t /home/theaven/scratch/manuscript4/P112020_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f Podosphaera_leucotricha_P112020_fix_pseudo_ECfixed_fix2.tbl \
  -a s -V vb -M n -euk -c w -j "[organism=Podosphaera leucotricha] [isolate=P112020] [moltype=genomic DNA] [gcode=1]" -Z \
  -o Podosphaera_leucotricha_P112020_table2asn4.sqn > table2asn.log 2>&1

#######################################################################################

awk '
BEGIN { inserted=0 }
/^>Feature/ { 
  inserted=0; 
  print; 
  next 
}
# First real feature line after >Feature: e.g. "1\t56736\tREFERENCE"
(!inserted && $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ && $3 != "") {
  len=$2
  print "1\t" len "\tsource"
  print "\t\ttransl_table\t1\n"
  inserted=1
}
{ print }
' /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_SCOTT2021_fix.tbl \
> /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_SCOTT2021_fix2.tbl



conda activate table2asn
table2asn \
  -t SCOTT2020_template.sbt \
  -i genome.fsa \
  -f Podosphaera_aphanis_SCOTT2021_fix.tbl \
  -a s -V vb -M n -euk -c w -j "[organism=Podosphaera aphanis] [isolate=SCOTT2021] [moltype=genomic DNA]" -Z \
  -o Podosphaera_aphanis_SCOTT2020_table2asn.sqn > table2asn.log 2>&1

grep "SEQ_FEAT.InternalStop" Podosphaera_aphanis_SCOTT2020_table2asn.val \
| sed -n 's/.*-> \[\(gnl|ncbi|[^]]*\)\].*/\1/p' \
| sort -u > internal_stop_protein_ids.txt

awk -F'\t' '
$4=="locus_tag"  { lt=$5 }
$4=="protein_id" { if (lt!="") print $5 "\t" lt }
' Podosphaera_aphanis_SCOTT2021.tbl \
| sort -t $'\t' -k1,1 -u > protein_to_locus.tsv

join -t $'\t' internal_stop_protein_ids.txt protein_to_locus.tsv \
| cut -f2 \
| sort -u > internal_stop_locus_tags.txt


awk -f /home/theaven/git_repos/Scripts/gruffalo/mark_pseudogenes.awk Podosphaera_aphanis_SCOTT2021.tbl \
> Podosphaera_aphanis_SCOTT2021_pseudo.tbl

awk -F'\t' 'BEGIN{OFS="\t"}
{
  if ($(NF-1)=="EC_number") {
    v=$NF
    gsub(/[[:space:]]+$/,"",v)
    n=split(v,a,".")
    if (n==1)      $NF=v".-.-.-"
    else if (n==2) $NF=v".-.-"
    else if (n==3) $NF=v".-"
  }
  print
}' Podosphaera_aphanis_SCOTT2021_pseudo.tbl > Podosphaera_aphanis_SCOTT2021_pseudo_ECfixed.tbl

awk -F'\t' 'BEGIN{OFS="\t"}
{
  print
  if ($3=="CDS") {
    print "\t\t\ttransl_table\t1"
  }
}' Podosphaera_aphanis_SCOTT2021_pseudo_ECfixed.tbl \
> Podosphaera_aphanis_SCOTT2021_final.tbl


table2asn \
  -t /home/theaven/scratch/manuscript4/SCOTT2020_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f Podosphaera_aphanis_SCOTT2021_pseudo_ECfixed.tbl \
  -a s -V vb -M n -euk -c w -j "[organism=Podosphaera aphanis] [isolate=SCOTT2021] [moltype=genomic DNA] [gcode=1]" -Z \
  -o Podosphaera_aphanis_SCOTT2020_table2asn2.sqn > table2asn.log 2>&1

grep "SEQ_FEAT.InternalStop" Podosphaera_aphanis_SCOTT2020_table2asn2.val \
| sed -n 's/.*-> \[\(gnl|ncbi|[^]]*\)\].*/\1/p' \
| sort -u 

gnl|ncbi|g12516.t1
gnl|ncbi|g13318.t1
gnl|ncbi|g9415.t1

grep "SEQ_INST.StopInProtein" Podosphaera_aphanis_SCOTT2020_table2asn2.val \
| sed -n 's/.*-> \[\(gnl|ncbi|[^]]*\)\].*/\1/p' \
| sort -u 

gnl|ncbi|g12516.t1
gnl|ncbi|g13318.t1
gnl|ncbi|g9415.t1

grep -n "gnl|ncbi|g67.t1" -A40 Podosphaera_aphanis_SCOTT2021_pseudo_ECfixed.tbl | grep EC_number

awk '
BEGIN{inCDS=0; hypo=0}

# feature header line
/^[0-9]+[[:space:]]+[0-9]+[[:space:]]+CDS[[:space:]]*$/ {inCDS=1; hypo=0; print; next}
# leaving CDS block when next feature starts
/^[0-9]+[[:space:]]+[0-9]+[[:space:]]+(gene|mRNA|tRNA|rRNA|misc_feature|exon|intron)[[:space:]]*$/ {inCDS=0; hypo=0; print; next}

{
  if (inCDS) {
    if ($0 ~ /^[[:space:]]+product[[:space:]]+hypothetical protein/ ||
        $0 ~ /^[[:space:]]+product[[:space:]]+unknown protein/) hypo=1;

    if (hypo && $0 ~ /^[[:space:]]+EC_number[[:space:]]+/) next;
  }
  print
}
' Podosphaera_aphanis_SCOTT2021_pseudo_ECfixed.tbl > Podosphaera_aphanis_SCOTT2021_pseudo_ECfixed_noEC_on_hypo.tbl

grep -n "gnl|ncbi|g67.t1" -A40 Podosphaera_aphanis_SCOTT2021_pseudo_ECfixed_noEC_on_hypo.tbl | grep EC_number

table2asn \
  -t /home/theaven/scratch/manuscript4/SCOTT2020_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f Podosphaera_aphanis_SCOTT2021_pseudo_ECfixed_noEC_on_hypo.tbl \
  -a s -V vb -M n -euk -c w -j "[organism=Podosphaera aphanis] [isolate=SCOTT2021] [moltype=genomic DNA] [gcode=1]" -Z \
  -o Podosphaera_aphanis_SCOTT2020_table2asn3.sqn > table2asn.log 2>&1

TBL=/home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_SCOTT2021_pseudo_ECfixed.tbl
egg=/home/theaven/scratch/manuscript4/eggnog/Pod_aph_SCOTT2020-2/Pod_aph_SCOTT2020-2.emapper.annotations
singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/python3.sif python3 ~/git_repos/Scripts/gruffalo/ec_unknown_from_tbl2.py "$TBL" --emapper "$egg" --join-on protein_id -o "$(dirname $TBL)"


Fun_Dir=$(dirname "$TBL")
singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/python3.sif python3 ~/git_repos/Scripts/gruffalo/remove_ec_from_problem_tbl.py \
  --tbl "$TBL" \
  --hits /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_SCOTT2021_EC_UNKNOWN_with_eggnog.tsv \
  --out "$Fun_Dir"/"$(basename $TBL | sed 's@.tbl@_fix.tbl@g')"

table2asn \
  -t /home/theaven/scratch/manuscript4/SCOTT2020_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f Podosphaera_aphanis_SCOTT2021_pseudo_ECfixed_fix.tbl \
  -a s -V vb -M n -euk -c w -j "[organism=Podosphaera aphanis] [isolate=SCOTT2021] [moltype=genomic DNA] [gcode=1]" -Z \
  -o Podosphaera_aphanis_SCOTT2020_table2asn4.sqn > table2asn.log 2>&1

grep -n "MissingGeneXref" -n Podosphaera_aphanis_SCOTT2020_table2asn4.val | head -n 20

cp Podosphaera_aphanis_SCOTT2021_pseudo_ECfixed_fix.tbl Podosphaera_aphanis_SCOTT2021_pseudo_ECfixed_fix_2.tbl
nano Podosphaera_aphanis_SCOTT2021_pseudo_ECfixed_fix_2.tbl

table2asn \
  -t /home/theaven/scratch/manuscript4/SCOTT2020_template.sbt \
  -i /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/funannotate3/annotate_misc/tbl2asn/genome.fsa \
  -f Podosphaera_aphanis_SCOTT2021_pseudo_ECfixed_fix_2.tbl \
  -a s -V vb -M n -euk -c w -j "[organism=Podosphaera aphanis] [isolate=SCOTT2021] [moltype=genomic DNA] [gcode=1]" -Z \
  -o Podosphaera_aphanis_SCOTT2020_table2asn5.sqn > table2asn.log 2>&1

###########################################################################################################################
cp /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_SCOTT2020_tbl2asn.sqn down_20260204/.
cp /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_P112020_tbl2asn.sqn down_20260204/. 
cp /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_OGB2021_tbl2asn.sqn down_20260204/. 
cp /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_OGB2019_tbl2asn.sqn down_20260204/. 
cp /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_DRCT72020_tbl2asn.sqn down_20260204/. 
cp /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_DRCT72021_tbl2asn.sqn down_20260204/. 


grep -c '^>Feature ' /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_SCOTT2021_fix.tbl
grep -c '^>Feature ' /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_SCOTT2021.tbl

# IDs from FASTA
grep '^>' Podosphaera_aphanis_SCOTT2021.contigs.fsa | sed 's/^>//' | awk '{print $1}' | sort -u > fsa.ids

# IDs referenced by tbl (original)
grep '^>Feature ' Podosphaera_aphanis_SCOTT2021.tbl | awk '{print $2}' | sort -u > tbl_orig.ids

# IDs referenced by tbl (edited)
grep '^>Feature ' Podosphaera_aphanis_SCOTT2021_fix.tbl | awk '{print $2}' | sort -u > tbl_fix.ids

echo "orig tbl missing from fasta:"
comm -23 tbl_orig.ids fsa.ids | head

echo "fix tbl missing from fasta:"
comm -23 tbl_fix.ids fsa.ids | head


# 1) IDs referenced by the feature table
grep '^>Feature' /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_SCOTT2021_fix.tbl | awk '{print $2}' | sort -u > tbl.ids

# 2) IDs present in the FASTA
grep '^>' /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_SCOTT2021.contigs.fsa | sed 's/^>//' | awk '{print $1}' | sort -u > fsa.ids

# 3) IDs in tbl but missing from fasta (this triggers your error)
comm -23 tbl.ids fsa.ids | head

# 4) IDs in fasta but not referenced by tbl (usually fine)
comm -13 tbl.ids fsa.ids | head


singularity exec --bind /mnt/shared/scratch/theaven,/mnt/shared/projects/niab/theaven /home/theaven/git_repos/Containers/python3.sif python3 - <<'PY'
import pandas as pd
p = "/home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_DRCT72020_EC_UNKNOWN_with_eggnog.tsv"
df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
cols = [c for c in ["protein_id","protein_id_norm","eggnog_query_norm","query","transcript_id","transcript_id_norm"] if c in df.columns]
print("rows:", len(df))
for c in cols:
    print(c, "nonempty:", (df[c]!="").sum(), "unique_nonempty:", df.loc[df[c]!="", c].nunique())
PY



srun --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=16G \
     --nodelist=n24-64-384-wesley \
     --pty bash

for TBL in $(ls /home/theaven/projects/niab/theaven/gene_pred/*/*/codingquarry/rep_modeling/final/funannotate3/annotate_results/*.tbl | grep -v 'part'); do
Fun_Dir=$(dirname "$TBL")
awk -F'\t' '
BEGIN{cur=""; inCDS=0; ec=""; hyp=0}
function flush() {
  if(inCDS && ec!="" && hyp && cur!="") print cur "\t" ec
}
# New feature starts on lines like: start end featurekey
NF>=3 && $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ && $3 ~ /^[A-Za-z_]+$/ {
  flush()
  inCDS = ($3=="CDS")
  ec=""; hyp=0
  next
}
# locus_tag lines look like: \t\t\tlocus_tag\tK3495_...
$4=="locus_tag" {cur=$5; next}

# Only look at product/EC_number within CDS features
inCDS && $4=="EC_number" {ec=$5; next}
inCDS && $4=="product" {
  p=tolower($5)
  if(p=="hypothetical protein" || p=="unknown protein" || p=="uncharacterized protein" || p ~ /hypothetical|unknown|uncharacterized/) hyp=1
  next
}
END{flush()}
' "$TBL" | sort -u > "$Fun_Dir"/ec_hyp.locus_tags.txt

awk -F'\t' '
BEGIN{
  OFS="\t"
}
{
  split($2,a,".")
  c=a[1]+0
  prod="putative enzyme"
  if(c==1) prod="putative oxidoreductase"
  else if(c==2) prod="putative transferase"
  else if(c==3) prod="putative hydrolase"
  else if(c==4) prod="putative lyase"
  else if(c==5) prod="putative isomerase"
  else if(c==6) prod="putative ligase"
  else if(c==7) prod="putative translocase"
  print $1,".",prod
}
' "$Fun_Dir"/ec_hyp.locus_tags.txt > "$Fun_Dir"/fix_ec_putative_enzyme.tsv
done

funannotate annotate \
  --gff /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3 \
  --fasta /home/theaven/projects/niab/theaven/genomes/Pod_leu_OGB2021-2_clean_renamed_fix_renamed2.fasta \
  --species "Podosphaera leucotricha" --strain OGB2021 \
  --iprscan /home/theaven/projects/niab/theaven/Pod_leu_OGB2021-2_clean_renamed.pep.fasta.interproscan.xml \
  --eggnog /home/theaven/scratch/manuscript4/eggnog/Pod_leu_OGB2021-2/Pod_leu_OGB2021-2.emapper.annotations \
  --busco_db ascomycota \
  --cpus 16 \
  -d /home/theaven/scratch/manuscript4/db/funannotate \
  -o /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/funannotate3 \
  --sbt /home/theaven/scratch/manuscript4/OGB2021_template.sbt \
  --rename MP843_ --force --fix /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2021_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/fix_ec_putative_enzyme.tsv

funannotate annotate \
  --gff /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3 \
  --fasta /home/theaven/projects/niab/theaven/genomes/Pod_leu_OGB2019-2_clean_renamed_fix_renamed2.fasta \
  --species "Podosphaera leucotricha" --strain OGB2019 \
  --iprscan /home/theaven/projects/niab/theaven/Pod_leu_OGB2019-2_clean_renamed.pep.fasta.interproscan.xml \
  --eggnog /home/theaven/scratch/manuscript4/eggnog/Pod_leu_OGB2019-2/Pod_leu_OGB2019-2.emapper.annotations \
  --busco_db ascomycota \
  --cpus 16 \
  -d /home/theaven/scratch/manuscript4/db/funannotate \
  -o /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/funannotate3 \
  --sbt /home/theaven/scratch/manuscript4/OGB2019_template.sbt \
  --rename MP842_ --force --fix /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenpOGB2019_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/fix_ec_putative_enzyme.tsv

funannotate annotate \
  --gff /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3 \
  --fasta /home/theaven/projects/niab/theaven/genomes/Pod_leu_OGBp112020_clean_renamed_renamed2.fasta \
  --species "Podosphaera leucotricha" --strain P112020 \
  --iprscan /home/theaven/projects/niab/theaven/Pod_leu_OGBp112020_clean_renamed.pep.fasta.interproscan.xml \
  --eggnog /home/theaven/scratch/manuscript4/eggnog/Pod_leu_OGBp112020/Pod_leu_OGBp112020.emapper.annotations \
  --busco_db ascomycota \
  --cpus 16 \
  -d /home/theaven/scratch/manuscript4/db/funannotate \
  -o /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/funannotate3 \
  --sbt /home/theaven/scratch/manuscript4/P112020_template.sbt \
  --rename K3496_ --force --fix /home/theaven/projects/niab/theaven/gene_pred/P_leucotricha/THeavenp11_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/fix_ec_putative_enzyme.tsv

funannotate annotate \
  --gff /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3 \
  --fasta /home/theaven/projects/niab/theaven/genomes/Pod_aph_SCOTT2020-2_clean_renamed_fix_renamed2.fasta \
  --species "Podosphaera aphanis" --strain SCOTT2021 \
  --iprscan /home/theaven/projects/niab/theaven/Pod_aph_SCOTT2020-2_clean_renamed.pep.fasta.interproscan.xml \
  --eggnog /home/theaven/scratch/manuscript4/eggnog/Pod_aph_SCOTT2020-2/Pod_aph_SCOTT2020-2.emapper.annotations \
  --busco_db ascomycota \
  --cpus 16 \
  -d /home/theaven/scratch/manuscript4/db/funannotate \
  -o /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/funannotate3 \
  --sbt /home/theaven/scratch/manuscript4/SCOTT2020_template.sbt \
  --rename LBW10_ --force --fix /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenSCOTT2020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/fix_ec_putative_enzyme.tsv

funannotate annotate \
  --gff /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3 \
  --fasta /home/theaven/projects/niab/theaven/genomes/Pod_aph_DRT72021-2_clean_renamed_fix_renamed2.fasta \
  --species "Podosphaera aphanis" --strain DRCT72021 \
  --iprscan /home/theaven/projects/niab/theaven/Pod_aph_DRT72021-2_clean_renamed.pep.fasta.interproscan.xml \
  --eggnog /home/theaven/scratch/manuscript4/eggnog/Pod_aph_DRT72021-2/Pod_aph_DRT72021-2.emapper.annotations \
  --busco_db ascomycota \
  --cpus 16 \
  -d /home/theaven/scratch/manuscript4/db/funannotate \
  -o /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/funannotate3 \
  --sbt /home/theaven/scratch/manuscript4/DRCT2021_template.sbt \
  --rename MP841_ --fix /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72021_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/fix_ec_putative_enzyme.tsv --force

funannotate annotate \
  --gff /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3 \
  --fasta /home/theaven/projects/niab/theaven/genomes/Pod_aph_DRT72020_clean_renamed_renamed2.fasta \
  --species "Podosphaera aphanis" --strain DRCT72020 \
  --iprscan /home/theaven/projects/niab/theaven/Pod_aph_DRT72020_clean_renamed.pep.fasta.interproscan.xml \
  --eggnog /home/theaven/scratch/manuscript4/eggnog/Pod_aph_DRT72020/Pod_aph_DRT72020.emapper.annotations \
  --busco_db ascomycota \
  --cpus 16 \
  -d /home/theaven/scratch/manuscript4/db/funannotate \
  -o /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/funannotate4 \
  --sbt /home/theaven/scratch/manuscript4/DRCT72020_template.sbt \
  --rename K3495_ --fix /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/fix_ec_putative_enzyme.tsv --force 


ID=$(head -n 1 /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/fix_ec_putative_enzyme.tsv | cut -f1)
echo "$ID"
grep -n "$ID" /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_DRCT72020.tbl



funannotate annotate \
  --gff /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/frozen_genes_appended_renamed.gff3 \
  --fasta /home/theaven/projects/niab/theaven/genomes/Pod_aph_DRT72020_clean_renamed_renamed2.fasta \
  --species "Podosphaera aphanis" --strain DRCT72020 \
  --iprscan /home/theaven/projects/niab/theaven/Pod_aph_DRT72020_clean_renamed.pep.fasta.interproscan.xml \
  --eggnog /home/theaven/scratch/manuscript4/eggnog/Pod_aph_DRT72020/Pod_aph_DRT72020.emapper.annotations \
  --busco_db ascomycota \
  --cpus 16 \
  -d /home/theaven/scratch/manuscript4/db/funannotate \
  -o /home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/funannotate3 \
  --sbt /home/theaven/scratch/manuscript4/DRCT72020_template.sbt \
  --rename K3495_ \
  --fix /home/theaven/scratch/manuscript4/fix_ec_byclass.current.protein.tsv \
  --force


```
```bash
TBL=/home/theaven/projects/niab/theaven/gene_pred/P_aphanis/THeavenDRCT72020_1/codingquarry/rep_modeling/final/funannotate3/annotate_results/Podosphaera_aphanis_DRCT72020.tbl
OUTTBL=/path/to/Podosphaera_aphanis_DRCT72020.fixed.tbl

awk '
BEGIN{inCDS=0; hasEC=0; ecClass=0}

function classprod(c){
  if(c==1) return "putative oxidoreductase"
  if(c==2) return "putative transferase"
  if(c==3) return "putative hydrolase"
  if(c==4) return "putative lyase"
  if(c==5) return "putative isomerase"
  if(c==6) return "putative ligase"
  if(c==7) return "putative translocase"
  return "putative enzyme"
}

{
  gsub(/\r$/,"")

  # Start of a new feature: "start end key"
  if($0 ~ /^[0-9]+[ \t]+[0-9]+[ \t]+/){
    inCDS = ($0 ~ /^[0-9]+[ \t]+[0-9]+[ \t]+CDS([ \t]|$)/)
    hasEC=0; ecClass=0
    print
    next
  }

  # EC_number inside CDS → record class
  if(inCDS && $0 ~ /^[ \t]+EC_number[ \t]+/){
    hasEC=1
    ec=$0
    sub(/^[ \t]+EC_number[ \t]+/,"",ec)
    split(ec,a,".")
    ecClass=a[1]+0
    print
    next
  }

  # If this CDS has an EC number, rewrite bad product names
  if(inCDS && hasEC && $0 ~ /^[ \t]+product[ \t]+/){
    prod=$0
    sub(/^[ \t]+product[ \t]+/,"",prod)
    p=tolower(prod)
    if(p ~ /(hypothetical|unknown|uncharacterized|protein of unknown|putative protein)/){
      # preserve indentation and keyword; replace only the value
      sub(/product[ \t]+.*/,"product\t" classprod(ecClass))
      print
      next
    }
  }

  print
}
' "$TBL" > "$OUTTBL"

```
