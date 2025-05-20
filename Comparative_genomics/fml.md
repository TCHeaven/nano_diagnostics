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

