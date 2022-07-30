
<p align="center"> 
<img src="images/IceR.jpg" style="width: 100%; height: 100%"/>
</p>

### Description
Label-free proteomics enables the unbiased quantification of thousands of proteins across large sample cohorts. Commonly used mass spectrometry-based proteomic workflows rely on data dependent acquisition (DDA). However, its stochastic selection of peptide features for fragmentation-based identification inevitably results in high rates of missing values, which prohibits the integration of larger cohorts as the number of recurrently detected peptides is a limiting factor. Peptide identity propagation (PIP) can mitigate this challenge, allowing to transfer sequencing information between samples. However, despite the promise of these approaches, current methods remain limited either in sensitivity or reliability and there is a lack of robust and widely applicable software. To address this, we here present IceR, an efficient and user-friendly quantification workflow introducing a hybrid PIP approach with superior quantification precision, accuracy, reliability and data completeness. IceR is available as an easy to-use R-package incorporating a graphical user interface and comprehensive quality control measures.

### Installation
#### Win10
The following installations are required for Windows 10:

 - [R](https://cran.r-project.org/bin/windows/base/) (Version between 3.6.3 and 4.2.1)
 - Optional: [Rtools](https://cran.r-project.org/bin/windows/Rtools/history.html) (select required version)
 - Optional: [ProteoWizard](http://proteowizard.sourceforge.net/download.html)
 - Optional: [RStudio](https://rstudio.com/products/rstudio/download/)
 
During installation, please keep default settings and follow respective instructions.

If an error occurs during installation of the package rJava, it could indicate that Java is not properly installed. Please follow the installation instructions for Java and rJava (installing suitable JDK version should solve the issue) 

Next, we install the IceR package from GitHub (development version)
```r
install.packages("devtools")
devtools::install_github("mathiaskalxdorf/IceR",ref = "develop")
```
If everything wents fine, IceR should be installed and can be used (e.g. with the GUI) with the following lines:
```r
library(IceR)
runIceR()
```

#### Ubuntu
Sorry for potentially unclear and/or amateurish installation instructions for Ubuntu.

The following steps have to be done using the terminal:

```t
sudo apt-get update
sudo apt-get install r-base
sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev libgit2-dev libv8-dev
sudo apt-get install r-base-dev default-jdk
```

Open R in administrator mode:

```t
sudo -i R
```

In the R console trigger installation of R packages:

```r
install.packages('devtools')
devtools::install_github("mathiaskalxdorf/IceR",ref = "develop")
```

If successful, close and restart R (not in administrator mode):

```t
R
```

... and try running IceR:

```r
library(IceR)
runIceR()
```

ProteoWizard and its related msCovert is currently (to my knowledge) not available for Ubuntu. Hence, processing of TIMS-ToF pro data is currently not supported on Ubuntu. In case of Thermo MS data, the conversion into required mzXML files using msConvert can not be triggered automatically but have to be manually converted. Resulting mzXML files have to be manually located in the folder "mzXML" within the same folder which contains the raw files.

#### macOS
Sorry for potentially unclear and/or amateurish installation instructions for macOS.

The following installations are required for macOS:

 - [R](https://cran.r-project.org/bin/macosx/) (Version 4.0.4)
 - [RStudio](https://rstudio.com/products/rstudio/download/)
 
In case of macOS (tested for macOS Big Sur 11.1.0) it is important to install and use RStudio as currently a bug in tctlk package prevents IceR from running properly when started directly from R console (tested for R 4.0.4).

Furthermore, please check that a latest version of [Java](https://www.java.com/de/download/) is installed.

In case of macOS, we additionally have to install [XQuartz](https://www.xquartz.org/)

Start RStudio in administrator/superuser mode and install IceR from GitHub (development version) using devtools:

```r
install.packages('devtools')
devtools::install_github("mathiaskalxdorf/IceR",ref = "develop")
```

If an error occurs while installing package utf8, please follow the displayed instructions.

If installation worked successfully, try running IceR:

```r
library(IceR)
runIceR()
```

ProteoWizard and its related msCovert is currently (to my knowledge) not available for macOS. Hence, processing of TIMS-ToF pro data is currently not supported on macOS. In case of Thermo MS data, the conversion into required mzXML files using msConvert can not be triggered automatically but have to be manually converted. Resulting mzXML files have to be manually located in the folder "mzXML" within the same folder which contains the raw files.

### Prerequisites
IceR was so far tested on Windows 10, Ubuntu 20.04.1, and macOS 11.1.0.

The current version of IceR requires raw MS files (from Thermo Mass Spectrometers or from Bruker TIMS-ToF pro) to be preprocessed with MaxQuant (tested for Versions 1.5.1.2, 1.6.12, 1.6.14, and 2.0.3.0, versions in between should work as well). Alternatively, the user can supply preprocessed data from any other pipeline by manually providing the required input data in a predefined format. The data has to be supplied in the following 3 tab-separated txt files:

 - features.txt table which containes all information about detected features in the following columns:
    - Raw.file - Name of the Raw file
    - Charge - Charge state (e.g. 1, 2, 3)
    - Mass - Monoisotopic mass of the feature
    - m.z - Mass to charge ratio of feature
    - Uncalibrated.m.z - (Optional) Uncalibrated mass to charge ratio of feature
    - Max.intensity.m.z.0 - (Optional) Mass to charge ratio at which highest intensity of ions was detected, typically m/z which was selected as precursor subsequent MS2
    - Retention.time - Chromatographic retention time in minutes of feature
    - Calibrated.retention.time - (Optional) Corrected chromatographic retention time in minutes of feature
    - Retention.Length - Elution peak width in minutes of feature
    - Sequence - Observed peptide sequence if feature was sequenced
    - Modifications - Observed post-translational modification of feature. If unmodified, should be `Unmodified`
    - Score - Any score value for the spectra to peptide matching with normal distributed scores and higher score representing higher significance
    - Proteins - Matching protein identifier (e.g. Uniprot ID) if feature was sequenced
    - MSMS.Scan.Numbers - (Optional) Indices of MSMS spectra which matched to the peptide sequence of feature
    - Intensity - Extracted feature intensity (raw, unloged) by preprocessing pipeline. Should not be `NA`

 - preprocess_quant.txt table which contains peptide-level quantifications observed by preprocessing pipeline in the following columns:
    - Sequence - Observed peptide sequence
    - ID - Protein identifier of matching peptide sequence (e.g. Uniprot ID)
    - Intensity - At least 2 columns containing intensities (in log2) of samples. Column names should start with `Intensity.` followed by the sample name e.g. `Intensity.E3_R1` and `Intensity.E3_R2`. Intensity columns should be ordered by name.
        
 - id_mapping.txt table which contains information for mapping from protein IDs to gene names in the following columns:
    - ID - Protein identifier (e.g. Uniprot ID)
    - Gene_Name - Corresponding gene name


Important general notes: 

 - File paths and files should, if possible, not contain blank spaces.
 - File paths and files should not contain special letters like e.g. any of the following: !"§$%&/()=?#*~+-,.
 - (Preprocess by MaxQuant) Fasta file should be parsed correctly. In this case the column "Gene names" can be found in the proteinGroups.txt.
 - IceR was yet only tested for the analysis of single-shot sample data but not for fractionated sample data.

Raw Orbitrap files are required in the mzXML format and raw TIMS-ToF data has to be converted into a readable format. If msconvert from ProteoWizard is installed (Windows-only), IceR triggers conversion automatically. Otherwise, the user has to place the converted files in the folder "mzXML" in case of Orbitrap data within the folder containing the raw files. Conversion of Bruker TIMS-ToF pro data is currently only possible if msConvert is installed. 

Generally, the more ressources are available the faster the data analysis can be performed. Still, IceR can run reasonably well even on a normal PC/Laptop. However, in case of TIMS-ToF data a potent machine will be required due to the huge file sizes. Here it is recommended to at least allocate 128 gb of RAM to enable at maximum 3 samples to be processed in parallel. Furthermore, TIMS-ToF raw data currently has to be converted (automatically triggered by IceR, requires msConvert to be installed) into a readable format which requires several hours per sample and at least 50 gb of space on the machines home drive. A future version of IceR will implement a faster conversion method.

### Example (Windows 10)
For this example, we will require [ProteoWizard](http://proteowizard.sourceforge.net/download.html) to be installed as raw files will have to be converted.

The example data set is stored [here](https://drive.google.com/drive/folders/1te8awhyliY4vKCCxjUdwJNJD_YlYjHJf?usp=sharing).

Please download the raw files and (to speed things up) the corresponding MaxQuant results. This example data set consists of 4 tool MS files representing 2 replicates of human lysate spiked with 3 % E. coli lysate and 2 replicates of human lysate spiked with 9 % E. coli lysate. This data was acquired on a Q-Exactive HF machine. Supplied MaxQuant results were generated using MaxQuant V1.5.1.2. Used search database from UniProt is supplied.

After running ...
```r
library(IceR)
runIceR()
```
... the GUI opens.
<p align="center"> 
<img src="images/IceR_gui.jpg" style="width: 100%; height: 100%"/>
</p>

It allows setting up the IceR run. Among others, the following parameters can be modified:

 - Paths to raw files, MaxQuant results files and IceR output folder
 - MassSpec Mode switching between Orbitrap and TIMS-ToF data (with or without using TIMS dimension)
 - Multiplicity Mode switching between label-free (1 - LFQ) and SILAC (2 - SILAC or 3 - SILAC) and adjusting corresponding label isotopes
 - Minimal retention (RT) and m/z feature alignment windows during feature-based PIP
 - Kernel density estimation (KDE) resolution (grid points per dimension, increasing resolution increases computation workload but also increases resolution of peak detections) 
 - Number of threads to be used during IceR workflow (please check on your system specifications. Should not be higher than numbers of samples to be analyzed. Number of threads should be also adapted to available RAM. In case of TIMS-ToF data, a single 2h gradient acquisition run will require ~ 40 gb of RAM, hence, number of threads will be automatically limited in this case to 3 parallel threads)

We can keep all settings at default (changes of some parameters are disabled at the moment) but reduce number of threads to 4. 

Please specify the path to the downloaded raw files and MaxQuant results by clicking on the respective "Choose directory" buttons. Similarly, please specify a folder where IceR results should be stored.

After clicking on "Run", the IceR workflow will start. A detailed description of the individual steps can be found in the original publication [IceR](https://www.biorxiv.org/content/10.1101/2020.11.01.363101v1.full)

A recent computational system should be able to complete the IceR workflow for the example data set within 2 - 3 hours.

When IceR finished, some quality control plots for the alignment and quantification can be visualized by clicking on the tap "QC" and subsequently on the respective button:

<p align="center"> 
<img src="images/IceR_gui_QC_alignment.jpg" style="width: 100%; height: 100%"/>
</p>

<p align="center"> 
<img src="images/IceR_gui_QC_quantification.jpg" style="width: 100%; height: 100%"/>
</p>

Subsequently, we can investigate results in more detail and compare results to MaxQuant outputs.

We need 3 additional R-packages to be installed for the following steps.
```{r}
install.packages("BiocManager")
BiocManager::install(c("PECA","limma","vsn"))
```
Now we load IceR library.
```{r}
library(IceR)
```

Next, we load IceR and MaxQ data using the respective functions supplyed by the IceR package. A file-choose box will open to specify the location of the respective data. 

```{r}
IceR <- load_Requant_data()
MaxQ <- load_MaxQ_data()
```
Add annotation information to both data sets. The first two samples correspond to two replicates of 3 % E.coli lysate spiked into constant human background. The second two samples correspond to 9 % E.coli spike-in conditions.
```{r}
anno <- data.frame(Spike = c(3,3,9,9))
IceR <- add_annotations(IceR,anno)
MaxQ <- add_annotations(MaxQ,anno)
```
Next, perform median normalization of protein- and peptide-level data based on background proteome (human proteins).
```{r}
IceR$Protein_level$Quant_data_norm <- normalize_data(IceR$Protein_level$Quant_data,method = "median",main = "IceR - Protein-level data",norm_on_subset = which(IceR$Protein_level$Meta_data$Organism == "Homo sapiens"))
IceR$Peptide_level$Quant_data_norm <- normalize_data(IceR$Peptide_level$Quant_data,method = "median",main = "IceR - Peptide-level data",norm_on_subset = which(IceR$Peptide_level$Meta_data$Organism == "Homo sapiens"))
MaxQ$Protein_level$Quant_data_norm <- normalize_data(MaxQ$Protein_level$Quant_data,method = "median",main = "MaxQ - Protein-level data",norm_on_subset = which(MaxQ$Protein_level$Meta_data$Organism == "Homo sapiens"))
MaxQ$Peptide_level$Quant_data_norm <- normalize_data(MaxQ$Peptide_level$Quant_data,method = "median",main = "MaxQ - Peptide-level data",norm_on_subset = which(MaxQ$Peptide_level$Meta_data$Organism == "Homo sapiens"))
```
<p align="center"> 
<img src="images/Normalization.jpg" style="width: 100%; height: 100%"/>
</p>

We continue with a look on general numbers per data set.
```{r}
MaxQ <- determine_general_numbers(MaxQ)
IceR <- determine_general_numbers(IceR)

compare_general_numbers(list(MaxQ=MaxQ,IceR=IceR),colors = c("darkgrey","chocolate2"))
```
<p align="center"> 
<img src="images/General_numbers.jpg" style="width: 100%; height: 100%"/>
</p>

We see a clear reduction of missing values.

Check CVs of quantification in both data sets.
```{r}
plot_accuracy(list(MaxQ=MaxQ,IceR=IceR),inset = c(0,0),Legendpos = "topright",colors = c("chocolate2","darkgrey"))
```
<p align="center"> 
<img src="images/CV_quant.jpg" style="width: 100%; height: 100%"/>
</p>

We see similar CVs in MaxQuant and IceR data, however, more data points are available especially on peptide level in case of IceR data.

Finally, we perform differential expression (DE) analysis. Perform DE on protein-level in case of MaxQ data and on peptide-level (using PECA) in case of IceR data to make use of the highly increased amount of available data for IceR results. As IceR robustly infers protein abundances e.g. using the MaxLFQ algorithm, DE analyses could be of course also performed on protein-level.

```{r}
DE_MaxQ <- LIMMA_analysis(MaxQ$Protein_level$Quant_data_norm,assignments = MaxQ$Annotations$Spike,contrast = "9_vs_3")
DE_IceR <- PECA_analysis(IceR$Peptide_level$Quant_data_norm,ids = IceR$Peptide_level$Meta_data$Gene_name,anno = IceR$Annotations$Spike,group1_name = "9",group2_name = "3")
```

Now plot Volcano plots with significance cutoffs (dashed black lines) at adj.pval < 0.05 and abs. logfc >= 1. Shape data points based on organism. Indicate true ratio (dashed red line). Add barchart representing number of detected true positives and false positives.

```{r}
plot(DE_MaxQ$logFC,-log10(DE_MaxQ$P.Value),xlab="Ratio, log2",ylab="-log10 pval",main="MaxQ - 9 % vs 3 % spike-in",col=ifelse(abs(DE_MaxQ$logFC)>=1 & DE_MaxQ$adj.P.Val < 0.05,"red","black"),pch=ifelse(rownames(DE_MaxQ) != toupper(rownames(DE_MaxQ)),19,1),xlim=c(-3,3))
abline(h=-log10(0.05),lty=2)
abline(v=1,lty=2)
abline(v=-1,lty=2)
abline(v=log2(9/3),lty=3,col="red")
legend("top",legend = c("significant","not significant","Spiked","Background"),title = "Legend",col=c("red","black","black","black"),pch=c(15,15,19,1))

plot(DE_IceR$logFC,-log10(DE_IceR$P.Value),xlab="Ratio, log2",ylab="-log10 pval",main="IceR - 9 % vs 3 % spike-in",col=ifelse(abs(DE_IceR$logFC)>=1 & DE_IceR$adj.P.Val < 0.05,"red","black"),pch=ifelse(rownames(DE_IceR) != toupper(rownames(DE_IceR)),19,1),xlim=c(-4,4))
abline(h=-log10(0.05),lty=2)
abline(v=1,lty=2)
abline(v=-1,lty=2)
abline(v=log2(9/3),lty=3,col="red")
legend("top",legend = c("significant","not significant","Spiked","Background"),title = "Legend",col=c("red","black","black","black"),pch=c(15,15,19,1))

MaxQ_TP <- length(which(DE_MaxQ$logFC>=1 & DE_MaxQ$adj.P.Val < 0.05 & rownames(DE_MaxQ) != toupper(rownames(DE_MaxQ))))
MaxQ_FP <- length(which(abs(DE_MaxQ$logFC)>=1 & DE_MaxQ$adj.P.Val < 0.05 & rownames(DE_MaxQ) == toupper(rownames(DE_MaxQ))))
IceR_TP <- length(which(DE_IceR$logFC>=1 & DE_IceR$adj.P.Val < 0.05 & rownames(DE_IceR) != toupper(rownames(DE_IceR))))
IceR_FP <- length(which(abs(DE_IceR$logFC)>=1 & DE_IceR$adj.P.Val < 0.05 & rownames(DE_IceR) == toupper(rownames(DE_IceR))))
plot_Data <- data.frame(MaxQ=c(MaxQ_TP,MaxQ_FP),IceR=c(IceR_TP,IceR_FP))
p <- Barplotsstacked(plot_Data,AvgLine = F,col=c("lightblue","grey"),margins = c(4,4,4,12),Legends = c("true positives","false positives"),Legendpos = "top",inset = c(0,0),main="Comparison of DE results",ylab="Count")
```

<p align="center"> 
<img src="images/DE_results.jpg" style="width: 100%; height: 100%"/>
</p>

We see an increase in true positives by +17 % in case of the IceR analysis. Number of false positives is comparable.

