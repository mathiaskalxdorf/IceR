# IceR - Quantitative Label-Free Proteomics Workflow

<p align="center"> 
<img src="images/IceR.jpg" style="width: 50%; height: 50%"/>​
</p>

### Description
Label-free proteomics enables the unbiased quantification of thousands of proteins across large sample cohorts. Commonly used mass spectrometry-based proteomic workflows rely on data dependent acquisition (DDA). However, its stochastic selection of peptide features for fragmentation-based identification inevitably results in high rates of missing values, which prohibits the integration of larger cohorts as the number of recurrently detected peptides is a limiting factor. Peptide identity propagation (PIP) can mitigate this challenge, allowing to transfer sequencing information between samples. However, despite the promise of these approaches, current methods remain limited either in sensitivity or reliability and there is a lack of robust and widely applicable software. To address this, we here present IceR, an efficient and user-friendly quantification workflow introducing a hybrid PIP approach with superior quantification precision, accuracy, reliability and data completeness. IceR is available as an easy to-use R-package incorporating a graphical user interface and comprehensive quality control measures.

### Installation
The following installations are required:

 - [R](https://cran.r-project.org/bin/windows/base/) (Version 3.63 or above)
 - [Rtools](https://cran.r-project.org/bin/windows/Rtools/history.html) (select required version)
 - Optional: [ProteoWizard](http://proteowizard.sourceforge.net/download.html)
 - Optional: [RStudio](https://rstudio.com/products/rstudio/download/)
 
During installation, please keep default settings.

Next, we install the IceR package from GitHub
```r
install.packages("devtools")
devtools::install_github("mathiaskalxdorf/IceR")
```
If everything wents fine, IceR should be installed and can be used (e.g. with the GUI) with the following lines:
```r
library(IceR)
runIceR()
```
### Prerequisites
The current version of IceR requires raw MS files (from Thermo Mass Spectrometers) to be preprocessed with MaxQuant (tested for Versions 1.5.1.2, 1.6.12 and 1.6.14, versions in between should work as well). Important note: Fasta file should be parsed correctly. In this case the column "Gene names" can be found in the proteinGroups.txt. Furthermore, raw files are required in the mzXML format. If msconvert from ProteoWizard is installed, IceR triggers conversion automatically. Otherwise, the user has to place the converted files in the folder "mzXML" within the folder containing the raw files. 

### Example
An example data set is stored [here](http://proteomecentral.proteomexchange.org/cgi/GetDataset). After running ...
```r
library(IceR)
runIceR()
```
... the GUI opens.
<p align="center"> 
<img src="images/IceR_gui.jpg" style="width: 50%; height: 50%"/>​
</p>

It allows setting up the IceR run. Among others, the following parameters can be modified:

 - Paths to raw files, MaxQuant results files and IceR output folder
 - Analysis name
 - Retention (RT) and m/z feature alignment windows
 - Kernel density estimation (KDE) resolution (increasing resolution increases computation workload but also increases resolution of peak detections) 
 - Number of peaks to be stored per quantification
 - Diagnostic and statistical significance cutoffs for quantifications
 - Number of threads to be used during IceR workflow (please check on your system specifications. Should not be higher than numbers of samples to be analyzed. Number of threads should be also adapted to available RAM)

We can keep all settings at default but reduce number of threads to 4. After clicking on "Total process", the IceR workflow will start. A detailed description of the individual steps can be found in the original publication [IceR](https://pubmed.ncbi.nlm.nih.gov/)

A recent computational system should be able to complete the IceR workflow for the example data set within 2 - 3 hours.

Afterwards, some quality control plots for the alignment and quantification can be visualized by clicking on the tap "QC" and subsequently on the respective button:

<p align="center"> 
<img src="images/IceR_gui_QC_alignment.jpg" style="width: 50%; height: 50%"/>​
</p>

<p align="center"> 
<img src="images/IceR_gui_QC_quantification.jpg" style="width: 50%; height: 50%"/>​
</p>
 
