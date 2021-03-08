[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/BiodataAnalysisGroup/godel-numbering/main?filepath=Biological-Sequences-and-Godel-numbers.ipynb)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Godel numbers and DNA

This material is made available under the [MIT license](https://opensource.org/licenses/MIT). Please see [LICENSE](LICENSE.md) for more details.

## Citation

- (_arXiv preprint_) Argyris Nicolaidis and Fotis Psomopoulos, "DNA coding and Gödel numbering", arXiv:1909.13574, 2019 (_[link](https://arxiv.org/abs/1909.13574)_)

## Structure
godel-numbering repository consists of four folders:
- `data`: This is where all the input files are stored.
- `R`: This is where all the R scripts are stored.
- `plots`: This is where all the output files are stored. This folder is created automatically while executing the project. There are also some examples of outputs in zip format.
- `lit`: some literature stuff


## Getting started
### Dependencies
Execute the following line to install the required packages:
- from CRAN:

```
install.packages(c("seqinr", "MASS", "ggplot2", "ggpubr", "stringr", "combinat" , "microseq"))
```

or run `install.R` script.

### Setting up
The project can be downloaded using git:
```
git clone https://github.com/BiodataAnalysisGroup/godel-numbering.git
```

### Extra files
- Download **S1_L001** and **S10_L001** samples from [E-MTAB-6962 - RNA-seq dataset](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6962/samples/), which correspond to ERR2681749_1.fastq.gz and ERR2681763_1.fastq.gz files respectively. 
- **Store them in the same directoy.**.

### Running the project
The framework consists of six scripts:
- ```script_mirs.R```
- ```qc.R``` 
- ```diff_analysis.R```
- ```fa.R```
- ```save_as_pdf.R``` 
- ```save_image.R```

In order to run the project:
1. Set the [R](https://github.com/BiodataAnalysisGroup/miRNAtool/tree/main/R) folder as your working directory
2. Place the required input files inside the `data` directory and set the parameters inside the .yaml file. See instructions [here](https://github.com/BiodataAnalysisGroup/miRNAtool/tree/main/data). 
3. Use the following command:
```
source("script_mirs.R")
```
