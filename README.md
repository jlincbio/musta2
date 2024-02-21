# MuSTA2: Multi-Sample Transcript Assembly, Version 2



Current release: 
[![DOI](https://zenodo.org/badge/748499811.svg)](https://zenodo.org/doi/10.5281/zenodo.10685676)

This is a pipeline and associated R codes which take long-read isoform sequencing results as input and genarate an assembled transcriptome.
Documentation for the original version of MuSTA as well as a quick-start guide can be found at https://musta.readthedocs.io.

A preview for MuSTA2 was unveiled during the annual meeting of the Japanese Cancer Association in 2022: [J. Lin, T. Morinaga, M. Kawazu. Characterization of cancer-associated transcriptional splicing variants from massive parallel sequencing by MuSTA (abstract). In: Abstracts of the 81st Annual Meeting of the Japanese Cancer Association; 2022 29 September - 1 October; Yokohama, Japan. _Cancer Sci._ 2023: **114** (Suppl 1); 1785](https://doi.org/10.1111/cas.15742)

Last updated: 1/18/2024

## Prerequisites
* Python 3.7 via a Miniconda or Anaconda environment<sup>1</sup>
* Perl 5.18 or later
* R 4.1 or later
* [SQANTI3](https://github.com/ConesaLab/SQANTI3.git)<sup>2</sup>
* [GMAP](http://research-pub.gene.com/gmap/)<sup>3</sup>
* [Minimap2](https://github.com/lh3/minimap2)
* [Salmon](https://github.com/COMBINE-lab/salmon)
* [Samtools](https://github.com/samtools/samtools)
* [STAR](https://github.com/alexdobin/STAR)
* [DTUrtle](https://github.com/TobiTekath/DTUrtle)

<sup>1</sup>The supplied Conda environment configuration file will automatically resolve and install all dependencies aside from DTUrtle.

<sup>2</sup>SQANTI3 versions after Jun. 14 2022 should theoretically work with MuSTA2.

<sup>3</sup>Only GMAP release dated 2017.11.15 is supported as newer versions contain fundamentally major rewrites that are hugely incompatible.

## Installation

MuSTA2 requires the presence of a Miniconda or Anaconda environment; this is mostly to accomodate companion tools such as SQANTI 3; if you alread have a SQANTI 3 python environment present then MuSTA can be installed within that environment as well.

```
mkdir musta2
conda env create -f musta.conda_env.yml
conda activate musta2
cd musta2 && git clone https://github.com/ConesaLab/SQANTI3.git
mv SQANTI3/sqanti3*.py .
mv SQANTI3/utilities .
rm -rf SQANTI3
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred -P ./utilities/
chmod +x ./utilities/gtfToGenePred 
git clone https://github.com/Magdoll/cDNA_Cupcake.git
CONDA_BIN=`dirname $(which gmap)`
python ./cDNA_Cupcake/setup.py build
python ./cDNA_Cupcake/setup.py install
cpanm Sort::Key::Natural
```

## Usage

`musta2 -h`

## Citation
Citation for the original version of MuSTA: [Namba et. al, _Commun Biol_ **4**: 1320, 2021](https://www.nature.com/articles/s42003-021-02833-4)

