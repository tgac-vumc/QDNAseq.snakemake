
To run this pipeline Snakemake is required.

for easy installation you need (mini)conda.

Miniconda installation from folder where you want to install miniconda:

```
cd </path/to/files/dir/>
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

follow the instructions of the installation process, give the location where you want miniconda to be installed and answer YES to add miniconda to your path.

go to the directory where the analysis need to be performed

```
cd </path/to/analysis/dir>
git clone https://github.com/tgac-vumc/QDNAseq.snakemake/
cd QDNAseq.snakemake
git checkout Tjitske
```

install mamba as drop-in replacement for conda with mamba's improved installation-performance:
```
conda install -c conda-forge mamba
```

install the environment using mamba:
```
mamba env create --name  QDNAseq-snakemake --file environment.yaml
```

activate the environment by:
```
conda activate QDNAseq-snakemake
```

Not all packages required by R can be installed by conda, these need to be installed manually from R. first download CGHtest in the QDNAseq.snakemake directory:
source code at: http://www.few.vu.nl/~mavdwiel/CGHtest.html

```
cd QDNAseq.snakemake
wget http://www.few.vu.nl/~mavdwiel/CGHtest/CGHtest_1.1.tar.gz
```

Then run the following R-script in the terminal to install the non-conda r-dependencies:

```
Rscript r-dependencies.R
```

go to analysis dir and prepare analysis by copy or create links to fastq.gz files:

```
cd </path/to/analysis/dir>

mkdir fastq
cd fastq

```
to link a single file:
```
ln -s <path/to/file>
```

to link a files from a folder:

```
for file in <path/to/fastq/files>/*.fastq.gz
do ln -s $file
done
```

go to QDNAseq.snakemake directory and optionally change configfile, than start snakemake.
One of the options in the configfile is dewaving, if set to true QNDAseq objects will be dewaved before segmentation. 

```
cd ../QDNAseq.snakemake

snakemake

```
Useful snakemake options
```
-j , --cores, --jobs : Use at most N cores in parallel (default: 1). If N is omitted, the limit is set to the number of available cores.
-n , --dryrun : Do not execute anything. but show rules which are planned to be performed.
-k , --keep-going : Go on with independent jobs if a job fails.
-f , --force : Force the execution of the selected target or the first rule regardless of already created output.
-U , --until : Runs the pipeline until it reaches the specified rules or files. Only runs jobs that are dependencies of the specified rule or files, does not run sibling DAGs.
-T , --timestamp : Add a timestamp to all logging output
```
for all options go to http://snakemake.readthedocs.io/en/stable/executable.html#all-options
