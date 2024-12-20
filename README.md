<p align="center">
  <img width="100%" height="100%" src="https://github.com/tgac-vumc/QDNAseq.snakemake/blob/master/DAG_all.svg">
</p>

## Installation

For the installation of this pipeline any Python install compatable Conda is required.

The pipeline itself will run on Python 3.8.5 and R 3.6.3. For exact dependencies view `environment.yaml` and `r-dependencies.R`.

### Using Conda/Mamba

for easy installation you need (Mini)Conda.

Miniconda installation from folder where you want to install Miniconda:

```
cd </path/to/files/dir/>
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

follow the instructions of the installation process, give the location where you want Miniconda to be installed and answer YES to add Miniconda to your path.

go to the directory where the analysis need to be performed

```
cd </path/to/analysis/dir>
git clone https://github.com/tgac-vumc/QDNAseq.snakemake/
cd QDNAseq.snakemake
```

install Mamba as drop-in replacement for Conda with Mamba's improved installation-performance:

```
conda install -c conda-forge mamba
```

create  the environment using Mamba:

```
mamba env create --name QDNAseq-snakemake --file environment.yaml
```

activate the environment by:

```
conda activate QDNAseq-snakemake
```

<!---
Not all packages required by R can be installed by conda, these need to be installed manually from R. first download CGHtest in the QDNAseq.snakemake directory:
source code at: http://www.few.vu.nl/~mavdwiel/CGHtest.html

```
cd QDNAseq.snakemake
wget http://www.few.vu.nl/~mavdwiel/CGHtest/CGHtest_1.1.tar.gz
```
-->

Then run the R-script r-dependencies.R in the terminal to install the non-conda R dependencies in the environment:

```
Rscript r-dependencies.R
```


### Using Singularity

Under development

<!---
If you have Singularity, you can simply pull the singularity container with all dependency resolved (in few minutes, depends on the network speed).

```
singularity pull shub://tgac-vumc/QDNAseq-snakemake
```
-->
## Preparing analysis

### Prepare the data

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

to link all files from a folder:

```
for file in <path/to/fastq/files>/*.fastq.gz
do ln -s $file
done
```

### Prepare the snakemake settings

Open the configuration file `config.yaml` to check the settings that snakemake will use and change according to your needs.
For providing service-analysis, set `setting` to `'service'`. For research purposes, set `setting` to `'research'`. For all settings set `setting` to `'all'`.

One of the options in the configfile is `dewaving`, if set to `'true'` QNDAseq objects will be dewaved before segmentation. 

These options change the rules performed in the pipeline, see the rule-graph in the next section.




## Running analysis

Make sure that snakemake is able to find the excecutive file Snakefile by performing a dry-run:

```
cd ../QDNAseq.snakemake
snakemake -n
```

Check the rules that are planned to be performed, conform the rule-graph.

An visualization of the order of rules to be performed can be viewed by running the following command and opening the DAG-file

```
snakemake --forceall --rulegraph | dot -Tsvg > DAG.svg
```

Rulegraphs for the intial settings `'service'`, `'research'` and `'all'` are commited to this repro in the files `DAG_<setting>.svg`.


When ready, run the analysis

```
snakemake
```

Useful snakemake options

`-j , --cores, --jobs` : Use at most N cores in parallel (default: 1). If N is omitted, the limit is set to the number of available cores.

`-n , --dryrun` : Do not execute anything. but show rules which are planned to be performed.

`-k , --keep-going` : Go on with independent jobs if a job fails.

`-f , --force` : Force the execution of the selected target or the first rule regardless of already created output.

`-R , --forcerun` : Force the re-execution or creation of the given rules or files. Use this option if you changed a rule and want to have all its output in your workflow updated.

`-U , --until` : Runs the pipeline until it reaches the specified rules or files. Only runs jobs that are dependencies of the specified rule or files, does not run sibling DAGs.

for all options go to https://snakemake.readthedocs.io/en/v5.31.1/executing/cli.html#all-options
