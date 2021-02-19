Bootstrap: yum
OSVersion: 7
MirrorURL: http://mirror.centos.org/centos-7/7/os/x86_64/
Include: yum

%environment    
    export PATH=/usr/local/bin:$PATH

%post
    ./environment
    
    yum -y update
    yum -qq -y install curl tar bzip2 git zip wget
    curl -sSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh
    bash /tmp/miniconda.sh -bfp /usr/local
    conda update conda -y

    git clone https://github.com/tgac-vumc/QDNAseq.snakemake/
    cd QDNAseq.snakemake

    wget http://www.few.vu.nl/~mavdwiel/CGHtest/CGHtest_1.1.tar.gz
    
    conda install mamba -c conda-forge -y
    mamba install -c conda-forge -c bioconda -c r -c defaults python=3.8.5 snakemake=5.31 pyyaml bwa samtools picard fastqc perl-bioperl perl-dbd-mysql bioconductor-qdnaseq.hg19 bioconductor-biobase bioconductor-cghcall bioconductor-cghregions bioconductor-limma bioconductor-impute r-matrixstats r-r.cache r-denstrip r-flexmix r-gtools r-mass r-devtools r-ggplot2 r-survival
    Rscript r-dependencies.R

