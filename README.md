crisprCHOPOFF: R wrapper for CHOPOFF, a CRISPR OFF target search algorithm
==============================================================================


This package is under heavy development to include more features.

#### About


Efficient search of CRISPR OFF targets


#### Installation

Package is available here on github (Experimental branch, R version >= 4.1.0)
```r
if (!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes")
remotes::install_github("JokingHero/crisprCHOPOFF)
```  

#### System requirements
OS support: Unix (linux & mac)
For windows, run through WSL (windows subsystem for linux)

Normally, first time running build_index(), will install the backend for you, 
through the function install_CHOPOFF()
```r
# The script ran, can be found from:
system.file("bash_script", "install_julia_and_CHOPOFF.sh", package = "crisprCHOPOFF")

install_CHOPOFF()
```  

For direct install into ~/bin, copy and run code below, it will install:
- Julia 1.8.5
- CHOPOFF.jl
The script content is this:
```sh
# In Shell
# Get julia 1.8.5 (required version)
mkdir -p ~/bin/julia
cd ~/bin/julia
wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.5-linux-x86_64.tar.gz
tar -xvzf julia-1.8.5-linux-x86_64.tar.gz
# Update PATH and source .bashrc
echo 'export PATH="$PATH:~/bin/julia/julia-1.8.5/bin"' >> ~/.bashrc
echo 'export JULIA_NUM_THREADS=11' >> ~/.bashrc
. ~/.bashrc # resource shell
# Clone julia CHOPOFF library (select folder where you have your github clones)
cd ~/bin # update path to your path if needed
git clone git@github.com:JokingHero/CHOPOFF.jl.git
cd CHOPOFF.jl
./build_standalone.sh
# In .bashrc (update to your custom path!)
echo 'export PATH="$PATH:~/bin/CHOPOFF.jl/build/bin/"' >> ~/.bashrc
# In Shell
. ~/.bashrc # resource shell
echo CHOPOFF=~/bin/CHOPOFF.jl/build/bin/CHOPOFF >> ~/.Renviron
# <- To allow RSTUDIO to see the CHOPOFF path, then restart RSTUDIO
```  

CHOPOFF should now be ready to use in the R wrapper

#### More information

After installation run:
```r
library(crisprCHOPOFF)

# Build a genome index to search against
?build_index

# Search the index with selected guides
?search_index
```  

Quick walk-through:

```r
name <- "CAS9"
genome <- system.file("extdata/sample_genome", "semirandom.fa", package = "crisprCHOPOFF")
## Note: a fasta index ".fai" file must exist in directory of genome.
# You can make it with:
#if (!file.exists(paste0(genome, ".fai"))) {
#Rsamtools::indexFa(genome)
#}
out_dir_index <- file.path(tempdir(), "CHOPOFF_sample_genome")
build_index(name, genome, out_dir_index, validate = FALSE, distance = 2)

# Now search some guides
guides <- system.file("extdata/sample_genome", "guides.txt", package = "crisprCHOPOFF")
# Quick preview in guides:
guide_candidates <- read.table(guides, col.names = "guides")
unique(nchar(unlist(guide_candidates))) # Unique lengths of guides
guide_hits <- search_index(guides, out_dir_index, validate = FALSE, distance = 2)
guide_hits_table <- read.table(guide_hits, sep = ",", header = TRUE)
# use data.table::fread for reading in large list
# Subset to 0 distance hits
dist0 <- guide_hits_table[guide_hits_table$distance == 0,]
head(dist0)
# Which chromosomes is a specific guide found on with 0 distance hits?
unique(dist0[dist0$guide == "TCCGGCCTGGTTATCGAAGG",]$chromosome) # 2 chromosomes
```

Please read Bioconductor vignettes for detailed tutorials and examples.

#### Feedback

Please feel free to provide feedback or desired functionality by creating a new issue on our github page.
