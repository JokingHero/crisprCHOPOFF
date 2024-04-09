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
export PATH="$PATH:~/bin/CHOPOFF.jl/build/bin/"
# In Shell
. ~/.bashrc # resource shell
echo CHOPOFF=/home/roler/Desktop/forks/CHOPOFF.jl/build/bin/CHOPOFF >> ~/.Renviron
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
Please read Bioconductor vignettes for detailed tutorials and examples.

#### Feedback

Please feel free to provide feedback or desired functionality by creating a new issue on our github page.
