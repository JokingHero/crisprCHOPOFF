#!/bin/bash

#Hakon Tjeldnes 09/04/2024
# Script to install Julia 1.8.5 and CHOPOFF.jl

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
echo CHOPOFF=~/bin/CHOPOFF.jl/build/bin/CHOPOFF >> ~/.Renviron
# <- To allow RSTUDIO to see the CHOPOFF path, then restart RSTUDIO
echo "INSTALL DONE"
