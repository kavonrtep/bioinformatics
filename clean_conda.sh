#!/bin/bash
conda install --revision 0
conda install -c  conda-forge -c bioconda fastx_toolkit
conda create -n busco -c conda-forge -c bioconda busco

