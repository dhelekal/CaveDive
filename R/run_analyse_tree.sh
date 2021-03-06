#!/bin/bash
cp analyse_tree.rmd $1;
cd $1;
R -e "rmarkdown::render('analyse_tree.rmd')";