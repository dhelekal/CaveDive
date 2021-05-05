#!/bin/bash
cp analyse_tree_s2.rmd $1;
cd $1;
R -e "rmarkdown::render('analyse_tree_s2.rmd')";