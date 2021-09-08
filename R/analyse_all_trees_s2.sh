#!/bin/bash
ls $1 | parallel -j 8 sh run_analyse_tree_s2.sh $1/{}