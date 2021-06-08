#!/bin/bash
ls $1 | parallel -j 6 sh run_analyse_tree.sh $1/{}