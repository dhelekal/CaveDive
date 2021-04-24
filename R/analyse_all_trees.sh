#!/bin/bash
ls ./Paper/Trees | parallel -j 8 sh run_analyse_tree.sh ./Paper/Trees/{}