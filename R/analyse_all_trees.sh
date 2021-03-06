#!/bin/bash
ls ./trees | parallel -j 20 --verbose --halon-error 2 eval sh "run_analyse_tree.sh ./trees/{}"