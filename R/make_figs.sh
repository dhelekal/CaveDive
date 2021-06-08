#!/bin/bash
cd ./Paper/Figures;
ls | grep .*.R | xargs -L 1 R -f;