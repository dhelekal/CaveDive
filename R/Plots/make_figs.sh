#!/bin/bash
ls | grep .*.R | xargs -L 1 R -f