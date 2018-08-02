#!/bin/sh
cd `git rev-parse --show-toplevel`
if [ 0 -lt `git status -s algorithms.math | wc -l` ]
then
    echo "You modified algorithms.math so I install the required python packages to compile them. If you want to prevent this, execute git checkout algorithms.math."
    python -m pip install --user -r requirements.txt
    python ./parser/math_parser.py
fi
