#! /bin/bash

if compgen -G "*.log" > /dev/null; then
    mv *.log log/
fi
if compgen -G "*.o*" > /dev/null; then
    mv *.o* log/
fi
if compgen -G "*.po*" > /dev/null; then
    mv *.po* log/
fi
if compgen -G "*report.*" > /dev/null; then
    mv *report.* report/
fi
if compgen -G "*.params" > /dev/null; then
    mv *.params report/
fi
if compgen -G "*.pdf" > /dev/null; then
    mv *.pdf plots/
fi

# http://stackoverflow.com/a/6364244
if compgen -G "log/fail.*" > /dev/null; then
    for f in log/fail.*; do
        if [ -e "$f" ] ; then
    	echo "Some jobs failed; check log/ directory"
        fi
        exit 1
    done
fi
