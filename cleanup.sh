#! /bin/bash

mv *.log log/
mv *.o* log/
mv *.po* log/
mv fail.* log/

mv *report* report/
mv *.params report/

# http://stackoverflow.com/a/6364244
for f in log/fail.*; do
    if [ -e "$f" ] ; then
	echo "Some jobs failed; check log/ directory"
    fi
    exit 1
done
