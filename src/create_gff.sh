#!/bin/bash -l

mkdir -p $2
tar -xf $1 -C $2
ls $2/*/*.gz | while read file
              do lines=$(gunzip -cd $file | grep "$3" | grep "$4")
                name=$(basename $file)
                echo -e "$name\t$lines" >> $5
              done
rm -r $2
