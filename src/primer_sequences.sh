#!/bin/bash -l

counter=0
grep -v 'Pos' $1 | cut -f1,7 | while read line
  do
    counter=$((counter+1))
    pos=$(echo $line | cut -d ' ' -f 1)
    check=$(echo -n "$pos" | wc -c)
    if [ $check -eq 0 ]; then
        echo "there is no primer that fullfil coverage and GC criteria" >> $2
    else
      seq=$(echo $line | cut -d ' ' -f 2)
      echo -e ">Primer_"$counter"_at_"$pos"\n"$seq"" >> $2
    fi
  done
