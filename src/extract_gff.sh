#!/bin/bash -l

if [[ ! -d $1 ]]; then
  mkdir -p tem_gf
  mkdir -p $1
  tar -xf $2 -C tem_gf
  mv tem_gf/*/*.gz $1/.
  rm -r tem_gf
fi
