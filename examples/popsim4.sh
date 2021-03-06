#!/bin/bash

../bin/rescape --model=segal --steps=200000 --length=50 --seed=6 \
  --pwm=pwm_pair.pwms --tfs=2 --lambda=-1,2,-3 --gamma=2,-0.5,-0.5,2 \
  --ffunc=f4 --coef=0.6,0.041,0.786,0.999 \
  --condition=0.02,0.3 --condition=0.3,0.3 --condition=0.3,0.02 \
  --popsim -N 1000 -g 1000 -u 1e-06 \
  --pstat=all_alleles,250 \
  --seq=catttcggtcttgtttttggcgggaagatgttctcgactgtgtgccgcgt

