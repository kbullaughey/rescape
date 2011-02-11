#!/bin/bash

../bin/regevscape --model=segal --steps=200000 --length=20 --seed=5 \
  --pwm=pwm_pair.pwms --tfs=2 --lambda=-1,2,-3 --gamma=2,-0.5,-0.5,2 \
  --ffunc=f4 --coef=0.6,0.041,0.786,0.999 \
  --condition=0.02,0.3 --condition=0.3,0.3 --condition=0.3,0.02 \
  --neighbors=1 --lstat=expression,1 --lstat=sequence,1
