#!/bin/sh

set echo
source ~/.bashrc

python convGRUNet_train.py | tee convGRUNet_train.log