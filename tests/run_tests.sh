#!/bin/bash

for i in *_test.py
do
    python $i
    echo -e "\n"
done

