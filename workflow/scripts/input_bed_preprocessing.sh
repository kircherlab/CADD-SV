#!/bin/bash

(grep -o '^[^#]*' $1) |\
awk 'BEGIN{{OFS="\t"}}{{
    print $1,$2,$3,$4,""}}'


