#!/bin/bash

(grep -o '^[^#]*' $1) |\
awk 'BEGIN{{OFS="\t"}}{{ 
    split($8, a, ";")
    for (i = 1; i <= length(a); i++) {
        
        if (a[i] ~ /^SEQ=/) {
            split(a[i], seq, "=");
            seq_value = seq[2];
        }
        if (a[i] ~ /^SVTYPE=/) {
            split(a[i], svtype, "=");
            svtype_value = svtype[2];
        }
        if (a[i] ~ /^END=/) {
            split(a[i], end, "=");
            end_value = end[2];
        }
    }
  }}
  {{if (end_value==NULL){
    end_value=$2
  }}}
  {{print $1,($2-1),end_value,svtype_value,seq_value}}
  {{end_value=""
  svtype_value=""
  seq_value=""}}' 


