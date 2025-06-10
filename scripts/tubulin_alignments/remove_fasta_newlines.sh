#!/bin/bash

#script to remove newlines from fasta files where sequences are across multiple lines and put them on a single line while keeping the header
#usage: bash remove_fasta_newlines.sh input.fasta output.fasta

input=$1
output=$2

#remove newlines from fasta files where sequences are across multiple lines and put them on a single line while keeping the header
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' $input > $output