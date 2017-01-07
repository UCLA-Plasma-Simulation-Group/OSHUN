#!/bin/sh

mkdir bin
mkdir source/tmp
# rm -r source/tmp/*

printf "OSHUN version control: \n\n" > version.txt 
git log --pretty=format:'%h : %s' --graph >> version.txt
printf "\n\n (PICKSC made this)" >> version.txt

cd source
make
#make debug
