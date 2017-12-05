#!/bin/sh

mkdir bin
mkdir bin/tmp
# rm -r source/tmp/*

printf "OSHUN version control: \n\n" > version.txt 
git log --pretty=format:'%h : %s' --graph >> version.txt
printf "\n\n (PICKSC made this)" >> version.txt

cd source
cp makefiles/makefile makefile
make
#make debug
