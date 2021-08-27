#!/bin/sh

# path to inner shell script
prog=sh/send-sortcons-commands.sh

mem=2G
tm=72::

f_log=sendbkgd.sortcons.log
qsub -l mem=${mem},time=${tm} -cwd -j y -o $f_log $prog
