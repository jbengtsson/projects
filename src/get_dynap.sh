#!/bin/sh

queue="prime_bd.q"

t1="10:00:00"
t2="24:00:00"

#dir=`pwd`
dir=$HOME/git_repos/projects/src

\rm dynap.cmd.o*

qsub -l s_rt=$t1 -l h_rt=$t2 -q $queue $dir/dynap.cmd
