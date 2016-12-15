#!/bin/sh

DIR=$HOME/git_repos/projects/src

echo $DIR

NOHUP $DIR/leac param.dat 1 >& leac_1.log &
NOHUP $DIR/leac param.dat 2 >& leac_2.log &

wait
