#!/bin/sh

DIR=$HOME/git_repos/projects/src

nohup $DIR/dnu      flat_file.dat &
nohup $DIR/ptc/fmap flat_file.dat 1 &
nohup $DIR/ptc/fmap flat_file.dat 2 &

wait
