#!/bin/sh

DIR=$HOME/git_repos/projects/src

echo $DIR

$DIR/dnu      flat_file.dat &
$DIR/ptc/fmap flat_file.dat 1 &
$DIR/ptc/fmap flat_file.dat 2 &

wait
