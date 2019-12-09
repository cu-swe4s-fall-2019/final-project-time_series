#!/bin/sh
export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.1.0:$PATH
output=`basename $1 _promoters.bed`

ame -o $output --control --shuffle-- $1 $2
