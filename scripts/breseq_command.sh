#!/bin/bash
# intended to be run by the breseq process, taking variables from the environment

if [[ -z $ADDITIONAL ]]; then
  #command="breseq $REQUIRED"
  breseq $REQUIRED
else
  #command="breseq $ADDITIONAL $REQUIRED"
  breseq $ADDITIONAL $REQUIRED
fi

#mkdir -p $(dirname $SAVETO)
#echo $command > $SAVETO