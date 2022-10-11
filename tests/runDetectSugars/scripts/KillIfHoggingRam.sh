#!/bin/sh

# $1 is process name

if [ $# -ne 1 ];
then
    echo "Invalid number of arguments"
    exit 0
fi
while true;
do
    pgrep "$1" | while read -r procId;
    do
       SIZE=$(pmap $procId | grep total | grep -o "[0-9]*")
       SIZE=${SIZE%%K*}
       SIZEMB=$((SIZE/1024))
       # echo "Process id = $procId Size = $SIZEMB MB"
       if [ $SIZEMB -gt 10000 ]; then
          sleep 120
	  COMMAND=$(ps -o cmd fp $procId)
	  CMD="CMD"
          if [ "$COMMAND" != "$CMD" ]; then
	    kill -9 "$procId"
            echo "Killed $COMMAND">>LoopKilled.txt
            exit 0
	  fi
       fi
   done
   sleep 10
done
