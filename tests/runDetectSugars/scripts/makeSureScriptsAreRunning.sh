#!/bin/sh

while true; do
	if ! pgrep -fx "KillIfHoggingRam.sh" > /dev/null; then
		./KillIfHoggingRam.sh detect_sugars &
	fi
	if ! pgrep -fx "ontology.sh" > /dev/null; then
		nohup ./controlCPU.sh /programs/repos/PDB/all/ 2 >> 12-7-20.log 2>&1
	fi
	sleep 60
done

