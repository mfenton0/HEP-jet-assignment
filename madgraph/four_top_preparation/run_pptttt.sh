#!/bin/bash
source ~/.bash_profile
for i in {1..9};do
	if [ $i == 1 ]
       	then
		echo $i
        	nohup ./bin/mg5_aMC pptttt_first_round.txt >> /home/david/log.txt
	else
 		echo $i
		nohup ./bin/mg5_aMC pptttt.txt >> /home/david/log.txt
	fi
done
for i in {10..11};do
	echo $i
       	nohup ./bin/mg5_aMC pptttt.txt >> /home/david/log.txt
done
