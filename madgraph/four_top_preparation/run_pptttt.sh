#!/bin/bash
source ~/.bash_profile
for i in {1..9};do
	if [ $i == 1 ]
       	then
		echo $i
        	nohup ./bin/mg5_aMC pptttt_first_round.txt >> /home/david/log.txt
            nohup python3 /home/david/workplace/HEP-jet-assignment/analysis_script/main.py -u parse -m four_top -i /home/david/pptttt/Events/run_0$i/tag_1_delphes_events.root -o event_record_four_top_$i.h5 -s 1 
            rm /home/david/pptttt/Events/run_0$i/*
	else
 		echo $i
		nohup ./bin/mg5_aMC pptttt.txt >> /home/david/log.txt
        nohup python3 /home/david/workplace/HEP-jet-assignment/analysis_script/main.py -u parse -m four_top -i /home/david/pptttt/Events/run_0$i/tag_1_delphes_events.root -o event_record_four_top_$i.h5 -s 1 
        rm /home/david/pptttt/Events/run_0$i/*
	fi
done
for i in {10..11};do
	echo $i
       	nohup ./bin/mg5_aMC pptttt.txt >> /home/david/log.txt
        nohup python3 /home/david/workplace/HEP-jet-assignment/analysis_script/main.py -u parse -m four_top -i /home/david/pptttt/Events/run_$i/tag_1_delphes_events.root -o event_record_four_top_$i.h5 -s 1 
        rm /home/david/pptttt/Events/run_$i/*
done
