#!/bin/bash
source ~/.bash_profile
for i in {1..9};do
	if [ $i == 1 ]
       	then
		echo $i 
        	nohup ./bin/mg5_aMC ttH_first_round.txt >> /home/david/log_ttH.txt 
		nohup python3 /home/david/workplace/HEP-jet-assignment/analysis_script/main.py -u parse -m ttH -i /home/david/ppttH/Events/run_0$i/tag_1_delphes_events.root -o event_record_ttH_$i.h5 -s 1
		rm /home/david/ppttH/Events/run_0$i/*
	else
 		echo $i
		nohup ./bin/mg5_aMC ttH.txt >> /home/david/log_ttH.txt
		nohup python3 /home/david/workplace/HEP-jet-assignment/analysis_script/main.py -u parse -m ttH -i /home/david/ppttH/Events/run_0$i/tag_1_delphes_events.root -o event_record_ttH_$i.h5 -s 1
		rm /home/david/ppttH/Events/run_0$i/*
	fi 		
done
for i in {10..100};do
	echo $i
	nohup ./bin/mg5_aMC ttH.txt >> /home/david/log_ttH.txt
	nohup python3 /home/david/workplace/HEP-jet-assignment/analysis_script/main.py -u parse -m ttH -i /home/david/ppttH/Events/run_$i/tag_1_delphes_events.root -o event_record_ttH_$i.h5 -s 1
	rm /home/david/ppttH/Events/run_$i/*
done
