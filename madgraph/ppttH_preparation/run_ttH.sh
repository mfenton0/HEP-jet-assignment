#!/bin/bash
source ~/.bash_profile
for i in {1..9};do
	if [ $i == 1 ]
       	then
		echo $i 
        	nohup ./bin/mg5_aMC ttH_first_run.txt >> /home/david/log_ttH.txt 
		nohup python3 parse_uproot_release_ttH.py /home/david/ppttH/Events/run_0$i/tag_1_delphes_events.root /home/david/mass_generation/event_record_0$i.h5 >> /home/david/log.txt
		rm /home/david/ppttH/Events/run_0$i/*
	else
 		echo $i
		nohup ./bin/mg5_aMC ttH.txt >> /home/david/log_ttH.txt
		nohup python3 parse_uproot_release_ttH.py /home/david/pptt/Events/run_0$i/tag_1_delphes_events.root /home/david/mass_generation/event_record_0$i.h5 >> /home/david/log.txt
		rm /home/david/ppttH/Events/run_0$i/*
	fi 		
done
for i in {10..100};do
	echo $i
	nohup ./bin/mg5_aMC ttH.txt >> /home/david/log_ttH.txt
	nohup python3 parse_uproot_release_ttH.py /home/david/pptt/Events/run_$i/tag_1_delphes_events.root /home/david/mass_generation/event_record_$i.h5 >> /home/david/log.txt
	rm /home/david/ppttH/Events/run_$i/*
done
