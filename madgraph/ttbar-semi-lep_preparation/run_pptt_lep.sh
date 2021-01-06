#!/bin/bash
source ~/.bash_profile
for i in {1..9};do
	if [ $i == 1 ]
       	then
		echo $i 
        	nohup ./bin/mg5_aMC pptt_lep_first_round.txt >> /home/david/log_lep.txt 
		#nohup python3 parse_uproot_release.py /home/david/pptt_lep/Events/run_0$i/tag_1_delphes_events.root /home/david/mass_generation/event_record_0$i.h5 >> /home/david/log_lep.txt
		#rm /home/david/pptt_lep/Events/run_0$i/*
	else
 		echo $i
		nohup ./bin/mg5_aMC pptt_lep.txt >> /home/david/log_lep.txt
		#nohup python3 parse_uproot_release.py /home/david/pptt_lep/Events/run_0$i/tag_1_delphes_events.root /home/david/mass_generation/event_record_0$i.h5 >> /home/david/log_lep.txt
		#rm /home/david/pptt_lep/Events/run_0$i/*
	fi 		
done
for i in {10..11};do
	echo $i
       	nohup ./bin/mg5_aMC pptt_lep.txt >> /home/david/log_lep.txt
        #nohup python3 parse_uproot_release.py /home/david/pptt_lep/Events/run_$i/tag_1_delphes_events.root /home/david/mass_generation/event_record_$i.h5 >> /home/david/log_lep.txt
	#rm /home/david/pptt_lep/Events/run_$i/*
done
