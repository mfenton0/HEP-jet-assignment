#!/bin/bash
source ~/.bash_profile
for i in {1..9};do
	if [ $i == 1 ]
       	then
		echo $i 
        	nohup ./bin/mg5_aMC pptt_lep_first_round.txt >> /home/david/log_lep.txt 
		nohup python3 /home/david/workplace/HEP-jet-assignment/analysis_script/main.py -u parse -m ttbar -i /home/david/pptt_lep/Events/run_0$i/tag_1_delphes_events.root -o event_record_top_lep_$i.h5 -s 1 >> /home/david/log_lep.txt
		rm /home/david/pptt_lep/Events/run_0$i/*
	else
 		echo $i
		nohup ./bin/mg5_aMC pptt_lep.txt >> /home/david/log_lep.txt
		nohup python3 /home/david/workplace/HEP-jet-assignment/analysis_script/main.py -u parse -m ttbar -i /home/david/pptt_lep/Events/run_0$i/tag_1_delphes_events.root -o event_record_top_lep_$i.h5 -s 1 >> /home/david/log_lep.txt
		rm /home/david/pptt_lep/Events/run_0$i/*
	fi 		
done
for i in {10..11};do
	echo $i
       	nohup ./bin/mg5_aMC pptt_lep.txt >> /home/david/log_lep.txt
	nohup python3 /home/david/workplace/HEP-jet-assignment/analysis_script/main.py -u parse -m ttbar -i /home/david/pptt_lep/Events/run_$i/tag_1_delphes_events.root -o event_record_top_lep_$i.h5 -s 1 >> /home/david/log_lep.txt
	rm /home/david/pptt_lep/Events/run_$i/*
done
