#!/bin/bash
source ~/.bash_profile
for i in {1..9};do
	if [ $i == 1 ]
       	then
		echo $i 
        	nohup ./bin/mg5_aMC bkg_first_run.txt >> /home/david/log_bkg.txt 
		#nohup python3 parse_uproot_release.py /home/david/pptt/Events/run_0$i/tag_1_delphes_events.root /home/david/mass_generation/event_record_0$i.h5 >> /home/david/log.txt
		#rm /home/david/pptt/Events/run_0$i/*
	else
 		echo $i
		nohup ./bin/mg5_aMC bkg.txt >> /home/david/log_bkg.txt
		#nohup python3 parse_uproot_release.py /home/david/pptt/Events/run_0$i/tag_1_delphes_events.root /home/david/mass_generation/event_record_0$i.h5 >> /home/david/log.txt
		#rm /home/david/pptt/Events/run_0$i/*
	fi 		
done
#for i in {10..11};do
#	echo $i
#       	nohup ./bin/mg5_aMC pptt.txt >> /home/david/log.txt
#        nohup python3 parse_uproot_release.py /home/david/pptt/Events/run_$i/tag_1_delphes_events.root /home/david/mass_generation/event_record_$i.h5 >> /home/david/log.txt
	#rm /home/david/pptt/Events/run_$i/*
#done
