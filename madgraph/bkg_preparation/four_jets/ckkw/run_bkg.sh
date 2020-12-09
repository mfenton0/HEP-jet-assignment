#!/bin/bash
source ~/.bash_profile
mkdir /home/david/four_jets_ckkw_root

for i in {1..10};do
	if [ $i == 1 ]
       	then
		echo $i
        	nohup ./bin/mg5_aMC ppjjjj_ckkw_first_run.txt >> /home/david/log_four_ckkw_jets.txt
		mv /home/david/ppjjjj/Events/run_0$i/tag_1_delphes_events.root /home/david/four_jets_ckkw_root/tag_1_delphes_events_$i.root
		rm /home/david/bkg/Events/run_0$i/*
	else
 		echo $i
		nohup ./bin/mg5_aMC ppjjjj_ckkw.txt >> /home/david/log_four_ckkw_jets.txt
		mv /home/david/ppjjjj/Events/run_0$i/tag_1_delphes_events.root /home/david/four_jets_ckkw_root/tag_1_delphes_events_$i.root
                rm /home/david/ppjjjj/Events/run_0$i/*
	fi
done
for i in {10..30};do
	echo $i
       	nohup ./bin/mg5_aMC ppjjjj_ckkw.txt >> /home/david/log_four_ckkw_jets.txt
	mv /home/david/ppjjjj/Events/run_$i/tag_1_delphes_events.root /home/david/four_jets_ckkw_root/tag_1_delphes_events_$i.root
        rm /home/david/ppjjjj/Events/run_$i/*
done
