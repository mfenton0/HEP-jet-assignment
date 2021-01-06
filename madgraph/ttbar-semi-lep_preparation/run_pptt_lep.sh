#!/bin/bash
source ~/.bash_profile
ANALYSIS_SCRIPT_PATH='/home/david/workplace/HEP-jet-assignment/analysis_script'
ROOT_FILE_PATH='/home/david/pptt_lep/Events'
LOG_FILE_PATH='/home/david'
OUTPUT_FILE_PATH='/home/david'
echo "Path of script: $ANALYSIS_SCRIPT_PATH"
echo "path of root file: $ROOT_FILE_PATH"
echo "path of log file: $LOG_FILE_PATH"
echo "path of output file: $OUTPUT_FILE_PATH"
for i in {1..9};do
	if [ $i == 1 ]
       	then
		echo $i 
        	nohup ./bin/mg5_aMC pptt_lep_first_round.txt >> $LOG_FILE_PATH/log_lep.txt 
		nohup python3 $ANALYSIS_SCRIPT_PATH/main.py -u parse -m ttbar -i $ROOT_FILE_PATH/run_0$i/tag_1_delphes_events.root -o $OUTPUT_FILE_PATH/event_record_top_lep_$i.h5 -s 1 >> $LOG_FILE_PATH/log_lep.txt
		rm $ROOT_FILE_PATH/run_0$i/*
	else
 		echo $i
		nohup ./bin/mg5_aMC pptt_lep.txt >> $LOG_FILE_PATH/log_lep.txt
		nohup python3 $ANALYSIS_SCRIPT_PATH/main.py -u parse -m ttbar -i $ROOT_FILE_PATH/run_0$i/tag_1_delphes_events.root -o $OUTPUT_FILE_PATH/event_record_top_lep_$i.h5 -s 1 >> $LOG_FILE_PATH/log_lep.txt
		rm $ROOT_FILE_PATH/run_0$i/*
	fi 		
done
for i in {10..11};do
	echo $i
       	nohup ./bin/mg5_aMC pptt_lep.txt >> $LOG_FILE_PATH/log_lep.txt
	nohup python3 $ANALYSIS_SCRIPT_PATH/main.py -u parse -m ttbar -i $ROOT_FILE_PATH/run_$i/tag_1_delphes_events.root -o $OUTPUT_FILE_PATH/event_record_top_lep_$i.h5 -s 1 >> $LOG_FILE_PATH/log_lep.txt
	rm $ROOT_FILE_PATH/run_$i/*
done