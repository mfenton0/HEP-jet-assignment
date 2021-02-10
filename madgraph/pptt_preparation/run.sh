#!/bin/bash
source ~/.bash_profile

export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

ANALYSIS_SCRIPT_PATH='/home/workplace/HEP-jet-assignment/analysis_script'
ROOT_FILE_PATH='/home/david/pptt/Events'
LOG_FILE_PATH='/home/david'
echo "Path of script: $ANALYSIS_SCRIPT_PATH"
echo "path of root file: $ROOT_FILE_PATH"
echo "path of log file: $LOG_FILE_PATH"
for i in {1..9};do
	if [ $i == 1 ]
       	then
		echo $i
        ./bin/mg5_aMC pptt_first_round.txt >> $LOG_FILE_PATH/log.txt
		python3 $ANALYSIS_SCRIPT_PATH/main.py -p 1 -u parse -m ttbar -i $ROOT_FILE_PATH/run_0$i/tag_1_delphes_events.root -o /home/david/mass_generation/event_record_top_FHD_$i.npz -s 1 >> $LOG_FILE_PATH/log.txt
		rm $ROOT_FILE_PATH/run_0$i/*
	else
 		echo $i
		awk -F"=" 'BEGIN{OFS=FS} $1=="set iseed "{$2=" "$2+1}1' pptt.txt > pptt2.txt && mv -f pptt2.txt pptt.txt
		./bin/mg5_aMC pptt.txt >> $LOG_FILE_PATH/log.txt
		python3 $ANALYSIS_SCRIPT_PATH/main.py -p 1 -u parse -m ttbar -i $ROOT_FILE_PATH/run_0$i/tag_1_delphes_events.root -o /home/david/mass_generation/event_record_top_FHD_$i.npz -s 1 >> $LOG_FILE_PATH/log.txt
		rm $ROOT_FILE_PATH/run_0$i/*
	fi 		
done
for i in {10..100};do
	echo $i
	awk -F"=" 'BEGIN{OFS=FS} $1=="set iseed "{$2=" "$2+1}1' pptt.txt > pptt2.txt && mv -f pptt2.txt pptt.txt
    ./bin/mg5_aMC pptt.txt >> $LOG_FILE_PATH/log.txt
	python3 $ANALYSIS_SCRIPT_PATH/main.py -p 1 -u parse -m ttbar -i $ROOT_FILE_PATH/run_$i/tag_1_delphes_events.root -o /home/david/mass_generation/event_record_top_FHD_$i.npz -s 1 >> $LOG_FILE_PATH/log.txt
	rm $ROOT_FILE_PATH/run_$i/*
done
