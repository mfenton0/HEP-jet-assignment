#!/bin/bash

#source ~/.bashrc
#source /ROOT/bin/thisroot.sh
#source /herwig/bin/activate

export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

#mkdir -p /home/david/mass_generation

ANALYSIS_SCRIPT_PATH='###YOUR AREA###/HEP-jet-assignment/analysis_script'
ROOT_FILE_PATH='###YOUR AREA###/HEP-jet-assignment/output/eett_had/Events'
LOG_FILE_PATH='###YOUR AREA###/HEP-jet-assignment/logs'
MG_PATH='###YOUR MG5 AREA###'
echo "Path of script: $ANALYSIS_SCRIPT_PATH"
echo "path of root file: $ROOT_FILE_PATH"
echo "path of log file: $LOG_FILE_PATH"
for i in {1..9};do
	if [ $i == 1 ]
       	then
		echo "Running round: $i, SEED: $SEED"
		echo "Running MG5."
	        $MG_PATH/bin/mg5_aMC eett_first_round.txt >> $LOG_FILE_PATH/mg5_eett_had.log
		echo "Parsing root file generate by pythia8 showering."
		python3 $ANALYSIS_SCRIPT_PATH/main.py -p 1 -u parse -m ttbar -i $ROOT_FILE_PATH/run_0$i/tag_1_delphes_events.root -o /afs/cern.ch/user/h/hokawa/workspace/private/mana/mass_generation/event_record_eett_had_$i -s 1 >> $LOG_FILE_PATH/log.txt
		echo "Removing generate files."
		rm $ROOT_FILE_PATH/run_0$i/*
	else
 		echo "Running round: $i, SEED: $SEED"
		awk -F"=" 'BEGIN{OFS=FS} $1=="set iseed "{$2=" "$2+1}1' eett.txt > eett2.txt && mv -f eett2.txt eett.txt
		echo "Running MG5."
		$MG_PATH/bin/mg5_aMC eett.txt >> $LOG_FILE_PATH/log.txt
		echo "Parsing root file generate by pythia8 showering."
		python3 $ANALYSIS_SCRIPT_PATH/main.py -p 1 -u parse -m ttbar -i $ROOT_FILE_PATH/run_0$i/tag_1_delphes_events.root -o /afs/cern.ch/user/h/hokawa/workspace/private/mana/mass_generation/event_record_eett_had_$i -s 1 >> $LOG_FILE_PATH/log.txt
		echo "Removing generate files."
		rm $ROOT_FILE_PATH/run_0$i/*
	fi 		
done

for i in {10..10};do
 	echo "Running round: $i, SEED: $SEED"
	awk -F"=" 'BEGIN{OFS=FS} $1=="set iseed "{$2=" "$2+1}1' eett.txt > eett2.txt && mv -f eett2.txt eett.txt
	echo "Running MG5."
	$MG_PATH/bin/mg5_aMC eett.txt >> $LOG_FILE_PATH/log.txt
	echo "Parsing root file generate by pythia8 showering."
	python3 $ANALYSIS_SCRIPT_PATH/main.py -p 1 -u parse -m ttbar -i $ROOT_FILE_PATH/run_$i/tag_1_delphes_events.root -o /afs/cern.ch/user/h/hokawa/workspace/private/mana/mass_generation/event_record_eett_had_$i -s 1 >> $LOG_FILE_PATH/log.txt
	echo "Removing generate files."
	rm $ROOT_FILE_PATH/run_$i/*
done
