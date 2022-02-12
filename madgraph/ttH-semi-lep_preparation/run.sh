#!/bin/bash

#source ~/.bashrc
#source /ROOT/bin/thisroot.sh
#source /herwig/bin/activate

#conda activate hepjet
#cd /afs/ihep.ac.cn/users/o/okawa/data/mana/HEP-jet-assignment-4/madgraph/pptt_preparation/default 

export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

#mkdir -p /home/david/mass_generation

ANALYSIS_SCRIPT_PATH='###YOUR_AREA###/HEP-jet-assignment/analysis_script'
ROOT_FILE_PATH='###YOUR_AREA###/HEP-jet-assignment/output/ppttH_lep/Events'
LOG_FILE_PATH='###YOUR_AREA###/HEP-jet-assignment/logs'
MG_PATH='###YOUR_AREA###/MG5_aMC_v3_2_0'
OUTPUT_PATH='###YOUR_OUTPUTAREA###'
echo "Path of script: $ANALYSIS_SCRIPT_PATH"
echo "path of root file: $ROOT_FILE_PATH"
echo "path of log file: $LOG_FILE_PATH"
for i in {1..10};do
	if [ $i == 1 ]
       	then
		echo "Running round: $i, SEED: $SEED"
		echo "Running MG5."
	        $MG_PATH/bin/mg5_aMC ttH_first_round.txt > $LOG_FILE_PATH/mg5_ppttH_lep.log
		echo "Parsing root file generate by pythia8 showering."
		python3 $ANALYSIS_SCRIPT_PATH/main.py -p 1 -u parse -m ttH_lep -i $ROOT_FILE_PATH/run_0$i/tag_1_delphes_events.root -o $OUTPUT_PATH/mass_generation/event_record_ttH_semilep_$i -s 1 >> $LOG_FILE_PATH/log.txt
		echo "Removing generate files."
		rm $ROOT_FILE_PATH/run_0$i/*
	else
 		echo "Running round: $i, SEED: $SEED"
		awk -F"=" 'BEGIN{OFS=FS} $1=="set iseed "{$2=" "$2+1}1' ttH.txt > ttH2.txt && mv -f ttH2.txt ttH.txt
		echo "Running MG5."
		$MG_PATH/bin/mg5_aMC ttH.txt >> $LOG_FILE_PATH/mg5_ttH_lep.log
		echo "Parsing root file generate by pythia8 showering."
		python3 $ANALYSIS_SCRIPT_PATH/main.py -p 1 -u parse -m ttH_lep -i $ROOT_FILE_PATH/run_0$i/tag_1_delphes_events.root -o $OUTPUT_PATH/mass_generation/event_record_ttH_semilep_$i -s 1 >> $LOG_FILE_PATH/log.txt
		echo "Removing generate files."
	        rm $ROOT_FILE_PATH/run_0$i/*
	fi 		
done

