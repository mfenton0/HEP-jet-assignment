#!/bin/bash

#source ~/.bashrc
#source /ROOT/bin/thisroot.sh
#source /herwig/bin/activate

export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

#mkdir -p ###YOUR_AREA###/mass_generation

YOUR_AREA='' # define by yourself
ANALYSIS_SCRIPT_PATH='###YOUR_AREA###/HEP-jet-assignment/analysis_script'
ROOT_FILE_PATH='###YOUR_AREA###/HEP-jet-assignment/output/pptt_lep/Events'
LOG_FILE_PATH='###YOUR_AREA###/HEP-jet-assignment/logs'
MG_PATH='###YOUR_AREA###/MG5_aMC_v3_2_0'
echo "Path of script: $ANALYSIS_SCRIPT_PATH"
echo "path of root file: $ROOT_FILE_PATH"
echo "path of log file: $LOG_FILE_PATH"
echo "path of MG5 dir: $MG_PATH"
for i in {1..10};do
	if [ $i == 1 ]
       	then
		echo "Running round: $i, SEED: $SEED"
		echo "Running MG5."
	        $MG_PATH/bin/mg5_aMC pptt_lep_first_round.txt >> $LOG_FILE_PATH/mg5_pptt_lep.log
		echo "Parsing root file generate by pythia8 showering."
		python3 $ANALYSIS_SCRIPT_PATH/main.py -p 1 -u parse -m ttbar_lep -i $ROOT_FILE_PATH/run_0$i/tag_1_delphes_events.root -o $YOUR_AREA/mass_generation/event_record_top_semilep_$i -s 1 >> $LOG_FILE_PATH/log.txt
		#echo "Removing generate files."
		#rm $ROOT_FILE_PATH/run_0$i/*
	else
 		echo "Running round: $i, SEED: $SEED"
		awk -F"=" 'BEGIN{OFS=FS} $1=="set iseed "{$2=" "$2+1}1' pptt_lep.txt > pptt_lep2.txt && mv -f pptt_lep2.txt pptt_lep.txt
		echo "Running MG5."
		$MG_PATH/bin/mg5_aMC pptt_lep.txt >> $LOG_FILE_PATH/log.txt
		echo "Parsing root file generate by pythia8 showering."
		python3 $ANALYSIS_SCRIPT_PATH/main.py -p 1 -u parse -m ttbar_lep -i $ROOT_FILE_PATH/run_0$i/tag_1_delphes_events.root -o $YOUR_AREA/mass_generation/event_record_top_semilep_$i -s 1 >> $LOG_FILE_PATH/log.txt
		#echo "Removing generate files."
		#rm $ROOT_FILE_PATH/run_0$i/*
	fi 		
done

for i in {10..10};do
 	echo "Running round: $i, SEED: $SEED"
	awk -F"=" 'BEGIN{OFS=FS} $1=="set iseed "{$2=" "$2+1}1' pptt_lep.txt > pptt_lep2.txt && mv -f pptt_lep2.txt pptt_lep.txt
	echo "Running MG5."
	$MG_PATH/bin/mg5_aMC pptt_lep.txt >> $LOG_FILE_PATH/log.txt
	echo "Parsing root file generate by pythia8 showering."
	python3 $ANALYSIS_SCRIPT_PATH/main.py -p 1 -u parse -m ttbar_lep -i $ROOT_FILE_PATH/run_$i/tag_1_delphes_events.root -o $YOUR_AREA/mass_generation/event_record_top_semilep_$i -s 1 >> $LOG_FILE_PATH/log.txt
	#echo "Removing generate files."
	#rm $ROOT_FILE_PATH/run_$i/*
done
