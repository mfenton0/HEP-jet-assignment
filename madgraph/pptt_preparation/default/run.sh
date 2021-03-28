#!/bin/bash
source ~/.bashrc
source /ROOT/bin/thisroot.sh
source /herwig/bin/activate
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

mkdir -p /mass_generation

ANALYSIS_SCRIPT_PATH='/workplace/HEP-jet-assignment/analysis_script'
ROOT_FILE_PATH='/MAP_VOLUME_FOR_CONTAINER/pptt/Events'
LOG_FILE_PATH='/MAP_VOLUME_FOR_CONTAINER'
echo "Path of script: $ANALYSIS_SCRIPT_PATH"
echo "path of root file: $ROOT_FILE_PATH"
echo "path of log file: $LOG_FILE_PATH"
for i in {1..9};do
	if [ $i == 1 ]
       	then
		echo $i
	        ./bin/mg5_aMC pptt_first_round.txt >> $LOG_FILE_PATH/log.txt
		gunzip $ROOT_FILE_PATH/run_0$i/unweighted_events.lhe.gz >> $LOG_FILE_PATH/log.txt
		sed -E -i "s/run_[0-9]{2,}/run_0${i}/g" pptt_mg_herwig.in 
		#grep 'run_' pptt_mg_herwig.in
		Herwig read pptt_mg_herwig.in >> $LOG_FILE_PATH/log.txt
		Herwig run pptt_mg_herwig.run -j 1 -N 10000 -s $SEED -d 1 >> $LOG_FILE_PATH/log.txt
		/delphes/DelphesHepMC /MG5_aMC_v2_7_3/Delphes/cards/delphes_card_ATLAS.tcl $ROOT_FILE_PATH/run_0$i/herwig_run_0$i.root  $ROOT_FILE_PATH/run_0$i/pptt_angular_run_0$i.hepmc
		python3 $ANALYSIS_SCRIPT_PATH/main.py -p 1 -u parse -m ttbar -i $ROOT_FILE_PATH/run_0$i/tag_1_delphes_events.root -o /mass_generation/event_record_top_FHD_$i.npz -s 1 >> $LOG_FILE_PATH/log.txt
		python3 $ANALYSIS_SCRIPT_PATH/main.py -p 1 -u parse -g herwig7 -m ttbar -i $ROOT_FILE_PATH/run_0$i/herwig_run_0$i.root -o /mass_generation/event_record_top_FHD_herwig_$i.npz -s 1 >> $LOG_FILE_PATH/log.txt
		#rm $ROOT_FILE_PATH/run_0$i/*
	else
 		echo $i
		#echo "seed = " $(( $SEED + $i ))
		awk -F"=" 'BEGIN{OFS=FS} $1=="set iseed "{$2=" "$2+1}1' pptt.txt > pptt2.txt && mv -f pptt2.txt pptt.txt
		./bin/mg5_aMC pptt.txt >> $LOG_FILE_PATH/log.txt
		gunzip  $ROOT_FILE_PATH/run_0$i/unweighted_events.lhe.gz >> $LOG_FILE_PATH/log.txt
		sed -E -i "s/run_[0-9]{2,}/run_0${i}/g" pptt_mg_herwig.in
		#grep 'run_' pptt_mg_herwig.in
		Herwig read pptt_mg_herwig.in >> $LOG_FILE_PATH/log.txt
		Herwig run pptt_mg_herwig.run -j 1 -N 10000 -s $(( $SEED + $i )) -d 1 >> $LOG_FILE_PATH/log.txt
		/delphes/DelphesHepMC /MG5_aMC_v2_7_3/Delphes/cards/delphes_card_ATLAS.tcl $ROOT_FILE_PATH/run_0$i/herwig_run_0$i.root  $ROOT_FILE_PATH/run_0$i/pptt_angular_run_0$i.hepmc
		python3 $ANALYSIS_SCRIPT_PATH/main.py -p 1 -u parse -m ttbar -i $ROOT_FILE_PATH/run_0$i/tag_1_delphes_events.root -o /mass_generation/event_record_top_FHD_$i.npz -s 1 >> $LOG_FILE_PATH/log.txt
		python3 $ANALYSIS_SCRIPT_PATH/main.py -p 1 -u parse -g herwig7 -m ttbar -i $ROOT_FILE_PATH/run_0$i/herwig_run_0$i.root -o /mass_generation/event_record_top_FHD_herwig_$i.npz -s 1 >> $LOG_FILE_PATH/log.txt
		#rm $ROOT_FILE_PATH/run_0$i/*
	fi 		
done
for i in {10..100};do
	echo $i
	awk -F"=" 'BEGIN{OFS=FS} $1=="set iseed "{$2=" "$2+1}1' pptt.txt > pptt2.txt && mv -f pptt2.txt pptt.txt
	./bin/mg5_aMC pptt.txt >> $LOG_FILE_PATH/log.txt
	gunzip  $ROOT_FILE_PATH/run_$i/unweighted_events.lhe.gz >> $LOG_FILE_PATH/log.txt
        sed -E -i "s/run_[0-9]{2,}/run_${i}/g" pptt_mg_herwig.in
	#grep 'run_' pptt_mg_herwig.in
        Herwig read pptt_mg_herwig.in >> $LOG_FILE_PATH/log.txt
        Herwig run pptt_mg_herwig.run -N 10000 -s $(( $SEED + $i )) -d 1 >> $LOG_FILE_PATH/log.txt
	/delphes/DelphesHepMC /MG5_aMC_v2_7_3/Delphes/cards/delphes_card_ATLAS.tcl run_0$i.root  $ROOT_FILE_PATH/run_$i/pptt_angular_run_0$i.hepmc
	python3 $ANALYSIS_SCRIPT_PATH/main.py -p 1 -u parse -m ttbar -i $ROOT_FILE_PATH/run_$i/tag_1_delphes_events.root -o /mass_generation/event_record_top_FHD_$i.npz -s 1 >> $LOG_FILE_PATH/log.txt
	python3 $ANALYSIS_SCRIPT_PATH/main.py -p 1 -u parse -g herwig7 -m ttbar -i $ROOT_FILE_PATH/run_$i/herwig_run_0$i.root -o /mass_generation/event_record_top_FHD_herwig_$i.npz -s 1 >> $LOG_FILE_PATH/log.txt
	#rm $ROOT_FILE_PATH/run_$i/*
done

