#!/bin/bash
source ~/.bash_profile
for i in {1..9};do
	if [ $i == 1 ]
       	then
		echo $i 
        	nohup ./bin/mg5_aMC ppjjjj_ckkw_first_run.txt >> /home/david/log_ppjjjj_ckkw.txt 
	else
 		echo $i
		nohup ./bin/mg5_aMC ppjjjj_ckkw.txt >> /home/david/log_ppjjjj_ckkw.txt
	fi 		
done
for i in {10..11};do
	echo $i
       	nohup ./bin/mg5_aMC ppjjjj_ckkw.txt >> /home/david/log_ppjjjj_ckkw.txt
done
