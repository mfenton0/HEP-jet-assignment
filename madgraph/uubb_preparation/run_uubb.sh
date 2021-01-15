#!/bin/bash
source ~/.bash_profile
for i in {1..10};do
	if [ $i == 1 ]
       	then
		echo $i 
        	nohup ./bin/mg5_aMC uubb_first_run.txt >> /home/david/log_uubb.txt 
	else
 		echo $i
		nohup ./bin/mg5_aMC uubb.txt >> /home/david/log_uubb.txt
	fi 		
done
for i in {10..100};do
	echo $i
       	nohup ./bin/mg5_aMC bkg.txt >> /home/david/log_bkg.txt
done
