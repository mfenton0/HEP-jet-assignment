if [ $# -lt 3 ]; then
	echo "Usage: create.sh INITIAL_SEED NUM_RUNS NUM_PARALLEL"
	exit 1
fi

INITIALSEED="$1"
NUM_RUNS="$2"
PARALLEL="$3"

HOME=`pwd`

for (( i=0; i<$PARALLEL; i++ ))
do
	j=`expr $NUM_RUNS \* $i + $INITIALSEED + 1`
	echo $j

	mkdir docker_$j
	echo $j > docker_$j/seed.txt

	docker run \
		-d \
		-v $HOME/workplace:/workplace \
		-v $HOME/docker_$j:/home/david \
		alan200276/centos:SVJsimulation \
		/bin/bash -c /workplace/HEP-jet-assignment/docker_script/pptttt/four_top_FHD/four_top.sh
done
