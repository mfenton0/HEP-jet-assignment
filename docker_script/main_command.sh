docker run \
    -it \
    --cpus="1.0" \
    -v /home/alex/Research/FeynmanData/vagrant1:/home/david \
    -v /home/alex/Research/FeynmanData/workplace:/home/workplace \
    alan200276/centos:SVJsimulation \
    /bin/bash -c /home/workplace/HEP-jet-assignment/docker_script/ttbar.sh
    