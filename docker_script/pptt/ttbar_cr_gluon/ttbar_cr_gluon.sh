# Set up environment
source ~/.bash_profile;
cd ~/MG5_aMC_v2_7_2;
pip3 install h5py tqdm;

# Copy over all of the necessary files to run the script
/bin/cp -f /home/workplace/HEP-jet-assignment/madgraph/delphes_card_ATLAS.tcl ~/MG5_aMC_v2_7_2/Delphes/cards/;
/bin/cp -f /home/workplace/HEP-jet-assignment/madgraph/pptt_preparation/color_reconnection/cr_gluon/pythia8_card_gluon_move.dat ~/MG5_aMC_v2_7_2/;
/bin/cp -avr /home/workplace/HEP-jet-assignment/madgraph/pptt_preparation/color_reconnection/cr_gluon/* ./;

# Set Seed
SEED=$(cat /home/david/seed.txt)
awk -v seedval=$SEED -F"=" 'BEGIN{OFS=FS} $1=="set iseed "{$2=" "seedval}1' CR_gluon_move.txt > CR_gluon_move2.txt && mv -f CR_gluon_move2.txt CR_gluon_move.txt
awk -v seedval=$SEED -F"=" 'BEGIN{OFS=FS} $1=="set iseed "{$2=" "seedval}1' CR_gluon_move_first_round.txt > CR_gluon_move_first_round2.txt && mv -f CR_gluon_move_first_round2.txt CR_gluon_move_first_round.txt

# Create output directory and log our config
mkdir -p /home/david/mass_generation_cr_gluon;
/bin/cp pptt.txt /home/david/mass_generation_cr_gluon/

# Start the simulation
./run.sh
