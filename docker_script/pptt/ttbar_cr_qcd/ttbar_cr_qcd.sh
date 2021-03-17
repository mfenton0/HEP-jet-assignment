# Set up environment
source ~/.bash_profile;
cd ~/MG5_aMC_v2_7_2;
pip3 install h5py tqdm;

# Copy over all of the necessary files to run the script
/bin/cp -f /home/workplace/HEP-jet-assignment/madgraph/delphes_card_ATLAS.tcl ~/MG5_aMC_v2_7_2/Delphes/cards/;
/bin/cp -f /home/workplace/HEP-jet-assignment/madgraph/pptt_preparation/color_reconnection/cr_qcd/pythia8_card_qcd_inspire.dat ~/MG5_aMC_v2_7_2/;
/bin/cp -avr /home/workplace/HEP-jet-assignment/madgraph/pptt_preparation/color_reconnection/cr_qcd/* ./;

# Set Seed
SEED=$(cat /home/david/seed.txt)
awk -v seedval=$SEED -F"=" 'BEGIN{OFS=FS} $1=="set iseed "{$2=" "seedval}1' CR_qcd_inspire.txt > CR_qcd_inspire2.txt && mv -f CR_qcd_inspire2.txt CR_qcd_inspire.txt
awk -v seedval=$SEED -F"=" 'BEGIN{OFS=FS} $1=="set iseed "{$2=" "seedval}1' CR_qcd_inspire_first_round.txt > CR_qcd_inspire_first_round2.txt && mv -f CR_qcd_inspire_first_round2.txt CR_qcd_inspire_first_round.txt

# Create output directory and log our config
mkdir -p /home/david/mass_generation_cr_qcd;
/bin/cp pptt.txt /home/david/mass_generation_cr_qcd/

# Start the simulation
./run.sh
