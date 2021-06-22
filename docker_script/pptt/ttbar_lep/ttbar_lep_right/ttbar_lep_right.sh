# Set up environment
source ~/.bash_profile;
cd ~/MG5_aMC_v2_7_2;
pip3 install h5py tqdm;

# Copy over all of the necessary files to run the script
/bin/cp -f /workplace/HEP-jet-assignment/madgraph/delphes_card_ATLAS.tcl ~/MG5_aMC_v2_7_2/Delphes/cards/;
/bin/cp -avr /workplace/HEP-jet-assignment/madgraph/ttbar-semi-lep_preparation/right/* ./;

# Set Seed
SEED=$(cat /home/david/seed.txt)
awk -v seedval=$SEED -F"=" 'BEGIN{OFS=FS} $1=="set iseed "{$2=" "seedval}1' pptt_lep.txt > pptt_lep2.txt && mv -f pptt_lep2.txt pptt_lep.txt
awk -v seedval=$SEED -F"=" 'BEGIN{OFS=FS} $1=="set iseed "{$2=" "seedval}1' pptt_lep_first_round.txt > pptt_lep_first_round2.txt && mv -f pptt_lep_first_round2.txt pptt_lep_first_round.txt

# Create output directory and log our config
mkdir -p /home/david/mass_generation;
/bin/cp pptt_lep.txt /home/david/mass_generation/

# Start the simulation
./run.sh
