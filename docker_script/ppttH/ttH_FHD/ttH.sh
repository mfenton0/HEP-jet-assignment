# Set up environment
source ~/.bash_profile;
cd ~/MG5_aMC_v2_7_2;
pip3 install h5py tqdm;

# Copy over all of the necessary files to run the script
/bin/cp -f /workplace/HEP-jet-assignment/madgraph/delphes_card_ATLAS.tcl ~/MG5_aMC_v2_7_2/Delphes/cards/;
/bin/cp -avr /workplace/HEP-jet-assignment/madgraph/ppttH_preparation/* ./;

# Set Seed
SEED=$(cat /home/david/seed.txt)
awk -v seedval=$SEED -F"=" 'BEGIN{OFS=FS} $1=="set iseed "{$2=" "seedval}1' ttH.txt > ttH2.txt && mv -f ttH2.txt ttH.txt
awk -v seedval=$SEED -F"=" 'BEGIN{OFS=FS} $1=="set iseed "{$2=" "seedval}1' ttH_first_round.txt > ttH_first_round2.txt && mv -f ttH_first_round2.txt ttH_first_round.txt

# Create output directory and log our config
mkdir -p /home/david/mass_generation;
/bin/cp ttH.txt /home/david/mass_generation/

# Start the simulation
./run_ttH.sh
