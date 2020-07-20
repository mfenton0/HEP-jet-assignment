import uproot

root_file = './tag_1_delphes_events.root'

data = uproot.open(root_file)['Delphes']
data.keys()
data.show()