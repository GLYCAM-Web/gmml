//Instructions for making the ontology

In home directory,
git clone -b gems-dev https://github.com/GLYCAM-Web/gems.git
cd gems
git clone -b gmml-dev https://github.com/GLYCAM-Web/gmml.git
./make.sh
cd ~/OntologyScripts
nohup ./controlCPU.sh /programs/repos/PDB/all 1 &

//1 CPU used because the detect_sugars_runArchive.cc file uses 20 CPUS to build the structure.
