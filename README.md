# Quantum-Espresso-Tools
Collection of Quantum Espresso tools to analyze phonon calculations. It consists of 4 files which are used to organize and produce a projection heat map of phonon calculations onto the Gamma eigenvectors orthogonal basis. The first code called eigenvalprocess.sh split the Quantum Espresso output file of eigenvectors into individual q-points subfiles automatically. 
The second file called projection.py performs a projection of any two q-points from the previous code. The third code called 
gen_modes.py generates an arrow representation of the modes computed by Quantum Espresso.


