WASA Fast Simulation
Andre Nepomuceno - january 2025

Observations
1. You need Geant4 version 4.11 and ROOT 6.x to run this application
2. At the moment, only the electromagnetic calorimeter response is implemented
3. The source files are located in WASA_FAST/src directory
4. The head files are located in WASA_FAST/include directory

HOW TO COMPILE

1. Copy the entire directory WASA_FAST_V2 to some location. For example, let us say that this location is /home/you/.
2. In the /home/you/ directory, alongside  WASA_FAST folder, create the folder WASA-build:

$ mkdir WASA-build
$ ls
WASA-build  WASA_FAST_V2

3. Inside the folder WASA-build, run CMake:

$cmake -DCMAKE_PREFIX_PATH=<path_to_geant4-v11-install>  /home/you/WASA_FAST_V2

where <path_to_geant4-v11-install> is the path to where Geant4.11 is installed. 

4. Compile:

$ make

If all goes well, the executable file wasa_main will be created in your WASA-build directory.

HOW TO  RUN

1. Visualize the geometry
   
To visualize the geometry, we must have Qt5 and OpenGL libraries and header in your system, and Geant4 must have been compile with the option GEANT4_USE_QT=ON.
To  visualize the geometry, first run the executable:

$ ./wasa_main

On the visualization window, type after “Session”: /control/execute vis_wasa.mac

2. Run a simulation
   
The macro wasa_simulation.in is the input file to simulate a neutral pion (pi0) event decaying into photons that will hit the EM calorimeter. This simulation shoots 10000 neutral pions of kinetic energy of 200 MeV,  from the center of the detector, in random directions. Check the file to see how to change the primary particle, kinematics and number of events.
It is possible to generate different energies and angle distributions. For details, see the Geant4 General Particle Source (GPS) documentation: 
https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/GettingStarted/generalParticleSource.html#g4gps

To run the simulation, type on a Linux terminal:

$./wasa_main wasa_simulation.in

The output of the simulation (energy recorded in the EM calorimeter and hit position) will be store in four ROOT files WASAFastOutput_t0.root. 

Run the macro pi0_analysis.C to analyze the output:

$ root -l 

root [0] .x pi0_analysis_v2.C
