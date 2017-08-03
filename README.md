# PyFoamSetup
A python library for setting up CFD simulations with OpenFOAM, with special focus on marine applications.

There is a general library that act as python interface to many of the settings one can adjust in OpenFOAM simulations. This library is used by several specialized classes that can set up simulations for specific cases. Currently the specialized classes consists of the following:
- TowingTank - A class that sets up simulations intended for finding ship resistance at different speeds, drift and heel angles. 
- HydrofoilTowingTank - A class that inherits from the TowingTank class, but that is altered to set up a simulation of a hydrofoil. The main difference is the length variables that is used to set up dimensions
- FoilSimulation - A class that sets up simulations for two-dimensional foils. Used to simulate both steady cases such as different angles of attack, but also unsteady cases, such as oscillating foils. 
- WingSimulation - A class that sets up simulations for three-dimensional wings
- PropelleSimulations - A class that sets up simulations of propellers, using sliding interface between and a rotating mesh with the propeller geometry.

## Dependencies
- Numpy
- Scipy
- My Mesh library [https://github.com/jarlekramer/polymesh]
- OpenFOAM. In theory it should work with both normal OpenFOAM and OpenFOAM+, but OpenFOAM+ is the version I am using, and therefore the only verison that is propely tested. The version of OpenFOAM+ that I currently use is 1612+. Compatebility with older versions will not be prioritised.  

## Install instructions
Download, cd into folder and execute python setup.py install

In order to install extensions to OpenFoam that comes with this library execute the "installOpenFoamExtenstensions.py" script in the "scripts" folder. This requires that OpenFoam is allready installed, and that the environment is set up (the FOAM_RUN environmental variable must be defined and the run folder must be created)
