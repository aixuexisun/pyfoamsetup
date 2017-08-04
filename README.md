# PyFoamSetup
A python library for setting up CFD simulations with OpenFOAM, with special focus on marine applications. That is, it is written by-and-for marine engineers, using CFD as an engineering tool. This is means that the supported features are somewhat limited to what is needed in the daily work of a marine enginner. The library is written to be as general as possible, so that additional features could be added with relative ease at a later stage if needed. Things you can currently do with this library is as follows:

- Simulate ship resistance. Both model scale and full scale is possible, but the library is mostly tuned to be efficient in model scale (roughly 3-10 m models). The ship can have a rudder, appendage, model a propeller using an actuator disk, and be either free to move in heave and pitch, or fixed. The simulation can either be symmetric about the ships center line, or not. The later is necessary when reistance as a function of drift and heel angle is of interest. Free surface modeling can be turned on or off.
- Simulate three-dimensional hydrofoils. Both with and without free surface modelling.
- Simulate air resistance on the superstructure of a ship
- Simulate forces on 2D foils. Both steady state simulations (for small angles of attack), unsteady simulations of fixed foils (large angles of attack), several foils together, multi element foils, and foils with arbitrary movement (such as oscillating foils for propulsion).
- Simulate forces on 3D wings. Typical use could be finding forces on a rudder, or a wing sail. Actuatur disks could be used to model the influence of propellers
- Simulate forces on rotating cylinders, such as a Flettner rotor
- Simulate a propeller using sliding interface and rotating mesh

Below comes a short explanation of how the library works, assuming basic knowledge of how OpenFoam works. There are examples scripts of how to use the library located in the "Examples" folder. More detailed explanation will come in the Wiki (work in progress).

The library is only intended for external flow simulations, allthough it should be possible to extend it other things as well. The OpenFOAM solvers that can be used with this library is simpleFoam (steady state incompressible), pimpleFoam (unsteady incompressible), and interFoam (unsteady, two-phase incompressible). The library supports several turbulence models; k-omega, k-omega SST, k-epsilon, realizable k-epsilon and Spalart-Allmaras. There is some limited support for LES/DES turubulence models as well, but as I mostly do RANS simulations, these features are not yet tested properly. 

There is a general library, located in the "pyfoamsetup/coreLibrary" folder, that act as python interface to many of the settings one can adjust in an OpenFOAM simulation. This is done by having python classes for most of the settings-files that is needed. For instance, one can set up a snappyHexMesh-settings-file by using the "Dict" class in "SnappyHexMesh.py". In the same way, one can set up a "fvSchemes" file using the "Dict" class in the "fvScemes.py". A "Dict" class is simply a collection of dictionaries containg keywords needed by OpenFOAM, together with values for each keyword, and a method for writing the settings to a dictionary file that OpenFoam understands.

In the core library, there is also a file called "CaseSetup.py", which contains a general class for setting up an OpenFOAM simulation. This class needs a representaive length dimension and velocity, a reference area, kinematic viscosity, density and an application name. The length dimensions and velocity is used to calculate the Reynolds number, and thereby decide the necessary cell length close to the walls in the simulation. The reference area is only used for setting up force coefficients. The input viscosity and density is the viscosity and density that used in the actual simulation. The application name is a used to set up some folder paths that is different depending on which "application" that is currently using the CaseSetup class. An application is a class that inherits from the CaseSetup class, but that specialize on something more specific. The current applications are as follows:
- TowingTank          - A class that sets up simulations of ships
- HydrofoilTowingTank - A class that inherits from the TowingTank class, but sets up simulations of hydrofoils
- FoilSimulation      - A class that sets up simulations of two-dimensional foils
- WingSimulation      - A class that sets up simulations of three-dimensional wings
- PropelleSimulations - A class that sets up simulations of propellers

Each application, together with the general CaseSetup class, does several things:
- They decide how the mesh should be created, based on the type of simulation and dimensions (length and velocity) of the problem
- They decide which solver, and which solver setings to be used
- They create the entire case folder for the simulation, and puts it in the ```$FOAM_RUN``` folder
- They set up all the necessary dictionary files, inluding the boundary conditions
- They set up inlet values, including corret values for the turbulence model.

Every thing is set up based on "best practices" determined by me, the developer of this software. As part of my research, I spend a lot of time trying different settings in OpenFOAM simulations, for the things that I work with (research into ships, ship resistance, hydrofoils and modern sail technology). I try to formalize my experience with CFD and OpenFOAM into a set of rules that is used by the applications in this library. My experience is off-course limited; there is alwasy a chance that I have done something wrong, or that the settings I think is appropriate should be something else. I will publish justification for at least some of the settings in the Wiki at a later stage (work in progress). I have most experience with ship resistance and 2D foils simulations. I do very little propeller simulations, so the propeller application is not very well developed at the moment (I hope to improve it in the future). My goal is to set up a good balance between accuracy and speed. What that actually means is a huge topic of discussion, but the goal is to set up simulations that are practical to run on a modern desktop work station, and do not need a super computer. For instance, finding the ship reistance of a 7 m ship model, with a default settings from this library on a 16-core work station takes around 12-hours.  

## Dependencies
- Numpy
- Scipy
- My Mesh library [https://github.com/jarlekramer/polymesh]
- OpenFOAM. In theory it should work with both normal OpenFOAM and OpenFOAM+, but OpenFOAM+ is the version I am using, and therefore the only verison that is propely tested. The version of OpenFOAM+ that I currently use is 1612+. Compatebility with older versions will not be prioritised.  

## Install instructions
This package can be installed in the same way as most other python packages by executing the following in the source folder:

```
python setup.py install
```

It only works for python 3, so you might need to use ```python3``` rather than ```python```, depending on your setup

In order to install extensions to OpenFoam that comes with this library execute the "installOpenFoamExtenstensions.py" script in the "scripts" folder. This requires that OpenFoam is allready installed, and that the environment is set up (the FOAM_RUN environmental variable must be defined and the run folder must be created)
