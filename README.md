# ML Force Field Optimisation (BETA)

**FFOPTI** is a python package for optimising flexible force fields of structures automatically using machine learning methods paired with GROMACS or LAMMPS software. Currently, FFOPTI only support GROMACS.
**FFOPTI** also have built in functions that can be used to extract connectivity in atoms and assign UFF force fields for classical MD without optimisation (see: class About_Structure)

# To Start
- To extract all the informations required a cif, mdf and car files are required. Making sure that all of these files have the same name. Stop here if FF optimisation is not required.

- To continue with FF optimisation, DFT reference data needs to be provided. We currently only support CP2K simulations. Please make sure that the energy, forces and the coordinate files (md-1.ener, md-frc-1.xyz and md-pos-1.xyz or md-pos-1.pdb (preferred)) are all included in the same $DFT_path within the structure directory (i.e $Structure_name/DFT_path)

- It is preferred if a few intial parameters are provided within the directory $Structure_name/classical/your_param for GROMACs, please provide a force_field.itp file. If not the initial structures will be assigned UFF, OPLS and CVFF force fields are an initial 'guesstimate' parameters.

- UPDATE: It is compulsory to provide an inital primary force_field.itp (at least 1). This file will be used to create limits for the parameters to ensure that the parameters are wirhin a physically reasonable range.

# How it works
- A single point energy calculations will be done with the coordinates provided by DFT reference data to compute the loss of the forces.
- Bayesian Optimisation is used to optimise the parameters with the X_train data being the parameters being optimised and the Y_train is the computed loss of the forces.
- The optimiser will then be able to predict a new parameter, which will the be used to calculate the loss of the forces on each atom.
- The Bayesian Optimisation loop will stop once it reaches the set *max_iterations* or it has not seen any improvements in the last *loss_limit* iterations.

