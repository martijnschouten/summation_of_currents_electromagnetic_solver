# Summation of currents based FEM solver

This solver can be used to predict the self-inductance and mutual inductance of coils, as well as the inductance change due to the presence of a conductive object near the coil. To do so the solver can create object that only have mesh elements on the surface in order to be able to model the skin effect.

The code has been developed to predict the influence of a 3D nozzle on the inductance of a coil in order to be able to calibration the nozzle offsets in x and y. The paper describing this method can be found in:
doi to actual published paper
open source author version

The code for calibrating a printer can be found in:

# Typical usage
The induced current density example demonstrates how to use the solver to calculate the current density induced in a 3D printer nozzle by the coil

The change impedance example demonstrates how to calculate the change impedance due to the presence of a 3D printer nozzle.

# Requirements
This program uses the gpu support in Matlab's parallel processing toolbox to perform some big matrix invertions. This means that in order to run the script you will need a NVIDIA CUDA enabled GPU with a couple of gigs of memory.

# Acknowledgement
This work was developed within the Wearable Robotics programme, funded by the Dutch Research Council (NWO)

<img src="https://user-images.githubusercontent.com/6079002/124443163-bd35c400-dd7d-11eb-9fe5-53c3def86459.jpg" width="62" height="100"><img src="https://user-images.githubusercontent.com/6079002/124443273-d3dc1b00-dd7d-11eb-9282-54c56e0f42db.png" width="165" height="100">
