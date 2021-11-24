# Summation of currents electromagnetic solver

This solver can be used to predict the self-inductance and mutual inductance of coils, as well as the inductance change due to the presence of a conductive object near the coil. To do so the solver can create object that only have mesh elements on the surface in order to be able to model the skin effect.

The code has been developed to predict the influence of a 3D nozzle on the inductance of a coil in order to be able to calibration the nozzle offsets in x and y. The method has been described in the follow [paper]() (open access version [available]())

A GUI for calibrating a printer in combination with a LDC1101EVM can be found [here](https://github.com/martijnschouten/inductive_calibration_GUI).

# Typical usage
1. Copy SOC_object.m into your matlab work folder
1. Make an coil using object using the SOC_object class `coil = SOC_object`
1. Make an nozzle using object using the SOC_object class `nozzle = SOC_object`
1. Set the parameters of the geometry and mesh of both the coil and the nozzle
1. Build the mesh
1. Permform calculations

Run `doc SOC_object` in your matlab terminal to get more info on the different functions to set the parameters, mesh and do calculations. Also make sure to check out the example scripts.

# Requirements
This program by default uses the gpu support in Matlab's parallel processing toolbox to perform some big matrix invertions. This means that in order to run at a reasonable speed the script you will need a NVIDIA CUDA enabled GPU with a couple of gigs of memory. If you don't have such a GPU the matrix inversions can also be run on the cpu by using
```
coil.use_gpu = false
```

# Examples


# Acknowledgement
This work was developed within the Wearable Robotics programme, funded by the Dutch Research Council (NWO)

<img src="https://user-images.githubusercontent.com/6079002/124443163-bd35c400-dd7d-11eb-9fe5-53c3def86459.jpg" width="62" height="100"><img src="https://user-images.githubusercontent.com/6079002/124443273-d3dc1b00-dd7d-11eb-9282-54c56e0f42db.png" width="165" height="100">
