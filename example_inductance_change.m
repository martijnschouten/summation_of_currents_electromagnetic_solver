RhoC = 1.68e-8;%ohm meter
NL = 55;%number of windings of coil 
HL = 1.4e-3;%m height of coil
DL1 = 1.84e-3;%m inner diameter of coil
DL2 = 2e-3;%m outer diameter of coil

nxL = 3;%number of elements on the coil in the x direction
nzL = 15;%number of elements on the coil in the z direction
nphiL = 30;%number of elements of the coil in radial direction

RhoB = 1/0.28*1.68e-8;%ohm meter
NN = 1;%number of windings of nozzle
HN = 1e-3;%m height of nozzle
DN1 = 0.5e-3;%m inner diameter of nozzle
DN2 = 3e-3;%m outer diameter of nozzle

zpos = 2e-3;%Location of the object in z in meter

nxN = 14;%Number of elements on the bottom/top boundary mesh in the x direction
nzN = 11;%Number of elements on the left/right/center boundary mesh in the z direction
nphiN = 20;%Number of elements of the coil in radial direction
nxbl = 3;%Number of elements on the left side in the x direction
nxbr = 10;%Number of elements on the right side in the x direction
nxc = 7;%Number of elements on the center in the x direction
nzbb = 10;%Number of elements on the bottom side in the z direction
nzbt = 0;%Number of elements on the top side in the z direction

xpos = 0;%Location of the object in x in meter
f = 2e6;%Frequency at which the calculation will take place
nmesh = 2;%Number of times the skind depth that the boundary mesh should extend

taper_width = 0.75e-3;%Horizontal distance of the taper of the nozzle in meter
taper_angle = pi/4;%Angle of the taper of the nozzle in radians
hole_offset = 0e-3;%How much the hole of the nozzle is offset from the center in meter
rotation = 0;%Rotation of the nozzle in degrees

sense_coil = SOC_object;%create an object for the sense coil
sense_coil = sense_coil.set_geometry(HL,DL1,DL2,NL,RhoC,0,0);%set the geometry of the coil
sense_coil = sense_coil.set_mesh(nxL,nzL,nphiL);%set the mesh of the coil
sense_coil = sense_coil.build_coil();%build the coil

%create a nozzle geometery in the 
nozzle = SOC_object;%create an object for the nozzle
nozzle = nozzle.set_nozzle_geometry(HN,DN1,DN2,RhoB,xpos,zpos,taper_width,taper_angle,hole_offset,rotation); %define the nozzle geometry
nozzle = nozzle.set_boundary_mesh(nxN,nzN,nphiN,nxbl,nxbr,nxc,nzbb,nzbt);%configure a mesh that only has points at the boundary of the coil
skin_depth = nozzle.skin_depth(f);%calculate the skin depth, so that we know approximately how deep the mesh should extent into the nozzle.
nozzle = nozzle.build_nozzle_surface(nmesh*skin_depth);%build the boundary mesh

%calculate the change in impedance due to the presence of the nozzle
dZ = sense_coil.calculate_impedance_change(nozzle,f);
disp('Change in impdance: ')
disp(dZ)