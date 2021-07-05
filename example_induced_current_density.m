clear all
close all
close all force

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
HN = 2e-3;%m height of nozzle
DN1 = 0.5e-3;%m inner diameter of nozzle
DN2 = 3e-3;%m outer diameter of nozzle

zpos = 2.5e-3; %distance between the bottom of the coil and the center bottom of the nozzle

nxN = 10;%number of elements on the coil in the x direction
nzN = 10;%number of elements on the coil in the z direction
nphiN = 15;%number of elements of the coil in radial direction
nxb = 10;%number of elements of the vertical boundaries
nzb = 15;%number of elements of the horizonal boundaries

xpos = 0;%position of the nozzle relative to the coil
f = 2e6;%frequency at which the measurement is performed
nmesh = 3;%number of times the skin depth that the boundary mesh will extend into the nozzle


taper_width = 0.75e-3; %distance between the outer radius of the nozzle and the point where the taper starts
taper_angle = pi/4; %the angle of the taper of the nozzle
hole_offset = 0e-6; %offset of the hole inside the nozzle
rotation = 2; %rotation of the nozzle in degree

sense_coil = SOC_object;%create an object for the sense coil
sense_coil = sense_coil.set_coil_geometry(HL,DL1,DL2,NL,RhoC,0,0);%set the geometry of the coil
sense_coil = sense_coil.set_mesh(nxL,nzL,nphiL);%set the mesh of the coil
sense_coil = sense_coil.build_coil();%build the coil

%calculate the self inductance and series resistance of the coil.
L_self = sense_coil.calculate_self_inductance();
disp("Inductance: ")
disp(L_self);
R = sense_coil.resistance();
disp("Resistance: ")
disp(R);


%create a nozzle geometery in the 
nozzle = SOC_object;
nozzle = nozzle.set_nozzle_geometry(HN,DN1,DN2,NN,RhoB,xpos,zpos,taper_width,taper_angle,hole_offset,rotation);
nozzle = nozzle.set_boundary_mesh(nxN,nzN,nphiN,nxb,nzb);
skin_depth = nozzle.skin_depth(f);
disp("Skin depth: ")
disp(skin_depth)
nozzle = nozzle.build_nozzle_surface(nmesh*skin_depth);

%plot the geometry of the nozzle
figure
sense_coil.plot_geometry(3,[1e3,1e3,1e3]);
xlabel('x position (mm)')
ylabel('y position (mm)')
zlabel('z position (mm)')

%plot the geometry of the nozzle
figure
nozzle.plot_geometry(3,[1e3,1e3,1e3]);
xlabel('x position (mm)')
ylabel('y position (mm)')
zlabel('z position (mm)')

%plot the geometry of the nozzle and the coil in one plot
figure
sense_coil.plot_geometry(3,[1e3,1e3,1e3]);
xlabel('x position (mm)')
ylabel('y position (mm)')
zlabel('z position (mm)')
hold on
nozzle.plot_geometry(3,[1e3,1e3,1e3]);
xlabel('x position (mm)')
ylabel('y position (mm)')
zlabel('z position (mm)')

%plot the surface vectors used for 
figure
sense_coil.plot_surface_vectors([1e3,1e3,1e3])
 xlabel('x position (mm)')
ylabel('y position (mm)')
zlabel('z position (mm)')
view(2)

%plot the current density induced in the nozzle
figure
sense_coil.plot_induced_current_density_2d(nozzle,1,f,0,[1e3,1e3,1e-6],[0,-HL*1e3,0])
xlim(1e3*[DN1/2,DN2/2-taper_width+0.2e-3])
ylim(1e3*[zpos-HL,zpos-HL+0.2e-3])
caxis([0 9e7*1e-6])
ylabel('z position (mm)')
xlabel('x position (mm)')






