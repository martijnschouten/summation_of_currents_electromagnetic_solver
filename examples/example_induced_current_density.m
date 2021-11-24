clear all
close all
close all force

addpath("../")

RhoC = 1.68e-8;%ohm meter
NL = 55;%number of windings of coil 
HL = 1.4e-3;%m height of coil
DL1 = 1.84e-3;%m inner diameter of coil
DL2 = 2e-3;%m outer diameter of coil

nxL = 3;%number of elements on the coil in the x direction
nzL = 20;%number of elements on the coil in the z direction
nphiL = 120;%number of elements of the coil in radial direction

RhoB = 1/0.28*1.68e-8;%ohm meter
HN = 1e-3;%m height of nozzle
DN1 = 0.4e-3;%m inner diameter of nozzle
DN2 = 3e-3;%m outer diameter of nozzle

zpos = 2.5e-3;%Location of the object in z in meter

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

sense_coil = SOC_object;
sense_coil = sense_coil.set_geometry(HL,DL1,DL2,NL,RhoC,0,0);
sense_coil = sense_coil.set_mesh(nxL,nzL,nphiL);
sense_coil = sense_coil.build_coil();
Zself = sense_coil.calculate_impedance(f)
Lp = (imag(Zself).^2+real(Zself).^2)./imag(Zself)./f/2/pi
Rp = (imag(Zself).^2+real(Zself).^2)./real(Zself)



nozzle = SOC_object;
nozzle = nozzle.set_nozzle_geometry(HN,DN1,DN2,RhoB,xpos,zpos,taper_width,taper_angle,hole_offset,rotation);
nozzle = nozzle.set_boundary_mesh(nxN,nzN,nphiN,nxbl,nxbr,nxc,nzbb,nzbt);
skin_depth = nozzle.skin_depth(f);
disp("Skin depth: ")
disp(skin_depth)
nozzle = nozzle.build_nozzle_surface(nmesh*skin_depth);

figure
sense_coil.plot_geometry(3,[1e3,1e3,1e3]);
set(gcf,'Position',[0,300,350,300])
xlabel('x position (mm)')
ylabel('y position (mm)')
zlabel('z position (mm)')

figure
nozzle.plot_geometry(2,[1e3,1e3,1e3]);
xlabel('x position (mm)')
ylabel('y position (mm)')
zlabel('z position (mm)')
set(gcf,'Position',[350,300,350,300])
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
set(gcf,'Position',[300,300,350,500])
xlim([-1.2,1.2])
ylim([-1.2,1.2])
zlim([0,4])

figure
sense_coil.plot_surface_vectors([1e3,1e3,1e3])
 xlabel('x position (mm)')
ylabel('y position (mm)')
zlabel('z position (mm)')
view(2)
set(gcf,'Position',[700,300,350,300])

    
figure
sense_coil.plot_induced_current_density_2d(nozzle,1,f,0,[1e3,1e3,1e-6],[0,-HL*1e3,0])
xlim(1e3*[0,DN2/2])
daspect([1 1 1])
ylim(1e3*[zpos-HL,zpos-HL+HN])
caxis([0 9e7*1e-6])
ylabel('z position (mm)')
xlabel('x position (mm)')
set(gcf,'Position',[1050,300,450,300])
c = colorbar;
c.Label.String = 'Induced current density (MAm^{-2})';



figure
T = table2array(readtable('I_induced.csv'));
sense_coil.plot_trisurf(real(T(:,1))*1e3,real(T(:,2)-HL)*1e3,abs(T(:,3))*1e-6)
xlim(1e3*[0,DN2/2])
daspect([1 1 1])
ylim(1e3*[zpos-HL,zpos-HL+HN])
caxis([0 9e7*1e-6])
ylabel('z position(mm)')
xlabel('x position(mm)')
c = colorbar;
c.Label.String = 'Induced current density (MAm^{-2})';



figure
nxc = 131;
x_comsol = reshape(real(T(:,1)),nxc,length(T(:,1))/nxc);
y_comsol = reshape(real(T(:,2)),nxc,length(T(:,1))/nxc);
J_comsol = reshape(abs(T(:,3)),nxc,length(T(:,1))/nxc);
sense_coil.plot_current_density_difference(nozzle,1,f,y_comsol,x_comsol,J_comsol,0,[1e3,1e3,1e-6],[0,-HL*1e3,0])
daspect([1 1 1])
xlim(1e3*[0,DN2/2])
ylim(1e3*[zpos-HL,zpos-HL+HN])
ylabel('z position(mm)')
xlabel('x position(mm)')
c = colorbar;
c.Label.String = 'Induced current density (MAm^{-2})';


figure
sense_coil.plot_current_density_difference_ratio(nozzle,1,f,y_comsol,x_comsol,J_comsol,0,[1e3,1e3,1],[0,-HL*1e3,0])
daspect([1 1 1])
xlim(1e3*[0,DN2/2])
ylim(1e3*[zpos-HL,zpos-HL+HN])
ylabel('z position(mm)')
xlabel('x position(mm)')
caxis([0 50])

c = colorbar;
c.Label.String = 'Error (%)';








