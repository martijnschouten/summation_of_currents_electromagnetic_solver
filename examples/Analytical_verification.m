clear all
close all

u0 = 4*pi*1e-7;

A = 1e-3;%radius first loop
a = 1e-3;%radius second loop
b = 2e-3;%distance between the loops
f = 2e6;%frequency of the current through the first lop
t = 10e-6;%thickness of the round wire
RhoB = 1/0.28*1.68e-8;%resistivity of the wire

%calculate mutual inductance using Maxwell's method
c = 2*sqrt(A*a)/sqrt((A+a)^2+b^2);
K = ellipticK(c^2);
E = ellipticE(c^2);
LMm = -1*u0*sqrt(A*a)*((c-2/c)*K+2/c*E)

%calculate the self inductance
LNm2 = u0*A*((3/16/exp(0.5)*t^2/a^2+1)*log(8*a/t)-1/64/exp(0.5)*t^2/a^2-7/4)

%calculate the resistance
Rm = 2*a*RhoB/t^2

%calculate the change in impedance
jw = 1i*f*2*pi;
dZm = -(f*2*pi)^2*LMm^2/(Rm-LNm2*jw)

%calculate the change in parallel inductance
dLpm = (imag(dZm).^2+real(dZm).^2)./imag(dZm)./f/2/pi

%calculate the change in parallel resistance
dRpm = (imag(dZm).^2+real(dZm).^2)./real(dZm)


%parameters for the summation of currents simulation
RhoC = 1.68e-8;%ohm meter
NL = 1;%number of windings of coil 
HL = 2*t;%m height of coil
DL1 = 2*a-2*t;%m inner diameter of coil
DL2 = 2*a+2*t;%m outer diameter of coil

nxL = 4;%number of elements on the coil in the x direction
nzL = 4;%number of elements on the coil in the z direction
nphiL = 400;%number of elements of the coil in radial direction


NN = 1;%number of windings of nozzle
HN = 2*t;%m height of nozzle
DN1 = 2*A-2*t;%m inner diameter of nozzle
DN2 = 2*A+2*t;%m outer diameter of nozzle

zpos = b;

nxN = 4;%number of elements on the coil in the x direction
nzN = 4;%number of elements on the coil in the z direction
nphiN = nphiL;%number of elements of the coil in radial direction

xpos = 0;

for i1 = 1:length(nphiL)
    sense_coil = summation_of_currents;
    sense_coil = sense_coil.set_geometry(HL,DL1,DL2,NL,RhoC,0,0);
    sense_coil = sense_coil.set_mesh(nxL,nzL,nphiL(i1));
    sense_coil = sense_coil.build_single_loop();

    nozzle = summation_of_currents;
    nozzle = nozzle.set_geometry(HN,DN1,DN2,NN,RhoB,xpos,zpos);
    nozzle = nozzle.set_mesh(nxN,nzN,nphiN(i1));
    nozzle = nozzle.build_single_loop();

    if i1 == length(nphiL)
        figure
        sense_coil.plot_geometry(3,[1,1,1]);
        hold on
        nozzle.plot_geometry(3,[1,1,1]);

        figure
        nozzle.plot_geometry(3,[1,1,1]);
    end

    
    dZs(i1) = sense_coil.calculate_impedance_change(nozzle,f);

    LMs(i1) = sense_coil.calculate_mutual_inductance(nozzle);

    LNs(i1) = nozzle.calculate_self_inductance();

    Rms(i1) = nozzle.resistance();

    dLps(i1) = (imag(dZs(i1)).^2+real(dZs(i1)).^2)./imag(dZs(i1))./f/2/pi;
    dRps(i1) = (imag(dZs(i1)).^2+real(dZs(i1)).^2)./real(dZs(i1));
end
