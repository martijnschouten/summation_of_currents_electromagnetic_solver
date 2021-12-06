clear all
close all

load("sim_result_0mu_offset_hr2.mat")

dLp = (imag(dZ).^2+real(dZ).^2)./imag(dZ)./f/2/pi;
dRp = (imag(dZ).^2+real(dZ).^2)./real(dZ);
dltot = zeros(3,13);
figure
ax1 = subplot(2,1,1);
plot(xpos,dLp)
ylabel('\Delta Inductance (H)')
xlabel('Nozzle position (m)')
hold on
ax2 = subplot(2,1,2);
plot(xpos,dRp)
ylabel('\Delta Parallel resistance (\Omega)')
xlabel('Nozzle position (m)')
hold on
dLptot(1,:) = dLp;
dRptot(1,:) = dRp;

load("sim_result_1deg_rotation_hr2.mat")
dLp = (imag(dZ).^2+real(dZ).^2)./imag(dZ)./f/2/pi;
dRp = (imag(dZ).^2+real(dZ).^2)./real(dZ);
subplot(ax1)
plot(xpos,dLp)
ylabel('\Delta Parallel Inductance (H)')
xlabel('Nozzle position (m)')
hold on
subplot(ax2)
plot(xpos,dRp)
ylabel('\Delta Parallel resistance (\Omega)')
xlabel('Nozzle position (m)')
hold on
dLptot(2,:) = dLp;
dRptot(2,:) = dRp;

load("sim_result_2deg_rotation_hr2.mat")
dLp = (imag(dZ).^2+real(dZ).^2)./imag(dZ)./f/2/pi;
dRp = (imag(dZ).^2+real(dZ).^2)./real(dZ);
subplot(ax1)
plot(xpos,dLp)
ylabel('\Delta Inductance (H)')
xlabel('Nozzle position (m)')
hold on
%legend('no offset','50mu offset', '100mu offset')
subplot(ax2)
plot(xpos,dRp)
ylabel('\Delta Parallel resistance (\Omega)')
xlabel('Nozzle position (m)')
hold on
legend('0 degree','1 degree', '2 degree')
dLptot(3,:) = dLp;
dRptot(3,:) = dRp;
set(gcf,'Position',[0,100,450,600])


style = '- .'

ddLp1 = dLptot(2,:)-dLptot(1,:)
ddRp1 = dRptot(2,:)-dRptot(1,:)
figure
ax1 = subplot(2,1,1);
plot(xpos*1e3,ddLp1*1e12,style)
ylabel('\Delta \Delta Parallel inductance (pH)')
xlabel('Nozzle position (mm)')
hold on
xline(0,':k')
ax2 = subplot(2,1,2);
plot(xpos*1e3,ddRp1*1e3,style)
ylabel('\Delta \Delta Parallel resistance (m\Omega)')
xlabel('Nozzle position (mm)')

ddLp2 = dLptot(3,:)-dLptot(1,:)
ddRp2 = dRptot(3,:)-dRptot(1,:)


hold on
subplot(ax1)
plot(xpos*1e3,ddLp2*1e12,style)
%ylabel('\Delta\Delta Inductance (H)')
%xlabel('Nozzle position (m)')
hold on
subplot(ax2)
plot(xpos*1e3,ddRp2*1e3,style)
%ylabel('\Delta\Delta Parallel resistance (\Omega)')
%xlabel('Nozzle position (m)')
xline(0,':k')

leg = legend('1 degree', '2 degree')
set(gcf,'Position',[0,100,450,600])
%leg.Position = [0.62794503703308,0.44938889569189,0.277135235219664,0.060833334604899]

%export_fig('../Images/inductance_vs_rotation.png', '-dpng', '-transparent', '-r600');

xpos = double(xpos*1000);
dLptot = double(dLptot*1e9);
for i1 = 1:3
    y_min = min(real(dLptot(i1,:)));
    y_max = max(real(dLptot(i1,:)));
    x_avg = mean(xpos);
    x_min = min(xpos);
    b0 = (y_max-y_min)/(x_min-x_avg)^2;
    p0 = double([0 x_avg y_min b0 0 0]);
    
    fun = @(p)(real(dLptot(i1,:))-(p(2) + p(3)*(xpos-p(1)).^2 + p(4)*(xpos-p(1)).^4 + p(5)*(xpos-p(1)).^6 + p(6)*(xpos-p(1)).^8));
    fun(p0)
    pfit = lsqnonlin(fun,p0);
    posL(i1) = pfit(1);
end

%xpos = double(xpos*1000);
dRptot = double(dRptot*1e3);
for i1 = 1:3
    y_min = min(real(dRptot(i1,:)));
    y_max = max(real(dRptot(i1,:)));
    x_avg = mean(xpos);
    x_min = min(xpos);
    b0 = (y_max-y_min)/(x_min-x_avg)^2;
    p0 = double([0 x_avg y_min b0 0 0]);
    
    fun = @(p)(real(dRptot(i1,:))-(p(2) + p(3)*(xpos-p(1)).^2 + p(4)*(xpos-p(1)).^4 + p(5)*(xpos-p(1)).^6 + p(6)*(xpos-p(1)).^8));
    fun(p0)
    pfit = lsqnonlin(fun,p0);
    posR(i1) = pfit(1);
end

