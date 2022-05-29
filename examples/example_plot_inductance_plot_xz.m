clear all
close all

addpath('../')

zpos = 2e-3:0.25e-3:3.5e-3;
load("sim_result_xz2.mat")
zpos = zpos - sense_coil.H;


figure
ax1 = subplot(2,1,1);
s = surf(xpos*1e3,zpos*1e3,real(dLp)*1e9);
ylabel('\Delta Inductance (H)')
xlabel('Nozzle position (mm)')
view(2)
s.EdgeColor = 'none';
s.FaceColor = 'interp';

xlabel('x (mm)')
ylabel('z (mm)')
c = colorbar;
c.Label.String = '\Delta Parallel inductance (nH)';
xlimits = [min(xpos)*1e3,max(xpos)*1e3];
ylimits = [10e-1,max(zpos*1e3)];
xlim(xlimits)
ylim(ylimits)
climits = 1e9*[min(min(real(dLp(3:end,:))));max(max(real(dLp(3:end,:))))];
caxis(climits);

hold on
dLlvec = zeros(length(zpos)*length(xpos),1);
xposvec = zeros(length(zpos)*length(xpos),1);
zposvec = zeros(length(zpos)*length(xpos),1);
for i1=1:length(zpos)
    loc = (i1-1)*length(xpos)+1:i1*length(xpos);
    zposvec(loc) = zpos(i1);
    xposvec(loc) = xpos;
end
    
scatter3(xposvec*1e3,zposvec*1e3,1e9*max(max(dLp))*ones(length(zposvec),1),8,'r','filled')

ax2 = subplot(2,1,2);
s = surf(xpos*1e3,zpos*1e3,dRp);
ylabel('\Delta Parallel resistance (\Omega)')
xlabel('Nozzle position (mm)')
view(2)
s.EdgeColor = 'none';
s.FaceColor = 'interp';

xlabel('x (mm)')
ylabel('z (mm)')
c = colorbar;
c.Label.String = '\Delta Parallel resistance (\Omega)';
xlim(xlimits)
ylim(ylimits)

climits2 = [min(min(real(dRp(3:end,:))));max(max(real(dRp(3:end,:))))];
caxis(climits2);

hold on
scatter3(xposvec*1e3,zposvec*1e3,max(max(dRp))*ones(length(zposvec),1),8,'r','filled')

set(gcf,'Position',[0,300,450,600])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
exportgraphics(gcf,'../../inductive 3d printer calibration 2/Images/inductance_xz.pdf','ContentType','vector');

zuse = 1.1e-3;
[~,nuse] = min(abs(zpos-zuse));
dLcomp_sim = dLp(nuse,:);
xcomp_sim = xpos;

load("diabase_xz_1.mat")

x = [];
z = [];
L = [];



for i1 = 1:length(zpos)
    use = xpos(:,i1)~=0;
    x = [x;xpos(use,i1)];
    z = [z;zpos(i1)*ones(sum(use),1)];
    L = [L;data(use,i1)];
end
x = (x-mean(x))*1e-3;
z = (z-6.1)*1e-3;

z = z +6e-4;
figure
nozzle.plot_trisurf(x*1e3,z*1e3,1e9*(L-max(L)))
xlabel('x (mm)')
ylabel('z (mm)')
c = colorbar;
c.Label.String = '\Delta Parallel inductance (nH)';
xlim(xlimits)
ylim(ylimits)
caxis(climits);



set(gcf,'Position',[0,300,450,300])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
exportgraphics(gcf,'../../inductive 3d printer calibration 2/Images/inductance_xz_meas.pdf','ContentType','vector');

[~,nuse] = min(abs(zpos-6.1+0.6-zuse));
use = xpos(:,nuse)~=0;
dLcomp_meas = data(use,nuse);
xcomp_meas = xpos(use,nuse)-mean(xpos(use,nuse));

disp('inductance:')
disp(min(dLcomp_meas))

xcomp_sim = xcomp_sim*1e3;
dLcomp_sim = dLcomp_sim*1e9;
dLcomp_meas = (dLcomp_meas-max(dLcomp_meas))*1e9;
figure
plot(xcomp_sim,dLcomp_sim);
hold on
plot(xcomp_meas,dLcomp_meas);
ylabel('\Delta Inductance (nH)')
xlabel('Position (mm)')
xlim([-3,3])
leg = legend('simulated','measured');


save('inductance_comparison_data.mat','xcomp_sim', 'dLcomp_sim', 'xcomp_meas', 'dLcomp_meas')

for i1 = 1:length(zpos)
    use = xpos(:,i1)~=0;
    
    xtemp = xpos(use,i1);
    Ltemp = data(use,i1)*1e9;
    
    y_min = min(real(Ltemp));
    y_max = max(real(Ltemp));
    x_avg = mean(xtemp);
    x_min = min(xtemp);
    b0 = (y_max-y_min)/(x_min-x_avg)^2;
    p0 = double([x_avg y_min b0 0 0 0]);
    
    fun = @(p)abs(Ltemp-(p(2) + p(3)*(xtemp-p(1)).^2 + p(4)*(xtemp-p(1)).^4 + p(5)*(xtemp-p(1)).^6 + p(6)*(xtemp-p(1)).^8));
    fun(p0);

    pfit = lsqnonlin(fun,p0);
    pos(i1) = pfit(1);
   
end

dist = xtemp(end)-xtemp(1)
ttemp =  time(use,i1)
time_per_meas = ttemp(end)-ttemp(1)
v = dist/time_per_meas

figure
zpos1 = zpos(1:2:end)-6.1+0.6-zuse;
zpos2 = zpos(2:2:end)-6.1+0.6-zuse;
plot(zpos1,pos(1:2:end),zpos2,pos(2:2:end))
xlim([min(zpos1) max(zpos2)])


yline(min(pos(1:2:end)),':','Color',[0, 0.4470, 0.7410])
yline(max(pos(1:2:end)),':','Color',[0, 0.4470, 0.7410])
yline(min(pos(2:2:end)),':','Color',[0.8500, 0.3250, 0.0980])
yline(max(pos(2:2:end)),':','Color',[0.8500, 0.3250, 0.0980])

xlabel('z position (mm)')
ylabel('Point of symmetry (mm)')
legend('Moving up','Moving down')

set(gcf,'Position',[0,300,450,300])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
exportgraphics(gcf,'../../inductive 3d printer calibration 2/Images/inductance_xz_meas_pos.pdf','ContentType','vector');



