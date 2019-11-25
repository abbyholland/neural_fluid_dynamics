%%% 492 Cerebrospinal Fluid modelled using Womersley's Pulsation Flow
clc; close all; clear all;
%% Global Variables

rho = 1060; % fluid density of blood [kg/m^3]
mu = 3.5e-3; % dynamic viscosity of blood [Ns/m^2]
grid=64;
timestep=28; 
freq=1.0;
dt = 0.1; % s

%% Relevant Constants from the ADAN Model
% Basilar Artery
ba_diameter = 0.003448; %m %[0.003 0.004 0.005 0.006 0.007]; % artery diameter varies
ba_pressure_norm = [9332.565574 15998.683842 14665.460188 14665.460188 12665.624708 13332.236535 12665.624708 11332.401054 10665.789228]; % [Pa] 
ba_pressure_norm = (ba_pressure_norm - mean(ba_pressure_norm))%/1.0e+03;
% ba_pressure_norm =[70 120 110 110 95 100 95 85 80];
ba_pressure_hyp = [15332.072015 23331.413936 23998.025762 23331.413936 21331.578455 19998.354802 18665.131149 17331.907495 16665.295668]; % [Pa]
ba_pressure_hyp = (ba_pressure_hyp - mean(ba_pressure_hyp))%/1.0e+03;

% Posterior Cerebral Artery
pca_diameter = 0.001633; %m %[0.0005 0.001 0.0015 0.002 0.0025 0.003 0.0035 0.004 0.0045 0.005]; % m
pca_pressure_norm = [9332.565574 14665.460188 13998.848361 13998.848361 12532.302343 12798.947073 11999.012881 11332.401054 10665.789228]; % [Pa] left
pca_pressure_norm = (pca_pressure_norm - mean(pca_pressure_norm))%/1.0e+03;
% pca_pressure_norm = [75 70 110 105 105 94 96 90 85 80 75];
pca_pressure_hyp = [14665.460188 22664.802109 23331.413936 22664.802109 21331.578455 19998.354802 18665.131149 16665.295668 15998.683842]; % [Pa]
pca_pressure_hyp = (pca_pressure_hyp - mean(ba_pressure_hyp))%/1.0e+03;

% Hippocampal Artery
ha_diameter = [2e-4 3e-4 4e-4 5e-4 6e-4 7e-4 8e-4]; %artery diameter; varies
ha_pressure = [] ;

%% Basilar Artery
% local variables
rad = linspace(-ba_diameter,ba_diameter, grid);
t = 2;

% normal
[p0,pn,phi] = FourierSeries(ba_pressure_norm,dt,4,1/freq);
[ba_n_u,ba_n_p,ba_n_ta,ba_n_q,dq,alpha]=PulsatileFlow(ba_diameter,rho,mu,freq,p0,pn,phi,timestep,grid);
% velocity feild definitions
n_u1 = permute(ba_n_u, [2 3 1]); % x y t
n_u1 = n_u1(:,:,t);
n_u2 = permute(ba_n_u, [3 2 1]); % y x t 
n_u2 = n_u2(:,:,t);

% hypertensive
[p0,pn,phi] = FourierSeries(ba_pressure_hyp,dt,4,1/freq);
[ba_h_u,ba_h_p,ba_h_ta,ba_h_q,dq,alpha]=PulsatileFlow(ba_diameter,rho,mu,freq,p0,pn,phi,timestep,grid);
% velocity feild definitions
h_u1 = permute(ba_h_u, [2 3 1]);
h_u1 = h_u1(:,:,1);
h_u2 = permute(ba_h_u, [3 2 1]);
h_u2 = h_u2(:,:,1);

% plot the volumetric flow over time
figure(1)
plot(ba_n_q)
hold on 
plot(ba_h_q)
title('Basilar Artery: Volumetric Flow over Time')
xlabel('Timestep')
ylabel('Volumetric Flow [m/s]')
legend('Normal Conditions', 'Hypertensive Conditions')

% plot the pressure over time
figure(2)
plot(ba_n_p)
hold on
plot(ba_h_p)
title('Basilar Artery: Pressure over Time')
xlabel('Timestep')
ylabel('Pressure')
legend('Normal Conditions', 'Hypertensive Conditions')

% plot the shear stress
figure(3)
plot(ba_n_ta)
hold on
plot(ba_h_ta)
title('Shear Stress Acting on the Wall over Time')
xlabel('Timestep')
ylabel('Shear Stress')
legend('Normal Conditions', 'Hypertensive Conditions')

% present velocity at given time, x , y
%   2D plots for changing x, constant y
figure(4) % normal
for t = 1:timestep
%     subplot(1,2,1)
%     title('Normal Conditions: X For Every Y Value')
%     for y = 1:grid
%         s = permute(ba_n_u, [2 1 3]);
%         s = s(:,t,y);
%         plot(rad,s,'g');
%         hold on
%     end
%     subplot(1,2,2)
    title('Normal Conditions')
    for x = 1:grid
        s = permute(ba_n_u, [3 1 2]);
        s = s(:,t,x);
        plot(s, rad,'r');
        hold on
    end
end

figure(5) % hypertensive
for t = 1:timestep
%     subplot(1,2,1)
%     title('Hypertensive Conditions: X For Y')
%     for y = 1:grid
%         s = permute(ba_h_u, [2 1 3]);
%         s = s(:,t,y);
%         plot(s, rad,'g');
%         hold on
%     end
%     subplot(1,2,2)
    title('Hypertensive Conditions')
    for x = 1:grid
        s = permute(ba_h_u, [3 1 2]);
        s = s(:,t,x);
        plot(s, rad,'r');
        hold on
    end
end

% Change in Vorticity as a Function of Pulsation

[vortex, info] = IDvortex(rad*1000, rad*1000, n_u1, n_u2) % x,y coordinate vectors; Vx, Vy fields

% Plot the U field with the found vortex for Normal conditions
figure(6); surface('ZData',n_u1,'YData',rad,'XData',rad,'CData',n_u1,'FaceColor','interp','EdgeColor','none'); 
title('Velocity Feild at Singular Point in Time for Normal Conditions')
hold on; 
plot3(vortex(:,1),vortex(:,2),9e9,'g+','MarkerSize',20)

figure(7)
for t = 1:timestep
    n_u1 = permute(ba_n_u, [2 3 1]);
    n_u1 = n_u1(:,:,t);
    surface('ZData',n_u1,'YData',rad,'XData',rad,'CData',n_u1,'FaceColor','interp','EdgeColor','none')
    title('Velocity Feild Over Time for Normal Conditions')
    pause(0.5)
    saveas(gcf,sprintf('ba_norm/%0.f.png', t))
end

[vortex, info] = IDvortex(rad*1000, rad*1000, h_u1, h_u2) %x,y coordinate vectors; Vx, Vy fields

% Plot the U field with the found vortex for Hypertensive conditions
figure(8); 
surface('ZData',h_u1,'YData',rad,'XData',rad,'CData',h_u1,'FaceColor','interp','EdgeColor','none'); 
title('Velocity Feild at Singular Point in Time for Hypertensive Conditions')
hold on; plot3(vortex(:,1),vortex(:,2),9e9,'g+','MarkerSize',20)

figure(9)
for t = 1:timestep
    h_u1 = permute(ba_h_u, [2 3 1]);
    h_u1 = h_u1(:,:,t);
    surface('ZData',h_u1,'YData',rad,'XData',rad,'CData',h_u1,'FaceColor','interp','EdgeColor','none')
    title('Velocity Feild Over Time for Hypertensive Conditions')
    pause(0.5)
    saveas(gcf,sprintf('ba_hyp/%0.f.png', t))
end


%% Posterior Cerebral Artery
% local variables
rad = linspace(-pca_diameter,pca_diameter, grid);
t = 2;

% normal
[p0,pn,phi] = FourierSeries(pca_pressure_norm,dt,4,1/freq);
[pca_n_u,pca_n_p,pca_n_ta,pca_n_q,dq,alpha]=PulsatileFlow(pca_diameter,rho,mu,freq,p0,pn,phi,timestep,grid);
% velocity feild definitions
n_u1 = permute(pca_n_u, [2 3 1]); % x y t
n_u1 = n_u1(:,:,t);
n_u2 = permute(pca_n_u, [3 2 1]); % y x t 
n_u2 = n_u2(:,:,t);

% hypertensive
[p0,pn,phi] = FourierSeries(ba_pressure_hyp,dt,4,1/freq);
[pca_h_u,pca_h_p,pca_h_ta,pca_h_q,dq,alpha]=PulsatileFlow(pca_diameter,rho,mu,freq,p0,pn,phi,timestep,grid);
% velocity feild definitions
h_u1 = permute(pca_h_u, [2 3 1]);
h_u1 = h_u1(:,:,1);
h_u2 = permute(pca_h_u, [3 2 1]);
h_u2 = h_u2(:,:,1);

% plot the volumetric flow over time
figure(10)
plot(pca_n_q)
hold on 
plot(pca_h_q)
title('Posterior Cerebral Artery: Volumetric Flow over Time')
xlabel('Timestep')
ylabel('Volumetric Flow [m/s]')
legend('Normal Conditions', 'Hypertensive Conditions')

% plot the pressure over time
figure(11)
plot(pca_n_p)
hold on
plot(pca_h_p)
title('Posterior Cerebral Artery: Pressure over Time')
xlabel('Timestep')
ylabel('Pressure')
legend('Normal Conditions', 'Hypertensive Conditions')

% plot the shear stress
figure(12)
plot(pca_n_ta)
hold on
plot(pca_h_ta)
title('Shear Stress Acting on the Wall over Time')
xlabel('Timestep')
ylabel('Shear Stress')
legend('Normal Conditions', 'Hypertensive Conditions')

% present velocity at given time, x , y
%   2D plots for changing x, constant y
figure(13) % normal
for t = 1:timestep
%     subplot(1,2,1)
%     title('Normal Conditions: X For Y')
%     for y = 1:grid
%         s = permute(pca_n_u, [2 1 3]);
%         s = s(:,t,y);
%         plot(s,rad,'g');
%         hold on
%     end
%     subplot(1,2,2)
    title('Normal Conditions')
    for x = 1:grid
        s = permute(pca_n_u, [3 1 2]);
        s = s(:,t,x);
        plot(s, rad, 'r');
        hold on
    end
end

figure(14) % hypertensive
for t = 1:timestep
%     subplot(1,2,1)
%     title('Hypertensive Conditions: X for Y')
%     for y = 1:grid
%         s = permute(pca_h_u, [2 1 3]);
%         s = s(:,t,y);
%         plot(s, rad,'g');
%         hold on
%     end
%     subplot(1,2,2)
    title('Hypertensive Conditions')
    for x = 1:grid
        s = permute(pca_h_u, [3 1 2]);
        s = s(:,t,x);
        plot(s, rad, 'r');
        hold on
    end
end

% Change in Vorticity as a Function of Pulsation

[vortex, info] = IDvortex(rad*1000, rad*1000, n_u1, n_u2) % x,y coordinate vectors; Vx, Vy fields

% Plot the U field with the found vortex for Normal conditions
figure(15); surface('ZData',n_u1,'YData',rad,'XData',rad,'CData',n_u1,'FaceColor','interp','EdgeColor','none'); 
title('Velocity Feild at Singular Point in Time for Normal Conditions')
hold on; 
plot3(vortex(:,1),vortex(:,2),9e9,'g+','MarkerSize',20)

figure(16)
for t = 1:timestep
    n_u1 = permute(pca_n_u, [2 3 1]);
    n_u1 = n_u1(:,:,t);
    surface('ZData',n_u1,'YData',rad,'XData',rad,'CData',n_u1,'FaceColor','interp','EdgeColor','none')
    title('Velocity Feild Over Time for Normal Conditions')
    pause(0.5)
    saveas(gcf,sprintf('pca_norm/%0.f.png', t))
end

[vortex, info] = IDvortex(rad*1000, rad*1000, h_u1, h_u2) %x,y coordinate vectors; Vx, Vy fields

% Plot the U field with the found vortex for Hypertensive conditions
figure(17); 
surface('ZData',h_u1,'YData',rad,'XData',rad,'CData',h_u1,'FaceColor','interp','EdgeColor','none'); 
title('Velocity Feild at Singular Point in Time for Hypertensive Conditions')
hold on; plot3(vortex(:,1),vortex(:,2),9e9,'g+','MarkerSize',20)

figure(18)
for t = 1:timestep
    h_u1 = permute(pca_h_u, [2 3 1]);
    h_u1 = h_u1(:,:,t);
    surface('ZData',h_u1,'YData',rad,'XData',rad,'CData',h_u1,'FaceColor','interp','EdgeColor','none')
    title('Velocity Feild Over Time for Hypertensive Conditions')
    pause(0.5)
    saveas(gcf,sprintf('pca_hyp/%0.f.png', t))
end


%% other potential functions / things that can be used or done
% print the angle as a movie over time
%figure(3)
%SurfMovie(u)
% save the produced values as a DICOM image
%folder='ba_normal';
%date = 20191123;
%[prefix,num,venc,pixel_size] = FlowDICOM(u,ba_diameter,folder)
% read and present the DICOM image over time
%present_dicom(folder, date)


% %   3D representation of velocity over time
% figure(4) % normal
% for t = 1:timestep
%     for g = 1:grid
%         sx = permute(ba_n_u, [3 1 2]);
%         sx = sx(:,t,g);
%         sy = permute(ba_n_u, [2 1 3]);
%         sy = sy(:,t,g);
%         plot3(rad, sx, sy,'g');
%         hold on
%     end
% end
% figure(5) % hypertensive
% for t = 1:timestep
%     for g = 1:grid
%         sx = permute(ba_h_u, [3 1 2]);
%         sx = sx(:,t,g);
%         sy = permute(ba_h_u, [2 1 3]);
%         sy = sy(:,t,g);
%         plot3(rad, sx, sy,'g');
%         hold on
%     end
% end


