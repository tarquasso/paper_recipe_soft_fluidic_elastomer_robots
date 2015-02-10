function plotCharacterizations_colored()
% cyl = cylindrical
% rib = ribbed
% ple = pleated
close all

KILOPASCAL_PER_VOLT = (15.0/4.0)*6.89475729; % [kPa/V] for the pressure transducer
CENTIMETERS_PER_BIT = 5.0/4096; % [cm/'number'] for the jrk positional feedback
PISTON_DIAMETER = 3.81; % [cm]
PISTON_AREA = pi*(PISTON_DIAMETER/2)^2; % [cm^2]

cyl_voltage = [0.51, 0.660, 0.832, 0.980, 1.14, 1.32, 1.48, 1.68, 1.76, 1.77, 1.80, 1.85, 1.92, 2.03, 2.10, 2.21, 2.33, 2.42, 2.53, 2.62];
cyl_pressure = (cyl_voltage - cyl_voltage(1))*KILOPASCAL_PER_VOLT; % [kPa]
cyl_curvature = [0.0, 0.014, 0.035, 0.071, 0.132, 0.249, 0.472, 1.34, 4.74, 7.25, 8.97, 10.83, 12.56, 14.23, 15.95, 17.42, 18.83, 20.18, 21.40, 22.58]; % [1/m]
cyl_volume = (0:200:3800)*CENTIMETERS_PER_BIT*PISTON_AREA; % [cm^3] or [mL]
cyl_energy = calculateEnergy(cyl_volume, cyl_pressure);
cyl_denergy = calculateEnergyError(cyl_volume, cyl_pressure);

% the curvature plot
p1 = figure(1);
a1 = axes('Parent',p1,'FontSize',14);
grid(a1,'on');
hold(a1,'all');
ylabel('Curvature [m^{-1}]', 'FontSize',14) % label left y-axis
xlabel('Volume [mL]', 'FontSize',14) % label x-axis
h1 = plot(cyl_volume, cyl_curvature);
set(h1, 'Color', 'k');
set(h1, 'LineStyle', '-');
set(h1, 'LineWidth', 2.0);
set(h1, 'Marker', 'o');

% the pressure plot
p2 = figure(2);
a2 = axes('Parent',p2,'FontSize',14);
grid(a2,'on');
hold(a2,'all');
ylabel('Pressure [kPa]', 'FontSize',25) % label right y-axis
xlabel('Volume [mL]', 'FontSize',25) % label x-axis
h2 = plot(cyl_volume, cyl_pressure);
set(h2, 'Color', 'k');
set(h2, 'LineStyle', '-');
set(h2, 'LineWidth', 2.0);
set(h2, 'Marker', 'o');

% the energy plot
p3 = figure(3);
a3 = axes('Parent',p3,'FontSize',14);
grid(a3,'on');
hold(a3,'all');
xlabel('Fluid Energy [J]', 'FontSize',14) % label right y-axis
ylabel('Curvature [m^{-1}]', 'FontSize',14) % label x-axis
h3 = plot(cyl_energy, cyl_curvature);
set(h3, 'Color', 'k');
set(h3, 'LineStyle', '-');
set(h3, 'LineWidth', 2.0);
set(h3, 'Marker', 'o');

% the error_energy plot
p4 = figure(4);
a4 = axes('Parent',p4,'FontSize',14);
grid(a4,'on');
hold(a4,'all');
ylabel('Fluid Energy Error [J]', 'FontSize',14) % label right y-axis
xlabel('Fluid Energy  J]', 'FontSize',14) % label x-axis
h4 = plot(cyl_energy,cyl_denergy);
set(h4, 'Color', 'k');
set(h4, 'LineStyle', '-');
set(h4, 'LineWidth', 2.0);
set(h4, 'Marker', 'o');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rib_voltage = [0.515, 0.65, 0.805, 0.95, 1.07, 1.16, 1.34, 1.40, 1.50, 1.50, 1.60, 1.65, 1.67, 1.74, 1.77];
rib_pressure = (rib_voltage - rib_voltage(1))*KILOPASCAL_PER_VOLT; % [kPa]
rib_curvature = [0.0, 0.468, 1.024, 1.661, 2.46, 3.60, 5.12, 7.29, 10.24, 13.24, 16.47, 20.33, 23.18, 26.75, 29.30]; % [1/m]
rib_volume = (0:200:2800)*CENTIMETERS_PER_BIT*PISTON_AREA; % [cm^3] or [mL]
rib_energy = calculateEnergy(rib_volume, rib_pressure);
rib_denergy = calculateEnergyError(rib_volume, rib_pressure);

figure(1);
h1 = plot(rib_volume, rib_curvature);
set(h1, 'Color', 'r');
set(h1, 'LineStyle', '-');
set(h1, 'LineWidth', 2.0);
set(h1, 'Marker', '+');
set(h1, 'MarkerSize', 10);

figure(2);
h2 = plot(rib_volume, rib_pressure);
set(h2, 'Color', 'r');
set(h2, 'LineStyle', '-');
set(h2, 'LineWidth', 2.0);
set(h2, 'Marker', '+');
set(h2, 'MarkerSize', 10);


figure(3);
h3 = plot(rib_energy, rib_curvature);
set(h3, 'Color', 'r');
set(h3, 'LineStyle', '-');
set(h3, 'LineWidth', 2.0);
set(h3, 'Marker', '+');
set(h3, 'MarkerSize', 10);

figure(4);
h4 = plot(rib_energy, rib_denergy);
set(h4, 'Color', 'r');
set(h4, 'LineStyle', '-');
set(h4, 'LineWidth', 2.0);
set(h4, 'Marker', '+');
set(h4, 'MarkerSize', 10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ple_voltage = [0.52, 0.675, 0.770, 0.920, 1.05, 1.19, 1.31, 1.47, 1.65, 1.82, 2.0, 2.16, 2.32, 2.45, 2.60, 2.67, 2.81, 2.90];
ple_pressure = (ple_voltage - ple_voltage(1))*KILOPASCAL_PER_VOLT; % [kPa]
ple_curvature = [0, 0.554, 1.146, 1.712, 2.234, 2.98, 3.73, 4.67, 5.586, 6.982, 8.339, 9.911, 11.55, 13.49, 15.39, 16.22, 17.10, 17.91]; % [1/m]
ple_volume = (0:200:3400)*CENTIMETERS_PER_BIT*PISTON_AREA; % [cm^3] or [mL]
ple_energy = calculateEnergy(ple_volume, ple_pressure);
ple_denergy = calculateEnergyError(ple_volume, ple_pressure);

figure(1);
h1 = plot(ple_volume, ple_curvature);
set(h1, 'Color', 'b');
set(h1, 'LineStyle', '-');
set(h1, 'LineWidth', 2.0);
set(h1, 'Marker', 's');
legend({'Cylindrical', 'Ribbed', 'Pleated'}, 'Location','northwest', 'FontSize',12)

figure(2);
h2 = plot(ple_volume, ple_pressure);
set(h2, 'Color', 'b');
set(h2, 'LineStyle', '-');
set(h2, 'LineWidth', 2.0);
set(h2, 'Marker', 's');
legend({'Cylindrical', 'Ribbed', 'Pleated'}, 'Location','northwest','FontSize',22)

figure(3);
h3 = plot(ple_energy, ple_curvature);
set(h3, 'Color', 'b');
set(h3, 'LineStyle', '-');
set(h3, 'LineWidth', 2.0);
set(h3, 'Marker', 's');
legend({'Cylindrical', 'Ribbed', 'Pleated'}, 'Location','northwest','FontSize',12)

figure(4);
h4 = plot(ple_energy, ple_denergy);
set(h4, 'Color', 'b');
set(h4, 'LineStyle', '-');
set(h4, 'LineWidth', 2.0);
set(h4, 'Marker', 's');
legend({'Cylindrical', 'Ribbed', 'Pleated'}, 'Location','northwest','FontSize',12)

end

function energy = calculateEnergy(volume, pressure)
volume = volume*1.0e-6; % [m^3]
pressure = pressure*1000;
max_index = length(volume);
energy = zeros(1,max_index);

for i=1:length(volume)
    if(i==1)
        energy = 0;
    else
        delta_volume = volume(i)-volume(i-1);
        delta_pressure = pressure(i)-pressure(i-1);
        energy(i) = energy(i-1) + delta_volume*pressure(i-1) + 0.5*delta_pressure*delta_volume;
    end
end
end


function denergy = calculateEnergyError(volume, pressure)
volume = volume*1.0e-6; % [m^3]
pressure = pressure*1000; %[Pa]
max_index = length(volume);
denergy = zeros(1,max_index);
perr = 0.02*15*6894.78; % [Pa] ASDXACX015PAAA5 Accuracy 2% of 15psi full scale
verr = 0.3/1000*(1.5*0.0254)^2; % [m^3] http://www.firgelli.com/Uploads/L16_datasheet.pdf 3mm accuracy

for i=1:length(volume)
    if(i==1)
        denergy(i)=0;
        
    else
        delta_volume = volume(i)-volume(i-1);
        delta_pressure = pressure(i)-pressure(i-1);
        deltap_err = additiveError(perr,perr,0);
        deltav_err = additiveError(verr,verr,0);
        
        dvdelpm1 = multiplicativeError(delta_volume,deltav_err,pressure(i-1),perr);
        if(isnan(dvdelpm1))
            dvdelpm1=0;
        end
        
        dpdelvdel = multiplicativeError(delta_pressure,deltap_err,delta_volume,deltav_err);     
        if(isnan(dpdelvdel))
            dpdelvdel=0;
        end
        
        denergy(i) = additiveError(denergy(i-1),dvdelpm1,0.5*dpdelvdel);
    end
end
end

function error = additiveError(xerr,yerr,zerr)
error = sqrt(xerr^2 + yerr^2 + zerr^2);        
end

function error = multiplicativeError(x,xerr,y,yerr)
R = abs(x*y);
error = R * sqrt((xerr/x)^2 + (yerr/y)^2);        
end

