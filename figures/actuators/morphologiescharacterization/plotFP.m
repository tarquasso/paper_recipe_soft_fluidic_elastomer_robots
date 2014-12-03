function plotFP()
% cyl = cylindrical
% rib = ribbed
% ple = pleated

KILOPASCAL_PER_VOLT = (15.0/4.0)*6.89475729; % [kPa/V] for the pressure transducer
CENTIMETERS_PER_BIT = 10.16/4096; % [cm/'number'] for the jrk positional feedback
PISTON_DIAMETER = 6.35; % [cm] 
PISTON_AREA = pi*(PISTON_DIAMETER/2)^2; % [cm^2] 
FORCE_PER_VOLT = 0.6642; %[N/V]

cyl_voltage_pressure = [0.545, 0.600, 0.670, 0.730, 0.80, 0.90, 1.00, 1.07, 1.19, 1.30, 1.40, 1.50, 1.60, 1.68, 1.76, 1.76, 1.78, 1.87, 1.92, 2.00, 2.08, 2.18, 2.30, 2.39, 2.48];
cyl_pressure = (cyl_voltage_pressure - cyl_voltage_pressure(1))*KILOPASCAL_PER_VOLT; % [kPa]
cyl_volume = (0:50:1200)*CENTIMETERS_PER_BIT*PISTON_AREA; % [cm^3] or [mL]
cyl_voltage_force = [0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0296, 0.0283, 0.0268, 0.0224, 0.017, 0.008, -0.007, -0.038, -0.076, -0.320, -0.773, -1.013, -1.200, -1.366, -1.536, -1.694, -1.835, -1.970, -2.106, -2.256 ];
cyl_force = abs(cyl_voltage_force-cyl_voltage_force(1))*FORCE_PER_VOLT;
cyl_energy = calculateEnergy(cyl_volume, cyl_pressure);

% the curvature plot
p1 = figure(4);
a1 = axes('Parent',p1,'FontSize',14);
grid(a1,'on');
hold(a1,'all');
ylabel('Tip Force [N]', 'FontSize',14) % label left y-axis
xlabel('Fluid Energy [J]', 'FontSize',14) % label x-axis
h1 = plot(cyl_energy, cyl_force);
set(h1, 'Color', 'k');
set(h1, 'LineStyle', '-');
set(h1, 'LineWidth', 2.0);
set(h1, 'Marker', 'o');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rib_voltage_pressure = [0.545, 0.615, 0.680, 0.725, 0.807, 0.90, 1.01, 1.04, 1.14, 1.24, 1.35, 1.44, 1.53, 1.60, 1.66, 1.68, 1.70, 1.76, 1.81, 1.84, 1.92, 1.94];
rib_pressure = (rib_voltage_pressure - rib_voltage_pressure(1))*KILOPASCAL_PER_VOLT; % [kPa]
rib_volume = (0:50:1050)*CENTIMETERS_PER_BIT*PISTON_AREA; % [cm^3] or [mL]
rib_voltage_force = [-0.1073, -0.1242, -0.1385, -0.157, -0.175, -0.1955, -0.218, -0.246, -0.2703, -0.305, -0.347, -0.410, -0.494, -0.600, -0.727, -0.855, -0.967, -1.072, -1.172, -1.279, -1.355, -1.441 ];
rib_force = abs(rib_voltage_force-rib_voltage_force(1))*FORCE_PER_VOLT;
rib_energy = calculateEnergy(rib_volume, rib_pressure);

figure(4);
h1 = plot(rib_energy, rib_force);
set(h1, 'Color', 'k');
set(h1, 'LineStyle', '-');
set(h1, 'LineWidth', 2.0);
set(h1, 'Marker', '+');
set(h1, 'MarkerSize', 10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ple_voltage_pressure = [0.540, 0.660, 0.800, 0.980, 1.12, 1.32, 1.52, 1.70, 1.92, 2.16, 2.40, 2.64, 2.82, 3.05, 3.20, 3.45, 3.68, 3.93];
ple_pressure = (ple_voltage_pressure - ple_voltage_pressure(1))*KILOPASCAL_PER_VOLT; % [kPa]
ple_volume = (0:100:1700)*CENTIMETERS_PER_BIT*PISTON_AREA; % [cm^3] or [mL]
ple_voltage_force =[ 0.04521, 0.01235, -0.01099, -0.05355, -0.1055, -0.1513, -0.1935, -0.2545, -0.3344, -0.4469, -0.5764, -0.7590, -0.9706, -1.258, -1.589, -2.00, -2.437, -2.887 ];
ple_force = abs(ple_voltage_force-ple_voltage_force(1))*FORCE_PER_VOLT;
ple_energy = calculateEnergy(ple_volume, ple_pressure);

figure(4);
h1 = plot(ple_energy, ple_force);
set(h1, 'Color', 'k');
set(h1, 'LineStyle', '-');
set(h1, 'LineWidth', 2.0);
set(h1, 'Marker', 's');
legend({'Cylindrical', 'Ribbed', 'Pleated'}, 'Location','northwest', 'FontSize',12)

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

