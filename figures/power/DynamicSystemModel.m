clc;
clear all;

MeasuredFlow = importdata('MeasuredFlow.mat');
MeasuredTime = importdata('MeasuredTime.mat');

dt = 0.01;
Tfinal = 2.25;
ein = 12.0;
toPSI = 0.000145037738; % Pascals -> psi
toLPM = 60000; % meters^3 / second to Liters/minute

t = 0:dt:Tfinal;
t_off = 0.02;
off_index = find(t >= t_off,1,'first');
%u = (square(t,30.3)+1)*ein/2 ;
u = [zeros(1,off_index), (square(t(off_index+1:end),30.3)+1)*ein/2];
x0 = [0; 0; 0];

g = 232; %Newtons / Amp
p = 475; %Volts / meter / sec
Mp = 0.185; % kiligrams
Rm = 18.5; % Ohms
Ap = 0.00079; % meters^2
Cf = -2.0*10^(-10); % meters^3 / Pascal
Cc = 1.25*10^(-10); %5.176*10^(-10); % meters^3 / Pascal
mu = 1.983*10^(-5); %kg / (meter*sec)
L = 0.4572; %meters
D = 0.001016; %meters
Rt = 128*mu*L/(pi*D^4); % Pascal / meter^3 / sec

a = [-1*g*p/(Mp*Rm), 0, -Ap/Mp; 0, ((Cf + Cc)*Rt)^-1, -((Cf + Cc)*Rt)^-1; -Ap/Cf, -1/(Cf*Rt), 1/(Cf*Rt)];
b = [g/(Mp*Rm); 0; 0];
c = [1, 0, 0; 0, 1, 0; 0, 0, 1];
d = [0];


sys = ss(a,b,c,d);
[y,t,x] = lsim(sys,u,t,x0);

flow = (x(:,3)-x(:,2))./(Rt); 

%index at compliance switch
p_switch = 2.25*10^(4);
switch_index = find(x(:,2)>=p_switch,1,'first');

figure;
hold on;
plot(t(1:switch_index), flow(1:switch_index)*toLPM, '-r')
plot(t(1:end), zeros(length(t),1), '--k')
plot(MeasuredTime'-(981/1000), MeasuredFlow)
axis([0 Tfinal -1 2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_after = t(switch_index):dt:Tfinal;
u = (square(t_after,30.3)+1)*ein/2 ;
x0 = [x((switch_index),1); x((switch_index),2); x((switch_index),3)];

Cc = 4.15*10^(-9); %5.176*10^(-10); % meters^3 / Pascal

a = [-1*g*p/(Mp*Rm), 0, -Ap/Mp; 0, ((Cf + Cc)*Rt)^-1, -((Cf + Cc)*Rt)^-1; -Ap/Cf, -1/(Cf*Rt), 1/(Cf*Rt)];
b = [g/(Mp*Rm); 0; 0];

sys = ss(a,b,c,d);
[y,t,x] = lsim(sys,u,t_after,x0);

flow = (x(:,3)-x(:,2))./(Rt); 

plot(t(1:end), flow(1:end)*toLPM, '-r')


% figure;
% subplot(3,1,1), plot(t,y(:,3)*toPSI,t,y(:,2)*toPSI, t, u)
% axis([0 Tfinal -3 13])
% subplot(3,1,2),plot(t, flow*toLPM, '-r', t, zeros(length(t),1), '--k', MeasuredTime'-(981/1000), MeasuredFlow)
% axis([0 Tfinal -1 2])
% 
% for n = 1:length(t)
%     displacement(n) = trapz(y((1:n),1))*dt;
%     Volume(n) = trapz(flow(1:n))* 0.000166667;
% end
% 
% subplot(3,1,3), plot(t,displacement*1000)
% Volume(end)





