function reactor1(g)
% REACTOR1 Zero-dimensional kinetics: adiabatic, constant pressure.
%
%
%    This example illustrates how to use class 'Reactor' for
%    zero-dimensional kinetics simulations. Here the parameters are
%    set so that the reactor is adiabatic and very close to constant
%    pressure.
%

help reactor1

P = 101325; %in Pascal
T = 800; %in Kelvin
dt = 1e-4; %time step
TotalTime = 5; % in seconds - includes autoignition phase

nSteps = ceil(TotalTime/dt); %number of steps. Total time = nSteps*dt


if strcmp(g,'kerosene')
   gas = Solution('kerosene.yaml', 'gas');
else
   gas = GRI30('None');
end

% set the initial conditions. All values below are for stoichiometric
% conditions. Divide the values of O2 and N2 by the equivalence ratio to
% get the lean condition.
% set(gas,'T',T,'P',P,'X','H2:1,O2:0.5,N2:1.88'); %H2
%set(gas,'T',T,'P',P,'X','CH4:1,O2:2,N2:7.52');  % methane
%set(gas,'T',T,'P',P,'X','NC10H22:1,O2:15.5,N2:58.28') % decane
%set(gas,'T',T,'P',P,'X','NC10H22:0.74,PHC3H7:0.15,CYC9H18:0.11,O2:14.76,N2:55.45') %kerosene
set(gas,'T',T,'P',P,'X','NC10H22:0.74,PHC3H7:0.15,CYC9H18:0.11,H2:1,O2:15.26,N2:57.38') % 50% H2 in volume  
%set(gas,'T',T,'P',P,'X','NC10H22:0.74,PHC3H7:0.15,CYC9H18:0.11,H2:60,O2:44.76,N2:168.3') % 50% H2 in volume  


% create a reactor, and insert the gas
r = IdealGasReactor(gas);

% create a reservoir to represent the environment
a = Solution('air.yaml','air','None');
set(a,'P',P)
env = Reservoir(a);

% Define a wall between the reactor and the environment and
% make it flexible, so that the pressure in the reactor is held
% at the environment pressure.
w = Wall;
install(w,r,env);

% set expansion parameter. dV/dt = KA(P_1 - P_2)
setExpansionRateCoeff(w, 1.0e6);

% set wall area
setArea(w, 1.0);

% create a reactor network and insert the reactor:
network = ReactorNet({r});




tim(nSteps) = 0;
temp(nSteps) = 0;
x(nSteps,7) = 0;
kero(nSteps,3) = 0;
t = 0.0;
t0 = cputime;
for n = 1:nSteps
  t = t + dt;
  advance(network, t);
  tim(n) = time(network);
  temp(n) = temperature(r);
  x(n,1:7) = massFraction(gas,{'CH4','CO','CO2','H2O','NO','NO2','H2'});
  kero(n,1:3) = massFraction(gas,{'NC10H22','PHC3H7','CYC9H18'});
end
disp(['CPU time = ' num2str(cputime - t0)]);

clf;
subplot(2,2,1);
plot(tim,temp,'r','LineWidth',2);
xlabel('Time (s)');
ylabel('Temperature (K)');
subplot(2,2,2)
plot(tim,sum(kero,2),'k',tim,x(:,3),'r',tim,x(:,4),'b',tim,x(:,end),'g');
xlabel('Time (s)');
ylabel('Mass fraction');
legend('Kerosene','CO2','H2O','H2');
subplot(2,2,3)
plot(tim,x(:,2));
xlabel('Time (s)');
ylabel('CO Mass Fraction');
subplot(2,2,4)
plot(tim,(x(:,5)+x(:,6))*1e6);
xlabel('Time (s)');
ylabel('NOX Mass Fraction (ppm)');
clear all
cleanup
