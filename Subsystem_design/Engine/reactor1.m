function [TPZ, MF_emis, name_emis,COf, NOxf] = reactor1(g, P_input, T_input, eqr_input,n_h2,n_ker,n_O2,n_N2)

%    REACTOR1 Zero-dimensional kinetics: adiabatic, constant pressure.
%
%    This example illustrates how to use class 'Reactor' for
%    zero-dimensional kinetics simulations. Here the parameters are
%    set so that the reactor is adiabatic and very close to constant
%    pressure.

    function [time_idx1, time_idx2] = time_res(t,dt,p)

            slope = (t(2:end) - t(1:end-1))/dt;
            perc = p*max(slope); %5 percent time
            %result = find(slope == 0.7*max(slope));

            i = 1;
            for c = 1:(length(slope)-1)
                % disp(slope(c))
                if (perc >= slope(c)) && (perc <= slope(c+1))
                    i = c+1;
                    break;
                end
            end

            j = 1;
            for c = 1:(length(t)-1)
                if t(c) >= 0.99*t(end)
                    j = c;
                    break;
                end
            end

            time_idx1 = i;
            time_idx2 = j;
    end
    
    function [rest, tim1, tim2] = residence_time(T,time,rlimit)
        rlim1 = 0.08;
        rlim2 = rlimit;
        T1 = T(1)+rlim1*(T(end)-T(1));
        idx1 = find(T >= T1);
        tim1 = time(idx1(1));
        T2 = T(1)+rlim2*(T(end)-T(1));
        idx2 =find(T >= T2);
        tim2 = time(idx2(1));
        rest = 1000*(tim2 - tim1);
    end

% help reactor1

%-----------------
% Phases - based on data from GSP 11

% Cruise:
% P = 11.70102*100000; %101325; %in Pascal
% T = 724.90893; %800; in Kelvin

% (Take-Off:)
% Inputs from Python
P = P_input; %101325; %in Pascal
T = T_input; %800; in Kelvin
eqr = eqr_input; %0.3;

%------------------

<<<<<<<<< Temporary merge branch 1
dt = 1e-1; %time step; -4 originally
TotalTime = 20; % in seconds - includes autoignition phase

=========
dt = 1e-6; %time step; -4 originally
TotalTime = 0.8; % in seconds - includes autoignition phase
>>>>>>>>> Temporary merge branch 2
nSteps = ceil(TotalTime/dt); %number of steps. Total time = nSteps*dt
rlim = 0.995;

if strcmp(g,'neo') %   compare string
   gas = Solution('kerosene.yaml', 'gas');
   p_o2 = n_O2; %14.76;
   p_n2 = n_N2; %55.45;
   p_ker = n_ker;
   p_o2_new = p_o2/eqr;
   p_n2_new = p_n2/eqr;

   str_ker = convertStringsToChars(join(['NC10H22:',string(0.74*p_ker),',PHC3H7:',string(0.15*p_ker),',CYC9H18:',string(0.11*p_ker),',O2:',string(p_o2_new),',N2:',string(p_n2_new)],"")); % only kerosene
   set(gas,'T',T,'P',P,'X',str_ker) % only kerosene
   %gas = Solution('nDodecane_Reitz.yaml','nDodecane_IG');

elseif strcmp(g,'hack_mix')
   gas = Solution('kerosene.yaml', 'gas');
   p_o2 = n_O2; %15.26;%44.76;
   p_n2 = n_N2; %57.38;%168.3;
   p_ker = n_ker;
   p_h2 = n_h2;
   p_o2_new = p_o2/eqr;
   p_n2_new = p_n2/eqr;
   %str_kerosene = convertStringsToChars(join(['NC10H22:0.74,PHC3H7:0.15,CYC9H18:0.11,O2:',string(p_o2_new),',N2:',string(p_n2_new)],""))
   str_ker_h2 = convertStringsToChars(join(['NC10H22:',string(0.74 * p_ker),',PHC3H7:',string(0.15*p_ker),',CYC9H18:',string(0.11*p_ker),',H2:',string(p_h2),',O2:',string(p_o2_new),',N2:',string(p_n2_new)],"")); % kerosene and H2, 50% in volume
   set(gas,'T',T,'P',P,'X',str_ker_h2) % 50% H2 in volume - change this (automatically based on mass ratio)

elseif strcmp(g,'hack_h2')
   gas = GRI30('None');
   p_o2 = n_O2; %0.5;
   p_n2 = n_N2; %1.88;
   p_h2 = n_h2;
   p_o2_new = p_o2/eqr;
   p_n2_new = p_n2/eqr;
   str_h2 = convertStringsToChars(join(['H2:',string(p_h2),',O2:',string(p_o2_new),',N2:',string(p_n2_new)],""));
   set(gas,'T',T,'P',P,'X',str_h2); %H2
end

% set the initial conditions. All values below are for stoichiometric
% conditions. Divide the values of O2 and N2 by the equivalence ratio to
% get the lean condition.

%set(gas,'T',T,'P',P,'X','H2:1,O2:0.5,N2:1.88'); %H2
%set(gas,'T',T,'P',P,'X','CH4:1,O2:2,N2:7.52');  % methane
%set(gas,'T',T,'P',P,'X','NC10H22:1,O2:15.5,N2:58.28') % decane

% Case: Kerosene Only

%p_o2 = 14.76;
%p_n2 = 55.45;
%eqr = eqr_input; %0.3;

%p_o2_new = p_o2/eqr;
%p_n2_new = p_n2/eqr;

%str_kerosene = convertStringsToChars(join(['NC10H22:0.74,PHC3H7:0.15,CYC9H18:0.11,O2:',string(p_o2_new),',N2:',string(p_n2_new)],""));
%set(gas,'T',T,'P',P,'X',str_kerosene) %kerosene

%set(gas,'T',T,'P',P,'X','NC10H22:0.74,PHC3H7:0.15,CYC9H18:0.11,O2:14.76,N2:55.45')% kerosene, stoich

%set(gas,'T',T,'P',P,'X','NC10H22:0.74,PHC3H7:0.15,CYC9H18:0.11,H2:1,O2:15.26,N2:57.38') % 50% H2 in volume
%set(gas,'T',T,'P',P,'X','NC10H22:0.74,PHC3H7:0.15,CYC9H18:0.11,H2:60,O2:44.76,N2:168.3') % 50% H2 in volume
%set(gas,'T',T,'P',P,'X','c12h26:1,O2:18.5,N2:69.56'); %For soot


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

% Is it here that I can change the reaction time?
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

name_emis = ['CH4, ','CO, ','CO2',', H2O',', NO',', NO2',', H2']
MF_emis =  x(nSteps,1:7)
disp(['CPU time = ' num2str(cputime - t0)]);

%p = 0.005;
%[t_five_i, t_com_i] = time_res(temp,dt,p);
%t_five = tim(t_five_i) + dt/2;
%t_com = tim(t_com_i);

[rt, tim1, tim2] = residence_time(temp,tim,rlim);
res_time = rt

%t_res = (t_com - t_five)*1000;
% disp(['Residence time = ', t_res, ' ms']);

T_tot = temp(length(temp));

idx = find(temp >= rlim*max(temp));
i = idx(1);

clf; %  clear figure
subplot(2,2,1);
hold on;
plot(tim,temp,'r','LineWidth',2)%,
plot([tim1 tim1],[0 3000],'k-');
plot([tim2 tim2],[0 3000],'k-');
hold off;
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
hold on;
plot(tim,(x(:,5)+x(:,6))*1e6);
plot([tim(i) tim(i)],[0 max((x(:,5)+x(:,6))*1e6)],'k-');
hold off;
xlabel('Time (s)');
ylabel('NOX Mass Fraction (ppm)');

COf = (x(end-1,2));
NOxf =(x(i,5)+x(i,6))*1e6;
% disp(['CO fraction = ', x(nSteps,2)]);
% disp(['NOx fraction = ', (x(nSteps,5)+x(nSteps,6))*1e6]);

TPZ = temp(length(temp))

% clear all
% cleanup
% Add a calculation of 5% steep angle
end

