%Coevo - finding GamESS for all spectrum values
clc;clear;
%% setting constants
n = 2;
r = zeros(1,n+1);
R = 1;
r(:,:) = R;             %growth rate vector, all set to 1
tend = 100000; 
 %max time, ensures steady state reached

Nu1 = 1;                %nutrient starting abundnace    
              
kn1 = 5;                %bateria nutrient affinity            


km = 0.05;              %high affinity toxin binding


HCE = 1;
E = 5;
EP = 2;


%% Looping values
         %Gam_ip and Gam_ir prodction steps

GamI0 = 0;              %production rate range
GamIF = 0.2;
GamR0 = GamI0;
GamRF = GamIF;
GamP = 0.2;

AbI = 0.0001;
AbR = 0.0001;          %abundances at T0
AbP = AbI;
AbTI = 0;
AbTR = 0;
AbTP = 0;
       
kir = 0.05;
kri = kir;
kip = 0.05;
krp = 0.05;
kpi = 0.05;
kpr = 0.05;

Degr = 1;


