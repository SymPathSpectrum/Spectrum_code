%Coevo - finding GamESS for all spectrum values
clc;clear;
%% setting constants
n = 2;
r = zeros(1,n+1);
R = 1;
r(:,:) = R;             %growth rate vector, all set to 1
r(1,3) = 1;             %for altering path growth rate
tend = 100000; 
 %max time, ensures steady state reached

Nu1 = 1;                %nutrient starting abundnace    
              
kn1 = 5;                %bacteria nutrient affinity            


km = 0.05;              %high affinity toxin binding
     
HCE = 1;        %hill coeff
E = 5;          %toxin efficiencies
EP = 2;


%% other values
         %Gam_ip and Gam_ir prodction steps


GamP = 0.02;     %pathogen production rate

AbI = 0.0001;
AbR = 0.0001;          %abundances at T0
AbP = 0.000001;
AbTI = 0;
AbTR = 0;
AbTP = 0;
       
kir = 0.05;         %affinities (will be ranged later)
kri = kir;
kip = 0.05;
krp = 0.05;
kpi = 0.05;
kpr = 0.05;

Degr = 1;




