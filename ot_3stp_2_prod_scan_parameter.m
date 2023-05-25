%% finding ESSs and sweeping an extra parameter

%%
stpP = 10; %how many steps for parameter

P0 = 0; %parameter start value
PF = 0.1;  %parameter end value

stpG = 99;

stpK = 49;
kip0 = 3;
kipF = 0.05;


%% looping control
tspan = [0 tend];
stpszG = (GamIF-GamI0)/stpG;      %calc step sizes
stpszP = (PF-P0)/stpP;
stpszK = (kipF-kip0)/stpK;

GamIndex = GamI0:stpszG:GamIF; %record of all tested gamma values
PIndex = P0:stpszP:PF;
KipIndex = kip0:stpszK:kipF;

igi = 0; iP = 0; igr = 0; ik = 1;
Gam_finalabI = zeros(stpG+1, stpG+1);
Gam_finalabP = zeros(stpG+1, stpG+1);
Gam_MasterCell = cell(stpK+1, stpP+1);
Gam_MasterCellP = cell(stpK+1, stpP+1);

%% production scan

countdown = (stpK+1)*(stpP+1);
ik = 0;
for kip = kip0:stpszK:kipF
    ik = ik + 1;
    krp = kip;
    iP = 0;
   for P = P0:stpszP:PF   %over all possible parameter values
        iP = iP + 1;
        igi = 0;
        tic
            for GamI = GamI0:stpszG:GamIF   %production of IR toxin by inv
                igi = igi + 1;
                igri = 0;
                    for GamR = GamR0:stpszG:GamRF      %production of RI toxin by res
                        igri = igri + 1;
                            eventfunc = @(t,y) ot_3stp_NutrSteadyState(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE);
                            optionsode=odeset('Events',eventfunc,'NonNegative', 1:7);
                            y0 = [AbI AbR AbP AbTI AbTR AbTP Nu1];
                            [t,y,te,ye,ie] = ode45(@(t,y)ot_3stp_patch_mod(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE),tspan,y0,optionsode);
                            G = [t,y];
                            Gam_finalabI(igi,igri) = G(end,2); %generate a gamxgam matrix of final abI
                            Gam_finalabP(igi,igri) = G(end,4); %generate a gamxgam matrix of final abP
                        
                    end
            end    
            Gam_MasterCell{ik,iP} = Gam_finalabI;
            Gam_MasterCellP{ik,iP} = Gam_finalabP;
            toc 
            countdown = countdown - 1
   end
end
 


load splat
sound(y,Fs)




