
%% running comps for all gamma values at every spectrum

GamI0 = 0.0;              %production rate range
GamIF = 0.2;
GamR0 = GamI0;
GamRF = GamIF;
stpG = 50;


stpK = 49;
kip0 = 3;
kipF = 0.05;


%% looping control
tspan = [0 tend];
stpszG = (GamIF-GamI0)/stpG;      %calc step sizes
stpszK = (kipF-kip0)/stpK;

GamIndex = GamI0:stpszG:GamIF; %record of all tested values
KipIndex = kip0:stpszK:kipF;

igi = 0; igr = 0; ik = 1;
Gam_finalabI = zeros(stpG+1, stpG+1);   %all gamma comp results for one spectrum (invader abundance (strain A))
Gam_finalabP = zeros(stpG+1, stpG+1);   %same for pathogen (strain C) final abundance

Gam_MasterCellI = cell(stpK+1,1);       %cell with a results matrix for each spectrum
Gam_MasterCellP = cell(stpK+1,1);       %same but for pathogen abundance


%% spectrum
countdown = (stpK+1);
ik = 0;
for kip = kip0:stpszK:kipF
    ik = ik + 1;
    krp = kip;
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
                            Gam_finalabP(igi,igri) = G(end,4);
                        
                    end
            end    
            Gam_MasterCellI{ik,1} = Gam_finalabI;
            Gam_MasterCellP{ik,1} = Gam_finalabP;
            toc 
            countdown = countdown - 1
end

 


load splat
sound(y,Fs)




