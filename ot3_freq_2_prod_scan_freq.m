% %Coevo - finding GamESS for all spectrum values

%%
 stpF = 49;     %how many frequency steps (slows down spec_scan later!)
 
 Freq0 = 0.0;
 FreqF = 1.0;  

stpG = 99;

stpK = 1;
kip0 = 3;
kipF = 0.05;


%% looping control
tspan = [0 tend];
stpszG = (GamIF-GamI0)/stpG;      %calc step sizes
stpszF = (FreqF-Freq0)/stpF;
stpszK = (kipF-kip0)/stpK;

GamIndex = GamI0:stpszG:GamIF; %record of all tested gamma values
FreqIndex = Freq0:stpszF:FreqF;
KipIndex = kip0:stpszK:kipF;

igi = 0; iF = 0; igr = 0; ik = 1;
Gam_finalabI = zeros(stpG+1, stpG+1);
G_none_finalabI = zeros(stpG+1, stpG+1);
G_path_finalabI = zeros(stpG+1, stpG+1);

Gam_MasterCell = cell(stpK+1, stpF+1);
G_none_MasterCell = cell(stpK+1, stpF+1);
G_path_MasterCell = cell(stpK+1, stpF+1);

%% spectrum
countdown = (stpK+1);
ik = 0;
for kip = kip0:stpszK:kipF
    ik = ik + 1;
    krp = kip;
    iF = 0;
    tic
            for GamI = GamI0:stpszG:GamIF   %production of IR toxin by inv
                igi = igi + 1;
                igri = 0;
                    for GamR = GamR0:stpszG:GamRF      %production of RI toxin by res
                        igri = igri + 1;
                            AbP = 0;
                            eventfunc = @(t,y) ot3_freq_NutrSteadyState(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE);
                            optionsode=odeset('Events',eventfunc,'NonNegative', 1:7);
                            y0 = [AbI AbR AbP AbTI AbTR AbTP Nu1];
                            [t,y,te,ye,ie] = ode45(@(t,y)ot3_freq_patch_mod(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE),tspan,y0,optionsode);
                            G_none = [t,y];
                             if G_none(end,2) > 1E-9
                                G_none_finalabI(igi,igri) = G_none(end,2); %contest with no pathogen
                            else
                                G_none_finalabI(igi,igri) = 0; %zeros values that are too small
                            end  

                            AbP = AbI; %contest with pathogen       
                            eventfunc = @(t,y) ot3_freq_NutrSteadyState(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE);
                            optionsode=odeset('Events',eventfunc,'NonNegative', 1:7);
                            y0 = [AbI AbR AbP AbTI AbTR AbTP Nu1];
                            [t,y,te,ye,ie] = ode45(@(t,y)ot3_freq_patch_mod(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE),tspan,y0,optionsode);
                            G_path = [t,y];
                            if G_path(end,2) > 1E-9
                                G_path_finalabI(igi,igri) = G_path(end,2); 
                            else
                                G_path_finalabI(igi,igri) = 0;
                            end  
                    end
            end    
            for freq = Freq0:stpszF:FreqF   %now calculate average based on presence and absence
                iF = iF + 1;
                igi = 0;
                Gam_finalabI = (1-freq)*G_none_finalabI + (freq)*G_path_finalabI;
                Gam_MasterCell{ik,iF} = Gam_finalabI;
                G_none_MasterCell{ik,iF} = G_none_finalabI;
                G_path_MasterCell{ik,iF} = G_path_finalabI;
            end    
            toc 
            countdown = countdown - 1

end
 
load splat
sound(y,Fs)




