%% Probio vs. pathogen

%Aiming to find the production rate which minimises pathogen abundance 
%if there are mutliple, aim to maximise probiotic abundance too

%range
% pathogen starting abundance
% pathogen toxin production
% pathogen toxicity

%sweep
% E
% Nutrient abundance


%% setting constants
clc;clear;

n = 2;
r = zeros(1,n+1);
R = 1;
r(:,:) = R;             %growth rate vector, all set to 1
tend = 100000; 
 %max time, ensures steady state reached

Nu1 = 1;                %nutrient starting abundnace    
              
kn1 = 5;                %bateria nutrient affinity            
nut = 1;                %pathogen nutrient binding modifier

km = 0.05;              %high affinity toxin binding

HCE = 1;
E = 5;
EP = 2;

GamI = 0.0;
GamR = 0.16;
GamP = 0.0;
AbI = 0.0001;
AbR = 0.0001;          %abundances at T0
AbP = AbI/1000;
AbTI = 0;
AbTR = 0;
AbTP = 0;
       
kir = 0.05;
kri = kir;
kip = 0.05;
krp = 300;  %conspecific toxoin is only narrow
kpi = 0.05;
kpr = 0.05;

Degr = 1;

%% looping

GamI0 = 0.0;
GamIF = 0.99;
stpG = 99;

GamR0 = 0.16;     %rwhether conspecific produces
GamRF = 0.0;
stpGamR = 1;

AbP0 = AbI/1000;
AbPF = AbI;
stpA = 1;

kip0 = 3;
kipF = 0.05;
stpKip = 24;

kir0 = 0.05;
kirF = 3;
stpKir = 24;

tspan = [0 tend];
stpszG = (GamIF-GamI0)/stpG;      %calc step sizes
AbIndex = [AbI/1000 AbI];
stpszKip = (kipF-kip0)/stpKip;
stpszKir = (kirF-kir0)/stpKir;
stpszGamR = (GamRF-GamR0)/stpGamR;

GamIndex = GamI0:stpszG:GamIF; %record of all tested gamma values
KipIndex = zeros(1,49); % allowing the toxin to target anyone
KirIndex = zeros(1,49);
KipIndex(1,:) = 0.05;
KirIndex(1,:) = 0.05;
KipIndex(1,1:25) = kip0:stpszKip:kipF;
KirIndex(1,25:49) = kir0:stpszKir:kirF;
GamRIndex = GamR0:stpszGamR:GamRF;

igi = 0; iA = 0; 
Gam_finalabsI = zeros((stpKip*2)+1, stpG+1);
Gam_finalabsP = zeros((stpKip*2)+1, stpG+1);
Gam_finalabsR = zeros((stpKip*2)+1, stpG+1);

%% sweeping different pathogen types

GamP0 = 0.0;
GamPF = 0.2;
stpGP = 1;

stpszGP = (GamPF-GamP0)/stpGP; 
GamPIndex = GamP0:stpszGP:GamPF;

Gam_MasterCellI = cell(4,stpGP+1);
Gam_MasterCellP = cell(4,stpGP+1);
Gam_MasterCellR = cell(4,stpGP+1);

%% compete
countdown = (stpA + 1)*(stpGP+1)*((stpKip*2)+1)*(stpGamR+1);
 
igamr = 0 ;
for GamR = GamR0:stpszGamR:GamRF
    igamr = igamr + 1;
    igp = 0;
    for GamP = GamP0:stpszGP:GamPF
        igp = igp + 1;
        iA = 0;
        for AbPStp = 1:1:stpA+1
            AbP = AbIndex(AbPStp);
                iA = iA + 1;
                ikip = 0;     
               for i = 1:1:(stpKip*2)+1
                   kip = KipIndex(i);
                   kir = KirIndex(i);
                    ikip = ikip+1;
                    igi = 0;
                    tic
                    for GamI = GamI0:stpszG:GamIF
                        igi = igi +1;
                            eventfunc = @(t,y) probio_NutrSteadyState(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE,nut);
                            optionsode=odeset('Events',eventfunc,'NonNegative', 1:7);
                            y0 = [AbI AbR AbP AbTI AbTR AbTP Nu1];
                            [t,y,te,ye,ie] = ode45(@(t,y)probio_patch_mod(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE,nut),tspan,y0,optionsode);
                            Vals = [t,y];
                            if Vals(end,2) > 1E-09 
                               Gam_finalabsI(ikip,igi) = Vals(end,2);
                            else 
                               Gam_finalabsI(ikip,igi) = 0;  
                            end
                           if Vals (end,4) > 1E-09
                               Gam_finalabsP(ikip,igi) = Vals(end,4);
                           else 
                               Gam_finalabsP(ikip,igi) = 0; 
                           end
                           if Vals (end,3) > 1E-09
                               Gam_finalabsR(ikip,igi) = Vals(end,3);
                           else 
                               Gam_finalabsR(ikip,igi) = 0; 
                           end
                    end  
            toc 
            countdown = countdown - 1
                end
           if igamr == 1
            Gam_MasterCellI{iA,igp} = Gam_finalabsI;
            Gam_MasterCellP{iA,igp} = Gam_finalabsP;
            Gam_MasterCellR{iA,igp} = Gam_finalabsR;
           elseif igamr == 2
            Gam_MasterCellI{iA + 2,igp} = Gam_finalabsI;
            Gam_MasterCellP{iA + 2,igp} = Gam_finalabsP;
            Gam_MasterCellR{iA + 2,igp} = Gam_finalabsR;
           end
        end
        
    end
end


%% finding minimum P val and its position

for krpcount = 1:1:stpGamR+1
    for abcount = 1:1:stpA+1
        for gampcount = 1:1:stpGP+1
            if krpcount == 1
                Gam_finalabsI = Gam_MasterCellI{abcount,gampcount};
                Gam_finalabsP = Gam_MasterCellP{abcount,gampcount};
            elseif krpcount == 2
                Gam_finalabsI = Gam_MasterCellI{abcount + 2,gampcount};
                Gam_finalabsP = Gam_MasterCellP{abcount + 2,gampcount};
            end
            if krpcount == 1
               Gam_minP(abcount,gampcount) = min(Gam_finalabsP,[],'all');
            elseif krpcount == 2
              Gam_minP(abcount + 2,gampcount) = min(Gam_finalabsP,[],'all');
            end
        end
    end
end 
%% finding values when P in minimal 
%if multiple min P, find the one where probiotic abundance is highest

bestAbsI = zeros(stpA+2,stpGP+1);
bestAbsP = zeros(stpA+2,stpGP+1);
idealGam = zeros(stpA+2,stpGP+1); 
idealKip = zeros(stpA+2,stpGP+1);




for krpcount = 1:1:stpGamR+1
    for abcount = 1:1:stpA+1
        for gampcount = 1:1:stpGP+1
            if krpcount == 1
                Gam_finalabsI = Gam_MasterCellI{abcount,gampcount};
                Gam_finalabsP = Gam_MasterCellP{abcount,gampcount};
            elseif krpcount == 2
                Gam_finalabsI = Gam_MasterCellI{abcount + 2,gampcount};
                Gam_finalabsP = Gam_MasterCellP{abcount + 2,gampcount};
            end
            for kipcount = 1:1:stpKip+1
                for gamcount = 1:1:stpG+1
                    pos = find(Gam_finalabsP == min(Gam_finalabsP,[],"all"));
                    amount = size(pos);
                    if amount(1,1) > 1 
                        [minprows,minpcolls] = find(Gam_finalabsP == 0);
                        Poss_final_Is = 0;
                        for trying = 1:1:size(minpcolls)
                            Poss_final_Is(1,trying) = Gam_finalabsI(minprows(trying), minpcolls(trying));
                        end
                          best_final_I = max(Poss_final_Is);
                          [bestrow,bestcol] = find(Gam_finalabsI == best_final_I);
                    else
                        [bestrow,bestcol] = find(Gam_finalabsP == min(Gam_finalabsP,[],"all"));
                    end
                    if krpcount == 1
                        bestAbsI(abcount,gampcount) = Gam_MasterCellI{abcount,gampcount}(bestrow(1,1),bestcol(1,1));
                        bestAbsP(abcount,gampcount) = Gam_MasterCellP{abcount,gampcount}(bestrow(1,1),bestcol(1,1));
                        idealGam(abcount,gampcount) = GamIndex(bestcol(1,1));
                        idealKip(abcount,gampcount) = KipIndex(bestrow(1,1));
                        idealKir(abcount,gampcount) = KirIndex(bestrow(1,1));
                    elseif krpcount == 2
                         bestAbsI(abcount + 2,gampcount) = Gam_MasterCellI{abcount + 2,gampcount}(bestrow(1,1),bestcol(1,1));
                        bestAbsP(abcount + 2,gampcount) = Gam_MasterCellP{abcount + 2,gampcount}(bestrow(1,1),bestcol(1,1));
                        idealGam(abcount + 2,gampcount) = GamIndex(bestcol(1,1));
                        idealKip(abcount + 2,gampcount) = KipIndex(bestrow(1,1));
                        idealKir(abcount + 2,gampcount) = KirIndex(bestrow(1,1));
                        
                    end
                end
                
            end
        end
    end
end

%% heatmaps of Path removal

%bestAbsP is a matrix showing the lowest final Path ratio for each condition
%                            Passive path   |  toxic path
% Low AbPath + narrow con
% Low AbPath + broad con
% mid AbPath + narrow con
% mid AbPath + broad con


    
    V1a = 1/255*[ 255, 255, 255;  242, 126, 126];
    
    X1a = [0 1];
    
    Xqa = linspace(0, 1, 255);
    Vq1a = interp1(X1a, V1a, Xqa, 'pchip');
    
    upper =max(bestAbsP,[],'all');

%heatmap of ratios of a low abundance pathogen 
figure(1)
    lowAbAbs = [bestAbsP(1,1) bestAbsP(1,2); bestAbsP(3,1) bestAbsP(3,2)];
    lowAbKips = [idealKip(1,1) idealKip(1,2); idealKip(3,1) idealKip(3,2)];
       lowAbKirs = [idealKir(1,1) idealKir(1,2); idealKir(3,1) idealKir(3,2)];
    lowAbGams = [idealGam(1,1) idealGam(1,2); idealGam(3,1) idealGam(3,2)];
heatmap(lowAbAbs)
colormap(Vq1a)
clim([0 upper])

figure(2)
    midAbAbs = [bestAbsP(2,1) bestAbsP(2,2); bestAbsP(4,1) bestAbsP(4,2)];
     midAbKips = [idealKip(2,1) idealKip(2,2); idealKip(4,1) idealKip(4,2)];
      lowAbKirs = [idealKir(2,1) idealKir(2,2); idealKir(4,1) idealKir(4,2)];
      midAbGams = [idealGam(2,1) idealGam(2,2); idealGam(4,1) idealGam(4,2)];
heatmap(midAbAbs)
colormap(Vq1a)
clim([0 upper])



