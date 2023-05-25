%% competeing all spectrum values using calculated gamma ESS

Kip_finalabI = zeros(stpKip+1, stpKip+1);
Kip_MasterCellI = cell(1, stpP+1);
Kip_finalabR = zeros(stpKip+1, stpKip+1);
Kip_MasterCellR = cell(1, stpP+1);
Kip_finalabP = zeros(stpKip+1, stpKip+1);
Kip_MasterCellP = cell(1, stpP+1);

krp0 = kip0;
krpF = kipF;

%%
countdown = stpP + 1;
tic
iP = 0;

for P = P0:stpszP:PF   %over all possible GamIP (= GamIP)
    iP = iP + 1;
    ikrp = 0;
    tic
        for krp = krp0:stpszK:kipF      %vary resident spectrum
        ikrp = ikrp + 1;
        iikip = 0;
            for kip = kip0:stpszK:kipF  %vary invader spectrum
                iikip = iikip + 1;
                if invalid_G_ESS(ikrp,iP) == 1
                    Kip_finalabI(iikip,ikrp) = 1; 
                    Kip_finalabR(iikip,ikrp) = 1; 
                    Kip_finalabP(iikip,ikrp) = 1; 
                else
                GamI = Gam_ESS(ikrp,iP);
                GamR = Gam_ESS(ikrp,iP);
                    eventfunc = @(t,y) ot_3stp_NutrSteadyState(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE);
                    optionsode=odeset('Events',eventfunc,'NonNegative', 1:7);
                    y0 = [AbI AbR AbP AbTI AbTR AbTP Nu1];
                    [t,y,te,ye,ie] = ode45(@(t,y)ot_3stp_patch_mod(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE),tspan,y0,optionsode);
                    K = [t,y];
                    Kip_finalabI(iikip,ikrp) = K(end,2); %generate a matrix of kipxkip of final abI
                    Kip_finalabR(iikip,ikrp) = K(end,3); %matrix of resident final abundances
                    Kip_finalabP(iikip,ikrp) = K(end,4); %matrix of pathogen final abundances
                end
            end
    end
    Kip_MasterCellI{1,iP} = Kip_finalabI;
    Kip_MasterCellR{1,iP} = Kip_finalabR;
    Kip_MasterCellP{1,iP} = Kip_finalabP;
    toc
    countdown = countdown - 1
end

toc
load splat
sound(y,Fs)