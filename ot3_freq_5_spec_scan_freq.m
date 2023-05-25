
Kip_finalabI = zeros(stpK+1, stpK+1);
Kip_MasterCellI = cell(1, stpF+1);
Kip_finalabR = zeros(stpK+1, stpK+1);
Kip_MasterCellR = cell(1, stpF+1);
Kip_finalabP = zeros(stpK+1, stpK+1);
Kip_MasterCellP = cell(1, stpF+1);

Kip_path_finalabI = zeros(stpK+1,stpK+1); %contest with no pathogen
Kip_path_finalabR = zeros(stpK+1,stpK+1);
Kip_path_finalabP = zeros(stpK+1,stpK+1);
Kip_none_finalabI = zeros(stpK+1,stpK+1); %contest with no pathogen
Kip_none_finalabR = zeros(stpK+1,stpK+1);
Kip_none_finalabP = zeros(stpK+1,stpK+1);
Kip_none_MasterCellI = cell(1, stpF+1);
Kip_path_MasterCellI = cell(1, stpF+1);

krp0 = kip0;
krpF = kipF;

%%
countdown = (stpF + 1);
tic
iF = 0;

for freq = Freq0:stpszF:FreqF   %over all possible GamIP (= GamIP)
        iF = iF + 1;
        ikrp = 0;
        tic
        for krp = krp0:stpszK:kipF      %vary resident spectrum
        ikrp = ikrp + 1;
        iikip = 0;
            for kip = kip0:stpszK:kipF  %vary invader spectrum
                iikip = iikip + 1;
                if invalid_G_ESS(ikrp,iF) == 1      %not checked for currently
                    Kip_finalabI(iikip,ikrp) = 1; 
                    Kip_finalabR(iikip,ikrp) = 1; 
                    Kip_finalabP(iikip,ikrp) = 1; 
                else
                GamI = new_Gam_ESS(ikrp,iF);
                GamR = new_Gam_ESS(ikrp,iF);
                           AbP = 0;
                            eventfunc = @(t,y) ot3_freq_NutrSteadyState(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE);
                            optionsode=odeset('Events',eventfunc,'NonNegative', 1:7);
                            y0 = [AbI AbR AbP AbTI AbTR AbTP Nu1];
                            [t,y,te,ye,ie] = ode45(@(t,y)ot3_freq_patch_mod(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE),tspan,y0,optionsode);
                            K_none = [t,y];
                            Kip_none_finalabI(iikip,ikrp) = K_none(end,2); %contest with no pathogen
                            Kip_none_finalabR(iikip,ikrp)= K_none(end,3);
                            Kip_none_finalabP(iikip,ikrp) = K_none(end,4);

                            AbP = AbI;
                            eventfunc = @(t,y) ot3_freq_NutrSteadyState(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE);
                            optionsode=odeset('Events',eventfunc,'NonNegative', 1:7);
                            y0 = [AbI AbR AbP AbTI AbTR AbTP Nu1];
                            [t,y,te,ye,ie] = ode45(@(t,y)ot3_freq_patch_mod(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE),tspan,y0,optionsode);
                            K_path = [t,y];
                            Kip_path_finalabI(iikip,ikrp) = K_path(end,2); %contest with no pathogen
                            Kip_path_finalabR(iikip,ikrp) = K_path(end,3);
                            Kip_path_finalabP(iikip,ikrp) = K_path(end,4);

                           if Kip_path_finalabI(iikip,ikrp) < 1E-9
                                 Kip_path_finalabI(iikip,ikrp) = 0;
                           end
                           if Kip_none_finalabI(iikip,ikrp) < 1E-9
                             Kip_none_finalabI(iikip,ikrp) = 0;
                           end
                             Kip_none_MasterCellI{1,iF} = Kip_none_finalabI;
                             Kip_path_MasterCellI{1,iF} = Kip_path_finalabI;
                             %generally keep a record of what happens here,
                             %helps with troubleshooting
                                       
                end %then do averging
                            Kip_finalabI(iikip,ikrp) = (1-freq)*Kip_none_finalabI(iikip,ikrp) + (freq)*Kip_path_finalabI(iikip,ikrp);
                            Kip_finalabR(iikip,ikrp) = (1-freq)*Kip_none_finalabR(iikip,ikrp) + (freq)*Kip_path_finalabR(iikip,ikrp);
                            Kip_finalabP(iikip,ikrp) = (1-freq)*Kip_none_finalabP(iikip,ikrp) + (freq)*Kip_path_finalabP(iikip,ikrp);

                
            end
    end
    Kip_MasterCellI{1,iF} = Kip_finalabI;
    Kip_MasterCellR{1,iF} = Kip_finalabR;
    Kip_MasterCellP{1,iF} = Kip_finalabP;
    toc
    countdown = countdown - 1
end

toc
load splat
sound(y,Fs)