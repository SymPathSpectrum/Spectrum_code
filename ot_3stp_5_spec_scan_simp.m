%% competitions for alll spectrum values at the calculated ESS

Kip_finalabI = zeros(stpK+1, stpK+1);
Kip_finalabR = zeros(stpK+1, stpK+1);
Kip_finalabP = zeros(stpK+1, stpK+1);

krp0 = kip0;
krpF = kipF;


%%

ikrp = 0;
    tic
    for krp = krp0:stpszK:krpF      %vary resident spectrum
    ikrp = ikrp + 1;
    iikip = 0;
         for kip = kip0:stpszK:kipF  %vary invader spectrum
             iikip = iikip + 1;
            GamI = Gam_ESS(1, ikrp);
            GamR = Gam_ESS(1, ikrp);
           if invalid_G_ESS(1,ikrp) == 1        %if GESS invalid, PIP is white for this spectrum value
                Kip_finalabI(iikip,ikrp) = 1;
           else
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


toc
load splat
sound(y,Fs)