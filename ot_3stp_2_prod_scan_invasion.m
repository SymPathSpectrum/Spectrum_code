%% finding invasion outcomes for all production rates
%% test values
 stpG = 100;

 GamI0 = 0;
 GamIF = 1;  

 GamR0 = GamI0;
 GamRF = GamIF;
 GamP = 0;

AbP = AbI/1000;

kip = 3;        %spectrum kept constant
krp = 3;

%% looping control
tspan = [0 tend];
    
%calc step sizes
stpszG = (GamIF-GamI0)/stpG;


%record of all tested gamma values
GamIndex = GamI0:stpszG:GamIF;

iG = 0; 
Gam_finalabI = zeros(1, stpG+1);
Gam_finalabR = zeros(1, stpG+1);
Gam_finalabP = zeros(1, stpG+1);




%Ab_Pips = [AbIndex(3),AbIndex(25),AbIndex(48)];

%% spectrum

    iG = 0;
   for GamI = GamI0:stpszG:GamIF  %over all possible GamIP (= GamIP)
        iG = iG + 1;
        GamR = GamI;

            eventfunc = @(t,y) ot_3stp_NutrSteadyState(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE);
            optionsode=odeset('Events',eventfunc,'NonNegative', 1:7);
            y0 = [AbI AbR AbP AbTI AbTR AbTP Nu1];
            [t,y,te,ye,ie] = ode45(@(t,y)ot_3stp_patch_mod(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE),tspan,y0,optionsode);
            G = [t,y];
            Gam_finalabI(1,iG) = G(end,2); %generate a vector of final abundances
            Gam_finalabR(1,iG) = G(end,3);
            Gam_finalabP(1,iG) = G(end,4);
         
    end    

%% find the relative final abundance (normalised to 1)

maxes = [max(Gam_finalabP), max(Gam_finalabR), max(Gam_finalabI)];
the_max = max(maxes);

RFA_I = Gam_finalabI./the_max;
RFA_R = Gam_finalabR./the_max;
RFA_P = Gam_finalabP./the_max;
% RFA_F = Gam_finalabF./the_max;



%%
figure(1)
hold on
%title("Impact of toxin production rate on invasion by a passive invader")

%plot(G_ESS_line, '.');


plot((RFA_I(1,:)+0.01), 'Color', 1/255*[204 105 102], 'lineWidth', 6)
plot((RFA_R(1,:)), 'Color', 1/255*[100 158 185], 'lineWidth', 4)
plot((RFA_P(1,:)), 'Color', 1/255*[232 191 59], 'lineWidth', 6)
ylabel("Relative Final Abundance")
xlabel("Broad Toxin Production Rate")
    xticks([0 11 21 31 41 51 61 71 81 91 101])
    xticklabels({'0','',GamIndex(21),'',GamIndex(41),'',GamIndex(61),'',GamIndex(81),'',GamIndex(101)})
    yticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
    yticklabels({'0', '', '0.2', '', '0.4', '', '0.6', '', '0.8', '', '1'})
    xlim([1 stpG+1])
 %      ylim([-0.005 1.05])
       set(gca,'Layer','top')
 legend("Resident strain A", "Resident strain B", "Invading strain C", 'location', 'east')
    hold off



%%
load splat
sound(y,Fs)




