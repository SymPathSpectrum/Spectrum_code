%% testing all spectrum values at each freq using gamESS

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
                if invalid_G_ESS(ikrp,iF) == 1
                    Kip_finalabI(iikip,ikrp) = 1; 
                    Kip_finalabR(iikip,ikrp) = 1; 
                    Kip_finalabP(iikip,ikrp) = 1; 
                else
                GamI = new_Gam_ESS(ikrp,iF);
                GamR = new_Gam_ESS(ikrp,iF);
                           AbP = 0;
                            eventfunc = @(t,y)abc_NutrSteadyState(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE);
                            optionsode=odeset('Events',eventfunc,'NonNegative', 1:7);
                            y0 = [AbI AbR AbP AbTI AbTR AbTP Nu1];
                            [t,y,te,ye,ie] = ode45(@(t,y)abc_patch_mod(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE),tspan,y0,optionsode);
                            K_none = [t,y];
                            Kip_none_finalabI(iikip,ikrp) = K_none(end,2); %contest with no pathogen
                            Kip_none_finalabR(iikip,ikrp)= K_none(end,3);
                            Kip_none_finalabP(iikip,ikrp) = K_none(end,4);

                            AbP = AbI;
                            eventfunc = @(t,y)abc_NutrSteadyState(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE);
                            optionsode=odeset('Events',eventfunc,'NonNegative', 1:7);
                            y0 = [AbI AbR AbP AbTI AbTR AbTP Nu1];
                            [t,y,te,ye,ie] = ode45(@(t,y)abc_patch_mod(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE),tspan,y0,optionsode);
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
                                       
                end
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

%% find kipESS make sure to see what production rates look like first!

%currently only using usual method (lowest invasion number) as issues with
%mutual invasibility not seen

xi2 = 0;
Kip_ESS = zeros(1,stpF+1); %matrix, Kess for each freq
total = 0;
toad = 0;
k_value = 0;
double_count = zeros(1, stpK+1);

%%
c = 0; c2 = 0; c3 = 0; c4 = 0;
for c = 1:1:stpF+1
    for c6 = 1:1:stpK+1
        for c7 = 1:1:stpK+1
           if Kip_MasterCellI{1,c}(c6,c7) < (1E-9)*(FreqIndex(c))
             Kip_MasterCellI{1,c}(c6,c7) = 0;
           end
        end
    end
end

ai = 0; a = 0; k = 0; m = 0; M = 0;
for a = 1:1:stpF+1   %for each freq
    ai = ai + 1;
    Kip_finalabI = Kip_MasterCellI{1,a};
        for k = 1:1:stpK+1      %for every tested kip value
        if Kip_finalabI(k,k) == 0 %if diagonal has a zero, zeros vector contains 1
           diag_zeros = diag_zeros + 1;
        end
        for m = k:1:stpK          %and every value 'down'
            if (k == 1)  
                if (Kip_finalabI(k,k) > Kip_finalabI(k+m,k))
                    total = total + 1;
                end
            end
                if (k > 1) && (Kip_finalabI(k,k) > Kip_finalabI(m+1,k))
                    total = total + 1;
                end
        end

        if k > 1           
            for M = k:-1:2          % and every value up
                    if (Kip_finalabI(k,k) > Kip_finalabI(M-1,k))
                        total = total + 1;
                    end
            end
        end
        if toad < total
            if Kip_finalabI(k,k) ~= 0
                Kip_ESS(1,a) = KipIndex(k); 
                k_value_history(1,k) = k_value;
                k_value = k;
                inv_tot(1,a) = stpK - total;
                toad = total;
                %if invasion unsuccessful, ESS rate is the resident find max
                %number of score increases
                %as long as res|res isnt 0 
            end
        elseif (toad == total) && (k > 1)
            if Kip_finalabI(k,k) ~= 0
%                if Kip_finalabI(k-1,k-1) ~= 0
                    if k_value == k-1    
                        double_count(1,toad+1) = double_count(1,toad+1) + 1;        %counts the amount of times that two consecuitive totals have occured
                        same_tot(toad+1,double_count(1,toad+1)) = k_value; %makes note of the k value each time a double occurs at this specific total value
                        same_tot(toad+1,1 + double_count(1,toad+1)) = k;
                        low_k = (same_tot(toad+1,1));
                        Kip_ESS(1,a) = (KipIndex(k) + KipIndex(low_k)) / 2; %finds the gamESS  by averaging the two values either end of the consecutive string
                        k_value = k;    
                        inv_tot(1,a) = stpK - total;
                        toad = total;
 %                   end
                    %if two consecutive columns have same value of unsuccessful invasions, ESS found 
                    %if they are not consecutive apart, lower (old) value is selected
                        %this is inverted to production as bias is seen
                    %as long as either value doesnnt give a 0 res|res
                end
            end
        end
    total = 0;
    end
    toad = 0;
end 
sigma = 1 - ((Kip_ESS - 0.05)/(3-0.05));
sigmaIndex = 1 - ((KipIndex - 0.05)/(3-0.05));



%% plot kip ESS with frequency

    figure(1)
    subplot(3,2,3:6)
   plot((sigma(1,:)),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
    xlabel('Pathogen frequency')
    ylabel('\sigma _{ESS}')
    xticks([0 1 (stpF+1)*0.25 (stpF+1)*0.5 (stpF+1)*0.75 (stpF+1)])
    xticklabels({'','','','','',FreqF})
    yticks([0 0.5 1])
    yticklabels({'0' '' '1'})
    xlim([2 stpF+1])
    ylim([-0.05 1.05])
    title('Effect of Pathogen Frequency on Spectrum ESS')
    hold on
    %plot(Kip_ESS_AbsP(1,:))
    %plot(MaxAbs)
    %legend('GamIR_{ESS}','Abundance at ESS')

    %subplot(3,2,1:2)
   % plot(Kip_ESS_AbsI(1,:))
    %xticks([0 1 (stpF+1)*0.25 (stpF+1)*0.5 (stpF+1)*0.75 (stpF+1)])
    %xticklabels({'','','','','',FreqF})
   % title('Final Abundance of focal strain at \sigma _{ESS}')
    hold off

    %% spectrum PIP

% make the PIP
%set colour to use by making a colormap
    hexvals1 = ['#3B4244';'#FFFFFF';'#EA696D']; %black - white - red
    hexvals2 = ['#FFFFFF';'#FCF3F1';'#E11C1C']; %white - red
    hexvals3 = ['#FFFFFF';'#F4E1FF';'#401858']; %white - purple
    
    V1 = hex2rgb(hexvals1);
    V2 = hex2rgb(hexvals2);
    V3 = hex2rgb(hexvals3);
    
    X1 = [0 0.5 1];
    X2 = [0 0.1 1];
    X3 = [0 0.1 1];
    
    Xq = linspace(0, 1, 255);
    Vq1 = interp1(X1, V1, Xq, 'pchip');
    Vq2 = interp1(X2, V2, Xq, 'pchip');
    Vq3 = interp1(X3, V3, Xq, 'pchip');
    
  
F = 1;      %pick freq

IA = Kip_MasterCellI{1,F}./diag(Kip_MasterCellI{1,F})';

    figure(F)
    %subplot(3,2,3:6)
    imagesc(IA)
    caxis([0.999999995 1.000000005])      %colour limits around 1 - creates central diagonal
    colormap(Vq1)
    set(gca, 'YDir','normal','XDir','normal')
    xlabel('Resident Spectrum (\it\sigma_{res})')
    ylabel('Mutant Spectrum (\it\sigma_{mut})')
  
    xticks([ 0 0.5 (stpK)*0.25 (stpK)*0.5 (stpK)*0.75 stpK+1.5])
    xticklabels({'','0','','','','1'})
    xlim([0.5 stpK+1.5])
    yticks([0 0.5 (stpK)*0.25 (stpK)*0.5 (stpK)*0.75 stpK+1.5])
    yticklabels({'','0','','','','1'})
    ylim([0.5 stpK+1.5])
 
 





