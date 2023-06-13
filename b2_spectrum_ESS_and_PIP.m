%% competeing all spectrum values using calculated gamma ESS at each param value

Kip_finalabI = zeros(stpK+1, stpK+1);
Kip_MasterCellI = cell(1, stpP+1);
Kip_finalabR = zeros(stpK+1, stpK+1);
Kip_MasterCellR = cell(1, stpP+1);
Kip_finalabP = zeros(stpK+1, stpK+1);
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
                    eventfunc = @(t,y)abc_NutrSteadyState(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE);
                    optionsode=odeset('Events',eventfunc,'NonNegative', 1:7);
                    y0 = [AbI AbR AbP AbTI AbTR AbTP Nu1];
                    [t,y,te,ye,ie] = ode45(@(t,y)abc_patch_mod(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE),tspan,y0,optionsode);
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

%% calculate kipESS

xi2 = 0;
Kip_ESS = zeros(1,stpP+1); %matrix, Kess for each abP
total2 = 0;
toad2 = 0;
toadscore2 = 0;
toadalert2 = 0;

c = 0; c2 = 0; c3 = 0; c4 = 0;
for c = 1:1:stpP+1
    for c6 = 1:1:stpK+1
        for c7 = 1:1:stpK+1
           if Kip_MasterCellI{1,c}(c6,c7) < 1E-9
             Kip_MasterCellI{1,c}(c6,c7) = 0;
           end
        end
    end
end

ai = 0; a = 0; k2 = 0; m2 = 0; M2 = 0;
for a = 1:1:stpP+1   %for each parameter value tested  
    ai = ai + 1;
    Kip_finalabI = Kip_MasterCellI{1,a};
        for k2 = 1:1:stpK+1     %for all tested spectrum values
            for m2 = k2:1:stpK
                if (k2 == 1)
                    if (Kip_finalabI(k2,k2) > Kip_finalabI(k2+m2,k2))
                        total2 = total2 + 1;
                    end
                end
                if (k2 > 1) && (Kip_finalabI(k2,k2) > Kip_finalabI(m2+1,k2))
                    total2 = total2 + 1;
                end
            end
        
            if k2 > 1
                for M2 = k2:-1:2
                    if (Kip_finalabI(k2,k2) > Kip_finalabI(M2-1,k2))
                        total2 = total2 + 1;
                    end
                end
            end
            
            if toad2 < total2
                if Kip_finalabI(k2,k2) ~= 0
                    toad2 = total2;
                    Kip_ESS(1,a) = KipIndex(k2);
                end
            elseif (toad2 == total2) && (k2 > 1)
                if Kip_finalabI(k2,k2) ~= 0
                    if Kip_finalabI(k2-1,k2-1) ~= 0
                        Kip_ESS(1,a) = (KipIndex(k2) + KipIndex(k2-1))/2;
                        toad2 = total2;
                        toadscore2 = toadscore2 + 1;
                        if toadscore2 > 1
                            toadalert2 = k2
                        end
                    end
                end
            end
         total2 = 0;
        end
        toad2 = 0;
end 
sigma = 1 - ((Kip_ESS - 0.05)/(3-0.05));
sigmaIndex = 1 - ((KipIndex - 0.05)/(3-0.05));


Kip_ESS_AbsI = zeros(1,stpP+1);

%% plot sigma ess against parameter value

    figure(1)
   plot((sigma(1,:)),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
    xlabel('parameter')
    ylabel('\sigma _{ESS}')
    xticks([0 1 (stpP+1)*0.25 (stpP+1)*0.5 (stpP+1)*0.75 (stpP+1)])
    xticklabels({P0,'','','','',PF})
    yticks([0 0.5 1])
    yticklabels({'0' '' '1'})
    xlim([1 stpP+1])
    ylim([-0.05 1.05])

    %% plot spectrum PIPs

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
    
  
p = 1;      %choose paramater value sweeping a parameter

IA = Kip_MasterCellI{1,p}./diag(Kip_MasterCellI{1,p})';

    figure(p)
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
 
 


    


