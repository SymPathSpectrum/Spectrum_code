%% finding ESSs and sweeping an extra parameter

%%
% P is a placeholder for whatever parameter you change 

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
Gam_MasterCellI = cell(stpK+1, stpP+1);
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
                            eventfunc = @(t,y) abc_NutrSteadyState(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE);
                            optionsode=odeset('Events',eventfunc,'NonNegative', 1:7);
                            y0 = [AbI AbR AbP AbTI AbTR AbTP Nu1];
                            [t,y,te,ye,ie] = ode45(@(t,y)abc_patch_mod(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE),tspan,y0,optionsode);
                            G = [t,y];
                            Gam_finalabI(igi,igri) = G(end,2); %generate a gamxgam matrix of final abI
                            Gam_finalabP(igi,igri) = G(end,4); %generate a gamxgam matrix of final abP
                        
                    end
            end    
            Gam_MasterCellI{ik,iP} = Gam_finalabI;
            Gam_MasterCellP{ik,iP} = Gam_finalabP;
            toc 
            countdown = countdown - 1
   end
end
 
load splat
sound(y,Fs)

%% calc GammaESS for all kip and paramater combinations

ai = 0;
xi = 0;
total = 0;
Gam_ESS = zeros(stpK+1,stpP + 1); 
toad = 0;
toadscore = 0;
toadalert =0;
k_value = zeros(stpK+1,stpP+1); %gives the column/row value of the G_ESS
inv_tot = zeros(stpK+1,stpP+1); %gives the total number of successful invasions
invalid_G_ESS = zeros(stpK+1,stpP+1); %is 1 if the G_ESS is considered invalid


for c = 1:1:stpK+1
for c2 = 1:1:stpP + 1
    for c3 = 1:1:stpG+1
        for c4 = 1:1:stpG+1
           if Gam_MasterCellI{c,c2}(c3,c4) < 1E-9
             Gam_MasterCellI{c,c2}(c3,c4) = 0;
           end
        end
    end
end
end

a = 0;
for y = 1:1:stpK+1
    a = a + 1;
    xi = 0;
for x = 1:1:stpP+1    %go through each spectrum and parameter value matrix
    Gam_finalabI = Gam_MasterCellI{a,x};
    xi = xi + 1;
    for k = 1:1:stpG+1      %for every tested prod value
        for m = k:1:stpG          %and every value 'down'
            if (k == 1)  
                if (Gam_finalabI(k,k) > Gam_finalabI(k+m,k))
                    total = total + 1;
                end
            end
                if (k > 1) && (Gam_finalabI(k,k) > Gam_finalabI(m+1,k))
                    total = total + 1;
                end
        end

        if k > 1           
            for M = k:-1:2          % and every value up
                    if (Gam_finalabI(k,k) > Gam_finalabI(M-1,k))
                        total = total + 1;
                    end
            end
        end
        if toad < total
            if Gam_finalabI(k,k) ~= 0
                Gam_ESS(a,x) = GamIndex(k); 
                k_value(a,x) = k;
                inv_tot(a,x) = stpG - total;
                toad = total;
                %if invasion unsuccessful, ESS rate is the resident find max
                %number of score increases
                %as long as res|res isnt 0 
            end
        elseif (toad == total) && (k > 1)
            if Gam_finalabI(k,k) ~= 0
                if Gam_finalabI(k-1,k-1) ~= 0
                    Gam_ESS(a,x) = (GamIndex(k) + GamIndex(k-1)) / 2; 
                    k_value(a,x) = (k + k-1)/2;
                    inv_tot(a,x) = stpG - total;
                    toad = total;
                    %if two columns have same value of unsuccessful invasions, ESS found 
                    %ESS given as average, prod rate at which success changed
                    %as long as either value doesnnt give a 0 res|res
                    toadscore = toadscore + 1; % count the amount of times the max total is reached after once
                    if toadscore > 1
                        toadalert = x,k; %release alert showing which columns have repeating totals for manual inspection
                    end
                end
            end
        end
    total = 0;
    end
    toad = 0;
end
end
 
%% check for invalid ESS

a = 0;
for y = 1:1:stpK+1
    a = a + 1;
    for x = 1:1:stpP+1      
        if inv_tot(a,x) > 1
            invalid_G_ESS(a,x) = 1;
        else 
            invalid_G_ESS(a,x) = 0;
        end
    end
end


%% plotting gamma ess over parameter (xaxis), each line is a different spectrum value

figure(1)
hold on
plot((Gam_ESS(1,:)),Color=[0 0.55 0.9], LineWidth=5)
plot((Gam_ESS(2,:)),Color=[0.8 0.35 0.21],LineWidth=5) %if only testing broad and narrow
ylabel('\gamma_{ESS}')
xlabel('parameter')
xticks([1 50])
xticklabels({'0' PIndex(stpP+1)})
yticks([0 0.1 0.2 0.3 0.4 0.5 0.6  0.7 0.8 0.9 1.0])
yticklabels({'0' '0.1' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7' '0.8' '0.9' '1'})
xlim([1 stpE+1])
ylim([0 0.95])

title('Effect of nutrient abundance on \gamma_{ESS}')
%legend('Narrow toxin','abundance using ESS', 'Broad toxin','abundance using ESS')
legend('Narrow toxin', 'Broad toxin')


%% Plotting production rate PIPs


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
   
p = 1;   %pick parameter value
Kip = 1;    %pick spectrum value
 
  IA = Gam_MasterCellI{Kip,p}./diag(Gam_MasterCellI{Kip,p})';

    figure(p)
    imagesc(IA)
    clim([0.9999999999999995 1.00000000000000005])      %colour limits around 1 - creates central diagonal
    colormap(Vq1)
    set(gca, 'YDir','normal','XDir','normal')
    xlabel('Resident production rate (\gamma_{res})')
    ylabel('Mutant production rate (\gamma_{mut})')
    xticks([0 0.5 (stpG+1)*0.25 (stpG+1)*0.5 (stpG+1)*0.75 (stpG+1.5)*1])
    xlim([0.5 (stpG+1.5)])
    yticks([0 0.5 (stpG+1)*0.25 (stpG+1)*0.5 (stpG+1)*0.75 (stpG+1.5)*1])
    ylim([0.5 (stpG+1.5)])

  




