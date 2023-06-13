
%% Set the range of tested values for gamma and kip 

GamI0 = 0.0;              %production rate range
GamIF = 0.2;
GamR0 = GamI0;
GamRF = GamIF;
stpG = 99;


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
                            eventfunc = @(t,y) abc_NutrSteadyState(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE);
                            optionsode=odeset('Events',eventfunc,'NonNegative', 1:7);
                            y0 = [AbI AbR AbP AbTI AbTR AbTP Nu1];
                            [t,y,te,ye,ie] = ode45(@(t,y)abc_patch_mod(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE),tspan,y0,optionsode);
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

%% Now finds the ESS value of gamma

ai = 0;
xi = 0;
total = 0;
Gam_ESS = zeros(1, stpK+1); 
toad = 0;
toadscore = 0;
toadalert =0;

k_value = zeros(1,stpK+1); %gives the column/row value of the G_ESS
inv_tot = zeros(1,stpK+1); %gives the total number of successful invasions
invalid_G_ESS = zeros(1,stpK+1); %is 1 if the G_ESS is considered invalid

for c = 1:1:stpK+1
    for c3 = 1:1:stpG+1
        for c4 = 1:1:stpG+1
           if Gam_MasterCellI{c,1}(c3,c4) < 1E-9
             Gam_MasterCellI{c,1}(c3,c4) = 0;
           end
        end
    end
end

x = 1;
a = 0;
for y = 1:1:stpK+1
    a = a + 1; 
    Gam_finalabI = Gam_MasterCellI{a,1};
     for k = 1:1:stpG+1      %for every tested prod value
        for m = k:1:stpG          %and every value 'down'
            if (k == 1)  
                if (Gam_finalabI(k,k) > Gam_finalabI(k+m,k))
                    total = total + 1;  %when invasion unsuccessful, increase total
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
                Gam_ESS(1,a) = GamIndex(k); 
                k_value(1,a) = k;
                inv_tot(1,a) = stpG - total;
                toad = total;
                %if unsuccessful invasions is greater than current max, ESS rate is the current resident
                %max number of unsuccessful invasions increases
                %as long as res|res isnt 0 
            end
        elseif (toad == total) && (k > 1)
            if Gam_finalabI(k,k) ~= 0
                if Gam_finalabI(k-1,k-1) ~= 0
                    Gam_ESS(1,a) = (GamIndex(k) + GamIndex(k-1)) / 2; 
                    k_value(1,a) = (k + k-1)/2;
                    inv_tot(1,a) = stpG - total;
                    toad = total;
                    %if two columns have same value of unsuccessful invasions, 
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


%%

a = 0;
for y = 1:1:stpK+1
    a = a + 1;
      if  inv_tot(1,a) > 1
          invalid_G_ESS(1,a) = 1;
      end
end
%labels a GamESS invalid if more than one invasion is successful
%an arbitrary number but mainly to use for checking validity (can arise if
%population goes extinct before a stable value is selected)

%% Plot Production over spectrum

figure(1)
hold on
plot((Gam_ESS(1,:)),Color=[0 0.55 0.9], LineWidth=2) 
%plot(Gam_ESS_Abs(1,:),'--',Color=[0 0.55 0.9])
ylabel('Evolutionarily stable production rate \gamma_{ESS}')
xlabel('Toxin spectrum (\sigma)')
xticks([0 1 (stpK+1)*0.2 (stpK+1)*0.4 (stpK+1)*0.6 (stpK+1)*0.8 stpK+1])
xticklabels({'','0', '','', '', '', '1'})                         
xlim([1 (stpK+1)+0.5])
ylim([GamI0 GamIF])
hold off

%% Plot a PIP of production for each spectrum value


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
   

Kip = 1;    %pick a spectrum value to plot (position in KipIndex)
 
  IA = Gam_MasterCellI{Kip,1}./diag(Gam_MasterCellI{Kip,1});

    figure(1)
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
        
