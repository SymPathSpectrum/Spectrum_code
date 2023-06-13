
%% set looping values

 stpF = 10;     %how many frequency steps (slows down spec_scan later!)
 
 Freq0 = 0.0;
 FreqF = 1.0;  

stpG = 6;
GamI0 = 0.0;
GamR0 = GamI0;
GamIF = 0.2;
GamRF = GamIF;

stpK = 4;  %can be 1 if you just want to test broad and narrow
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

Gam_MasterCellI = cell(stpK+1, stpF+1);
G_none_MasterCellI = cell(stpK+1, stpF+1);
G_path_MasterCellI = cell(stpK+1, stpF+1);

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
                            AbP = 0;        %contest wihtout pathogen
                            eventfunc = @(t,y)abc_NutrSteadyState(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE);
                            optionsode=odeset('Events',eventfunc,'NonNegative', 1:7);
                            y0 = [AbI AbR AbP AbTI AbTR AbTP Nu1];
                            [t,y,te,ye,ie] = ode45(@(t,y)abc_patch_mod(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE),tspan,y0,optionsode);
                            G_none = [t,y];
                             if G_none(end,2) > 1E-9
                                G_none_finalabI(igi,igri) = G_none(end,2); %contest with no pathogen
                            else
                                G_none_finalabI(igi,igri) = 0; %zeros values that are too small
                            end  

                            AbP = AbI;  %contest with pathogen       
                            eventfunc = @(t,y)abc_NutrSteadyState(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE);
                            optionsode=odeset('Events',eventfunc,'NonNegative', 1:7);
                            y0 = [AbI AbR AbP AbTI AbTR AbTP Nu1];
                            [t,y,te,ye,ie] = ode45(@(t,y)abc_patch_mod(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE),tspan,y0,optionsode);
                            G_path = [t,y];
                            if G_path(end,2) > 1E-9
                                G_path_finalabI(igi,igri) = G_path(end,2); 
                            else
                                G_path_finalabI(igi,igri) = 0;
                            end  
                    end
            end    
            for freq = Freq0:stpszF:FreqF  %now calculate average based on presence and absence
                iF = iF + 1;
                igi = 0;
                Gam_finalabI = (1-freq)*G_none_finalabI + (freq)*G_path_finalabI;
                Gam_MasterCellI{ik,iF} = Gam_finalabI;
                G_none_MasterCellI{ik,iF} = G_none_finalabI;
                G_path_MasterCellI{ik,iF} = G_path_finalabI;
            end    
            toc 
            countdown = countdown - 1

end
 
load splat
sound(y,Fs)

%% Gamma ESS calc - alternative to normal
%this finds the point with the fewest invasions 'below'
%and the fewest invasions 'above'
%and averages these
%useful for when there are areas of mutual invasibility (often seen)

ai = 0;
xi = 0;
totalup = 0;
totaldown = 0;
new_Gam_ESS = zeros(stpK+1,stpF+1); 
toad = 0;
inv_toad_up = stpG;
inv_toad_down = stpG;
toadscore = 0;
toadalert =0;
diag_zeros = zeros(stpK+1,stpF+1); %shows how many zeros on the diagonal there are per invasion test
k_value_up = 0;
k_value_down = 0;
inv_tot = zeros(stpK+1,stpF+1); %gives the total number of successful invasions
invalid_G_ESS = zeros(stpK+1,stpF+1); %is 1 if the G_ESS is considered invalid
inv_tot_down = 0;
inv_tot_up = 0;

%%
% 
for c = 1:1:stpK+1
for c2 = 1:1:stpF + 1
    for c3 = 1:1:stpG+1
        for c4 = 1:1:stpG+1
           if Gam_MasterCellI{c,c2}(c3,c4) < (1E-9)*(FreqIndex(c2))   %zeros any value that could not arise from the frequency calculation 
             Gam_MasterCellI{c,c2}(c3,c4) = 0;
           end
        end
    end
end
end
%%

a = 0;
for y = 1:1:stpK+1
    a = a + 1;
    xi = 0;
for x = 1:1:stpF+1    %go through each spectrum value matrix
    Gam_finalabI = Gam_MasterCellI{a,x};
    xi = xi + 1;
     for k = 1:1:stpG+1      %for every tested prod value
        if Gam_finalabI(k,k) == 0 %if diagonal has a zero, zeros vector contains 1
           diag_zeros = diag_zeros + 1;
        end
        for m = k:1:stpG          % every value 'down'
            if (k == 1)  
                if (Gam_finalabI(k,k) > Gam_finalabI(k+m,k))
                    totaldown = totaldown + 1;
                end
            end
                if (k > 1) && (Gam_finalabI(k,k) > Gam_finalabI(m+1,k))
                    totaldown = totaldown + 1;
                end
        end
        inv_tot_down(1,k) = ((stpG+1)-k) - totaldown;   %number of successful invasions at values 'below' in abundance table

        if k > 1           
            for M = k:-1:2          % and every value up
                    if (Gam_finalabI(k,k) > Gam_finalabI(M-1,k))
                        totalup = totalup + 1;
                    end
            end
        end

        inv_tot_up(1,k) = (k-1) - totalup;  %number of successful invasions at values 'above' in abundance table
        
            if Gam_finalabI(k,k) ~= 0
                if k ~= 1
                if inv_tot_up(1,k) < inv_toad_up     %select lowest invasion total
                    k_value_up = k;         
                    inv_toad_up = inv_tot_up(1,k);   %gives new minimum total
                elseif inv_tot_up(1,k) == inv_toad_up
                    k_value_up = k;
                    inv_toad_up = inv_tot_up(1,k);
                end
                end
            end

            if Gam_finalabI(k,k) ~= 0   %same for invasions 'down'
                if k ~= stpG+1
                    if inv_tot_down(1,k) < inv_toad_down
                        k_value_down = k;
                        inv_toad_down = inv_tot_down(1,k);
                    end
                    if inv_tot_down(1,k) == inv_toad_down
                        if k - k_value_down > 1     %if two maxima are not consecutive, replace the old with the new
                            if inv_tot_down(1,k) ~= inv_tot_down(1,k-1) %this upsets things on the larger scale :(
                               k_value_down = k;
                               inv_toad_down = inv_tot_down(1,k);
                            end
                        end
                    end
                end
            end
        
        if k_value_down - k_value_up > 3
            check_please(a,x) = 1;
        elseif k_value_down - k_value_up < -3
            check_please(a,x) = 1;      %if two points with lowest invasion totals are far apart (arbitrary) add +1 to check_please, in position which describes the spectrum value and freq value
        else 
            check_please(a,x) = 0;
        end
    totalup = 0;
    totaldown = 0;
     end
    inv_toad_up = stpG;
    inv_toad_down = stpG;   
    inv_tot_up(1,:) = stpG;
    inv_tot_down(1,:) = stpG;
    new_Gam_ESS(a,x) = (GamIndex(1,k_value_up) + GamIndex(1,k_value_down))/2;  
    k_value_up = 0;
    k_value_down = 0;
end

end

%% currently unsure how to scan for invalid

%% Plotting ESS production rate against spectrum, lines per frequency value

darkBLUE = [0.6, 0.7176, 0.98686];
lightBLUE = [0.00, 0.00, 0.17058];

blueGRADflex = @(i,N) lightBLUE + (darkBLUE-lightBLUE)*((i-1)/(N-1));
N = stpF+1;

figure(1)
hold on

for i = 1:N
    plot(new_Gam_ESS(:,i),'color',blueGRADflex(i,N),'LineWidth',3.0);
    hold on
end

%plot((Gam_ESS(:,:)))
ylabel('\gamma_{ESS}')
xlabel('Spectrum')
xticks([0 1 (stpG+1)*0.25 (stpG+1)*0.5 (stpG+1)*0.75 stpG+1])
xticklabels({'','narrow','','broad','narrow','narrow'})
yticks([0 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2 0.3 0.4 0.5 0.6  0.7 0.8 0.9 1.0])
yticklabels({'0' '0.02' '0.04' '0.06' '0.08' '0.1' '0.12' '0.14' '0.16' '0.18' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7' '0.8' '0.9' '1'})
xlim([0 50.05])
ylim([-0.005 0.2])
hold on
title('Effect of pathogen frequency on \gamma_{ESS}')

%% Production rate PIPs for spectrum and freq

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
   

F = 1;  %which frequency
Kip = 1;    %which kip
 
  IA = Gam_MasterCellI{Kip,F}./diag(Gam_MasterCellI{Kip,F})';

    figure(1)
    imagesc(IA)
    clim([0.9999999999999995 1.00000000000000005])      %colour limits around 1 - creates central diagonal
    colormap(Vq1)
    set(gca, 'YDir','normal','XDir','normal')
    xlabel('Resident production rate (\gamma_{res})')
    ylabel('Mutant production rate (\gamma_{mut})')
    xticks([0 0.5 (stpG+1)*0.25 (stpG+1)*0.5 (stpG+1)*0.75 (stpG+1.5)*1])
    xticklabels({'','0','','','',GamIF})
    xlim([0.5 (stpG+1.5)])
    yticks([0 0.5 (stpG+1)*0.25 (stpG+1)*0.5 (stpG+1)*0.75 (stpG+1.5)*1])
    yticklabels({'','0','','','',GamIF})
    ylim([0.5 (stpG+1.5)])
    %end









