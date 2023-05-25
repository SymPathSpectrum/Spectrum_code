%% Gamma ESS calc - gets complex. check PIPs to see whether to use this or
%'new'
%this finds point with fewest successful invasions overall.
%also replaces any lower ESS values with higher ones

ai = 0;
xi = 0;
total = 0;
Gam_ESS = zeros(stpK+1,stpF+1); 
toad = 0;
toadscore = 0;
toadalert =0;
diag_zeros = zeros(stpK+1,stpF+1); %shows how many zeros on the diagonal there are per invasion test
k_value = 0;
inv_tot = zeros(stpK+1,stpF+1); %gives the total number of successful invasions
invalid_G_ESS = zeros(stpK+1,stpF+1); %is 1 if the G_ESS is considered invalid

%%
% 
for c = 1:1:stpK+1
for c2 = 1:1:stpF + 1
    for c3 = 1:1:stpG+1
        for c4 = 1:1:stpG+1
           if Gam_MasterCell{c,c2}(c3,c4) < (1E-9)*(FreqIndex(2))   %this zeros any value that could not arise from the frequency calculation 
             Gam_MasterCell{c,c2}(c3,c4) = 0;
           end
        end
    end
end
end

a = 0;
for y = 1:1:stpK+1
    a = a + 1;
    xi = 0;
for x = 1:1:stpF+1   
    Gam_finalabI = Gam_MasterCell{a,x};
    double_count = zeros(1,stpG);
    xi = xi + 1;
    for k = 1:1:stpG+1      %for every tested prod value
        if Gam_finalabI(k,k) == 0 %if diagonal has a zero, zeros vector contains 1
           diag_zeros = diag_zeros + 1;
        end
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
                k_value_history(1,k) = k_value;
                k_value = k;
                inv_tot(a,x) = stpG - total;
                toad = total;
                %if invasion unsuccessful, ESS rate is the resident find max
                %number of score increases
                %as long as res|res isnt 0 
            end
        elseif (toad == total) && (k > 1)
            if Gam_finalabI(k,k) ~= 0
                if Gam_finalabI(k-1,k-1) ~= 0
                    if k_value == k-1      
                        double_count(1,toad) = double_count(1,toad) + 1;        %counts the amount of times that two consecuitive totals have occured
                        same_tot(toad,double_count(1,toad)) = k_value; %makes note of the k value each time a double occurs at this specific total value
                        same_tot(toad,1 + double_count(1,toad)) = k;
                        low_k = (same_tot(toad,1));
                        Gam_ESS(a,x) = (GamIndex(k) + GamIndex(low_k)) / 2; %finds the gamESS  by averaging the two values either end of the consecutive string
                        k_value = k;    
                        inv_tot(a,x) = stpG - total;
                        toad = total;
                    else
                        Gam_ESS(a,x) = GamIndex(k); 
                        k_value = k;
                        inv_tot(a,x) = stpG - total;
                        toad = total;
                    end
                    %if two consecutive columns have same value of unsuccessful invasions, ESS found as average 
                    %if they are not consecutive apart, higher (new) value is selected
                    %as long as either value doesnt have a 0 res|res
                end
            end
        end
    total = 0;
    end
    toad = 0;
end
end
%%
a = 0;
for y = 1:1:stpK+1
    a = a + 1;
    xi = 0;
for x = 1:1:stpF+1      %go through each spectrum value matrix
    Gam_finalabI = Gam_MasterCell{a,x};
    xi = xi + 1;
        if diag_zeros(a,x) ~= 0       %if there are any zeros on the diag
            if inv_tot(a,x) > 1
                invalid_G_ESS(a,x) = 1;
            else 
                invalid_G_ESS(a,x) = 0;
            end
        end
    end
end


%% Plots gramma ESS over spectrum, each line a different frequency

lightBLUE = [0.35686, 0.81176, 0.95686];
darkBLUE = [0.00, 0.00, 0.57058];

blueGRADflex = @(i,N) lightBLUE + (darkBLUE-lightBLUE)*((i-1)/(N-1));
N = 50;

figure(1)
hold on

for i = 1:N
    plot(Gam_ESS(:,i),'color',blueGRADflex(i,N));
    hold on
end

ylabel('\gamma_{ESS}')
xlabel('Spectrum')
xticks([0 1 (stpG+1)*0.25 (stpG+1)*0.5 (stpG+1)*0.75 stpG+1])
xticklabels({'','narrow','','broad','narrow','narrow'})
yticks([0 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2 0.3 0.4 0.5 0.6  0.7 0.8 0.9 1.0])
yticklabels({'0' '0.02' '0.04' '0.06' '0.08' '0.1' '0.12' '0.14' '0.16' '0.18' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7' '0.8' '0.9' '1'})
xlim([0 50.05])
ylim([-0.005 0.2])
hold on
