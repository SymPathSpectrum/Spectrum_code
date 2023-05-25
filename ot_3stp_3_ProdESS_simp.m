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
%population goes extinct before a stable value is selected
%% Plot!

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

