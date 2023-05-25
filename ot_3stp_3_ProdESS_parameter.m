%% finding gamma ESS for each spectrum and parameter value

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



%%

for c = 1:1:stpK+1
for c2 = 1:1:stpP + 1
    for c3 = 1:1:stpG+1
        for c4 = 1:1:stpG+1
           if Gam_MasterCell{c,c2}(c3,c4) < 1E-9
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
for x = 1:1:stpP+1    %go through each spectrum and parameter value matrix
    Gam_finalabI = Gam_MasterCell{a,x};
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
%%
a = 0;
for y = 1:1:stpK+1
    a = a + 1;
    for x = 1:1:stpP+1      %go through each spectrum value matrix
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
xlabel('Symbiont toxin efficiency (E)')
xticks([1 25.5 50])
xticklabels({'0' EIndex(25) EIndex(50)})
yticks([0 0.1 0.2 0.3 0.4 0.5 0.6  0.7 0.8 0.9 1.0])
yticklabels({'0' '0.1' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7' '0.8' '0.9' '1'})
xlim([1 stpE+1])
ylim([0 0.95])

%plot(Gam_ESS_Abs(2,:),'--',Color=[0.8 0.35 0.21])
%plot(MaxAbs)
title('Effect of nutrient abundance on \gamma_{ESS}')
%legend('Narrow toxin','abundance using ESS', 'Broad toxin','abundance using ESS')
legend('Narrow toxin', 'Broad toxin')

