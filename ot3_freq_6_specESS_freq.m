%% finding kip ESS
%currently only using usual method (lowest invasion number) as issues with
%mutual invasibility not seen
%I'll be honest i didn't really get this far since it went downhill before
%getting here

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
           if Kip_MasterCellI{1,c}(c6,c7) < (1E-9)*(FreqIndex(2))
             Kip_MasterCellI{1,c}(c6,c7) = 0;
           end
        end
    end
end

ai = 0; a = 0; k = 0; m = 0; M = 0;
for a = 1:1:stpF+1   %for each fre
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
                    if k_value == k-1    
                        double_count(1,toad+1) = double_count(1,toad+1) + 1;        %counts the amount of times that two consecuitive totals have occured
                        same_tot(toad+1,double_count(1,toad+1)) = k_value; %makes note of the k value each time a double occurs at this specific total value
                        same_tot(toad+1,1 + double_count(1,toad+1)) = k;
                        low_k = (same_tot(toad+1,1));
                        Kip_ESS(1,a) = (KipIndex(k) + KipIndex(low_k)) / 2; %finds the gamESS  by averaging the two values either end of the consecutive string
                        k_value = k;    
                        inv_tot(1,a) = stpK - total;
                        toad = total;
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



%% plot signma over frequency

    figure(1)
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
  

