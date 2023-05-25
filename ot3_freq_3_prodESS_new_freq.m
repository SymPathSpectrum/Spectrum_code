%% alternative ESS calc
%this finds the point with the fewest invasions below
%and the fewest invasions above
%and averages these
%useful for when there are areas of mutual invasibility

ai = 0;
xi = 0;
totalup = 0;
totaldown = 0;
new_Gam_ESS = zeros(stpKip+1,stpF+1); 
toad = 0;
inv_toad_up = stpG;
inv_toad_down = stpG;
toadscore = 0;
toadalert =0;
diag_zeros = zeros(stpKip+1,stpF+1); %shows how many zeros on the diagonal there are per invasion test
k_value_up = 0;
k_value_down = 0;
inv_tot = zeros(stpKip+1,stpF+1); %gives the total number of successful invasions
invalid_G_ESS = zeros(stpKip+1,stpF+1); %is 1 if the G_ESS is considered invalid
inv_tot_down = 0;
inv_tot_up = 0;

%%
% 
for c = 1:1:stpKip+1
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
for y = 1:1:stpKip+1
    a = a + 1;
    xi = 0;
for x = 1:1:stpF+1    %go through each spectrum value matrix
    Gam_finalabI = Gam_MasterCell{a,x};
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
% a = 0;
% for y = 1:1:stpKip+1
%     a = a + 1;
%     xi = 0;
% for x = 1:1:stpF+1      %go through each spectrum value matrix
%     Gam_finalabI = Gam_MasterCell{a,x};
%     xi = xi + 1;
%         if diag_zeros(a,x) ~= 0       %if there are any zeros on the diag
%             if inv_tot_up(a,x) > 1
%                 invalid_G_ESS(a,x) = 1;
%             else 
%                 invalid_G_ESS(a,x) = 0;
%             end
%             if inv_tot_down(a,x) > 1
%                 invalid_G_ESS(a,x) = 1;
%             else 
%                 invalid_G_ESS(a,x) = 0;
%             end
%         end
%     end
% end




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
title('Effect of Third Strain Frequency on \gamma_{ESS}')


