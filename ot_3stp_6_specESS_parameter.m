
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
    xlabel('GamP')
    ylabel('\sigma _{ESS}')
    xticks([0 1 (stpP+1)*0.25 (stpP+1)*0.5 (stpP+1)*0.75 (stpP+1)])
    xticklabels({P0,'','','','',PF})
    yticks([0 0.5 1])
    yticklabels({'0' '' '1'})
    xlim([1 stpP+1])
    ylim([-0.05 1.05])
    


