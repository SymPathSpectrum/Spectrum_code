
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
   

param = 1;  %if doing a parameter sweep, use this
Kip = 1;    %pick spectrum value
 
  IA = Gam_MasterCellI{Kip,ab}./diag(Gam_MasterCellI{Kip,ab})';

    figure(ab)
    imagesc(IA)
    clim([0.9999999999999995 1.00000000000000005])      %colour limits around 1 - creates central diagonal
    colormap(Vq1)
    set(gca, 'YDir','normal','XDir','normal', 'fontname', 'arial', 'fontsize', 25)
    xlabel('Resident production rate (\gamma_{res})')
    ylabel('Mutant production rate (\gamma_{mut})')
    xticks([0 0.5 (stpG+1)*0.25 (stpG+1)*0.5 (stpG+1)*0.75 (stpG+1.5)*1])
    xticklabels({'','0','','','',GamIF})
    xlim([0.5 (stpG+1.5)])
    yticks([0 0.5 (stpG+1)*0.25 (stpG+1)*0.5 (stpG+1)*0.75 (stpG+1.5)*1])
    yticklabels({'','0','','','',GamIF})
    ylim([0.5 (stpG+1.5)])
    %end
   

