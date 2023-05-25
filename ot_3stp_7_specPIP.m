
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
    
  
i = 1      %only need to change if sweeping a parameter
IA = Kip_finalabI./diag(Kip_finalabI)';
%IA = Kip_MasterCellI{1,i}./diag(Kip_MasterCellI{1,i})';

    figure(i)
    %subplot(3,2,3:6)
    imagesc(IA)
    caxis([0.999999995 1.000000005])      %colour limits around 1 - creates central diagonal
    colormap(Vq1)
    set(gca, 'YDir','normal','XDir','normal')
    xlabel('Resident Spectrum (\it\sigma_{res})')
    ylabel('Mutant Spectrum (\it\sigma_{mut})')
  
    xticks([ 0 0.5 (stpK)*0.25 (stpK)*0.5 (stpK)*0.75 stpK+1.5])
    xticklabels({'','0','','','','1'})
    xlim([0.5 stpK+1.5])
    yticks([0 0.5 (stpK)*0.25 (stpK)*0.5 (stpK)*0.75 stpK+1.5])
    yticklabels({'','0','','','','1'})
    ylim([0.5 stpK+1.5])
 
 

