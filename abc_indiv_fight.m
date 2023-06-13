tspan = [0 tend];
%   running individual fight 
%usually only changing gamma and kip 
GamI =  Gam_ESS(i);
GamR = Gam_ESS(i);
kip = KipIndex(i);
krp = kip;

% E = 5;
% EP = 0;
% GamP = 0.2;
% r(1,3) = 1;
% Nu1 = 1;
% krp = 0.05;
% kpi = 0.05;
% kpr = 0.05;
% kir = 0.05;
% kri = kir;
% AbR = AbI;
% AbI = AbR;
% AbP = AbI;

        invSigma = 1 - ((kip - 0.05)/(3-0.05)); %gives the spectrum values used
        resSigma = 1 - ((krp - 0.05)/(3-0.05));
eventfunc = @(t,y)abc_NutrSteadyState(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE,nut);
optionsode=odeset('Events',eventfunc,'NonNegative', 1:7);
y0 = [AbI AbR AbP AbTI AbTR AbTP Nu1];
[t,y] = ode45(@(t,y)abc_patch_mod(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE,nut),tspan,y0,optionsode);
K = [t,y];
absI = y(:,1)';
absR = y(:,2)';
absP = y(:,3)';

%%normalising to a max of 1
r_tot = [max(absI) max(absR) max(absP)];
r_tot_max = max(r_tot);

relI = (absI./r_tot_max);
relR = (absR./r_tot_max);
relP = (absP./r_tot_max);

toxI = y(:,4)';
toxR = y(:,5)';
toxP = y(:,6)';
r_tot_tox = [max(toxI + toxR + toxP)];
r_tot_tox_max = max(r_tot_tox);

rtoxI = toxI./r_tot_tox_max;
rtoxR = (toxR./r_tot_tox_max);
rtoxP = toxP./r_tot_tox_max;

Nu = y(:,7)';
r_Nu = Nu./max(Nu);
%% plot

figure(1)

a = area((r_Nu), 'lineWidth',3);
hold on
a.FaceColor = 1/255*[245 254 251];
a.EdgeColor = 1/255*[186 215 211];
plot(relR,'-', 'Color', 1/255*[100 158 185], 'lineWidth',6.0)
plot(relI,'-', 'Color', 1/255*[204 105 102], 'lineWidth',6.0)
plot(relP, 'Color', 1/255*[232 191 59], 'lineWidth',6.0)
plot(rtoxR,'--', 'Color', 1/255*[100 158 185], 'lineWidth',3)
plot(rtoxI, '--', 'Color', 1/255*[204 105 102], 'lineWidth',3.5)
plot(rtoxP,'--', 'Color', 1/255*[232 191 59], 'lineWidth',3.0)

yticks([0 0.2 0.4 0.6 0.8 1])
ylabel('Relative abundance')
xlabel('Time (AU)')
tfinal = size(t);
xlim([1 tfinal(1)])
%xlim([1 400])
ylim([0 1.05])
legend(' Nutrient',' strain A',' strain B',' strain C',' A toxin',' B toxin','C toxin','Location','southoutside','Orientation', 'vertical' ,LineWidth=1.0)
hold off

endvalI = K(end,2)
endvalR = K(end,3)
endlavP = K(end,4)
maxtoxI = max(toxI);
maxtoxR = max(toxR);
maxtoxP = max(toxP);



            