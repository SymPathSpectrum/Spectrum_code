tspan = [0 tend];
% EP = 5;
% E = 5;  
% AbI = A;
% GamI = GamIndex(35);
% GamR = GamIndex(1);
% GamP = 0
% kip = KipIndex(1);
% krp = krpIndex(1);
% EP = 0;
% kip = KipIndex(41);
% krp = kip;
% kpi = 2;
% kpr = 2;
% kir = 3.05 - kip;
% kri = kir;
% AbR = 1e-04;
% AbI = AbR;
% AbP = AbR/1000;
% Nu1 = 1;
eventfunc = @(t,y)probio_NutrSteadyState(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE,nut);
optionsode=odeset('Events',eventfunc,'NonNegative', 1:7);
y0 = [AbI AbR AbP AbTI AbTR AbTP Nu1];
[t,y] = ode45(@(t,y)probio_patch_mod(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE,nut),tspan,y0,optionsode);
K = [t,y];
absI = y(:,1)';
absR = y(:,2)';
absP = y(:,3)';
 r_tot = [max(absI) max(absR) max(absP)];
 r_tot_max = max(r_tot);

  relI = absI./r_tot_max;
  relR = (absR./r_tot_max);
  relP = absP./r_tot_max;

%%

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
%%
figure(1)
a = area((r_Nu), 'lineWidth',1.5);
hold on
a.FaceColor = 1/255*[242 255 238];
a.EdgeColor = 1/255*[160 232 134];
plot(relI,'-', 'Color', 1/255*[255 100 10], 'lineWidth',2.0)
plot(relR,'-', 'Color', 1/255*[120 202 255], 'lineWidth',2.0)
plot(relP, 'Color', 1/255*[255 216 34], 'lineWidth',2.0)
plot(rtoxI, '--', 'Color', 1/255*[255 100 10], 'lineWidth',1.0)
plot(rtoxR,'--', 'Color', 1/255*[120 202 255], 'lineWidth',1.0)
plot(rtoxP,'--', 'Color', 1/255*[255 216 34], 'lineWidth',1.0)

yticks([0 0.2 0.4 0.6 0.8 1])
ylabel('Relative abundance')
xlabel('Time (AU)')
tfinal = size(t);
xlim([1 tfinal(1)])
%xlim([1 400])
ylim([0 1.05])
hold off

%legend('Nutrient','focal','sibling','third strain', 'toxF','toxS','Location','east')

endvalI = K(end,2)
endvalR = K(end,3)
endlavP = K(end,4)

            