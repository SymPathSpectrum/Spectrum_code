tspan = [0 tend];
% E = 5;
% Nu2 = 1;
% Ov = 0.2;
% kn2 = 5;
% AbP = AbIndex(1)
%  GamR = 0.2;
%  GamI = 0.05;
% GamP = 0.0;
% kip = 0.05;
% krp = krpIndex(1);
% kip = 1.6479;
% kir = 30000;
% kpi = 2;
% kpr = 2;
% kir = 3.05 - kip;
% kri = 0.05;
% AbR = 1e-04;
% AbI = AbR;
% AbP = AbI;


eventfunc = @(t,y)N2_probio_NutrSteadyState(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE,kn2,Ov);
optionsode=odeset('Events',eventfunc,'NonNegative', 1:8);
y0 = [AbI AbR AbP AbTI AbTR AbTP Nu1 Nu2];
[t,y] = ode45(@(t,y)N2_probio_patch_mod(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE,kn2,Ov),tspan,y0,optionsode);
K = [t,y];
absI = y(:,1)';
absR = y(:,2)';
absP = y(:,3)';
 r_tot = [max(absI) max(absR) max(absP)];
 r_tot_max = max(r_tot);

  relI = absI./r_tot_max;
  relR = (absR./r_tot_max);
  relP = (absP./r_tot_max)+0.015;

%%

toxI = y(:,4)';
toxR = y(:,5)';
toxP = y(:,6)';
    r_tot_tox = [max(toxI + toxR + toxP)];
  r_tot_tox_max = max(r_tot_tox);

rtoxI = toxI./r_tot_tox_max;
rtoxR = (toxR./r_tot_tox_max)+0.05;
rtoxP = toxP./r_tot_tox_max;

Nu = y(:,7)';
r_Nu = Nu./max(Nu);
Nu2 = y(:,8)';
r_Nu2 = Nu2./max(Nu2);
%%
figure(2)


a = area((r_Nu), 'lineWidth',2);
hold on
a2 = area((r_Nu2), 'lineWidth',2);

a.FaceColor = 1/255*[242 255 238];
a.EdgeColor = 1/255*[160 232 134];
a2.FaceColor = 1/255*[254 245 251];
a2.EdgeColor = 1/255*[226 203 223];
plot(relR,'-', 'Color', 1/255*[100 158 185], 'lineWidth',5.0)
plot(relI,'-', 'Color', 1/255*[204 105 102], 'lineWidth',5.0)
%plot(relR+0.005,'-', 'Color', 1/255*[100 158 185], 'lineWidth',2.0)
plot(relP, 'Color', 1/255*[232 191 59], 'lineWidth',5.0)
plot(rtoxR,'--', 'Color', 1/255*[100 158 185], 'lineWidth',2)
plot(rtoxI, '--', 'Color', 1/255*[204 105 102], 'lineWidth',2)

%plot(rtoxP,'--', 'Color', 1/255*[232 191 59], 'lineWidth',4.0)

yticks([0 0.2 0.4 0.6 0.8 1])
ylabel('Relative abundance')
xlabel('Time (AU)')
tfinal = size(t);
xlim([1 tfinal(1)])
xlim([1 1000])
ylim([0 1.05])
hold off

%legend('Nutrient','focal','sibling','third strain', 'toxF','toxS','Location','east')

endvalI = K(end,2)
endvalR = K(end,3)
endlavP = K(end,4)

            