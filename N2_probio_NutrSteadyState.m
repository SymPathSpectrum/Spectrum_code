function [value,isterminal,direction] = N2_probio_NutrSteadyState(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE,kn2,Ov)

dydt = N2_probio_patch_mod(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE,kn2,Ov);

SS = [abs(dydt(7) + dydt(8))];
SS1 = max(SS);
val1 = SS1 - 0.00000005;

value = val1;
isterminal = 1;
direction = -1;

end

%function outputs a negative value when the dN/dt < 5e-8

