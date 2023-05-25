%% probiotic 

function dydt = probio_patch_mod(t,y,r,GamI,GamR,GamP,Degr,kn1,kpi,kpr,kip,krp,kri,kir,E,EP,HCE)
%GamI/R/P = toxin investment
%Ovrlp = niche overlap
%Degr = toxin degredation/absorption constant
%E = toxin potency/efficacy
%HCE = hill coefficient
%kn1/2 = binding affinity for nutrient 1/2
% Kir = kri  - max affinity of toxin binding
%kpi/r = pathogen toxin affinity for invader/resident
%ki/rp = invader/resident toxin affinity for pathogen

dydt = zeros(1,7);
dydt = dydt';
%row vector to contain the equations

%set values to 0 if sufficiently low
if y(1) < 1E-9
    y(1) = 0;
end

if y(2) < 1E-9
    y(2) = 0;
end

if y(3) < 1E-9
    y(3) = 0;
end

if y(4) < 1E-9
    y(4) = 0;
end

if y(5) < 1E-9
    y(5) = 0;
end

if y(6) < 1E-9
    y(6) = 0;
end

if y(7) < 1E-9
    y(7) = 0;
end


%monod equations for simplicity 

monoN1 = (y(7))/((y(7))+kn1);           %Nutrient 1
monoTir = y(4)^HCE/((y(4)^HCE)+kir^HCE);     %invader toxin on resident
monoTri = y(5)^HCE/((y(5)^HCE)+kri^HCE);     %resident toxin on invader
monoTpi = y(6)^HCE/((y(6)^HCE)+kpi^HCE);    %pathogen toxin on invader
monoTpr = y(6)^HCE/((y(6)^HCE)+kpr^HCE);    %pathogen toxin on resident
monoTip = y(4)^HCE/((y(4)^HCE)+kip^HCE);    %invader toxin on pathogen
monoTrp = y(5)^HCE/((y(5)^HCE)+krp^HCE);    %resident toxin on pathogen  

%Probiotic (= focal strain)
dydt(1) = y(1)*(monoN1*r(1)*(1-GamI) - E*monoTri - EP*monoTpi);

%Resident (= coevolving strain)
dydt(2) = y(2)*(monoN1*r(2)*(1-GamR) - E*monoTir - EP*monoTpr);

%Pathogen 
dydt(3) = y(3)*(monoN1*r(3)*(1-GamP) - E*monoTip - E*monoTrp);

%I toxin
dydt(4) = GamI*y(1)*monoN1 - Degr*(y(2)*monoTir + y(3)*monoTip);

%R toxin 
dydt(5) = GamR*y(2)*monoN1 - Degr*(y(1)*monoTri + y(3)*monoTrp);

%P toxin 
dydt(6) = GamP*y(3)*monoN1 - Degr*(y(1)*monoTpi + y(2)*monoTpr);

%Nutrient 1
dydt(7) = -(monoN1)*(y(1) + y(2) + y(3));



