%KipIndex = KirIndex;

xi2 = 0;
total2 = 0;
toad2 = 0;
toadscore2 = 0;
toadalert2 = 0;

c = 0; c2 = 0; c3 = 0; c4 = 0;
    for c6 = 1:1:stpK+1
        for c7 = 1:1:stpK+1
           if Kip_finalabI(c6,c7) < 1E-9
             Kip_finalabI(c6,c7) = 0;
           end
           if Kip_finalabR(c6,c7) < 1E-9
             Kip_finalabR(c6,c7) = 0;
           end
           if Kip_finalabP(c6,c7) < 1E-9
             Kip_finalabP(c6,c7) = 0;
           end
        end
    end
    
k2 = 0; m2 = 0; M2 = 0;
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
                    Kip_ESS = KipIndex(k2);
                end
            elseif (toad2 == total2) && (k2 > 1)
                if Kip_finalabI(k2,k2) ~= 0
                    if Kip_finalabI(k2-1,k2-1) ~= 0
                        Kip_ESS = (KipIndex(k2) + KipIndex(k2-1))/2;
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
   %%     
sigma = 1 - ((Kip_ESS - 0.05)/(3-0.05));
sigmaIndex = 1 - ((KipIndex - 0.05)/(3-0.05));



Kip_ESS




