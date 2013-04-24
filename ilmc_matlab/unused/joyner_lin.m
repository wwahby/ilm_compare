function a = joyner_lin(l,Nt)
    if l < 0
        a = 0;
    elseif l==0
        a = 1/Nt;
    elseif l < sqrt(Nt)
        a = 1/Nt^2 * (4*Nt*l - 4*sqrt(Nt)*l^2 + 2/3*l^3);
    elseif l < 2*sqrt(Nt) - 1
        a = 2/3/Nt^2 * (2*sqrt(Nt) - l)^3;
    else
        a = 0;
    end