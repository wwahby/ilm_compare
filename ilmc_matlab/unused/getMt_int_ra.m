function Mt_int_ra = getMt_int_ra(Ns,l,Mt_intra_corr,use_corrected)

if use_corrected == 0
    if (l == 0)
        Mt_int_ra = Ns;
    elseif ( (0 < l) && (l < sqrt(Ns) ) )
        Mt_int_ra = ( 2*Ns*l - 2*sqrt(Ns)*(l^2) + 1/3*l^3 );
    elseif ( (sqrt(Ns) <= l) && (l < 2*sqrt(Ns)-1) )
        Mt_int_ra = 1/3 * (2*sqrt(Ns) - l)^3;
    else
        Mt_int_ra = 0;
    end
else
    if l < 0
        Mt_int_ra = 0;
    elseif l >= 2*sqrt(Ns)
        Mt_int_ra = 0;
    else
        Mt_int_ra = Mt_intra_corr(l+1);
    end
end
