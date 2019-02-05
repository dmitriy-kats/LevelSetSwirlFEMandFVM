function kappa=getkappa(s,eq_numb)
if eq_numb==20
    kappa=1/2;
elseif eq_numb==21
    if s>0
        kappa=1;
    elseif s<0
        kappa=0;
    elseif s==0
        kappa=1/2;
    end
elseif eq_numb==22
    if s>0
        kappa=0;
    elseif s<0
         kappa=1;
    else
        kappa=1/2;
    end   
end
end