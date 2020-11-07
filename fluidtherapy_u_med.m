function [] = fluidtherapy_u_med(vect_vol, vect_serum, vect_drug, vect_foodintake, vect_extraord, vect_insens)

    %Out
    out = sum(vect_vol) + sum(vect_extraord) + sum(vect_insens);

    %In
    in = sum(vect_serum) + sum(vect_drug) + sum(vect_foodintake);

    %Balance
    balance = in - out;

    if balance > -300 && balance < 0
        msgbox('CORRECT BALANCE')
    end

    if balance < -300 
        msgbox('DEHIDRATION','Altert', 'error');
        msgbox('Drink or increase serum');
    end

    if balance > 0
        msgbox('HYDRIC EQUILIBRIUM')
        if sum(vect_serum)~=0
            msgbox('Remove Serum','Altert', 'error');
        end
    end

end