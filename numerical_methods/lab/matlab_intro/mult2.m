% Conditional statement

lat0 = -70;

if lat0>=0   
    if (lat0<60)
        disp('Hot')
    elseif (lat0>=60)
        disp('cold')
    else
        disp('Moderate')
    end  
end