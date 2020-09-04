iex =1;
if iex ==1
    for a=10:20
        a = a + 1;
        if (a < 13)
            fprintf('value of a: %d %f \n', a,a);
        elseif a<15
            fprintf('value of a: %d %f \n', a,a);
        elseif a<17
            fprintf('value of a: %d %f \n', a,a);
        else
            fprintf('Nothing \n')
        end
    end
    
elseif iex==2
    a=10;
    while (a<20)
        if (a<13)
            fprintf('Yes \n')
        elseif (a >=13 || a<16)
            fprintf('No \n')
        else
            fprintf('Dont know \n')
        end
        a=a+1;
    end
end