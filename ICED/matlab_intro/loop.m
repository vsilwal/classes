stn = {'dehradun', 'roorkee', 'haridwar'} ;
val = [3 6 8];

xvec = 1:1:20;

out = zeros(1,length(xvec));
for ii=1:length(xvec)
    out(ii) = 0;
    for jj = 1:ii
        out(ii) = out(ii) + xvec(jj);
    end
end

%----
for ii = 1:length(xvec)
    out2(ii) = sum(xvec(1:ii));
end
