%
% Bootstrap Method
% This function takes N bootstrap samples with replacement of the
% original data and stores the calculation of the desired statistic
% in the returned array.
%
% bootstrap(data,N)
% data is the original array of data to be sampled.
% N is the number of bootstrap samples to be taken.
%

function [stat]=bootstrap(data,N)
rand('state', sum(100*clock));
%reset random generator to different state at each time

n=length(data);
stat=zeros(1,N); %initialize array to be returned of bootstrap estimations
sample=zeros(1,n); %initialize array to hold bootstrap sample

for i=1:N
    choose=round(((n-1)*rand(1,n))+1);
    %choose an array length n of random #'s between [1,n]
    %the sampled data will come from the values at these random indices.
    
    %fill the sample array with values at randomized indices
    for j=1:n
        sample(j)=data(choose(j));
    end;
    stat(i)=mean(sample); %fill stat array with bootstrap estimations
end;

%% to be ran before executing the function:
if 0
clear all;
close all;
%data=[141,156.5,162,159,157,143.5,158,140,142,150,148.5,138.5,161,153,145,147,158.5,160.5,167.5,155,137];
    
data = [41    49    51    38    41    16    41    52    49    43    38    33, ...
    43     8    28    45    28    25    17    28    33    25    39    45,...
    45    14    43    33    21    33    48     7    27    34     8];
stat = bootstrap(data,10000)
hist(stat,30)
end