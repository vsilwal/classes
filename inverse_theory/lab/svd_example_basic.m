close all, clear all
%% Homogeneous figure
X = ones(657,1800);

figure(1)
subplot(2,2,1)
imagesc(X)
title('homogeneous image')

[U,S,V] = svd(X);

% Informaiton contained in major singular values
figure(2)
subplot(2,2,1)
semilogy(diag(S))

%% Slightly change the figure
X(1:323,1:900) = 0;

figure(1)
subplot(2,2,2)
imagesc(X)
title('two independent vectors');

[U,S,V] = svd(X);

% Informaiton contained in major singular values
figure(2)
subplot(2,2,2)
semilogy(diag(S),'.')
title('two independent vectors');

%% Slightly change the figure
X(:,1:900) = 0;

figure(1)
subplot(2,2,3)
imagesc(X)
title('(basically) one independent vector')

[U,S,V] = svd(X);

% Informaiton contained in major singular values
figure(2)
subplot(2,2,3)
semilogy(diag(S))
title('(basically) one independent vector')

%% Slightly change the figure
X(:,901:1800) = rand(657,900);

figure(1)
subplot(2,2,4)
imagesc(X)
title('random vector appended')

[U,S,V] = svd(X);

% Informaiton contained in major singular values
figure(2)
subplot(2,2,4)
semilogy(diag(S))
title('random vector appended')