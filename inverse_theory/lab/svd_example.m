close all
clear all

fname = '/Users/vipul/Downloads/IITR.jpeg';

% load image
s = imread(fname);
X = double(rgb2gray(s));
figure(1); subplot(2,2,1)
imagesc(X);
colormap('gray')

[U,S,V] = svd(X);

% Informaiton contained in major singular values
figure(2)
subplot(1,2,1)
semilogy(diag(S))
subplot(1,2,2)
plot(cumsum(diag(S))/sum(diag(S)))

% truncated svd
r = [2 5 10];
for ii=1:length(r)
    Xapprox = U(:,1:r(ii))*S(1:r(ii),1:r(ii))*V(:,1:r(ii))';
    figure(1); subplot(2,2,ii+1);
    imagesc(Xapprox);
end

[U,S,V]=svd(X,'econ');
Xecon=U*S*V';
imagesc(Xecon)
