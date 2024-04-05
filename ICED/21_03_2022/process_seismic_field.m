clear all
close all

% select example
iex = 1 ;

% Define filename
ddir='/Users/vipul/Desktop/REPOSITORIES/classes/ICED/21_03_2022/';
fname = 'MC405';

% read data
fid=fopen(strcat(ddir,fname));
C=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','Headerlines',1);

% plot seismogram
figure
scale_factor = 3;
for ii=2:25
    c_max = max(C{ii});
    y_normalized = scale_factor * (C{ii}/c_max);
    plot(C{1},y_normalized - (2*(ii-2)),'b')
    hold on, axis tight
end

% plot Envelope
figure
jj=1;
for ii=2:25
    if ii>13; jj=2 ; end
    [up,lo] = envelope(C{ii});
    c_max = max(up);
    y_normalized = scale_factor *(up/c_max);
    subplot(1,2,jj)
    plot(C{1},y_normalized - (2*(ii-2)),'b')
    hold on, axis tight
end

if 0
% plot spectrogram
figure(3); 
for ii=1:9
subplot(3,3,ii)
colormap('jet')
sig = C{ii+1};
nwin = 100;
wind = kaiser(nwin,20);
nlap = nwin-10;
nfft = 512;
Fs = 512;
spectrogram(sig,wind,nlap,nfft,Fs,'yaxis')
end
end