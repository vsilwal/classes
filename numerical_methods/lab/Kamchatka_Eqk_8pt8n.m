clc
clear all
%close all

Data=load('2025-07-30-0200-00Snew.TXT'); % loading sgdata file
X=Data(:,1); % Theoretical load data from sg data file
% C=Dat(121931:294730); % data range for 24 hours
% plot(X)

ndays = 2.25;

Cn=zeros(86400*ndays,1);

figure(100)
for n=1
    m=3600*0;
    Cn(:)=X(1+(n-1)*m:86400*ndays+(n-1)*m);
    % storeCn(:,n)=Dat(1+(n-1)*m:864000+(n-1)*m);
    l= length(Cn);
    fs=1;
    % hann(n)
    window = 0.5 - 0.5*cos(2*pi*linspace(0, 1, l));

    Xcc=fft(Cn.*window');        %calculation of fft points
    Fs=1;                %defining sampling rate
    Pcc2=abs(Xcc).^2/Fs/length(Cn);
    Pcc=Pcc2(1:l/2);    %calculating power spectral density (psd)
    f=0:Fs/(l-1):Fs/2;%defining frequencies
    f=f';               %transposing frequency vector
    
    area(f,Pcc)
    hold on
    plot(f,Pcc,'k','LineWidth',1.5),grid    %plotting frequency vs. psd

    Dn=[f,Pcc];          %Forming a matrix for freq and psd.
    %plot(v(1:end),Xcc(1:end)),grid
    fname=num2str(n);
    save(fname,'Dn','-ascii')
    %store(:,n)=Dn(:);
    %new_data=Dn(1:600000,:);
    %new_data=store(1:2000000,1:n);
end

%-----------------------
% Plotting
ymin = 0;
ymax = 2e7;
xmin = 0.00018;
xmax = .005;

ylim([ymin ymax]);
xlim([xmin xmax]);
xlabel('freq (Hz)');
ylabel('Power Amplitude')
%figure
%plot(X)

% Read the Excel file into a table
T = readtable('normal_modes_data2.xlsx');

Nmodes = 100;
figure(100)
hold on
for ii=1:Nmodes
    label = T.Mode{ii};
    plot([T.PREM_freq(ii) T.PREM_freq(ii)], [ymin ymax],'r--')
    t=text(T.PREM_freq(ii), 1.8*10^7, label, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    t.Rotation=90;
end

% save ('Combined_attnuation.txt','store','-ASCII')
% plot(residual)
