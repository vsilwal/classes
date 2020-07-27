%
% lab_svdgeom.m
% Carl Tape, GEOS 627, Inverse Problems and Parameter Estimation
%
% Misc concepts associated with the singular value decomposition
%

clc, clear, close all

bprint = false;
pdir = "./";
bdagger = true;
bdaggersol = true;
  
% plot settings
f = 0.5;    % radial shift for text label
fsize = 10;    
axmax = 3;
ax0 = axmax*[-1 1 -1 1];
mtest = [-1 -1]';

btest = true;
if bdaggersol
    ftag = '_sol';
    be1e2 = true;
    bu1u2 = true;
else
    ftag = '';
    be1e2 = false;
    bu1u2 = false;
end

G = [ 2  1 ; -1  1]
% G = [ 0 -2 ;  2  0]
% G = [ 2 -2 ;  2 -1]
% G = [ 2  0 ;  2  0]
% G = [ 1 -1 ; -2  0]       % default
% G = [ 0  1 ; -1  1]
% G = [-2  0 ;  0  1]
% G = [ 1 -1 ;  1  1]
% G = randi([-2 2],2,2)

% singular value decomposition
[U,S,V] = svd(G)
Gcheck = U*S*V'

swapm = [0 1 ; 1 0 ];
rotm = @(thetad) ( [ cosd(thetad) -sind(thetad) ; sind(thetad) cosd(thetad) ] );

% singular values
s1 = S(1,1);
s2 = S(2,2);
% basis vectors
v1 = V(:,1);
v2 = V(:,2);
u1 = U(:,1);
u2 = U(:,2);

% unit circle for plotting
R = 1;
n = 100;
t = linspace(0,2*pi,n);
x = R*cos(t);
y = R*sin(t);
Cx = [x ; y];

% unit circles representing basis U and basis V
Ux = Cx;
Vx = Cx;

% standard basis
e1 = [1 0]';
e2 = [0 1]';

figure; nr=2; nc=2;

subplot(nr,nc,1); hold on; grid on;
set(gca,'xtick',[-axmax:axmax],'ytick',[-axmax:axmax]);
plot(Vx(1,:),Vx(2,:),'b-'); axis equal, axis(ax0);
plot([0 v1(1)],[0 v1(2)],'b',[0 v2(1)],[0 v2(2)],'b','linewidth',2);
if be1e2
    plot([0 e1(1)],[0 e1(2)],'c',[0 e2(1)],[0 e2(2)],'c','linewidth',1);
    %plot_text(e1,'e_1',f,fsize);
    %plot_text(e2,'e_2',f,fsize); 
end

GVx = G*Vx;     % ellipse for plotting
su1 = s1*u1;
su2 = s2*u2;
Ge1 = G*e1;
Ge2 = G*e2;
dtest = G*mtest;

subplot(nr,nc,2); hold on; grid on;
set(gca,'xtick',[-axmax:axmax],'ytick',[-axmax:axmax]);
plot(GVx(1,:),GVx(2,:),'r'); axis equal, axis(ax0);
plot(Ux(1,:),Ux(2,:),'r--');
plot([0 su1(1)],[0 su1(2)],'r',[0 su2(1)],[0 su2(2)],'r','linewidth',2);
if be1e2
    plot([0 Ge1(1)],[0 Ge1(2)],'c',[0 Ge2(1)],[0 Ge2(2)],'c','linewidth',1);
    %plot_text(Ge1,'G e_1',f,fsize);
    %plot_text(Ge2,'G e_2',f,fsize); 
end

if bu1u2
    [Veig,Deig] = eig(G)
    if isreal(Deig)
        GVeig = G*Veig;
        subplot(nr,nc,1);
        plot([0 Veig(1,1)],[0 Veig(2,1)],'k',...
             [0 Veig(1,2)],[0 Veig(2,2)],'k','linewidth',1);
        plot_text(Veig(:,1),'h_1',f,fsize);
        plot_text(Veig(:,2),'h_2',f,fsize); 
        subplot(nr,nc,2);
        plot([0 GVeig(1,1)],[0 GVeig(2,1)],'k',...
             [0 GVeig(1,2)],[0 GVeig(2,2)],'k','linewidth',1);
        plot_text(GVeig(:,1),'G h_1',f,fsize);
        plot_text(GVeig(:,2),'G h_2',f,fsize); 
    end
end

if btest
    subplot(nr,nc,1); plot(mtest(1),mtest(2),...
        'kp','markersize',12,'markerfacecolor','c','linewidth',1);
    %plot_text(mtest,'m',f,fsize);
    subplot(nr,nc,2); plot(dtest(1),dtest(2),...
                'kp','markersize',12,'markerfacecolor','r','linewidth',1);
end

% text labels
subplot(nr,nc,1); %plot_text(v1,'v_1',f,fsize); plot_text(v2,'v_2',f,fsize);
subplot(nr,nc,2); %plot_text(su1,'s u_1',f,fsize); plot_text(su2,'s u_2',f,fsize);
% print figure
if bprint, print(gcf,'-depsc',sprintf('%ssvd_2D_%i%s',pdir,kk,ftag)); end

%------------------------

if bdagger
    subplot(nr,nc,3); hold on; grid on;
    set(gca,'xtick',[-axmax:axmax],'ytick',[-axmax:axmax]);
    plot(Ux(1,:),Ux(2,:),'r-'); axis equal, axis(ax0);
    subplot(nr,nc,4); hold on; grid on;
    set(gca,'xtick',[-axmax:axmax],'ytick',[-axmax:axmax]);
    plot(Vx(1,:),Vx(2,:),'b--'); axis equal, axis(ax0);

    if bdaggersol
        Gdagger = XXX;
        inv(G)

        subplot(nr,nc,2);
        title(sprintf('G = [%i %i ; %i %i]',G),'fontsize',9);
        
        subplot(nr,nc,3);
        plot([0 u1(1)],[0 u1(2)],'r',[0 u2(1)],[0 u2(2)],'r','linewidth',2);
        if be1e2
            plot([0 e1(1)],[0 e1(2)],'m',[0 e2(1)],[0 e2(2)],'m','linewidth',1);
            %plot_text(e1,'e_1',f,fsize); plot_text(e2,'e_2',f,fsize); 
        end

        % ellipse for plotting
        GdUx = Gdagger*Ux;
        sv1 = (1/s1)*v1;
        sv2 = (1/s2)*v2;
        Gde1 = Gdagger*e1;
        Gde2 = Gdagger*e2;
        mfinal = XXX

        subplot(nr,nc,4);
        plot(GdUx(1,:),GdUx(2,:),'b'); axis equal, axis(ax0);
        plot([0 sv1(1)],[0 sv1(2)],'b',[0 sv2(1)],[0 sv2(2)],'b','linewidth',2);
        if be1e2
            plot([0 Gde1(1)],[0 Gde1(2)],'m',[0 Gde2(1)],[0 Gde2(2)],'m','linewidth',1);
            %plot_text(Gde1,'Gd e_1',f,fsize); plot_text(Gde2,'Gd e_2',f,fsize); 
        end
        title(sprintf('Gd = [%.2f %.2f ; %.2f %.2f]',Gdagger),'fontsize',9);
        
        if btest
            subplot(nr,nc,3); plot(dtest(1),dtest(2),...
                'kp','markersize',12,'markerfacecolor','r','linewidth',1);
           % plot_text(dtest,'d',f,fsize);
            subplot(nr,nc,4); plot(mtest(1),mtest(2),...
                'kp','markersize',12,'markerfacecolor','c','linewidth',1);
           % plot_text(mfinal,'m',f,fsize);
        end
        
        % text labels
        subplot(nr,nc,3);% plot_text(sv1,'u_1',f,fsize); plot_text(sv2,'u_2',f,fsize); 
        subplot(nr,nc,4);% plot_text(sv1,'(1/s_1) v_1',f,fsize); plot_text(sv2,'(1/s_2) v_2',f,fsize); 
    end

    % print figure
    if bprint, print(gcf,'-depsc',sprintf('%ssvd_2D_both_%i%s',pdir,kk,ftag)); end
end

    % for kk
