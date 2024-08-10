close all
clear all

% initial vector
origin = [ 0; 0];
vec = [3 ; 5];
figure(1); subplot(2,2,1)
plot([origin(1) vec(1)],[origin(2) vec(2)],'r')

% G input matrix 
G = [-11 4; 11 9];

% prepare all vectors (0 to 360)
for theta=0:5:360
    % R is the rotation matrix
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    Rv = R*vec;
    figure(1); subplot(2,2,1)
    axis equal
    hold on
    plot([origin(1) Rv(1)],[origin(2) Rv(2)],'b')
    
    % Multiple vector v with G (Gv)
    Gv = G*Rv;   % Rv is rotated vector
    subplot(2,2,2)
    hold on
    plot([origin(1) Gv(1)],[origin(2) Gv(2)],'r')
end

% obtain eigenvalues and eigenvectors
[V,D] = eig(G);
v1 = V(:,1); v2 = V(:,2);

% plot
figure(1)
subplot(2,2,2)
plot([origin(1) V(1)],[origin(2) V(2)],'k')
plot([origin(1) V(3)],[origin(2) V(4)],'k')

% plot eigenvector
subplot(2,2,3)
hold on
plot([origin(1) v1(1)],[origin(2) v1(2)],'r','LineWidth',4)
plot([origin(1) v2(1)],[origin(2) v2(2)],'r','LineWidth',4)
axis equal

% plot eigenvector multiplied with G
subplot(2,2,3)
Gv1 = G*v1; Gv2 = G*v2;
hold on
plot([origin(1) Gv1(1)],[origin(2) Gv1(2)],'k','LineWidth',2)
plot([origin(1) Gv2(1)],[origin(2) Gv2(2)],'k','LineWidth',2)
axis equal
