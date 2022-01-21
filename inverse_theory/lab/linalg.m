close all
clear all
clc

% Multiplication matrix
% Each column represents seperate vector
G = [2 1 5; 
    -1 3 1; 
    5 -1 4];
% scaling factor for each vector
% Hence, number of scaling factors must be equal to number of vector in G
m = [3 2 7]';

% Extract vectors
u = G(:,1);
v = G(:,2);
w = G(:,3);

% Peform matrix multiplication for comparison
p = G*m;
orig = zeros(3,3);

% Plot
% figure
% hold on
% quiver3(0,0,0,u(1),u(2),u(3),'b','Linewidth',2), text(u(1),u(2),u(3),'u')
% quiver3(0,0,0,v(1),v(2),v(3),'b','Linewidth',2), text(v(1),v(2),v(3),'v')
% quiver3(0,0,0,w(1),w(2),w(3),'b','Linewidth',2), text(w(1),w(2),w(3),'w')
% quiver3(0,0,0,p(1),p(2),p(3),'r','Linewidth',5)
% view(3); grid on

%-----------------
% Alternate approach
% scale each vector by corresponding model parameter
u1 = m(1) * u;
v1 = m(2) * v;
w1 = m(3) * w;

% Sum of scaled vectors
p1= u1+ v1 + w1;

% Plot
figure; hold on 
quiver3(0,0,0,u1(1),u1(2),u1(3),'g','Linewidth',2), text(u1(1),u1(2),u1(3),'u1')
quiver3(0,0,0,v1(1),v1(2),v1(3),'g','Linewidth',2), text(v1(1),v1(2),v1(3),'v1')
quiver3(0,0,0,w1(1),w1(2),w1(3),'g','Linewidth',2), text(w1(1),w1(2),w1(3),'w1')
quiver3(0,0,0,p1(1),p1(2),p1(3),'r','Linewidth',5)

quiver3(0,0,0,u(1),u(2),u(3),'b','Linewidth',2), text(u(1),u(2),u(3),'u')
quiver3(0,0,0,v(1),v(2),v(3),'b','Linewidth',2), text(v(1),v(2),v(3),'v')
quiver3(0,0,0,w(1),w(2),w(3),'b','Linewidth',2), text(w(1),w(2),w(3),'w')
quiver3(0,0,0,p(1),p(2),p(3),'k','Linewidth',2)
view(3); grid on


% OUTCOME
% The resultant vector in black (p1) matches with previous result in red
% (p)

%========================================================================
% Eigenvalues and eigenvectors

[V,D] = eig(G)

% Get eigenvectors
v1= V(:,1);
v2= V(:,2);
v3= V(:,3);

% Check by multiplyin eigenvectors to the matrix
vg1 = G*v1;
vg2 = G*v2;
vg3 = G*v3;

% PLot
figure; hold on
quiver3(0,0,0,vg1(1),vg1(2),vg1(3),'r','Linewidth',2), text(vg1(1),vg1(2),vg1(3),'vg1')
quiver3(0,0,0,vg2(1),vg2(2),vg2(3),'r','Linewidth',2), text(vg2(1),vg2(2),vg2(3),'vg2')
quiver3(0,0,0,vg3(1),vg3(2),vg3(3),'r','Linewidth',2), text(vg3(1),vg3(2),vg3(3),'vg3')
quiver3(0,0,0,v1(1),v1(2),v1(3),'b','Linewidth',2), text(v1(1),v1(2),v1(3),'v1')
quiver3(0,0,0,v2(1),v2(2),v2(3),'b','Linewidth',2), text(v2(1),v2(2),v2(3),'v2')
quiver3(0,0,0,v3(1),v3(2),v3(3),'b','Linewidth',2), text(v3(1),v3(2),v3(3),'v3')
view(3); grid on

%====================================================
% Problem 1: Check if eigenvectors of G form orthogonal basis

% Problem 2: Find eigenvectors of modified matrix A= (G'G). Check if eigenvectors of A form orthogonal basis
