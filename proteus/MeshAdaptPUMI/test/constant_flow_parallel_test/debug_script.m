%%Debug Script

clear 
clc
close all

A = load('stiffness.csv');
b = load('force.csv');

cond(A)

x=inv(A)*b

%% Debugging ERM.cpp

qpt=[1.381966011250110e-01, 1.381966011250110e-01, 1.381966011250110e-01;
    5.854101966249690e-01, 1.381966011250110e-01, 1.381966011250110e-01;
    1.381966011250110e-01, 5.854101966249690e-01, 1.381966011250110e-01;
    1.381966011250110e-01, 1.381966011250110e-01, 5.854101966249690e-01];

%eID=0
% verts = [0.000000000000000e+00, 3.533333333333336e-02, 2.500000000000000e-02;
%        0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00;
%        0.000000000000000e+00, 4.999999999999999e-02, 0.000000000000000e+00;
%        4.043948768719352e-02, 4.037218712569610e-02, 0.000000000000000e+00];

% verts = [0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00;
%        5.000000000000000e-02, 0.000000000000000e+00, 0.000000000000000e+00;
%        4.043948768719352e-02, 4.037218712569610e-02, 0.000000000000000e+00;
%        3.463448691282167e-02, 0.000000000000000e+00, 2.498723421202763e-02];

%eID=4
%verts = [0.05, 0, 0; 0.0404395, 0.0403722, 0; 0.0346345, 0, 0.0249872; 0.0710784, 0.0312237, 0];

%eID=860 %check 258 as bottom side
verts = [0.0400052, 0.08, 0.025; 0.05, 0.06, 0.05; 0.05, 0.08, 0.05; 0.0400052, 0.08, 0.05];

plot3(verts(:,1),verts(:,2),verts(:,3),'o');    
     
a=verts(1,:);
b=verts(2,:);
c=verts(3,:);
d=verts(4,:);

J = [(b-a)',(c-a)',(d-a)'];
Jdet=det(J);
Vol = Jdet/6; %abs(dot(a-d,cross(b-d,c-d)))/6;

invJ=inv(J);

nu = 1.004e-6;
%s=1,t=1, term1 
% val = 40*40*Vol*nu;
% 
% % Shape function Gradients
% 
% Ni(1,:) = [-1,-1,-1];
% Ni(2,:) = [1,0,0];
% Ni(3,:) = [0,1,0];
% Ni(4,:) = [0,0,1];
% 
% Ngi(1,:) = -sum(invJ);
% Ngi(2,:) = invJ(1,:);
% Ngi(3,:) = invJ(2,:);
% Ngi(4,:) = invJ(3,:);
% 
% 
% %Construct single block stiffness:
% 
% Ablock = zeros(4,4);
% 
% for s=1:4
%     for t=1:4
%         for j=1:3
%             Ablock(s,t) = Ablock(s,t)+Ngi(s,j)*Ngi(t,j);
%         end
%         %Ablock(s,t) = Ablock(s,t)+Ngi(s,1)*Ngi(t,1);
%     end
% end
% 
% Ablock = Ablock*nu*Vol;
% 
% invAblock=inv(Ablock)
% cond(invAblock)
% %Term 2
% 
% %%RHS
% 
% %a(u,v)
% %rhs_a = (N1g_i(3))*nu*20*Vol;

%% Hierarchic "Bubble" Functions

syms r s t u;
u = 1-r-s-t;

% Nbub(1) = 27*r*s*u;
% Nbub(2) = 27*s*t*u;
% Nbub(3) = 27*r*t*u;
% Nbub(4) = 27*r*s*t;

c = -2.44948974278318;
Nbub(1) = c*(u*r);
Nbub(2) = c*(r*s);
Nbub(3) = c*(s*u);
Nbub(4) = c*(u*t);
Nbub(5) = c*(r*t);
Nbub(6) = c*(s*t);

dNbub(:,1) =  diff(Nbub,r);
dNbub(:,2) =  diff(Nbub,s);
dNbub(:,3) =  diff(Nbub,t);

dNgbub = dNbub*invJ;

%Ng(1,:) = -sum(invJ);
%Ng(2,:) = invJ(1,:);
%Ng(3,:) = invJ(2,:);
%Ng(4,:) = invJ(3,:);
%Ng = Ng*invJ;

Ng = [-1 -1 -1; 1 0 0; 0 1 0; 0 0 1];
Ng = Ng*invJ;

a4=0.5854101966249685;
b4=0.1381966011250150;
w4 = 4.166666666666666e-02;
qpt = [a4,b4,b4,b4;b4,a4,b4,b4;b4,b4,a4,b4;b4,b4,b4,a4;];
w =[w4;w4;w4;w4];
nshl = 6;
nsd = 3;
Atest = zeros(nshl*nsd,nshl*nsd);
for k=1:length(qpt)
    dNgqpt = subs(dNgbub,{r,s,t},qpt(k,2:4)); %to meet the definition in PUMI
    Ablock = zeros(nshl,nshl);
    for a=1:nshl
        for b=1:nshl
            for j=1:3
                Ablock(a,b) = Ablock(a,b)+dNgqpt(a,j)*dNgqpt(b,j)*w(k);
            end
        end
    end
    for i=1:nsd
        offset = 1+(i-1)*nshl;
        Atest(offset:offset+nshl-1,offset:offset+nshl-1)=Atest(offset:offset+nshl-1,offset:offset+nshl-1)+Ablock;
    end
    for i=1:nsd
        for j=1:nsd
            for a=1:nshl
                idx1 = a+(i-1)*nshl;
                for b=1:nshl
                    idx2 = b+(j-1)*nshl;
                    Atest(idx1,idx2) = Atest(idx1,idx2)+dNgqpt(a,j)*dNgqpt(b,i)*w(k);
                end
            end
        end
    end
end
Atest= Atest*Jdet*nu;

% invAblock=inv(Ablock)
% cond(invAblock)

%RHS 

% Ly = 0.2;
% Lz = 0.05;
% qpt_g = qpt*verts; %Linear interpolation of points to real space
% figure(1)
% hold on
% plot3(qpt_g(:,1),qpt_g(:,2),qpt_g(:,3),'xr');    
% hold off
%     
% syms x y z;
% mu = 0.0010021928;
% rho = 998.2;
% dpdy=-1;
% u = 0;
% %v = z/Lz;
% v = 0.5/mu*dpdy*(z^2-Lz*z); %Poiseuille Flow
% wz = 0;
% vel = [u,v,wz];
% p = 1 + dpdy*y/Ly;
% 
% rhs = zeros(nshl*nsd,1);
% force = [0,0,0];
% 
% vellin = zeros(4,3);
% for i=1:4
%     vellin(i,:) = subs(vel,{x,y,z},verts(i,:));
% end
% 
% for k=1:length(qpt)
%     Ngqpt = subs(Nbub,{r,s,t},qpt(k,2:4)); %to meet the definition in PUMI
%     dNgqpt = subs(dNgbub,{r,s,t},qpt(k,2:4)); %to meet the definition in PUMI
%     velqpt = qpt*vellin;
%     velgradqpt = vellin'*Ng;
%     %velqpt = subs(vel,{x,y,z},qpt_g(k,1:3));
%     %velgrad = [diff(vel,x);diff(vel,y);diff(vel,z)]';
%     %velgradqpt=subs(velgrad,{x,y,z},qpt_g(k,1:3));
%     p_qpt=subs(p,{x,y,z},qpt_g(k,1:3));
%     
%     for i=1:nsd
%         for a =1:nshl
%             idx = (i-1)*nshl+a;
%             f_k = force(i)*Ngqpt(a);
%             a_k = nu*(velgradqpt+velgradqpt')*(dNgqpt(a,:))';
%             a_k =a_k(i); %only need the ith component
%             b_k = p_qpt/rho*dNgqpt(a,i);
%             c_k = velgradqpt*velqpt'*Ngqpt(a);
%             c_k = c_k(i);
%             rhs(idx) =  rhs(idx) +(f_k-a_k+b_k-c_k)*w(k); 
%         end
%     end
% end
% rhs=rhs*Jdet;
% 
% 
% %Boundary Terms
% 
% w = [1/6];
% 
% adjvert1=[0, 0.0353333, 0.025;
%             0, 0, 0; 
%             0.0404395, 0.0403722, 0; 
%             0.0346345, 0, 0.0249872];
% normal1=[0.499531, -0.500364, 0.707181];
% A1= 0.00101025;
% bqptshp1 =  [2/3,0,1/6;1/6,0,2/3;1/6,0,1/6];
% qptlocal11= [1/6,2/3,1/6,0;1/6,1/6,2/3,0;2/3,1/6,1/6,0];
% qptlocal12= [2/3,0,1/6;1/6,0,2/3;1/6,0,1/6];
% 
% vellinother = zeros(4,3);
% for i=1:4
%     vellinother(i,:) = subs(vel,{x,y,z},adjvert1(i,:));
% end
% 
% qpt_g = qptlocal11*adjvert1; %Linear interpolation of points to real space
% bflux = zeros(nshl*nsd,1);
% for k=1:3
%     Ngqpt = subs(Nbub,{r,s,t},bqptshp1(k,:)); %to meet the definition in PUMI
%     dNgqpt = subs(dNgbub,{r,s,t},qptlocal11(k,2:4)); %to meet the definition in PUMI
%     velqpt = qpt*vellinother;
%     velgradqpt = vellinother'*Ng;
%     p_qpt=subs(p,{x,y,z},qpt_g(k,1:3));
% 
%   for i=1:nsd
%     for a=1:nshl
%       idx = (i-1)*nshl+a;
%       sigma_k = nu*(velgradqpt+velgradqpt')-p_qpt/rho*eye(nsd,nsd);
%       sigma_k = sigma_k*normal1';
%       sigma_k=sigma_k*0.5; %fluxweight
%       bflux(idx) = bflux(idx)+sigma_k(i)*w*Ngqpt(a);
%     end
%   end
% end
% bflux = bflux*2*A1;
% 
% 
% adjvert2=[0.0404395, 0.0403722, 0;
%           0, 0.0353333, 0.025;
%           0, 0.0764907, 0.025;
%           0, 0.05, 0];
% 
% normal2=[0.201152, 0.844895, 0.495672];
% A2= 0.000598292;
% 
% bqptshp2 =  [0,1/6,1/6;0,1/6,2/3;0,2/3,1/6];
% qptlocal21= [2/3,0,1/6;1/6,0,1/6;1/6,0,2/3];
% qptlocal22= [0,1/6,1/6;0,1/6,2/3;0,2/3,1/6];

