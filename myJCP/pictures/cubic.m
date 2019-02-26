clc,clear,clf;
%theta(0,pi)%phi(0,2*pi)
alpha=0.8;theta=0.354*pi;phi=1.35*pi;
nx=sin(theta)*cos(phi);
ny=sin(theta)*sin(phi);
nz=cos(phi);
n0=[nx ny nz];
n=0.5*n0;
d=1;
a=d;b=0.9*d;c=1.1*d;
center=[a/2 b/2 c/2];
p =[0 0 0;
 a 0 0;
 a b 0;
 0 b 0;
 0 0 c;
 a 0 c;
 a b c;
 0 b c];
points=[1 2 3 4 5 6 7 8];
edges=[1 2 3 4 5 6 7 8 9 10 11 12];
e=[1 2;
   2 3;
   3 4;
   4 1;
   1 5;
   5 8;
   8 4;
   8 7;
   7 3;
   5 6;
   6 7;
   6 2;
   7 3];
% x=[0 a a 0 0 0 0 0 0 0 0 a a a a a a a 0];
% y=[0 0 0 0 0 0 b b 0 b b b b 0 0 b b b b];
% z=[0 0 c c 0 c c 0 0 0 c c 0 0 c c 0 0 0];
figure(1)
set(0,'defaultfigurecolor','w'),
for i=1:12
    x=[p(e(i,1),1) p(e(i,2),1)];
    y=[p(e(i,1),2) p(e(i,2),2)];
    z=[p(e(i,1),3) p(e(i,2),3)];
    plot3(x,y,z,'.-b','LineWidth',2);hold on,
end
axis off,
hold on,
% x0=[center(1) center(1)+n(1)];
% y0=[center(2) center(2)+n(2)];
% z0=[center(3) center(3)+n(3)];
% plot3(x0,y0,z0,'.-k'),
%quiver3(center(1),center(2),center(3),n(1),n(2),n(3),'k','LineWidth',2); 
cc=[center;center;center;center;center;center;center;center];
pn=p-cc;
dd=zeros(1,8);
for i=1:8
    dd(i)=dot(pn(i,:),n0);
end
PII=[];
%  dmax=max(dd);
%  dmin=min(dd);
[dds,index]=sort(dd); 
V=zeros(1,8);
V(8)=1;
for k=2:7
%k=4;
f0=dds(k);
xc=nx*f0+center(1);
yc=ny*f0+center(2);
zc=nz*f0+center(3);
PII=[];
PIJ=[];
%xc,yc,zc,
for i=1:12
    j=(dd(e(i,1))-f0)*(dd(e(i,2))-f0);
    if(j<0)
        f1 = dd(e(i,1));
        f2 = dd(e(i,2));
        x1 = p(e(i,1),1);
        x2 = p(e(i,2),1);
        x0 = x1 - (f1-f0)*((x1-x2)/(f1-f2));
        y1 = p(e(i,1),2);
        y2 = p(e(i,2),2);
        y0 = y1 - (f1-f0)*((y1-y2)/(f1-f2));
        z1 = p(e(i,1),3);
        z2 = p(e(i,2),3);
        z0 = z1 - (f1-f0)*((z1-z2)/(f1-f2));
        p0=[x0 y0 z0];
        PII=[PII;p0];
    end
end
%dot(PII(1,:)-PII(4,:),n),
PII=[PII;p(index(k),:)];
s=size(PII);
pcenter=sum(PII)./s(1);
%PII([3 4],:)=PII([4 3],:);
%plot3(PII(:,1),PII(:,2),PII(:,3),'.k'),
%plot3([PII(:,1);pcenter(1);xc],[PII(:,2);pcenter(2);yc],[PII(:,3);pcenter(3);zc],'.k'),
PIJ=PII(1:s(1),:);
PIJ(:,4)=zeros(s(1),1);
test1=PII(1,:)-pcenter;
for ii=2:s(1)
    figure(1),
    test=PII(ii,:)-pcenter;
    ssin=dot(cross(test,test1),n0);
    ccos=dot(test,test1);
    if ccos==0
%     ttan=ssin/(ccos+1e-12);
        if ssin >0
            Theta=0.5*pi;
            PIJ(ii,4)= Theta;
        else
            Theta=-0.5*pi;
            PIJ(ii,4)= Theta;
        end
    else
        ttan=ssin/ccos;
        Theta=atan(ttan);
        if ttan>0
            if ssin>0
                PIJ(ii,4)= Theta;
            else
                PIJ(ii,4)= Theta - pi;
            end
        else
            if ssin>0
                PIJ(ii,4)= Theta+pi;
            else
                PIJ(ii,4)= Theta;
            end
        end
    end
end
PIJ=sortrows(PIJ,4);
PII(:,[1 2 3])=PIJ(:,[1 2 3]);
%plot3([PII(:,1);pcenter(1);xc],[PII(:,2);pcenter(2);yc],[PII(:,3);pcenter(3);zc],'.-k'),
%plot3(PII(:,1),PII(:,2),PII(:,3),'.-k'),
fill3(PII(:,1),PII(:,2),PII(:,3),'r'),
PIK=PII;
for kk=1:8
    if dd(kk)<f0
        PIK=[PIK;p(kk,:)];
    end
end
dt=delaunayTriangulation(PIK);
[ch V(k)] = convexHull(dt);
figure(k),
set(0,'defaultfigurecolor','w'),
for i=1:12
    x=[p(e(i,1),1) p(e(i,2),1)];
    y=[p(e(i,1),2) p(e(i,2),2)];
    z=[p(e(i,1),3) p(e(i,2),3)];
    plot3(x,y,z,'.-b','LineWidth',2);hold on,
end
axis off,
trisurf(ch, dt.Points(:,1),dt.Points(:,2),dt.Points(:,3), 'FaceColor', 'cyan'),
end
figure(8),
plot(dds,V,'.-k'),