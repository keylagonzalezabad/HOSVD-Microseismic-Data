%Created by Maxim Yatsenko and Milan Brankovic

%THIS CODE CREATES A NUMBER OF SYNTHETIC EVENTS ON WHICH WE TRAIN OUR CODE.

%VARIABLES i1, i2 and i3 DETERMINE NORMAL DISTANCE FROM SOURCE TO RECEIVER
%LINE, VERTICAL POSITION, AND ANGLE OF FRACTURE OPENING.

%FOR EACH EVENT, POSITION AND ANGLE ARE RECORDED IN VARIABLE dm, WHILE
%OTHER PROPERTIES OF CORE TENSOR ARE RECORDED IN VARIABLES rez1, and rez2.

%THE PROPERTIES OF rez1 ARE MOSTLY AFFECTED BY POSITION, WHILE rez2 IS
%MOSTLY AFFECTED BY THE ANGLE

clear all
h=[0 10 20 30 40 50 60 70 80 90];
c=1;
er0=100;

i1u=21;
i2u=21;
i3u=21;

for i1=1:i1u
    for i2=1:i2u
        for i3=1:i3u
            er=er0+rand*1000;
            xp=-2*i1-8;
            zp=i2+14;
            ang=sin(3.14*(i3-1)*3/180);
            rn=1+0.2*rand;
            spr=0.5+0.1*rand;
            ra=1+5*rand;
            
            for i=1:10
                vec=[-xp 0 h(i)-zp];
                r(i)=norm(vec);
                p(i,1:3)=vec/r(i);
                vecn=[0 h(i)-zp xp*ang/((1-ang^2)^0.5)];
                if norm(vecn)==0
                    ps=[0 0 0];
                else
                    ps=cross(p(i,1:3),vecn);
                    ps=ps/norm(ps);
                end
                
                angv=dot(p(i,1:3),[1 0 0]);
                rp=angv*ra;
                rs=rp*spr*sqrt(1-angv^2);
                for t=1:1000
                    trc(t,i,1)=rp*pwf(t/10-r(i)/2,rn)*p(i,1)*(r(i)/10)^(-2)+rs*swf(t/10-r(i),rn)*ps(1)*(r(i)/10)^(-2)+(rand-0.5)/er;
                    trc(t,i,2)=rp*pwf(t/10-r(i)/2,rn)*p(i,2)*(r(i)/10)^(-2)+rs*swf(t/10-r(i),rn)*ps(2)*(r(i)/10)^(-2)+(rand-0.5)/er;
                    trc(t,i,3)=rp*pwf(t/10-r(i)/2,rn)*p(i,3)*(r(i)/10)^(-2)+rs*swf(t/10-r(i),rn)*ps(3)*(r(i)/10)^(-2)+(rand-0.5)/er;
                end
            end
            trcd((i1-1)*i3u*i2u+(i2-1)*i3u+i3,:,:,:)=trc;

%   THE CODE BELOW IS USED FOR THE ANALSYS OF THE SYNTHETIC DATA
%             [U,S]=mlsvd(trc,[10,10,3]);
% 
%             for i=1:10
%                 plt1(i)=sqrt(sum(sum((S(i,:,:).^2))));
%                 plt2(i)=sqrt(sum(sum((S(:,i,:).^2))));
%             end
%             for i=1:3
%                 plt3(i)=sqrt(sum(sum(S(:,:,i).^2)));
%             end
%             p1(i1,i2,i3)=log(gsum(plt1(1:3))/gsum(plt1(8:10)));
%             p2(i1,i2,i3)=log(gsum(plt2(1:3))/gsum(plt2(8:10)));
%             p6(i1,i2,i3)=log(plt3(2)/plt3(3));
%             rez1((i1-1)*i3u*i2u+(i2-1)*i3u+i3)=0.5*p1(i1,i2,i3)+0.5*p2(i1,i2,i3);
%             rez2((i1-1)*i3u*i2u+(i2-1)*i3u+i3)=p6(i1,i2,i3);
%             dm(1,(i1-1)*i3u*i2u+(i2-1)*i3u+i3)=-xp;
%             dm(2,(i1-1)*i3u*i2u+(i2-1)*i3u+i3)=zp;
%             dm(3,(i1-1)*i3u*i2u+(i2-1)*i3u+i3)=ang;
        end
    end
    i1
end

%LINES BELOW ARE JUST USED TO PLOT ONE OF THE DATA SET EXAMPLES
s=size(trcd);
for i=1:s(2)
    for j=1:s(3)
        plt(i,j)=abs(trcd(120,i,j,3))^0.5*sign(trcd(120,i,j,3));
    end
end
surface(plt,'EdgeColor','none')