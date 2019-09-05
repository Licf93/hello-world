%***********************************************
%异常向下延拓_2D
%***********************************************
clc
clear all
close all
fclose('all');
tic%保存当前时间

B=load('Paralines_Data_2D.mat');
Y=B.y;
d=B.d;
Bx=B.Bx;
By=B.By;
Bz=B.Bz;

H=1.5;
N=1000;
s=1;

Ny=length(Y);

Ydelt=(Y(Ny)-Y(1))/(Ny-1);%测点间距

n1=0:Ny-1;

n=zeros(Ny,1);

NyOddOrEvent=(mod(Ny,2)==0);

switch NyOddOrEvent
    case 0
        for i=1:Ny
            if n1(i)<=fix(Ny/2)
                n(i)=n1(i);
            else
                n(i)=n1(i)-Ny;
            end
        end
        
    case 1
        for i=1:Ny
            if n1(i)<=Ny/2-1
                n(i)=n1(i);
            else
                n(i)=n1(i)-Ny;
            end
        end
end

Ly=Ny*Ydelt;

Ky=n/(Ly);

A=zeros(1,Ny);

for iy=1:Ny
    A(iy)=exp(-2*pi*H*sqrt(Ky(iy)^2));
end

BX=Bx;
for qq=1:1:N %迭代次数
    BF=fft2(BX);
    BF=BF.*A;
    B2=ifft2(BF);
    err=Bx-B2;
    if err<1e-4
        break;
    end
    BX=BX+s*err;
end

BY=By;
for qq=1:1:N %迭代次数
    BF=fft2(BY);
    BF=BF.*A;
    B2=ifft2(BF);
    err1=By-B2;
    if err1<1e-4
        break;
    end
    BY=BY+s*err1;
end

BZ=Bz;
for qq=1:1:N %迭代次数
    BF=fft2(BZ);
    BF=BF.*A;
    B2=ifft2(BF);
    err2=Bz-B2;
    if err2<1e-4
        break;
    end
    BZ=BZ+s*err2;
end

save('E:\学术研究\我的成果\科技论文\科技论文\Estimation of the location and depth information of underground parallel pipeline using magnetic tilt angle\matlab\Paralines_Data_downward_2D_0.8.mat','BX','BY','BZ','Bx','By','Bz','Y','d')