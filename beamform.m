clear
clc
format long;
v=1;
M=3;
N=1000;%%%%%%%快拍数
f0=16*10^3;%%%%%%%%%%%信号与干扰的频率
f1=16*10^3;
f2=16*10^3;
omiga0=2*pi*f0;%%%%%%%信号与干扰的角频率
omiga1=2*pi*f1;
omiga2=2*pi*f2;
sita0=9.5941; %信号方向
sita1=77.8985; %干扰方向1
sita2= 48.1114; %干扰方向2
for t=1:N           %%%%%%%%%%%%信号
    adt(t)=sin(omiga0*t/(N*f0));
    a1t(t)=sin(omiga1*t/(N*f1));
    a2t(t)=sin(omiga2*t/(N*f2));
end
for i=1:M    %%%%%%%%%%%%信号的导向矢量：线阵的形式
    ad(i,1)=exp(j*(i-1)*pi*sin(sita0));
    a1(i,1)=exp(j*(i-1)*pi*sin(sita1));
    a2(i,1)=exp(j*(i-1)*pi*sin(sita2));
end
R=zeros(M,M);
for t=1:N
 x1=adt(t)*ad;
 x=adt(t)*ad+a1t(t)*a1+a2t(t)*a2; %阵列对信号的完整响应
     R=R+x*x';%信号的协方差矩阵
end
R=R/N;%%%%%%%%%协方差矩阵，所有快拍数的平均
miu=1/(ad'*inv(R)*ad);%%%%%%这个貌似是LMS算法的公式，具体我记不太清，这里是求最优权值，根据这个公式求出，然后加权
w=miu*inv(R)*ad;
%%%%%%形成波束%%%%%%%%%%%%%%%%%%%
for sita=0:pi/100:pi
    for i=1:M
        x_(i,1)=exp(j*(i-1)*pi*sin(sita));
    end
    y(1,v)=w'*x_;%%%%%%%对信号进行加权，消除干扰
    v=v+1;
end
y_max=max(y(:));%%%%%%%%%%%%%%%归一化
y_1=y/y_max;
y_db=20*log(y_1);
 
sita=0:pi/100:pi;
plot(sita,y)
xlabel('sita');
ylabel('天线增益db');
