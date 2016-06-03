fprintf('此程序由hone编写，谢谢使用！！\n');
fprintf('\n');
y=input('请输入要进行的实验（比如：1）');
if(y<4)
    switch(y)
    %选择实验（1·6）
       case 1
            k1=input('请输入要进行的实验内容（比如：1）');        
            if k1==1
                %1
                t=-2:0.01:2;
                y=rectpuls(t);
                plot(t,y);title('方波');%画出方波波形，并添加标题'方波'；
        k2=input('请输入要进行的实验内容（比如：1）');
              elseif k1==2
                %2
                t=-2*pi:0.03:2*pi;
                y2=sawtooth(t,0.5);
                plot(t,y2);
                xlabel('幅值');ylabel('时间');
                title('三角波');
    k3=input('请输入要进行的实验内容（比如：1）');
              elseif k1==3
                %31
                t=0:0.0003:0.2;
                y=square(2*pi*30*t);
                plot(t,y);
                xlabel('幅值');ylabel('时间');
                title('周期方波信号');
    k4=input('请输入要进行的实验内容（比如：1）');
            elseif k1==4
                %4
                t=0:0.0003:0.2;
                y=square(2*pi*30*t);
                plot(t,y);
                xlabel('幅值');ylabel('时间');
                title('周期方波信号');
            end
         case 2
            k2=input('请输入要进行的实验内容（比如：1）');
            if k2==1
                %实验二
                %1
                t1=1:0.01:2 
                f1=ones(size(t1)).*(t1>1);% 表示一个高度为1的门函数，时间从t=1到t=2
                t2=2:0.01:3;
                f2=ones(size(t2)).*(t2>2);% 表示另一个高度为1的门函数，时间从t=2到t=)
                c=conv(f1,f2);% 表示卷积
                t3=3:0.01:5;
                subplot(3,1,1),plot(t1,f1);
                subplot(3,1,2),plot(t2,f2);
                subplot(3,1,3),plot(t3,c);t1=1:0.01:2;



            elseif k2==2
                %2
                t1=0:0.01:1;
                f1=t1.*(t1>0);
                t2=-1:0.01:2;
                f2=t2.*exp(-t2).*(t2>0)+exp(t2).*(t2<0);
                c=conv(f1,f2);
                t3=-1:0.01:3;
                subplot(3,1,1),plot(t1,f1);title('f1(t)');
                subplot(3,1,2),plot(t2,f2);title('f2(t)');
                subplot(3,1,3),plot(t3,c);title('c(t)');


            elseif k2==3
                %3
                x=[2 -1 3];h=[1 2 2 3];    %信号x[n]和h[n]
                y=conv(x,h);                %计算卷积y[n]
                nx=-1:1;nh=-2:1;           % x[n]和y[n]的序号
                ns=nx(1)+nh(1);            % y[n]的开始序号
                ne=nx(length(nx))+nh(length(nh)); % y[n]的结果序号
                ny=ns:ne;                  % y[n]的序号
                figure;
                subplot(311);
                stem(nx,x);  title('x[n]'); %绘制信号x，离散波形的画图用stem
                subplot(312);
                stem(nh,h); title('h[n]');    %绘制信号h
                subplot(313);
                stem(ny,y); title('y[n]');     %绘制信号y
            end
         case 3
            %实验三
            %1
            k3=input('请输入要进行的实验内容（比如：1）');
            if k3==1
                %function fangbofj(N)
                %
                N=input('谐波次数：N=');
                N=double(N);
                t=-1:0.01:1;
                w=2*pi;
                y=square((t+1/4)*2*pi,50);

                %循环画前N次谐波函数体
                for(k=1:2:1+N)%循环画前N次谐波
                    n=1:2:k;
                    b=4./(pi*n);
                    x=b*sin(w*n'*(t+1/4));
                    figure;
                    plot(t,y);
                    hold on;
                    plot(t,x);
                    hold off;
                    xlabel('t');ylabel('部分波形');
                    axis([-1 1 -1.5 1.5]);
                    grid on;
                title(['最大谐波数',num2str(k)]);
                 end
            elseif k3==2
                %2
                %function fangbofj1(N)
                %
                N=input('谐波次数：N=');
                N=double(N);
                t=-1:0.01:1;
                w=2*pi;
                y=square((t+1/4)*2*pi,50);

                    n=1:2:N;
                    b=4./(pi*n);
                    x=b*sin(w*n'*(t+1/4));
                    figure;
                    plot(t,y);
                    hold on;
                    plot(t,x);
                    hold off;
                    xlabel('t');ylabel('部分波形');
                    axis([-1 1 -1.5 1.5]);
                    grid on;
                    title(['最大谐波数=',num2str(N)]);
            end
        case 4
            %实验四

            %function tiaozhi(n)
            n=input('选择调制方式（1、正弦波sin 2、方波square 3、三角波sawtooth）n=');
            n=double(n);
            t=0:pi/100:50*pi;  

            switch(n)
                case 1
                    Carrier=sin(t);
                case 2
                    Carrier=square(0.5*t,50);
                case 3
                    Carrier=sawtooth(t,0.5); 
                otherwise
                    disp('请输入正确的数字！');
            end

            Mod_Sig=sin(t/20);
            Dsb_am=Carrier.*(1+Mod_Sig);
            figure;
            subplot(221);
            plot(t,Carrier);
            mm=max(Carrier);
            title('载波');
            %axis([0,20*pi,-1.2*Carrier,1.2*Carrier);
            subplot(222);
            plot(t,Mod_Sig);title('调制信号');
            subplot(212);
            plot(t,Dsb_am);
            title('调制后的信号');
            xlabel('time');
            ylabel('voltage');
            grid on;
            %legend('carrier','baseband','modulated signal') 
        end
elseif y>4
    if y==5
           
                        %通过输入零极点矩阵画出系统零极图，以及时域，频域波形；
                    %B为系统函数分子系数；A为分母系数；
                    %eg：B=[2]     A=[1 0 4],
                    %系统时域波形即为正弦信号；
                    %注意：输入系数时加方括号，因A、B均为矩阵；

                    %B=[2,30];A=[1,10,50];
                    B=input('输入分子系数（注意为矩阵[2,30]）：');
                    A=input('输入分母系数（注意为矩阵[1,10,50]）：');
                    %{
                    也可以替换为下列程序段实现动态输入变量：
                    B=input('输入分子系数（注意为矩阵）：');
                    A=input('输入分母系数（注意为矩阵）：');
                    %}
                    sys1=tf(B,A);
                    figure;
                    subplot(2,2,1);
                    pzmap(sys1);


                    subplot(2,2,2);
                    grid on;
                    impulse(B,A);


                    w=-8*pi:0.01:8*pi;
                    H=freqs(B,A,w);
                    subplot(2,2,3);
                    plot(w,abs(H));grid on;
                    xlabel('\omega(rad/s)');ylabel('|H(\omega)|');
                    title('幅频特性');

                    subplot(2,2,4);
                    plot(w,angle(H));grid on;
                    xlabel('\omega(rad/s)'),ylabel('\phi(\omega)');
                    title('H(s)的相频特性');
            
    elseif y==6
                f1=input('输入信号频率(比如：1000)：f1=');%1000;
                f2=input('输入信号频率（比如：2500）：f2=');%2500;
                f3=input('输入信号频率（比如：3800）：f3=');%3800;
                p=input('输入归一化截止频率[0~1]之间（比如：0.5）：w=');%3800;
                %}
                %f1=double(f1);f2=double(f2);f3=double(f3);p=double(p);
                s=p+0.05;
                f=[f1 f2 f3];
                ff=max(f);

                N=1024;      %这个是你举得信号的点数
                fs=8000;       %这个是抽样频率，要高于信号中最高频率的2倍
                t=(0:N-1)/fs;     %信号时域横轴向量
                f=fs*(0:1023)/2048;
                f=(0:N-1)*fs/2048;    %信号频域横轴向量，不过待会要减半处理一下
                A1=0.05;A2=0.10;A3=0.08;
                x=A1*sin(2*pi*f1*t)+A2*sin(2*pi*f2*t)+A3*sin(2*pi*f3*t);   %想分析的函数

                wp=p*pi;
                ws=s*pi;
                Rp=1;
                Rs=15;
                Fs=8000;
                Ts=1/Fs;
                wp1=2/Ts*tan(wp/2);                 %将模拟指标转换成数字指标
                ws1=2/Ts*tan(ws/2); 
                [N,Wn]=buttord(wp1,ws1,Rp,Rs,'s');  %选择滤波器的最小阶数
                [Z,P,K]=buttap(N);                  %创建butterworth模拟滤波器
                [Bap,Aap]=zp2tf(Z,P,K);
                [b,a]=lp2lp(Bap,Aap,Wn);   
                [bz,az]=bilinear(b,a,Fs);           %用双线性变换法实现模拟滤波器到数字滤波器的转换
                [H,W]=freqz(bz,az);                 %绘制频率响应曲线
                figure(1)
                plot(W*Fs/(2*pi),abs(H))
                grid on;axis tight;
                xlabel('频率（Hz）')
                ylabel('频率响应')
                title('Butterworth')

                f1=filter(bz,az,x);
                figure(2)
                subplot(2,1,1)
                plot(t,x)                          %画出滤波前的时域图
                grid on;axis tight;
                title('滤波前的时域波形');
                subplot(2,1,2)
                plot(t,f1);                         %画出滤波后的时域图
                grid on;axis tight;
                title('滤波后的时域波形');
                %}

                y3=fft(f1,2048);
                figure(3)
                y2=fft(x,2048);
                subplot(2,1,1);
                plot(f,abs(y2(1:1024)));             %画出滤波前的频谱图
                grid on;axis tight;
                title('滤波前的频谱')
                xlabel('Hz');
                ylabel('幅度');
                subplot(2,1,2)
                plot(f,abs(y3(1:1024)));          %画出滤波后的频谱图
                grid on;axis tight;
                title('滤波后的频谱')
                xlabel('Hz');
                ylabel('幅度');
                warning off；


    end
end
