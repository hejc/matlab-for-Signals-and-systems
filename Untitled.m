fprintf('�˳�����HeJuncheng��д��ллʹ�ã���\n');
fprintf('\n');
y=input('������Ҫ���е�ʵ�飨���磺1��');
if(y<4)
    switch(y)
    %ѡ��ʵ�飨1��6��
       case 1
            k1=input('������Ҫ���е�ʵ�����ݣ����磺1��');        
            if k1==1
                %1
                t=-2:0.01:2;
                y=rectpuls(t);
                plot(t,y);title('����');%�����������Σ�����ӱ���'����'��
        k2=input('������Ҫ���е�ʵ�����ݣ����磺1��');
              elseif k1==2
                %2
                t=-2*pi:0.03:2*pi;
                y2=sawtooth(t,0.5);
                plot(t,y2);
                xlabel('��ֵ');ylabel('ʱ��');
                title('���ǲ�');
    k3=input('������Ҫ���е�ʵ�����ݣ����磺1��');
              elseif k1==3
                %31
                t=0:0.0003:0.2;
                y=square(2*pi*30*t);
                plot(t,y);
                xlabel('��ֵ');ylabel('ʱ��');
                title('���ڷ����ź�');
    k4=input('������Ҫ���е�ʵ�����ݣ����磺1��');
            elseif k1==4
                %4
                t=0:0.0003:0.2;
                y=square(2*pi*30*t);
                plot(t,y);
                xlabel('��ֵ');ylabel('ʱ��');
                title('���ڷ����ź�');
            end
         case 2
            k2=input('������Ҫ���е�ʵ�����ݣ����磺1��');
            if k2==1
                %ʵ���
                %1
                t1=1:0.01:2 
                f1=ones(size(t1)).*(t1>1);% ��ʾһ���߶�Ϊ1���ź�����ʱ���t=1��t=2
                t2=2:0.01:3;
                f2=ones(size(t2)).*(t2>2);% ��ʾ��һ���߶�Ϊ1���ź�����ʱ���t=2��t=)
                c=conv(f1,f2);% ��ʾ���
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
                x=[2 -1 3];h=[1 2 2 3];    %�ź�x[n]��h[n]
                y=conv(x,h);                %������y[n]
                nx=-1:1;nh=-2:1;           % x[n]��y[n]�����
                ns=nx(1)+nh(1);            % y[n]�Ŀ�ʼ���
                ne=nx(length(nx))+nh(length(nh)); % y[n]�Ľ�����
                ny=ns:ne;                  % y[n]�����
                figure;
                subplot(311);
                stem(nx,x);  title('x[n]'); %�����ź�x����ɢ���εĻ�ͼ��stem
                subplot(312);
                stem(nh,h); title('h[n]');    %�����ź�h
                subplot(313);
                stem(ny,y); title('y[n]');     %�����ź�y
            end
         case 3
            %ʵ����
            %1
            k3=input('������Ҫ���е�ʵ�����ݣ����磺1��');
            if k3==1
                %function fangbofj(N)
                %
                N=input('г��������N=');
                N=double(N);
                t=-1:0.01:1;
                w=2*pi;
                y=square((t+1/4)*2*pi,50);

                %ѭ����ǰN��г��������
                for(k=1:2:1+N)%ѭ����ǰN��г��
                    n=1:2:k;
                    b=4./(pi*n);
                    x=b*sin(w*n'*(t+1/4));
                    figure;
                    plot(t,y);
                    hold on;
                    plot(t,x);
                    hold off;
                    xlabel('t');ylabel('���ֲ���');
                    axis([-1 1 -1.5 1.5]);
                    grid on;
                title(['���г����',num2str(k)]);
                 end
            elseif k3==2
                %2
                %function fangbofj1(N)
                %
                N=input('г��������N=');
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
                    xlabel('t');ylabel('���ֲ���');
                    axis([-1 1 -1.5 1.5]);
                    grid on;
                    title(['���г����=',num2str(N)]);
            end
        case 4
            %ʵ����

            %function tiaozhi(n)
            n=input('ѡ����Ʒ�ʽ��1�����Ҳ�sin 2������square 3�����ǲ�sawtooth��n=');
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
                    disp('��������ȷ�����֣�');
            end

            Mod_Sig=sin(t/20);
            Dsb_am=Carrier.*(1+Mod_Sig);
            figure;
            subplot(221);
            plot(t,Carrier);
            mm=max(Carrier);
            title('�ز�');
            %axis([0,20*pi,-1.2*Carrier,1.2*Carrier);
            subplot(222);
            plot(t,Mod_Sig);title('�����ź�');
            subplot(212);
            plot(t,Dsb_am);
            title('���ƺ���ź�');
            xlabel('time');
            ylabel('voltage');
            grid on;
            %legend('carrier','baseband','modulated signal') 
        end
elseif y>4
    if y==5
           
                        %ͨ�������㼫����󻭳�ϵͳ�㼫ͼ���Լ�ʱ��Ƶ���Σ�
                    %BΪϵͳ��������ϵ����AΪ��ĸϵ����
                    %eg��B=[2]     A=[1 0 4],
                    %ϵͳʱ���μ�Ϊ�����źţ�
                    %ע�⣺����ϵ��ʱ�ӷ����ţ���A��B��Ϊ����

                    %B=[2,30];A=[1,10,50];
                    B=input('�������ϵ����ע��Ϊ����[2,30]����');
                    A=input('�����ĸϵ����ע��Ϊ����[1,10,50]����');
                    %{
                    Ҳ�����滻Ϊ���г����ʵ�ֶ�̬���������
                    B=input('�������ϵ����ע��Ϊ���󣩣�');
                    A=input('�����ĸϵ����ע��Ϊ���󣩣�');
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
                    title('��Ƶ����');

                    subplot(2,2,4);
                    plot(w,angle(H));grid on;
                    xlabel('\omega(rad/s)'),ylabel('\phi(\omega)');
                    title('H(s)����Ƶ����');
            
    elseif y==6
                f1=input('�����ź�Ƶ��(���磺1000)��f1=');%1000;
                f2=input('�����ź�Ƶ�ʣ����磺2500����f2=');%2500;
                f3=input('�����ź�Ƶ�ʣ����磺3800����f3=');%3800;
                p=input('�����һ����ֹƵ��[0~1]֮�䣨���磺0.5����w=');%3800;
                %}
                %f1=double(f1);f2=double(f2);f3=double(f3);p=double(p);
                s=p+0.05;
                f=[f1 f2 f3];
                ff=max(f);

                N=1024;      %�������ٵ��źŵĵ���
                fs=8000;       %����ǳ���Ƶ�ʣ�Ҫ�����ź������Ƶ�ʵ�2��
                t=(0:N-1)/fs;     %�ź�ʱ���������
                f=fs*(0:1023)/2048;
                f=(0:N-1)*fs/2048;    %�ź�Ƶ�������������������Ҫ���봦��һ��
                A1=0.05;A2=0.10;A3=0.08;
                x=A1*sin(2*pi*f1*t)+A2*sin(2*pi*f2*t)+A3*sin(2*pi*f3*t);   %������ĺ���

                wp=p*pi;
                ws=s*pi;
                Rp=1;
                Rs=15;
                Fs=8000;
                Ts=1/Fs;
                wp1=2/Ts*tan(wp/2);                 %��ģ��ָ��ת��������ָ��
                ws1=2/Ts*tan(ws/2); 
                [N,Wn]=buttord(wp1,ws1,Rp,Rs,'s');  %ѡ���˲�������С����
                [Z,P,K]=buttap(N);                  %����butterworthģ���˲���
                [Bap,Aap]=zp2tf(Z,P,K);
                [b,a]=lp2lp(Bap,Aap,Wn);   
                [bz,az]=bilinear(b,a,Fs);           %��˫���Ա任��ʵ��ģ���˲����������˲�����ת��
                [H,W]=freqz(bz,az);                 %����Ƶ����Ӧ����
                figure(1)
                plot(W*Fs/(2*pi),abs(H))
                grid on;axis tight;
                xlabel('Ƶ�ʣ�Hz��')
                ylabel('Ƶ����Ӧ')
                title('Butterworth')

                f1=filter(bz,az,x);
                figure(2)
                subplot(2,1,1)
                plot(t,x)                          %�����˲�ǰ��ʱ��ͼ
                grid on;axis tight;
                title('�˲�ǰ��ʱ����');
                subplot(2,1,2)
                plot(t,f1);                         %�����˲����ʱ��ͼ
                grid on;axis tight;
                title('�˲����ʱ����');
                %}

                y3=fft(f1,2048);
                figure(3)
                y2=fft(x,2048);
                subplot(2,1,1);
                plot(f,abs(y2(1:1024)));             %�����˲�ǰ��Ƶ��ͼ
                grid on;axis tight;
                title('�˲�ǰ��Ƶ��')
                xlabel('Hz');
                ylabel('����');
                subplot(2,1,2)
                plot(f,abs(y3(1:1024)));          %�����˲����Ƶ��ͼ
                grid on;axis tight;
                title('�˲����Ƶ��')
                xlabel('Hz');
                ylabel('����');
                warning off��


    end
end