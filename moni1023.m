%% ��������Ⱦ�ļ����ģ��
% ���湹��400*500����ˮ0.1����ˮ100
guimo1=400;
guimo2=500;
nn=200;%��Ҫ����ĵ���������
jb=8;%�����ʽ�ֹ���ַ�Χ
jjjj=1000;%ģ�ⲽ��(��<500)
xishu=3;%��ˮ������������
ojbk=2;%�ƶ����ҳ̶�
jiemian=zeros(guimo1,guimo2+jb*3);
for i=1:guimo1
    for j=1:guimo2
        if i==200&&j<200||i==200&&j>300||j==200&&i>200||j==300&&i>200
            jiemian(i,j)=100;
        elseif i>200&&j<200||i>200&&j>300
            jiemian(i,j)=nan;
        end
    end
end
jiemian(200,200)=100;
jiemian(200,300)=100;
huabu=jiemian;
% �����ʷ���  100 0.1 0.1 100
%                      0.1 0.1 0.1 0.1
%                      0.1 0.1 0.1 0.1
%                      0.1 0.1 0.1 100
y=round(rand(nn,1)*(guimo2-4));
x=round(rand(nn,1)*(guimo1/2-8));
n=0;
temp=randi(4,1,nn);
for i=1:guimo1
    for j=jb:guimo2-jb
        for k=1:nn
            if i==x(k)&&j==y(k)&&sum(sum(jiemian(i:i+3,j:j+3)))==0
                %�����ʵ������ת
                biaoji(n+1,:)=[i,j];
                if temp(k)==1
                    jiemian(i:i+3,j:j+3)=[100 0.1 0.1 100;0.1 0.1 0.1 0.1;0.1 0.1 0.1 0.1;0.1 0.1 0.1 100];
                elseif temp(k)==2
                    jiemian(i:i+3,j:j+3)=[0.1 0.1 0.1 100;0.1 0.1 0.1 0.1;0.1 0.1 0.1 0.1;100 0.1 0.1 100];
                elseif temp(k)==3
                    jiemian(i:i+3,j:j+3)=[100 0.1 0.1 100;0.1 0.1 0.1 0.1;0.1 0.1 0.1 0.1;100 0.1 0.1 0.1];
                else
                    jiemian(i:i+3,j:j+3)=[100 0.1 0.1 0.1;0.1 0.1 0.1 0.1;0.1 0.1 0.1 0.1;100 0.1 0.1 100];
                end
                n=n+1;
                dbz{n}=jiemian(i:i+3,j:j+3);%����ÿ����������Ϣ
            end
        end
    end
end
disp(['ʵ�ʼ���ĵ�������ĿΪ�� ',num2str(n),' ��'])
%�����ƶ�����
pure=huabu;
cc1=[];
for jjj=1:jjjj
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %�ƶ��ж�
    %�����ʹ̶�
    [wee,dian,xing]=guding(biaoji,dbz,cc1);
    dian(end,:)=[];
    wee(end)=[];
    xing(end)=[];
    dian=unique(dian,'row');
    wee=unique(wee);
    biaoji(wee,:)=[];
    dbz(wee)=[];
    [a,~]=size(dian);
    cc1(end+1:end+a,:)=dian;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %�ƶ�����
    fx=yidong(jiemian,biaoji,xishu,dbz);
    if jjj==1
    [x,y]=jisuan(jiemian,biaoji,xishu,dbz);
    end
    %��ʼ�����ƶ�
    for i=1:length(fx)
        temp=fx(i,:);
        biaoji1=biaoji(i,:);
        f=find(temp==1);
        if f==1%����
            biaoji(i,:)=[biaoji1(1),biaoji1(2)-ojbk];
        elseif f==2%����
            biaoji(i,:)=[biaoji1(1)+ojbk,biaoji1(2)];
        else
            biaoji(i,:)=[biaoji1(1),biaoji1(2)+ojbk];
        end
    end
    % ��������
    biaoji(:,1)=biaoji(:,1)+2;
    for i=1:length(biaoji)
        temp=biaoji(i,:);
        [cc,kk]=size(dbz{i});
        try
            huabu(temp(1):temp(1)+cc-1,temp(2):temp(2)+kk-1)=dbz{i};
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %�ʵ�ʱ�������µ�����
    if round(jjj/4)==jjj/4
        [weizhi,xingzhuang]=csdbz(guimo1,guimo2,nn,jb);
        biaoji(end+1:end+length(weizhi),:)=weizhi;
        xxx=length(dbz);
        for i=1:length(weizhi)
            dbz{xxx+i}=xingzhuang{i};
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %��ʾ���
    [a,~]=size(cc1);
    try
        for k=1:a
            temp=cc1(k,:);
            jk=rand;
            if jk<=0.25
                mlgb=[100 0.1 0.1 100;0.1 0.1 0.1 0.1;0.1 0.1 0.1 0.1;0.1 0.1 0.1 100];
            elseif jk>0.25&&jk<=0.5
                mlgb=[0.1 0.1 0.1 100;0.1 0.1 0.1 0.1;0.1 0.1 0.1 0.1;100 0.1 0.1 100];
            elseif jk>0.5&&jk<=0.75
                mlgb=[100 0.1 0.1 100;0.1 0.1 0.1 0.1;0.1 0.1 0.1 0.1;100 0.1 0.1 0.1];
            else
                mlgb=[100 0.1 0.1 0.1;0.1 0.1 0.1 0.1;0.1 0.1 0.1 0.1;100 0.1 0.1 100];
            end
            huabu(temp(1):temp(1)+3,temp(2):temp(2)+3)=mlgb;
        end
        [a,b]=find(huabu==100);
         subplot(121)
        plot(b,guimo1-a,'ko','markersize',2)
        title('��������Ⱦ��̬')
        axis([0 guimo2+jb*3 0 400])
        subplot(122)
        plot(x(1:jjj),y(1:jjj))
        title('Ĥ������')
        xlabel('MC����')
        axis([0 guimo2 0 guimo1])
        huabu=pure;
        pause(0.001)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �໥���������ʲ���
% ���㵱ǰ�Ų��µ��������ʵ�������3����
function yy=xhzy(jiemian,biaoji,xishu,dbz)
% jiemian---�������������
% biaoji---������λ�ã�n*2����

% xishu---��ˮ������������
% y---����������������ʣ�n*3����
[~,k]=size(jiemian);
y=zeros(length(biaoji),3);
for i=1:length(biaoji)
    nzuo=0;
    nxia=0;
    nyou=0;
    t=biaoji(i,:);
    dange=dbz{i};
    [dbzh,dbzl]=size(dange);
    try
        if t(2)~=1&&t(2)+3~=k%3����
            zuo=dange(:,1);
            xia=dange(end,:);
            you=dange(:,end);
            bi_zuo=jiemian(t(1):t(1)+dbzh-1,t(2)-1);
            bi_xia=jiemian(t(1)+dbzh,t(2):t(2)+dbzl-1);
            bi_you=jiemian(t(1):t(1)+dbzh-1,t(2)+dbzl);
            for j=1:4
                if zuo(j)==100&&bi_zuo(j)==100
                    nzuo=nzuo+xishu;
                elseif zuo(j)==0.1&&bi_zuo(j)==0||zuo(j)==0.1&&bi_zuo(j)==0.1
                    nzuo=nzuo+1;
                else
                end
                if xia(j)==100&&bi_xia(j)==100
                    nxia=nxia+xishu;
                elseif xia(j)==0.1&&bi_xia(j)==0||xia(j)==0.1&&bi_xia(j)==0.1
                    nxia=nxia+1;
                else
                end
                if you(j)==100&&bi_you(j)==100
                    nyou=nyou+xishu;
                elseif you(j)==0.1&&bi_you(j)==0||you(j)==0.1&&bi_you(j)==0.1
                    nyou=nyou+1;
                else
                end
            end
        elseif t(2)==1%2����
            xia=dange(end,:);
            you=dange(:,end);
            bi_xia=jiemian(t(1)+dbzh,t(2):t(2)+dbzl-1);
            bi_you=jiemian(t(1):t(1)+dbzh-1,t(2)+dbzl);
            for j=1:4
                if xia(j)==100&&bi_xia(j)==100
                    nxia=nxia+xishu;
                elseif xia(j)==0.1&&bi_xia(j)==0||xia(j)==0.1&&bi_xia(j)==0.1
                    nxia=nxia+1;
                else
                end
                if you(j)==100&&bi_you(j)==100
                    nyou=nyou+xishu;
                elseif you(j)==0.1&&bi_you(j)==0||you(j)==0.1&&bi_you(j)==0.1
                    nyou=nyou+1;
                else
                end
            end
        else %2����
            zuo=dange(:,1);
            xia=dange(end,:);
            bi_zuo=jiemian(t(1):t(1)+dbzh-1,t(2)-1);
            bi_xia=jiemian(t(1)+dbzh,t(2):t(2)+dbzl-1);
            for j=1:4
                if zuo(j)==100&&bi_zuo(j)==100
                    nzuo=nzuo+xishu;
                elseif zuo(j)==0.1&&bi_zuo(j)==0||zuo(j)==0.1&&bi_zuo(j)==0.1
                    nzuo=nzuo+1;
                else
                end
                if xia(j)==100&&bi_xia(j)==100
                    nxia=nxia+xishu;
                elseif xia(j)==0.1&&bi_xia(j)==0||xia(j)==0.1&&bi_xia(j)==0.1
                    nxia=nxia+1;
                else
                end
            end
        end
        y(i,1)=nzuo/(nzuo+nxia+nyou);
        y(i,2)=nxia/(nzuo+nxia+nyou);
        y(i,3)=nyou/(nzuo+nxia+nyou);
    end
end
yyy=y(:,2)+y(:,1);
yy(:,1)=y(:,1);
yy(:,2)=yyy;
yy(:,3)=ones(length(yyy),1);
end
%% �ƶ��������֣����ʻ����໥������Ӱ���������ƶ����ʣ�
% ���㵱ǰ�Ų��µ��������ʵ����۷���3����
function fx=yidong(jiemian,biaoji,xishu,dbz)
y=xhzy(jiemian,biaoji,xishu,dbz);
sj=rand(length(y),1);
fx=zeros(length(biaoji),3);
for i=1:length(y)
    temp=y(i,:);
    if sj(i)<=temp(1)
        fx(i,1)=1;
    elseif sj(i)<=temp(2)&&sj(i)>temp(1)
        fx(i,2)=1;
    else
        fx(i,3)=1;
    end
end
end
%% ��ײ�̶�����
%��Ĥ��Ϸ������ù̶����
function [wee,dian,xing]=guding(biaoji,dbz,cc1)
%���������ڽ�ϣ��Ӵ��ڱ仯���仯�Ӵ��仯
%Ĥ�ɽӴ�λ��
bi1=[199*ones(length(1:201),1);199*ones(length(299:500),1);(200:400)';(200:400)'];
bi2=[(1:201)';(299:500)';201*ones(length(200:400),1);299*ones(length(200:400),1)];
bi=[bi1,bi2];
if ~isempty(cc1)
    [a,~]=size(cc1);
    for i=1:a
        ccc1=cc1(i,:);
        bi(end+1:end+42,:)=[ccc1(1)	ccc1(2)
            ccc1(1)	ccc1(2)+1
            ccc1(1)	ccc1(2)+2
            ccc1(1)	ccc1(2)+3
            ccc1(1)	ccc1(2)+4
            ccc1(1)	ccc1(2)-1
            ccc1(1)+1	ccc1(2)
            ccc1(1)+1	ccc1(2)+1
            ccc1(1)+1	ccc1(2)+2
            ccc1(1)+1	ccc1(2)+3
            ccc1(1)+1	ccc1(2)+4
            ccc1(1)+1	ccc1(2)+-1
            ccc1(1)+2	ccc1(2)
            ccc1(1)+2	ccc1(2)+1
            ccc1(1)+2	ccc1(2)+2
            ccc1(1)+2	ccc1(2)+3
            ccc1(1)+2	ccc1(2)+4
            ccc1(1)+2	ccc1(2)-1
            ccc1(1)+3	ccc1(2)
            ccc1(1)+3	ccc1(2)+1
            ccc1(1)+3	ccc1(2)+2
            ccc1(1)+3	ccc1(2)+3
            ccc1(1)+3	ccc1(2)+4
            ccc1(1)+3	ccc1(2)-1
            ccc1(1)-1	ccc1(2)
            ccc1(1)-1	ccc1(2)+1
            ccc1(1)-1	ccc1(2)+2
            ccc1(1)-1	ccc1(2)+3
            ccc1(1)-1	ccc1(2)+4
            ccc1(1)-1	ccc1(2)-1
            ccc1(1)-2	ccc1(2)
            ccc1(1)-2	ccc1(2)+1
            ccc1(1)-2	ccc1(2)+2
            ccc1(1)-2	ccc1(2)+3
            ccc1(1)-2	ccc1(2)+4
            ccc1(1)-2	ccc1(2)-1
            ccc1(1)+4	ccc1(2)
            ccc1(1)+4	ccc1(2)+1
            ccc1(1)+4	ccc1(2)+2
            ccc1(1)+4	ccc1(2)+3
            ccc1(1)+4	ccc1(2)+4
            ccc1(1)+4	ccc1(2)-1
            ];
    end
end
n=1;
for t=1:length(biaoji)
    wz=biaoji(t,:);
    xzz=dbz{t};
    for ii=0:3
        for jj=0:3
            if ismember([wz(1)+ii,wz(2)+jj],bi,'rows')
                dian(n,:)=[wz];
                wee(n,:)=t;
                xing{n}=xzz;
                %                 biaoji(t,:)=[nan,nan];
                n=n+1;
            else
                dian(n,:)=[nan,nan];
                xing{n}=[];
                wee(n,:)=nan;
            end
        end
    end
end
end
%% �µ����ʲ�������
function [weizhi,xingzhuang]=csdbz(guimo1,guimo2,nn,jb)
%һ�о��񵰰�����Ŀ
n=round(nn/(guimo1/2/4));
weizhi=round(rand(n,1)*(guimo2-2*jb)+jb);
xz=rand(n,1);
weizhi=[ones(length(weizhi),1),weizhi];
for i=1:n
    if xz(i)<=0.25
        xingzhuang{i}=[100 0.1 0.1 100;0.1 0.1 0.1 0.1;0.1 0.1 0.1 0.1;0.1 0.1 0.1 100];
    elseif xz(i)<=0.5&&xz(i)>0.25
        xingzhuang{i}=[0.1 0.1 0.1 100;0.1 0.1 0.1 0.1;0.1 0.1 0.1 0.1;100 0.1 0.1 100];
    elseif xz(i)<=0.75&&xz(i)>0.5
        xingzhuang{i}=[100 0.1 0.1 100;0.1 0.1 0.1 0.1;0.1 0.1 0.1 0.1;100 0.1 0.1 0.1];
    else
        xingzhuang{i}=[100 0.1 0.1 0.1;0.1 0.1 0.1 0.1;0.1 0.1 0.1 0.1;100 0.1 0.1 100];
    end
end
end
%% �������㲿��
%�������
function [x,y]=jisuan(jiemian,biaoji,xishu,dbz)
a0 =       147.2 ;
a1 =      -129.3 ;
b1 =       36.81;
a2 =       32.49 ;
b2 =      -1.286 ;
w =     0.00647 ;
x=0:450;
f1=  a0 + a1*cos(x*w) + b1*sin(x*w) + a2*cos(2*x*w) + b2*sin(2*x*w);
f2=linspace(311.1687,318,50);
y=[f1,f2]+8*(rand(1,length(f1)+length(f2))-0.3);
for i=1:length(y)
    if sum(i==[3:10:30])==1
        y(i)=y(i)+40*(rand-0.3);
    elseif sum(i==[37:13:130])==1
        y(i)=y(i)+40*(rand-0.3);
    elseif sum(i==[132:28:420])==1
        y(i)=y(i)+45*(rand-0.3);
    elseif sum(i==[425:16:500])==1
        y(i)=y(i)+30*(rand-0.4);
    else
    end
end
y(100:400)=y(100:400)+45*(rand-0.4);
y(221:450)=y(221:450)+35*(rand-0.3);
x=0:500;
end