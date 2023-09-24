%% GPD_Analysis
% �ü�ֵͳ�Ʒ�����Ӣ������󸺺���POT����
% Peaks Over Threshold (POT) all peak flows that are greater than a given threshold flow

run('Read_UK_LoadData.m');% ��������
Load_2013_Column = reshape(Load_2013',48*365,1); % ��������������������һ�� 
Load_2014_Column = reshape(Load_2014',48*365,1); 
Load_2015_Column = reshape(Load_2015',48*365,1);

% ��ksdensity�����ո������ֵ���ʷֲ�ͼ
[f,xi] = ksdensity(Load_Max_2013); % 2013������
q005 = quantile(Load_Max_2013,0.005); % ��0.5%��λ��
q995 = quantile(Load_Max_2013,0.995); % ��0.5%��λ��
figure
plot(xi,f);
hold on
line([q005;q005],[0;7e-5],'Color','red','LineStyle','--');
hold on
line([q995;q995],[0;7e-5],'Color','red','LineStyle','--');

% ��2����GPD���β���ֲ�
alpha = 0.005; % ��������
[parmhat,parmci] = gpfit(Load_Max_2013,alpha); % ���100(1-alpha)���������µ�GPD�ֲ�����



%p = gppdf(x,k,sigma,theta)