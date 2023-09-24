%% Read_UK_LoadData 
% ��Ӣ���������ݼ��ж�ȡ����
% �˹��̺�ʱ�ϳ�������ÿ�ζ�ִ��
% uncertainty2015data.xlsx��uncertainty2014data.xlsx��uncertainty2013data.xlsx
% �зֱ𱣴�Ӣ��ĳ��ȫ��365�죬ÿ��Сʱһ��ĸ������ݣ�����ȡ���е��ո������ֵ����

% ��ȡ2015����ո��ɼ���ֵ
[Load_2015,] = xlsread('uncertainty2015data.xlsx','365-48','B2:AW366');
[Load_Max_2015,] = xlsread('uncertainty2015data.xlsx','365-48','AX2:AX366');
[row,col]=find(Load_2015==0); % Ѱ�������е�0ֵ
for i=1:length(row)
    Load_2015(row(i),col(i)) = 0.5*(Load_2015(row(i),col(i)+1)+Load_2015(row(i),col(i)-1)); % ��ÿһ��0ֵ�滻Ϊ����������ֵ֮��ľ�ֵ
end

% ��ȡ2014����ո��ɼ���ֵ
[Load_2014,] = xlsread('uncertainty2014data.xlsx','365-48','B2:AW366');
[Load_Max_2014,] = xlsread('uncertainty2015data.xlsx','365-48','AX2:AX366');
[row,col]=find(Load_2014==0); % Ѱ�������е�0ֵ
for i=1:length(row)
    Load_2014(row(i),col(i)) = 0.5*(Load_2014(row(i),col(i)+1)+Load_2014(row(i),col(i)-1)); % ��ÿһ��0ֵ�滻Ϊ����������ֵ֮��ľ�ֵ
end

% ��ȡ2013����ո��ɼ���ֵ
[Load_2013,] = xlsread('uncertainty2013data.xlsx','365-48','B2:AW366');
[Load_Max_2013,] = xlsread('uncertainty2013data.xlsx','365-48','AX2:AX366');
[row,col]=find(Load_2013==0); % Ѱ�������е�0ֵ
for i=1:length(row)
    Load_2013(row(i),col(i)) = 0.5*(Load_2013(row(i),col(i)+1)+Load_2013(row(i),col(i)-1)); % ��ÿһ��0ֵ�滻Ϊ����������ֵ֮��ľ�ֵ
end