%% Read_UK_LoadData 
% 从英国负荷数据集中读取数据
% 此过程耗时较长不建议每次都执行
% uncertainty2015data.xlsx；uncertainty2014data.xlsx；uncertainty2013data.xlsx
% 中分别保存英国某点全年365天，每半小时一点的负荷数据，这里取其中的日负荷最大值数据

% 读取2015年的日负荷及峰值
[Load_2015,] = xlsread('uncertainty2015data.xlsx','365-48','B2:AW366');
[Load_Max_2015,] = xlsread('uncertainty2015data.xlsx','365-48','AX2:AX366');
[row,col]=find(Load_2015==0); % 寻找数据中的0值
for i=1:length(row)
    Load_2015(row(i),col(i)) = 0.5*(Load_2015(row(i),col(i)+1)+Load_2015(row(i),col(i)-1)); % 将每一个0值替换为其左右两个值之间的均值
end

% 读取2014年的日负荷及峰值
[Load_2014,] = xlsread('uncertainty2014data.xlsx','365-48','B2:AW366');
[Load_Max_2014,] = xlsread('uncertainty2015data.xlsx','365-48','AX2:AX366');
[row,col]=find(Load_2014==0); % 寻找数据中的0值
for i=1:length(row)
    Load_2014(row(i),col(i)) = 0.5*(Load_2014(row(i),col(i)+1)+Load_2014(row(i),col(i)-1)); % 将每一个0值替换为其左右两个值之间的均值
end

% 读取2013年的日负荷及峰值
[Load_2013,] = xlsread('uncertainty2013data.xlsx','365-48','B2:AW366');
[Load_Max_2013,] = xlsread('uncertainty2013data.xlsx','365-48','AX2:AX366');
[row,col]=find(Load_2013==0); % 寻找数据中的0值
for i=1:length(row)
    Load_2013(row(i),col(i)) = 0.5*(Load_2013(row(i),col(i)+1)+Load_2013(row(i),col(i)-1)); % 将每一个0值替换为其左右两个值之间的均值
end