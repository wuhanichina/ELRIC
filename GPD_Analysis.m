%% GPD_Analysis
% 用极值统计方法对英国日最大负荷做POT分析
% Peaks Over Threshold (POT) all peak flows that are greater than a given threshold flow

run('Read_UK_LoadData.m');% 加载数据
Load_2013_Column = reshape(Load_2013',48*365,1); % 将负荷数据重新排列至一列 
Load_2014_Column = reshape(Load_2014',48*365,1); 
Load_2015_Column = reshape(Load_2015',48*365,1);

% 用ksdensity绘制日负荷最大值概率分布图
[f,xi] = ksdensity(Load_Max_2013); % 2013年数据
q005 = quantile(Load_Max_2013,0.005); % 下0.5%分位点
q995 = quantile(Load_Max_2013,0.995); % 上0.5%分位点
figure
plot(xi,f);
hold on
line([q005;q005],[0;7e-5],'Color','red','LineStyle','--');
hold on
line([q995;q995],[0;7e-5],'Color','red','LineStyle','--');

% 用2参数GPD拟合尾部分布
alpha = 0.005; % 置信区间
[parmhat,parmci] = gpfit(Load_Max_2013,alpha); % 拟合100(1-alpha)置信区间下的GPD分布参数



%p = gppdf(x,k,sigma,theta)