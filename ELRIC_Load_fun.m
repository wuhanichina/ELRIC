function [n_old_ELRIC, n_new_ELRIC, ELRIC] = ELRIC_Load_fun(case_name,load,eta,r,d,Asset,C_max)
%ELRIC_FUN 计算扩展长期增量成本
% case_name Matpower中的case
% load 负荷数据
% eta 单位阻塞管理 dollar/MW
% r 负荷增长率
% d 折现率
% Asset 设备成本
% C_max 线路容量
MPC_case = ext2int(case_name);
H = makePTDF(MPC_case); % DCPTDF矩阵
delta_Pi = 5; % 节点注入功率增量
P_i = MPC_case.bus(:,3); % 节点负荷Pd 第三列
P_l = abs(H*P_i); % 原始线路功率

q_lo = 0.1;
q_up = 0.9;
q_up_load = quantile(load,q_up); % 取上分位点

full_load_pd = paretotails(load,q_lo,q_up); % 完整负荷概率模型
pd_load_up_k = full_load_pd.UpperParameters(1); pd_load_up_s = full_load_pd.UpperParameters(2);

n_old_ELRIC = zeros(length(P_l),length(P_i)); % 各条线路的初始投资年限
n_new_ELRIC = zeros(length(P_l),length(P_i)); % 节点i新增负荷各条线路的投资年限
ELRIC = zeros(length(P_i),1); % 扩展长期增量成本定义
PV_eLRIC_old = zeros(length(P_l),length(P_i));% 节点i新增负荷前线路l投资的净现值
PV_eLRIC_new = zeros(length(P_l),length(P_i));% 节点i新增负荷后线路l投资的净现值
AF_eLRIC = zeros(length(P_l),1); % 年金系数定义
IC_eLRIC = zeros(length(P_l),1); % 现值定义
alpha = zeros(length(P_i),1); % 节点i的系数
beta = zeros(length(P_i),1);% i节点需要加上的系数
PCC_Voltage= 1;

for i = 1:length(P_i)
    % 第一步，计算n_old,即计算上面这个方程的解
    % 功率从节点i到线路j的转移
    for l=setdiff(1:length(P_l),find(P_l==0)) % 
        alpha(l) = abs(H(l,i));% 注意每换一个注入节点，这两个参数都要重新计算
        beta(l) = abs(H(l,setdiff(1:length(P_i),i))*MPC_case.bus(setdiff(1:length(P_i),i),3)); % Pbus就不需要调整节点编号
        ELRIC_get_n = @(n)cvar_equal(n, d, C_max, Asset, eta, r, alpha(l), beta(l), q_up_load, pd_load_up_s, pd_load_up_k);
        interval = [1  100];
        n_old_ELRIC(l,i) = fzero(ELRIC_get_n,interval); % 求解ELRIC_get_n=0
        % 第二步，计算原始设备投资成本的净现值
        % 为了减少循环次数，就放在这里一起做了
        PV_eLRIC_old(l,i) = Asset/(1+d)^n_old_ELRIC(l,i);
    end
    % 第三步，计算节点新增负荷之后的年限和净现值，公式在上面
    for l=setdiff(1:length(P_l),find(P_l==0)) % 线路14上没有功率。。。所以这里就略掉吧
        alpha(l) = abs(H(l,i));% 注意每换一个注入节点，这两个参数都要重新计算
        beta(l) = abs(H(l,setdiff(1:length(P_i),i))*MPC_case.bus(setdiff(1:length(P_i),i),3)); % Pbus就不需要调整节点编号
        ELRIC_get_n = @(n)new_cvar_equal(n, delta_Pi,d, C_max, Asset, eta, r, alpha(l), beta(l), q_up_load, pd_load_up_s, pd_load_up_k);
        interval = [1  100];
        n_new_ELRIC(l,i) = fzero(ELRIC_get_n,interval); % 求解ELRIC_get_n=0
        % 第二步，计算原始设备投资成本的净现值
        % 为了减少循环次数，就放在这里一起做了
        PV_eLRIC_new(l,i) = Asset/((1+d)^n_new_ELRIC(l,i));
    end
    % 第四步，计算净现值变化量、IC
    for l=setdiff(1:length(P_l),find(P_l==0))
        AF_eLRIC(l) = 1/r-1/(1+r)^n_new_ELRIC(l,i); % 年金系数
        IC_eLRIC(l) = (PV_eLRIC_new(l,i)-PV_eLRIC_old(l,i))*AF_eLRIC(l); % 现值
    end
    % 第五步，计算扩展长期增量成本
    ELRIC(i) = sum(IC_eLRIC)/delta_Pi;
end

end

