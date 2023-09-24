function [n_old_LRIC, n_new_LRIC, LRIC] = LRIC_Gen_fun(case_name,r,d,Asset,C_max)
%LRIC_FUN 计算经典长期增量成本
% r 负荷增长率
% d 折现率
% Asset 设备成本
% C_max 线路容量

MPC_case = ext2int(case_name);
H = makePTDF(MPC_case); % DCPTDF矩阵

delta_Pi = 5; % 节点注入功率增量

P_i = MPC_case.bus(:,3); % 节点负荷Pd 第三列
for i = 1:length(P_i) % 每个节点注入一次功率计算LRIC
    P_i_old = P_i;
    P_i_old(i) = -P_i(i);
    P_l = abs(H*P_i_old); % 原始线路功率
    n_old_LRIC = zeros(length(P_l),1); % 旧的年限定义
    PV_LRIC_old = zeros(length(P_l),1); % 线路功率定义
    AF_LRIC_old = zeros(length(P_l),1); % 年金因子的定义
    for l=1:length(P_l)
        n_old_LRIC(l) = (log(C_max)-log(P_l(l)))/log(1+r); % 现有的年限
        PV_LRIC_old(l) = Asset/((1+d).^n_old_LRIC(l)); % 当前净现值
        AF_LRIC_old(l) = 1/r-1/(1+r)^n_old_LRIC(l); % 年金系数
    end
end

n_new_LRIC = zeros(length(P_l),length(P_i)); % 新的年限定义
PV_LRIC_new = zeros(length(P_l),length(P_i)); % 线路净现值定义
AF_LRIC = zeros(length(P_l),1); % 年金系数定义
IC_LRIC = zeros(length(P_l),1); % 现值定义
LRIC = zeros(length(P_i),1); % 长期增量成本

for i = 1:length(P_i) % 每个节点注入一次功率计算LRIC
    P_i_new = P_i;
    P_i_new(i) = -P_i_new(i)-delta_Pi;
    P_l_new = abs(H*P_i_new);
    for l = 1:length(P_l_new)
        n_new_LRIC(l,i) = (log(C_max)-log(P_l_new(l)))/log(1+r); % 节点i新增注入功率后的年限
        PV_LRIC_new(l,i) = Asset/((1+d).^n_new_LRIC(l,i)); % 节点i新增注入功率后的净现值
        AF_LRIC(l) = 1/r-1/(1+r)^n_new_LRIC(l,i); % 年金系数
        IC_LRIC(l) = (PV_LRIC_new(l,i)-PV_LRIC_old(l))*AF_LRIC(l); % 现值
    end
LRIC(i) = sum(IC_LRIC)/delta_Pi;
end
end

