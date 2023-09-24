function [n_old_LRIC, n_new_LRIC, LRIC] = LRIC_Load_fun_Dis(case_name,r,d,Asset,C_max)
%LRIC_FUN 计算经典长期增量成本
% r 负荷增长率
% d 折现率
% Asset 设备成本
% C_max 线路容量

MPC_case = ext2int(case_name);
delta_Pi = 0.1; % 节点注入功率增量

P_i = MPC_case.bus(:,3); % 节点负荷Pd 第三列
Q_i = MPC_case.bus(:,4); % 节点负荷Pd 第三列
nbus = length(P_i);
nline = nbus-1;
n_old_LRIC = zeros(nline,nbus); % 旧的年限定义
PV_LRIC_old = zeros(nline,nbus); % 线路功率定义
AF_LRIC_old = zeros(nline,nbus); % 年金因子的定义
PCC_Voltage = 1;
for i = 1:nbus % 每个节点注入一次功率计算LRIC
    [~,~,~,P_l,~,~,~] = Compute_Energy_Loss('case59China',[P_i';Q_i']'); % 单位分别为W,A,p.u.
    P_l = abs(P_l);
    for l=1:nline
        n_old_LRIC(l,i) = (log(C_max)-log(P_l(l)))/log(1+r); % 现有的年限
        PV_LRIC_old(l,i) = Asset/((1+d).^n_old_LRIC(l,i)); % 当前净现值
        AF_LRIC_old(l,i) = 1/r-1/(1+r)^n_old_LRIC(l,i); % 年金系数
    end
end

n_new_LRIC = zeros(nline,nbus); % 新的年限定义
PV_LRIC_new = zeros(nline,nbus); % 线路净现值定义
AF_LRIC = zeros(nline,nbus); % 年金系数定义
IC_LRIC = zeros(nline,nbus); % 现值定义
LRIC = zeros(nbus,1); % 长期增量成本

for i = 1:nbus % 每个节点注入一次功率计算LRIC
    delta_p = zeros(nbus,1);
    delta_p(i) = delta_Pi;
    P_i_new = P_i + delta_p;
    Q_i_new = Q_i + delta_p*tan(0.95);
    [~,~,~,P_l_new,~,~,~] = Compute_Energy_Loss('case59China',[P_i_new';Q_i_new']'); % 单位分别为W,A,p.u.
    P_l_new = abs(P_l_new);
    for l = 1:nline
        n_new_LRIC(l,i) = (log(C_max)-log(P_l_new(l)))/log(1+r); % 节点i新增注入功率后的年限
        PV_LRIC_new(l,i) = Asset/((1+d).^n_new_LRIC(l,i)); % 节点i新增注入功率后的净现值
        AF_LRIC(l,i) = 1/r-1/(1+r)^n_new_LRIC(l,i); % 年金系数
        IC_LRIC(l,i) = (PV_LRIC_new(l,i)-PV_LRIC_old(l,i))*AF_LRIC(l,i); % 现值
    end
LRIC(i) = sum(IC_LRIC(:,i))/delta_Pi;
end
end