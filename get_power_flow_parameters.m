function [para_aij,para_beta,para_cij,para_delta,gamma_ij,tau_ij] = get_power_flow_parameters(casename,inbus,V_R)
%GET_POWER_FLOW_PARAMETERS 计算线性化潮流参数
%   casename Matpower case
%   inbus 注入节点编号
%   V_R = 1.05; % 平衡节点电压
Distribution_Network = ext2int(casename);
[Ybus,~,~] = makeYbus(casename); % case59China不存在shunt elements
G = real(Ybus);
B = imag(Ybus);

Distribution_Network = ext2int(casename);
% makeSbus这里对节点注入功率做了标幺化，除以-100
Sbus = makeSbus(Distribution_Network.baseMVA, Distribution_Network.bus, Distribution_Network.gen); 
Pbus = real(Sbus);
Qbus = tan(acos(0.95))*real(Sbus); % imag(Sbus);
nbus = length(Distribution_Network.bus);
nline = length(Distribution_Network.branch);

% 尝试用公式(9)将平衡节点、PV节点、PQ节点分开,后缀R、S、L分别代表平衡节点、PV节点、PQ节点
% 相应的节点导纳矩阵Y也要修改一下
VT_BUS = Distribution_Network.bus((Distribution_Network.bus(:,2)==3),1); % 平衡节点节点编号
PV_BUS = Distribution_Network.bus((Distribution_Network.bus(:,2)==2),1); % PV节点节点编号
PQ_bus = Distribution_Network.bus((Distribution_Network.bus(:,2)==1),1); % PQ节点节点编号
Y_rr = Ybus(VT_BUS,VT_BUS); G_rr = real(Y_rr); B_rr = imag(Y_rr);
Y_rs = Ybus(VT_BUS,PV_BUS); G_rs = real(Y_rs); B_rs = imag(Y_rs);
Y_rl = Ybus(VT_BUS,PQ_bus); G_rl = real(Y_rl); B_rl = imag(Y_rl);
Y_sr = Ybus(PV_BUS,VT_BUS); G_sr = real(Y_sr); B_sr = imag(Y_sr);
Y_ss = Ybus(PV_BUS,PV_BUS); G_ss = real(Y_ss); B_ss = imag(Y_ss);
Y_sl = Ybus(PV_BUS,PQ_bus); G_sl = real(Y_sl); B_sl = imag(Y_sl);
Y_lr = Ybus(PQ_bus,VT_BUS); G_lr = real(Y_lr); B_lr = imag(Y_lr);
Y_ls = Ybus(PQ_bus,PV_BUS); G_ls = real(Y_ls); B_ls = imag(Y_ls);
Y_ll = Ybus(PQ_bus,PQ_bus); G_ll = real(Y_ll); B_ll = imag(Y_ll);
Ybus_RSL = [Y_rr,Y_rs,Y_rl;Y_sr,Y_ss,Y_sl;Y_lr,Y_ls,Y_ll]; % 式(10)，重写的导纳矩阵
H = -[B_ss,B_sl;B_ls,B_ll]; % 式(13)
N = -[-G_sl;-G_ll];
M = -[G_ls,G_ll];
L = -B_ll;
H_dash = H-N*inv(L)*M; % 式(18)
L_dash = L-M*inv(H)*N; % 式(19)

% % 对于平衡节点的V、theta值和PV节点的P、V值，我们是已知的。于是利用式(12)计算P_dash和Q_dash
theta_R = 0; % 平衡节点相角
V_S = []; % PV节点电压
% S_dash = [Pbus(PV_BUS);Pbus(PQ_bus);Qbus(PQ_bus)] + ...
%     [B_sr,-G_sr,-G_ss; B_lr,-G_lr,-G_ls;G_lr,B_lr,B_ls]*[theta_R;V_R;V_S];  % 式(12)

% P_dash = S_dash(1:(length(PV_BUS)+length(PQ_bus)));  % 已知节点的有功功率
% Q_dash = S_dash((length(PV_BUS)+length(PQ_bus)+1):end) ;  % 已知节点的无功功率
% % 于是对公式(8)进行一下改造,得到式(20)和(21)

% theta_dash = inv(H_dash)*P_dash - inv(H_dash)*N*inv(L)*Q_dash; % 式(20) theta_s,theta_l PV节点、PQ节点相角
% V_dash = inv(L_dash)*Q_dash - inv(L_dash)*M*inv(H)*P_dash; % 式(21), V_l, PQ节点电压
% 
% V_final(VT_BUS) = V_R; V_final(PV_BUS) = V_S; V_final(PQ_bus) = V_dash; % 最终计算得到的节点电压
% theta_final(VT_BUS) = theta_R; theta_final([PV_BUS,PQ_bus]) = theta_dash; % 最终计算得到的节点电压相角

from = Distribution_Network.branch(1:nline,1);
to = Distribution_Network.branch(1:nline,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 首先我们要把上面的矩阵表达式拆分成单个元素，以便对每个节点的电压单独建模
R_Ps = [B_sr,-G_sr,-G_ss]*[theta_R;V_R;V_S];
R_Pl = [B_lr,-G_lr,-G_ls]*[theta_R;V_R;V_S];
R_Ql = [G_lr,B_lr,B_ls]*[theta_R;V_R;V_S];
A = zeros(nbus);
Beta = ones(nbus,1);
A(2:end,2:end) = tan(acos(0.95))*inv(L_dash)-inv(L_dash)*M*inv(H); % 这个A矩阵是个32*32矩阵
Beta(2:end) = inv(L_dash)*R_Ql-inv(L_dash)*M*inv(H)*R_Pl;% 这里全是1

para_aij = zeros(length(Pbus),1); % 节点i的系数
para_beta = zeros(length(Pbus),1);% i节点需要加上的系数
for j=1:length(Pbus) % 节点编号从1开始,与上一版code区别
    para_aij(j) = A(j,inbus);% 注意每换一个注入节点，这两个参数都要重新计算
    para_beta(j) = A(j,setdiff(2:nbus,inbus))*Pbus(setdiff(PQ_bus,inbus))+Beta(j); % Pbus就不需要调整节点编号
end

% 计算相角Va
C = zeros(length(Pbus));
D = zeros(length(Pbus),1);
% 注意C、D矩阵是32维的，计算前要把属于松弛节点的对应项补上
C(2:end,2:end) = inv(H_dash)-tan(acos(0.95))*inv(H_dash)*N*inv(L); % 原始C矩阵是个32*32矩阵
D(2:end) = inv(H_dash)*R_Pl-inv(H_dash)*N*inv(L)*R_Ql;% 

para_cij = zeros(length(Pbus),1); % 节点i的系数
para_delta = zeros(length(Pbus),1);% i节点需要加上的系数
for j=2:length(Pbus) % 这里C矩阵里节点序号-1，是因为C矩阵没有考虑平衡节点1，节点编号从2开始
    para_cij(j) = C(j,inbus);% 注意每换一个注入节点，这两个参数都要重新计算
    para_delta(j) = C(j,setdiff(2:nbus,inbus))*Pbus(setdiff(PQ_bus,inbus))+D(j); % Pbus就不需要调整节点编号
end

gamma_ij = zeros(nline,1); % 线路Pij的末端系数
tau_ij = zeros(nline,1); % 线路Pij的sum系数

for l=1:nline
    gamma_ij(l) = real(Ybus(from(l),to(l)))*(para_aij(from(l))-para_aij(to(l)))-imag(Ybus(from(l),to(l)))*(para_cij(from(l))-para_cij(to(l)));
    tau_ij(l) =  real(Ybus(from(l),to(l)))*(para_beta(from(l))-para_beta(to(l)))-imag(Ybus(from(l),to(l)))*(para_delta(from(l))-para_delta(to(l)));
end
end

