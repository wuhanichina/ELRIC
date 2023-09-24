function y = new_cvar_equal(n, delta_Pi, d, C_max, Asset, L, r, alpha, beta, mu, sigma, xi)
% n 年限; delta_Pi 注入功率变化量; d 折旧率; C_max 线路容量; Asset 设备价值; L 单位切负荷损失 dollar/MW; r 利率;
% kappa 线路功率系数; beta 线路功率系数;
% mu 位置参数; sigma 尺度参数; xi 形状参数;
y = C_max*L+L*((1+r)^n*alpha*sigma+xi*(C_max-(1+r)^n*alpha*mu-(1+r)^n*alpha*delta_Pi-beta))/(1-xi)-Asset/(1+d)^n;
end