clc
clear all
% 打开 fig 文件
fig = openfig('opt_radar1.fig');

% 找到图形对象，根据需求自动调整参数
lines = findobj(fig, 'type', 'line');

% 获取数据
xdata = get(lines, 'XData');
ydata = get(lines, 'YData');



fig2 = openfig('opt_radar2.fig');
% 找到图形对象，根据需求自动调整参数
lines2 = findobj(fig2, 'type', 'line');

% 获取数据
xdata2 = get(lines2, 'XData');
ydata2 = get(lines2, 'YData');
figure(3)
plot(xdata2{1,1},ydata2{1,1},'y--',xdata2{3,1},ydata2{3,1},'blue');

% plot(xdata{5,1},ydata{5,1},'r--',xdata{6,1},ydata{6,1},'black--');
% xlim([-90,120]);
% ylim([-20,25]);