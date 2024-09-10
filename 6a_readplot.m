clear all
% 打开.fig文件并提取数据
fig = openfig('opt_radar1.fig','invisible');  % 'invisible'选项可在不显示窗口的情况下打开.fig文件

% 获取图像标题
title = get(gcf,'Name');

% 获取图像大小和位置
pos = get(gcf,'Position');

% 获取图像的颜色和样式
color = get(gca,'Color');
style = get(gca,'LineStyleOrder');

% 获取图像的标注和注释
x_label = get(gca,'XLabel');
y_label = get(gca,'YLabel');
title_label = get(gca,'Title');
legend_label = get(gca,'Legend');

% 获取图像数据
x_data = get(findobj(gcf,'Type','line'),'XData');
y_data = get(findobj(gcf,'Type','line'),'YData');

% 获取图像注释
annots = findobj(gcf,'Type','annotation');
annotations = cell(size(annots));
for i = 1:length(annots)
    annotations{i} = get(annots(i));
end

% 关闭.fig文件
close(fig);
X1=x_data{1,1};%propose best
Y1=y_data{1,1};
X2=x_data{2,1};%propose cnn
Y2=y_data{2,1};
X3=x_data{3,1};%desired
Y3=y_data{3,1};


fig2 = openfig('opt_radar2.fig','invisible');  % 'invisible'选项可在不显示窗口的情况下打开.fig文件
x2_data = get(findobj(gcf,'Type','line'),'XData');
y2_data = get(findobj(gcf,'Type','line'),'YData');
close(fig2);
X4=x_data{1,1};%propose OMP
Y4=y_data{1,1};
X5=x_data{2,1};%propose PEMO
Y5=y_data{2,1};

figure(2)
plot(X3,Y3,'Color',[0 0 0],'Linewidth',1.5);
axis([-90 130 -20 25])
hold on
plot(X1,Y1,'Color',[1 128/255 0],'Linewidth',1.5);
axis([-90 130 -20 25])
hold on
plot(X2,Y2,'Color',[0 1 1],'Linewidth',1.5);
axis([-90 130 -20 25])
hold on
plot(X4,Y4,'--','Color',[1 1 0],'Linewidth',1);
axis([-90 130 -20 25])
hold on
plot(X5,Y5,'--','Color',[65/255 105/255 225/255],'Linewidth',1);
axis([-90 130 -20 25])
hold on
legend('Desired','Propose Best','Propose CNN','OMP','PE-MO');