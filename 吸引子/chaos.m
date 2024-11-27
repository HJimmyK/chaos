%%%%%%%%%
% 吸引子可视化
% 


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%1-Chen-Lee Attractor 

clc,clear
chen_lee_ab = @(t,yi)( [ ...
    5.*yi(1) - yi(2).*yi(3) , ...
    -10.*yi(2) + yi(1).*yi(3) , ...
    -0.38.*yi(3)+yi(1).*yi(2)./3 ...
    ]);


p_len = 6;
p0 = zeros(p_len.*p_len,3);
% 生成点
A = 10;
for i=1:p_len
    for j=1:p_len
        p0((i-1)*p_len+j,:) = A.*[cos(2.*pi.*(i./p_len)),sin(2.*pi.*(i./p_len)).*sin(2.*pi.*(j./p_len)),sin(2.*pi.*(i./p_len)).*cos(2.*pi.*(j./p_len))];
    end
end
p_len = length(p0);

% 计算的点的数量
len = 100000;

x=p0(:,1)' .* ones(len,p_len);
y=p0(:,2)' .* ones(len,p_len);
z=p0(:,3)' .* ones(len,p_len);
step0 = 10^(-3);


for j = 1:p_len
    step = step0;
    for i = 2:len
        [x(i,j),y(i,j),z(i,j),step] = myfun2(x(i-1,j),y(i-1,j),z(i-1,j),step);
    end
end

% 彼此两点间开始计算的差值点数
pnum = 1000;
x_ = zeros(len+(p_len-1)*pnum,p_len);
y_ = zeros(len+(p_len-1)*pnum,p_len);
z_ = zeros(len+(p_len-1)*pnum,p_len);
for i=1:p_len
    tempx = (x(1,i).*ones(1,(i-1).*pnum))';
    tempy = (y(1,i).*ones(1,(i-1).*pnum))';
    tempz = (z(1,i).*ones(1,(i-1).*pnum))';
    tempx_ = (x(1,i).*ones(1,(p_len-i).*pnum))';
    tempy_ = (y(1,i).*ones(1,(p_len-i).*pnum))';
    tempz_ = (z(1,i).*ones(1,(p_len-i).*pnum))';
    x_(:,i) = [tempx;x(:,i);tempx_];
    y_(:,i) = [tempy;y(:,i);tempy_];
    z_(:,i) = [tempz;z(:,i);tempz_];
end


%% 动画

clc
figure(Color='k',Position=[100, 100, 1000, 800])
scatter3(0,0,0,5,'filled','MarkerEdgeColor',"red","MarkerFaceColor","red")
hold on
scatter3(p0(:,1),p0(:,2),p0(:,3),5,'filled','MarkerEdgeColor',"k","MarkerFaceColor","k")
axis([-35 35 -35 35 0 40])
ax = gca;
grid off
ax.Color = 'k';
ax.XColor = 'k';
ax.YColor = 'k';
ax.ZColor = 'k';



startColor = [0.8, 0.5, 0.2]; % 橙色 RGB  
endColor = [0, 0.5, 1];   % 青色 RGB  
colors = zeros(p_len,3);

for i=1:p_len
    temh(i).h = animatedline(x_(1,i),y_(1,i),z_(1,i),'MaximumNumPoints',7000);
    tcolor = (i - 1) / (p_len - 1); % 计算插值比例   
    temh(i).h.Color = (1 - tcolor) * startColor + tcolor * endColor; % 线性插值 
end



% % 创建 VideoWriter 对象  
% v = VideoWriter('temp.mp4', 'MPEG-4'); % 设置视频文件名和格式  
% open(v); % 打开 VideoWriter 对象以开始写入  

for i=1:len
    for j=1:p_len
        temh(j).h.addpoints(x_(i,j),y_(i,j),z_(i,j))
        % 刷新图形  
        drawnow limitrate; % 使用 limitrate 限制绘制速率，提高效率  
        
    end
    % % 将当前帧写入视频  
    % if mod(i,2)
    %     continue;
    % end
    % if mod(i,3)
    %     continue;
    % end
    % frame = getframe(gcf); % 捕获当前图形窗口  
    % writeVideo(v, frame);  % 将帧写入视频  
end
% close(v)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rene-thomas 吸引子
clc,clear
chen_lee_ab = @(t,yi)( [ ...
    5.*yi(1) - yi(2).*yi(3) , ...
    -10.*yi(2) + yi(1).*yi(3) , ...
    -0.38.*yi(3)+yi(1).*yi(2)./3 ...
    ]);
b = 0.19;
rene_thomas_ab = @(t,yi)( [ ...
    sin(yi(2)) - b.*sin(yi(1)) , ...
    sin(yi(3)) - b.*sin(yi(2)) , ...
    sin(yi(1)) - b.*sin(yi(3)) ...
    ]);

p_len = 3;
p0 = zeros(p_len.*p_len,3);
% 生成点
A = 0.1;
for i=1:p_len
    for j=1:p_len
        p0((i-1)*p_len+j,:) = A.*[cos(2.*pi.*(i./p_len)),sin(2.*pi.*(i./p_len)).*sin(2.*pi.*(j./p_len)),sin(2.*pi.*(i./p_len)).*cos(2.*pi.*(j./p_len))];
    end
end
p_len = length(p0);

% 计算的点的数量
len = 100000;
x = zeros(len,p_len);
y = zeros(len,p_len);
z = zeros(len,p_len);
step0 = 10^(-3);

for j=1:p_len
    temp = zeros(len,3);
    temp(1,:) = p0(j,:);
    step = step0;
    for i=2:len
        [temp(i,:), step] = runge_kutta_chage(rene_thomas_ab, 1, temp(i-1,:), step);
    end
    x(:,j) = temp(:,1);
    y(:,j) = temp(:,2);
    z(:,j) = temp(:,3);
end



%% 绘图
clc
figure(Color='k')
scatter3(0,0,0,5,'filled','MarkerEdgeColor',"re","MarkerFaceColor","red")
hold on
%axis([-80 80 -75 75 0 60])
ax = gca;
grid off
ax.Color = 'k';
ax.XColor = 'w';
ax.YColor = 'w';
ax.ZColor = 'w';


maincolor = "#009efa";
pcolor  = "#00d2fc";
maincolors = generate_gradient_colors("#ffc75f",maincolor,  p_len);
pcolors = generate_gradient_colors( "#ffc75f",pcolor, p_len);
for j=1:p_len
scatter3(x(:,j),y(:,j),z(:,j),1,'filled','MarkerEdgeColor',maincolors(j,:),'MarkerFaceColor',maincolors(j,:));
hold on
plot3(x(:,j),y(:,j),z(:,j),'LineWidth',0.15,'Color',pcolors(j,:));
hold on
end

%% 动图
figure(Color='k',Position=[100, 100, 1000, 800])
scatter3(0,0,0,5,'filled','MarkerEdgeColor',"re","MarkerFaceColor","red")
hold on
axis([-150 150 -190 190 -50 250])
ax = gca;
grid off
ax.Color = 'k';
ax.XColor = 'k';
ax.YColor = 'k';
ax.ZColor = 'k';

startColor = [0.8, 0.5, 0.2]; % 橙色 RGB  
endColor = [0, 0.5, 1];   % 青色 RGB  
colors = zeros(p_len,3);

for i=1:p_len
    temh(i).h = animatedline(x(1,i),y(1,i),z(1,i),'MaximumNumPoints',10000);
    tcolor = (i - 1) / (p_len - 1); % 计算插值比例   
    temh(i).h.Color = (1 - tcolor) * startColor + tcolor * endColor; % 线性插值 
end
for i=1:len
    for j=1:p_len
        temh(j).h.addpoints(x(i,j),y(i,j),z(i,j))
        % 刷新图形  
        drawnow limitrate; % 使用 limitrate 限制绘制速率，提高效率  
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tsucs2 吸引子
clc,clear

p_len = 6;
p0 = zeros(p_len.*p_len,3);
% 生成点
A = 1;
for i=1:p_len
    for j=1:p_len
        p0((i-1)*p_len+j,:) = A.*[cos(2.*pi.*(i./p_len)),sin(2.*pi.*(i./p_len)).*sin(2.*pi.*(j./p_len)),sin(2.*pi.*(i./p_len)).*cos(2.*pi.*(j./p_len))];
    end
end
p_len = length(p0);


tsucs_ab = @(t,yi)( [ ...
    40.*(yi(2)-yi(1)) + 0.16.*yi(1).*yi(3) , ...
    55.*yi(1) - yi(1).*yi(3) + 20.*yi(2) , ...
    -0.65.*(yi(1))^2 + yi(1).*yi(2) + 1.833.*yi(3) ...
    ]);

% 计算的点的数量
len = 100000;
x = zeros(len,p_len);
y = zeros(len,p_len);
z = zeros(len,p_len);
step0 = 10^(-3);

for j=1:p_len
    temp = zeros(len,3);
    temp(1,:) = p0(j,:);
    step = step0;
    for i=2:len
        [temp(i,:), step] = runge_kutta_chage(tsucs_ab, 1, temp(i-1,:), step);
    end
    x(:,j) = temp(:,1);
    y(:,j) = temp(:,2);
    z(:,j) = temp(:,3);
end

%% 绘图
clc
figure(Color='k')
scatter3(0,0,0,5,'filled','MarkerEdgeColor',"re","MarkerFaceColor","red")
hold on
%axis([-80 80 -75 75 0 60])
ax = gca;
grid off
ax.Color = 'k';
ax.XColor = 'w';
ax.YColor = 'w';
ax.ZColor = 'w';


maincolor = "#009efa";
pcolor  = "#00d2fc";
maincolors = generate_gradient_colors("#ffc75f",maincolor,  p_len);
pcolors = generate_gradient_colors( "#ffc75f",pcolor, p_len);
for j=1:p_len
scatter3(x(1:4:end,j),y(1:4:end,j),z(1:4:end,j),1,'filled','MarkerEdgeColor',maincolors(j,:),'MarkerFaceColor',maincolors(j,:));
hold on
plot3(x(1:2:end,j),y(1:2:end,j),z(1:2:end,j),'LineWidth',0.15,'Color',pcolors(j,:));
hold on
end


%% 动图

figure(Color='k',Position=[100, 100, 1000, 800])
scatter3(0,0,0,5,'filled','MarkerEdgeColor',"re","MarkerFaceColor","red")
hold on
axis([-150 150 -190 190 -50 250])
ax = gca;
grid off
ax.Color = 'k';
ax.XColor = 'k';
ax.YColor = 'k';
ax.ZColor = 'k';

startColor = [0.8, 0.5, 0.2]; % 橙色 RGB  
endColor = [0, 0.5, 1];   % 青色 RGB  
colors = zeros(p_len,3);

for i=1:p_len
    temh(i).h = animatedline(x(1,i),y(1,i),z(1,i),'MaximumNumPoints',10000);
    tcolor = (i - 1) / (p_len - 1); % 计算插值比例   
    temh(i).h.Color = (1 - tcolor) * startColor + tcolor * endColor; % 线性插值 
end
for i=1:len
    for j=1:p_len
        temh(j).h.addpoints(x(i,j),y(i,j),z(i,j))
        % 刷新图形  
        drawnow limitrate; % 使用 limitrate 限制绘制速率，提高效率  
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rucklidge 吸引子
clc,clear

p_len = 6;
p0 = zeros(p_len.*p_len,3);
% 生成点
A = 1;
for i=1:p_len
    for j=1:p_len
        p0((i-1)*p_len+j,:) = A.*[cos(2.*pi.*(i./p_len)),sin(2.*pi.*(i./p_len)).*sin(2.*pi.*(j./p_len)),sin(2.*pi.*(i./p_len)).*cos(2.*pi.*(j./p_len))];
    end
end
p_len = length(p0);


tsucs_ab = @(t,yi)( [ ...
    -2.*yi(1) + 6.7.*yi(2) - yi(2).*yi(3) , ...
    yi(1) , ...
    -yi(3) + yi(2)^2 ...
    ]);

% 计算的点的数量
len = 100000;
x = zeros(len,p_len);
y = zeros(len,p_len);
z = zeros(len,p_len);
step0 = 10^(-3);

for j=1:p_len
    temp = zeros(len,3);
    temp(1,:) = p0(j,:);
    step = step0;
    for i=2:len
        [temp(i,:), step] = runge_kutta_chage(tsucs_ab, 1, temp(i-1,:), step);
    end
    x(:,j) = temp(:,1);
    y(:,j) = temp(:,2);
    z(:,j) = temp(:,3);
end

%% 绘图
clc
figure(Color='k')
scatter3(0,0,0,5,'filled','MarkerEdgeColor',"re","MarkerFaceColor","red")
hold on
%axis([-80 80 -75 75 0 60])
ax = gca;
grid off
ax.Color = 'k';
ax.XColor = 'w';
ax.YColor = 'w';
ax.ZColor = 'w';


maincolor = "#009efa";
pcolor  = "#00d2fc";
maincolors = generate_gradient_colors("#ffc75f",maincolor,  p_len);
pcolors = generate_gradient_colors( "#ffc75f",pcolor, p_len);
for j=1:p_len
scatter3(x(1:4:end,j),y(1:4:end,j),z(1:4:end,j),1,'filled','MarkerEdgeColor',maincolors(j,:),'MarkerFaceColor',maincolors(j,:));
hold on
plot3(x(1:2:end,j),y(1:2:end,j),z(1:2:end,j),'LineWidth',0.15,'Color',pcolors(j,:));
hold on
end


%% 动图


figure(Color='k',Position=[100, 100, 1000, 800])
scatter3(0,0,0,5,'filled','MarkerEdgeColor',"re","MarkerFaceColor","red")
hold on
axis([-14 14 -8 8 -2 16])
ax = gca;
grid off
ax.Color = 'k';
ax.XColor = 'k';
ax.YColor = 'k';
ax.ZColor = 'k';
view([-2,-2,-5])
startColor = [0.8, 0.5, 0.2]; % 橙色 RGB  
endColor = [0, 0.5, 1];   % 青色 RGB  
colors = zeros(p_len,3);
pause(20);
for i=1:p_len
    temh(i).h = animatedline(x(1,i),y(1,i),z(1,i),'MaximumNumPoints',10000);
    tcolor = (i - 1) / (p_len - 1); % 计算插值比例   
    temh(i).h.Color = (1 - tcolor) * startColor + tcolor * endColor; % 线性插值 
end
for i=1:len
    for j=1:p_len
        temh(j).h.addpoints(x(i,j),y(i,j),z(i,j))
        % 刷新图形  
        drawnow limitrate; % 使用 limitrate 限制绘制速率，提高效率  
    end
end




function gradient_colors = generate_gradient_colors(maincolor, pcolor, num)  
    % 将十六进制颜色转换为RGB  
    main_rgb = hex2rgb(maincolor);  
    p_rgb = hex2rgb(pcolor);  
    
    % 生成主颜色与目标颜色之间的RGB渐变  
    gradient_colors = zeros(num, 3); % 初始化颜色数组  
    for i = 1:num  
        % 计算插值  
        t = (i - 1) / (num - 1); % 范围从0到1  
        gradient_colors(i, :) = main_rgb * (1 - t) + p_rgb * t; % 线性插值  
    end  
end  

function rgb = hex2rgb(hex)  
    % 将十六进制颜色转换为RGB格式  
    hex = char(hex);  
    rgb = sscanf(hex(2:end), '%2x') / 255; % 从十六进制获取RGB并归一化  
end  


function [yi_, step_] = runge_kutta_chage(odefun, ti, yi, step)
    k1 = odefun(ti, yi);  
    k2 = odefun(ti + step/2, yi + step/2 * k1);  
    k3 = odefun(ti + step/2, yi + step/2 * k2);  
    k4 = odefun(ti + step, yi + step * k3); 
    yi_ = yi + step/6 * (k1 + 2*k2 + 2*k3 + k4); 
    step_ = step;
    no = sum(abs((yi-yi_)./yi));
    if no < 0.01
        step_ = step .* 1.05;
    end
    if no > 0.015
        step_ = step .* 0.99;
    end

end

function [t, y] = runge_kutta4(odefun, tspan, init_conditions, dt)  
    % 初始化时间和状态  
    t = tspan(1):dt:tspan(2);  
    n = length(t);  
    y = zeros(length(init_conditions), n);  
    y(:, 1) = init_conditions;  

    for i = 1:n-1  
        ti = t(i);  
        yi = y(:, i);  

        k1 = odefun(ti, yi);  
        k2 = odefun(ti + dt/2, yi + dt/2 * k1);  
        k3 = odefun(ti + dt/2, yi + dt/2 * k2);  
        k4 = odefun(ti + dt, yi + dt * k3);  

        y(:, i+1) = yi + dt/6 * (k1 + 2*k2 + 2*k3 + k4);  
    end  
end

function [x, y, z, step] = myfun2(x_, y_, z_, step_t)  
    a = 5;  
    b = -10;  
    c = -0.38;  

    % 当前状态导数  
    k1x = a .* x_ - y_ .* z_;         
    k1y = b .* y_ + x_ .* z_;        
    k1z = c .* z_ + (x_ .* y_) / 3;  

    % 计算中间值  
    x_temp = x_ + 0.5 .* k1x .* step_t;  
    y_temp = y_ + 0.5 .* k1y .* step_t;  
    z_temp = z_ + 0.5 .* k1z .* step_t;  
    
    % 中间导数  
    k2x = a .* x_temp - y_temp .* z_temp;         
    k2y = b .* y_temp + x_temp .* z_temp;        
    k2z = c .* z_temp + (x_temp .* y_temp) / 3;  

    % 计算中间值  
    x_temp = x_ + 0.5 .* k2x .* step_t;  
    y_temp = y_ + 0.5 .* k2y .* step_t;  
    z_temp = z_ + 0.5 .* k2z .* step_t;  
    
    % 中间导数  
    k3x = a .* x_temp - y_temp .* z_temp;         
    k3y = b .* y_temp + x_temp .* z_temp;        
    k3z = c .* z_temp + (x_temp .* y_temp) / 3;  

    % 计算最后的中间值  
    x_temp = x_ + k3x .* step_t;  
    y_temp = y_ + k3y .* step_t;  
    z_temp = z_ + k3z .* step_t;  
    
    % 最后的导数  
    k4x = a .* x_temp - y_temp .* z_temp;         
    k4y = b .* y_temp + x_temp .* z_temp;        
    k4z = c .* z_temp + (x_temp .* y_temp) / 3;  

    % 更新新的值  
    x = x_ + (k1x + 2.*k2x + 2.*k3x + k4x) .* step_t / 6;         
    y = y_ + (k1y + 2.*k2y + 2.*k3y + k4y) .* step_t / 6;        
    z = z_ + (k1z + 2.*k2z + 2.*k3z + k4z) .* step_t / 6;  
    
    step = step_t;
    
   
    if abs((x(1)-x_(1))./x_(1))+abs((y(1)-y_(1))./y_(1))+abs((z(1)-z_(1))./z_(1)) < 0.01
        step = 1.005.*step_t;
    end
    if abs((x(1)-x_(1))./x_(1))+abs((y(1)-y_(1))./y_(1))+abs((z(1)-z_(1))./z_(1)) > 0.015
        step = 0.99.*step_t;
    end
end



