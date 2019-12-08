function myplot(number, value, time, xname, yname, legend1, legend2, legend3, png_name)
%MYPLOT Summary of this function goes here
%   Detailed explanation goes here
syms t real
value_time = double(subs(value, {t}, {time}));
fig = figure(number);
set(gcf,'color','w','Position',[0 0 800 600]);
plot(time, value_time, 'LineWidth', 2)
grid on
xlabel(xname)
ylabel(yname)
legend(legend1, legend2, legend3, 'interpreter', 'latex', 'Fontsize', 14)
set(gca,'FontSize',16)
frame = getframe(fig);
im = frame2im(frame);
[img,map] = rgb2ind(im,256);
imwrite(img,map,png_name,'png');
close;
end