function [] = plot_moments(M_array_global,dt,moments_in_figure)

[M,N,T ]= size(M_array_global);
step= ceil(N/moments_in_figure);

T_Global= linspace(0,dt*(T-1),T);
figure('Position',[600 0 1000 190*moments_in_figure])
tiledlayout(N/step,1);
for i = 1:step:N
    nexttile;
    plot(T_Global,squeeze(M_array_global(1,i,:)),'linewidth',2,'DisplayName','M_x');
    hold on
    plot(T_Global,squeeze(M_array_global(2,i,:)),'linewidth',2,'DisplayName','M_y');
    plot(T_Global,squeeze(M_array_global(3,i,:)),'linewidth',4,'DisplayName','M_z')
    plot(T_Global,squeeze(((M_array_global(1,i,:)).^2 + (M_array_global(2,i,:)).^2 + (M_array_global(3,i,:)).^2).^0.5),'linewidth',1.5,'DisplayName','Magnetization Absolute Value');
    hold off
    title(['Moment num: ' num2str(i)])
    xlabel('time [sec]'), ylabel('magentizasion value [A/m]')
    legend()
end
