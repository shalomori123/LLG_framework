function []=time_plot_2D_fields(E,dz,z1, z2,speed,varargin)
%plots a matrix in time in 3D
%ussage: with saving AVI file:time_plot_2D_fields(E,dz,medium,speed,'testing123') 
%        without saving AVI file: time_plot_2D_fields(E,dz,medium,speed)
name=[varargin{:}];

lim=max(max(max(E)));
%Y_axis=5;

[D ,I, N]=size(E);
z=0:dz:dz*(I-1);

[x_s,y_s] = meshgrid(-lim:lim/3:lim);
z1_s= z1*ones(size(x_s));
z2_s= z2* ones(size(x_s));

if isempty (name)
    for n=1:speed:N
        
        plot3(z,E(1,:,n),E(2,:,n),'linewidth',2);
        hold on
        %axis([0 z(end) -1.1*lim 1.1*lim -1.1*lim 1.1*lim]);
        surf(z1_s, x_s,y_s, 'FaceAlpha',0.3 , 'EdgeColor', 'none', 'FaceColor', 'r');
        surf(z2_s,x_s,y_s, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'r');

               
           
            
            
           
        hold off
        drawnow;
        if mod(n,300)==0
            xlabel('z')
            ylabel('E_x)')
            zlabel('E_y')
        end
    end
else
    v= VideoWriter(name, 'Uncompressed AVI');
    open(v)
    for n=1:speed:N
        plot3(z,E(1,:,n),E(2,:,n),'linewidth',1.5);
        hold on
        axis([ 0 z(end)   -1.1*lim 1.1*lim -1.1*lim 1.1*lim])
        surf(z1_s, x_s,y_s, 'FaceAlpha',0.3, 'EdgeColor', 'none', 'FaceColor', 'r');
        surf(z2_s,x_s,y_s, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'r');
        hold off
        drawnow;
        frame = getframe;
        writeVideo(v,frame)
    end
    close(v)
       hold on;

end