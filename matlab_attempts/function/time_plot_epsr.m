function []=time_plot_epsr(E,dz,medium,speed,varargin)
%plots a matrix in time : F(time,length)
%ussage: with saving AVI file:time_plot_epsr(E,dz,medium,speed,'testing123') 
%        without saving AVI file: time_plot_epsr(E,dz,medium,speed)
name=[varargin{:}];

Y_axis=max(max(E));
%Y_axis=5;

[N ,I]=size(E);
z=0:dz:dz*(I-1);

if isempty (name)
    for n=1:speed:N
        plot(z,E(n,:),z,medium,'linewidth',2);
        axis([ 0 z(end)   -Y_axis Y_axis])
        drawnow;
    end;
else
    clear M;
    v= VideoWriter(name, 'Uncompressed AVI');
    open(v)
    for n=1:speed:N
        plot(z,E(n,:),z,medium,'linewidth',2);
        axis([ 0 z(end)   -1.15*Y_axis 1.15*Y_axis])
        drawnow;
        frame = getframe;
        writeVideo(v,frame)
    end;
    close(v)
end