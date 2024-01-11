function [t,y] = odefunc(func,y0,tspan,step,sel)
    n=floor((tspan(2)-tspan(1))/step)+1;
    y=zeros(n,length(y0));
    y(1,:)=y0;
    t=0;
    if strcmp(sel,'euler')
        for i=2:n
            y(i,:)=y(i-1,:)+step*func(t,y(i-1,:))';
        end
    elseif strcmp(sel,'rk4')
        for i=2:n
            k1=func(t,y(i-1,:))';
            k2=func(t+step/2,y(i-1,:)+step*k1/2)';
            k3=func(t+step/2,y(i-1,:)+step*k2/2)';
            k4=func(t+step,y(i-1,:)+step*k3)';
            y(i,:)=y(i-1,:)+step*(k1+2*k2+2*k3+k4)/6;
        end
    else
        error('Invalid method');
    end
end