function [x, y, data]=oshunplot(directory,quantity,time)



if (time < 1e1 && time >= 0)
    timestr = ['0000' num2str(time)];
elseif (time < 1e2 && time >= 1e1)
    timestr = ['000' num2str(time)];
elseif (time < 1e3 && time >= 1e2)
    timestr = ['00' num2str(time)];
elseif (time < 1e4 && time >= 1e3)
    timestr = ['0' num2str(time)];
else
    timestr = num2str(time);
end
    
[x, y, info, data, typeofquantity] = oshun_getdata(directory,quantity,timestr);

if strcmp(typeofquantity,'FLD')
    
    plot(x,data);
    ylabel([info.Datasets.Attributes(1).Value '/' info.Datasets.Attributes(2).Value]);
    xlabel('x (c/\omega_p)');
    y = 0;
    
elseif strcmp(typeofquantity,'f')
    
    
    subplot(2,2,1);
    contourf(x,y,real(log10(data(:,:,1))),64,'LineColor','none');
    xlabel('x (c/\omega_p)');
    ylabel('p (m_e c)');
    colorbar;
    title(['real(' quantity ')']);
    
    subplot(2,2,2);
    contourf(x,y,real(log10(data(:,:,2))),64,'LineColor','none');
    xlabel('x (c/\omega_p)');
    ylabel('p (m_e c)');
    title(['imag(' quantity ')']);
    colorbar;
    
    subplot(2,2,3);
    contourf(x,y,data(:,:,1),64,'LineColor','none');
    xlabel('x (c/\omega_p)');
    ylabel('p (m_e c)');
    title(['real(' quantity ')']);
    colorbar;
    
    subplot(2,2,4);
    contourf(x,y,data(:,:,2),64,'LineColor','none');
    xlabel('x (c/\omega_p)');
    ylabel('p (m_e c)');
    title(['imag(' quantity ')']);
    colorbar;
    
elseif strcmp(typeofquantity,'fulldist')
    
    contourf(x,y,real((data)),64,'LineColor','none');
    xlabel('x (c/\omega_p)');
    ylabel('p (m_e c)');
    title(quantity);
    colorbar;
    
    
end





   
end



