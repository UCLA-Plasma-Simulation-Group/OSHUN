function [t, data]=oshunxtplot(directory,quantity,timerange)

i=1;

for time = timerange

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

timeaxis(i) = info.Attributes(2).Value;
datamax(i) = max(max(abs(data)));
i=i+1;
end


semilogy(timeaxis,datamax);
   
end



