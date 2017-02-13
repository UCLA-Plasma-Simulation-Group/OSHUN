function [x,y,info, data,typeofquantity]= oshun_getdata(directory,quantity,timestr)

if (strcmp(quantity,'Ex') || strcmp(quantity,'Ey') || strcmp(quantity,'Ez') || ...
        strcmp(quantity,'Bx') || strcmp(quantity,'By') || strcmp(quantity,'Bz'))

    subdir = ['OUTPUT/FLD/' quantity];
    typeofquantity = 'FLD';
    filename = [directory '/' subdir '/' quantity '_' timestr '.h5'];
    info = h5info(filename);field = info.Datasets.Name;data = h5read(filename,['/' field]);
    
elseif (strcmp(quantity,'ne') || strcmp(quantity,'Jx') )

    subdir = ['OUTPUT/MOM/' quantity];
    typeofquantity = 'FLD';
    filename = [directory '/' subdir '/' quantity '_s0_' timestr '.h5'];
    info = h5info(filename);field = info.Datasets.Name;data = h5read(filename,['/' field]);    
    
elseif (strcmp(quantity,'T') )

    subdir = ['OUTPUT/MOM/' quantity];
    typeofquantity = 'FLD';
    filename = [directory '/' subdir '/' quantity '_eV_s0_' timestr '.h5'];
    info = h5info(filename);field = info.Datasets.Name;data = h5read(filename,['/' field]);        
    
elseif (strcmp(quantity,'f0-x') || strcmp(quantity,'f10-x') )
    
    subdir = ['OUTPUT/DISTR/' quantity];
    typeofquantity = 'f';
    filename = [directory '/' subdir '/' quantity '_s0_' timestr '.h5'];
    info = h5info(filename);field = info.Datasets.Name;data = h5read(filename,['/' field]);

elseif (strcmp(quantity,'px-x') )
    
    subdir = ['OUTPUT/DISTR/' quantity];
    typeofquantity = 'fulldist';
    filename = [directory '/' subdir '/' quantity '_s0_' timestr '.h5'];
    info = h5info(filename);field = info.Datasets.Name;data = h5read(filename,['/' field]);
    
    
elseif strcmp(quantity,'df')
    subdir = ['OUTPUT/DISTR/'];
    typeofquantity = 'fulldist';
    filename = [directory '/' subdir '/px-x/px-x_s0_' timestr '.h5'];
    info = h5info(filename);field = info.Datasets.Name;px = h5read(filename,['/' field]);
    
    filename1 = [directory '/' subdir '/f0-x/f0-x_s0_' timestr '.h5'];
    info = h5info(filename1);field = info.Datasets.Name;f00 = h5read(filename1,['/' field]);
    
    f00 = [flipud(f00); f00(1,:,:); f00];
    
    data = px-f00(:,:,1);
    data((size(data,1)-1)/2+1,:,:) = 0.0;
   
end


if strcmp(typeofquantity,'FLD')
    numaxis1 = size(data,1);
    axis1lims = h5read(filename,['/AXIS/AXIS1']);
    x = linspace(axis1lims(1),axis1lims(2),numaxis1);
    y = 0;
    
elseif strcmp(typeofquantity,'f')
    numaxis1 = size(data,2);
    axis1lims = h5read(filename,['/AXIS/AXIS1']);
    x = linspace(axis1lims(1),axis1lims(2),numaxis1);
    
    numaxis2 = size(data,1);
    axis2lims = h5read(filename,['/AXIS/AXIS2']);
    y = linspace(axis2lims(1),axis2lims(2),numaxis2);    
    
elseif strcmp(typeofquantity,'fulldist')
    numaxis1 = size(data,2);
    axis1lims = h5read(filename,['/AXIS/AXIS1']);
    x = linspace(axis1lims(1),axis1lims(2),numaxis1);
    
    numaxis2 = size(data,1);
    axis2lims = h5read(filename,['/AXIS/AXIS2']);
    y = linspace(-axis2lims(2),axis2lims(2),numaxis2);        
    
end



end
