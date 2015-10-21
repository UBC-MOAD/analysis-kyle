clear all
close all
load ORCA2_deptht
deptht=double(deptht);
d18O=ncread('calculated_d18O_v1_1.nc', 'd18o');
depth=double(ncread('calculated_d18O_v1_1.nc', 'depth'));
load Exchange/coord_ANHA4.mat
test = permute(d18O, [2, 1, 3]);
d18O=double(test);
d18O(d18O<-1e+10)=nan;
clear test

x=ncread('calculated_d18O_v1_1.nc', 'lon');
y=ncread('calculated_d18O_v1_1.nc', 'lat');
[lon, lat]=meshgrid(x, y);
clear x y
for i=1:33
    layer=inpaint_nans(d18O(:, :, i), 3);
    layer(layer<=-7)=-7.0;
    layer(layer>5)=5;
    d18O(:, :, i)=layer;
end

d18O_ORCA2=zeros([size(d18O, 1), size(d18O, 2), 31]);
for i=1:size(d18O, 1)
    for j=1:size(d18O, 2)
        d18O_ORCA2(i, j, :)=interp1(depth, squeeze(d18O(i, j, :)), deptht, 'linear');
    end
end

d18O_ANHA4=zeros([size(d18O, 1), size(d18O, 2), 50]);
for i=1:size(d18O, 1)
    for j=1:size(d18O, 2)
        d18O_ANHA4(i, j, :)=interp1(depth, squeeze(d18O(i, j, :)), nav_lev, 'linear');
    end
end
d18O_ANHA4(:, :, 50)=d18O_ANHA4(:, :, 49);

save d18O d18O_ORCA2 d18O_ANHA4        
%test = interp2(lon, lat, d18O(:, :, 1), lon, lat, 'nearest');