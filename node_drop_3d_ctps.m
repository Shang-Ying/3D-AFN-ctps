function xyz = node_drop_3d_ctps (b_box, ninit, dotmax, radius, ctps, vargin)
% Function: advancing_front_node_generation
% Author: Shang-Ying Chen
% Date: April, 2024
% Description: This MATLAB function is based on the code from Van der Sande (2019), modified to include the functionality of control points (ctps).
% It generates advancing front nodes within a three-dimensional bounding box (b_box) using a moving window search.
% Parameters ninit and radius represent the initial PDPs and exclusion radius, respectively.
% The output is generated node coordinates (xyz), suitable for meshless methods.

if nargin < 6
    vargin = [];
end

dotnr   = 0;
rng(0);
xyz      = zeros(dotmax,3);
excess_height = 0.1;
dx = (b_box(2)-b_box(1))/(ninit(1)-1);
dy = (b_box(4)-b_box(3))/(ninit(2)-1);
r = radius;
pdp = b_box(5)+0.01*min(r)*rand(ninit);
nodeindices = zeros(size(pdp));
[zm,idx]  = min(pdp(:));
[i1,i2] = ind2sub(size(pdp),idx);

xyz_new = [b_box(1)+dx*(i1-1),b_box(3)+dy*(i2-1),pdp(i1,i2)];
mdis = 0;

while zm <= (1+excess_height)*b_box(6) && dotnr < dotmax
    dotnr = dotnr + 1;
    xyz(dotnr,:) = xyz_new;
    nodeindices(i1,i2) = dotnr;
    r = radius_3d_ctps(radius, mdis);
    ileft  = max(1,i1 - floor(r/dx));
    iright = min(ninit(1),i1 + floor(r/dx));
    ibottom = max(1,i2 - floor(r/dy));
    itop = min(ninit(2),i2 + floor(r/dy));
    xx = ileft:iright;
    yy = ibottom:itop;
    [X,Y] = ndgrid(xx,yy);
    
    height =  sqrt(r^2 - (dx*(X-i1)).^2 - (dy*(Y-i2)).^2);
    pdp(ileft:iright,ibottom:itop) = (max(pdp(ileft:iright,ibottom:itop),pdp(i1,i2) + real(height)));
    [zm,ix] = min(pdp(ileft:iright,ibottom:itop));
    [zm,iy] = min(zm);
    i1 = ileft+ix(iy)-1;
    i2 = ibottom+iy-1;
    
    searchr = min(2*ceil(r/dx),floor(ninit(1)/2)-1);
    
    while 1
        if i1-searchr < 1
            xsearch = [ninit(1)+i1-searchr:ninit(1),1:i1+searchr];
        elseif i1+searchr > ninit(1)
            xsearch = [i1-searchr:ninit(1),1:i1+searchr-ninit(1)];
        else
            xsearch = i1-searchr:i1+searchr;
        end
        if i2-searchr < 1
            ysearch = [ninit(2)+i2-searchr:ninit(2),1:i2+searchr];
        elseif i2+searchr > ninit(2)
            ysearch = [i2-searchr:ninit(2),1:i2+searchr-ninit(2)];
        else
            ysearch = i2-searchr:i2+searchr;
        end
        [zm,ix] = min(pdp(xsearch,ysearch));
        [~,iy] = min(zm);
        ix = ix(iy);
        i1 = xsearch(ix);
        i2 = ysearch(iy);
        zm = pdp(i1,i2);
        
        if ix > searchr/2 && ix < length(xsearch)-searchr/2 && ...
                iy > searchr/2 && iy < length(ysearch)-searchr/2
            break
        end
    end

    dis_ctps = pdist2([b_box(1)+dx*(i1-1),b_box(3)+dy*(i2-1),pdp(i1,i2)], ctps);
    [mdis, mindex] = min(dis_ctps);
    r_ctps = radius_3d_ctps(0.5*radius, 0);
    if (mdis < r_ctps) == 1
        xyz_new = ctps(mindex,:);
    else
        xyz_new = [b_box(1)+dx*(i1-1),b_box(3)+dy*(i2-1),pdp(i1,i2)];
    end
end

xyz = xyz(1:dotnr,:);
xyz = xyz(xyz(:,3)<=b_box(6),:);
catchup = ismember(ctps,xyz,'row');
xyz = [xyz; ctps(~catchup,:)];
xyz = unique(xyz,'rows');