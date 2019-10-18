function [ r, c ] = imgRectangle( b, h, c)

rot = [cosd(c) -sind(c); sind(c) cosd(c)];
bh = round((b-1)/2);
hh = round((h-1)/2);
c = (-bh:0.5:bh);
r = (-hh:0.5:hh)';
c = c(ones(1,4*hh+1),:);
c = c(:);
r = r(:,ones(4*bh+1,1));
r = r(:);
a = [c r]*rot;
a = floor(a);
a = unique(a,'rows');
c = a(:,1);
r = a(:,2);
end

