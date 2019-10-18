function [ImgOrg, nc] = FillCircle(ImgOrg, icolor, d, ncir, gap)

n=size(ImgOrg,1);
low = floor((0.5+gap)*d); high = ceil(n-(0.5+gap)*d);
cnt = zeros(ncir,2);
cnt(1,:) = round(random('Uniform',low, high, [1,2]));
nc = 1;
ntry = 1;
d2 = ((1+gap)*d)^2;
while ((ntry<10*ncir)&(nc<ncir))
    a = round(random('Uniform',low, high, [1,2]));
    dmin = min((cnt(:,1)-a(1)).^2 + (cnt(:,2)-a(2)).^2);
    if dmin > d2
        nc = nc+1;
        cnt(nc,:) = a;
    end
    ntry = ntry + 1;
end
cnt = cnt(1:nc,:);

[r,c]=imgCircle(d);
for ii = 1:nc;
    ind = cnt(ii,1)+round(r)+(cnt(ii,2)+round(c)-1)*n;
    ind = ind(ind>0&ind<n*n); ind = unique(ind);
    ImgOrg(ind) = icolor;
end

end