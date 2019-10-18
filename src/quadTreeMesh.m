function [coord, ele, eleQT, eleColor, eleSize, eleCentre, eleDof] = quadTreeMesh(S,I,voidColor)
[M,N] = size(S);
maxDim = full(max(max(S))); minDim = full(min(min(S(S>0))));
if nargin < 3; voidColor = []; end

%% enforce 2:1 ratio
nSplit = 1;
J = repmat(maxDim,size(S));
dim = maxDim/2;
while dim >= minDim;
    blocks = find(S==dim)';
    numBlocks = length(blocks);
    if numBlocks>0
        % Compute block indices for a dim-by-dim block.
        ind = blockIndex(M, dim, blocks);
        J(ind) = dim;
    end
    dim = dim/2;
end

[sr, sc, sdim] = find(S);
sLen = length(sr);
sLenIncr = sLen;
sr = [sr; zeros(sLenIncr,1)];
sc = [sc; zeros(sLenIncr,1)];
sdim = [sdim; zeros(sLenIncr,1)];
sLenMax = 2*sLen;
while nSplit>0
    dim = maxDim;
    nSplit = 0;
    while (dim > 2*minDim)
        % Find all the blocks at the current size.
        a = find(sdim==dim);
        r = sr(a); c = sc(a);
        
        dim2 = dim/2;
        doSplit = zeros(length(r),1);
        
        bindx = find(r>1);
        if ~isempty(bindx)
            Sind = r(bindx)-1+M*(c(bindx)-1);
            ind = (0:M:(dim-1)*M)';
            ind = bsxfun(@plus, ind(:), Sind');
%            Jblks = reshape(J(ind), dim, []);
            Jblks = J(ind);
            doSplit(bindx(min(Jblks)<dim2)) = 1;
        end
        bindx = find((r+dim) < M);
        if ~isempty(bindx)
            Sind = r(bindx) + dim + M*(c(bindx)-1);
            ind = (0:M:(dim-1)*M)';
            ind = bsxfun(@plus, ind(:), Sind');
%            Jblks = reshape(J(ind), dim, []);
            Jblks = J(ind);
            doSplit(bindx(min(Jblks)<dim2)) = 1;
        end
        bindx = find(c>1);
        if ~isempty(bindx)
            Sind = r(bindx) + M*(c(bindx)-2);
            ind = (0:(dim-1))';
            ind = bsxfun(@plus, ind(:), Sind');
%            Jblks = reshape(J(ind), dim, []);
            Jblks = J(ind);
            doSplit(bindx(min(Jblks)<dim2)) = 1;
        end
        bindx = find((c+dim) < N);
        if ~isempty(bindx)
            Sind = r(bindx) + M*(c(bindx)+dim-1);
            ind = (0:(dim-1))';
            ind = bsxfun(@plus, ind(:), Sind');
%            Jblks = reshape(J(ind), dim, []);
            Jblks = J(ind);
            doSplit(bindx(min(Jblks)<dim2)) = 1;
        end
        
        sDoSplit = find(doSplit>0);
        if sLen+length(sDoSplit) > sLenMax
            sr = [sr; zeros(sLenIncr,1)];
            sc = [sc; zeros(sLenIncr,1)];
            sdim = [sdim; zeros(sLenIncr,1)];
            sLenMax = sLenMax + sLenIncr;
        end
        % Record results in output matrix.
%         blocks =  (r(sDoSplit) + M*(c(sDoSplit)-1))';
%         ind = blockIndex(M, dim, blocks);
%         J(ind(:)) = dim2;
        for ii = sDoSplit'
            J(r(ii):r(ii)+dim-1,c(ii):c(ii)+dim-1) = dim2;
            sdim( a(ii) ) = dim2;
            
            sLen = sLen + 1;
            sr(sLen) = r(ii)+dim2; sc(sLen) = c(ii);  sdim(sLen) = dim2;
            sLen = sLen + 1;
            sr(sLen) = r(ii); sc(sLen) = c(ii)+dim2; sdim(sLen) = dim2;
            sLen = sLen + 1;
            sr(sLen) = r(ii)+dim2; sc(sLen) = c(ii)+dim2; sdim(sLen) = dim2;
        end
        dim = dim2;
        nSplit = nSplit + length(sDoSplit);
    end
end

S = sparse(sr(sr>0),sc(sr>0), sdim(sr>0), M, N);

%% get color in elements
maxDim = full(max(max(S))); minDim = full(min(min(S(S>0))));
disp(['    Maximum element size = ',num2str(maxDim)]);
disp(['    Minimum element size = ',num2str(minDim)]);

nLevel = 1+nextpow2(maxDim/minDim);
r = cell(nLevel,1);
c = cell(nLevel,1);
eleColor = cell(nLevel,1);
dim = maxDim;
ii = 1;
while dim >= minDim;
%    [eleColor{ii},r{ii},c{ii}] = qtgetblkavg(I,S,dim);
    [r{ii}, c{ii}] = find(S == dim);
    if isempty(r{ii});  continue;  end
    Sind = (r{ii} + M*(c{ii}-1))';
    ind = blockIndex(M, dim, Sind); 
    eleColor{ii} = (sum(I(ind),1))'/(dim*dim);
    ii = ii+1;
    dim = dim/2;
end
r=cell2mat(r); c=cell2mat(c);
eleColor = cell2mat(eleColor);
if ~isempty(voidColor) %remove elements in voids
    if length(voidColor) == 1
        a = eleColor~=voidColor;
    else
        a = (eleColor<voidColor(1) | eleColor>voidColor(2));
    end
    r = r(a); c = c(a);
    eleColor = uint8(round(eleColor(a)));
end

%% get nodal coordinates and construct element connectivity
nEle = length(r);

Sind = r+M*(c-1);
eleSize = J(Sind);
x1 = c-0.5;  y1=r-0.5; x2=x1+eleSize; y2=y1+eleSize;
eleCentre = [(x1+x2)/2 (y1+y2)/2];
eCd = [ x1  x2  x2 x1 zeros(nEle,4) y2 y2 y1 y1 zeros(nEle,4) ]';
eCd = reshape(eCd,8,2*nEle);

Jnorth = maxDim*ones(nEle,1); Jsouth = Jnorth;  Jwest = Jnorth;  Jeast = Jnorth;
b = r > 1; Sind = r(b)-1 + M*(c(b)-1); 
Jnorth(b) = J(Sind);
b = (r+eleSize) < M; Sind = r(b) + eleSize(b) + M*(c(b)-1);
Jsouth(b) = J(Sind);
b = c > 1; Sind = r(b) + M*(c(b)-2);
Jwest(b) = J(Sind);
b = (c+eleSize) < N; Sind = r(b) + M*(c(b)+eleSize(b)-1);
Jeast(b) = J(Sind);
clear J;

ele = cell(nEle,1);
eleQT = zeros(nEle,1);
eleNode = zeros(nEle,1);
eleCoord = cell(nEle,1);
a = [1:4 1:4];
for ii = 1:nEle
    dim = eleSize(ii); dim2 = dim/2;
    exy = eCd(:,2*ii-1:2*ii);
    inode = 4;  mNode = zeros(1,4);
    if dim > minDim
        if Jsouth(ii) < dim
            inode = inode + 1; mNode(1) = inode;
            exy(inode,:) = [ exy(1,1)+dim2 exy(1,2)];
        end
        if Jeast(ii) < dim
            inode = inode + 1; mNode(2) = inode;
            exy(inode,:) = [exy(2,1) exy(3,2)+dim2];
        end
        if Jnorth(ii) < dim
            inode = inode + 1; mNode(3) = inode;
            exy(inode,:) = [ exy(1,1)+dim2 exy(3,2)];
        end
        if Jwest(ii) < dim
            inode = inode + 1; mNode(4) = inode;
            exy(inode,:) = [ exy(1,1) exy(3,2)+dim2];
        end
    end
    switch inode
        case 4
            qtType = 1;
        case 5
            qtType = 2; k = find(mNode>0);
            exy = exy([k mNode(k) a(k+1:k+3)],:);
        case 6
            k = find(mNode>0);
            if k(2)-k(1) == 1
                qtType = 3; k = k(1);
                exy = exy([k mNode(k) a(k+1) mNode(k+1) a(k+2) a(k+3)],:);
            elseif k(2)-k(1) == 2
                qtType = 4; k = k(1);
                exy = exy([k mNode(k) a(k+1) a(k+2) mNode(k+2) a(k+3)],:);
            else
                qtType = 3; exy = exy([4 6 1 5 2 3],:);
            end
        case 7
            qtType = 5; k = 1+find(mNode==0);
            exy = exy([a(k) mNode(a(k)) a(k+1) mNode(a(k+1)) a(k+2) mNode(a(k+2)) a(k+3)],:);
        case 8
            qtType = 6; exy = exy([ 1 5 2 6 3 7 4 8],:);
    end
    ele{ii} = 1:inode;
    eleNode(ii) = inode;
    eleCoord{ii} = exy(1:inode,:);
    eleQT(ii,1) = qtType;
    dx = exy(2,1)-exy(1,1); dy = exy(2,2)-exy(1,2);
    if abs(dx) > abs(dy)
        if dx > 0; eleQT(ii,2) = 1;  else eleQT(ii,2) = 3; end
    else
        if dy > 0; eleQT(ii,2) = 4;  else eleQT(ii,2) = 2; end %note: y is positive when pointing downward
    end
end

%% merge nodes and update element connectivity
coord = cell2mat(eleCoord);
coord(:,2) = M+1 - coord(:,2); %transform y to positive when pointing upwards
eleCentre(:,2)= M+1 - eleCentre(:,2);
coord = round(2*coord(:,:));
[coord, ~, ic] = unique(coord,'rows');
coord = coord/2;
%ele = cellfun(@(x) ic(x)',ele,'UniformOutput',false); %slower than the following loop
cumNode = 0;
eleDof = cell(nEle,1);
nNode = size(coord,1);
gDof = [1:2:2*nNode;2:2:2*nNode];
for ii = 1:nEle
    ele{ii} = ic(ele{ii}+cumNode)';
    eDof = gDof(:,ele{ii});
    eleDof{ii} = eDof(:);
    cumNode =  cumNode+eleNode(ii);
end

end

