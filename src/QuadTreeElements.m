function [ QTEle ] = QuadTreeElements( mat, g )
% compute quadtree elements for all materials
QTEle{6}.xy = [-1 -1; 0 -1; 1 -1; 1 0; 1 1; 0 1; -1 1; -1  0 ];%8-node
QTEle{1}.xy = [-1 -1; 1 -1; 1 1; -1 1];%4-node
QTEle{2}.xy = [-1 -1; 0 -1; 1 -1; 1 1; -1 1]; %5-node
QTEle{3}.xy = [-1 -1; 0 -1; 1 -1; 1 0; 1 1; -1 1]; %6-node with adjacent mid-nodes
QTEle{4}.xy = [-1 -1; 0 -1; 1 -1; 1 1; 0 1; -1 1]; %6-node with oppositing mid-nodes
QTEle{5}.xy = [-1 -1; 0 -1; 1 -1; 1 0; 1 1; 0 1; -1 1];%7-node
nMat = length(mat);
for jj = 1:6
    xy = QTEle{jj}.xy(:,1:2); nNode = size(xy,1); nDof = nNode+nNode;
    Conn = [1:nNode; 2:nNode 1];
    stff = zeros(nDof,nDof,4,nMat);
    mass = zeros(nDof,nDof,4,nMat);
    lambda = zeros(nDof,4,nMat);
    v = zeros(nDof,nDof,4,nMat);
    vinv = zeros(nDof,nDof,4,nMat);
    selfw = zeros(nDof,4,nMat);
    strnMode = zeros(3*size(Conn,2),nDof-2,4,nMat);
    sp = [ones(nNode,1) (xy+[xy(2:end,:); xy(1,:)])/2];
    QTEle{jj}.strLSFit = (sp'*sp)\sp';
    for kk = 1:4
        for ii = 1:nMat;
            [stff(:,:,kk,ii), lambda(:,kk,ii), v(:,:,kk,ii), vinv(:,:,kk,ii),...
                mass(:,:,kk,ii), selfw(:,kk,ii), strnMode(:,:,kk,ii) ] ...
                =   getSBFEMStiffMat_2NodeEle( xy, Conn, mat{ii}.D, mat{ii}.density, g);
        end
        if mat{ii}.phantom == 1; 
            stff(:,:,kk,ii) = (1.d-12)*stff(:,:,kk,ii); 
            mass(:,:,kk,ii) = (1.d-10)*mass(:,:,kk,ii); 
            selfw(:,kk,ii) = 0;
        end;

        xy = xy*[0 1; -1 0]; %rotate the element by 90deg
    end
    QTEle{jj}.stff = stff;
    QTEle{jj}.lambda = lambda;
    QTEle{jj}.v = v;
    QTEle{jj}.vinv = vinv;
    QTEle{jj}.mass = mass;
    QTEle{jj}.selfw = selfw;
    QTEle{jj}.strnMode = strnMode;
end


end

