function [ K, rhs, M ] = AssembleQuadtreeMesh(QTEle, coord, ele, ...
    eleQT, eleSize, eleDof, eleMat )
% assemble quadtree mesh

nNode = size(coord,1);
rhs = zeros(2*nNode,1);  %external load vector
eNDof = 2*cellfun(@length,ele); % # of DOFs per element
femi = zeros(sum(eNDof.^2),1); femj=femi; femk=femi; femm = femi;
index = 0;
%for el = 1:length(ele)
nMat = max(eleMat);
for jj = 1:6
    cEle = find(eleQT(:,1)==jj);
    for ii = 1:nMat
        mcEle = cEle(eleMat(cEle)==ii);
        mcQTEleStff = QTEle{jj}.stff(:,:,:,ii);
        mcQTEleMass = QTEle{jj}.mass(:,:,:,ii);
        mcQTEleSelfw = QTEle{jj}.selfw(:,:,ii);

        for el = mcEle'
            NDof = eNDof(el);
            eDof = eleDof{el};
            I = eDof(:,ones(1,NDof)); J=I';
            ncoe = NDof*NDof;
            kk = eleQT(el,2);
            femi(index+1:index+ncoe) = I(:);
            femj(index+1:index+ncoe) = J(:);
            femk(index+1:index+ncoe) = mcQTEleStff(:,:,kk);
            s = (eleSize(el)/2)^2;
            femm(index+1:index+ncoe) = s*mcQTEleMass(:,:,kk);
            rhs(eDof) = rhs(eDof) + s*mcQTEleSelfw(:,kk);

            index = index + ncoe;
        end
    end
end

K = sparse(femi,femj, femk);
K = (K+K')/2;

if nargout >2
    M = sparse(femi,femj, femm);
    M = (M+M')/2;
end

end

