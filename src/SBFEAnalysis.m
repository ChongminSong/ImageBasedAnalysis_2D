function [U, varargout] = SBFEAnalysis(cntrl, QTEle, coord, ele, eleQT, eleSize, eleDof, ...
    eleMat, mat, bc_force, bc_disp)


[ K, rhs ] = AssembleQuadtreeMesh(QTEle,...
    coord, ele,  eleQT, eleSize, eleDof, eleMat);

if strncmpi(cntrl.prbtdm.type, 'statics', 3)
    [U, varargout{1}] = static_solver(2, bc_force, bc_disp, K, rhs);
    % elseif strncmpi(cntrl.prbtdm.type, 'modal', 3)
    %     [varargout{1} varargout{2}] = modes_solver(cntrl.prbtdm.modalPara(1), ndn, bc_disp, K, M);
elseif strncmpi(cntrl.prbtdm.type, 'time', 3)
    [varargout{1}] = tm_newmark(cntrl.prbtdm.TIMEPara,cntrl.prbtdm.forceHistory ...
        , ndn, bc_force, bc_disp, K, M, F);
else
end

end
