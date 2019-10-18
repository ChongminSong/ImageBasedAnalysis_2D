%function [ ] = test ()
close all
clearvars
dbstop if error

%% set path
% restoredefaultpath;
sourceFolder = './src';
addpath(genpath(sourceFolder));
imageFolder = './images';
myImread = @(x) imread([imageFolder,'/',x]);

matD = @(E, p) E/(1-p^2)*[1 p 0; p 1 0;0 0 (1-p)/2];


ticTotal = tic;
ticImage = tic;

cntrl.phyprb = 'ELASTICITY';    % type of physical problem: 'ELASTICITY'; 'DIFFUSION'
%cntrl.typfunc = 'SCALAR';  % type of function: 'VECTOR'; 'SCALAR'
cntrl.prbtgeo.type='TWO-DIM';  %'AXISYMMETRIC'; %'TWO-DIM'
cntrl.prbtdm.type='statics';  %'FREQUENCY'; %'STATICS':'TIME'; 'MODAL';
cntrl.prbtdm.modalPara = [120 200];
%cntrl.prbtdm.frqPara = [100 0.05 0.1];  
cntrl.prbtdm.TIMEPara = [1000 0.025 0.5 0.25];
cntrl.ForceHistory = [ 0 0; 0.2d0 1.d0; 0.4d0 0.d0;  900 0.d0];

QTthreshold = 0.05;
Resolution = 1;
SBFEMAnalysis = 1;

% n = 128*8; ImgOrg = zeros(n,'uint8'); d=40*6; [r,c]=imgCircle(d); ImgOrg(n/2+r+(n/2+c-1)*n) = 200;

%  n = 1024; ImgOrg = zeros(n,'uint8'); 
% [ImgOrg, nc] = FillCircle(ImgOrg, 155, 150, 20, 0.1);%FillCircle(ImgOrg, icolor, d, ncir, gap)
% % [ImgOrg, nc] = FillCircle(ImgOrg, 255, 100, 3, 0.01);%FillCircle(ImgOrg, icolor, d, ncir, gap)
% minDim = 2^0;
% matColor = [0 50 1; 51 220 2; 221 255 3]; voidColor = [254 255];
% mat{1}.D = matD(70,0.2); mat{1}.density = 0.01; mat{1}.phantom = 0;
% mat{2}.D = matD(25,0.2); mat{2}.density = 0.01; mat{2}.phantom = 0;
% mat{3}.D = matD(0.01,0.2); mat{3}.density = 0; mat{3}.phantom = 1;


% % ImgOrg = myImread('a123.bmp');ImgOrg(ImgOrg>250)=200; ImgOrg(:,:,1)=ImgOrg(:,:,3);ImgOrg(:,:,2)=ImgOrg(:,:,3);
ImgOrg = myImread('concreteBW512.bmp');
% ImgOrg = myImread('ZhenJunConcrete1.bmp');
minDim = 2^0;
matColor = [0 50 1; 51 230 2; 231 255 3]; voidColor = [230 255];
mat{1}.D = matD(70,0.2); mat{1}.density = 0.01; mat{1}.phantom = 0;
mat{2}.D = matD(25,0.2); mat{2}.density = 0.01; mat{2}.phantom = 0;
mat{3}.D = matD(11,0.2); mat{3}.density = 0; mat{3}.phantom = 0;


% ImgOrg = myImread('cementpaste.bmp');
% matColor = [0 50 1; 51 200 2; 201 251 3; 252 255 4]; voidColor = [252 255];
% mat{1}.D = matD(80,0.2); mat{1}.density = 0.01; mat{1}.phantom = 0;
% mat{2}.D = matD(40,0.2); mat{2}.density = 0.01; mat{2}.phantom = 0;
% mat{3}.D = matD(20,0.2); mat{3}.density = 0.01; mat{3}.phantom = 0;
% mat{4}.D = matD(1,0.2); mat{4}.density = 0; mat{4}.phantom = 1;


% ImgOrg = myImread('eng_z_2.jpg');
% ImgOrg = myImread('ZhenJunConcrete1.bmp');
% matColor = [0 50 1; 51 220 2; 221 255 3]; voidColor = [240 255];
% mat{1}.D = matD(70,0.2); mat{1}.density = 0.01; mat{1}.phantom = 0;
% mat{2}.D = matD(25,0.2); mat{2}.density = 0.01; mat{2}.phantom = 0;
% mat{3}.D = matD(0.1,0.2); mat{3}.density = 0; mat{3}.phantom = 1;

%  ImgOrg = myImread('tabinterior0249.bmp');
% ImgOrg = myImread('myTab0249.bmp'); ImgOrg = 100*ImgOrg;
%  ImgOrg = myImread(['concrete3D\tmp\',num2str(1),'.jpg']);
% ImgOrg = myImread('ceramiccomposite.bmp');
% ImgOrg = myImread('micro512.jpg');
% ImgOrg = myImread('EX2514C14c.bmp'); %another concrete

% ImgOrg = myImread('multiscale1.bmp');
% Img = uint8(round(sum(ImgOrg,3)/size(ImgOrg,3)));

% ImgOrg = myImread('cirHole.bmp');

% n = 2*1024; ImgOrg = zeros(n,'uint8'); 
% d = 200; nc = 10; gap = 0.1;%
% [ImgOrg, nc] = FillCircle(ImgOrg, 150, d, nc, gap);
% VolumeFraction = nc*0.25*pi*d^2/n^2


% ImgOrg = myImread('test.bmp');


%ImgOrg = myImread('Concrete_Yvonnet_IJNME2012.bmp');

%ImgOrg = myImread('Metallic–ceramic_kikuchi.bmp');


 Img = uint8(round(sum(ImgOrg,3)/size(ImgOrg,3)));
% I = imresize(Img,2,'nearest');
% c=contour(seg,3,'r'); hold on;
% ib = 1;
% nend = size(c,2);
% while ib<nend
% n = c(2,ib);
% xy = c(1:2,ib+1:ib+n)';
% plot(xy(:,1),xy(:,2),'b');
% sx = smooth(xy(:,1));
% sy = smooth(xy(:,2));
% plot(sx,sy,'r');
% %%in = inpolygon(xq,yq,sx,sy);
% % dxy = xy(2:end,:)-xy(1:end-1,:);
% % s=[0;cumsum(sqrt(sum(dxy.^2,2)))];
% % fx=fit(s,xy(:,1),'smoothingspline');
% % fy=fit(s,xy(:,2),'smoothingspline');
% % plot(fx(s),fy(s))
% ib = ib+n+1;
% end
% % imwrite(Img,[imageFolder,'/','multiscale1.bmp'],'bmp')
% matColor = [0 50 1; 51 220 2; 221 255 3]; voidColor = [240 255];
% mat{1}.D = matD(70,0.2); mat{1}.density = 0.01; mat{1}.phantom = 0;
% mat{2}.D = matD(25,0.2); mat{2}.density = 0.01; mat{2}.phantom = 0;
% mat{3}.D = matD(0.1,0.2); mat{3}.density = 0; mat{3}.phantom = 1;
% matColor = [0 70 1; 71 179 2; 180 250 3; 251 255 4]; voidColor = [250 255];
% mat{1}.D = matD(70,0.2); mat{1}.density = 0.01; mat{1}.phantom = 0;
% mat{2}.D = matD(25,0.2); mat{2}.density = 0.01; mat{2}.phantom = 0;
% mat{3}.D = matD(55,0.2); mat{3}.density = 0.01; mat{3}.phantom = 0;
% mat{4}.D = matD(0.1,0.2); mat{4}.density = 0; mat{4}.phantom = 0;
g = 0;% -9.81;
 
nImg = size(Img);
n = min(1024*512,2^nextpow2(max(nImg)));

I = repmat( feval(class(Img),255),n);
I(1:min(n,nImg(1)), 1:min(n,nImg(2))) = Img(1:min(n,nImg(1)), 1:min(n,nImg(2)));

% figure(1)
% %subplot(121)
% imshow(ImgOrg)
% clear ImgOrg;
figure
% subplot(122)
 imshow(Img)
%clear Img;

disp(['*** Input image completed. time = ',num2str(toc(ticImage))]);
disp(['    Image size = ',num2str(nImg(1:2))]);

ticMesh = tic;

maxDim = n/4;

S = qtdecomp(I, QTthreshold, [minDim, maxDim]);
[coord, ele, eleQT, eleColor, eleSize, eleCentre, eleDof] = ...
    quadTreeMesh(S,I,voidColor);
nNode = size(coord,1);
nEle = length(ele);

%%
if SBFEMAnalysis
    nMat = length(mat);
    eleMat = zeros(nEle,1);
    eleColor = round(eleColor);
    phantomMat = cellfun(@(x) x.phantom, mat);
    phantomMat = find(phantomMat);
    for ii =  1:nMat
        eleMat( eleColor>=matColor(ii,1) & eleColor<=matColor(ii,2) ) = ii;
    end
    vele = find( eleMat==0 );
    if ~isempty(vele)
        disp('elements without assigned material exist')
        pause;
    end
end
if exist('Resolution','var'); 
    coord = Resolution*coord; 
    eleSize = Resolution*eleSize; 
    eleCentre = Resolution*eleCentre; 
end;
minEleSize = min(eleSize);
Tolerance = 0.5*minEleSize;

disp(['*** Mesh generation completed. time = ',num2str(toc(ticMesh))]);
disp(['    Number of nodes = ',num2str(nNode)]);
disp(['    Number of elements = ',num2str(nEle)]);
disp(['    Ratio of number of elements to number of pixel = ',num2str(nEle/(nImg(1)*nImg(2)))]);
%%
if SBFEMAnalysis
    ticBC = tic;
    
    eleTrue = 1:nEle;
    if exist('phantomMat','var')
        if ~isempty(phantomMat)
            eleTrue(eleMat==phantomMat)=[];
        end
    end

    [ QTedge, QTedgeCentre, eleEdge, edge2Ele ] = findElementEdges( ele(eleTrue), coord );
    
    bottomNodes = find(abs(coord(:,2)-min(coord(:,2)))<Tolerance);
    topNodes = find(abs(coord(:,2)-max(coord(:,2)))<Tolerance);
    leftNodes = find(abs(coord(:,1)-min(coord(:,1)))<Tolerance);
    rightNodes = find(abs(coord(:,1)-max(coord(:,1)))<Tolerance);

% uniform extension
    dNode = length(topNodes)+length(bottomNodes)+1;
    bc_disp = zeros(dNode,3);
    bc_disp(:,1) = [topNodes; bottomNodes; bottomNodes(1)];
    bc_disp(:,2) = 2; bc_disp(end,2) = 1;
    bc_disp(1:length(topNodes),3) = 0.01*nImg(1);
% %uniform tension
%     dNode = length(bottomNodes)+1;
%     bc_disp = zeros(dNode,3);
%     bc_disp(:,1) = [bottomNodes; bottomNodes(1)];
%     bc_disp(:,2) = 2; bc_disp(end,2) = 1;


    topEdges = find(abs(QTedgeCentre(:,2)-max(coord(:,2)))<Tolerance);
    topEdgeLen = abs(sum(coord(QTedge(topEdges,2),:)-coord(QTedge(topEdges,1),:),2));
    bc_force = zeros(2*nNode,3); 
    bc_force(:,1) = reshape([1:nNode;1:nNode],[],1);
    bc_force(:,2) = reshape([ones(1,nNode);2*ones(1,nNode)],[],1);
    minx = min(min(coord(QTedge(topEdges,2),1)), min(coord(QTedge(topEdges,1))));
    maxx = max(max(coord(QTedge(topEdges,2),1)), max(coord(QTedge(topEdges,1))));
    x = QTedgeCentre(topEdges,1);
%    pressure = 2*(x-minx)/(maxx-minx)-1;
    pressure = 0;
    bc_force(QTedge(topEdges,1)*2,3) = bc_force(QTedge(topEdges,1)*2,3)+ pressure.*topEdgeLen/2;
    bc_force(QTedge(topEdges,2)*2,3) = bc_force(QTedge(topEdges,2)*2,3)+ pressure.*topEdgeLen/2;
    disp(['*** Enforcing boundary condition completed. time = ',num2str(toc(ticBC))]);
end


ticPlotmesh = tic;
figure
axis equal; 
PolyMshr_PlotMsh(coord, ele);
axis on; % axis tight;


% figure('Color','white');
% %PolyMshr_PlotMsh(coord, ele);
% %PlotElement(coord, ele, eleColor, newmap); 
% PlotElement(coord, ele, eleColor); 
% caxis([0 255])

% hold on
% ib = 1;
% nend = size(c,2);
% while ib<nend
% n = c(2,ib);
% xy = c(1:2,ib+1:ib+n)';
% plot(xy(:,1),xy(:,2),'b');
% sx = smooth(xy(:,1));
% sy = smooth(xy(:,2));
% plot(sx,512-sy,'r');
% % dxy = xy(2:end,:)-xy(1:end-1,:);
% % s=[0;cumsum(sqrt(sum(dxy.^2,2)))];
% % fx=fit(s,xy(:,1),'smoothingspline');
% % fy=fit(s,xy(:,2),'smoothingspline');
% % plot(fx(s),fy(s))
% ib = ib+n+1;
% end
% 
% 

disp(['*** Plotting mesh completed. time = ',num2str(toc(ticPlotmesh))]);

if ~SBFEMAnalysis
    return
end

%%
ticSolution = tic;

QTEle  = QuadTreeElements( mat, g );

[U, ReFrc] = SBFEAnalysis(cntrl, QTEle, coord, ele, eleQT, eleSize, eleDof, eleMat, ...
    mat, bc_force, bc_disp);
ReFrc = reshape(ReFrc,2,[]); 
sum(ReFrc(2,topNodes))
disp(['*** Solution completed. time = : ',num2str(toc(ticSolution))]);

ticPost = tic;

U = reshape(U,2,[]);
[ eleStrs, eleResult, eleStrsNode] = ElementStress(U, QTEle, coord, ele, eleQT, eleSize, eleMat, mat );

maxU = max(abs(U(:,[ele{eleTrue}])), [], 2);
maxU = max(maxU);
maxD = max(max(coord)-min(coord));
deformed = 0.06*maxD/maxU*U' + coord;

figure
PolyMshr_PlotMsh(deformed, ele(eleTrue));
axis equal; axis on

cmax = 0.6;

figure('Color','white')
myjet = jet; myjet = myjet(1:6:end,:);
PlotResult(coord, U, ele, eleResult, eleStrsNode,  eleMat, eleCentre, ...
    'Component', 'SP1', 'Element',eleTrue, 'Material', ([1]), ...
    'Average', 'YES', 'ColorMap', myjet, 'DeformFactor',0.05);
caxis([0 cmax])
colorbar

figure('Color','white')
myjet = jet; myjet = myjet(1:6:end,:);
PlotResult(coord, U, ele, eleResult, eleStrsNode,  eleMat, eleCentre, ...
    'Component', 'SP1', 'Element',eleTrue, 'Material', ([2]), ...
    'Average', 'YES', 'ColorMap', myjet, 'DeformFactor',0.05);
caxis([0 cmax])
colorbar

figure('Color','white')
myjet = jet; myjet = myjet(1:6:end,:);
PlotResult(coord, U, ele, eleResult, eleStrsNode,  eleMat, eleCentre, ...
    'Component', 'SP1', 'Element',eleTrue, 'Material', ([1:2]), ...
    'Average', 'YES', 'ColorMap', myjet, 'DeformFactor',0.05);
caxis([0 cmax])
colorbar
% 
disp(['*** Post-processing completed. time = : ',num2str(toc(ticPost))]);

disp(['Total time = : ',num2str(toc(ticTotal))]);

