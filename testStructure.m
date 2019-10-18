%function [ ] = test ()
close all
clearvars
dbstop if error

matD = @(E, p) E/(1-p^2)*[1 p 0; p 1 0;0 0 (1-p)/2];
imageFolder = 'images';
myImread = @(x) imread([imageFolder,'/',x]);

%% set path
sourceFolder = 'D:\sbfem_matlab\src\solution';
addpath(genpath(sourceFolder));
sourceFolder = 'D:\sbfem_matlab\imagQuadtree\src';
addpath(genpath(sourceFolder));

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

ImgOrg = myImread('mytunnel.bmp');%
%ImgOrg = myImread('submodel.bmp');

% ImgOrg = myImread('a.bmp');
% ImgOrg = myImread('canvas.bmp');
% ImgOrg = myImread('bigImage8192.bmp');
% ImgOrg = myImread('hugeImage16384.bmp');
% ImgOrg = myImread('multiscale1.bmp');


%ImgOrg = myImread('multiscale1.bmp');
% Img = uint8(round(sum(ImgOrg,3)/size(ImgOrg,3)));

% ImgOrg = myImread('cirHole.bmp');
% 
%  n = 1024; ImgOrg = zeros(n,'uint8'); 
%  [ImgOrg, nc] = FillCircle(ImgOrg, 150, 50, 10, 0.1);%FillCircle(ImgOrg, icolor, d, ncir, gap)
% 

%n = 1024; ImgOrg = zeros(n,'uint8'); % d=200; [r,c]=imgCircle(d); ImgOrg(n/2+r+(n/2+c-1)*n) = 200;

%n = 1024*1; ImgOrg = zeros(n,'uint8'); [r,c]= imgRectangle( 400, 100, 0); ImgOrg(n/2+r+(n/2+c-1)*n) = 200; 

%%[ImgOrg] = ScatteredEllipse(ImgOrg, icolor, dmax, dmin, dmid, ncir, dratioMin)
% n=1024*4; ImgOrg = zeros(n,'uint8'); ImgOrg = ScatteredEllipse(ImgOrg, 150, 400, 10, 50, 200, 1/4);

%%ImgOrg = ScatteredCrack(ImgOrg, icolor, dmax, dmin, thickness, ncrack)
%n=1024*4; ImgOrg = zeros(n,'uint8');ImgOrg =  ScatteredCrack(ImgOrg, 230, 400, 100, 5, 100);

%ImgOrg = myImread('ZhenJunConcrete1.bmp');

% ImgOrg = myImread('working.bmp');
% ImgOrg = myImread('test.bmp');

% Img( Img<=50 ) = 0; Img( Img>50 & Img<200 ) = 100; Img( Img >= 200& Img < 250) = 240;
Img = uint8(round(sum(ImgOrg,3)/size(ImgOrg,3)))-2;

% imwrite(Img,[imageFolder,'/','multiscale1.bmp'],'bmp')
QTthreshold = 0.15;
Resolution = 1;
SBFEMAnalysis = 1;
matColor = [0 50 1; 51 230 2; 231 255 3]; voidColor = [231 255];
mat{1}.D = matD(10,0.2); mat{1}.density = 0; mat{1}.phantom = 0;
mat{2}.D = matD(100,0.2); mat{2}.density = 1; mat{2}.phantom = 0;
mat{3}.D = matD(1,0.2); mat{3}.density = 0; mat{3}.phantom = 0;
g = -9.81;
 
nImg = size(Img);
n = min(1024*16,2^nextpow2(max(nImg)));

I = repmat( feval(class(Img),255),n);
I(1:min(n,nImg(1)), 1:min(n,nImg(2))) = Img(1:min(n,nImg(1)), 1:min(n,nImg(2)));
%I(1:1024*4, 1:1024*8) = Img(1:1024*4, 1:1024*8);

% figure(1)
% %subplot(121)
% imshow(ImgOrg)
% clear ImgOrg;
figure
% subplot(122)
 imshow(Img)
clear Img;

disp(['*** Input image completed. time = ',num2str(toc(ticImage))]);
disp(['    Image size = ',num2str(nImg(1:2))]);

ticMesh = tic;

minDim = 2^1;
maxDim = n/32;

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
Tolerance = 0.2*minEleSize;

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

    dNode = length(leftNodes)+length(rightNodes)+length(bottomNodes);
    bc_disp = zeros(dNode,3);
    bc_disp(:,1) = [bottomNodes; leftNodes; rightNodes];
    bc_disp(1:length(bottomNodes),2) = 2; 
    bc_disp(length(bottomNodes)+1:end,2) = 1;
%     dNode = length(bottomNodes)+1;
%     bc_disp = zeros(dNode,3);
%     bc_disp(:,1) = [bottomNodes; bottomNodes(1)];
%     bc_disp(:,2) = 2; bc_disp(end,2) = 1;


    topEdges = find(abs(QTedgeCentre(:,2)-max(coord(:,2)))<Tolerance);
    topEdgeLen = abs(sum(coord(QTedge(topEdges,2),:)-coord(QTedge(topEdges,1),:),2));
    Load = zeros(nNode,3);
    Load(:,1) = 1:nNode;
    bc_force = zeros(2*nNode,3); 
    bc_force(:,1) = reshape([1:nNode;1:nNode],[],1);
    bc_force(:,2) = reshape([ones(1,nNode);2*ones(1,nNode)],[],1);
    minx = min(min(coord(QTedge(topEdges,2),1)), min(coord(QTedge(topEdges,1))));
    maxx = max(max(coord(QTedge(topEdges,2),1)), max(coord(QTedge(topEdges,1))));
    x = QTedgeCentre(topEdges,1);
%    pressure = 2*(x-minx)/(maxx-minx)-1;
    pressure = 1;
    Load(QTedge(topEdges,1),3) = Load(QTedge(topEdges,1),3) + pressure.*topEdgeLen/2;
    Load(QTedge(topEdges,2),3) = Load(QTedge(topEdges,2),3) + pressure.*topEdgeLen/2;
    Load = Load(Load(:,3)~=0,:);
    bc_force(QTedge(topEdges,1)*2,3) = bc_force(QTedge(topEdges,1)*2,3)+ pressure.*topEdgeLen/2;
    bc_force(QTedge(topEdges,2)*2,3) = bc_force(QTedge(topEdges,2)*2,3)+ pressure.*topEdgeLen/2;
    disp(['*** Enforcing boundary condition completed. time = ',num2str(toc(ticBC))]);
end


ticPlotmesh = tic;
figure
%subplot(122)
axis equal; 
% if exist('BCDisp','var')
%     PolyMshr_PlotMsh(coord, ele, BCDisp, Load);
% else
    PolyMshr_PlotMsh(coord, ele);
% end
axis on; % axis tight;

% figure
% a = find(abs(eleCentre(:,1)-200)<120&abs(eleCentre(:,2)-300)<50);
% PolyMshr_PlotMsh(coord, ele(a));
% axis equal
 
% % figure
%  PlotElement(coord, ele, eleColor, gray);
%  axis on
disp(['*** Plotting mesh completed. time = ',num2str(toc(ticPlotmesh))]);

if ~SBFEMAnalysis
    return
end

%%
ticSolution = tic;

QTEle  = QuadTreeElements( mat, g );

% [U, QTEle, ReFrc] = SBFEAnalysis(coord, ele, eleQT, eleSize, eleDof, eleMat, ...
%     mat, g, bc_force, bc_disp);
[U, ReFrc] = SBFEAnalysis(cntrl, QTEle, coord, ele, eleQT, eleSize, eleDof, eleMat, ...
    mat, bc_force, bc_disp);
ReFrc = reshape(ReFrc,2,[]); sum(ReFrc(2,topNodes));
disp(['*** Solution completed. time = : ',num2str(toc(ticSolution))]);

ticPost = tic;

U = reshape(U,2,[]);
[ eleStrs, eleResult, eleStrsNode] = ElementStress(U, QTEle, coord, ele, eleQT, eleSize, eleMat, mat );

maxU = max(abs(U(:,[ele{eleTrue}])), [], 2);
maxU = max(maxU);
maxD = max(max(coord)-min(coord));
%deformed = 0.05*maxD/maxU*U' + coord;
deformed = 0.0003*U' + coord;

figure
PolyMshr_PlotMsh(deformed, ele(eleTrue));
axis equal; axis on

% figure
% PlotElement(coord, ele(eleTrue), eleStrs(5,eleTrue)', jet);
% caxis([-20000 2000]);
% colorbar
% %axis on
% 
% figure
% myjet = jet; myjet = myjet(1:4:end,:);
% PlotResult(coord, U, ele, eleResult, eleStrsNode,  eleMat, eleCentre, ...
%     'Component', 'SYY', 'Element',eleTrue, 'Material', ([1:2]), ...
%     'Average', 'YES', 'ColorMap', myjet, 'DeformFactor',0.0001);
% %
% colorbar
% 

% figure
% myjet = jet; myjet = myjet(1:6:end,:);
% PlotResult(coord, U, ele, eleResult, eleStrsNode,  eleMat, eleCentre, ...
%     'Component', 'SP1', 'Element',eleTrue, 'Material', ([1:2 ]), ...
%     'Average', 'YES', 'ColorMap', myjet, 'DeformFactor',0.05);
% caxis([0 cmax])
% colorbar
% 
% 

disp(['*** Post-processing completed. time = : ',num2str(toc(ticPost))]);

disp(['Total time = : ',num2str(toc(ticTotal))]);

fclose('all');
