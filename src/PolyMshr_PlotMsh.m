function PolyMshr_PlotMsh(coord,ele, Supp, Load)
hold on;
ne  = length(ele);
elePoint = NaN(ne,8);
for ii = 1:ne
    nodes = ele{ii};
    elePoint(ii,1:length(nodes)) = nodes;
end
patch('Faces',elePoint,'Vertices',coord,'FaceColor','w');
if exist('Supp','var')&&~isempty(Supp)&&~isempty(Load)%Plot BC if specified
    plot(coord(Supp(:,1),1),coord(Supp(:,1),2),'b>','MarkerSize',8);
    plot(coord(Load(:,1),1),coord(Load(:,1),2),'m^','MarkerSize',8); hold off;
end
end

