function pctOverUnity = tangleScatter( tangle_1, tangle_2, labels, varargin )
% marker = 'o'; color = 'k';
% plotDataOnly = false;

% myTitle = sprintf('%s  vs  %s', data1{:}, data2{:});
% if nargin == 1
%    myTitle = '';
% elseif nargin == 2
%    myTitle = varargin{1};
% elseif nargin >= 3 
%    myTitle = varargin{1};
%    color = varargin{2};
% end
params = inputParser;

% Set default options
params.addParameter('marker',  'o');
params.addParameter('color', 'k');
params.addParameter('plotDataOnly', false, @islogical);


% Parse command line for parameter settings
params.parse(varargin{:});

marker = params.Results.marker;
color = params.Results.color;
plotDataOnly = params.Results.plotDataOnly;

if ~plotDataOnly
   figure;
end
%%
% if any(strcmp(varargin,'plotDataOnly')); plotDataOnly = true; end

% scale = 1e3;
% scale = 1e3;


 gray = [.2 .2 .2]; 
% D1 = P2D( P, data1 );


[h,p] = ttest(tangle_1, tangle_2);
scale = floor(log10(max([tangle_2; tangle_1])));
%* (10^(scale-1));
% step = 5 * (10^(scale-1));
% for simplicity :
% tangle_1 = tangle_1./(10^scale); tangle_2 = tangle_2./(10^scale);
Q90_m1 = prctile(tangle_1,90); Q90_emg = prctile(tangle_2,90);

% if ~plotDataOnly; blankFigure([-1 1 -1 1]); end

maxAxis = max([tangle_2; tangle_1]); % get max point
% maxAxis = 1000;
step = 10^(scale-1);
maxAxis = ceil(maxAxis/step)*step; 

plot(tangle_2, tangle_1, marker,'MarkerFaceColor', color,'MarkerEdgeColor', color,'MarkerSize',4);  hold on;
plot([Q90_emg Q90_emg], [-.01*maxAxis .01*maxAxis],'Color',[1 .5 0],'LineWidth',3)
plot( [-.01*maxAxis .01*maxAxis], [Q90_m1 Q90_m1],'Color',color,'LineWidth',3)

% general figure prettying
if ~plotDataOnly
h = gca;
daspect([1 1 1])

% text(maxAxis/2, maxAxis*.95, myTitle, 'FontSize', 24,'Color',gray, 'HorizontalAl','center')


% Axis Params
genParams = struct();
genParams.fontSize = 18;
genParams.axisOffset = -.01*maxAxis;
genParams.color = gray;
genParams.lineThickness = 2;
genParams.tickLocations = [0 maxAxis];
genParams.longTicks = [0 maxAxis];
genParams.tickLabels = arrayfun(@(n) {num2str(n)}, genParams.longTicks);
genParams.tickLabelLocations = genParams.longTicks;
% genParams.extraLength = 0.75;
genParams.tickLength = maxAxis/95;

% X Axis
xParams = genParams;
xParams.axisLabel = labels{2};
% xParams.axisLabel = sprintf('%s (x 10 %1g)',labels{2}, scale);
AxisMMC(0, maxAxis, xParams);

% Y Axis
yParams = genParams;
yParams.axisLabel = labels{1};
% yParams.axisLabel = sprintf('%s (x 10 %1g)',labels{1}, scale);
yParams.axisOrientation = 'v';
AxisMMC(0, maxAxis, yParams);

% axis([-1 maxAxis*1.25 -1 maxAxis*1.25])

% 
% % X axis
% % set(h, 'XLim',[0 maxAxis], 'XTick',[0:step*2:maxAxis]);%,'XTickLabel',[0:step*2:maxAxis])
% h.XAxis.MinorTickValues = [0:step:maxAxis];
% h.XAxis.MinorTick = 'on';
% h.XAxis.Color = gray;
% h.XAxis.Exponent = 3;
% xlabel('Q_{emg}', 'interpreter', 'tex','FontSize',20);
% 
% % Y axis
% set(h, 'YLim',[0 maxAxis], 'YTick',[0:step*2:maxAxis]);%,'YTickLabel',[0:step*2:maxAxis])
% h.YAxis.MinorTickValues = [0:step:maxAxis];
% h.YAxis.MinorTick = 'on';
% h.YAxis.Color = gray;
% h.YAxis.Exponent = 3;
% ylabel('Q_{neural}', 'interpreter', 'tex','FontSize',20);


   
plot([0 maxAxis],[0 maxAxis],'b','LineWidth',1.5)

end
pctOverUnity = sum((tangle_1./tangle_2)>1)/length(tangle_1);
pctOverUnity = pctOverUnity*100;
% text(step/4, maxAxis/2, sprintf('%1.1f%% points > 1', pctOverUnity),'FontSize',14,'Color', gray);
set(gca,'Visible','off')
set(gcf,'Color','w')
end

