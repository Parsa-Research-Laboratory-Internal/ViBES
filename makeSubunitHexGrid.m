function [posX, posY] = makeSubunitHexGrid(stimSize, subunit_spacing)
maxDim = max(stimSize);
gridSize = ceil((maxDim./subunit_spacing));
if rem(gridSize,2) == 1, gridSize = gridSize+1; end

Rad3Over2 = sqrt(3) / 2;
[X, Y] = meshgrid(1:1:gridSize);
n = size(X,1);
X = Rad3Over2 * X;
Y = Y + repmat([0 0.5],[n,n/2]);

%set spacing
X = X * subunit_spacing;
Y = Y * subunit_spacing;

Ind = (X <= stimSize(2) & Y <= stimSize(1));
posX = X(Ind);
posY = Y(Ind);

% Ind = (X-subunit_spacing*4 <= stimSize(1) & Y-subunit_spacing*4 <= stimSize(2));
% posX = X(Ind) + subunit_spacing*2;
% posY = Y(Ind) + subunit_spacing*2;