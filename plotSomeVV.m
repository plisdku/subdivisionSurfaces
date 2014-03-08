function plotSomeVV(VV, vertices, selected, varargin)
% plotSomeVV(VV, vertices, selected, plotStyle...)
% plotSomeVV(VV, vertices, selected, radius, plotStyle...)
%

if nargin > 3 && isnumeric(varargin{1})
    % plotSomeVV(VV, vertices, selected, neighborhood, varargin)
    
    neighborhood = varargin{1};
    
    A = vv2adjacency(VV);
    
    flags = zeros(size(VV,1),1);
    flags(selected) = 1;
    
    neighborhood = (A^neighborhood)*flags | flags;
    
    [VV vertices] = truncateVV(VV, vertices, neighborhood);
    
    plotOpts = varargin(2:end);
else
    [VV, vertices] = truncateVV(VV, vertices, selected);
    
    plotOpts = varargin(1:end);
end

plotVV(VV, vertices, plotOpts{:});