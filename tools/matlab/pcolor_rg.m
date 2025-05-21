%==========================================================================
% Function: pcolor_rg
%
% Description:
%   Checkerboard plot of matrix C using pcolor with flat shading.
%   All elements of C are used.
%
% Input:
%   X - x vector
%   Y - y vector
%   C - matrix to be plotted
%   V - two-element vector [cmin cmax]
%
% Output:
%   H - returned handle
%
% Author: Ralf Greve
% Date:   2018-12-26
%==========================================================================

function H = pcolor_rg(X, Y, C, V)

SX = size(X); if SX(2) == 1; X = X'; end   % ensure X is row vector
SY = size(Y); if SY(2) == 1; Y = Y'; end   % ensure Y is row vector

dX = X(2)-X(1) ;
XX = [X-0.5*dX X(end)+0.5*dX] ;

dY = Y(2)-Y(1) ;
YY = [Y-0.5*dY Y(end)+0.5*dY] ;

CC = [[C nan*zeros(size(C,1),1)] ; nan*zeros(1,size(C,2)+1)] ;

H = pcolor(XX, YY, CC); shading('flat')
caxis(V)

end % function pcolor_rg

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
