%==========================================================================
% Function: frac_area
%
% Description:
%   Computes the fraction of the area of a small rectangle (xx,yy)
%   within a larger rectangle (x,y).
%
% Author: Ralf Greve
% Date:   2024-03-13
%==========================================================================

function frac_area_r = frac_area(x, y, xx, yy)

%-------- Ensure correct order --------

if x(1)>x(2); x=flip(x); end
if y(1)>y(2); y=flip(y); end
if xx(1)>xx(2); xx=flip(xx); end
if yy(1)>yy(2); yy=flip(yy); end

%-------- Mapping of large rectangle to unit square --------

% xm = [0 1];   % large rectangle -> unit rectangle (values not needed)
% ym = [0 1];   % large rectangle -> unit rectangle (values not needed)

%-------- Mapping of small rectangle --------

diff_x_inv = 1.0/(x(2)-x(1));
diff_y_inv = 1.0/(y(2)-y(1));

xxm = [xx(1)-x(1) xx(2)-x(1)] * diff_x_inv;
yym = [yy(1)-y(1) yy(2)-y(1)] * diff_y_inv;

if xxm(2)-xxm(1) > 1 || yym(2)-yym(1) > 1
    msg = 'frac_area: Small rectangle larger than large rectangle!';
    error(msg)
end

%-------- Computation of the fraction of area --------

if (xxm(1)>=0 && xxm(2)<=1) && (yym(1)>=0 && yym(2)<=1)
    frac_area_r = 1.0;
    return
elseif xxm(1)>=0 && xxm(2)<=1
    if yym(1)<0 && yym(2)>=0
        frac_area_r = yym(2)/(yym(2)-yym(1));
        return
    elseif yym(1)<=1 && yym(2)>1
        frac_area_r = (1-yym(1))/(yym(2)-yym(1));
        return
    else
        frac_area_r = 0.0;
        return
    end
elseif yym(1)>=0 && yym(2)<=1
    if xxm(1)<0 && xxm(2)>=0
        frac_area_r = xxm(2)/(xxm(2)-xxm(1));
        return
    elseif xxm(1)<=1 && xxm(2)>1
        frac_area_r = (1-xxm(1))/(xxm(2)-xxm(1));
        return
    else
        frac_area_r = 0.0;
        return
    end
elseif (xxm(1)<0 && xxm(2)>=0) && (yym(1)<0 && yym(2)>=0)
        frac_area_r = xxm(2)*yym(2) ...
                      /((xxm(2)-xxm(1))*(yym(2)-yym(1)));
        return
elseif (xxm(1)<=1 && xxm(2)>1) && (yym(1)<0 && yym(2)>=0)
        frac_area_r = (1-xxm(1))*yym(2) ...
                      /((xxm(2)-xxm(1))*(yym(2)-yym(1)));
        return
elseif (xxm(1)<0 && xxm(2)>=0) && (yym(1)<=1 && yym(2)>1)
        frac_area_r = xxm(2)*(1-yym(1)) ...
                      /((xxm(2)-xxm(1))*(yym(2)-yym(1)));
        return
elseif (xxm(1)<=1 && xxm(2)>1) && (yym(1)<=1 && yym(2)>1)
        frac_area_r = (1-xxm(1))*(1-yym(1)) ...
                      /((xxm(2)-xxm(1))*(yym(2)-yym(1)));
        return
else
    frac_area_r = 0.0;
end

end % function frac_area

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
