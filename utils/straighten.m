
function [IM2] = straighten(IM,pts,width,showFigs)
%STRAIGHTEN Straightens an image based on the spline interpolation between
%the input points.
% based on the straightenLine function in the Straightener plugin from ImageJ1
% https://github.com/imagej/imagej1/blob/master/ij/plugin/Straightener.java#L101
% https://github.com/imagej/imagej1/blob/master/ij/gui/PolygonRoi.java
% https://ch.mathworks.com/matlabcentral/fileexchange/72266-straighten
%
% [IM2] = straighten3D(IM,pts,width)
%
% INPUT
% IM       input image
% pts      points
% width    extract width
% showFigs figures for debugging
%
% OUTPUT
% IM2      straightened image

if ndims(IM)==3
    for ct = 1:size(IM,3)
        % process each slice independently
        disp(num2str(ct))
        IM2(:,:,ct)=straighten(IM(:,:,ct),pts,width);
    end
    return;
end
% must make spline with one pixel segment lengths
[x,y]=fitSplineForStraightening(pts(:,1),pts(:,2));
if showFigs
    figure(99);clf;
    imagesc(IM);axis image
    hold on;
    plot(x,y,'o');
    plot(pts(:,1),pts(:,2),'r.-');
    hold off
    title('spline')
    colormap('gray')
    axis off
end
for ct = 1:length(x)
    if ct == 1
        dx = x(2)-x(1); %extrapolate first direction
        dy = y(1)-y(2);
    else
        dx = x(ct)-x(ct-1);
        dy = y(ct-1)-y(ct);
    end
    %want to move 1 pixel when we move the pointer dx dy
    L = sqrt(dx^2+dy^2);
    dx = dx/L;
    dy = dy/L;
    
    xStart = x(ct)-(dy*width)/2; %start half the width away
    yStart = y(ct)-(dx*width)/2;
    if showFigs && mod(ct,100)==1
        % show every 100th
        figure(98);
        if ct==1
            imagesc(IM)%;axis image;
            hold on
            plot(x,y,'-');
            colormap('gray')
            axis equal, axis tight, axis off
        end
        plot(x(ct),y(ct),'o-');
        plot([xStart,xStart+dy*width],[yStart,yStart+dx*width],'r--')
        pause(0.1)
    end
    for ct2 = 1:width
        try
            IM2(ct,ct2) = getInterpolatedValue(IM,xStart,yStart);
        catch
            if showFigs
                figure(888)
                subplot(1,2,1)
                imagesc(IM2)
                hold on
                plot(ct2,ct,'rx')
                hold off
                title(['IM2(' num2str(ct) ',' num2str(ct2) ')'])
                axis equal tight
                subplot(1,2,2)
                imagesc(IM)
                hold on
                plot(xStart,yStart,'ro')
                hold off
                title(['IM(' num2str(xStart) ',' num2str(yStart) ')'])
                axis equal tight
            end
            keyboard
        end
        xStart=xStart+dy;
        yStart=yStart+dx;
    end
end
end
function I = getInterpolatedValue(IM,y,x)
% interpolation of pixelvalues one-indexed
% IM = [0,1,2;2,3,4;4,5,6];
% getInterpolatedValue(IM,1,1);
% > 0
x = [floor(x),ceil(x),rem(x,1)];
y = [floor(y),ceil(y),rem(y,1)];
bgVal=0; % background value if outside
try
    I(1,1) = IM(x(1),y(1))*(1-x(3))*(1-y(3));
catch
    I(1,1) = bgVal;
end
try
    I(2,1) = IM(x(1),y(2))*(1-x(3))*(0+y(3));
catch
    I(2,1) = bgVal;
end
try
    I(1,2) = IM(x(2),y(1))*(0+x(3))*(1-y(3));
catch
    I(1,2) = bgVal;
end
try
    I(2,2) = IM(x(2),y(2))*(0+x(3))*(0+y(3));
catch
    I(2,2) = bgVal;
end
I = sum(I(:));
end
function [xspl,yspl] = fitSplineForStraightening(xVal,yVal)
%return a spline fit with distances of 1 between the points
%based on https://github.com/imagej/imagej1/blob/master/ij/gui/PolygonRoi.java#L1006
%generate intermediate spline with approximately half pixel steps
% check if unique values in y, in x already done
% [yVal,ia]=unique(yValIn);
% if length(yVal)==length(yValIn)
%     xVal=xValIn;
% else
%     xVal=xValIn(ia);
% end
pp = spline(xVal,yVal);
%pp = spline(xVal,yVal);
xIspl = interpolate(xVal,yVal,0.5);
yIspl = ppval(pp,xIspl);
L=0;%measure spline distance
%generate spline
xspl(1) = xVal(1);
yspl(1) = yVal(1);
L=0; %keep track of length of spline
ptw = 1; %points written
for ct = 2:length(xIspl)
    dx=xIspl(ct)-xIspl(ct-1);
    dy=yIspl(ct)-yIspl(ct-1);
    d = sqrt(dx^2+dy^2);
    L=L+d;
    overshoot = L-ptw;
    if overshoot>0 %we went over the length, must add a point on the spline
        ptw=ptw+1;
        frac = overshoot/d; %fractional overshoot for the last step
        %move back in a straight line towards the previous point
        xspl(ptw) = xIspl(ct) - frac*dx;
        yspl(ptw) = yIspl(ct) - frac*dy;
    end
end
end
function [Xi] = interpolate(xVal,yVal,d)
%interpolate line with aproximate distance d between points
Xi = [];
dists = sqrt(diff(xVal).^2+diff(yVal).^2); %segment distances 
for ct = 1:length(xVal)-1
    x = linspace(xVal(ct),xVal(ct+1),round(dists(ct)/d)+1);x(end)=[];
    Xi = [Xi,x];
end
end
