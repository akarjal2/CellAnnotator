function ellipse_t=antti_fit_ellipse_to_polygon(x,y)

% see documentation from antti_fit_ellipse

x=round(x);
y=round(y);

if x(end)~=x(1) & y(end)~=y(1)
    x(end+1)=x(1);
    y(end+1)=y(1);
end

min_x=min(x);
min_y=min(y);

x=x-min_x+1;
y=y-min_y+1;

s_img=[max(y) max(x)];

img=zeros(s_img,'uint8');
img=antti_plot_line([x(1:(end-1)); x(2:end)],[y(1:(end-1)); y(2:end)],img,1);
img=imfill(img,'holes');

info=regionprops(img,'Centroid','Orientation','Majoraxis','minoraxis');
param=struct2array(info);

A=polyarea(x,y);

% param(1)=param(1)+min_x-1;
% param(2)=param(2)+min_y-1;
% param(3:4)=[sqrt(param(3)/param(4)*4/pi*A) sqrt(param(4)/param(3)*4/pi*A)];
% param(5)=-param(5);

ellipse_t.X0=param(1)+min_x-1;
ellipse_t.Y0=param(2)+min_y-1;
param(3:4)=[sqrt(param(3)/param(4)*4/pi*A) sqrt(param(4)/param(3)*4/pi*A)];
ellipse_t.a=param(3);
ellipse_t.b=param(4);
ellipse_t.phi=-param(5)*pi/180;


