% merge red, green, blue and gray channels
% only one of the inputs have to be non-empty
% all four inputs need to be supplied

% example: comb = antti_merge_RGBG(Red,[],[],Grey);

function comb=antti_merge_RGBG(red,green,blue,grey)

% handle grey
if ~isempty(grey)
    comb1=grey;
    comb2=grey;
    comb3=grey;
else
    comb1=zeros(max([size(red);size(green);size(blue)]));
    comb2=comb1;
    comb3=comb1;    
end

% handle red
if ~isempty(red)
    comb1=comb1+red;
end

% handle green
if ~isempty(green)
    comb2=comb2+green;
end

% handle green
if ~isempty(blue)
    comb3=comb3+blue;
end

comb=zeros([size(comb1) 3]);
comb(:,:,1)=comb1;
comb(:,:,2)=comb2;
comb(:,:,3)=comb3;

% if all input data types are uint8, use uint8 for output as well 
if (isa(red,'uint8') | isempty(red)) & (isa(green,'uint8') | isempty(green)) & (isa(blue,'uint8') | isempty(blue)) & (isa(grey,'uint8') | isempty(grey))
    comb=uint8(comb);
end







