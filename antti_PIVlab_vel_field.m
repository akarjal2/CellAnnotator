% input image1 and image2
% output x y u_filt v_filt

function [x,y,u_filtered1,v_filtered1]=antti_PIVlab_vel_field(image1,image2,filt_off)

if isempty(filt_off)
    filt_off=0;
end

% % % % Example script how to use PIVlab from the commandline
% % % % You can adjust the settings in "s" and "p", specify a mask and a region of interest
% % % clc; clear all
% % % 
% % % %% Create list of images inside specified directory
% % % directory=uigetdir; %directory containing the images you want to analyze
% % % suffix='*.tif'; %*.bmp or *.tif or *.jpg
% % % direc = dir([directory,filesep,suffix]); filenames={};
% % % [filenames{1:length(direc),1}] = deal(direc.name);
% % % filenames = sortrows(filenames); %sort all image files
% % % amount = length(filenames);

%% Standard PIV Settings
s = cell(10,2); % To make it more readable, let's create a "settings table"
%Parameter                       %Setting           %Options
s{1,1}= 'Int. area 1';           s{1,2}=64;         % window size of first pass
s{2,1}= 'Step size 1';           s{2,2}=32;         % step of first pass
s{3,1}= 'Subpix. finder';        s{3,2}=1;          % 1 = 3point Gauss, 2 = 2D Gauss
s{4,1}= 'Mask';                  s{4,2}=[];         % If needed, generate via: imagesc(image); [temp,Mask{1,1},Mask{1,2}]=roipoly;
s{5,1}= 'ROI';                   s{5,2}=[];         % Region of interest: [x,y,width,height] in pixels, may be left empty
s{6,1}= 'Nr. of passes';         s{6,2}=1;          % 1-4 nr. of passes
s{7,1}= 'Int. area 2';           s{7,2}=32;         % second pass window size
s{8,1}= 'Int. area 3';           s{8,2}=16;         % third pass window size
s{9,1}= 'Int. area 4';           s{9,2}=16;         % fourth pass window size
s{10,1}='Window deformation';    s{10,2}='*linear'; % '*spline' is more accurate, but slower

%% Standard image preprocessing settings
p = cell(8,1);
%Parameter                       %Setting           %Options
p{1,1}= 'ROI';                   p{1,2}=s{5,2};     % same as in PIV settings
p{2,1}= 'CLAHE';                 p{2,2}=1;          % 1 = enable CLAHE (contrast enhancement), 0 = disable
p{3,1}= 'CLAHE size';            p{3,2}=50;         % CLAHE window size
p{4,1}= 'Highpass';              p{4,2}=0;          % 1 = enable highpass, 0 = disable
p{5,1}= 'Highpass size';         p{5,2}=15;         % highpass size
p{6,1}= 'Clipping';              p{6,2}=0;          % 1 = enable clipping, 0 = disable
p{7,1}= 'Clipping thresh.';      p{7,2}=0;          % 0-255 clipping threshold
p{8,1}= 'Intensity Capping';     p{8,2}=0;          % 1 = enable intensity capping, 0 = disable

% % % %% PIV analysis loop
% % % if mod(amount,2) == 1 %Uneven number of images?
% % %     disp('Image folder should contain an even number of images.')
% % %     %remove last image from list
% % %     amount=amount-1;
% % %     filenames(size(filenames,1))=[];
% % % end
% % % x=cell(amount/2,1);
% % % y=x;
% % % u=x;
% % % v=x;
% % % typevector=x; %typevector will be 1 for regular vectors, 0 for masked areas
% % % counter=0;
%% PIV analysis loop:
% % % for i=1:2:amount
% % %     counter=counter+1;
% % %     image1=imread(fullfile(directory, filenames{i})); % read images    
% % %     image2=imread(fullfile(directory, filenames{i+1}));

max_val=max([image1(:);image2(:)]);
min_val=min([image1(:);image2(:)]);

image1=uint8((image1-min_val)*(255/double(max_val-min_val+1)));
image2=uint8((image2-min_val)*(255/double(max_val-min_val+1)));

image1 = PIVlab_preproc (image1,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2}); %preprocess images
image2 = PIVlab_preproc (image2,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2});
% % %     tic
% % %     [x{counter} y{counter} u{counter} v{counter} typevector{counter}] = piv_FFTmulti (image1,image2,s{1,2},s{2,2},s{3,2},s{4,2},s{5,2},s{6,2},s{7,2},s{8,2},s{9,2},s{10,2});
[x, y, u, v, ~] = antti_piv_FFTmulti (image1,image2,s{1,2},s{2,2},s{3,2},s{4,2},s{5,2},s{6,2},s{7,2},s{8,2},s{9,2},s{10,2});
% % %     toc
% % % %     clc
% % %     disp([int2str(i/amount*100) ' %']);
% % %     
% % %     % Graphical output (disable to improve speed)
% % %     %%{
% % %     figure(1)
% % %     imagesc(double(image1)+double(image2));colormap('gray');
% % %     hold on
% % %     quiver(x{counter},y{counter},u{counter},v{counter},'g','AutoScaleFactor', 1.5);
% % % %     quiver(x{counter},y{counter},u_filt{counter},v_filt{counter},'r','AutoScaleFactor', 1.5);
% % % %     quiver(x{counter},y{counter},u{counter}-u_filt{counter},v{counter}-v_filt{counter},'g','AutoScaleFactor', 1.5);
% % % %     quiver(x{counter},y{counter},u{counter}-mean_u,v{counter}-mean_v,'g','AutoScaleFactor', 1.5);
% % %     hold off;
% % %     axis image;
% % %     title(filenames{i},'interpreter','none')
% % %     set(gca,'xtick',[],'ytick',[])
% % %     drawnow;
% % %     %%}
% % % end

%% PIV postprocessing loop
% Settings
umin = -20; % minimum allowed u velocity
umax = 20; % maximum allowed u velocity
vmin = -20; % minimum allowed v velocity
vmax = 20; % maximum allowed v velocity
stdthresh=6; % threshold for standard deviation check
epsilon=0.15; % epsilon for normalized median test
thresh=3; % threshold for normalized median test

% % % u_filt=cell(amount/2,1);
% % % v_filt=u_filt;

% % for PIVresult=1:size(x,1)
u_filtered=u;
v_filtered=v;

if ~filt_off

    % % % typevector_filtered=typevector;
    % vellimit check
    u_filtered(u_filtered<umin)=NaN;
    u_filtered(u_filtered>umax)=NaN;
    v_filtered(v_filtered<vmin)=NaN;
    v_filtered(v_filtered>vmax)=NaN;
    % stddev check
    meanu=nanmean(nanmean(u_filtered));
    meanv=nanmean(nanmean(v_filtered));
    std2u=nanstd(reshape(u_filtered,size(u_filtered,1)*size(u_filtered,2),1));
    std2v=nanstd(reshape(v_filtered,size(v_filtered,1)*size(v_filtered,2),1));
    minvalu=meanu-stdthresh*std2u;
    maxvalu=meanu+stdthresh*std2u;
    minvalv=meanv-stdthresh*std2v;
    maxvalv=meanv+stdthresh*std2v;
    u_filtered(u_filtered<minvalu)=NaN;
    u_filtered(u_filtered>maxvalu)=NaN;
    v_filtered(v_filtered<minvalv)=NaN;
    v_filtered(v_filtered>maxvalv)=NaN;
    % normalized median check
    % Westerweel & Scarano (2005): Universal Outlier detection for PIV data
    [J,I]=size(u_filtered);
    % medianres=zeros(J,I);
    normfluct=zeros(J,I,2);
    b=1;
    for c=1:2
        if c==1; velcomp=u_filtered;else;velcomp=v_filtered;end %#ok<*NOSEM>
        for i=1+b:I-b
            for j=1+b:J-b
                neigh=velcomp(j-b:j+b,i-b:i+b);
                neighcol=neigh(:);
                neighcol2=[neighcol(1:(2*b+1)*b+b);neighcol((2*b+1)*b+b+2:end)];
                med=median(neighcol2);
                fluct=velcomp(j,i)-med;
                res=neighcol2-med;
                medianres=median(abs(res));
                normfluct(j,i,c)=abs(fluct/(medianres+epsilon));
            end
        end
    end
    info1=(sqrt(normfluct(:,:,1).^2+normfluct(:,:,2).^2)>thresh);
    u_filtered(info1==1)=NaN;
    v_filtered(info1==1)=NaN;

    % % % typevector_filtered(isnan(u_filtered))=2;
    % % % typevector_filtered(isnan(v_filtered))=2;
    % % % typevector_filtered(typevector==0)=0; %restores typevector for mask
end
    % Interpolate missing data
u_filtered=inpaint_nans(u_filtered,4);
v_filtered=inpaint_nans(v_filtered,4);

    % x=repmat([x(1,1)-(x(1,2)-x(1,1)) x(1,:) x(1,end)+(x(1,2)-x(1,1))],[size(x,1)+2 1]);
    % y=repmat([y(1,1)-(y(2,1)-y(1,1)) ; y(:,1) ; y(end,1)+(y(2,1)-y(1,1))],[1 size(y,2)+2]);


x=repmat([0 x(1,:) size(image1,2)+1],[size(x,1)+2 1]);
y=repmat([0 ; y(:,1) ; size(image1,1)],[1 size(y,2)+2]);

u_filtered1=zeros(size(u_filtered)+[2 2]);
u_filtered1(2:end-1,2:end-1)=u_filtered;
u_filtered1([1 end],2:end-1)=u_filtered1([2 end-1],2:end-1);
u_filtered1(1:end,[1 end])=u_filtered1(1:end,[2 end-1]);

v_filtered1=zeros(size(v_filtered)+[2 2]);
v_filtered1(2:end-1,2:end-1)=v_filtered;
v_filtered1([1 end],2:end-1)=v_filtered1([2 end-1],2:end-1);
v_filtered1(1:end,[1 end])=v_filtered1(1:end,[2 end-1]);


% % end

% % % i=1;
% % % counter=1;
% % % figure(2)
% % % imagesc(double(image1)+double(image2));colormap('gray');
% % % hold on
% % % % quiver(x{counter},y{counter},u{counter},v{counter},'g','AutoScaleFactor', 1.5);
% % % quiver(x{counter},y{counter},u_filt{counter},v_filt{counter},'r','AutoScaleFactor', 1.5);
% % % %     quiver(x{counter},y{counter},u{counter}-u_filt{counter},v{counter}-v_filt{counter},'g','AutoScaleFactor', 1.5);
% % % %     quiver(x{counter},y{counter},u{counter}-mean_u,v{counter}-mean_v,'g','AutoScaleFactor', 1.5);
% % % hold off;
% % % axis image;
% % % title(filenames{i},'interpreter','none')
% % % set(gca,'xtick',[],'ytick',[])
% % % drawnow;
% % % 
% % % clearvars -except p s x y u v typevector directory filenames u_filt v_filt typevector_filt