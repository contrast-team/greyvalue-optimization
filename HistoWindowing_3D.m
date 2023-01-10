
clear all
close all
%%
% author initial version : Arne Maes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modif by Piens Catherine's father as part of her Thesis work
%   - allow for asymetrical boundaries
%   - simplification of question structure
%   - automatic selection of the representative image
%   - introduction of preset value in specific session
%   - loop until values LL/UL are OK
%   - introduction of "mask concept" to visualize boundaries impact on
%     image
%   - introduction of an enhanced contrast picture to optimize constrasts
%     and improve visualisation (modified Histogram equalizer)
%   - for BMP/percentage creates specific directory per percentage
%
% 1Q21
%   - adding 2.5D and 3D histo windowing
%      2.5 Histo windowng based on frequence of colors in all the slices
%  
%      3D - Same as histo 2.5 but allows to create view on the 3 axes
%         - Files with images summaries are produced in directory-  just before creating output 
%          - selection UL and LL  made in 3 steps (value can be changed at each step
%           1. via histogram
%           2. histogram based value are presented in %
%           3. % values are translated in LL UL and can still be changed
%         -retained value indicate the % of pixels within UL and LL (an
%         indirect indication of the volume of the selected object)
%   - Flags
%   - Canny Edge drawings (borders) - the threshold value can be
%   automatically selected or set as an initial value and changed during
%   processing  (if CANNY flag =1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FLAGS Code tailoring variables  
TEST=0;   %set to 1 to reduce size for tests
TIM=0;    %set to O for simplified sequence
CANNY=1;  %indicate Canny Edge is requested (O to avoid to run Canny optimizer
ZIP=0;    %1 enable the auto creation of a ZIP file for 2D BMP (to easy data transfer)


%  PRESET DEFAULT VALUES

% Default data directory
%shortp='/Volumes/mpfast/Users/dossier cath/CPI_3_WF_CA4+2_63/CPI_3_WF_CA4+2_63_HBM_4d_/recon/rec/*.tif'; %
%shortp='/Volumes/mpfast/Users/dossier cath/*.tif';
shortp='*.tif';


% Default values (corresponding)
Voxelsize= 6 * power(10,-6)
Voxelvol=  power(Voxelsize,3)
Thre=0.2;    % Starting threshold value for Canny Edge  (0 will take the Canny default)
oThre=0;
a=0.75;      %lower limit default percent (remove a percent of point of darkest value          
au=0.001;    %higher limit default percent (remove au percent of points with lightest values
buffer=0;    %manual correction over the lower border
bufferu=0;   %manual correction over the upper border
section=50;  %representative image percent (1-99) (the higher, the higher the number of the represnetative image  
UL=0    ;   % preset UL (for manual method - 65535 is the default, will be reset to hist size if left to 0
LL=1;       % preset LL for manual method
metho='Sobel';    % Sobel or Canny

%get specs of the display to adapt size of the cuts
screen_size = get(0,'ScreenSize');
pc_width  = screen_size(3);
pc_height = screen_size(4);

% Step 1 SELECT DATASAT (take any image)
temp=msgbox('Select directory of the *.TIF and click on any file');
[repr_image, repr_folder] = uigetfile(shortp, 'Select an image of the dataset (.tif).');
if isequal(repr_image,0)             %if no selection
   disp('User selected Cancel');
   return
else
    fileN=repr_image(1:7);   
    file_names = dir(fullfile(repr_folder, '*.tif')); %reads all names of Tif
    number_images = numel(file_names);                %number of images
end
delete(temp);            % delete warning msgbox

% Step 2 SELECT REFERENCE IMAGE
asw=inputdlg('Enter slice ref if you want a specific one, 0 will select the default one, a value 0<x1 woukld be handled as a %');
asw=str2double(asw);
if  isnan(asw)           % if no input go to default (based on section)
    asw=0;
end
if asw>number_images     % if image doesn't exist
    asw=0;
end
if asw>1                 % if valid select the requested image
    mid_img=asw;
    %I1_repr = imread(fullfile(repr_folder,repr_image)); % use selected image
else
    if asw==0            % take the preselected section
        %section=section;
    else
        section=asw*100;% take the number between 0 and 1 AS percentage section
    end
   mid_img = uint16(number_images*section/100);  %takes image at selected section
end
sprintf('Selected image %4.0f ',mid_img)         % written to the console for info
   I1_repr = imread(fullfile(repr_folder,file_names(mid_img).name));
[m,n] = size(I1_repr);
hist3DT=0;

%Step 3 SELECT DESIRED Processing type
answer0=questdlg('Select type of work','2D 2.5D 3D selection','2D','2.5D','3D','2D');
 if(strcmp(answer0,'2.5D'))                                % if 2.5D read all slices (but just keep histo data 
  wait = waitbar(0,'Images are being read...', ...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
  setappdata(wait,'canceling',0);

   for ni=1:number_images        %read  dataset image by image
        I1 = imread(fullfile(repr_folder,file_names(ni).name));
        hist=imhist(I1,65535);
        hist3DT=hist3DT+hist;    %sum the histo
    waitbar(ni/number_images)    
   end
delete(wait);
     hist=hist3DT;                    % store the histo in...hist (to be compatible with 2D 
     answer0='2D';                    % all the rest is like for 2D
 else                                 % for 2 D (and 3D) do an histo on selected slide
     hist = imhist(I1_repr, 65535);
 end

  if (UL==0)    % if no UL preset then set to histo size ( 65535)
    UL = size(hist,1);
  end
%Chosen_method = 'canceled';
if TIM || ~(strcmp(answer0,'3D'))  % if TIM=O and 3D by pass this section (and go to 3D directly
    
%STEP 4 'START CALIBRATION PHASE (selection of LL UL) for 2D and 2.5D
% select LL UL input method (Percentage, Histogram, Manual)
rep = 'Yes';
while (strcmp(rep , 'Yes'))  % Loop until user is OK with UL/LL selection
answer1 = questdlg('Select calibration method', ...
	'Calibration', ...
	'Percentage','On histogram' ,'Enter numerically','Percentage');
        if isempty(answer1)
         disp('User has cancelled.')
         return
        end

 
switch answer1  %Select LL UL based on selected method   
    case 'Percentage'            % methode based on % of image non considered
        Chosen_method = 'Percentage';
        %Ask user for threshold value (percentage) LL
        prompt = {'Enter LL threshold value (in percentage)'};
        dlgtitle = 'Enter Lower threshold value ';
        a=a*100;
        definput = {num2str(a)};
        dims = [1 65];
        opts.Resize= 'on';
        answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
        a = str2double(answer{1})/100;
        aper=answer{1};
        
        %Ask user for threshold value (percentage) UL
        prompt = {'Enter UL threshold value (in percentage)'};
        dlgtitle = 'Enter Upper threshold value ';
        au=au*100;
        definput = {num2str(au)};
        dims = [1 65];
        opts.Resize= 'on';
        answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
        au = str2double(answer{1})/100;
        
        %Ask user if you want to add a buffer
        prompt = {'Enter value for lower buffer (in grey value 0 - 65535)'};
        dlgtitle = 'Enter Lower buffer value ';
        definput = {num2str(buffer)};
        dims = [1 65];
        opts.Resize= 'on';
        answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
        buffer = str2double(answer{1});
        %Ask user if you want to add a buffer
        prompt = {'Potentially, enter value for upper buffer (in grey value 0 - 65535)'};
        dlgtitle = 'Enter higher buffer value';
        definput = {num2str(bufferu)};
        dims = [1 65];
        opts.Resize= 'on';
        answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
        bufferu = str2double(answer{1});
        
        %Calculate LL based on percentages
        hist_sum = sum(hist);
        hist_size = size(hist,1);
        L_sum = 0;
        U_sum = 0;
        for n = 1:hist_size
            L_sum = L_sum + hist(n);
            if L_sum > a*hist_sum
                LL = n - buffer;
                break
            end
        end
        % Calculate UL
        for n = 1:(hist_size-1)
            i = hist_size-n;
            U_sum = U_sum + hist(i);
            if U_sum > au*hist_sum
                UL = hist_size-n+bufferu;
                break
            end
        end
      
       case 'On histogram'   % note ensure 1st click is lower than 2nd
                sprintf('current value:  LL= %d  UL = %d',LL,UL)
                Chosen_method = 'Enter manually on histogram';
                f0=figure('name','selection histo')
                %imhist(I1_repr, 65535)
                plot(hist)
                xline(LL,'-.b','LL'); xline(UL,'-.b','UL');
                title('First select lower limit, then select upper limit')
                [x,y] = ginput(2);       %collect coordinates of selection
                LL = x(1);               %select x value (ie shade of grey value)
                UL = x(2);
                close (gcf)
                delete(f0);
                
        case 'Enter numerically'      %get directly the values to be used
                sprintf('current value:  LL= %d  UL = %d',LL,UL)
                prompt = {'Enter value for lower limit (in grey value 0 - 65535)'};
                dlgtitle = 'Lower limit value';
                definput = {num2str(LL)};
                dims = [1 65];
                opts.Resize= 'on';
                answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
                LL = str2num(answer{1});
        
                prompt = {'Enter value for upper limit (in grey value 0 - 65535)'};
                dlgtitle = 'Upperlimit value';
                definput = {num2str(UL)};
                dims = [1 65];
                opts.Resize= 'on';
                answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
                UL = str2num(answer{1});        
               
end             
%  Calculations example image and histograms to let user decide to continue or not with selected values)
    sprintf('final values:  LL= %d  UL = %d',LL,UL)
    width = UL-LL;
    
% Step 5 IMAGE TRANSFORMATIONS
%Calculation of images
    I2_16= winhist(I1_repr,LL,UL,0);    %histo windowing function call 
    I2_ce= edge(I1_repr,metho)*250;       % Sobel filter
    I2_m=mmask(I2_16);                  % call mask function
    I2_ci=mhisteq(I2_16);               % call histo equalizer 
    
    [I2ce,Thres]=edge(I2_16,'Canny'); %getting Canny edge default Threshold
    %I2ce=I2ce*128;% (Calculated in case we are not refined canny
    if oThre>0
        Thre=oThre
    end
    if Thre==0
    Thre=Thres(2);
    end %Taking high threshold (as starting point)                      %Taking high threshold (as starting point) 
    
    
   if CANNY                             %will only run if CANNY flag is on 1 to optimize Thre(shold) value
    while(Thre>0)
       f9=figure('name','selection de canny threshold');
       I2ce = (edge(I2_16,'Canny',Thre))*128;
       %I2ce=I2ce*126; 
       oThre=Thre;
       imshow(I2ce);                    %diplay the canny pitcure with optimized value
       Thre=inputdlg(sprintf('Enter Threshold vlaue for Canny Edge 0 to confirm current one  ( %3.3f)',Thre));
       Thre=str2double(Thre)
       delete(f9);
    end
   end

% STEP 6 display selected image transformed  based of current LL/UL selection
f1 = figure('name','Initial image, histo windowed,mask, Histo equalized, Canny edge, Sobel');
  montage({I1_repr,I2_16,I2_m,I2_ci,I2_ce,I2ce});
    set_fig_position(f1,pc_height*0.3, pc_width*0.0, pc_height*0.5, pc_width*0.6);
    title(['Transformation of selected image nr ',num2str(mid_img),'/',num2str(number_images)])

f2 = figure('name','Histogram comparition');                      %display 2 histo
  subplot(2,1,1)             %histo non modified image
    %imhist(I1_repr,65535);
    plot(hist);
    xline(LL,'-.b','LL'); xline(UL,'-.b','UL');
    title(['Histogram Original image',num2str(mid_img),'/',num2str(number_images)])
    axis 'auto y'
  subplot(2,1,2)              %histo modified image
    imhist(I2_16,65535);
    title(['Histo windowed image','UL = ',num2str(uint16(UL)) , 'LL= ',num2str(uint16(LL))])
    set_fig_position(f2,pc_height*0.3, pc_width*0.65, pc_height*0.65, pc_width*0.35);

%ask user if he is OK with the selection or redo the while
rep = questdlg(sprintf('Current Values  LL= %5.0f  UL = %5.0f \n Do you want to change values ? ',LL,UL),'Possibility to change Upper and Lower Limits'); % is USER OK with current selection?
if isempty(rep)
         disp('User has cancelled1.')
         return
end
        delete(f1);
        delete(f2);
end
end
if(strcmp(answer0,'2D'))            % process the specific part of 2D (and 2.5D)
    
% STEP 7 PROCESSING OF THE IMAGES BASED ON UL/LL Selection in calibration
% user choice on export formats (BMP, JPG or TIF)
BMP = 0;
JPG = 0;
TIF = 0;
bit_depth = 255;          %BMP AND JPG
answer = questdlg('Export images as BMP or JPG?', ...
	'BMP, JPG or TIF?', ...
	'BMP','JPG','TIF','BMP');

if isempty(answer)
    disp('User has cancelled.')
    return
end
  if strcmp ('Percentage', answer1)
     fext = append(answer,aper);
  else
     fext=answer;
  end
mkdir(repr_folder, fext);
new_folder=fullfile(repr_folder,fext); %pathway to new folder
% Handle format as requested
    switch answer
        case 'BMP'
            BMP = 1;
        case 'JPG'
            JPG = 1;
        case 'TIF'
            TIF = 1;
            bit_depth = 65535;
         case 'Cancel'
             message = msgbox('Operation has been cancelled by the user', 'Cancelled');
             return
    end
close all
width = UL-LL;
wait = waitbar(0,'Images are being processed...', ...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(wait,'canceling',0);

%process all images of the directory
tic
nim=number_images;
for n=1:number_images        %Process dataset image by image
    nt=n;
    if getappdata(wait,'canceling')
        message = msgbox('Operation has been cancelled by the user', 'Cancelled');
        return
    end
    I1 = imread(fullfile(repr_folder,file_names(nt).name));
    I1_d = double(I1);
    I2_d = (I1_d-LL)/width*bit_depth;
    
    baseFilename= file_names(n).name;  %Naming of the image
    fullFilename=fullfile(repr_folder,baseFilename);
    [folder, baseFilename, extension] = fileparts(fullFilename);
    
    %Write modified images in right format
    if BMP
        I2 = uint8 (I2_d);
        imwrite(I2, fullfile(new_folder,strcat(baseFilename,'.bmp'))); %Writing of the image
    elseif JPG
        I2 = uint8 (I2_d);
        imwrite(I2, fullfile(new_folder,strcat(baseFilename,'.jpg'))); %Writing of the image
    elseif TIF
        I2 = uint16 (I2_d);
        imwrite(I2, fullfile(new_folder,strcat(baseFilename,'.tif'))); %Writing of the image
    end
    waitbar(n/number_images,wait)
end
delete(wait)

toc
if ZIP
%create a zip file to ease copies (
 if (BMP)
  cd BMP
  zip('tiffzip',{'*.bmp'});
  cd ../
 end
end

% Make figures to compare starting image and resulting image to prepare
% file creations
f3 = figure(3);
   montage({I1_repr,I2_16,I2_m,I2_ci,I2_ce,I2ce});
   set_fig_position(f3,0, 0, pc_height*.7, pc_width*0.65);

f4 = figure(4);
  subplot(2,1,1)
    plot(hist);
    xline(LL,'-.b','LL'); xline(UL,'-.b','UL');
    axis 'auto y'
    subplot(2,1,2)
    imhist(I2_16,65535);
    %axis 'auto y'
    set_fig_position(f4,0, pc_width*0.65, pc_height*.7, pc_width*0.35);



%%
% STEP 8 PRODUCTION OF LOG FILES

%create text file
file_name = 'Histogram_windowing.txt';
%create log folder
mkdir(new_folder, 'LOG files');
log_folder=fullfile(new_folder, 'LOG files');
out = fullfile(log_folder,file_name);

%Write parameters to log file
fileID = fopen(out,'w+');
fprintf(fileID, ' ===================================================================\r\n');
fprintf(fileID, '||                         LOG INFORMATION                        ||\r\n');
fprintf(fileID, ' ===================================================================\r\n\r\n\r\n');
fprintf(fileID, 'Dataset folder: \t%s\r\n', repr_folder);
fprintf(fileID, 'Representative image of dataset: \t%s\r\n', baseFilename); % "sObject" is a string.
fprintf(fileID, 'Number of images: \t%i\r\n', number_images);
fprintf(fileID, ' ===================================================================\r\n\r\n');
fprintf(fileID, 'Method of choice: \t%s\r\n', Chosen_method); 
if (strcmp(Chosen_method,'Percentage'))
fprintf(fileID, 'Lower Threshold percentage (in percents): \t%f\r\n', a*100); 
fprintf(fileID, 'Upper Threshold percentage (in percents): \t%f\r\n', au*100);
end
fprintf(fileID, 'lower Buffer: \t%f\r\n\r\n', buffer);
fprintf(fileID, 'Upper Buffer: \t%f\r\n\r\n', bufferu);
fprintf(fileID, 'Lower limit of grey value: \t%f\r\n', LL); 
fprintf(fileID, 'Upper limit of grey value: \t%f\r\n', UL); 
fprintf(fileID, ' ===================================================================\r\n\r\n');

%fprintf(fileID, 'Value\t%f\r\n', value); % "value" is a float.
fclose(fileID); % Close file.
%save figures on disk files
saveas(f3, fullfile(log_folder,'comparison_images.png')); %Writing of the image
savefig(f3, fullfile(log_folder,'comparison_images.fig')); %Writing of the image
saveas(f4, fullfile(log_folder,'comparison_histograms.png')); %Writing of the image 
savefig(f4, fullfile(log_folder,'comparison_histograms.fig')); %Writing of the image   


message = msgbox('Operation Completed!', 'Success');
return
%%%%%%%%%%%%%%%%%  END OF 2D AND 2.5D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Start of 3D specific - requires a lot of memory but allows 3D outputs
%(could require some setting changes in Mathlab to allow larger memory
%allocations
%%%%%%%%%%%%%%%%%
else %starting 3D
    
temp=msgbox('Wait for initalization of the 3D matrix(check you removed the memory limit if needed'); %msg written to console
%if 1==1  % dummy if ...   

niz=number_images;       %Zdimention
[m,n] = size(I1_repr);
nix=m;                   %Xdimention
niy=n;                   %Ydimention
six=0;                   %
siy=0;
siz=0;
if TEST                  %allow to work on a subset of 200 images for test on smaller machine
    siz=uint16(niz/2);
    niz=200;
end

selx=uint16(nix/2);      %default representative images (middle of the selection
sely=uint16(niy/2);
selz=uint16(niz/2);

IT=uint16(zeros(nix,niy,niz));  %creation of the 3D matrix in 1 step (performance)
delete(temp);
wait = waitbar(0,'Please wait. Images are being loaded in memory...', ...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(wait,'canceling',0);

%bit_depth=255;
bit_depth=65535;                %Tif 16bits (expectation to work on TIF             %
histo3=1;
rx=0;
ry=0;
rz=0;
repi=1;
for n=1:niz        %Process dataset image by image
    nt=n+siz;
    I1 = imread(fullfile(repr_folder,file_names(nt).name)); % read images and store in IT
    IT(:,:,n) = uint16(I1);%load with raw data    
    waitbar(n/niz)
end
delete(wait);
repi =1;

sprintf('initial values: LL= %d  UL = %d',LL,UL)  %write to console
                                % display histo and make a first selection
                                % on histo
rep = 'Yes';
while (strcmp(rep , 'Yes'))  % Loop until user is OK with UL/LL selection (based on Percerntages now)
 
                                Chosen_method = 'Enter manually on histogram';
   f0=figure('name','Initial selection on histogram'); 
        imhist(IT, bit_depth)
        xline(LL,'-.b','LL'); xline(UL,'-.b','UL');
        title('First select lower limit, then select upper limit')
        [x,y] = ginput(2);
        LL = x(1);
        UL = x(2);
        close (gcf);
        delete(f0);
               
  sprintf('new value after histo:  LL= %d UL = %d',LL,UL)
     hist3D = imhist(IT, bit_depth);     % compute 3D histo (on 3D matrix)
     LL=uint16(LL);
     UL=uint16(UL);
     T3D=sum(hist3D);                %Total number of pixels in 3D
     R3D=sum(hist3D(LL : UL));       %Pixels inside windowing
     Ra3D=R3D/T3D;                    
     L3D=sum(hist3D(1 : (LL-1)))/T3D;      % percent pixels below LL
     U3D=sum(hist3D(UL+1 : bit_depth))/T3D;% percent pixels above UL
  
  sprintf('Retained pixel share = %d  Under LL = %d above UL %d',Ra3D,L3D,U3D)
     a=L3D;     % to reuse initial code
     au=U3D;
     hist=hist3D;
   

        Chosen_method = 'Percentage';
        prompt = {'Enter LL threshold value (in percentage)'};
        dlgtitle = 'Enter Lower threshold value ';
        a=a*100;
        definput = {num2str(a)};
        dims = [1 65];
        opts.Resize= 'on';
        answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
        a = str2double(answer{1})/100;
        aper=answer{1};
        
        %Ask user for threshold value (percentage) UL
        prompt = {'Enter UL threshold value (in percentage)'};
        dlgtitle = 'Enter Upper threshold value ';
        au=au*100;
        definput = {num2str(au)};
        dims = [1 65];
        opts.Resize= 'on';
        answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
        au = str2double(answer{1})/100;
        
        buffer=0;                      %buffers no longer used but for compatibility
        bufferu=0;
        
        %Calculate LL
        hist_sum = sum(hist);
        hist_size = size(hist,1);
        L_sum = 0;
        U_sum = 0;
        for n = 1:hist_size
            L_sum = L_sum + hist(n);
            if L_sum > a*hist_sum
                LL = n - buffer;
                break
            end
        end
        
        % Calculate UL
        for n = 1:(hist_size-1)
            i = hist_size-n;
            U_sum = U_sum + hist(i);
            if U_sum > au*hist_sum
                UL = hist_size-n+bufferu;
                break
            end
        end

      
        
       
    
   LLx=sprintf('LL Value (1-65535) %5.0f',LL);    %request if a new XYZ view is wanted
   ULx=sprintf('UL Value (1-65535) %5.0f',UL);    %  all 3 at 0 (or nothing) means view is OK
   definput={num2str(LL),num2str(UL)};
   prompt={LLx,ULx};
   dlgtitle=sprintf('Current Values  LL= %5.0f  UL = %5.0f \n Do you want to change values ? ',LL,UL);
   dims=[1 50; 1 50];
   rt=inputdlg(prompt,dlgtitle,dims,definput);
    if isempty(rt)
         disp('User has cancelled.')
         return
    end
    LL = str2num(rt{1});
    UL = str2num(rt{2});   
     
    
     R3D=sum(hist3D(uint16(LL) : uint16(UL)));       %Pixels inside windowing
     Ra3D=R3D/T3D;                    
     L3D=sum(hist3D(1 : (LL-1)))/T3D;      % percent pixels below LL
     U3D=sum(hist3D(UL+1 : bit_depth))/T3D;% percent pixels above UL

 sprintf('value percentage based:  LL= %d  UL = %d',LL,UL)
 sprintf('Retained pixel share   = %d  Under LL = %d above UL %d',Ra3D,L3D,U3D)
 sprintf('Retained object volume = %d', R3D*Voxelvol')
 
while (repi < 3)                           % while another XYZ view is requested   
    
I3=uint16(reshape(IT(selx,:,:),niy,niz));  %raw image in X dimention (selx position
I33=winhist(I3,LL,UL,bit_depth);               %Histo windowed
I4=reshape(IT(:,sely,:),nix,niz);          %reshape to have a 2 dimention image
I44=winhist(I4,LL,UL,bit_depth);
I5=uint16(reshape(IT(:,:,selz),nix,niy));  %Z dimention
I55=winhist(I5,LL,UL,bit_depth);
I333=mhisteq(I33);                          %histo equalizer functions
I444=mhisteq(I44);
I555=mhisteq(I55);

f7 = figure('name','Selected image Window equalized and Histo equalized views, XYZ views');
  montage({I33,I44,I55,I333,I444,I555});       %display  initial ,mask of changes and changed picture
   title(['Transformation of selected images nr ',num2str(selx),'/',num2str(sely),'/',num2str(selz)]);
   set_fig_position(f7,0, 0, pc_height*.7, pc_width*0.65);

f6 = figure('name','Histogram');                              %display 2 histo -compares 2D based histo
  subplot(2,1,1)
    imhist(I1_repr,bit_depth);
    title("Histogram Selected 2D")
    axis 'auto y';
  subplot(2,1,2)
    plot(hist3D);   
    title("test") 
    xline(LL,'-.b','LL'); xline(UL,'-.b','UL');
    %title("Histogram 3D")
    axis 'auto y';
    title("Histogram 3D")
    set_fig_position(f6,0, pc_width*0.65, pc_height*.7, pc_width*0.35);
    repi = 0;                                    %leaves the while except if new images are requested


  if CANNY                             %will only run if CANNY flag is on 1
    if Thre==0                         % if no predefined Canny value (would depend on image type)
       [I33ce,Thres]=edge(I3,'Canny'); %    getting Canny edge default (high)Threshold   
       Thre=Thres(2)
    end                      
    
    while(Thre>0)                      % while not the right Canny threshold
       f9=figure('name','selection de canny threshold');
       I33ce = (edge(I33,'Canny',Thre))*128;
       oThre=Thre;
       imshow(I33ce);
       title (['Threshold value',num2str(oThre)]);
       Thre=inputdlg(sprintf('Enter Threshold value for Canny Edge 0 to confirm current one  ( %3.3f)',Thre));
       Thre=str2double(Thre);
       delete(f9);
    end
    Thre=oThre;
    I33ce= ((edge(I33,'Canny',oThre))*64);      % Calculate 3 view with Canny 
    I44ce= ((edge(I44,'Canny',oThre))*64);
    I55ce= ((edge(I55,'Canny',oThre))*64);               
    
    f8 = figure('name','Mask and canny views');% diplay Canny and Mask views
      montage({mmask(I33),mmask(I44),mmask(I55),I33ce,I44ce,I55ce});  %displays mask views and Sobel/Canny
        title(['Mask and Canny Edge views for  ',num2str(selx),'/',num2str(sely),'/',num2str(selz),'   CannyThreshold = ',num2str(oThre)]);
        set_fig_position(f8,0, pc_width*0.5, pc_height*.8, pc_width*0.5);
  end

   ttx=sprintf('X Slice (1- %d)',nix);    %request if a new XYZ view is wanted
   tty=sprintf('Y Slice (1- %d)',niy);    %  all 3 at 0 (or nothing) means view is OK
   ttz=sprintf('Z Slice (1- %d)',niz);    %  il only 1 change the 2 other dimentions remains the same value
   prompt={ttx,tty,ttz};
   dlgtitle='Select views to be displayed (all 0 to leave)';
   dims=[1 50; 1 50;1 50];
   rt=inputdlg(prompt,dlgtitle,dims);
    if isempty(rt)
         disp('User has cancelled.')
         return
        end
   rx = str2num(rt{1});
   ry = str2num(rt{2});
   rz = str2num(rt{3});
  
   if  isempty(rx)||(rx==0)           % if no  reponse or 0
    rx=selx;
    repi=repi+1;
   end
   if  isempty(ry)||(ry==0)          
    ry=sely;
    repi=repi+1;
   end
   if  isempty(rz)||(rz==0)            
    rz=selz;
    repi=repi+1;              %  if none of the 3 values introduced repi=3 (and while continues)
  end  
   selx=rx;                   % new selected views
   sely=ry;
   selz=rz;                          
         delete(f6);
         delete(f7);
         delete(f8);
end
rep = questdlg(sprintf('Current Values  LL= %5.0f  UL = %5.0f \n Do you want to change values ? ',LL,UL),'Possibility to change Upper and Lower Limits','No'); % is USER OK with current selection?
  if isempty(rep)
         disp('User has cancelled.')
         return

end
        repi=0;               %reset of repi to ensure that the 2nd while loop will be re asked
end

%Creation of result files                    % user choice on export formats (BMP, JPG or TIF)
BMP = 0;
JPG = 0;
TIF = 0;
bit_depth = 255;
ans1 = questdlg('Export images as BMP or JPG?', ...  % Selection of output format
	'BMP, JPG or TIF?', ...
	'BMP','JPG','TIF','BMP');
 if isempty(ans1)
    disp('User has cancelled.')
    return
 end
                                            %selection of view (XYZ)
ans2 = questdlg('Axial view selecte?', ...
	'X, Y or Z?', ...
	'X','Y','Z','Z');
if isempty(ans2)
    disp('User has cancelled.')
    return
end
% save combined figures
new_folder=repr_folder;
f11 = figure(11);
  montage({I33,I44,I55,I333,I444,I555});       %display  initial ,mask of changes and changed pictur
   title(['Transformation of selected images nr ',num2str(selx),'/',num2str(sely),'/',num2str(selz),'  dataset:' , fileN,'  LL= ',num2str(LL),'  UL= ',num2str(UL)]);
   set_fig_position(f11,pc_height*.3,0, pc_height*.5, pc_width*0.75);
   
 f12 = figure(12);                              %display 2 histo -compares 2D based histo
  subplot(2,1,1)
    imhist(I1_repr,bit_depth);
    xline(LL,'-.b','LL'); xline(UL,'-.b','UL');
    title("Histogram Selected 2D")
  subplot(2,1,2)
    plot(hist3D);   
    title("test") 
    %and 3-D based histo
    xline(LL,'-.b','LL'); xline(UL,'-.b','UL');
    title(["Histogram 3D - DataSet:  ",fileN,' LL= ',num2str(LL),'  UL= ',num2str(UL) ]);
    set_fig_position(f12,pc_height*.3, 0, pc_height*.7, pc_width);
          
f13 = figure(13);
    imshow(I33ce);
    title (['Threshold value',oThre]);
    
if CANNY       
    f14 = figure('name','Mask and canny views');% diplay Canny and Mask views
     %title(['Mask and Canny Edge views for  ',num2str(selx),'/',num2str(sely),'/',num2str(selz),oThre])
     montage({mmask(I33),mmask(I44),mmask(I55),I33ce,I44ce,I55ce});  %displays mask views and Sobel/Canny
     title(['Mask and Canny Edge views for  ',num2str(selx),'/',num2str(sely),'/',num2str(selz),'CannyThres. = ',num2str(oThre),' LL= ',num2str(LL),'  UL= ',num2str(UL) ]);
     set_fig_position(f14,pc_height*.3, 0, pc_height*.7, pc_width);
    savefig(f14,fullfile(new_folder,'a-XYZcanny.fig')); %Writing of the image    
end

savefig(f11, fullfile(new_folder,'a-XYZimages.fig')); %Writing of the image
saveas(f11,  fullfile(new_folder,'a-XYZimages.png')); %Writing of the image
savefig(f12, fullfile(new_folder,'a-histograms.fig')); %Writing of the image
saveas(f12,  fullfile(new_folder,'a-histograms.png')); %Writing of the image
saveas(f13,  fullfile(new_folder,'a-canny.png')); %Writing of the image

delete(f12);                % delete old views 
delete(f13);
delete(f14);

%Write parameters to log file
out = fullfile(new_folder,'a-3D-LogFile.txt');

fileID = fopen(out,'w');
fprintf(fileID, ' ===================================================================\r\n');
fprintf(fileID, '||                         LOG INFORMATION                        ||\r\n');
fprintf(fileID, ' ===================================================================\r\n\r\n\r\n');
fprintf(fileID, 'Dataset folder: \t%s\r\n', repr_folder);
fprintf(fileID, 'Representative image of dataset: \t%s\r\n', fileN); % "sObject" is a string.
fprintf(fileID, 'Matrix Dimentions (XYZ): \t%i \t%i \t%i \r\n', nix,niy,niz);
fprintf(fileID, ' ===================================================================\r\n\r\n');
fprintf(fileID, 'Method of choice: \t%s\r\n', Chosen_method); 
fprintf(fileID, 'Lower Threshold percentage (in percents): \t%f\r\n', a*100); 
fprintf(fileID, 'Upper Threshold percentage (in percents): \t%f\r\n', au*100);
fprintf(fileID, 'lower Buffer: \t%f\r\n\r\n', buffer);
fprintf(fileID, 'Upper Buffer: \t%f\r\n\r\n', bufferu);
fprintf(fileID, 'Lower limit of grey value: \t%f\r\n', LL); 
fprintf(fileID, 'Upper limit of grey value: \t%f\r\n', UL); 
fprintf(fileID, 'Canny Edge Upper Threshold: \t%f\r\n', oThre); 
fprintf(fileID, 'Retained voxel share = %3.3f  Under LL = %3.3f above UL %3.3f \r\n',Ra3D,L3D,U3D);
fprintf(fileID, '   (retained voxel share is percent or Voxels between LL and UL\r\n');
fprintf(fileID, '       Retained object volume = %d  cubic meter\r\n',R3D*Voxelvol);
fprintf(fileID, ' ===================================================================\r\n\r\n');
fclose(fileID); % Close Logfile.
sprintf('Retained object volume = %d', R3D*Voxelvol')


fext=append(ans1,ans2) ;                     %name of the directory (BMPX)

mkdir(repr_folder, fext);
new_folder=fullfile(repr_folder,fext); %pathway to new folder
% Handle format as requested
    switch ans1
        case 'BMP'
            BMP = 1;
        case 'JPG'
            JPG = 1;
        case 'TIF'
            TIF = 1;
            bit_depth = 65535;              %TIF have a different depth
         case 'Cancel'
             message = msgbox('Operation has been cancelled by the user', 'Cancelled');
             return
    end
close all
width = UL-LL;
wait = waitbar(0,'Images are being processed...', ...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(wait,'canceling',0);
  
switch ans2                                 %selects number of image to create
    case "X"
        X=1;
        number_images=nix;
    case "Y"
        Y=1;
        number_images=niy;
    case "Z"
        Z=1;
        number_images=niz; 
end
tic

for n=1:number_images        %Process dataset image by image
    if getappdata(wait,'canceling')
        message = msgbox('Operation has been cancelled by the user', 'Cancelled');
        return;
    end
    switch ans2
        case "X"
            I1=reshape(IT(n,:,:),niy,niz);    %reshape from 3D to 2D
        case "Y"
            I1=reshape(IT(:,n,:),nix,niz);
        case "Z"
            I1=reshape(IT(:,:,n),nix,niy);
     end
    
    I1_d = double(I1);                        %histo windowng
    I2_d = (I1_d-LL)/width*bit_depth;
    baseFilename= strcat(fileN,num2str(n));

    %Write modified images in right format
    if BMP
        I2 = uint8 (I2_d);
        imwrite(I2, fullfile(new_folder,strcat(baseFilename,'.bmp'))); %Writing of the image
    
    elseif JPG
        I2 = uint8 (I2_d);
        imwrite(I2, fullfile(new_folder,strcat(baseFilename,'.jpg'))); %Writing of the image
    
    elseif TIF
        I2 = uint16 (I2_d);
        imwrite(I2, fullfile(new_folder,strcat(baseFilename,'.tif'))); %Writing of the image
    end
    waitbar(n/number_images,wait)
end
delete(wait);
delete(f8) ;               % delete old views 
delete(f6);
delete(f7);
toc
end
% do we need to create log files ?
close all
message = msgbox('Operation Completed!', 'Success');

return
%
%======= winhist (Histo Windowing)==============
function winhist=winhist(Im,LL,UL,bitdep)
        if bitdep==0
            bitdep=65535;
       end
        width1=UL-LL;
        Im0 = double(Im);%load with windowed data
        Im1 = (Im0-LL)/width1*bitdep;
        winhist=uint16(Im1);
        
  end
%===== mhisteq Histo Equalizer special ======= 
 % contrast improvement (histo equalizer modified)
%    - remove excluded value for equalizer
%    - put a 1000 margin between excluded and retained values)
function mhisteq=mhisteq(Im)
I2_ci = Im;
[m,n] = size(Im);
U= imhist(Im, 65535);
hist_size2 = size(U,1);

U(1)=1000;              %allow to differentiate excluded values (back)
for i = 2:(hist_size2-1)
     U(i) = U(i-1) + U(i);        
end
U(hist_size2)=U(hist_size2-1)+1000; % differentiate excluded val (White)
U=U/U(hist_size2)*65535;

for i = 1:m
    for j=1:n
        I2_ci(i,j) = U(Im(i,j)+1);      % gray for all values in range (LL/UL)
    end
end
   mhisteq=I2_ci;
end
%===mmask (Mask of remaining values===============
function mmask = mmask(Im)
[m,n] = size(Im);
ni=n;
for i = 1:m
    for j=1:n
        if ((Im(i,j) < 65535) && (Im(i,j) > 1))
         Im(i,j) = 32000;% gray for all values in range (LL/UL)
 
        end
    end
end
mmask=Im;
end
%===== Set-fig=========
function set_fig_position(h, top, left, height, width)
% Matlab has a wierd way of positioning figures so this function
% simplifies the poisitioning scheme in a more conventional way.
%
% Usage:      SET_FIG_POSITION(h, top, left, height, width);
%
%             H is the handle to the figure.  This can be obtain in the 
%               following manner:  H = figure(1);
%             TOP is the "y" screen coordinate for the top of the figure
%             LEFT is the "x" screen coordinate for the left side of the figure
%             HEIGHT is how tall you want the figure
%             WIDTH is how wide you want the figure
%
% Author: sparafucile17

% PC's active screen size
screen_size = get(0,'ScreenSize');
pc_width  = screen_size(3);
pc_height = screen_size(4);

%Matlab also does not consider the height of the figure's toolbar...
%Or the width of the border... they only care about the contents!
toolbar_height = 5;
window_border  = 0;

% The Format of Matlab is this:
%[left, bottom, width, height]
m_left   = left + window_border;
m_bottom = pc_height - height - top - toolbar_height - 1;
m_height = height;
m_width  = width - 1;

%Set the correct position of the figure
set(h, 'Position', [m_left, m_bottom, m_width, m_height]);

%If you want it to print to file correctly, this must also be set
% and you must use the "-r72" scaling to get the proper format
% set(h, 'PaperUnits', 'points');
% set(h, 'PaperSize', [width, height]); 
% set(h, 'PaperPosition', [0, 0, width, height]); %[ left, bottom, width, height]
end
%create text file
% Function to print the log for 3D
