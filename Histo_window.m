
clear all
close all
%%

%dialog for representative image of the dataset

[repr_image, repr_folder] = uigetfile('*.tif', 'Select a representative image of the dataset (.tif).');
if isequal(repr_image,0)
   disp('User selected Cancel');
   return
else
file_names = dir(fullfile(repr_folder, '*.tif'));
number_images = numel(file_names);
end

I1_repr = imread(fullfile(repr_folder,repr_image));

hist = imhist(I1_repr, 65535);

%Initilization of parameter
LL = 1;
UL = size(hist,1);
buffer = 0;
Chosen_method = 'canceled';
a = nan;


%%
%Ask user to adapt dynamic range automatically, or use manual input.
answer1 = questdlg('adapt dynamic range automatically or enter dynamic range manually', ...
	'Automatic or manual', ...
	'Adapt automatically','Enter manually','Adapt automatically');
% Handle response
switch answer1
    case 'Adapt automatically'
        
        Chosen_method = 'Adapt automatically';
        %Ask user for threshold value (percentage)
        prompt = {'Enter threshold value (in percentage)'};
        dlgtitle = 'Enter threshold value (Default: 0.1%)';
        definput = {'0.1'};
        dims = [1 65];
        opts.Resize= 'on';
        answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
        a = str2num(answer{1})/100;
        
        %Ask user if you want to add a buffer
        prompt = {'Potentially, enter value for buffer (in grey value 0 - 65535)'};
        dlgtitle = 'Enter buffer value (Default: 0)';
        definput = {'0'};
        dims = [1 65];
        opts.Resize= 'on';
        answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
        buffer = str2num(answer{1});
        
        %Calculate LL
        hist_sum = sum(hist);
        hist_size = size(hist,1);
        L_sum = 0;
        U_sum = 0;
        for n = 1:hist_size
            L_sum = L_sum + hist(n);
            if L_sum > a*hist_sum
                LL = n-buffer;
                break
            end
        end

        % Calculate UL
        for n = 1:(hist_size-1)
            i = hist_size-n;
            U_sum = U_sum + hist(i);
            if U_sum > a*hist_sum
                UL = hist_size-n+buffer;
                break
            end
        end

        LL;
        UL;
       
    
    case 'Enter manually'
        
        answer2 = questdlg('Enter dynamic range on the histogram or enter numerically?', ...
        'On histogram or numerically', ...
        'On histogram','Enter numerically','On histogram');
        %Ask user to indicate lower and upper limit on histogram
        
        switch answer2
            case 'On histogram'
                Chosen_method = 'Enter manually on histogram';
                figure 
                imhist(I1_repr, 65535)
                title('First select lower limit, then select upper limit')
                [x,y] = ginput(2);
                
                LL = x(1);
                UL = x(2);
                close (gcf)
                
                
            case 'Enter numerically'
                Chosen_method = 'Enter manually - numeric values';
                prompt = {'Enter grey value of lower limit: ','Enter grey value of upper limit: '};
                dlgtitle = 'Numerical input of thresholds';
                dims = [1 65];
                answer = inputdlg(prompt,dlgtitle,dims);

                LL = str2num(answer{1});
                UL = str2num(answer{2});
        end
end

                
%%   Calculations (example image to let user decide to continue or not)
        

width = UL-LL;
I1_d = double(I1_repr);
I2_d = (I1_d-LL)/width*65535;
I2_16 = uint16(I2_d);

screen_size = get(0,'ScreenSize');
pc_width  = screen_size(3);
pc_height = screen_size(4);

f1 = figure(1);
montage({I1_repr,I2_16});
set_fig_position(f1,0, 0, pc_height, pc_width*0.65);

f2 = figure(2);
subplot(2,1,1)
imhist(I1_repr,65535);
xline(LL,'-.b','LL'); xline(UL,'-.b','UL');
axis 'auto y'
subplot(2,1,2)
imhist(I2_16,65535);
axis 'auto y'
set_fig_position(f2,0, pc_width*0.65, pc_height, pc_width*0.35);





%% Save images as 8-bit
% f3 = figure(3);
% imshow(I2)

% %Make a folder 'Normalized'
% mkdir(im_folder, 'Normalized')
% 
% norm_folder=fullfile(im_folder,'Normalized'); %pathway to Normalized

% Ask user to export as BMP or as JPG
BMP = 0;
JPG = 0;
TIF = 0;
answer = questdlg('Export images as BMP or JPG?', ...
	'BMP, JPG or TIF?', ...
	'BMP','JPG','TIF','BMP');

if isempty(answer)
    disp('User has cancelled.')
    return
end


% Handle response
    switch answer
        case 'BMP'
            BMP = 1;
            mkdir(repr_folder, 'BMP');
            new_folder=fullfile(repr_folder,'BMP'); %pathway to new folder
        case 'JPG'
            JPG = 1;
            mkdir(repr_folder,'JPG');
            new_folder=fullfile(repr_folder,'JPG'); %pathway to new folder
        case 'TIF'
            TIF = 1;
            mkdir(repr_folder,'TIF');
            new_folder=fullfile(repr_folder,'TIF'); %pathway to new folder
%         case 'Cancel'
%             message = msgbox('Operation has been cancelled by the user', 'Cancelled');
%             return
    end

close all

if BMP || JPG
    bit_depth = 255;
    
    elseif TIF
    	bit_depth = 65535;
end

width = UL-LL;

wait = waitbar(0,'Images are being processed...', ...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

setappdata(wait,'canceling',0);


tic
for n=1:number_images
    if getappdata(wait,'canceling')
        message = msgbox('Operation has been cancelled by the user', 'Cancelled');
        return
    end
    
    I1 = imread(fullfile(repr_folder,file_names(n).name));
    I1_d = double(I1);
    I2_d = (I1_d-LL)/width*bit_depth;
    
    
    baseFilename= file_names(n).name;  %Naming of the image
    fullFilename=fullfile(repr_folder,baseFilename);
    [folder, baseFilename, extension] = fileparts(fullFilename);
    

    %Process dataset image by image
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
    
    %Save processed image that was chosen as representative image by the
    %user
    if strcat(baseFilename,'.tif') == repr_image
        I2_repr = I2;
    end   
    
    waitbar(n/number_images,wait)
    
end
delete(wait)

toc

% Make figures to compare starting image and resulting image
f3 = figure(3);
montage({I1_repr,I2_repr});
set_fig_position(f3,0, 0, pc_height, pc_width*0.65);

f4 = figure(4);
subplot(2,1,1)
imhist(I1_repr,65535);
xline(LL,'-.b','LL'); xline(UL,'-.b','UL');
axis 'auto y'
subplot(2,1,2)
imhist(I2_repr,bit_depth);
axis 'auto y'
set_fig_position(f4,0, pc_width*0.65, pc_height, pc_width*0.35);



  


%%
%Creating log file

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
fprintf(fileID, 'Threshold percentage (in percents): \t%f\r\n', a*100); 
fprintf(fileID, 'Buffer: \t%f\r\n\r\n', buffer); 
fprintf(fileID, 'Lower limit of grey value: \t%f\r\n', LL); 
fprintf(fileID, 'Upper limit of grey value: \t%f\r\n', UL); 
fprintf(fileID, ' ===================================================================\r\n\r\n');



%fprintf(fileID, 'Value\t%f\r\n', value); % "value" is a float.
fclose(fileID); % Close file.

saveas(f3, fullfile(log_folder,'comparison_images.png')); %Writing of the image
savefig(f3, fullfile(log_folder,'comparison_images.fig')); %Writing of the image
saveas(f4, fullfile(log_folder,'comparison_histograms.png')); %Writing of the image 
savefig(f4, fullfile(log_folder,'comparison_histograms.fig')); %Writing of the image   


message = msgbox('Operation Completed!', 'Success');

  %% Functions
function set_fig_position(h, top, left, height, width)
% Matlab has a wierd way of positioning figures so this function
% simplifies the poisitioning scheme in a more conventional way.
%
% Usage:      SET_FIG_POISITION(h, top, left, height, width);
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
