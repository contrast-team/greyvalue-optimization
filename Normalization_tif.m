clear all
%% Inititalization and Dataset for reference
% Initialization
flag_manual_entry = 0;
rescale_tif_ref = 1;
rescale_tif_im = 255/65535;
ext_ref = 'Manual entry of reference data';
bit_depth_ref=255;
bit_depth_im=65535;
bit_depth_mod=0;
Mat1 = 'Material 1';
Mat2 = 'Material 2';

%Ask user which reference materials are used
prompt = {'Enter name of reference material 1:',...
          'Enter name of reference material 2:'};
        dlgtitle = 'Reference materials';
        dims = [1 45];
        str = inputdlg(prompt,dlgtitle,dims);
        
        Mat1 = str{1};
        Mat2 = str{2};
%%
%Ask user to insert values used previously, or select a reference dataset.
answer = questdlg('Enter grey values manually or choose reference dataset', ...
                  'Reference input', ...
                  'Choose reference dataset',...
                  'Enter manually', 'Choose reference dataset');
% Handle response
    switch answer
    case 'Enter manually'
        flag_manual_entry = 1;
        prompt = {'Enter name of reference dataset',...
                  ['Enter grey value of ', Mat1,': '],...
                  ['Enter grey value of ', Mat2,': ']};
        dlgtitle = 'Reference input (8-bit)';
        dims = [1 45];  
        answer = inputdlg(prompt,dlgtitle,dims);

        ref_folder = ['Manual entry of gray values from: ' answer{1}];
        a_ref = str2num(answer{2});
        a_ref_original = a_ref;
        b_ref = str2num(answer{3});
        b_ref_original = b_ref;
    
    case 'Choose reference dataset'
        
        [repr_image_ref_name, ref_folder] = uigetfile({'*.tif';'*.bmp';'*.jpg'}, ...
            'Select a representative image of the reference dataset');
        
        repr_image_ref_path = fullfile(ref_folder,repr_image_ref_name);
        [~,~,ext_ref] = fileparts(repr_image_ref_path);
        repr_image_ref = imread(repr_image_ref_path);
        
        if ext_ref == '.tif'
            rescale_tif_ref = 255/65535;
            bit_depth_ref = 65535;
        end   
        ref_file_names = dir(fullfile(ref_folder, ext_ref));

        % Initialiation of the grey values [reference]
        figure
        imshow(repr_image_ref)
        title (['Select 2 points in: ', Mat1,...
               '. A line will be plotted between these points.'])
        [xa1,ya1] = ginput(2); % choose 2 points in material 1.
        title (['Select 2 points in: ', Mat2,...
               '. A line will be plotted between these points.'])
        [xb1,yb1] = ginput(2); % choose 2 points in material 2.
        close all

            
        
        profile_a_ref = (improfile(repr_image_ref,xa1,ya1)); % calculate grey value profile between 2 points
        profile_b_ref = improfile (repr_image_ref,xb1,yb1); % calculate grey value profile between 2 points
        a_ref_original = mean(profile_a_ref,'all'); % average along the profile
        b_ref_original = mean(profile_b_ref,'all'); % average along the profile
        a_ref = a_ref_original*rescale_tif_ref; % Rescales to 8-bit (if ref = tif)
        b_ref = b_ref_original*rescale_tif_ref; % Rescales to 8-bit (if ref = tif)
        
end

%% Dataset for normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[repr_image_im_name, im_folder] = uigetfile({'*.tif'}, ...
            'Select a representative image of the dataset to be normalized (.tif).');

repr_image_im_path = fullfile(im_folder,repr_image_im_name);
[~,repr_image_im_noext,ext_im] = fileparts(repr_image_im_path);
repr_image_im = imread(repr_image_im_path);

im_file_names = dir(fullfile(im_folder, '*.tif'));

rescale_tif_im = 255/65535;
total_im = numel(im_file_names); % Count total number of images present in that folder


% Initialization of grey values [images before normalization]
figure
imshow(repr_image_im)
title (['Select 2 points in: ', Mat1,...
               '. A line will be plotted between these points.'])
[xa2,ya2] = ginput(2);  % choose 2 points in material A.
title (['Select 2 points in: ', Mat2,...
               '. A line will be plotted between these points.'])
[xb2,yb2] = ginput(2);  % choose 2 points in material B.
close all

%Note: 'im' is always tif format
profile_a_im = improfile(repr_image_im,xa2,ya2);   % calculate grey value profile between 2 points
profile_b_im = improfile (repr_image_im,xb2,yb2);  % calculate grey value profile between 2 points
a_im_16 = mean(profile_a_im,'all');  % average along the profile
b_im_16 = mean(profile_b_im,'all');  % average along the profile

Rescale = (255/65535);
a_im_8 = a_im_16 * rescale_tif_im;
b_im_8 = b_im_16 * rescale_tif_im;



%Make a folder 'Normalized'

mkdir(im_folder, 'Normalized');
norm_folder=fullfile(im_folder,'Normalized'); %pathway to Normalized

answer = questdlg('Export images as BMP, JPG or TIF?', ...
	'BMP, JPG or TIF?', ...
	'BMP','JPG','TIF','BMP');
% Handle response
    switch answer
        case 'BMP'
            ext_mod = '.bmp';
            bit_depth_mod = 255;
            
        case 'JPG'
            ext_mod = '.jpg';
            bit_depth_mod = 255;
            
        case 'TIF'
            ext_mod = '.tif';
            bit_depth_mod = 65535;
            
    end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Factor = abs(a_ref-b_ref)/abs(a_im_8-b_im_8); %Factor = difference (a,b)_ref / (a,b)_im

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate the repr image after normalization

repr_image_im_d = double(repr_image_im); %change from uint16 to double class
    
    if ext_mod == '.tif'
        repr_image_mod_d = (repr_image_im_d*Factor + a_ref/rescale_tif_im - a_im_16*Factor);  %apply formula
        repr_image_mod = uint16 (repr_image_mod_d);  %convert back from double to uint16
    else
        repr_image_mod_d = (repr_image_im_d*Factor*rescale_tif_im + a_ref - a_im_8*Factor);  %apply formula
        repr_image_mod = uint8 (repr_image_mod_d);  %convert back from double to uint8
    end
    
    
    
    
 %% Intermediate check
 
 %Check output
    %Calculate profiles on the normalized images
profile_a_mod = improfile(repr_image_mod,xa2,ya2);
profile_b_mod = improfile (repr_image_mod,xb2,yb2);
a_mod_original = mean(profile_a_mod,'all');
b_mod_original = mean(profile_b_mod,'all');

if ext_mod == '.tif'
    a_mod = a_mod_original*Rescale;
    b_mod = b_mod_original*Rescale;
else
    a_mod = a_mod_original;
    b_mod = b_mod_original;
end

    

% Plot showing reference , before normalization and after normalization.
% Together with a bar graph showing the grey values.

%Calculate screen size to get height
screen_size = get(0,'ScreenSize');
pc_width  = screen_size(3);
pc_height = screen_size(4);

%make figure with four subplots
f1=figure(1);
    subplot (2,2,1)
        if flag_manual_entry == 0 %reference dataset was chosen
            imshow(repr_image_ref, 'border', 'tight')
            title('Reference image')
            hold on
            line([xa1],[ya1],'LineWidth',2)
            line([xb1],[yb1],'LineWidth',2,'Color','r')
        else %display text that no image is available
            text(0.1,0.5,['No image available since reference'...
                 newline  'gray values were entered manually'],...
                 'Fontsize', 14); axis off
        end

    subplot(2,2,2)
        imshow(repr_image_mod, 'border', 'tight')
        title ('Image after normalization')
        hold on
        line([xa2],[ya2],'LineWidth',2)
        line([xb2],[yb2],'LineWidth',2,'Color','r')
    
    subplot(2,2,3)
        imshow(repr_image_im, 'border', 'tight')
        title ('Image before normalization')
        hold on
        line([xa2],[ya2],'LineWidth',2)
        line([xb2],[yb2],'LineWidth',2,'Color','r')
 
 
 %%
    
    
    
    
    
    
    
    
% Calculate other images

%Waitbar initializing
wait = waitbar(0,'Images are being processed...', ...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

setappdata(wait,'canceling',0);
flag_cancelled = 0;

tic
%Calculating normalized images       
for n = 1:total_im
    if getappdata(wait,'canceling')
        flag_cancelled=1;
        break
    end
    
    f = fullfile(im_folder, im_file_names(n).name);
    im_f = imread(f); %read current image
    im_f_d = double(im_f); %change from uint16 to double class
    
    if ext_mod == '.tif'
        im_f_d_mod = (im_f_d*Factor + a_ref/rescale_tif_im - a_im_16*Factor);  %apply formula
        im_f_mod = uint16 (im_f_d_mod);  %convert back from double to uint16
    
    else
        im_f_d_mod = (im_f_d*Factor*rescale_tif_im + a_ref - a_im_8*Factor);  %apply formula
        im_f_mod = uint8 (im_f_d_mod);  %convert back from double to uint8
    end
    
    [~,mod_noext,~]= fileparts(f);
    im_f_mod_path = fullfile(norm_folder,strcat(mod_noext,ext_mod));
    imwrite (im_f_mod, im_f_mod_path);
    
    waitbar(n/total_im,wait)

end

toc
delete(wait)

%%
%Check output
    %Calculate profiles on the normalized images
profile_a_mod = improfile(repr_image_mod,xa2,ya2);
profile_b_mod = improfile (repr_image_mod,xb2,yb2);
a_mod_original = mean(profile_a_mod,'all');
b_mod_original = mean(profile_b_mod,'all');

if ext_mod == '.tif'
    a_mod = a_mod_original*Rescale;
    b_mod = b_mod_original*Rescale;
else
    a_mod = a_mod_original;
    b_mod = b_mod_original;
end

    

% Plot showing reference , before normalization and after normalization.
% Together with a bar graph showing the grey values.

%Calculate screen size to get height
screen_size = get(0,'ScreenSize');
pc_width  = screen_size(3);
pc_height = screen_size(4);

%make figure with four subplots
f1=figure(1);
    subplot (2,2,1)
        if flag_manual_entry == 0 %reference dataset was chosen
            imshow(repr_image_ref, 'border', 'tight')
            title('Reference image')
            hold on
            line([xa1],[ya1],'LineWidth',2)
            line([xb1],[yb1],'LineWidth',2,'Color','r')
        else %display text that no image is available
            text(0.1,0.5,['No image available since reference'...
                 newline  'gray values were entered manually'],...
                 'Fontsize', 14); axis off
        end

    subplot(2,2,2)
        imshow(repr_image_mod, 'border', 'tight')
        title ('Image after normalization')
        hold on
        line([xa2],[ya2],'LineWidth',2)
        line([xb2],[yb2],'LineWidth',2,'Color','r')
    
    subplot(2,2,3)
        imshow(repr_image_im, 'border', 'tight')
        title ('Image before normalization')
        hold on
        line([xa2],[ya2],'LineWidth',2)
        line([xb2],[yb2],'LineWidth',2,'Color','r')

    subplot(2,2,4)
        X = categorical({Mat1, Mat2});
        X = reordercats(X,{Mat1, Mat2});
        Y = round([a_ref  a_im_8   a_mod; b_ref   b_im_8   b_mod],1);
        b = bar(X,Y);
        legend ({'reference', 'before normalization', 'after normalization'}, 'Location','northwest')

        xtips1 = b(1).XEndPoints;
        ytips1 = b(1).YEndPoints;
        labels1 = string(b(1).YData);
        text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')

        xtips2 = b(2).XEndPoints;
        ytips2 = b(2).YEndPoints;
        labels2 = string(b(2).YData);
        text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')

        xtips3 = b(3).XEndPoints;
        ytips3 = b(3).YEndPoints;
        labels3 = string(b(3).YData);
        text(xtips3,ytips3,labels3,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')

        title('Average grey value')
        set_fig_position(f1,0, 0, pc_height, pc_width*0.6);
        
  %%   Plot histograms   
        f2 = figure(2);

            subplot(3,1,1)
            if flag_manual_entry
                text(0.1,0.5,['No histogram available since reference'...
                 newline  'gray values were entered manually'...
                 newline ['Grey value of ', Mat1,': ', num2str(a_ref)]...
                 newline ['Grey value of ', Mat2,': ', num2str(b_ref)]],...
                 'Fontsize', 14); axis off
            else
            imhist(repr_image_ref, bit_depth_ref)
            axis 'auto y'
            title('reference image')
            xline(a_ref_original,'-.b',Mat1); xline(b_ref_original,'-.b',Mat2);
            end

            subplot(3,1,2)
            imhist(repr_image_im, bit_depth_im)
            axis 'auto y'
            title('Before normalization')
            xline(a_im_16,'-.b',Mat1); xline(b_im_16,'-.b',Mat2);
            
            subplot(3,1,3)
            imhist(repr_image_mod, bit_depth_mod)
            axis 'auto y'
            title('After normalization')
            xline(a_mod_original,'-.b',Mat1); xline(b_mod_original,'-.b',Mat2);
            
        set_fig_position(f2,0, pc_width*0.6, pc_height, pc_width*0.4);

%%
%Creating log file

%create folder and text file
file_name = 'Normalization.txt';
%create log folder
mkdir(norm_folder, 'LOG files');
log_folder=fullfile(norm_folder, 'LOG files');
out = fullfile(log_folder,file_name);

%Write parameters to log file

fileID = fopen(out,'w+');
fprintf(fileID, ' ===================================================================\r\n');
fprintf(fileID, '||                         LOG INFORMATION                        ||\r\n');
fprintf(fileID, ' ===================================================================\r\n\r\n\r\n');
fprintf(fileID, 'Dataset: \t%s\r\n', mod_noext); % "sObject" is a string.
fprintf(fileID, 'Number of images: \t%i\r\n', total_im);
fprintf(fileID, 'Reference dataset: \t%s\r\n', ref_folder); 
fprintf(fileID, 'Reference Material 1: \t%s\r\n', Mat1);
fprintf(fileID, 'Reference Material 2: \t%s\r\n', Mat2);
fprintf(fileID, ' ===================================================================\r\n\r\n');
fprintf(fileID, 'Reference image type: \t%s\r\n', ext_ref); 
fprintf(fileID, 'Before normalization image type: \t%s\r\n', ext_im);
fprintf(fileID, 'After normalization image type: \t%s\r\n', ext_mod);
fprintf(fileID, ' ===================================================================\r\n\r\n');
fprintf(fileID, ['Reference gray value of ', Mat1, ' (original): \t%f\r\n'], a_ref_original); 
fprintf(fileID, ['Reference gray value of ', Mat2, ' (original): \t%f\r\n\r\n'], b_ref_original);
fprintf(fileID, ['Before normalization gray value of ', Mat1, ' (original): \t%f\r\n'], a_im_16); 
fprintf(fileID, ['Before normalization gray value of ', Mat2, ' (original): \t%f\r\n\r\n'], b_im_16); 
fprintf(fileID, ['After normalization gray value of ', Mat1, ' (original): \t%f\r\n'], a_mod_original); 
fprintf(fileID, ['After normalization gray value of ', Mat2, ' (original): \t%f\r\n\r\n\r\n'], b_mod_original);

fprintf(fileID, ['Reference gray value of ', Mat1, ' (8-bit): \t%f\r\n'], a_ref); 
fprintf(fileID, ['Reference gray value of ', Mat2, ' (8-bit): \t%f\r\n\r\n'], b_ref);
fprintf(fileID, ['Before normalization gray value of ', Mat1, ' (8-bit): \t%f\r\n'], a_im_16*Rescale); 
fprintf(fileID, ['Before normalization gray value of ', Mat2, ' (8-bit): \t%f\r\n\r\n'], b_im_16*Rescale);
fprintf(fileID, ['After normalization gray value of ', Mat1, ' (8-bit): \t%f\r\n'], a_mod); 
fprintf(fileID, ['After normalization gray value of ', Mat2, ' (8-bit): \t%f\r\n'], b_mod);


fprintf(fileID, ' ===================================================================\r\n\r\n');

%fprintf(fileID, 'Value\t%f\r\n', value); % "value" is a float.
fclose(fileID); % Close file.

saveas(f1, fullfile(log_folder,'comparison_images.png')); %Writing of the image
savefig(f1, fullfile(log_folder,'comparison_images.fig')); %Writing of the image
saveas(f2, fullfile(log_folder,'comparison_histograms.png')); %Writing of the image 
savefig(f2, fullfile(log_folder,'comparison_histograms.fig')); %Writing of the image

%%

%Text box indicating if operation was completed or canceled

if flag_cancelled == 0
message = msgbox('Operation Completed!', 'Success');
disp('Finished!')
else
message = msgbox('Operation has been canceled', 'Canceled');
end

set(message, 'position', [pc_height/1.777 pc_height/3 pc_height/6 pc_height/18]); %makes box bigger

% figure
% diff_I = imshow(ref_middle-mod_middle);

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







    
 