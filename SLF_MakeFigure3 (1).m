% Make FIgure 3
% we are asking how similar the different tracking protocols are
% we are going to load track density for all subjects, for all methods
% then loop through, and for each combination of method, calculate a Dice
% overlap coefficient for each pair
% we will plot this in a boxplot, with three boxes, one for each subject

% Hecht: A
% Reijmer: B

addpath(genpath('/Volumes/schillkg/MATLAB/NIFTI_20130306/'))
addpath(genpath('/Volumes/schillkg/MATLAB/Mathworks/NotBoxPlot/'))

%% we are going to load track density for all subjects, for all methods

tracks = {'A';'B';'C';'D';'E';'F';'G';'H';'I'}

TDI = zeros(157,189,136,2,length(tracks));

for i = 2:3
    for j = 1:length(tracks)
        % load tracks
        try
        nii = load_untouch_nii_gz(['subject' num2str(i) filesep tracks{j} '.nii.gz']);
        catch
        nii = load_untouch_nii(['subject' num2str(i) filesep tracks{j} '.nii']);
        end
        img = nii.img;
        disp(size(img))
        TDI(:,:,:,i-1,j) = img;
    end
end


%% then loop through, and for each combination of method, calculate a Dice

D1=[]; M1=[];
for i=1:length(tracks)-1
    disp(i)
    base_scan=squeeze(TDI(:,:,:,1,i));
    for j = i+1:1:9
        disp(j)
        comparing_to=squeeze(TDI(:,:,:,1,j));
        A = base_scan; A(A>0)=1;
        B = comparing_to; B(B>0)=1;
        [msd hd] = surface_distance(A, B, [1 1 1]);
        intersection = (A & B);
        a = sum(intersection(:));
        b = sum(A(:));
        c = sum(B(:));
        Dice = 2*a/(b+c);
        D1 = [D1 Dice];
        M1 = [M1 msd];
    end
end


D2=[]; M2=[];
for i=1:length(tracks)-1
    disp(i)
    base_scan=squeeze(TDI(:,:,:,2,i));
    for j = i+1:1:9
        disp(j)
        comparing_to=squeeze(TDI(:,:,:,2,j));
        A = base_scan; A(A>0)=1;
        B = comparing_to; B(B>0)=1;
        [msd hd] = surface_distance(A, B, [1 1 1]);
        intersection = (A & B);
        a = sum(intersection(:));
        b = sum(A(:));
        c = sum(B(:));
        Dice = 2*a/(b+c);
        D2 = [D2 Dice];
        M2 = [M2 msd];
    end
end

%%

tracks = {'A';'B';'C';'D';'E';'F';'G';'H';'I'}

TDI = zeros(181,218,181,1,length(tracks));

for i = 1
    for j = 1:length(tracks)
        % load tracks
        try
        nii = load_untouch_nii_gz(['subject' num2str(i) filesep tracks{j} '.nii.gz']);
        catch
        nii = load_untouch_nii(['subject' num2str(i) filesep tracks{j} '.nii']);
        end
        img = nii.img;
        disp(size(img))
        TDI(:,:,:,i,j) = img;
    end
end

D3=[]; M3=[];
for i=1:length(tracks)-1
    disp(i)
    base_scan=squeeze(TDI(:,:,:,1,i));
    for j = i+1:1:9
        disp(j)
        comparing_to=squeeze(TDI(:,:,:,1,j));
        A = base_scan; A(A>0)=1;
        B = comparing_to; B(B>0)=1;
        [msd hd] = surface_distance(A, B, [1 1 1]);
        intersection = (A & B);
        a = sum(intersection(:));
        b = sum(A(:));
        c = sum(B(:));
        Dice = 2*a/(b+c);
        D3 = [D3 Dice];
        M3 = [M3 msd];
    end
end

%%
y = [D3; D1; D2]
figure;
notBoxPlot(y')
xlabel('Subject')      
ylabel('Dice Coefficient')
set(gca,'FontSize',24); %grid on;
ylim([0 1])
whitebg(gcf,[1 1 1]);         % background color (in+out plot)
set(gcf,'Color',[1,1,1]);
saveas(gcf,'DiceCoefficients.png')


y = [M3; M1; M2]
figure;
notBoxPlot(y')
xlabel('Subject')      
ylabel('Mean Surface Distance (mm)')
set(gca,'FontSize',24); %grid on;
whitebg(gcf,[1 1 1]);         % background color (in+out plot)
set(gcf,'Color',[1,1,1]);
saveas(gcf,'MeanSurfaceDistance.png')



