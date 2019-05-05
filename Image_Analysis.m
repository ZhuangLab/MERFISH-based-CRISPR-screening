%% Measure the MERFISH and cell identities
%% Tian LU
%% 3/1/2017
%% Setup path and parameters
% Define data path
dataPath='';
analysisSavePath = SetFigureSavePath([dataPath,'\FOV\'], ...
     'makeDir', true);
convertPath=SetFigureSavePath([dataPath,'\Convert1\'], ...
     'makeDir', true);
RawPath = ''; 
barcodePath=SetFigureSavePath([dataPath,'\barcode\'], ...
     'makeDir', true);
% Useful data structure for spotfinding
[~,~,all]=xlsread([RawPath 'Arrangement.xlsx']);
Gene=all(2:end,1);
color=all(2:end,2);
desiredframe=all(2:end,4);
threshold=all(2:end,3);
desiredround=all(2:end,5);
filetype=all(2:end,6);
colocalizationthreshold=3; %unit: pixel
bitnumber=36;
roundnumber=18;
testnum=200;
mapPath=SetFigureSavePath([analysisSavePath,'\map_info_with_MALAT1_2\'], ...
     'makeDir', true);
decodePath=SetFigureSavePath([analysisSavePath,'\decode_info\'], ...
     'makeDir', true);
%% ------------------------------------------------------------------------
% Start logging
%%-------------------------------------------------------------------------
if ~isempty(mfilename) % Only log if running as a script
    diaryFile = [analysisSavePath 'matlab_output.log']; % File name
    diary off; % Turn off diary if already in use
    if exist(diaryFile) % Delete existing file
        diary off;
        delete(diaryFile);
    end
    diary(diaryFile); % Set diary file
    diary on;

    % Display information
    PageBreak();
    display(['Running: ' mfilename('fullpath') '.m']);
    display(['Started: ' datestr(now)]);
    
    % Archive script
    copyfile( [mfilename('fullpath'),'.m'],[analysisSavePath,mfilename,'.m']);
    display('------------------------------------------------------------------');
    display(['Copied analysis script to ' analysisSavePath,mfilename,'.m']);
    
    % Start script timer
    scriptTimer = tic;
end
%% Correct the illumination of 647 and 750
    Path=[RawPath 'data\'];
    file=['Epi-750s1-650s1-560s1-488s2-405s1_(?<fov>[0-9]+)_00'];    
    tempFiles = BuildFileStructure(Path, ...
    'fileExt', 'dax', ...
    'regExp', file, ...
    'fieldNames', {'fov','round'}, ...
    'fieldConv', {@str2num});

    display(['Found ' num2str(length(tempFiles)) ' dax files']);

    % Run analysis of all fov
    tempsum=zeros(2048,2048);
    for f=1:length(tempFiles)
        data=ReadDax(tempFiles(f).filePath,'startFrame', 1, 'endFrame', 1); 
        tempsum=tempsum+double(data);
    end
    a=tempsum/length(tempFiles);
    amax=max(max(a));
    ratio750=a/amax;


    % Run analysis of all fov
    tempsum=zeros(2048,2048);
    for f=1:length(tempFiles)
        data=ReadDax(tempFiles(f).filePath,'startFrame', 2, 'endFrame', 2); 
        tempsum=tempsum+double(data);
    end
    a=tempsum/length(tempFiles);
    amax=max(max(a));
    ratio647=a/amax;


%% Create parallel pool
if isempty(gcp('nocreate'))
    p = parpool(15); % Set this number to control the number of parallel workers.
                     % Polite usage would suggest these maximum values: morgan = 20, cajal = 10
else
    p = gcp;
end

%% spotfinding
for i=1:length(Gene)
    analysisPath=[dataPath,'\',Gene{i},'\'];
     if ~exist(analysisPath)
     
        Path=[RawPath 'data\'];
        file=['Epi-750s1-650s1-560s1-488s2-405s1_(?<fov>[0-9]+)_',num2str(desiredround{i,1},'%02d')];    
       
        tempFiles = BuildFileStructure(Path, ...
        'fileExt', 'dax', ...
        'regExp', file, ...
        'fieldNames', {'fov','round'}, ...
        'fieldConv', {@str2num});
       display(['Found ' num2str(length(tempFiles)) ' dax files of ' Gene{i}]);
       tempinfos = BuildFileStructure(Path, ...
        'fileExt', 'inf', ...
        'regExp', file, ...
        'fieldNames', {'fov','round'}, ...
        'fieldConv', {@str2num});
        display(['Found ' num2str(length(tempinfos)) ' info files of ' Gene{i}]);
    
        %% Write a daoSTORM parameter
        analysisLabel = [num2str(color{i}) '_' num2str(threshold{i,1})];
        parametersName = ['daoSTORMParameters_' analysisLabel '.xml'];
        WriteDaoSTORMParameters([analysisPath parametersName], ...
            'start_frame', 0, ... % the desired 
            'max_frame', 1, ...
            'iterations', 1, ... % Do not fit overlapping molecules
            'threshold', threshold{i,1}, ...
            'sigma', 3, ... % The best guess for the PSF size in terms of pixels
            'pixel_size', 110);
        %% Run analysis of all fov
        parfor f=1: length(tempFiles)
            frameName = [analysisPath, '\','Bit_c1_' num2str(f,'%02d') '_' num2str(color{i}),'_',num2str(threshold{i,1}),'_mlist.bin'];
            
            temp=ReadDax(tempFiles(f).filePath, 'startFrame', desiredframe{i,1}, 'endFrame',desiredframe{i,1});
            if desiredframe{i,1}==1
              temp=int16(double(temp)./ratio750);
            else
              temp=int16(double(temp)./ratio647);
            end
            frameName=['Bit_c1_' num2str(f,'%02d')];
            infPath_FISH = ReadInfoFile(tempinfos(f).filePath);
            infPath_FISH.number_of_frames = 1;
            infPath_FISH.localName = [frameName '.inf'];
            infPath_FISH.localPath = [convertPath,'\'];  % there must be a '\' here
            infPath_FISH.file = [convertPath, '\', frameName '.dax'];
            WriteDAXFiles(temp, infPath_FISH);
            daoSTORM([convertPath, '\', frameName '.dax'], ... % Path to dax to analyze
                [analysisPath parametersName], ... % Path to daoSTORM configuration file
                'overwrite', true, ... % Overwrite any existing analysis
                'numParallel', 2, ... % Only process one file at a time
                'savePath', analysisPath, ... % Location of the mlist files
                'mListType', [analysisLabel '_mlist'], ... % A label to mark the analysis for each file
                'outputInMatlab', false); % Display the output of daoSTORM in matlab
            disp(f);
            
        end
     end
end

%% ------------------------------------------------------------------------
% Map image files for ease of loading
%%-------------------------------------------------------------------------
%% smFISH binary files
mlistfile={};
for i=1:length(Gene)
    file=['Bit_c1_(?<fov>[0-9]+)_',num2str(color{i}),'_',num2str(threshold{i,1}),'_mlist'];
    mlistfile{i} = BuildFileStructure([dataPath,'\',Gene{i},'\'], ...
    'regExp', file, ...
    'fileExt', 'bin', ...
    'fieldNames', {'fov'}, ...
    'fieldConv', {@str2num});
end
mlistfile561=mlistfile(bitnumber+1:end);
mlistfile=mlistfile(1:bitnumber);
%% ------------------------------------------------------------------------
% Load and display panels for each fov
%%-------------------------------------------------------------------------
%% decode the 561 spot in every FOV

map=[1,1,0,0,2,2,3,3,4,4,9,9,5,5,6,6,7,7,8,8,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17];
parfor fovID =1:length(mlistfile{1})-1
        %% count and plot the colocalization and record nearest distance from colocalization 
        %if exist(strcat(decodePath,['FOV_' num2str(fovID),'_decode.mat']),'file')==0
        spotlist={};
        spotC={};
        empty=0;
        for i=1:length(mlistfile)
            list=mlistfile{i};
            spotlist{i}=ReadMasterMoleculeList(list([list.fov]==fovID).filePath);
            spotC{i}=[spotlist{i}.xc,spotlist{i}.yc];
            if isempty(spotC{i})
                empty=1;
            end             
        end
        shift=zeros(2,roundnumber)
        spotlist561={};
        spotC561={};
        spotC561original={};
        for i=1:length(mlistfile561)
            list=mlistfile561{i};
            spotlist561{i}=ReadMasterMoleculeList(list([list.fov]==fovID).filePath);
            spotC561{i}=[spotlist561{i}.xc,spotlist561{i}.yc];
            spotC561original{i}=spotC561{i};             
        end   
        count=0;
        while isempty(spotC561{1}) & count<100
                spotC561{1}=spotC561{randi(roundnumber)};
                count=count+1;
        end
            %%  shift rounds
        decode=zeros(bitnumber,length(spotC561{1}),'int16');
        if count<100
            display(['FOV ',num2str(fovID),'Start shifting']);
            for rounds=1:roundnumber
                while  isempty(spotC561{rounds})
                    spotC561{rounds}=spotC561{randi(rounds+1)};
                end
                if rounds~=1
                    for x=1:10
                        dist=pdist2(spotC561{rounds},spotC561{1});
                        [M,I]=min(dist);
                        [temp,ids]=sort(M);
                        dev=diff(temp);
                        idd=(dev<median(dev));
                        temp1=spotC561{1}(ids(idd),:);
                        temp2=spotC561{rounds}(I(ids(idd)),:);
                        co=temp1-temp2;
                        coordi=[median(co(:,1)),median(co(:,2))];    
                        spotC561{rounds}=[spotC561{rounds}(:,1)+coordi(1,1),spotC561{rounds}(:,2)+coordi(1,2)];
                        for i=1:bitnumber                             
                                 if (map(i)==rounds-1)
                                   if  ~isempty(spotC{i})
                                   spotC{i}=[spotC{i}(:,1)+coordi(1,1),spotC{i}(:,2)+coordi(1,2)]; 
                                   end
                                end
                            
                        end
                    end
                    shift(1,rounds)=mean(spotC561{rounds}(:,1)-spotC561original{rounds}(:,1));
                    shift(2,rounds)=mean(spotC561{rounds}(:,2)-spotC561original{rounds}(:,2));
                end
            end
            %% get rid of overlapping 561 spots
           
            display(['FOV ',num2str(fovID),'Start removing the spots']);
            dist=pdist2(spotC561{1},spotC561{1});           
            dist(dist==0)=inf;
            [M,I]=min(dist);
            for ii=1:length(M)
                if M(ii)<5
                    decode(:,ii)=-1;
                end
            end
          
            %%  decode the smFISH spots    
            display('Start decoding');
            for j=1:bitnumber     
                if ~isempty(spotC{j})
                    dist=pdist2(spotC{j},spotC561{1});
                    [M,I]=min(dist);
                    for ii=1:length(M)
                        if (M(ii)<=colocalizationthreshold)
                            if (decode(j,ii)>-1)                            
                                decode(j,ii)=1; %decode the spot
                            end
                        end
                    end
                end
            end  
        end
        SpotCoordinates=spotlist561{1};
        display('Save Files');
        parsavedecode(strcat(decodePath,['FOV_' num2str(fovID),'_decode.mat']),decode,SpotCoordinates,shift);
        %end
end
%% Correct the illumination
    Path=[RawPath 'data\'];
    file=['Epi-750s1-650s1-560s1-488s2-405s1_(?<fov>[0-9]+)_00'];    
    tempFiles = BuildFileStructure(Path, ...
    'fileExt', 'dax', ...
    'regExp', file, ...
    'fieldNames', {'fov','round'}, ...
    'fieldConv', {@str2num});

    display(['Found ' num2str(length(tempFiles)) ' dax files']);

    % Run analysis of all fov
    tempsum=zeros(2048,2048);
    for f=1:length(tempFiles)
        data=ReadDax(tempFiles(f).filePath,'startFrame', 6, 'endFrame', 6); 
        tempsum=tempsum+double(data);
    end
    a=tempsum/length(tempFiles);
    amax=max(max(a));
    nucleusratio=a/amax;
    
    Path=[RawPath 'data\'];
    file=['Epi-750s1-650s1-560s1-488s2-405s1_(?<fov>[0-9]+)_08'];    
    tempFiles = BuildFileStructure(Path, ...
    'fileExt', 'dax', ...
    'regExp', file, ...
    'fieldNames', {'fov','round'}, ...
    'fieldConv', {@str2num});

    display(['Found ' num2str(length(tempFiles)) ' dax files']);

    % Run analysis of all fov
    tempsum=zeros(2048,2048);
    for f=1:length(tempFiles)
        data=ReadDax(tempFiles(f).filePath,'startFrame', 6, 'endFrame', 6); 
        tempsum=tempsum+double(data);
    end
    b=tempsum/length(tempFiles);
    %b=a-b;
    bmax=max(max(b));
    SONratio=b/bmax;

    Path=[RawPath 'data\'];
    file=['Epi-750s1-650s1-560s1-488s2-405s1_(?<fov>[0-9]+)_07'];    
    tempFiles = BuildFileStructure(Path, ...
    'fileExt', 'dax', ...
    'regExp', file, ...
    'fieldNames', {'fov','round'}, ...
    'fieldConv', {@str2num});

    display(['Found ' num2str(length(tempFiles)) ' dax files']);

    % Run analysis of all fov
    tempsum=zeros(2048,2048);
    for f=1:length(tempFiles)
        data=ReadDax(tempFiles(f).filePath,'startFrame', 5, 'endFrame', 5); 
        tempsum=tempsum+double(data);
    end
    c=tempsum/length(tempFiles);
    cmax=max(max(c));
    cytoratio=c/cmax;
    
save(([barcodePath 'illumination.mat']),'nucleusratio','cytoratio','SONratio'); 
%% Try to segment cells on each FOV
nucleusbwThresh1=0.6;
nucleusbwThresh2=0.4;
cytoplasmbwThresh=0.07;
mapPath=SetFigureSavePath([analysisSavePath,'\map_info_optimized\'], ...
     'makeDir', true);
parameterPath=[RawPath 'parameter1.xlsx'];
[~,~,all]=xlsread(parameterPath);
parameter=all(2:end,2:end); 
for fovID =1:length(mlistfile{1})-2
   % if ~exist(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']),'file')
        marker='o>h^<o>h^<o>h^<o>h^<o>h^<o>h^<o>h^<o>h^<o>h^<o>h^<o>h^<o>h^<o>h^<o>h^<o>h^<';
        % segmentation  
        display(['Start segmentation: FOV ',num2str(fovID)]);
        cyto=ReadDax([RawPath,'data\Epi-750s1-650s1-560s1-488s2-405s1_',num2str(fovID,'%03d'),'_07.dax'], 'startFrame', 5, 'endFrame', 5);
        Nucleus=ReadDax([RawPath,'data\Epi-750s1-650s1-560s1-488s2-405s1_',num2str(fovID,'%03d'),'_00.dax'], 'startFrame', 6, 'endFrame', 6);
       cyto=int16(double(cyto)./cytoratio);
        Nucleus=imadjust(int16(double(Nucleus)./nucleusratio));
        Nucleus=imadjust(Nucleus);
        figure;
        imshow(Nucleus);
%         figure;
%         imshow(cyto);
        
        nucleusImgsi = im2double(Nucleus);
        nucleusImgsi = nucleusImgsi - min(nucleusImgsi(:));
        nucleusImgsi = nucleusImgsi./max(nucleusImgsi(:));
        nucleusbw = (nucleusImgsi>nucleusbwThresh1);
        nucleusbw2 = imfill(nucleusbw, 'holes');
        nucleusbw3 = bwareaopen(nucleusbw2,1000);
        nucleusbw4 = imopen(nucleusbw3,strel( 'disk', 30 ) );
%         figure;
%         imshow(nucleusbw4);
        
        nucleusbw = (nucleusImgsi>nucleusbwThresh2);%& (nucleusImgsi<nucleusbwThresh1);
        nucleusbw2 = imfill(nucleusbw, 'holes');
        nucleusbw3 = nucleusbw2-bwareaopen(nucleusbw2,50000);
        nucleusbw5=bwareaopen(nucleusbw3,1000)+nucleusbw4;
        nucleusbw5 = imopen(nucleusbw5,strel( 'disk', 2 ) );
        nucleus_perim = bwperim(nucleusbw5);
%         
        figure;
        imshow(nucleusbw5);
        
        cytoplasmImgsi=im2double(cyto);
        cytoplasmImgsi = cytoplasmImgsi - min(cytoplasmImgsi(:));
        cytoplasmImgsi = cytoplasmImgsi./max(cytoplasmImgsi(:));
        % change cytoplasmImgsi to binary images
        cytoplasmbw = (cytoplasmImgsi >cytoplasmbwThresh);
        cytoplasmbw1 = imfill(cytoplasmbw,'holes');
        % get rid of the isolated dots and smooth the edge
        cytoplasmbw2 = bwareaopen(cytoplasmbw1, 2000);
        % erode the edges and get the eroded image dilated
        cytoplasmbw3 = imdilate(cytoplasmbw2,ones(30,30));
        cytoplasmbw3 = imfill(cytoplasmbw3, 'holes');
        % figure; imshowpair(cytoplasmImgsi, cytoplasmbw3, 'montage');
        cytoplasmbw3_perim = bwperim(cytoplasmbw3);
        
%         figure;
%         imshow(cytoplasmbw3);
       
        % find centroid for nucleus
        [labeledImage, numberOfBlobs] = bwlabel(nucleus_perim ~= 0);
        measurements = regionprops(labeledImage, 'Centroid');
        allCentroids = [measurements.Centroid];
        xCentroids = allCentroids(1:2:end);
        yCentroids = allCentroids(2:2:end);
        imageprop=regionprops(labeledImage, 'Image');
        indexID=regionprops(labeledImage, 'PixelList');

        mask_em = nucleusbw5;
        mask_em = imclose(mask_em, ones(5,5));
        mask_em = imfill(mask_em, 'holes');
        mask_em = bwareaopen(mask_em, 40);

        cytoplasmImgsi_c = imcomplement(imadjust(cytoplasmImgsi));
        I_mod = imimposemin(cytoplasmImgsi_c, ~cytoplasmbw3 | mask_em);
        L = watershed(I_mod);
        cellNum = max(L(:))+1;
        cellCandidates = cell(cellNum,1);
        for celli = 0:cellNum-1
            [cellCandidates{celli+1}(:,1), cellCandidates{celli+1}(:,2)] = ind2sub(size(L), find(L==celli));
        end

        nucleusxCentroids = floor(xCentroids); % convert from float to pixels
        nucleusyCentroids = floor(yCentroids);
        % Count the nucleus number within FOVi
        nucleusCount = length(nucleusyCentroids);
        cellIdx = [];
        nucleusIdx = [];
        nucleusL = bwlabel(nucleusbw5);
        for celli = 1:cellNum
            % -------------------------------------------------------------------------
            % Keep only identified cell regions with nucleus staining of FOVi in the
            % middle
            % -------------------------------------------------------------------------
            for j = 1:length(nucleusxCentroids)
                distMat =  sqrt((cellCandidates{celli}(:,1) - nucleusyCentroids(j)).^2 + (cellCandidates{celli}(:,2) - nucleusxCentroids(j)).^2);                
                if sum(distMat < 5)
                    cellIdx(end+1) = celli;
                    nucleusIdx(end+1) = j;
                    break;
                end
            end
        end
        cellInfo = cellCandidates(cellIdx);
        nucleusColCentroid = nucleusxCentroids(nucleusIdx);
        nucleusRowCentroid = nucleusyCentroids(nucleusIdx);
        % reconstruct color map for L to represent cell with >0, and empty regions
        % as 0
        cellMap = zeros(size(L));
        for celli = 1:length(cellInfo)
            cellMap(sub2ind(size(cellMap),cellInfo{celli}(:,1), cellInfo{celli}(:,2))) = celli;
        end
        %     runningTime = toc;
        %     fprintf('Time for selecting out potential cells in the FOVi: %f\n s.', runningTime); 
        % -------------------------------------------------------------------------
        % Plot and keep records of segregated cells
        % -------------------------------------------------------------------------
        figHandle = figure(...
            'Name', ['spots_FOV' num2str(fovID)], ...
            'visible', 'on');
        imshow(label2rgb(cellMap));
        hold on;
        cellCount=0;
        cellBoundaryStruct = struct(...
            'cytoplasmBoundary',{},...
            'nucleusBoundary',{},...
            'cellID',[]);

        for celli = 1:length(cellInfo)
            % give cell ID
            cellCount = cellCount + 1;
            cellBoundaryStruct(end+1).cellID = cellCount;
            % cytoplasm boundary
            boundaryPosi =  bwboundaries(cellMap == celli);
            boundaryPosi = boundaryPosi{1};
            cellBoundaryStruct(end).cytoplasmBoundary = boundaryPosi;
            % nucleus boundary
            boundaryPosi = bwboundaries(nucleusL == nucleusIdx(celli));%nucleusL(nucleusRowCentroid(celli), nucleusColCentroid(celli)));
            boundaryPosi = boundaryPosi{1};
            cellBoundaryStruct(end).nucleusBoundary =boundaryPosi;
        end
%         figure;
%         scatter(nucleusyCentroids,nucleusxCentroids);
%         for i=1:length(nucleusxCentroids)
%             text(nucleusyCentroids(i),nucleusxCentroids(i),num2str(i));
%         end
%% load 561 spots and decode information
        disp('load spots and decode information');        
        temp=load(strcat(decodePath,['FOV_' num2str(fovID+1),'_decode.mat']));
        B=[temp.SpotCoordinates.xc,temp.SpotCoordinates.yc];
        decode=temp.decode;
        shift=temp.shift;
      %% find the cells the nucleus is not close to edge
       ValidCell=ones(1,length(cellBoundaryStruct),'int16');
       for i=1:length(cellBoundaryStruct)
            temp=cellBoundaryStruct(i).nucleusBoundary;
            x=temp(:,1);
            y=temp(:,2);
            mask=poly2mask(cellBoundaryStruct(i).nucleusBoundary(:,2),cellBoundaryStruct(i).nucleusBoundary(:,1),2048,2048);
            b=sum(x==1)+sum(x==2048)+sum(y==1)+sum(y==2048);
            if b>20 && sum(mask(:))>1000 && sum(mask(:))<60000 
                ValidCell(i)=0;
            end
       end  

        %% Quantify phenotype in the nucleus
        disp('Start quantifying phenotype');
        Cellfeature = {};
        for i=[1 2 7 9]
            if i==9
                phenotype=ReadDax([RawPath,'data\Epi-750s1-650s1-560s1-488s2-405s1_',num2str(fovID,'%03d'),'_' num2str(i-1,'%02d'),'.dax'], 'startFrame', 6, 'endFrame', 6);
                SPoriginal=int64(double(phenotype)./SONratio);   
            else    
              phenotype=ReadDax([RawPath,'data\Epi-750s1-650s1-560s1-488s2-405s1_',num2str(fovID,'%03d'),'_' num2str(i-1,'%02d'),'.dax'], 'startFrame', 5, 'endFrame', 5); % load phenotype
              SPoriginal=int64(double(phenotype)./cytoratio);   
            end
            
            %SPoriginal=phenotype;
            figure;
            imshow(int16(double(phenotype)./cytoratio));
            if i==7 | i==9
            [feature,result]=SpeckleFinder(RawPath,i,fovID,cytoratio,shift,SPoriginal,cellBoundaryStruct,parameter{i,1},parameter{i,2},parameter{i,3},parameter{i,4},parameter{i,5},parameter{i,6},1);            
            else 
             [feature,result]=SpeckleFinder(RawPath,i,fovID,cytoratio,shift,SPoriginal,cellBoundaryStruct,parameter{i,1},parameter{i,2},parameter{i,3},parameter{i,4},parameter{i,5},parameter{i,6},0);
            end
            figure;
            imshow(result);
            Cellfeature{i}=feature;
        end
        
        %% map spots to cells
        disp('Start mapping the spots to cells');       
        ratio=zeros(bitnumber,length(cellBoundaryStruct));
        spotnum=zeros(1,length(cellBoundaryStruct));
        if ~isempty(B)
        for l=1:length(cellBoundaryStruct)    
               %assign spots
               
                index=(inpolygon(B(:,1), B(:,2),cellBoundaryStruct(l).cytoplasmBoundary(:,2),cellBoundaryStruct(l).cytoplasmBoundary(:,1))==1)& (inpolygon(B(:,1), B(:,2),cellBoundaryStruct(l).nucleusBoundary(:,2),cellBoundaryStruct(l).nucleusBoundary(:,1))~=1) & (sum(decode,1)~=-bitnumber)';
                
                index1= decode(:,index);
                
                num=sum(index1,2);
                spotnum(l)=spotnum(l)+sum(index);
                
                if spotnum(l)>0
                    ratio(:,l)=double(num(:,1))/double(spotnum(l));
                end
                
                %plot spots
                
                B1=B(index,:);
                if ~isempty(B1)
                    index2=sum(index1,1)>0;
                    scatter(B1(index2,1),B1(index2,2),'k',marker(l),'filled');
                    hold on;
                    index3=sum(index1,1)<0;
                    scatter(B1(index3,1),B1(index3,2),'k',marker(l)');
                    hold on;
                end
                
        end        
        end     
        %axis off;
        SaveFigure(figHandle, 'overwrite', true, ...
            'formats', {'fig', 'png'}, ...
            'savePath', mapPath);
        close(figHandle);
        disp('Saving data');
        %parsave(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']),ratio,spotnum,cellBoundaryStruct);
    parsave_SP(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']),ratio,spotnum,cellBoundaryStruct,ValidCell,Cellfeature);
     %end
end 
%close all;
%% for quantify the phenotype only:
mapPath=SetFigureSavePath([analysisSavePath,'\map_info_test\'], ...
     'makeDir', true);
parameterPath=[RawPath 'parameter1.xlsx'];
[~,~,all]=xlsread(parameterPath);
parameter=all(2:end,2:end); 
for fovID =212%:length(mlistfile{1})-2
     disp('Start quantifying phenotype');
     temp=load(['\\neptune\analysis\lncRNA\analysis\PS166_U2OS_7_broad\FOV\map_info_optimized\FOV_' num2str(fovID),'data.mat']);
     ratio=temp.ratio;
     spotnum=temp.spotnum;
     cellBoundaryStruct=temp.cellBoundaryStruct;
     ValidCell=temp.ValidCell;
     temp=load(strcat(decodePath,['FOV_' num2str(fovID+1),'_decode.mat']));
     shift=temp.shift;
        Cellfeature = {};
        for i=1:9
            if i==9
                phenotype=ReadDax([RawPath,'data\Epi-750s1-650s1-560s1-488s2-405s1_',num2str(fovID,'%03d'),'_' num2str(i-1,'%02d'),'.dax'], 'startFrame', 6, 'endFrame', 6);
                SPoriginal=int64(double(phenotype)./SONratio);   
            else    
              phenotype=ReadDax([RawPath,'data\Epi-750s1-650s1-560s1-488s2-405s1_',num2str(fovID,'%03d'),'_' num2str(i-1,'%02d'),'.dax'], 'startFrame', 5, 'endFrame', 5); % load phenotype
              SPoriginal=int64(double(phenotype)./cytoratio);   
            end
            
            %SPoriginal=phenotype;
            figure;
            imshow(int16(double(phenotype)./cytoratio));
            if i==7 | i==9
            [feature,result]=SpeckleFinder(RawPath,i,fovID,cytoratio,shift,SPoriginal,cellBoundaryStruct,parameter{i,1},parameter{i,2},parameter{i,3},parameter{i,4},parameter{i,5},parameter{i,6},1);            
            else 
             [feature,result]=SpeckleFinder(RawPath,i,fovID,cytoratio,shift,SPoriginal,cellBoundaryStruct,parameter{i,1},parameter{i,2},parameter{i,3},parameter{i,4},parameter{i,5},parameter{i,6},0);
            end
            figure; 
            imshow(result);
            Cellfeature{i}=feature;
        end
       parsave_SP(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']),ratio,spotnum,cellBoundaryStruct,ValidCell,Cellfeature);
%        boundary=cellBoundaryStruct; 
%        sample=zeros(2048,2048);
%         for j=1:length(boundary)
%             mask=poly2mask(boundary(j).nucleusBoundary(:,2),boundary(j).nucleusBoundary(:,1),2048,2048);
%             sample=sample+mask.*feature(j,6);
%         end
%         sample=sample./max(sample(:));
%         imagesc(sample,[0,1]);
%     %colorbar;
%     colorbar('southoutside');
end

%% load data
f1=[];
f2=[];
f3=[];
f4=[];
f5=[];
f6=[];
f7=[];
f8=[];
f9=[];

total=0;
ratio=[];
fov=[];
cellID=[];
valid=[];
numthreshold=30; %count the cells with spots number > threshold
for fovID=1:length(mlistfile{1})-1
    if exist(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']), 'file') == 2
        temp=load(strcat(mapPath,['FOV_' num2str(fovID),'data.mat']));
        total=total+length(temp.spotnum);
        id=temp.spotnum>numthreshold;
        ratio=[ratio temp.ratio(:,id)];
        cellID=[cellID find(id)];
        fov=[fov fovID*ones(1,length(find(id)))];
        f1=[f1 temp.Cellfeature{1}(id,1:9)'];
        f2=[f2 temp.Cellfeature{2}(id,1:9)'];
        f3=[f3 temp.Cellfeature{3}(id,1:9)'];
        f4=[f4 temp.Cellfeature{4}(id,1:9)'];
        f5=[f5 temp.Cellfeature{5}(id,1:9)'];
        f6=[f6 temp.Cellfeature{6}(id,1:9)'];
        f7=[f7 temp.Cellfeature{7}(id,:)'];
        f8=[f8 temp.Cellfeature{8}(id,1:9)'];
        f9=[f9 temp.Cellfeature{9}(id,:)'];
        valid=[valid temp.ValidCell(id)];
    end
end

% allfeature=[f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11]';
% variance=std(allfeature);
% allfeature=(allfeature-mean(allfeature))./variance;
% [coeff,score,latent,tsquared,explained]=pca(allfeature)
% figure;
% scatter3(score(:,1),score(:,2),score(:,3),'b.')
% axis equal
% xlabel('1st Principal Component')
% ylabel('2nd Principal Component')
% zlabel('3rd Principal Component')
% % for i=1:length(score(:,1))
% %  text(score(i,1),score(i,2),score(i,3),num2str(i));
% % end
% s1=score(:,1);
% s2=score(:,2);
% s3=score(:,3);
% s4=score(:,4);
Threshold=[0.1,0.1,0.1,0.1, 0.1,0.1,0.1,0.1,0.1, 0.1,0.1];
disp(['Have found ' num2str(length(ratio)) ' cells']);
% decode cells bit by bit
barcodePath=SetFigureSavePath([dataPath,'\barcode\'], ...
     'makeDir', true);
cut=0.04;
%crop the cells
for i=1:12
%     cut=0.1
%     if i==7 
%         cut=0.01;
%     end
    bits=ratio(i*3-2:i*3,:)';
    %throw away the cells on the diaganal and cell ratio all <0.1 
    idx=(bits(:,1)<cut)&(bits(:,2)<cut)&(bits(:,3)<cut);
    disp(length(find(idx)));
    idx=~idx;
    ratio=ratio(:,idx);
    cellID=cellID(idx);
    fov=fov(:,idx);
    f1=f1(:,idx);
    f2=f2(:,idx);
    f3=f3(:,idx);
    f4=f4(:,idx);
    f5=f5(:,idx);
    f6=f6(:,idx);
    f7=f7(:,idx);
    f8=f8(:,idx);
    f9=f9(:,idx);
%     
    valid=valid(idx);
%     s1=s1(idx);
%     s2=s2(idx);
%     s3=s3(idx);
%     s4=s4(idx);
%     disp('try');
%     bits=ratio(i*3-2:i*3,:)';
%     idx=(bits(:,1)>0.2& bits(:,2)>0.2)|(bits(:,3)>0.2& bits(:,2)>0.2)|(bits(:,1)>0.2& bits(:,3)>0.2);
%     disp(length(find(idx)));
%     idx=~idx;
%     ratio=ratio(:,idx);
%    cellID=cellID(idx);
%     fov=fov(idx);

end
disp(['After cropping, left ' num2str(length(ratio)) ' cells']);
decode=zeros(36,length(ratio));
for i=1:12
    %cluster three bits
    bits=ratio(i*3-2:i*3,:)';
    idx=kmeans(bits,3,'Distance','cosine');
    %find cluster number correlates to bit number
    fig = figure(...
                'Name', ['Clustering of bits',num2str(i*3-2),'-',num2str(i*3)], ...
                'visible', 'on');
    
    scatter3(bits(:,1),bits(:,2),bits(:,3),3,idx);
    xlabel(['Bit ',num2str(i*3-2)]);
    ylabel(['Bit ',num2str(i*3-1)]);
    zlabel(['Bit ',num2str(i*3)]);
    title([['Clustering of bits',num2str(i*3-2),'-',num2str(i*3)]]);
    SaveFigure(fig, 'overwrite', true, ...
                'formats', {'fig', 'png'}, ...
                'savePath', barcodePath);
    % start decoding
    temp=zeros(3,3);
    for j=1:3
        temp(j,:)=mean(bits(idx==j,:));
    end
    [~,match]=max(temp);
    for j=1:length(ratio)
        decode(i*3-3+find(match==idx(j)),j)=1;
    end
end
close all;
disp(['After cropping, left ' num2str(length(ratio)) ' cells']);

%% match barcode to the MiSeq
codebookPath='\\neptune\analysis\lncRNA\20180514_MiSeq_Broad\SG_codebook.xlsx';
[~,~,codebook]=xlsread(codebookPath);
Genename={codebook{:,1}};
barcode=zeros(length(Genename),36);
% convert barcode to binary barcode
err=zeros(1,length(codebook));
for i=1:length(Genename)
    temp=str2num(codebook{i,2});  
    for j=1:length(temp)
        barcode(i,temp(j))=1;
    end
    err(i)=12-length(temp);
end

barcodeT=barcode(:,1:36) ;
decodeT=decode(1:36,:);
Biterror=zeros(36,1);
zerotoone=zeros(36,1);
onetozero=zeros(36,1);
% match decoded barcode to the codebook
decodenum=zeros(1,length(codebook));
decodenumall=ones(1,length(codebook))*0.001;
errornum=zeros(1,5);
cellName={};
cellexp=zeros(1,length(ratio),'int64');
cellbarcode=zeros(1,length(ratio),'int64');
cellerror=zeros(1,length(ratio),'int64');
for i=1:length(ratio)
    num=-err;
    for j=1:length(codebook)
        num(j)=num(j)+sum(abs(barcodeT(j,:)-decodeT(:,i)'));
    end
    [error,barcodeid]=min(num);   
    if error/2<1&& length(find(num==error))==1
        decodenum(barcodeid)=decodenum(barcodeid)+1;
        cellName{i}=codebook{barcodeid,1};
        cellbarcode(i)=barcodeid;
        cellerror(i)=error/2;
    end
    %if error/2>=1 %
        idx=find(barcodeT(barcodeid,:)<decodeT(:,i)');
        zerotoone(idx)=zerotoone(idx)+1;
        idx=find(barcodeT(barcodeid,:)>decodeT(:,i)');
        onetozero(idx)=onetozero(idx)+1;
        decodenumall(barcodeid)=decodenumall(barcodeid)+1;
        errornum(int16(error/2+1))=errornum(int16(error/2+1))+1;
   % end 
    
end
close all;
fig = figure(...
                'Name', ['Number of cell matched to each barcode (with less than one mismatch)'], ...
                'visible', 'on');
bar(decodenum);
xlabel('Barcode');
ylabel('Cell number');
title(['Number of cell matched to each barcode']);
SaveFigure(fig, 'overwrite', true, ...
                'formats', {'fig', 'png'}, ...
                'savePath', barcodePath);    
fig = figure(...
                'Name', ['Error rate'], ...
                'visible', 'on');
bar(errornum);
xlabel('Error number');
ylabel('Cell number');
set(gca,'XTickLabel',0:1:length(errornum)-1);
title(['Exact ratio: ',num2str((errornum(1))/sum(errornum))]);
SaveFigure(fig, 'overwrite', true, ...
                'formats', {'fig', 'png'}, ...
                'savePath', barcodePath);    
figure;
bar(decodenum./decodenumall);
fig = figure(...
    'Name', ['zero to one error'], ...
    'visible', 'on');
bar(zerotoone);
SaveFigure(fig, 'overwrite', true, ...
                'formats', {'fig', 'png'}, ...
                'savePath', barcodePath);
fig = figure(...
    'Name', ['one to zero error'], ...
    'visible', 'on');
bar(onetozero);
SaveFigure(fig, 'overwrite', true, ...
                'formats', {'fig', 'png'}, ...
                'savePath', barcodePath);
save([barcodePath 'decoded_data.mat'],'cellName','fov','cellID','cellerror','cellbarcode','decodeT','cytoratio');            
%% count the gene distribution
idx=~cellfun(@isempty, cellName);
allname=cellName(idx);
valid1=valid(idx);
f11=f1(:,idx);
f21=f2(:,idx);
f31=f3(:,idx);
f41=f4(:,idx);

% idx2=~isnan(f91(8,:));
% f92=f91(:,idx2);
% valid2=valid1(idx2);
% allname1=allname(idx2);

valid2=valid1;
clear unique;
[unique,ic,ib]=unique(allname);
h2=histogram(ib);

%% plot the phenotype of MALAT1
PhenotypePath=SetFigureSavePath([dataPath,'\Phenotype_optimized\'], ...
     'makeDir', true);
Gene='HNRNPM_';
feature=7;
phenotype=2;
% allpheno={f92};
allpheno={f11,f21,f31,f41};
% allpheno={f11,f21,f71,f91};
%allblankmean={};
for pheno=1:4
    f12=allpheno{pheno};
    temp=size(f12);
    if pheno==7 || pheno==9
        temp=33;
    end
    meanPhenotype=zeros(temp(1),length(unique));
    SEMPhenotype=zeros(temp(1),length(unique));
    totalCell=zeros(1,length(unique));
    ValidCell=zeros(1,length(unique));
    idx=(ib>=115)' & (ib<=154)' & (valid2==1);
    null=f12(1:temp(1),idx);
    blankmean=mean(null');
    allblankmean{pheno}=blankmean;
    pvalue=zeros(temp(1),length(unique));
    for i=1:length(unique)
        idx=(ib==i)' & (valid2==1);
        for j=1:temp(1)
            meanPhenotype(j,i)=mean(f12(j,idx));
            SEMPhenotype(j,i)=std(f12(j,idx))/sqrt(length(f12(j,idx)));
             [~,pvalue(j,i)]=ttest2(f12(j,idx),null(j,:));
            %pvalue(j,i) = ranksum(f12(j,idx),null(j,:));
        end    
    end

    for i=1:temp(1)
        fig = figure(...
                        'Name', ['Phenotype ' num2str(pheno) ' feature ' num2str(i) ' distribution'], ...
                        'visible', 'on');
        X=1:length(unique);
        errorbar(X,meanPhenotype(i,:),SEMPhenotype(i,:),'o');
        set(gca,'XTick',X)
        set(gca,'XTickLabel',unique)
        set(gca,'XTickLabelRotation',45)
        SaveFigure(fig, 'overwrite', true, ...
                        'formats', {'fig', 'png'}, ...
                        'savePath', PhenotypePath); 
        fig = figure(...
                        'Name', ['volcano plot of Phenotype ' num2str(pheno) ' feature ' num2str(i) ' '], ...
                        'visible', 'on');
        scatter(log10(meanPhenotype(i,:)/blankmean(i)),-log10(pvalue(i,:)),'.');
        hold on;
        scatter(log10(meanPhenotype(i,161:165)/blankmean(i)),-log10(pvalue(i,161:165)),'bo');
        hold on;
        scatter(log10(meanPhenotype(i,66:68)/blankmean(i)),-log10(pvalue(i,66:68)),'ro');
        SaveFigure(fig, 'overwrite', true, ...
                        'formats', {'fig', 'png'}, ...
                        'savePath', PhenotypePath); 
    end
    close all;
end

save([barcodePath 'phenotype_data_optimized2.mat'],'allpheno','f12','f1','ib','allname','valid2','unique','pvalue','meanPhenotype','allblankmean');  
%% PCA
variance=std(f12(3:end,:)');
f13=(f12(3:end,:)'-mean(f12(3:end,:)'))./variance;
[coeff,score,latent,tsquared,explained]=pca(f13);
figure;
scatter3(score(:,1),score(:,2),score(:,3),'b.')
axis equal
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')

meanComponent=zeros(3,length(unique));
SEMComponent=zeros(3,length(unique));
idx=(ib>=163)' & (valid2==1);
null=score(idx,1:3);
blankmean=mean(null);
pvalue=zeros(3,length(unique));
for i=1:length(unique)
    idx=(ib==i) & (valid2==1)';
%     if i==100
%         idx=ib==122;
%     end
%     if i>121
%         idx=ib==i+1;
%     end
    for j=1:3
        meanComponent(j,i)=mean(score(idx,j));
        SEMComponent(j,i)=std(score(idx,j))/sqrt(length(score(idx,j)));
        [~,pvalue(j,i)]=ttest2(score(idx,j),null(:,j));
    end    
end

for i=1:3
    fig = figure(...
                    'Name', ['MALAT1 component ' num2str(i) ' distribution'], ...
                    'visible', 'on');
    X=1:length(unique);
    errorbar(X,meanComponent(i,:),SEMComponent(i,:),'o');
    set(gca,'XTick',X)
    set(gca,'XTickLabel',unique)
    set(gca,'XTickLabelRotation',45)
    SaveFigure(fig, 'overwrite', true, ...
                    'formats', {'fig', 'png'}, ...
                    'savePath', barcodePath); 
    fig = figure(...
                    'Name', ['volcano plot of MALAT1 component ' num2str(i) ' '], ...
                    'visible', 'on');
    scatter(meanComponent(i,:)/blankmean(i),-log(pvalue(i,:)),'o');
    hold on;
    scatter(meanComponent(i,163:167)/blankmean(i),-log(pvalue(i,163:167)),'o');
    SaveFigure(fig, 'overwrite', true, ...
                    'formats', {'fig', 'png'}, ...
                    'savePath', barcodePath); 
end

    
    
%%
% for i=1:length(score(:,1))
%  text(score(i,1),score(i,2),score(i,3),num2str(i));
% end
s1=score(:,1);
s2=score(:,2);
s3=score(:,3);
s4=score(:,4);


% fig = figure(...
%                 'Name', ['MALAT1 speckle intensity distribution'], ...
%                 'visible', 'on');
% X=1:length(unique);
% errorbar(X,meanSPint,SEMSPint,'o');
% set(gca,'XTick',X)
% set(gca,'XTickLabel',unique)
% set(gca,'XTickLabelRotation',45)
% SaveFigure(fig, 'overwrite', true, ...
%                 'formats', {'fig', 'png'}, ...
%                 'savePath', barcodePath);              
% fig = figure(...
%                 'Name', ['MALAT1 nucleus median Intensity distribution'], ...
%                 'visible', 'on');
% X=1:length(unique);
% errorbar(X,meanmedianIN,SEMmedianIN,'o');
% set(gca,'XTick',X)
% set(gca,'XTickLabel',unique)
% set(gca,'XTickLabelRotation',45)
% SaveFigure(fig, 'overwrite', true, ...
%                 'formats', {'fig', 'png'}, ...
%                 'savePath', barcodePath);  
% fig = figure(...
%                 'Name', ['MALAT1 speckle area Intensity distribution'], ...
%                 'visible', 'on');
% X=1:length(unique);
% bar(temp);
% set(gca,'XTick',X)
% set(gca,'XTickLabel',unique)
% set(gca,'XTickLabelRotation',45)
% SaveFigure(fig, 'overwrite', true, ...
%                 'formats', {'fig', 'png'}, ...
%                 'savePath', barcodePath);   
%% plot certain gene phenotype
phenotype=7;
codebookPath='\\neptune\analysis\lncRNA\20180514_MiSeq_Broad\SG_codebook.xlsx';
[~,~,codebook]=xlsread(codebookPath);
TargetGene='HNRNPK_2';
mosaic=[];
allmosaic=[];
num=0;
allfeature1=[];
allfeature7=[];
allfeature8=[];
allfeature6=[];
allfeature9=[];
for i=1:length(cellName) 
    if strcmp(cellName{i},TargetGene)
       image=ReadDax([RawPath,'data\Epi-750s1-650s1-560s1-488s2-405s1_',num2str(fov(i),'%03d'),'_',num2str(phenotype,'%02d'),'.dax'], 'startFrame', 5, 'endFrame', 5);
       %image1=int16(double(image)/cytoratio);
       image1=image;
        disp(['Index: ' num2str(i), ',FOV: ', num2str(fov(i)), ',ID: ' num2str(cellID(i)),',Error: ' num2str(cellerror(i))]);
        disp(['Barcode ID: ', codebook{cellbarcode(i),2},',Decoded Barcode: ',num2str(find(decodeT(:,i)'))]);
         
        temp=load(strcat(mapPath,['FOV_' num2str(fov(i)),'data.mat']));
        
        coordi=temp.cellBoundaryStruct(cellID(i)).nucleusBoundary;
       disp('hey');
        X=mean(coordi(:,1));
        Y=mean(coordi(:,2));
        allfeature1=[allfeature1 f1(1,i)];
        allfeature7=[allfeature7 f1(7,i)];
        allfeature8=[allfeature8 f1(8,i)];
        allfeature9=[allfeature9 f1(9,i)];
        allfeature6=[allfeature6 f1(6,i)];
        num=num+1;
        if num>=100
            break;
        end
        if Y-149<0 
            Y=150;
        end
        if X-149<0
            X=150;
        end     
        if Y+150>2048
            Y=2048-150;
        end
        if X+150>2048
            X=2048-150;
        end
        mosaic=[mosaic,image(X-149:X+150,Y-149:Y+150)];
        if mod(num,7)==0
            allmosaic=[allmosaic;mosaic];
            mosaic=[];
        end
            
    end
end
disp(num);
figure;
imshow(imadjust(allmosaic));
figure;
imshow(imadjust(mosaic));
title(['Find ',num2str(num),' cells']); 
disp(allfeature1);
disp(allfeature6);
disp(allfeature7);
disp(allfeature8);
disp(allfeature9);
disp(mean(allfeature1));
disp(mean(allfeature6));
disp(mean(allfeature7));
disp(mean(allfeature8));
disp(mean(allfeature9));
%%
bit=8;
cell=[];
for i=1:length(ratio)
    num=-err;
    for j=1:length(codebook)
        num(j)=num(j)+sum(abs(barcodeT(j,:)-decodeT(:,i)'));
    end
    [error,barcodeid]=min(num);   
    if error/2>0 && length(find(num==error))==1
        idx=find(barcodeT(barcodeid,:)<decodeT(:,i)');
        if ~isempty(find(idx==bit))
            cell=[cell ratio(bit-1:bit+1,i)];
        end
%         idx=find(barcodeT(barcodeid,:)>decodeT(:,i)');
%         if ~isempty(find(idx==bit))
%             cell=[cell ratio(bit-1:bit+1,i)];
%         end
    end
end
figure;
scatter3(cell(1,:),cell(2,:),cell(3,:));