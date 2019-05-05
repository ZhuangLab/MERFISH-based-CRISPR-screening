% Load the reads
experiment = ''; %% experiment name
%library = 'Kpn1_S1';
library = 'BstX1_S2'; % different library
addpath(genpath('\\neptune\analysis\htseq_fun\Matlab'));
addpath('\\neptune\analysis\htseq_fun\sequence_analysis\functions');
javaclasspath('\\neptune\analysis\htseq_fun\sequence_analysis\');

dataRoot = ['']; %Path where the raw sequencing data is
saveRoot = [''];
savePath = [saveRoot library '\']; %Path where the save Path is
mkdir(saveRoot);
mkdir(savePath);
display('Start reading fastq');
read1File = [dataRoot library '_L001_R1_001.fastq'];
read2File = [dataRoot library '_L001_R2_001.fastq'];

read1 = fastqread(read1File);
read2 = fastqread(read2File);
display('Finished reading fastq');
SetFigureSavePath(savePath, 'makeDir', true);

%% Get the protospacer and barcode sequence from list

PSPath='...\one-step_order.xlsx';  % put the XLSX files of name and sequence of the PS
BarcodePath='...\All_Sequence_Used.xlsx'; % the barcode sequence files


[~,~,temp1]=xlsread(BarcodePath);
[~,~,temp2]=xlsread(PSPath,1);

finalname=temp2(:,1);
finalPS=temp2(:,2);
for i=1:length(finalPS)
   finalPS{i}=finalPS{i}(25:45); 
end
Barcode=temp1(3:38,2);

%% Extract sequence components from set 1
display('Start extract the sequence');
%% Create parallel pool
if isempty(gcp('nocreate'))
    p = parpool(25); % Set this number to control the number of parallel workers.
                     % Polite usage would suggest these maximum values: morgan = 20, cajal = 10
else
    p = gcp;
end
%%
display('Start finding the barcodes');
umi5prime = 'tctagacatcgaaaacgagatgacggacggc';
temp=cell(1,length(read1));
temp(:)={''};
Arrange=struct('Header',temp,'Sequence',temp,'UMI',temp);
Bitnumber=zeros(1,length(read1));
parfor i=1:length(read1)
    currentSequence = read1(i).Sequence;
    currentSequence2 = read2(i).Sequence;
    umiAlign = localalign(currentSequence2(1:80), umi5prime);
    if umiAlign.Start(1)>20
        Arrange(i).UMI= currentSequence2(umiAlign.Start(1)-20:umiAlign.Start(1)-1);
    else
        Arrange(i).UMI='';
    end
    FS=currentSequence;
    RS=seqrcomplement(currentSequence2);
    score=zeros(1,length(finalPS));  
    Arrange(i).Header='';
    Arrange(i).Sequence='';
    score=zeros(1,length(Barcode));
%     ali=[];
    for j=1:length(Barcode)
        [score(j),idx]= max([swalign(FS,Barcode{j}),swalign(RS,Barcode{j})]);
    end
    for j=1:length(finalPS)
            if ~isempty(findstr(FS,finalPS{j}))
                Arrange(i).Header=finalname{j};
            end       
    end
    [Score,id]=sort(score,'descend');
    bitarrange=zeros(1,12);
    bitnum=zeros(1,12,'int16');
    for j=1:length(score)
       if bitarrange(ceil(id(j)/3))==0 & Score(j)*1.5>max(Score) & max(Score)>70
        bitnum(ceil(id(j)/3))=id(j);
        bitarrange(ceil(id(j)/3))=1;
       end
    end   
    for j=1:12
        if bitnum(j)>0
           Arrange(i).Sequence=[Arrange(i).Sequence ' ' num2str(bitnum(j))];
           Bitnumber(i)=Bitnumber(i)+1;
        end
    end
    if mod(i,10000)==0
        display([num2str(i/length(read1)*100) ' percent is completed']);
    end
end

%% Get rid of reads with bit number less than 10
display('Start analyzing sequence');
% idx=(Bitnumber==10) |(Bitnumber==11) | (Bitnumber==12); for KpnI library
idx=Bitnumber>0;
 New=Arrange(idx);
%New=Arrange;  for KpnI library
%fastawrite([savePath,'\BarcodeArrangement.fasta'],New);for KpnI library
%clear the memory
%% Plot the unique barcode distribution
[uniqueBarcode,ib,id]=unique({New.Sequence});
fig = figure(...
                'Name', ['Distribution of unique Barcode' ], ...
                'visible', 'on');
 h1=[];
 for i=1:length(uniqueBarcode)
     h1(i)=sum(id==i);
 end
 plot(sort(log10(h1),'descend'));
 tempid=h1>2;
 title(['Number of unique Barcode: ' num2str(length(ib)) ' number of Barcode>2 reads: ' num2str(length(find(tempid)))]);
 SaveFigure(fig, 'overwrite', true, ...
                'formats', {'fig', 'png'}, ...
                'savePath', savePath);

%% Plot the the distribution of reads number of unique Gene
[uniqueGene,ia,ic]=unique({New.Header});    
fig = figure(...
                'Name', ['Distribution of unique Gene' ], ...
                'visible', 'on');
 %h2=histogram(ic);
 h2=[];
 for i=1:length(uniqueGene)
     h2(i)=sum(ic==i);
 end
 
 plot(sort(log10(h2),'descend'));
 xlabel('Protospacer');
 ylabel('Counts (log10)');
 tempid=h2>2;
 title(['Number of unique Gene: ' num2str(length(ia)) ' number of Gene>2 reads: ' num2str(length(find(tempid)))]);
    SaveFigure(fig, 'overwrite', true, ...
                'formats', {'fig', 'png'}, ...
                'savePath', savePath);
            
%% Plot the the distribution of reads number of UMI
[uniqueUMI,ix,iy]=unique({New.UMI});    
fig = figure(...
                'Name', ['Distribution of unique Gene' ], ...
                'visible', 'on'); 
 h2=histogram(iy);
 plot(sort(h2,'descend'),'.');
 tempid=h2>2;
 title(['Number of unique UMI: ' num2str(length(uniqueUMI)) ' number of UMI>2 reads: ' num2str(length(find(tempid)))]);
    SaveFigure(fig, 'overwrite', true, ...
                'formats', {'fig', 'png'}, ...
                'savePath', savePath);

%% Number of unique combination
% fullid=~cellfun(@isempty,{New.Header});
% New1=New(fullid);
Combination=strcat({New.Header},{' '},{New.Sequence},{' '},{New.UMI});
[uniqueCombin,ie,idx]=unique(Combination);
fig = figure(...
                'Name', ['Distribution of unique Protospacer_Barcode_UMI combination' ], ...
                'visible', 'on');
 %h5=histogram(idx);
 h5=[];
 for i=1:length(uniqueCombin)
     h5(i)=sum(idx==i);
 end
 plot(sort(log10(h5),'descend'),'.');
 tempid=h5>2;
 title(['Number of unique Combination: ' num2str(length(ie)) ' number of combination>2 reads: ' num2str(length(find(tempid)))]);
    SaveFigure(fig, 'overwrite', true, ...
                'formats', {'fig', 'png'}, ...
                'savePath', savePath);
%% show the possible codebook
[sorted,id]=sort(h5,'descend');
tempid=sorted>2;
temp=uniqueCombin(id(tempid));
codebook=[];
for i=1:length(temp)
    tempstring=strsplit(temp{i});
    codebook(i).Header=tempstring{1};
    tempstring2=tempstring;
    for j=1:length(tempstring2)
        tempstring2{j}=[tempstring2{j} ' '];
    end
    codebook(i).Sequence=horzcat(tempstring2{2:end-1});
    codebook(i).UMI=tempstring{end};
end

%%
%fastawrite([savePath '\codebook.fasta'],codebook);
Gene={codebook.Header};
Code={codebook.Sequence};
UMI={codebook.UMI};
reads=sorted(tempid);
save([savePath '\codebook.mat'],'Gene','Code','UMI','reads','New');
disp('Saved all reads');

%%
codebook1=[];
index1=1;
percent=reads./sum(reads);
check=zeros(1,length(Gene));
rank=1:length(UMI);
for i=1:length(Gene)
    if check(i)==0 && ~isempty(Code{i}) && ~isempty(UMI{i}) && ~isempty(Gene{i})
        idx=i;
        for j=i+1:length(Gene)
            if  ~isempty(UMI{j})&& ~isempty(Code{j})
                [score,align]=swalign(UMI{i},UMI{j});
                if (check(j)==0) && ContainBarcode(Code{i},Code{j}) && sum(align(2,:)=='|')>14
                    check(j)=1;
                    if length(Code{j})>length(Code{i})
                        idx=j;
                    end
                end
            end
        end
        codebook1(index1).Header=Gene(idx);
        codebook1(index1).Sequence=Code{idx};
        codebook1(index1).UMI=UMI{idx};
        codebook1(index1).percentage=percent(i);
        codebook1(index1).rank1=rank(i);
        index1=index1+1;
    end   
    disp(i);
end

save([savePath '\codebook1.mat'],'codebook1');
disp('Saved codebook1');