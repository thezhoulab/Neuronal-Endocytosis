%% This file is used to calculate average 1D autocorrelation amplitude of MPS

clear all;
close all;
ipath = '\\E1-054867\Jinyu\20250801 a2KD towards a2b2\#4 DIV14 a2GFP MAP2_Cy3 a2_AF647\Analysis\diat_regions';
pixelSize = 167;
dbin = 10/1000;  % micr
bb=0;
average_c=[];
txt_file='autocorrelation_amplitude.txt';
fileID = fopen([ipath txt_file],'w');
fprintf(fileID,'%8s\t%8s\t%26s\t%8s\t%8s\r\n', 'foldername', 'filename','autocorrelation_amplitude','diameter','length');
d = dir(ipath);
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'..'})) = [];
aamplitude=[];
for jj=1%:length(nameFolds)
    filename = dir(fullfile(strcat([ipath nameFolds{jj,1} '\'],'*.bin')));
    
    for ii=1:length(filename)
        
        [MList, memoryMap] = ReadMasterMoleculeList([ipath nameFolds{jj,1} '\' filename(ii).name]);
        % get the y_corr value to do the FFT
        xcdata = MList.xc(find(MList.c==1)); %pixel -> nm
        ycdata = MList.yc(find(MList.c==1));
        scatter(xcdata,ycdata,'.')
        x = (xcdata-min(xcdata))*pixelSize/1000;    % micron ' dual objective 141 '148
        % x= X(:,7);
        
        bin=min(x):dbin:max(x);
        L=size(bin,2);
        nx=hist(x,bin);
        nx = (nx-bb); % important to keep the zero in the middle of account
        
        figure(2);
        hist(x,bin); % do FFT from -800 to 800?
        xlabel('in micron');
        
        ind = find(bin>-5 & bin<5); % to plot and calculat which part
        Nind = size(ind,2);
        
        NFFT = 2^nextpow2(L);
        nx1 = zeros(NFFT,1);
        nx1(1:Nind) = nx(ind);
        fx = fft (nx1,NFFT)/L;
        f = 1/dbin/2*linspace(0,1,NFFT/2+1);
        
        [c,i] = max(abs(fx));
        
        
        figure(3)
        zz=2*abs(fx(1:NFFT/2+1));
        plot(f,zz)
        title(sprintf('period = % .2f nm',1/f(i)*1000));
        xlabel('in 1/micron');
        
        figure(3)
        plot(bin(ind),nx(ind))
        xlabel('in micron');
        
        tlag = 80;
        tt=nx (ind)';
        c=autocorr(tt,tlag);% time lag 65
        % cc=xcorr(tt+bb,tt+bb);
        xxx=0:(tlag+1)*10/length(c):(tlag+1)*10-(tlag+1)*10/length(c);
        plot(xxx,c,'k','LineWidth',1.25);
        x2 = 600;
        y2 = 0.9;
        txt2 = ['amplitude=' num2str(-(c(11)+c(30))/2+c(20))];
        text(x2,y2,txt2)
        
        
        %   figHandle=figure;
        % idx = find(hist{i,j}.BinEdges<900);
        % plot(hist{i,j}.BinEdges(idx), acfValues{i, j}(idx),'k','LineWidth',1.25);
        xlabel('Lag (nm)','FontSize', 25)
        ylabel('Autocorrelation','FontSize', 25);
        set(gca,'XTick',0:200:900);
        set(gca,'YTick',-0.4:0.2:1);
        xlim([0 900])
        set(gca,'FontSize',20)
        %   saveas(gcf,[ipath nameFolds{jj,1} '\' filename(ii).name(1:end-4) '_autoCorrelation.png']);
        average_c=[average_c,c];
        %   low_peak_1=find(xxx>=90&xxx<=100);
        %   low_peak_2=find(xxx>=270&xxx<=310);
        %   high_peak=find(xxx>=180&xxx<=200);
        
        ycdata = MList.yc; %pixel -> nm
        y = (ycdata-min(ycdata))*pixelSize/1000;    % micron ' dual objective 141 '148
        % x= X(:,7);
        y(find(x>=1))=[];
        
        bin=min(y):dbin:max(y);
        L=size(bin,2);
        ny=hist(y,bin);
        threshold=0.005*sum(ny);
        consecutive_index=find(ny>threshold);
        consecutive_index=[0,consecutive_index,consecutive_index(end)+2];
        diameter=(max(diff(find(diff(consecutive_index)~=1)))-1)*dbin;
        aamplitude(ii)=-(c(11)+c(30))/2+c(20);
        fprintf(fileID,'%10s\t%20s\t%4.2f\t%4.2f\t%4.2f\r\n',nameFolds{jj,1},filename(ii).name,-(c(11)+c(30))/2+c(20),diameter, max(x));
        
    end
end

figure;
h=hist(aamplitude,-0.2:0.1:1.2);
average_total=mean(average_c,2);
readout=-(average_total(11)+average_total(30))/2+average_total(20)
  std(aamplitude)/sqrt(length(aamplitude))
  length(find(aamplitude>=0.27))/length(aamplitude)
'type1'
aamplitude_random=aamplitude;%(randperm(length(aamplitude)));
replicate=[mean(aamplitude_random(1:round(length(aamplitude_random)/3))),mean(aamplitude_random(round(length(aamplitude_random)/3)+1:round(2*length(aamplitude_random)/3))),mean(aamplitude_random(round(2*length(aamplitude_random)/3)+1:end))];
ste=std(replicate)/sqrt(3)
average=mean(replicate)
replicate

'type2'
replicate=[mean(aamplitude_random(1:3:end)),mean(aamplitude_random(2:3:end)),mean(aamplitude_random(3:3:end))];
ste=std(replicate)/sqrt(3)
replicate
average=mean(replicate)

figure;
plot(xxx, mean(average_c,2),'k','LineWidth',1.25);
readout_average=mean(average_c,2);
xlabel('Lag (nm)','FontSize', 25)
ylabel('Autocorrelation','FontSize', 25);
set(gca,'XTick',0:200:900);
set(gca,'YTick',-0.4:0.2:1);
xlim([0 900])
set(gca,'FontSize',20)
saveas(gcf,[ipath 'folder_average' '.png']);
savefig(gcf,[ipath 'folder_average' '.fig']);
average=mean(average_c,2);
figure;
errorbar(xxx, mean(average_c,2),std(average_c,0,2)/sqrt(size(average_c,2)),'k','LineWidth',1.25);
xlabel('Lag (nm)','FontSize', 25)
ylabel('Autocorrelation','FontSize', 25);
set(gca,'XTick',0:200:900);
set(gca,'YTick',-0.4:0.2:1);
xlim([0 900])
set(gca,'FontSize',20)
saveas(gcf,[ipath 'folder_average_with_error_bar' '.png']);
savefig(gcf,[ipath 'folder_average_with_error_bar' '.fig']);
%%
readout_group=-(average_c(11,:)+average_c(30,:))/2+average_c(20,:);
readout_group=readout_group';
readout_stderror=std(readout_group)/sqrt(size(readout_group,2));
