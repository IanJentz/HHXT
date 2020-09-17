function obj = readFiberFile(obj,file,subsections)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 2; subsections = []; end

%% read the text file
 fid = fopen(file, 'r') ;              % Open source file.
 for i = 1:4
 fgetl(fid) ;                                  % Read/discard line.
 end
 buffer = fread(fid, Inf) ;                    % Read rest of the file.
 fclose(fid) ;
 fid = fopen('workingFiber.txt', 'w')  ;   % Open destination file.
 fwrite(fid, buffer) ;                         % Save to file.
 fclose(fid) ;
 
 fid = fopen('workingFiber.txt', 'r')  ;   % Open destination file.
 pline =  fgetl(fid) ;
 fclose(fid) ;
 pos = str2double(strsplit(pline,'\t'));
 pos(1) = [];
 npos = length(pos);
 
  tdelim = tdfread('workingFiber.txt','\t');
  
%   load('tdelim.mat','tdelim');


 %% handle times
 fnames = fieldnames(tdelim);
 timestrings = string(tdelim.(fnames{1}));
 ntimes = length(timestrings);
 
 datetimes = cell(ntimes,1);
 for i = 1:ntimes
 split1 = strsplit(timestrings(i),'T');
 split2 = strsplit(split1(2),'.');
 split = [strsplit(split1(1),'-'),strsplit(split2(1),':'),split2(2)];
 YMDHMISMS = zeros(1,7);
 for j = 1:7
     YMDHMISMS(j) = str2double(split(j));
 end
 YMDHMISMS(7) = YMDHMISMS(7)*10^(2-floor(log10(YMDHMISMS(7))));
 datetimes{i} = datetime(YMDHMISMS(1),YMDHMISMS(2),YMDHMISMS(3),YMDHMISMS(4),YMDHMISMS(5),YMDHMISMS(6),YMDHMISMS(7));
 end
 tstart = datetimes{1};
 elapsedTime = cell(ntimes,1);
 elapsedSeconds = zeros(ntimes,1);
 for i = 2:ntimes
    elapsedTime{i} = (datetimes{i}-tstart);
    elapsedSeconds(i) = seconds(elapsedTime{i});
 end
 
 %% get freqency shift (GHz)
 FreqShift = zeros(ntimes,npos);
 for i = 2:length(fnames)
     FreqShift(:,i-1) = tdelim.(fnames{i});
 end
 
 %%apply filters
 if obj.FilterFiber

 load('fiberFilter.mat','lpFilt','PhaseShift','Flims');
 
 secfilt = elapsedSeconds'-PhaseShift;
 new_ind = secfilt>=0;
 timestrings(~new_ind) = [];
 elapsedTime(~new_ind) = [];
 elapsedSeconds(~new_ind) = [];
 ntimes = length(elapsedSeconds);
 Old_FreqShift = FreqShift;
 FreqShift(~new_ind,:) = [];
 for i = 1:npos
    xdata = Old_FreqShift(:,i)';
    xbound = xdata;
    xbound(xdata<Flims(1)) = nan;
    xbound(xdata>Flims(2)) = nan;
    xbound = movmean(xbound,5);
    xfill = fillmissing(xbound,'linear');
    xfilt = filter(lpFilt,xfill); 
    xfilt = xfilt(new_ind);
    FreqShift(:,i) = xfilt';
 end

%  load('fiberFilter.mat','Flims');
%  
%  for i = 1:npos
%     xdata = FreqShift(:,i)';
%     xbound = xdata;
%     xbound(xdata<Flims(1)) = nan;
%     xbound(xdata>Flims(2)) = nan;
%     xbound = movmean(xbound,5);
%     xfill = fillmissing(xbound,'linear');
%     FreqShift(:,i) = xfill';
%  end
 
 for i = 1:ntimes
    xdata = FreqShift(i,:);
    xbound = xdata;
    xbound(xdata<Flims(1)) = nan;
    xbound(xdata>Flims(2)) = nan;
    xbound = movmean(xbound,5);
    xfill = fillmissing(xbound,'linear','EndValues',0); 
    xfill = movmean(xfill,10); 
    FreqShift(i,:) = xfill;
 end
 
 end
 
 %% calculate temperature (C)
 Temperature = -1.33e-4*FreqShift.^2 - 0.748*FreqShift - 0.229 + 15;
 
 %% find subsections, if necessary
 if isempty(subsections)
     
     [~,locs] = findpeaks(mean(Temperature),pos,'MinPeakHeight',20,'MinPeakDistance',0.1);
     if mod(length(locs),2); locs(end+1) = locs(end); end
     locs = locs - (.25*25.4/1000)*( (-1).^(1:length(locs)) );
     subsections = sort(locs);   
 else
     if isa(subsections,'cell') && ~mod(length(subsections),2)
         subsect = [];
         for i = 1:length(subsections)
             subsect = [subsect,subsections{i}];
         end
         subsections = sort(subsect);
     else
         error('invalid subsections input');
     end
 end
 subsecInd = cell(length(subsections)/2,1);
 subsecPos = cell(length(subsections)/2,1);
 for i = 1:(length(subsections)/2)
     subsecPos{i} = [subsections(2*i-1),subsections(2*i)];
     subsecInd{i} = find( (pos >= subsections(2*i-1)) & (pos <= subsections(2*i)) );
 end
 
  %flip every other fiber reading
 for i = 1:length(subsecInd)
     if ~mod(i,2)
     subsecInd{i} = flip(subsecInd{i});
     end
 end
 
 %% build structure
 FiberReading = struct('npos',npos,'ntimes',ntimes,'pos',1,...
     'tstart',tstart,'timestrings',1,...
     'elapsedTime',1,'elapsedSeconds',1,...
     'subsecPos',1,'subsecInd',1,...
     'FreqShift',1,'Temperature',1);
 FiberReading.pos = pos;
 FiberReading.timestrings = timestrings;
 FiberReading.elapsedTime = elapsedTime;
 FiberReading.elapsedSeconds = elapsedSeconds;
 FiberReading.subsecPos = subsecPos;
 FiberReading.subsecInd = subsecInd;
 FiberReading.FreqShift = FreqShift;
 FiberReading.Temperature = Temperature;
 
 obj.Fiber = FiberReading;
 obj.sec_diff = seconds(obj.time_start-obj.Fiber.tstart)+( (5*60 + 5)*60 + 55);
 
end

