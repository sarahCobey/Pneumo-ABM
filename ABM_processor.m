function [ simpsonInfo prevalenceInfo meanThetas stdThetas hfluInfo rankCorrs fitCorrs meanCocols SDCocols simpsonAgeMean simpsonAgeSD shared kidReducts] = ABM_processor( sims, extThold, treatment, treatmentVal, numNeighborhoods, modelName, Hflu_simulated )

% ABM Processor
% Pneumo-ABM (c) 2012, S. Cobey


% Plotting options
cohortNames = {'<1 y';'1-2';'3-4';'5-20';'>20' };
durationNames = {'<1';'1';'2';'3';'4';'5'};
cohorts = [ 0 0; 1 2; 3 4; 5 20; 21 110 ];
if ( Hflu_simulated == 1 ) 
    serotypes = [ 0:25 ];
    pneumoSerotypes = [ 1:(length(serotypes)-1) ];
    serotypeNames = {'Strain 0';'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25'};
    numStrains = length(serotypes);
    Hflu_index = numStrains;
    numPneumoSerotypes = numStrains-1;
else
    serotypes = [ 0:24 ]; % Numbering in input files
    numStrains = length(serotypes);
	serotypeNames = {'Serotype 0';'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24'};
    numPneumoSerotypes = length(serotypes);
end
ages = [0:110]; 
timeInt = 20; % interval between time points on x-axes (used for epid results)
demInt = 100;
sampleDuration = 20; % duration (in years) of period in which to obtain sampling statistics; end of epid sampling
coinfectionBuffer = 10; % max number of coinfections to consider
numSims = length(sims);
minToPlot = min(numPneumoSerotypes,10);
yearsPostVaccine = 50;
postVaccineStart = 101; % epid times index
preVaccineYear = 99; % epid times index

if (length(serotypeNames)~=length(serotypes))
    warning('Number of serotype names should match number of serotypes')
end

% Set up cohort indices
cohortIndices = cohorts+1;
nDimension = 1; % neighborhood dimension (used for summing)
rowSP = 2 + min(minToPlot,length(serotypes));
if ( Hflu_simulated == 1 )
   rowSP = rowSP + 1; 
end
colSP = 3;
numSP = rowSP*colSP;

% For durations by age, establish lower bounds for age groups
durationCohorts = [ 0 1 2 3 4 5 ];

% Set up summary statistics for simulations
extinctions = zeros(numPneumoSerotypes+1,1); % holds number of simulations with each number of extinctions
finalTS = zeros(numStrains,sampleDuration,numSims); % holds final time series for each simulation
finalPrevs = zeros(numStrains,sampleDuration,numSims); % holds prevalences
totalCoinfections = zeros(sampleDuration,numSims);
finalPrevalences = zeros(numSims,1);
thetasByAge = zeros(length(durationCohorts),numSims);
clear effAgeDist
nrow = size(cohortIndices,1);
percentHflu = zeros(numPneumoSerotypes,100,numSims); % not robust (2nd index)
finalPercentHflu = zeros(numPneumoSerotypes,sampleDuration,numSims); % not robust (2nd index)
pneumoCocols = zeros(10,numSims); % assumes <10 co-colonizations, 2nd index is mean, SD
simpsonAges = zeros(length(cohortNames),numSims);
reductions = zeros(numSims,1);
prePostPrevs = zeros(numPneumoSerotypes,numSims);
carriageDiffsByYear = zeros(150,numSims); % should be length(epidTimes)
carriageDiffsByYearOver5 = zeros(150,numSims);
simpsonTSUnder5AllSims = zeros(150,numSims); % should be length(epidTimes)

for s = 1:numSims
    
	s
    figs(s) = figure;
    thisP = 1;
    thisSim = sims(s);
    numExtinctions = 0;

    % Load files that are simulation- but not neighborhood-specific
    demTimes = load( strcat('tr_',int2str(treatment),'_sim_',int2str(thisSim),'_dem_times'));
    hholdDist = load( strcat('tr_',int2str(treatment),'_sim_',int2str(thisSim),'_hh_dist'));
    demTimes = load( strcat('tr_',int2str(treatment),'_sim_',int2str(thisSim),'_dem_times'));
    epidTimes = load( strcat('tr_',int2str(treatment),'_sim_',int2str(thisSim),'_epid_times'));
    coinfections = load( strcat('tr_',int2str(treatment),'_sim_',int2str(thisSim),'_coinfection_dist_pneumo'));
    thetas = load( strcat('tr_',int2str(treatment),'_sim_',int2str(thisSim),'_theta'));
    numInfecteds = load( strcat('tr_',int2str(treatment),'_sim_',int2str(thisSim),'_totCarriage'));
    demTimes = demTimes./365; % convert to years
    epidTimes = epidTimes./365;
    infections = zeros(length(epidTimes),length(ages),numStrains,numNeighborhoods);
    infecteds = zeros(length(epidTimes),length(ages),numStrains,numNeighborhoods);
    effInfectedsEveryN = zeros(length(epidTimes),length(cohorts),numStrains,numNeighborhoods);
    effInfectionsEveryN = zeros(length(epidTimes),length(cohorts),numStrains,numNeighborhoods);
    effInfecteds = zeros(length(epidTimes),length(cohorts),numStrains);
    effInfections = zeros(length(epidTimes),length(cohorts),numStrains);
    effAgeDist = zeros(length(demTimes),length(cohorts),numNeighborhoods);
    ageDist = zeros(length(demTimes),length(cohorts));
    paredCoinfections = coinfections(size(coinfections,1)-sampleDuration+1:size(coinfections,1),1:coinfectionBuffer);
    totalCoinfections(:,s) = sum(paredCoinfections,2);
    
    if ( Hflu_simulated==1 )
        hflu_coinfections = load( strcat('tr_',int2str(treatment),'_sim_',int2str(thisSim),'_coinfection_dist_hflu-pneumo'));  
        paredHfluCoinfections = hflu_coinfections(size(hflu_coinfections,1)-sampleDuration+1:size(hflu_coinfections,1),1:numPneumoSerotypes);
        totalHfluPneumoCoinfections = mean(paredHfluCoinfections,1);
    end
    
    for n = 1:numNeighborhoods
        thisNhood = n-1;
        
        % Load files that are simulation- and neighborhood-specific
        for stype = 1:numStrains
            thisSType = serotypes(stype);
            theseInfections = load( strcat('tr_',int2str(treatment),'_sim_',int2str(thisSim),'_infections_',int2str( thisSType ),'_neighborhood_',int2str( thisNhood )) );
            theseInfecteds = load( strcat('tr_',int2str(treatment),'_sim_',int2str(thisSim),'_infecteds_',int2str( thisSType ),'_neighborhood_',int2str(thisNhood)) );
            infections(:,:,stype,n) = theseInfections;
            infecteds(:,:,stype,n ) = theseInfecteds;
            theseAbundances = sum(theseInfecteds,2);
            finalTS(stype,:,s) = finalTS(stype,:,s) + theseAbundances(length(theseAbundances)-sampleDuration+1:length(theseAbundances))';
        end
        
        % Convert age-associated data to cohort-associated data
        ageDist = load( strcat('tr_',int2str(treatment),'_sim_',int2str(thisSim),'_age_dist_neighborhood_',int2str(thisNhood)));
        for r = 1:nrow
            effAgeDist(:,r,n) = sum(ageDist(:,cohortIndices(r,1):cohortIndices(r,2)),2);
        end
        
        for ss = 1:numStrains
            tempInfecteds = infecteds(:,:,ss,n);
            tempInfections = infections(:,:,ss,n);
            for r = 1:nrow
                effInfectedsEveryN(:,r,ss,n) = sum( tempInfecteds(:,cohortIndices(r,1):cohortIndices(r,2)),2 );
                effInfectionsEveryN(:,r,ss,n) = sum( tempInfections(:,cohortIndices(r,1):cohortIndices(r,2)),2 );
                effInfecteds(:,r,ss) = effInfecteds(:,r,ss) + sum( tempInfecteds(:,cohortIndices(r,1):cohortIndices(r,2)),2 );
                effInfections(:,r,ss) = effInfections(:,r,ss) + sum( tempInfections(:,cohortIndices(r,1):cohortIndices(r,2)),2 );
            end
        end
    
    end
    ageDist = sum(effAgeDist,3);
    
    
    % Plot demographic dynamics
    epidStart = epidTimes(1);
    
    % Population size and age structure over time
    subplot(rowSP,colSP,thisP)
    h = bar(demTimes,ageDist,'stack');
    xlim([0 max(demTimes)])
    set(gca,'XTick',[min(demTimes):demInt:max(demTimes)],'FontName','Helvetica','FontSize',9) %'EdgeColor','none');
    xlabel('Time (y)','FontName','Helvetica','FontSize',9)
    ylabel('Number of people','FontName','Helvetica','FontSize',9)
    ylim([0 max(sum(ageDist,2))])
    hold on
    bar(epidStart,max(sum(ageDist,2)),'r') % consider proper line plot overlay
    legend(cohortNames)
    thisP = thisP + 1;
    
    % 	Final household sizes
    subplot(rowSP,colSP,thisP)
    thisMax = max(find(hholdDist(length(hholdDist),:) > 0 ));
    bar(hholdDist(length(hholdDist),:)./sum(hholdDist(length(hholdDist),:)));
    xlabel('Household size','FontName','Helvetica','FontSize',9)
    ylabel({'Fraction of households';'at end of simulation'},'FontName','Helvetica','FontSize',9)
    set(gca,'FontName','Helvetica','FontSize',9)
    xlim([0 thisMax+1])
    title(['Treatment value = ', num2str(treatmentVal),', simulation ',num2str(s),', ',modelName],'FontName','Helvetica','FontSize',12)
    thisP = thisP + 1;
    box on;
    grid on;
 
    
    % Epidemiological dynamics
     
    % Calculate means of theta for each cohort for this sim
    for c = 1:length(durationCohorts)
        firstFilter = find( thetas(:,1) > durationCohorts(c)*365 );
        if ( c ~= length(durationCohorts))
            secondFilter = find(thetas(:,1) < durationCohorts(c+1)*365);
            netFilter = intersect(firstFilter,secondFilter);
            thetasByAge(c,s) = mean(thetas(netFilter,2));
        else
            thetasByAge(c,s) = mean(thetas(firstFilter,2));
        end
    end
    
    
    % First find maxima to use in setting axis limits
    epidAgeDist = effAgeDist(find(demTimes>=epidStart),:,:);
    epidAgeDist = sum(epidAgeDist,3);
    maxPrev = 0;
    maxInf = 0;
    maxTotInf = 0;
    for t = 1:numPneumoSerotypes
        maxPrev = max(maxPrev, max(max(effInfecteds((length(epidTimes)-sampleDuration):length(epidTimes),:,t)./epidAgeDist((length(epidTimes)-sampleDuration):length(epidTimes),:,:))));
        maxInf = max(maxInf, max(max((effInfections((length(epidTimes)-sampleDuration):length(epidTimes),:,t)))) );
        maxTotInf = max(maxTotInf, max(sum(effInfecteds((length(epidTimes)-sampleDuration):length(epidTimes),:,t),2)));
        maxTotInf = max(maxTotInf, max(sum(effInfections((length(epidTimes)-sampleDuration):length(epidTimes),:,t),2)));
    end 
    
    % Record total fraction of people infected overall
    totCarriagePrevalence = numInfecteds./sum(epidAgeDist,2); % total population
    finalPrevalences(s) = mean(totCarriagePrevalence); % finalPrevalences in kids under 5
    under5TotCarriagePrevalence = sum(coinfections,2)./sum(epidAgeDist(:,1:3),2);
    over5TotCarriagePrevalence = (numInfecteds - sum(coinfections,2))./sum(epidAgeDist(:,4:5),2);
    
    % Plot fraction colonized with each pneumo strain (same plot)
    subplot(rowSP,colSP,thisP)
    totFracInfecteds = zeros(length(epidTimes),length(serotypes));
    under5FracInfecteds = zeros(length(epidTimes),length(serotypes));
    for i = 1:numPneumoSerotypes
        totFracInfecteds(:,i) = sum(effInfecteds(:,:,i),2); 
        totFracInfecteds(:,i) = totFracInfecteds(:,i)./sum(epidAgeDist,2);
        under5FracInfecteds(:,i) = sum(effInfecteds(:,1:3,i),2);
        under5FracInfecteds(:,i) = under5FracInfecteds(:,i)./sum(epidAgeDist(:,1:3),2);
        % recall finalPrevs = zeros(length(serotypes),sampleDuration,length(sims)); % holds prevalences
        finalPrevs(i,:,s) = totFracInfecteds((size(totFracInfecteds,1)-sampleDuration+1):size(totFracInfecteds,1),i);
    end
    thisMap = flipud(colormap(jet(numPneumoSerotypes))); % color map for serotypes -- add one more color for total prev
    [ax h1 h2] = plotyy(epidTimes,under5FracInfecteds,epidTimes,under5TotCarriagePrevalence);
    for p = 1:numPneumoSerotypes
        set(h1(p),'color',thisMap(p,:));
    end
    set(h2,'color','k');
    set(ax(2),'ycolor','k');
    hold on
    box off;
    set(get(ax(1),'ylabel'),'String',{'Fraction of people';'colonized with serotype'},'FontName','Helvetica','FontSize',9);
    set(get(ax(2),'ylabel'),'String','Total carriage prevalence','FontName','Helvetica','FontSize',9);
    xlabel('Time (y)','FontName','Helvetica','FontSize',9);
    set(ax(1),'ylim',[0 1.5*max(max(finalPrevs(:,:,s)))]);
    set(ax(2),'ylim',[0 0.2]);
    set(ax(1),'xlim',[min(epidTimes) max(epidTimes)]);
    set(ax(2),'xlim',[min(epidTimes) max(epidTimes)]);
    set(ax(1),'YTickMode','auto','FontName','Helvetica','FontSize',9)
    ylimits = get(ax(2),'YLim');
    yinc = (ylimits(2)-ylimits(1))/4;
	set(ax(2),'YTickMode','auto','FontName','Helvetica','FontSize',9)
    if ( numStrains < minToPlot )
        legend(serotypeNames)
        legend('boxoff')
    end

	carriageDiffsByYear(:,s) = under5TotCarriagePrevalence./under5TotCarriagePrevalence(preVaccineYear);
	carriageDiffsByYearOver5(:,s) = over5TotCarriagePrevalence./over5TotCarriagePrevalence(preVaccineYear);
    numExtinct = sum( totFracInfecteds(length(totFracInfecteds),1:numPneumoSerotypes) < extThold );
    extinctions(numExtinct+1) = extinctions(numExtinct+1) + 1; 
    thisP = thisP + 1;

    % Output totCarriagePrevalence, tot under 5, tot >=5
    filename = strcat('tr_',int2str(treatment),'_sim_',int2str(sims(s)),'_totFracInfecteds_',modelName,'.csv');
    csvwrite(filename,totFracInfecteds);
    filename = strcat('tr_',int2str(treatment),'_sim_',int2str(sims(s)),'_under5FracInfecteds_',modelName,'.csv');
    csvwrite(filename,under5FracInfecteds);
    filename = strcat('tr_',int2str(treatment),'_sim_',int2str(sims(s)),'_under5TotCarriagePrevalence_',modelName,'.csv');
    csvwrite(filename,under5TotCarriagePrevalence);
    filename = strcat('tr_',int2str(treatment),'_sim_',int2str(sims(s)),'_over5TotCarriagePrevalence_',modelName,'.csv');
    csvwrite(filename,over5TotCarriagePrevalence);    
    
    
    % Output Simpson Index over time for whole pop and young kids
    simpsonTS = zeros(length(epidTimes),1);
    simpsonTSUnder5 = zeros(length(epidTimes),1);
    for thisTime = 1:length(epidTimes)
        for thisS = 1:numPneumoSerotypes
            simpsonTS(thisTime) = simpsonTS(thisTime) + ( totFracInfecteds(thisTime,thisS)./sum(totFracInfecteds(thisTime,:)))^2;
            simpsonTSUnder5(thisTime) = simpsonTS(thisTime) + ( under5FracInfecteds(thisTime,thisS)./sum(under5FracInfecteds(thisTime,:)))^2;
        end
    end
    simpsonTS = 1 - simpsonTS;
    simpsonTSUnder5 = 1 - simpsonTSUnder5;
    simpsonTSUnder5AllSims(:,s) = simpsonTSUnder5;
    filename = strcat('tr_',int2str(treatment),'_sim_',int2str(sims(s)),'_simpsonTS_',modelName,'.csv');
    csvwrite(filename,simpsonTS)
    filename = strcat('tr_',int2str(treatment),'_sim_',int2str(sims(s)),'_simpsonTSUnder5_',modelName,'.csv');
    csvwrite(filename,simpsonTSUnder5)
    
 
    % Plot frequencies changing over time
    subplot(rowSP,colSP,thisP)
    theseFrequencies = zeros(numPneumoSerotypes,sampleDuration);
    for y = 1:sampleDuration
        theseFrequencies(:,y) = finalPrevs(1:numPneumoSerotypes,y,s)./sum(finalPrevs(1:numPneumoSerotypes,y,s));
    end
    theseFrequencies = theseFrequencies';
    area(theseFrequencies)
    colormap(flipud(colormap(jet(numPneumoSerotypes))))
    axis tight
    set(gca,'FontName','Helvetica','FontSize',9)
    ylabel({'Fraction of total';'pneumococcal prevalence'},'FontName','Helvetica','FontSize',9);
    xlabel('Time (y)','FontName','Helvetica','FontSize',9);
    set(gca,'XTickLabel',[(max(epidTimes)-sampleDuration+1):round(sampleDuration/5):max(epidTimes)],'xtick',[1:round(sampleDuration/5):sampleDuration],'FontName','Helvetica','FontSize',9) %'EdgeColor','none');
    thisP = thisP + 1;
    
    
	% Plot serotype diversity for major age groups
    % Currently plots MEAN OF SAMPLE DURATION FOR EACH COHORT
    subplot(rowSP,colSP,thisP)
    relevantInfecteds = effInfecteds((length(epidTimes)-sampleDuration+1):length(epidTimes),:,1:numPneumoSerotypes);
    meanFracInfecteds = zeros(sampleDuration,length(cohortNames),numPneumoSerotypes);
    for y = 1:sampleDuration
        for r = 1:length(cohortNames)
            totInfectedsThisRow = sum(relevantInfecteds(y,r,:));
            meanFracInfecteds(y,r,:) = relevantInfecteds(y,r,:)./totInfectedsThisRow;
        end
    end
    meanOverTime = mean(meanFracInfecteds,1);
    meanStuff = zeros(size(meanOverTime,2),size(meanOverTime,3));
    meanStuff(1:length(cohortNames),:) = meanOverTime(:,1:length(cohortNames),1:numPneumoSerotypes);
    bar(meanStuff,'stack');
    colormap(flipud(colormap(jet(numPneumoSerotypes))))
    axis tight
    ylabel('Final mean fraction of total carriage','FontName','Helvetica','FontSize',9)
    xlabel('Age (y)','FontName','Helvetica','FontSize',9)
    set(gca,'FontName','Helvetica','FontSize',9)
    set(gca,'XTickLabel',cohortNames)
    thisP = thisP + 1;
    
    % Calculate's Simpson Index for each cohort
    for r = 1:length(cohortNames)
       thisSimpsonVal = 0;
       for thisS = 1:numPneumoSerotypes
            % thisSimpsonVal = thisSimpsonVal + meanStuff(r,thisS)^2; % for
            % mean
            thisSimpsonVal = thisSimpsonVal + mean( meanFracInfecteds(1:sampleDuration,r,thisS).^2 );
        end
       simpsonAges(r,s) = 1-thisSimpsonVal;      
    end
    

    % Final coinfection distribution of pneumo
    subplot(rowSP,colSP,thisP)
    thisP = thisP + 1;
    thisMax = max(find(coinfections(size(coinfections,1),:) > 0 ));
    theseCoinfections = coinfections(size(coinfections,1),:)./sum(coinfections(size(coinfections,1),:));
    bar( theseCoinfections );
    pneumoCocols(:,s) = theseCoinfections(1:10);
    grid on;
    xlabel('Pneumococcal colonizations per person','FontName','Helvetica','FontSize',9)
    ylabel({'Fraction of colonized people';'at end of simulation'},'FontName','Helvetica','FontSize',9)
    set(gca,'XTickLabel',[1:thisMax],'xtick',[1:t],'FontName','Helvetica','FontSize',9) %'EdgeColor','none');
    xlim([ 0 thisMax+0.5])
    set(gca,'FontName','Helvetica','FontSize',9)
    
    
    % Plot H. flu stuff, if simulated
    if ( Hflu_simulated == 1 ) 
        
       % Plot prevalence by age
       subplot(rowSP,colSP,thisP)
       thisP = thisP + 1;
       plot(epidTimes,effInfecteds(:,:,Hflu_index)./epidAgeDist);
       xlim([min(epidTimes) max(epidTimes)]+1)
       set(gca,'XTick',[min(epidTimes):timeInt:max(epidTimes)],'FontName','Helvetica','FontSize',9) %'EdgeColor','none');
       ylabel('Fraction colonized with H. influenzae','FontName','Helvetica','FontSize',9);
       xlabel('Time (y)','FontName','Helvetica','FontSize',9);
       %ylim([0 maxPrev])
       legend(cohortNames)
       legend('boxoff')
      
       % Plot colonization rates by age
       subplot(rowSP,colSP,thisP)
       thisP = thisP + 1;
       plot(epidTimes,effInfections(:,:,Hflu_index));
       xlim([min(epidTimes) max(epidTimes)]+1)
       set(gca,'XTick',[min(epidTimes):timeInt:max(epidTimes)],'FontName','Helvetica','FontSize',9) %'EdgeColor','none');
       ylabel( 'Number of people colonized' );
       xlabel('Time (y)','FontName','Helvetica','FontSize',9);
       ylim([0 1.2*maxInf])
       legend('boxoff');
       
       % Plot co-colonizations with pneumo
       subplot(rowSP,colSP,thisP)
       thisP = thisP + 1;
       totalHfluPneumoCoinfections = mean(paredHfluCoinfections,1);
       bar(0:(numPneumoSerotypes-1),totalHfluPneumoCoinfections);
       xlim([-0.5 numPneumoSerotypes])
       set(gca,'XTick',0:(numPneumoSerotypes-1),'FontName','Helvetica','FontSize',9) %'EdgeColor','none');
       xlabel('Serotype','FontName','Helvetica','FontSize',9);
       ylabel({'Mean no. of co-colonizations'; 'with H. influenzae'},'FontName','Helvetica','FontSize',9);
 
    end
    
     
    for t = 1:min(minToPlot,numPneumoSerotypes) % only plots first ten     
        thisType = serotypes(t);

        %   Prevalence by age
        subplot(rowSP,colSP,thisP)
        thisP = thisP + 1;
        plot(epidTimes,effInfecteds(:,:,t)./epidAgeDist);
        xlim([min(epidTimes) max(epidTimes)]+1)
        set(gca,'XTick',[min(epidTimes):timeInt:max(epidTimes)],'FontName','Helvetica','FontSize',9) %'EdgeColor','none');
        someText = strcat('Serotype',{' '}, int2str( serotypes(t)),':' );
        ylabel({char(someText);'Fraction colonized' }); %,'FontName','Helvetica','FontSize',9);
        xlabel('Time (y)','FontName','Helvetica','FontSize',9);
        ylim([0 1.2*maxPrev])
        if ( Hflu_simulated == 0 && t == 1 )
            legend(cohortNames)
        end
        legend('boxoff')
        
        %   Infections by age
        subplot(rowSP,colSP,thisP)
        thisP = thisP + 1;
        plot(epidTimes,effInfections(:,:,t));
        xlim([min(epidTimes) max(epidTimes)]+1)
        set(gca,'XTick',[min(epidTimes):timeInt:max(epidTimes)],'FontName','Helvetica','FontSize',9) %'EdgeColor','none');
        ylabel( 'Number of people colonized' );
        xlabel('Time (y)','FontName','Helvetica','FontSize',9);
        ylim([0 1.2*maxInf])
        legend('boxoff');
        
        if ( Hflu_simulated == 0 )
            %   Total infecteds and infections
            subplot(rowSP,colSP,thisP)
            thisP = thisP + 1;
            plot(epidTimes,sum(effInfecteds(:,:,t),2),epidTimes,sum(effInfections(:,:,t),2));
            xlim([min(epidTimes) max(epidTimes)]+1)
            set(gca,'XTick',[min(epidTimes):timeInt:max(epidTimes)],'FontName','Helvetica','FontSize',9) %'EdgeColor','none');
            ylabel('Total number','FontName','Helvetica','FontSize',9)
            xlabel('Time (y)','FontName','Helvetica','FontSize',9);
            ylim([0 1.2*maxTotInf]);
            if ( t == 1 )
                legend('Colonized','Colonizations')
            end
            legend('boxoff');
        else
            % Show fraction of all colonizations in people with H. flu
            subplot(rowSP,colSP,thisP)
            thisP = thisP + 1;
            fracWithHflu = hflu_coinfections(:,t)./sum(effInfecteds(:,:,t),2);
            bar(epidTimes,fracWithHflu);
            %hold on
            %plot(epidTimes,hflu_coinfections(:,t),epidTimes,sum(effInfecteds(:,:,t),2));
            xlim([min(epidTimes) max(epidTimes)]+1)
            %ylim([0 1])
            set(gca,'XTick',[min(epidTimes):timeInt:max(epidTimes)],'FontName','Helvetica','FontSize',9) %'EdgeColor','none');
            ylabel({'Fraction co-colonized';'with H. flu'},'FontName','Helvetica','FontSize',9);
            xlabel('Time (y)','FontName','Helvetica','FontSize',9);
        end        
    end % end serotypes to plot 
    
    % Save co-colonization with H. flu data
    if ( Hflu_simulated == 1 )
        for t = 1:numPneumoSerotypes
            percentHflu(t,:,s) = hflu_coinfections(:,t)./sum(effInfecteds(:,:,t),2);
        end
        finalPercentHflu(:,:,s) = percentHflu(:,(max(epidTimes)-sampleDuration+1):max(epidTimes),s);
    end
    
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [12 5+3*min(minToPlot,length(serotypes))]);
    set(gcf, 'PaperPosition',[0,0,12,6+3*min(minToPlot,length(serotypes))])
    filename = strcat('tr_',int2str(treatment),'_sim_',int2str(sims(s)),'_results_',modelName);
    saveas(gcf,filename,'epsc');

	close all
    
    % Now calculate changes in carriage before and after vaccination
    % in kids
    preVaccRange = 90:99;
    postVaccRange = 110:119;
    preVaccPrevs = under5TotCarriagePrevalence(preVaccRange);
    postVaccPrevs = under5TotCarriagePrevalence(postVaccRange);
    percentReduction = (mean(preVaccPrevs)-mean(postVaccPrevs))/mean(preVaccPrevs);
    reductions(s) = percentReduction;
    for thisType = 1:numPneumoSerotypes
		prePostPrevs(thisType,s) = mean(under5FracInfecteds(postVaccRange,thisType))./mean(under5FracInfecteds(preVaccRange,thisType));
    end
    
end % end for this sim

% Output shift b/w pre- and post-PCV prevs
prePostOutput = [ mean(prePostPrevs,2) std(prePostPrevs,0,2)] ;
filename = strcat('tr_',int2str(treatment),'_prePostPrevs_',modelName,'.csv');
csvwrite(filename,prePostOutput);

% Output carriage diffs by year
carriageDiffsMeans = mean(carriageDiffsByYear,2);
carriageDiffsSTDs = std(carriageDiffsByYear,0,2);
carriageDiffsSummary = [ carriageDiffsMeans carriageDiffsSTDs ];
filename = strcat('tr_',int2str(treatment),'_totCarriageDiffs_',modelName,'.csv');
csvwrite(filename,carriageDiffsSummary);

carriageDiffsOver5Means = mean(carriageDiffsByYearOver5,2);
carriageDiffsOver5STDs = std(carriageDiffsByYearOver5,0,2);
carriageDiffsOver5Summary = [ carriageDiffsOver5Means carriageDiffsOver5STDs ];
filename = strcat('tr_',int2str(treatment),'_totCarriageDiffsOver5_',modelName,'.csv');
csvwrite(filename,carriageDiffsOver5Summary);

simpsonUnder5Means = mean(simpsonTSUnder5AllSims,2);
simpsonUnder5STDs = std(simpsonTSUnder5AllSims,0,2);
simpsonTSSummary = [ simpsonUnder5Means simpsonUnder5STDs ];
filename = strcat('tr_',int2str(treatment),'_simpsonTSUnder5_',modelName,'.csv');
csvwrite(filename,simpsonTSSummary);


% Summarize extinctions
figure
bar(0:numPneumoSerotypes,extinctions)
xlim([-0.5 numPneumoSerotypes+0.5])
xlabel('Number of serotypes extinct at simulation end','FontName','Helvetica','FontSize',9)
ylabel('Number of simulations','FontName','Helvetica','FontSize',9)
title(['Extinctions at ', num2str(extThold),' threshold'],'FontName','Helvetica','FontSize',9)
set(gca,'FontName','Helvetica','FontSize',9)
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [2.5 2.5]);
set(gcf, 'PaperPosition',[0,0,2.5,2.5])
filename = strcat( 'tr_',int2str(treatment),'_extinction_summary_',modelName);
saveas(gcf,filename,'eps');
    

% Recall finalTS = zeros(numStrains,sampleDuration,length(sims)); % holds final time series for each simulation
meanAbundances = zeros(numPneumoSerotypes,numSims);
meanPrevalences = zeros(numPneumoSerotypes,numSims);
for s = 1:numSims
    for type = 1:numPneumoSerotypes
       meanAbundances(type,s) = mean(finalTS(type,:,s));
       meanPrevalences(type,s) = mean(finalPrevs(type,:,s));
    end
end

simpsonIndex = zeros(numSims,1);
for m = 1:numSims
    for thisS = 1:numPneumoSerotypes
        simpsonIndex(m) = simpsonIndex(m) + mean( ( finalTS(thisS,:,m)./sum(finalTS(1:numPneumoSerotypes,:,m),1) ) .^ 2 );
    end
end


% Calculate correlations across simulations and plot rank abundance (mean +
% std)
endPrevs = zeros(numPneumoSerotypes,numSims);
for s = 1:numSims
    endPrevs(:,s) = finalPrevs(:,size(finalPrevs,2),s); % recall only contains pneumo serotypes and has dim (serotype,sim)
end
denoms = sum(endPrevs,1);
for ( s = 1:numSims )
    endPrevs(:,s) = endPrevs(:,s)./denoms(s); % converts to frequencies (proportional abundance)
end
ranks = zeros(numSims,numPneumoSerotypes);
fitnessRanks = 1:numPneumoSerotypes;
for ( s = 1:numSims )
    [ sortedEndPrevs thisStype ] = sort(endPrevs(:,s),'descend');
    for ( t = 1:length(thisStype))
       ranks(s,thisStype(t)) = t;
    end
end
% Output correlations between simulations (rankCorrs) and b/w sims and
% fitness (fitCorrs)
rankCorrs = corr(ranks','type','Spearman'); % ranks here indexed (stype,sim); show correlations across sims
offDiags = tril(rankCorrs,-1);
offDiags = nonzeros(offDiags);
rankCorrs = [ mean(offDiags) std(offDiags) ];
fitRanks = [ fitnessRanks; ranks ]; % makes new matrix w/ fitness on top
rankFitCorrs = corr(fitRanks','type','Spearman');
meanFitCorr = mean(rankFitCorrs(1,2:numSims+1));
stdFitCorr = std(rankFitCorrs(1,2:numSims+1));
fitCorrs = [ meanFitCorr stdFitCorr ];

% Produce plot showing ranks and mean abundances (from endPrevs)
[ sortedPrevs sortedix ] = sort(endPrevs,1,'descend'); % still indexed (stype,sim)
errorbar(1:numPneumoSerotypes,mean(sortedPrevs,2),std(sortedPrevs,0,2),'xk');%,'MarkerFaceColor','k');
xlabel('Rank','FontName','Helvetica','FontSize',9)
ylabel({'Frequency'},'FontName','Helvetica','FontSize',9)
set(gca,'FontName','Helvetica','FontSize',9)
axis tight
axis square
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [2 2]);
set(gcf, 'PaperPosition',[0,0,2,2])
filename = strcat('tr_',int2str(treatment),'_rankAbundance_',modelName);
print('-deps',filename)
filename = strcat(filename,'.csv');
csvwrite(filename,[ mean(sortedPrevs,2) std(sortedPrevs,0,2)]);

% Calculate mean shared serotypes in top 10
numShared = zeros((numSims^2-numSims)/2,1);
marker = 1;
for s1 = 1:numSims
    s1types = unique(sortedix(1:10,s1));
    for s2 = (s1+1):numSims
        s2types = unique(sortedix(1:10,s2));
        numShared(marker) = sum(intersect(s1types,s2types)>0);
        marker = marker + 1;
    end
end
shared = [ mean(numShared) std(numShared) ];


% Plot relative prevalence for all sims
        % recall finalPrevs =
        % zeros(length(serotypes),sampleDuration,length(sims)); % holds prevalences
meanFinalPrevs = mean(finalPrevs,2);
finalAbsRelativePrevs = zeros(length(serotypes),length(sims));
for ( s = 1:numSims )
    for ( t = 1:numPneumoSerotypes)
        plot((finalPrevs(t,:,s)-meanFinalPrevs(t,s))./meanFinalPrevs(t,s));
        finalAbsRelativePrevs(t,s) = mean(abs((finalPrevs(t,:,s)-meanFinalPrevs(t,s))./meanFinalPrevs(t,s)));
        set(gca,'ColorOrder',flipud(colormap(jet(numPneumoSerotypes))));
        hold all;
    end
    xlim([min(epidTimes) max(epidTimes)])
    axis tight
    axis square
    box on;
    set(gca,'XTickLabel',[(max(epidTimes)-sampleDuration+1):round(sampleDuration/5):max(epidTimes)],'xtick',[1:round(sampleDuration/5):sampleDuration],'FontName','Helvetica','FontSize',9) %'EdgeColor','none');
    ylabel({'Normalized prevalence'},'FontName','Helvetica','FontSize',9);
    xlabel('Time (y)','FontName','Helvetica','FontSize',9);
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [2.5 2.5]);
    set(gcf, 'PaperPosition',[0,0,2.5,2.5])
    filename = strcat('tr_',int2str(treatment),'_sim_',int2str(sims(s)),'_relativePrevalence_',modelName);
    saveas(gcf,filename,'epsc');
    close all
end
errorbar(mean(finalAbsRelativePrevs,2),std(finalAbsRelativePrevs,0,2),'.k')
axis tight
axis square
set(gca,'FontName','Helvetica','FontSize',9)
xlabel('Serotype','FontName','Helvetica','FontSize',9)
ylabel({'Mean absolute';'normalized prevalence'},'FontName','Helvetica','FontSize',9)    
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [6 6]);
set(gcf, 'PaperPosition',[0,0,6,6])
filename = strcat('tr_',int2str(treatment),'_relativePrevalence_',modelName);
saveas(gcf,filename,'eps');
csvwrite(strcat(filename,'.csv'),[mean(finalAbsRelativePrevs,2) std(finalAbsRelativePrevs,0,2)]);

% Simpson Index for major age groups
simpsonAgeMean = mean(simpsonAges,2);
simpsonAgeSD = std(simpsonAges,0,2);

% Calculate diversity index
simpsonDiversity = mean(simpsonIndex);
simpsonSD = std(simpsonIndex);
simpsonInfo = [ 1-simpsonDiversity simpsonSD ];

meanPrevalence = mean( finalPrevalences  );
prevalenceSD = std( finalPrevalences);
prevalenceInfo = [ meanPrevalence prevalenceSD ];

% Output pneumo cocolonization info
meanCocols = mean(pneumoCocols,2);
SDCocols = std(pneumoCocols,0,2);

% Calculate durations of carriage
minCarriage = 25;
intrinsicDuration = 220;
epsilon = 0.25;
thetaInfo = [ mean(thetasByAge,2) std(thetasByAge,0,2)];
meanThetas = mean(thetasByAge,2);
stdThetas = std(thetasByAge,0,2);
meanDurations = zeros(length(meanThetas),1);
stdDurations = zeros(length(meanThetas),1);
for ( d = 1:length(meanDurations) )
    meanDurations(d) = minCarriage + (intrinsicDuration - minCarriage )*exp(-epsilon*meanThetas(d));
    upperDuration = minCarriage + ( intrinsicDuration - minCarriage )*exp(-epsilon*( meanThetas(d) + stdThetas(d) ) );
    lowerDuration = minCarriage + ( intrinsicDuration - minCarriage )*exp(-epsilon*( meanThetas(d) - stdThetas(d) ) );
    stdDurations(d) = (upperDuration - lowerDuration)/2;
end
errorbar(meanDurations,stdDurations,'.k')
axis tight
axis square
set(gca,'FontName','Helvetica','FontSize',9)
xlabel('Age (y)','FontName','Helvetica','FontSize',9)
ylabel({'Mean carriage duration of';'intrisically fittest serotype (days)'},'FontName','Helvetica','FontSize',9) 
set(gca,'XTickLabel',durationNames,'xtick',[1:length(meanDurations)],'FontName','Helvetica','FontSize',9) %'EdgeColor','none');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [6 6]);
set(gcf, 'PaperPosition',[0,0,6,6])
filename = strcat('tr_',int2str(treatment),'_durationsByAge_',modelName);
saveas(gcf,filename,'eps');

% Hflu stuff
hfluInfo = zeros(numPneumoSerotypes,numSims);
hfluInfoMeans = zeros(numPneumoSerotypes,1);
hfluInfoSD = zeros(numPneumoSerotypes,1);

if ( Hflu_simulated == 1 )
    for s = 1:numSims
        hfluInfo(:,s) = nanmean(finalPercentHflu(:,:,s),2);
    end

    hfluInfoMeans = mean(hfluInfo,2);
    hfluInfoSD = std(hfluInfo');
    errorbar(1:numPneumoSerotypes,hfluInfoMeans,hfluInfoSD);
    xlabel('Serotype','FontName','Helvetica','FontSize',9)
    ylabel({'Mean percent co-colonized';'with H. influenzae'},'FontName','Helvetica','FontSize',9)
    axis tight
    axis square
    set(gca,'FontName','Helvetica','FontSize',9)
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [2 2]);
    set(gcf, 'PaperPosition',[0,0,2,2])
    filename = strcat('tr_',int2str(treatment),'_hflu-cocol_summary_',modelName);
    saveas(gcf,filename,'eps');
end

% Plot final proportion of each serotype (by average prevalence)
% and print diversity at top
subplot(2,1,1)
bar(0:(numPneumoSerotypes-1),nanmean(meanPrevalences,2));
xlabel('Serotype','FontName','Helvetica','FontSize',9)
ylabel(['Mean prevalence'],'FontName','Helvetica','FontSize',9)
title(['Simpson index = ',num2str(simpsonDiversity),', SD=', num2str(simpsonSD)],'FontName','Helvetica','FontSize',9)
set(gca,'FontName','Helvetica','FontSize',9)
axis tight
grid on

subplot(2,1,2)
bar(0:(numPneumoSerotypes-1),nanmean(meanPrevalences,2)./sum(nanmean(meanPrevalences,2)));
xlabel('Serotype','FontName','Helvetica','FontSize',9)
ylabel(['Frequency'],'FontName','Helvetica','FontSize',9)
set(gca,'FontName','Helvetica','FontSize',9)
axis tight
grid on

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [2.5 4]);
set(gcf, 'PaperPosition',[0,0,2.5,4])
filename = strcat('tr_',int2str(treatment),'_diversity_summary_',modelName);
saveas(gcf,filename,'epsc');

kidReducts = [ mean(reductions) std(reductions)];

close all
        
end








