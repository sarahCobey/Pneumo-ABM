% Treatment processor
% Pneumo-ABM - (c) 2012, S. Cobey


%% Processing treatments
modelName = 'sampleModel';
numSims = 10; % number of replicates for each model
numModels = 1; % models are named above; can here analyze results for multiple models (assuming same treatments and numSims for each)
Hflu_simulated = 0; % 0 = false, 1 = true
extThold = 0.001; % threshold prevalence for considering serotype extinct
numNeighborhoods = 1;
numPneumoSerotypes = 25;

treatmentVals = 0:0.1:0.9; % values of serotype-specific immunity
numTreatments = length(treatmentVals);
simpsonVals = zeros(numTreatments,2,numModels);
prevVals = zeros(numTreatments,2,numModels);
rankCorrelations = zeros(numTreatments,2,numModels);
rankFitCorrelations = zeros(numTreatments,2,numModels);
thetasByAgeMean = zeros(numTreatments,6,numModels);
thetasByAgeSD = zeros(numTreatments,6,numModels);
errors = zeros(numTreatments,1,numModels);
fittedBetas = zeros(numTreatments,1,numModels);
hfluFreqs = zeros(numTreatments,numPneumoSerotypes,numSims,numModels);
pneumoCocolsMean = zeros(numTreatments,10,numModels);
pneumoCocolsSD = zeros(numTreatments,10,numModels);
simpsonAgeMeans = zeros(numTreatments,5,numModels);
simpsonAgeSDs = zeros(numTreatments,5,numModels);
shared = zeros(numTreatments,2,numModels);
thisModel = 1;

betas = load('Betas_used.txt');

for t = 1:numTreatments
    treatmentVal = treatmentVals( t )
    [ simpsonStuff totalPrev thetaInfoMean thetaInfoSD hfluInfo error rankCorrs fitCorrs meanCocols SDCocols simpsonAgeMean simpsonAgeSD sharedStuff] = ABM_processor(1:numSims,extThold,t,treatmentVal,numNeighborhoods,modelName,Hflu_simulated);
    simpsonVals(t,1,thisModel) = simpsonStuff(1);
    simpsonVals(t,2,thisModel) = simpsonStuff(2);
    prevVals(t,1,thisModel) = totalPrev(1);
    prevVals(t,2,thisModel) = totalPrev(2);
	rankCorrelations(t,1,thisModel) = rankCorrs(1);
	rankCorrelations(t,2,thisModel) = rankCorrs(2);
	rankFitCorrelations(t,1,thisModel) = fitCorrs(1);
	rankFitCorrelations(t,2,thisModel) = fitCorrs(2);
	thetasByAgeMean(t,:,thisModel) = thetaInfoMean;
	thetasByAgeSD(t,:,thisModel) = thetaInfoSD;
    errors(t,1,thisModel) = error;
    beta = load(strcat('tr_',int2str(t),'_sim_',int2str(1),'_BETA'));
    fittedBetas(t,1,thisModel) = beta(1);
	pneumoCocolsMean(t,:,thisModel ) = meanCocols;
	pneumoCocolsSD(t,:,thisModel ) = SDCocols;
	simpsonAgeMeans(t,:,thisModel) = simpsonAgeMean;
	simpsonAgeSDs(t,:,thisModel)=simpsonAgeSD;
	shared(t,1,thisModel) = sharedStuff(1);
	shared(t,2,thisModel) = sharedStuff(2);
	if ( Hflu_simulated == 1 )
		hfluFreqs(t,:,:,thisModel) = hfluInfo;
	end
end

filename = strcat('simpsonVals_',modelName,'.csv');
csvwrite(filename,simpsonVals(:,:,thisModel))
filename = strcat('prevVals_',modelName,'.csv');
csvwrite(filename,prevVals(:,:,thisModel))
filename = strcat('rankCorrs_',modelName,'.csv');
csvwrite(filename,rankCorrelations(:,:,thisModel))
filename = strcat('rankFitCorrs_',modelName,'.csv');
csvwrite(filename,rankFitCorrelations(:,:,thisModel))
filename = strcat('errors_',modelName,'.csv');
csvwrite(filename,errors(:,1,thisModel))
filename = strcat('fittedBetas_',modelName,'.csv');
csvwrite(filename,fittedBetas(:,1,thisModel))
filename = strcat('thetasByAgeMean_',modelName,'.csv');
csvwrite(filename,thetasByAgeMean(:,:,thisModel))
filename = strcat('thetasByAgeSD_',modelName,'.csv');
csvwrite(filename,thetasByAgeSD(:,:,thisModel));
filename = strcat('pneumoCocolsMean_',modelName,'.csv');
csvwrite(filename,pneumoCocolsMean(:,:,thisModel));
filename = strcat('pneumoCocolsSD_',modelName,'.csv');
csvwrite(filename,pneumoCocolsSD(:,:,thisModel));
filename = strcat('simpsonAgeMeans_',modelName,'.csv');
csvwrite(filename,simpsonAgeMeans(:,:,thisModel));
filename = strcat('simpsonAgeSDs_',modelName,'.csv');
csvwrite(filename,simpsonAgeSDs(:,:,thisModel));
filename = strcat('numShared_',modelName,'.csv');
if ( Hflu_simulated == 1 )
	filename = strcat('hfluInfo_',modelName,'.csv');
	csvwrite(filename,hfluFreqs(:,:,thisModel))
end


%% Plot model summaries - currently shows one model
% but load data from files above and add others to same plots

thisModel = 1;
if ( numModels < 4 )
cmap = [0.0784    0.1686    0.5490;
		0.5392    0.0843    0.7745;
		1.0000         0    1.0000 ];
elseif ( numModels == 4 )
cmap = [ 0.0784    0.1686    0.5490;
		0.3856    0.1124    0.6993;
		0.6928    0.0562    0.8497;
		1.0000         0    1.0000 ];
end
shapes = {'o';'s';'^'};

% Simpson Index
figure
errorbar(treatmentVals,simpsonVals(:,1,thisModel),simpsonVals(:,2,thisModel),'o','color',cmap(thisModel,:)); 
legend('Location','SouthEast')
hold on
ylim([0 1])
xlim([0 1])
axis square
xlabel('Anticapsular immunity','FontName','Helvetica','FontSize',9)
ylabel('Simpson Index','FontName','Helvetica','FontSize',9);
set(gca,'FontName','Helvetica','FontSize',9)
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [2 2]);
set(gcf, 'PaperPosition',[0,0,2,2])
filename = strcat('Simpson_diversity_',modelName);
saveas(gcf,filename,'epsc')


% Fitted values of beta
figure		
scatter(treatmentVals(1:length(fittedBetas)),fittedBetas(:,1,thisModel));
grid on
box on
legend('Location','NorthWest')
ylim([0 0.5])
xlim([0 1])
axis square
hold on
xlabel('Serotype-specific immunity','FontName','Helvetica','FontSize',9)
ylabel('Fitted transmission rate (1/day)','FontName','Helvetica','FontSize',9)
set(gca,'FontName','Helvetica','FontSize',9)
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [2 2]);
set(gcf, 'PaperPosition',[0,0,2,2])
filename = strcat('Fitted_betas_',modelName);
saveas(gcf,filename,'epsc')


