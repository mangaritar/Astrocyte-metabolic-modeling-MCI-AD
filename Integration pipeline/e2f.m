% % % 
%         Simple interface for running Exp2Flux on Matlab
%
%
%          Bioinformatics and Systems Biology Laboratory
%     Institute for Genetics - National University of Colombia
%
%             created by: Maria Angarita
%             Last revision by: ampinzonv@unal.edu.co
%             last revised: February 11 2025
% % %


% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%                   USER PROVIDED INFO
%
% The following information has to be provided by the user.
% It assumes the simplest case of use, where the user has a
% metabolic model in .xml or .mat format as well as a 
% .csv file holding the expression data. Make sure this is a two columns file. 
% First colum has to be named: Gene 
% Second column has to be named Expression_Value
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%cobraModel  ='/Users/apinzon/dev/E2F/test_data/e_coli_core.xml';
cobraModel  = 'Astrocyte_Osorio2019 (2).xml';
csvFile    ='C:/Users/u0500/Documents/MATLAB/ast_flux_comparison.csv';
%csvFile     = '/Users/apinzon/dev/E2F/test_data/ast_flux_comparison.csv';
outputDir   = '/Users/apinzon/Documents/MATLAB/e2fOutdir';

% Are you aware of the type of ID used in GPR rules? It can be any of the
% following: kegg, entrez or NA (if not sure which type of ID you have).
typeOfRule  = 'NA';




% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%              !!! WARNING !!!!                        %
% DO NOT TOUCH ANYTHING FROM HERE ON UNLESS YOU ARE     %
% SURE WHAT YOUR DOING. ANY CHANGE COULD PREVENT THIS   %
% SCRIPT FROM FUNCTIONING OR WORST, ALTER THE RESULTS.  %
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%               MAKE SURE INPUT FILES EXIST
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if exist(csvFile, 'file') == 0
    fprintf('Error: File "%s" not found. Exiting.\n', csvFile);
    return;
end

%We can not do the same with cobraModel because if cobraModel is 
%not a file but for instance a variable in the environment such
% as ecoli_core_model.mat which exist only AFTER initCobraToolbox
%is called, then we will have an error.
% This will be checked right after initcobraToolBox is called.


% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%               MAKE SURE outputDir EXISTS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Verificar si el directorio ya existe
if isfolder(outputDir)
    error(['Output directory "', outputDir, '" exists. Make sure to rename it or delete it before proceeding.']);
else
    mkdir(outputDir); % Crear la carpeta si no existe
    disp(['Creating output directory: ', outputDir]);
end



% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%               LET'S BEGING
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
initCobraToolbox(false);

fprintf('== Initializing EXP2FLUX analysis ==\n\n')
fprintf('= Checking environment...\n\n')

% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%               LOAD MODEL
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
model = readCbModel(cobraModel);

% Check if the model was successfully loaded
if isempty(model)
    error('Error: Model "%s" failed to load. Stopping...', cobraModel);
end

fprintf('= Model "%s" loaded successfully. Proceeding...\n', cobraModel);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%               LOAD EXPRESSION DATA
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%Actually read csv file
dataTable = readtable(csvFile);
fprintf('= CSV file "%s" correctly loaded.\n', csvFile);

%TODO since this script assumes that expression file holds two columns:
%Gene and Expression_Value, makes sense to check for it before proceeding?
%Even bettertake the first colum and change it for Gene and
%Expression_Value if it has a different name.




% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%            Check if exp2flux is installed
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
scripts = {'plotDifferences.m', 'fluxdifferences.m', 'exp2flux.m'};

% Loop through each script and check if it exists
for i = 1:length(scripts)
    if exist(scripts{i}, 'file') == 0
        error('Error: %s is missing from the MATLAB path. Exiting MATLAB.', scripts{i});
    end
end

%Uncomment the following line just for testing
fprintf('= Exp2flux is correctly installed.\n');




% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%            Check if cobratoolbox is installed
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cobra_functions = {'initCobraToolbox', 'changeCobraSolver', 'optimizeCbModel'};
for i = 1:length(cobra_functions)
    if exist(cobra_functions{i}, 'file') == 0
        fprintf('Error: COBRA Toolbox function "%s" is missing. Exiting MATLAB.\n', cobra_functions{i});
        quit force; % Exit MATLAB
    end
end

fprintf('= COBRA Toolbox is correctly installed.\n');



% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%            Check model correctness
%
% Before running exp2flux we need to make sure model
% necessary characteristics are correct.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% TODO: check not for NaN or Inf but for any empty or
% non numeric values.
disp('Checking model boundaries...');
if any(isnan(model.lb)) || any(isinf(model.lb))
    warning('Some values are not valid in lb. Making necessary corrections.');
    model.lb(isnan(model.lb) | isinf(model.lb)) = -1000;
end
if any(isnan(model.ub)) || any(isinf(model.ub))
    warning('Some values are not valid in ub. Making necessary corrections.');
    model.ub(isnan(model.ub) | isinf(model.ub)) = 1000;
end

disp('Checking C vector correctness...');
if length(model.c) ~= length(model.rxns)
    error('Vector C length differs from reactions vector.');
end

%TODO: display possible objective functions, wait for user input and
% assign.Maybe check for something called "biomass" or so.
if all(model.c == 0)
    warning('Objective function not found. Assigning one.');
    model.c(find(model.ub > 0, 1)) = 1; % Asignar un objetivo genérico
end

% Agregar los campos faltantes para compatibilidad con exp2flux
model.lb_original = model.lb;
model.ub_original = model.ub;
model.rev = (model.lb < 0 & model.ub > 0); % Campo 'rev'


% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%            Prepare expression data
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Extract genes and expression values
disp('Extracting genes names and expression values')

%Remember that dataTable holds the expression file. It has two columns
%Gene and Expression_Value. 
genesFile = dataTable.Gene; 
expressionValues = dataTable.Expression_Value;

% Verificar la coincidencia de genes entre el modelo y el archivo
% CommonGenes: an array with the genes in model and expression file.
% idxModel: índices del modelo en commonGenes
% idxData: índices del archivo de genes en commonGenes
%
% Intersection cases:
% model.genes > genesFile (M>E):
% model holds more genes than expressed ones.
%
% model.genes < genesFile (M<E):
% model has less genes than expression file, this suggests an incomplete model.
%
% model.genes == genesFile (M=E)
% suggests that allgenes in model where expressed.

[commonGenes, idxModel, idxData] = intersect(model.genes, genesFile);

% Contar coincidencias y genes no encontrados
numCommon = length(commonGenes);

% IMPORTANT, If there are none coincidences between expression file genes
% and model genes, makes no sense to continue.
if numCommon == 0
    error('No common genes between model and expression file. Quiting.')
end


%If (M>E)numMissing will be a subset of model genes.
%If (M<E)numMissing will be a subset of expressed genes.
%If (M=E)numMissing will be "0" 
%
% setdiff checks for elements in A not present in B. So it assumes
% that always B is a subset of A. BUt we can have the case where A
% could be a subset of B, so it is necessary to check for both cases.
%
if length(genesFile) > length(model.genes) 

    numMissing = length(genesFile) - numCommon;
    missingGenes = setdiff(genesFile, model.genes);%(M<E)
    msg="Genes in expression file not present in the model."

else

    numMissing = length(model.genes) - numCommon;
    missingGenes = setdiff(model.genes, genesFile);%(M>E)
    msg="Genes in model not present in the expression file."

end
%

% Mostrar resultados de la comparación
disp(' ')
disp(msg);
disp(['Genes in model:             ', num2str(length(model.genes))]);
disp(['Genes in expression file:   ', num2str(length(genesFile))]);
disp(['Common genes:               ', num2str(numCommon)]);
disp(['Missing genes:              ', num2str(numMissing)]);

%save missing genes in a file
if numMissing > 0
    disp('')
    
    filePath = fullfile(outputDir,'missingGenes.txt');
    fprintf('Saving missing genes list to %s',filePath);
    T = table(missingGenes);
    writetable(T,filePath,"WriteVariableNames", false);

    disp('')
end

% VENN DIAGRAM
% TODO: find a better way to create this UGLY venn diagram
% Define the two gene sets
genesInModel = model.genes;  % Genes present in the metabolic model
genesInExpressionFile = genesFile;  % Genes from the expression dataset

% Compute common and missing genes
commonGenes = intersect(genesInModel, genesInExpressionFile);
missingGenes = setdiff(genesInModel, genesInExpressionFile);
onlyInExpression = setdiff(genesInExpressionFile, genesInModel);

% Get counts for Venn diagram
numModelGenes = length(genesInModel);
numExpressionGenes = length(genesInExpressionFile);
numCommon = length(commonGenes);
numMissing = length(missingGenes);
numOnlyExpression = length(onlyInExpression);

%plot the ven diagram
figure;
hold on;
r = 1; % Circle radius
theta = linspace(0, 2*pi, 100);

% Define circle positions
x1 = -0.5; y1 = 0;
x2 = 0.5; y2 = 0;

% Plot circles
fill(x1 + r*cos(theta), y1 + r*sin(theta), 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none'); % Model genes
fill(x2 + r*cos(theta), y2 + r*sin(theta), 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none'); % Expression file genes

% Add text labels
text(x1, y1, sprintf('Model Genes\n(%d)', numModelGenes), 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'w');
text(x2, y2, sprintf('Expression Genes\n(%d)', numExpressionGenes), 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'w');
text(0, 0, sprintf('Common\n(%d)', numCommon), 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');

hold off;
title('Genes in Model vs. Expression File');
axis equal;
axis off;

filePath = fullfile(outputDir, 'venn_diagram_shared_genes.png');
saveFigure(filePath);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%           HANDLING MISSING EXPRESSION VALUES
%
% TODO:This code here could be better managed in the exp2flux.m
% file itself (line 131) where the specific method is implemented.
% Somehow it makes no sense to have it here if this methods are hadled 
% by the exp2flux file.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% Create a vector with the same length of genes in model and then populate
% with the expression data present in the expressionValues array.
% If (M>E): there will be some zeros in the expressionProcessed array.
% If (M<E) || (M=E): there will be none zeros in expressionProcessed array.
% So the following is not necessary.

% Create a vector of zeros with same model.genes length
expressionProcessed = zeros(length(model.genes), 1);
% map expression values to correct positions in expressionProcessed. 
expressionProcessed(idxModel) = expressionValues(idxData);

% Guardar los datos procesados en la estructura de expresión
expression.genes = model.genes; % Mantener el orden del modelo
expression.data = expressionProcessed; % Asignar valores de expresión en el mismo orden

% Handle the zeros when (M>E). This is the only case when we will find
% zeros in the expressionProcessed array. Meaning that there will be some
% genes in the model that will be assigned the expression value of "0".
% Clearly this is not correct, so we have to assign a plausible value to
% this genes in the model that are not represented in the model.
% Here we have some options to fix those zero expressions:
% 1. Assign them the min expression value.
% 2. Assign them the max expression value.
% 3. Assign them the mean expression value.
% 4. Assign them an adjusted a distribution model value.
missingIdx = (expressionProcessed == 0);
if any(missingIdx)
    
    %inferExpValue = min(expressionProcessed(expressionProcessed > 0)); % get min value different from 0
    %inferExpValue = max(expressionProcessed(expressionProcessed > 0)); % get min value different from 0
    %inferExpValue = mean(expressionProcessed(expressionProcessed ~= 0))
    
    inferExpValue = median(expressionProcessed(expressionProcessed ~= 0))

   %Need to learn how to implement Gamma. THis code is not working.
   %inferExpValue = fitdist(expressionProcessed(expressionProcessed ~= 0), 'Gamma')
    


    %hereby we assign
    expression.data(missingIdx) = inferExpValue; % Asignar valor mínimo a genes faltantes
    disp(['It was assigned the value: (', num2str(inferExpValue), ') to all ', num2str(sum(missingIdx)), ' genes without expression data.']);
end


% Save new expression values data in file
filePath = fullfile(outputDir,'inputedExpression.csv');
fprintf('Saving inputed expression data to %s',filePath);
T = table(expression.genes, expression.data, 'VariableNames',{'Gene','New_Expression_value'});
writetable(T,filePath,"WriteVariableNames", true, "Delimiter", "\t");
disp('')


% save original boundaries
initial_lb = model.lb;
initial_ub = model.ub;


% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%       APPLY EXPRESSION RESTRICTIONS WITH EXP2FLUX
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
disp('Applying exp2flux method...');


% INPUTS:
%   model    - COBRA model structure (must include rules and bounds).
%   expression - A structure with two fields:
%                * genes: Cell array of gene IDs corresponding to expression data.
%                * data: Matrix where rows represent genes and columns represent samples.
%   organism - (Optional) KEGG organism identifier for pathway data.
%   typeID   - (Optional) Type of gene IDs in GPR rules ('entrez' or 'kegg').
%   missing  - (Optional) Method for handling missing gene expression values. %             
%                         Options: 'min', '1q', 'mean', 'median', '3q', 'max'. Default: 'mean'.
%   scale    - (Optional) Boolean indicating whether to scale expression data to
%              a maximum value of 1000. Default: false.

% Here missing='mean' has no effect since expression has no missing values.
updatedModel = exp2flux(model, expression, [], [], 'mean', true);



% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%            Compare boundaries pre and post exp2flux
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Save pre and post exp2flux boundaries
boundComparison = table(model.rxns, initial_lb, updatedModel.lb, initial_ub, updatedModel.ub, ...
    'VariableNames', {'Reaction', 'LB_model', 'LB_exp2flux', 'UB_model', 'UB_exp2flux'});


filePath = fullfile(outputDir,'boundaries.txt');
fprintf('Saving old and new model boundaries to %s',filePath);

writetable(boundComparison,filePath,"WriteVariableNames", true, "Delimiter","\t");
disp('')


% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%           PERFORM OPTIMIZATIONS
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -


disp('Optimizing original model...');
solution_original = optimizeCbModel(model);


%disp('Metabolic phenotype:');
%disp(solution_original.x);

disp('Optimizing exp2flux modified model...');
solution_updated = optimizeCbModel(updatedModel);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%            DISPLAY OBJECTIVE FUNCTIONS
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
disp('')
disp('--- Optimization results ---');
disp(['Original model:   ', num2str(solution_original.f)]);
disp(['Restricted model: ', num2str(solution_updated.f)]);
disp('')

% Also save solutions to file
% Format the optimization results as a string
resultsText = sprintf(['\n--- Optimization results ---\n', ...
    'Original model:   %f\n', ...
    'Restricted model: %f\n\n'], ...
    solution_original.f, solution_updated.f);

filePath = fullfile(outputDir, 'optimization_results.txt');
saveToFile(filePath, resultsText);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%            COMPARE PHENOTYPES
%
% Hereby we compare the differences between metabolic 
% phenotype in original model (A) vs. exp2flux model (B)
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%get the absolute difference between both phenotypes.
diff_abs = abs(solution_original.x - solution_updated.x);

% RATIO. Let's see how much increses or decreases B related to A.
% A/B=2, for instance, means that A doubles B.
ratioAB = solution_original.x ./ solution_updated.x;
ratioBA = solution_updated.x ./ solution_original.x;

%RELATIVE DIFFERENCE. |A-B|/max(A,B)
%measures the proportional change in flux, related (normalized) to max flux
%value. So we can comapre fluxex regarding the scale.
%Values close to 1 mean a drastic flux change.
%>Values close to 0 mean low variation.
diffRel = abs(solution_original.x - solution_updated.x) ./ max(abs(solution_original.x), abs(solution_updated.x));


%save table with this data
fluxDifference = table(model.rxns, solution_original.x, solution_updated.x,diff_abs,ratioAB,ratioBA, diffRel, ...
    'VariableNames', {'Reaction', 'Model_fluxes', 'exp2flux_fluxes', 'ABS (A-B)','Ratio_AB','Ratio_BA','Relative_difference'});
%disp(fluxDifference);

filePath = fullfile(outputDir,'phenotypeComparisson.csv');
writetable(fluxDifference,filePath,"WriteVariableNames", true, "Delimiter","\t");
disp('')


%                           OTHER GLOBAL METRICS
%MEAN SQUARE ERROR (MSE)
MSE = mean((solution_original.x - solution_updated.x).^2);
filePath = fullfile(outputDir,'MSE.txt');
saveToFile(filePath,MSE);

%ROOT MEAN SQUARE. A high RMSE could mean that there are important changes
% in various reactions.
RMSE = sqrt(MSE);
filePath = fullfile(outputDir,'RMSE.txt');
saveToFile(filePath,RMSE);

%EUCLIDIAN DISTANCE. Measures the total change in the metabolic profile.
% The higher the value the higher the difference between B (restricted model) and
% A (original model).
ED = norm(solution_original.x - solution_updated.x);
filePath = fullfile(outputDir,'euclidean.txt');
saveToFile(filePath,ED);

%PEARSON CORRELATION. Measures the general tendency between conditions.
% r=1 means no change at all.
pearson = corr(solution_original.x, solution_updated.x);
filePath = fullfile(outputDir,'pearson.txt');
saveToFile(filePath,pearson);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
%                   PLOTTING SECTION
%   
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Graficar comparación de flujos
figure;
bar([solution_original.x, solution_updated.x]);
legend({'Original model', 'exp2flux model'});
xlabel('Reaction index');
ylabel('Flux');
title('Non restricted vs. restricted model');

filePath = fullfile(outputDir, 'restricted_vs_nonrestricted.png');
saveFigure(filePath);



% Helps to identify if there are global shifts in flux distributions 
% (e.g., wider/narrower spread, higher/lower mean).
figure;
histogram(solution_original.x, 'FaceColor', 'b', 'FaceAlpha', 0.5);
hold on;
histogram(solution_updated.x, 'FaceColor', 'r', 'FaceAlpha', 0.5);
hold off;
legend({'Original model', 'exp2flux model'});
xlabel('Flux Value');
ylabel('Frequency');
title('Flux Distribution');

filePath = fullfile(outputDir, 'flux_distributions.png');
saveFigure(filePath);



% Shows correlation between fluxes in both models. Deviations from 
% the diagonal indicate significant metabolic shifts.
figure;
scatter(solution_original.x, solution_updated.x, 'filled');
hold on;
plot([min(solution_original.x), max(solution_original.x)], ...
     [min(solution_original.x), max(solution_original.x)], 'k--'); % Diagonal reference line
hold off;
xlabel('Original Model Flux');
ylabel('exp2flux Model Flux');
title('Flux Correlation');
grid on;

filePath = fullfile(outputDir, 'flux_correlation.png');
saveFigure(filePath);


% Identify reactions with large changes in fluxes (either positive or negative).
% Highlights which reactions are significantly altered, helping in 
% pathway-level interpretations.
fluxChanges = solution_updated.x - solution_original.x;
logFC = log2(abs(fluxChanges) + 1); % Log2 fold change for better scaling

figure;
scatter(1:length(fluxChanges), logFC, 15, fluxChanges, 'filled'); % Color by change magnitude
colorbar;
xlabel('Reaction Index');
ylabel('Log2 Flux Change');
title('Volcano Plot of Flux Changes');

filePath = fullfile(outputDir, 'volcano_plot.png');
saveFigure(filePath);

% Identify reactions with the largest absolute flux changes.
% Depicts the most significant metabolic shifts, helping identify key affected pathways.
fluxChanges = abs(solution_updated.x - solution_original.x);
[sortedChanges, idx] = sort(fluxChanges, 'descend'); % Sort by highest change

topN = 10; % Number of top reactions to plot
topReactions = model.rxns(idx(1:topN));
topChanges = sortedChanges(1:topN);

figure;
barh(topChanges, 'FaceColor', 'c');
set(gca, 'YTickLabel', topReactions, 'YDir', 'reverse'); % Set reaction names
xlabel('Absolute Flux Change');
title('Top 10 Most Changed Reactions');

filePath = fullfile(outputDir, 'top10_rxn_changes.png');
saveFigure(filePath);

%Visualize how entire pathways change between models. 
% Provides an overview of where major flux shifts occur.
fluxDiffMatrix = reshape(solution_updated.x - solution_original.x, [], 1); % Ensure it's a column
figure;
imagesc(fluxDiffMatrix);
colormap jet;
colorbar;
xlabel('Reaction Index');
ylabel('Flux Change');
title('Heatmap of Flux Differences');

filePath = fullfile(outputDir, 'heatMap.png');
saveFigure(filePath);



% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
%                  PATHWAY ANALYSIS
%
% This analysis needs that the "subSystems" field is present.
%   
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if isfield(model, 'subSystems')
    pathways = model.subSystems; % Cell array of pathways

    % Obtain absolute flux changes
    fluxChanges = abs(solution_updated.x - solution_original.x);

    %%Aggregate fluxchanges by pathways
    uniquePathways = unique(pathways);
    numPathways = length(uniquePathways);
    pathwayFluxChanges = zeros(numPathways, 1);

    % Sum flux changes for each pathway
    for i = 1:numPathways
        pathwayIdx = strcmp(pathways, uniquePathways{i}); % Find reactions in this pathway
        pathwayFluxChanges(i) = sum(fluxChanges(pathwayIdx)); % Sum flux changes
    end
    
    %%Visualize most affected pathways
    % Sort pathways by flux change
    [sortedFlux, idx] = sort(pathwayFluxChanges, 'descend');
    sortedPathways = uniquePathways(idx);
    
    % Plot top 10 affected pathways
    topN = 10;
    figure;
    barh(sortedFlux(1:topN));
    set(gca, 'YTickLabel', sortedPathways(1:topN), 'YDir', 'reverse');
    xlabel('Total Flux Change');
    title('Top 10 Most Affected Metabolic Pathways');

    filePath = fullfile(outputDir, 'top10_affected_pathways.png');
    saveFigure(filePath);

else
    warning('No pathway information available in model.subSystems');
end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
%             TOPOLOGICAL ANALYSIS (experimental)
%
%    TODO: add support for cytoscape
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Construct a reaction-reaction network based on stoichiometry matrix
% Basically this digraph will be used from now on for all reamining
% topological analysis.
%
% The idea to filter out cofactors is because they appear in many reactions
% and can skew the connectivity of the network, since we are going to
% create the graph based in shared metabolites (and connectivity
% directionality) so sharing atp or water is very common but not meaningful.
%
% TODO: figure out how to filter cofactors with non standard names such as
% those in models created in from other databases.
S = model.S(~ismember(model.mets, {'atp[c]', 'nad[c]', 'h2o[c]', 'pi[c]'}), :);


%Consider connection where the metabolite is produced and then consumed.
[reaction1_idx, reaction2_idx] = find((S' > 0) * (S < 0));

reactionNames = model.rxns;
reaction1 = reactionNames(reaction1_idx);
reaction2 = reactionNames(reaction2_idx);

%Maybe the flux difference could be a ratio or something like that? Instead
%of the absolute value.
%In this case high values represent that the change between models for that
%particular reaction was strong. Usually towards red color.
fluxDiffs = abs(solution_updated.x - solution_original.x);

% Create directed graph (nodes = reactions, edges = shared metabolites)
G = digraph(reaction1, reaction2, fluxDiffs(reaction1_idx));


%                  --- FLUX CHANGE NETWORK ---
%
% Plot flux change network (highlighting most altered reactions)
% Helps visualize connections between reactions.
% Pathways with high flux changes will appear as densely connected regions.

figure;
p = plot(G, 'Layout', 'circle', 'EdgeCData', G.Edges.Weight);
colormap jet;
colorbar;
title('Flux Change Network: Reactions Linked by Metabolites');
filePath = fullfile(outputDir, 'graph_flux_change.png');
saveFigure(filePath);



%Top ten of most affected connections after exp2flux
sortedEdges = sortrows(G.Edges, 'Weight', 'descend'); % sort
top10Edges = sortedEdges(1:10, :);

filePath = fullfile(outputDir,'top10_affected_connections.txt');
writetable(top10Edges, filePath, 'Delimiter', '\t', 'WriteVariableNames', true);





%    -- COMMUNITY (CONNECTED SUBGRAPHS) DETECTION --
% Use modularity analysis to detect clusters of 
% reactions that behave similarly.
% Identifies clusters of co-regulated metabolic reactions.
% Useful for detecting functional metabolic modules.

% Compute modularity-based community detection
% bins is a vector where each node gets assigned a community index.
%
% Example with two communities (1 and 2):
% Glycolisis  --> 1
% TCA_Cycle   --> 1
% Lipid_synth --> 2
bins = conncomp(G, 'Type', 'weak'); % Find communities
numCommunities = max(bins); %Check how many communities we have.
colors = lines(numCommunities); % 'lines' provides distinct colors

% Assign colors based on communities
figure;
p = plot(G, 'Layout', 'circle');
%highlight(p, find(bins == 1), 'NodeColor', 'r'); % Example: first community
%title(sprintf('Metabolic Reaction Communities (n=%d)', numCommunities));
nodeColors = colors(bins, :); % Assign each node a color based on community
p.NodeColor = nodeColors; % Apply the colors to the nodes

% Create a legend for communities
hold on;
legendHandles = gobjects(numCommunities, 1);
for i = 1:numCommunities
    % Create dummy scatter plots to act as legend items
    legendHandles(i) = scatter(nan, nan, 100, colors(i, :), 'filled'); 
end
hold off;

% Add the legend with community labels
legend(legendHandles, compose('Community %d', 1:numCommunities), 'Location', 'bestoutside');

% Set title
title(['Metabolic Network Communities (', num2str(numCommunities), ' found)']);


filePath = fullfile(outputDir, 'graph_community.pdf');
saveFigure(filePath);


%                   -- FLUX FLOW PROPAGATION --
%Use graph centrality measures to identify:
%Reactions acting as metabolic bottlenecks (high betweenness centrality).
%Highly flux-altered reactions affecting the entire network.

outdegreeCentrality = centrality(G, 'outdegree'); % Outgoing connections per reaction
indegreeCentrality = centrality(G, 'indegree'); % Incoming connections per reaction
betweennessCentrality = centrality(G, 'betweenness'); % Nodes that act as bottlenecks

% Highlight betweenness centrality
figure;
p = plot(G, 'Layout', 'circle', 'NodeCData', betweennessCentrality);
colormap hot;
colorbar;
title('Betweenness Centrality of Metabolic Reactions');

filePath = fullfile(outputDir, 'graph_fluxPropagation_betweenness.pdf');
saveFigure(filePath);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%               CLOSING UP
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
disp('')
disp('All done. quitting exp2flux analysis. Good bye.')
disp('')


% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%              GENERAL FUNCTIONS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function saveToFile(filePath, data)
    % Abrir el archivo en modo escritura
    fileID = fopen(filePath, 'w');
    
    % Verificar si el archivo se abrió correctamente
    if fileID == -1
        error(['No se pudo abrir el archivo: ', filePath]);
    end

    % Verificar el tipo de datos y escribirlos en el archivo
    if isnumeric(data)
        % Si es un solo número, escribirlo en formato flotante
        if isscalar(data)
            fprintf(fileID, '%f\n', data);
        % Si es un vector o matriz, escribirlo con formato adecuado
        elseif isvector(data)
            fprintf(fileID, '%f ', data);
            fprintf(fileID, '\n');
        else
            for i = 1:size(data, 1)
                fprintf(fileID, '%f ', data(i, :));
                fprintf(fileID, '\n');
            end
        end
    elseif ischar(data) || isstring(data)
        % Guardar cadenas de texto
        fprintf(fileID, '%s\n', data);
    elseif iscell(data)
        % Guardar celdas de texto o mixtas
        for i = 1:length(data)
            fprintf(fileID, '%s\n', data{i});
        end
    else
        % Tipo de dato no soportado
        fclose(fileID);
        error('El tipo de datos proporcionado no es compatible con esta función.');
    end

    % Cerrar el archivo
    fclose(fileID);

    % Mensaje de confirmación
    %disp(['✅ Archivo guardado en: ', filePath]);
end


%Function to save figures
function saveFigure(filePath)
    % Verificar si hay una figura activa
    if isempty(findall(0, 'Type', 'Figure'))
        error('No hay ninguna figura abierta para guardar.');
    end
    
    % Obtener la extensión del archivo
    [~, ~, ext] = fileparts(filePath);
    
    % Verificar formato compatible
    validFormats = {'.png', '.jpg', '.jpeg', '.pdf', '.eps', '.tiff', '.bmp', '.fig'};
    if ~ismember(lower(ext), validFormats)
        error('Formato de archivo no válido. Usa: .png, .jpg, .pdf, .eps, .tiff, .bmp o .fig');
    end

    % Guardar la figura en el formato deseado
    saveas(gcf, filePath);

    % Mensaje de confirmación
    %disp(['✅ Figura guardada en: ', filePath]);
end

