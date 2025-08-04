function fChange = fluxDifferences(model1, model2, foldReport)
    % FLUXDIFFERENCES Calculate fold change of fluxes between two models.
    %
    % This function computes the fold change for the fluxes of two given
    % metabolic models. The fold change is calculated as:
    %     (fluxModel2 / fluxModel1) - 1
    %
    % INPUTS:
    %   model1     - COBRA model structure for the first scenario.
    %   model2     - COBRA model structure for the second scenario. Must
    %                have the same reactions as model1 (same IDs and number).
    %   foldReport - Threshold for reporting fold changes. Only reactions
    %                with an absolute fold change >= foldReport are included
    %                in the output. Default: 2.
    %
    % OUTPUT:
    %   fChange    - Matrix containing flux values and fold changes for
    %                reactions meeting the threshold:
    %                Columns: {"fluxModel1", "fluxModel2", "foldChange"}
    %
    % DEPENDENCIES:
    %   COBRA Toolbox functions (optimizeCbModel, etc.)
    %
    %
    % VERSION:
    %   1.0 - January 2025

    if nargin < 3
        foldReport = 2; % Default threshold
    end

    % Get flux distributions for both models
    solution1 = optimizeCbModel(model1);
    solution2 = optimizeCbModel(model2);

    % Ensure both models solved successfully
    if isempty(solution1.x) || isempty(solution2.x)
        error('One or both models did not solve successfully.');
    end

    f_m1 = solution1.x; % Flux distribution for model1
    f_m2 = solution2.x; % Flux distribution for model2

    % Identify common reactions
    [commonReactions, idx1, idx2] = intersect(model1.rxns, model2.rxns);

    % Initialize results matrix
    fChange = zeros(length(commonReactions), 3);
    flux1 = f_m1(idx1); % Fluxes from model1 for common reactions
    flux2 = f_m2(idx2); % Fluxes from model2 for common reactions

    % Calculate fold changes
    fold = (flux2 - flux1) ./ abs(flux1); % Fold change formula

    % Handle special cases where flux1 is zero
    zeroFluxIdx = (flux1 == 0);
    fold(zeroFluxIdx) = flux2(zeroFluxIdx);

    % Replace NaN and Inf values with 0
    fold(isnan(fold)) = 0;
    fold(isinf(fold)) = 0;

    % Store results in fChange matrix
    fChange(:, 1) = flux1; % Flux values for model1
    fChange(:, 2) = flux2; % Flux values for model2
    fChange(:, 3) = fold;  % Fold changes

    % Filter results by the foldReport threshold
    significant = abs(fChange(:, 3)) >= foldReport;
    fChange = fChange(significant, :);
end
