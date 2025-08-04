function plotDifferences(model1, model2, varargin)
    % PLOTDIFFERENCES Plot the fold change of fluxes between two models as a bipartite graph.
    %
    % This function calculates the fold change for the fluxes of two given metabolic
    % models and visualizes it as a bipartite graph. Vertex size is proportional
    % to the absolute fold change. Vertex color indicates the direction of change:
    % green for positive and red for negative fold changes.
    %
    % INPUTS:
    %   model1 - COBRA model structure for the first scenario.
    %   model2 - COBRA model structure for the second scenario. Must have the same
    %            reactions (reaction identifiers) as model1 but with different restrictions.
    %   varargin - Additional arguments for customizing the plot.
    %
    % OUTPUT:
    %   A bipartite graph is displayed, showing flux changes between the two models.
    %
    % DEPENDENCIES:
    %   Requires a MATLAB graph library or manual implementation for plotting bipartite graphs.
    %
    % AUTHORS:
    %   Daniel Camilo Osorio <dcosorioh@unal.edu.co>
    %   Kelly Botero <kjboteroo@unal.edu.co>
    %
    % VERSION:
    %   1.0 - January 2025

    % Get stoichiometric matrix and set up reaction-matrix relationships
    S = model1.S';
    [numReactions, numMetabolites] = size(S);
    
    % Normalize stoichiometric matrix for visualization
    S(S < 0) = -1;
    S(S > 0) = 1;
    S(:, model1.rev) = 2;

    % Calculate flux differences using `fluxDifferences`
    fD = fluxDifferences(model1, model2, 0);
    nonZeroFluxIdx = find(model1.c ~= 0);
    S = S(nonZeroFluxIdx, :); % Filter reactions with non-zero flux
    S = S(:, sum(S ~= 0, 1) ~= 0); % Remove columns without contributions

    % Assign types to nodes (1 = reactions, 0 = metabolites)
    reactionNames = model1.rxns(nonZeroFluxIdx);
    metaboliteNames = model1.mets(any(S ~= 0, 1));
    allNodes = [reactionNames; metaboliteNames];
    types = [ones(length(reactionNames), 1); zeros(length(metaboliteNames), 1)];

    % Create edges for the bipartite graph
    edges = [];
    for i = 1:length(reactionNames)
        reactIdx = i;
        connectedMets = find(S(reactIdx, :) ~= 0);
        for j = connectedMets
            if S(reactIdx, j) == 2 % Reversible reactions
                edges = [edges; [j + length(reactionNames), reactIdx]]; % Metabolite -> Reaction
                edges = [edges; [reactIdx, j + length(reactionNames)]]; % Reaction -> Metabolite
            elseif S(reactIdx, j) == 1 % Positive flux
                edges = [edges; [reactIdx, j + length(reactionNames)]];
            else % Negative flux
                edges = [edges; [j + length(reactionNames), reactIdx]];
            end
        end
    end

    % Create the graph object
    G = digraph(edges(:, 1), edges(:, 2), [], allNodes);

    % Assign vertex properties: colors and sizes
    nodeColors = repmat([0.8, 0.8, 0.8], length(allNodes), 1); % Default: gray
    nodeSizes = ones(length(allNodes), 1) * 10; % Default size for metabolites

    for i = 1:length(reactionNames)
        rxn = reactionNames{i};
        if ismember(rxn, fD.Properties.RowNames)
            foldChange = fD{rxn, 'foldChange'};
            if foldChange > 0
                nodeColors(i, :) = [0, 1, 0]; % Green for positive
            elseif foldChange < 0
                nodeColors(i, :) = [1, 0, 0]; % Red for negative
            end
            nodeSizes(i) = abs(foldChange) / max(abs(fD.foldChange)) * 20;
        end
    end

    % Plot the graph
    figure;
    plot(G, 'NodeColor', nodeColors, 'MarkerSize', nodeSizes, ...
         'NodeLabel', allNodes, 'ArrowSize', 8);
    title('Fold Change of Fluxes Between Models');
end
