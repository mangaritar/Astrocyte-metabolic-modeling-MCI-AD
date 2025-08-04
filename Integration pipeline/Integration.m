%% 1. Load model
astrocito = readCbModel('Astrocyte_Mendoza2022_1.xml');
astrocito.rules = mapIndexToGeneRule(astrocito);
astrocito.rev = astrocito.lb < 0;
astrocito.lb_original = astrocito.lb;
astrocito.ub_original = astrocito.ub;

%% 2. Load and map expression (helper function)
expressionStructControl  = cargarExpresion('GEP_astrocitos_Control.csv',  'symbols_to_entrez_control.csv');
expressionStructIncip    = cargarExpresion('GEP_astrocitos_Incipient.csv','symbols_to_entrez_control.csv');
expressionStructModerate = cargarExpresion('GEP_astrocitos_Moderate.csv', 'symbols_to_entrez_control.csv');
expressionStructSevere   = cargarExpresion('GEP_astrocitos_Severe.csv',   'symbols_to_entrez_control.csv');

%% 3. Unlock lactate reaction
idxLac = find(strcmp(astrocito.rxns, 'EX_lac_L[e]'));
astrocito.lb(idxLac) = 0;
astrocito.ub(idxLac) = 10;
astrocito.lb_original(idxLac) = 0;
astrocito.ub_original(idxLac) = 10;

%% 4. Integrate expression
modelo_control_astrocito  = exp2flux(astrocito, expressionStructControl,  [], [], 'mean', true);
modelo_incip_astrocito    = exp2flux(astrocito, expressionStructIncip,    [], [], 'min', true);
modelo_moderate_astrocito = exp2flux(astrocito, expressionStructModerate, [], [], 'min', true);
modelo_severe_astrocito   = exp2flux(astrocito, expressionStructSevere,   [], [], 'min', true);

modelo_control_astrocito  = changeObjective(modelo_control_astrocito,  'biomass_maintenance');
modelo_incip_astrocito    = changeObjective(modelo_incip_astrocito,    'biomass_maintenance');
modelo_moderate_astrocito = changeObjective(modelo_moderate_astrocito, 'biomass_maintenance');
modelo_severe_astrocito   = changeObjective(modelo_severe_astrocito,   'biomass_maintenance');

% Optimize
sol_control  = optimizeCbModel(modelo_control_astrocito);
sol_incip    = optimizeCbModel(modelo_incip_astrocito);
sol_moderate = optimizeCbModel(modelo_moderate_astrocito);
sol_severe   = optimizeCbModel(modelo_severe_astrocito);

%% 6. Display results
fprintf('CONTROL   : %.4f\n', sol_control.f);
fprintf('INCIPIENT : %.4f\n', sol_incip.f);
fprintf('MODERATE  : %.4f\n', sol_moderate.f);
fprintf('SEVERE    : %.4f\n', sol_severe.f);

%% 7. Export complete metabolic phenotype
fenotipo = table(astrocito.rxns, astrocito.rxnNames, ...
    sol_control.x, sol_incip.x, sol_moderate.x, sol_severe.x, ...
    'VariableNames', {'RxnID', 'Name', 'CONTROL', 'INCIPIENT', 'MODERATE', 'SEVERE'});

writetable(fenotipo, 'FenotipoMetabolico_Lactato.csv');

%%% Compare flux bounds between models
% Initialize table
tablaCambiosLimites = table(astrocito.rxns, astrocito.rxnNames, ...
    modelo_control.lb, modelo_incip.lb, modelo_moderate.lb, modelo_severe.lb, ...
    modelo_control.ub, modelo_incip.ub, modelo_moderate.ub, modelo_severe.ub, ...
    'VariableNames', {'RxnID', 'Name', ...
    'lb_control', 'lb_incip', 'lb_moderate', 'lb_severe', ...
    'ub_control', 'ub_incip', 'ub_moderate', 'ub_severe'});

% Calculate absolute differences for lower bounds
tablaCambiosLimites.diff_lb_incip    = abs(tablaCambiosLimites.lb_incip    - tablaCambiosLimites.lb_control);
tablaCambiosLimites.diff_lb_moderate = abs(tablaCambiosLimites.lb_moderate - tablaCambiosLimites.lb_control);
tablaCambiosLimites.diff_lb_severe   = abs(tablaCambiosLimites.lb_severe   - tablaCambiosLimites.lb_control);

% Calculate absolute differences for upper bounds
tablaCambiosLimites.diff_ub_incip    = abs(tablaCambiosLimites.ub_incip    - tablaCambiosLimites.ub_control);
tablaCambiosLimites.diff_ub_moderate = abs(tablaCambiosLimites.ub_moderate - tablaCambiosLimites.ub_control);
tablaCambiosLimites.diff_ub_severe   = abs(tablaCambiosLimites.ub_severe   - tablaCambiosLimites.ub_control);

% Compute max change in any bound
tablaCambiosLimites.maxCambio = max([ ...
    tablaCambiosLimites.diff_lb_incip, tablaCambiosLimites.diff_lb_moderate, tablaCambiosLimites.diff_lb_severe, ...
    tablaCambiosLimites.diff_ub_incip, tablaCambiosLimites.diff_ub_moderate, tablaCambiosLimites.diff_ub_severe ...
], [], 2);

% Sort by max change
tablaLimitesOrdenada = sortrows(tablaCambiosLimites, 'maxCambio', 'descend');

% Select top 10
top10Limites = tablaLimitesOrdenada(1:10, :);

% Display result
disp(top10Limites)
writetable(top10Limites, 'top10_cambios_limites_reacciones.csv');

% Save integrated models
save('modelo_ANGARITA_control.mat',  'modelo_control_astrocito');
save('modelo_ANGARITA_incip.mat',    'modelo_incip_astrocito');
save('modelo_ANGARITA_moderate.mat', 'modelo_moderate_astrocito');
save('modelo_ANGARITA_severe.mat',   'modelo_severe_astrocito');

