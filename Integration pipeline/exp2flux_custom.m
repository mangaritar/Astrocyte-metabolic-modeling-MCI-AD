function modifiedModel = exp2flux_custom(model, expression)
    % Copiar el modelo original
    modifiedModel = model;

    % Copiar los límites originales
    lb = model.lb;
    ub = model.ub;

    % Aplicar los valores de expresión a los límites de flujo
    modifiedModel.lb = -1 * expression;
    modifiedModel.ub = expression;

    % Identificar reacciones irreversibles (donde lb es >= 0 originalmente)
    irreversible_rxns = lb >= 0;
    modifiedModel.lb(irreversible_rxns) = 0;

    % Mantener límites originales en reacciones de intercambio
    exchangeRxns = findExchangeReactions(model);
    modifiedModel.lb(exchangeRxns) = lb(exchangeRxns);
    modifiedModel.ub(exchangeRxns) = ub(exchangeRxns);
end
