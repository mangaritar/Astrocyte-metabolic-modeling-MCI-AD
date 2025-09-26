function expressionStruct = cargarExpresion(exprFile, mappingFile)
    % Cargar expresión y hacer el mapeo SYMBOL → ENTREZ → ENTREZ.1
    expr = readtable(exprFile);
    map = readtable(mappingFile);
    
    % Renombrar columna y estandarizar
    expr.Properties.VariableNames{1} = 'gene';
    expr.gene = upper(string(expr.gene));
    map.SYMBOL = upper(string(map.SYMBOL));
    map.ENTREZID = string(map.ENTREZID);

    % Mapear genes
    [found, idx] = ismember(expr.gene, map.SYMBOL);
    entrez = map.ENTREZID(idx(found));

    % Construir estructura
    expressionStruct.genes = entrez + ".1";  % tipo 8639.1
    expressionStruct.data  = mean(expr{found, 2:end}, 2, 'omitnan');
end
