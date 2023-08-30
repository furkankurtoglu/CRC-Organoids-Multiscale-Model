
function [Model, stat, MinimizedFlux, v, r, p, q] = upFBA_pipeline(Model, initvalue, met_IDs, foldchange_means, foldchange_sds, SPECIES_BOUND)

    %% estimate metabolite rates of change
    changeSlopes = zeros(length(met_IDs), 1);
    changeIntervals = zeros(length(met_IDs), 1);
    SPECIES_HARD_CONSTRAINT = [];
    SPECIES_COUPLED_CONSTRAINT = [];
    for i = 1:length(met_IDs)
        tmp1 = initvalue(i) * (foldchange_means{i, 1} - 1) / 24.0;
        changeSlopes(i, 1) = tmp1;
        changeIntervals(i, 1) = 0;
        tmp_met_ID = split(met_IDs{i}, ";");
        if numel(tmp_met_ID) == 1
            SPECIES_HARD_CONSTRAINT = [SPECIES_HARD_CONSTRAINT; met_IDs(i) tmp1];
        else
            SPECIES_COUPLED_CONSTRAINT = [SPECIES_COUPLED_CONSTRAINT; met_IDs(i) tmp1];
        end
    end
    
    measurable_mets = [];
    for i = 1:length(met_IDs)
        measurable_mets = [measurable_mets; split(met_IDs{i}, "; ")];
    end
    
    %% change species bounds
    Model = addSpeciesBounds_new(Model, SPECIES_BOUND);
    Model = addSpeciesHardConstraints_new(Model, SPECIES_HARD_CONSTRAINT);
    % optimize the model
    [MinimizedFlux, modelIrrev] = minimizeModelFlux(Model, 'min');
    if MinimizedFlux.stat == 1
        disp('Model is feasible. Nothing to do.');
        stat = MinimizedFlux.stat;
        v = zeros(size(Model.rxns));
        for i = 1:length(Model.rxns)
            if modelIrrev.match(i)
                v(i) = MinimizedFlux.x(i) - MinimizedFlux.x(modelIrrev.match(i));
            else
                v(i) = MinimizedFlux.x(i);
            end
        end
        r = zeros(size(Model.mets));
        p = zeros(size(Model.rxns));
        q = zeros(size(Model.rxns));
        writeCbModel(Model, 'format', 'SBML', 'fileName','TEEST.sbml')
    else
        disp('Model is infeasible');
        % Identify the exchange reactions and biomass reaction(s) heuristically
        Model = findSExRxnInd(Model,size(Model.S,1),0);
        relaxOption.internalRelax = 0;
        relaxOption.exchangeRelax = 0;
        relaxOption.steadyStateRelax = 1;
        relaxOption.epsilon = 10e-6;
        mets_to_exclude = false(size(Model.mets));
        idx = ismember(Model.mets, measurable_mets);
        idx2 = idx;
        mets_to_exclude(idx2) = true;
        mets_to_exclude(contains(Model.mets, "_greater")) = true;
        mets_to_exclude(contains(Model.mets, "_lower")) = true;
        relaxOption.excludedMetabolites = mets_to_exclude;
        solution = relaxedFBA(Model,relaxOption);
        [stat,v,r,p,q] = deal(solution.stat, solution.v, solution.r, solution.p, solution.q);
        if stat == 0
           v = zeros(size(Model.rxns));
           r = zeros(size(Model.mets));
           p = zeros(size(Model.rxns));
           q = zeros(size(Model.rxns));
        else
            b_nonzero_indices = find(r ~= 0);
            b_nonzero_values = r(b_nonzero_indices);
            for i=1:numel(b_nonzero_indices)
                Model.b(b_nonzero_indices(i)) = Model.b(b_nonzero_indices(i)) - b_nonzero_values(i);
            end    
            [MinimizedFlux, modelIrrev] = minimizeModelFlux(Model, 'min');
            stat = MinimizedFlux.stat;
            if stat == 1
                v = zeros(size(Model.rxns));
                for i = 1:length(Model.rxns)
                    if modelIrrev.match(i)
                        if endsWith(modelIrrev.rxns{i}, '_r')
                            v(i) = -(MinimizedFlux.x(i) - MinimizedFlux.x(modelIrrev.match(i)));
                        else
                            v(i) = MinimizedFlux.x(i) - MinimizedFlux.x(modelIrrev.match(i));
                        end
                    else
                        if endsWith(modelIrrev.rxns{i}, '_r')
                            v(i) = -MinimizedFlux.x(i);
                        else
                            v(i) = MinimizedFlux.x(i);
                        end
                    end
                end
            else
               v = zeros(size(Model.rxns));
               r = zeros(size(Model.mets));
               p = zeros(size(Model.rxns));
               q = zeros(size(Model.rxns));
            end
        end
    end
end
