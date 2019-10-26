"""
Genetic Algorithm
"""

"""
generate_chromosome(nRow, nFactor)

Returns a random chromosome.
"""
function generate_chromosome(nRow::Int64, nFactor::Int64)
    new_chromosome::Array{Int64, 2} = 2 .* (rand(nRow, nFactor) .> 0.5) .- 1
    return(new_chromosome)
end


"""
generate_chromosome_mixed(nRow, nFactor)

Returns a random chromosome. This version of the function can handle mixed 2 and
3 level factors.
"""
function generate_chromosome_mixed(nRow::Int64, nFactor::Array{Int64, 1})
    new_chromosome::Array{Int64, 2} = zeros(nRow, sum(nFactor) )
    # Generate random values for the 2 level factors
    for iFactor = 1:nFactor[1]
        for iRow = 1:nRow
            new_chromosome[iRow, iFactor] = 2 * (rand()>0.5) - 1
        end
    end
    # Generate random values for the 3 level factors
    for iFactor = (nFactor[1]+1):sum(nFactor)
        for iRow = 1:nRow
            new_chromosome[iRow, iFactor] = rand(0:2)
        end
    end
    return(new_chromosome)
end


"""
generate_X_matrix_mixed!()

Fills out the entire X matrix for a chromosome.
"""
# More general function to generate a chromosome for mixed level factor
# NOTE, THE FIRST NFACTORS COLUMNS OF X SHOUD ONLY BE FILLED OUT ONCE,
# AFTER THAT WE JUST RE-UPDATE THE INTERACTIONS!!!!!
function generate_X_matrix_mixed!(X::Array{Int64, 2}, chromosome::Array{Int64, 2}, iModel::Int64, interactions::Array{Int64, 2},
    nRow::Int64, nFactor::Array{Int64, 1}, intercept::Int64, k::Int64, factor_start_cols::Array{Int64, 1},
    int_settings::Array{Int64, 1})

    if intercept == 1
        X[:,1] .= 1
    end

    # Preallocate int term (dont actually need to do this)
    intTerm::Int64 = 1

    for iCol = (intercept+1):(sum(nFactor) + intercept + k)
        for iRow = 1:nRow

            # 2 level factors
            if iCol <= (nFactor[1] + intercept)
                X[iRow, iCol] = copy(chromosome[iRow, iCol - intercept])
            # 3 level factors
            elseif iCol <= (nFactor[1] + nFactor[2] + intercept)
                X[iRow, factor_start_cols[iCol]] = copy(chromosome[iRow, iCol - intercept]) - 1
                if X[iRow, factor_start_cols[iCol]] == 0
                    X[iRow, factor_start_cols[iCol]+1] = -2
                else
                    X[iRow, factor_start_cols[iCol]+1] = 1
                end
            # Fill out the interactions
            else
                # add the interaction term TODO figure this out for 3level
                intTerm = 1
                three_flag = 0 # if a three level factor was added, have to subtract 1
                two_flag = 0
                first_linear_or_quad = 0
                for iFactor = 1:sum(nFactor)

                    # Linear Terms
                    #if interactions[iModel, iFactor] == 1
                    if interactions[int_settings[iCol - sum(nFactor) - intercept], iFactor] == 1
                        if first_linear_or_quad == 0
                            first_linear_or_quad = 1;
                        end
                        # interaction with a two-level factor(multiplication)
                        if iFactor <= nFactor[1]
                            intTerm = intTerm * copy(chromosome[iRow, iFactor])
                        # interaction with a three-level factor(addition)
                        else
                            three_flag = 1
                            intTerm = intTerm + copy(chromosome[iRow, iFactor])
                        end

                    # Quadratic Terms
                    #elseif interactions[iModel, iFactor] == 2
                    elseif interactions[int_settings[iCol - sum(nFactor) - intercept], iFactor] == 2
                        if first_linear_or_quad == 0
                            first_linear_or_quad = 2
                        end
                        two_flag = two_flag + 1
                        intTerm = intTerm + 2*copy(chromosome[iRow, iFactor])
                    end

                end
                # Since intTerm is originally set up to be multiplied,
                # need this to keep the addition version from being 1 unit too big
                if two_flag == 0
                    # Subtract 1 if a three level factor was included
                    # this is NOT the subtraction in the formula,
                    # this is to account for starting intTerm at 1 for multiplication
                    if three_flag == 1
                        intTerm = intTerm - 1
                    end
                    if intTerm <= 0
                        X[iRow, factor_start_cols[iCol]] = intTerm
                    else
                        X[iRow, factor_start_cols[iCol]] = mod(intTerm, 3)
                    end
                #TODO
                elseif two_flag == 1
                    # next step is based on whether 1 came first or 2 came first
                    # Linear - Quadratic
                    X[iRow, factor_start_cols[iCol]] = mod(intTerm-1, 3)
                    # Quadratic - Linear
                    if (first_linear_or_quad == 2)
                        #print("caseacheived \n")
                        if X[iRow, factor_start_cols[iCol]] == 2
                            X[iRow, factor_start_cols[iCol]] = -2
                        else
                            X[iRow, factor_start_cols[iCol]] = 1
                        end
                    end
                #TODO
                elseif two_flag == 2
                    intTerm = intTerm - 1
                    X[iRow, factor_start_cols[iCol]] = mod(intTerm-1, 3)
                    if (X[iRow, factor_start_cols[iCol]]) == 1
                        (X[iRow, factor_start_cols[iCol]]) = -2
                    else
                        (X[iRow, factor_start_cols[iCol]]) = 1
                    end
                end
            end # end of if elseif end
        end # end of loop over rows
    end # end of loop over columns
end



"""
crossover!(TODO)

crossover mutation.
"""
# this should create 2 new chromosomes
# these chromosomes replace the two worst chromosomes (replace1, replace2)
# TODO random starting crossover point, let it loop around
function crossover!(population::Array{Int64, 3}, chromosome1::Int64, chromosome2::Int64,
    nFactor::Array{Int64, 1}, nRow::Int64, worst_chromosome_indices::Array{Int64, 1}, totalFactor::Int64,
    nReplace::Int64)

    # Crossover points
    crossover_point_row::Int64 = rand(1:nRow)
    crossover_point_column::Int64 = rand(1:totalFactor)

    # Exchange
    for iCol = 1:crossover_point_column
        for iRow = 1:crossover_point_row
            for i = 1:nReplace
                population[iRow, iCol, worst_chromosome_indices[i]] = copy(population[iRow, iCol, chromosome2])
            end
        end
    end

    # Now replace the rest of the weakest two with the new children
    for iCol = (crossover_point_column+1):totalFactor
        for iRow = (crossover_point_row+1):nRow
            for i = 1:nReplace
                population[iRow, iCol, worst_chromosome_indices[i]] = copy(population[iRow, iCol, chromosome1])
            end
        end
    end

end



"""
mutate!(TODO)

 mutation.

 # Function to perform mutation
 # TODO THIS IS CURRENTLY INEFFICIENT, ALSO IT NEEDS TO BE MODIFIED FOR MORE MUTATIONS
"""

function mutate!(population::Array{Int64, 3}, offspring_index::Int64, mutate_prob::Float64,
    nFactor::Array{Int64, 1}, nRow::Int64, totalFactor::Int64)

    # TODO check loop ordering
    for iCol = 1:nFactor[1]
        for iRow = 1:nRow
            if rand() < mutate_prob
                population[iRow, iCol, offspring_index] = -population[iRow, iCol, offspring_index]
            end
        end
    end

    # Mutation for three level factors
    for iCol = (nFactor[1]+1):totalFactor
        for iRow = 1:nRow
            if rand() < mutate_prob
                population[iRow, iCol, offspring_index] = sample(0:2)
            end
        end
    end

end



"""
calc_fitness(TODO)


calculate fitness. This needs a lot of re-writing. Switch to just the one
obj function? Also use try-catch instead of determinant to decide on the inverse of
    the fisher information.
"""
# TODO avoid doing the same work over and over again
# TODO move X outside and update in place - dont want to do all these allocations
function calc_fitness(chromosome::Array{Int64, 2}, nModel::Int64, k::Int64,
    nFactor::Array{Int64, 1}, nRow::Int64, interactions::Array{Int64, 2}, intercept::Int64,
    totalFactor::Int64, factor_start_cols::Array{Int64, 1}, slack::Float64,
    penalty::Float64, fit_crit::String)

    X::Array{Int64, 2} = zeros( nRow, intercept + nFactor[1] + 2*nFactor[2] + k )

    # Storage for results
    average_determinant::Float64 = 0.0
    all_fisher_var::Array{Float64, 1} = zeros(nModel)
    maxvar::Float64 = 0
    minvar::Float64 = Inf
    detfisher::Float64 = 0.0

    # Number of possible models
    #nPossibleModel::Int64 = binomial(nInteraction, k)
    nInteraction::Int64 = size(interactions)[1]
    # Track which interactions are under consideration
    int_settings::Array{Int64, 1} = zeros(k)
    for i = 1:k
        int_settings[i] = i
    end

    fisher::Array{Float64, 2} = zeros( intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k )
    inv_fisher::Array{Float64, 2} = zeros( intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k )

    failed_design::Bool = false

    # Block indices for lower right part of inverse fisher information
    startIndex::Int64 = intercept + nFactor[1] + 2*nFactor[2] + 1
    endIndex::Int64 = intercept + nFactor[1] + 2*nFactor[2] + k
    detvalue::Float64 = 0.0

    # Loop over possible models, calculate fitness for each model
    # TODO parallelize this
    for iModel = 1:nModel

        # generate the corresponding design matrix for X'X
        #TODO make this only do the full thing for iModel = 1, for the rest
        # only the interactions need to be updated!!
        generate_X_matrix_mixed!(X, chromosome, iModel, interactions,
            nRow, nFactor, intercept, k, factor_start_cols, int_settings)

        # Calculate determinant
        fisher = X' * X
        detfisher = det(fisher)
        average_determinant += 1.0/nModel * detfisher

        # inv_fisher[:, :] = try
        #     inv(fisher)
        # catch
        #     failed_design = true
        #     fisher .* -999.9
        # end
        #
        # if !failed_design
        #     if k == 1
        #         all_fisher_var[iModel] = inv_fisher[intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k]
        #         if inv_fisher[intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k] > maxvar
        #             maxvar = inv_fisher[intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k]
        #         end
        #         if inv_fisher[intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k] < minvar
        #             minvar = inv_fisher[intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k]
        #         end
        #     # using determinant of lower right sub-matrix for k > 1
        #     else
        #         detvalue = det(inv_fisher[startIndex:endIndex, startIndex:endIndex])
        #         all_fisher_var[iModel] = detvalue
        #         if detvalue > maxvar
        #             maxvar = detvalue
        #         end
        #         if detvalue < minvar
        #             minvar = detvalue
        #         end
        #     end
        # end


        if round(detfisher, digits = 10) == 0
            iModel = nModel
            failed_design = true
        else
            # calculate variance of interaction terms
            inv_fisher = inv(fisher)
            if k == 1
                all_fisher_var[iModel] = inv_fisher[intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k]
                if inv_fisher[intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k] > maxvar
                    maxvar = inv_fisher[intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k]
                end
                if inv_fisher[intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k] < minvar
                    minvar = inv_fisher[intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k]
                end
            # using determinant of lower right sub-matrix for k > 1
            else
                detvalue = det(inv_fisher[startIndex:endIndex, startIndex:endIndex])
                all_fisher_var[iModel] = detvalue
                if detvalue > maxvar
                    maxvar = detvalue
                end
                if detvalue < minvar
                    minvar = detvalue
                end
            end
        end

        # Update the interaction indicator
        #print("nint is: ", nInteraction)
        #display("Updating the interaction indicator. Started at: ")
        #display(int_settings)
        for i = 0:(k-1)
            if int_settings[k-i] < (nInteraction - i)
                int_settings[k-i] += 1
                if (i > 0)
                    for j = (k-i+1):k
                        int_settings[j] = int_settings[j-1] + 1
                    end
                end
                break
            end
        end
        #display(int_settings)

    end

    if failed_design == true
        #display("failed_design")
        #print("failed design\n")
        return(-Inf)
    else
        #display("Success!")
        #print("\nFitness will be: ", (1.0/minvar) / (1 + penalty*(maxvar - minvar)*(abs(maxvar - minvar)>slack) ))

        #print("Maxvar was: ", maxvar, " and minvar was: ", minvar, "with AD: ", average_determinant, "\n")
        # calculate the fitness function
        if fit_crit == "ACV"
            return( (1.0/minvar) / (1 + penalty*(maxvar - minvar)*(abs(maxvar - minvar)>slack) ) )
        elseif fit_crit == "SqDevMean"
            #display("using new technique")
            meanval = mean(all_fisher_var)
            ssdev = sum( (all_fisher_var .- meanval).^2 )
            return(  (1.0 / meanval) / (1 + penalty * ssdev)  )
        else
            display("ERROR")
        end
    end

end

"""
calculate_worst_chromosomes!(TODO)


Find the worst, the worst get replaced.
"""
function calculate_worst_chromosomes!(fitness::Array{Float64, 1}, nReplace::Int64,
    worst_chromosome_indices::Array{Int64, 1})

    for i = 1:nReplace
        worst_chromosome_indices[i] = argmin(fitness)
        fitness[worst_chromosome_indices[i]] = Inf
    end

end

# Function to remove duplicate design points from the chromosome
# ideally use hash table with all possible design points
#TODO redo so that not re-allocating at each iteration!!!
"""
remove_duplicates!(TODO)

Under the default settings all of the design points are required to be unique. This function
finds duplicated design points and replaces them with a random point drawn from the unused
design points.

"""
function remove_duplicates!(chromosome::Array{Int64, 2},
    design_point_bank::Array{Int64, 2})

    # Step 1: Reduce to just duplicated rows, TODO how to efficiently store
    uniquerows::Array{Int64, 2} = unique(chromosome; dims=1)
    nunique::Int64 = size(uniquerows, 1)

    availableRows::Array{Int64, 1} = zeros( size(design_point_bank,1) - nunique )

    if nunique < size(chromosome, 1)

        # find the available rows
        ind::Int64 = 0
        foundindex = 0
        for icand = 1:size(design_point_bank,1)
            #print(icand, " \n")
            foundindex = findfirst([ design_point_bank[icand,:]  == uniquerows[i,:] for i=1:nunique])
            #display(design_point_bank[icand,:])
            #display(foundindex)
            if foundindex == nothing
                ind = ind + 1
                availableRows[ind] = icand
            end
        end

        # find the replacement rows
        #display(availableRows)
        newrows::Array{Int64, 1} = sample(availableRows, size(chromosome, 1) - nunique; replace=false)
        #display(newrows)

        # TODO - reordering like this re-allocates. Want to plug in new rows at
        # replacement indices instead
        # chromosome update
        chromosome[1:nunique, :] = uniquerows
        chromosome[(nunique+1):end, :] = copy(design_point_bank[newrows, :])
    end

    # TODO this needs to work in place!!!!!!!!
    return(chromosome)

end


"""
remove_duplicates_index_based!(TODO)

Under the default settings all of the design points are required to be unique. This function
finds duplicated design points and replaces them with a random point drawn from the unused
design points. This index based version should be more efficient

TODO setting for if should check all points or just a specific one/set.

    bank_index - just the column of the 2d matrix corresponding to the current chromosome
"""
function remove_duplicates_index_based!(chromosome::SubArray{Int64, 2},
    design_point_bank::Array{Int64, 2}, bank_index::SubArray{Int64, 1})

    # modified from
    # https://stackoverflow.com/questions/50899973/indices-of-unique-elements-of-vector-in-julia
    uniqueset = Set{Int64}()
    ex = eachindex(bank_index)
    idxs = Vector{eltype(ex)}()
    for i in ex
        xi = bank_index[i]
        if !(xi in uniqueset)
            push!(uniqueset, xi)
        else
            push!(idxs, i) #indices to be replaced
        end
    end
    availablerows = setdiff(1:size(design_point_bank, 1), uniqueset)

    # Draw new rows and replace
    new_design_points = sample(availablerows, length(idxs), replace=false)
    ind::Int64 = 0
    for iRow in idxs
        ind += 1
        bank_index[iRow] = new_design_points[ind]
        chromosome[iRow, :] = (design_point_bank[new_design_points[ind], :])
    end

    nothing

end



"""
generate_breeding_partners(fitness)
"""

function generate_breeding_partners(fitness::Array{Float64, 1})

    rand1::Float64 = rand()
    rand2::Float64 = rand()
    rand1final = min(rand1, rand2)
    rand2final = max(rand1, rand2)

    cumsum::Float64 = 0.0
    totalfit = sum(fitness)

    #display(fitness ./ totalfit)

    chromosome1::Int64 = 0
    chromosome2::Int64 = 0

    popsize::Int64 = length(fitness)
    for i = 1:popsize
        cumsum += fitness[i] / totalfit
        if rand1final < cumsum
            chromosome1 = i
            break
        end
    end

    if chromosome1 == 0
        chromosome1 = sample(1:popsize, 1)[1]
    end

    for i = (chromosome1):popsize
        cumsum += fitness[i] / totalfit
        if rand2final < cumsum
            #display("worked")
            chromosome2 = i
            break
        end
    end

    if chromosome2 == 0
        chromosome2 = sample(1:popsize, 1)[1]
        #display("Had to fix")
    end

    return(chromosome1, chromosome2)

end




"""
acvGA(TODO)

main function
"""

function optGA(intercept::Int64, nFactor::Array{Int64, 1}, nRow::Int64, k::Int64;
    nReplace::Int64 = 2, maxgen::Int64 = 100, popsize::Int64 = 10, slack::Float64 = 0.0,
    penalty::Float64 = -999.9, mutation_prob::Float64 = -999.9, allow_replicates::Bool = false,
    fit_crit="ACV", DUPSTYPE::Int64 = 1)

    # Calculate the number of possible interactions (2 factor)
    totalFactor::Int64 = sum(nFactor)
    nInteraction::Int64 = binomial(totalFactor, 2) + nFactor[1]*nFactor[2] + 2*binomial(nFactor[2], 2) + binomial(nFactor[2], 2)
    #print("number of interactions: ", nInteraction, "\n")
    nModel::Int64 = binomial(nInteraction, k)
    factor_start_cols::Array{Int64, 1} = generate_factor_start_cols(totalFactor, nFactor,
        intercept, k)

    # Handle case where user did not specify penalty term or gave something invalid
    if penalty <= 0.0
        penalty = 10.0^( round( 2^nFactor[1] + 3^nFactor[2] ) )
    end

    # Handle case where user did not specify mutation probability
    if mutation_prob < 0.0
        mutation_prob = 1.0 / (nRow * sum(nFactor))
    end

    # create the bank of possible design points
    design_point_bank::Array{Int64, 2}, flip_indices = create_design_point_bank(nFactor)

    # Convergence criteria
    converged::Bool = false
    generation::Int64 = 0

    # Storage for all chromosomes
    population::Array{Int64, 3} = zeros(nRow, totalFactor, popsize)
    # Keeps track of which row of the bank each design point is for each chromosome
    bank_index::Array{Int64, 2} = zeros(nRow, popsize)
    # Storage for all fitness values
    fitness::Array{Float64, 1} = zeros(popsize)
    # chromosomerate the matrix of possible interactions
    interactions::Array{Int64, 2} = generate_interaction_matrix_mixed(nInteraction, nFactor)
    # Store the array of indices for the worst chromosomes
    worst_chromosome_indices::Array{Int64, 1} = zeros(nReplace)

    # Variables for storing running time
    generation_time::Float64 = 0.0
    fitness_time::Float64 = 0.0
    crossover_time::Float64 = 0.0
    mutation_time::Float64 = 0.0
    duprem_time::Float64 = 0.0
    lookup_time::Float64 = 0.0
    removal_time::Float64 = 0.0

    # chromosomerate starting values for all the chromosomes
    #print("generating initial population...\n")
    time0 = time()
    for ichromosome = 1:popsize

        population[:,:,ichromosome] = generate_chromosome_mixed(nRow, nFactor)

        if !allow_replicates
            if DUPSTYPE == 1

                for iRow = 1:nRow
                    bank_index[iRow, ichromosome] = lookup_bank_row(design_point_bank, flip_indices, population[iRow, :, ichromosome])
                end

                remove_duplicates_index_based!(
                    view(population, :, :, ichromosome),
                    design_point_bank, view(bank_index, :, ichromosome) )

            # Phase this out
            elseif DUPSTYPE == 0
                population[:,:,ichromosome] = remove_duplicates!(population[:,:,ichromosome], design_point_bank)
            end
        end

        fitness[ichromosome] = calc_fitness(population[:,:,ichromosome], nModel, k,
            nFactor, nRow, interactions, intercept, totalFactor, factor_start_cols, slack, penalty, fit_crit)

    end

    generation_time = time() - time0

    while !converged
        generation += 1

        ##### Perform crossover
        # Select two chromosomes for breeding
        # cumulative probability
        chromosome1, chromosome2 = generate_breeding_partners(fitness)

        # find the current worst chromosomes
        calculate_worst_chromosomes!(fitness, nReplace, worst_chromosome_indices)

        time0 = time()
        crossover!(population, chromosome1, chromosome2, nFactor, nRow,
            worst_chromosome_indices, totalFactor, nReplace)
        crossover_time += time() - time0

        ##### Calc Fitness
        for ichromosome = 1:nReplace

            time0 = time()
            if rand() < 0.95
                time0 = time()
                mutate!(population, worst_chromosome_indices[ichromosome], mutation_prob, nFactor, nRow, totalFactor)
                mutation_time += time() - time0
                time0 = time()

                if !allow_replicates
                    if DUPSTYPE == 1
                        tt0 = time()
                        for iRow = 1:nRow
                            bank_index[iRow, worst_chromosome_indices[ichromosome]] = lookup_bank_row(design_point_bank, flip_indices, population[iRow, :, worst_chromosome_indices[ichromosome]])
                        end
                        lookup_time += time() - tt0
                        tt0 = time()
                        remove_duplicates_index_based!(
                            view(population, :, :, worst_chromosome_indices[ichromosome]),
                            design_point_bank, view(bank_index, :, worst_chromosome_indices[ichromosome]) )
                        removal_time += time() - tt0
                    # Phase this out
                    elseif DUPSTYPE == 0
                        #print("here")
                        population[:,:,worst_chromosome_indices[ichromosome]] = remove_duplicates!(population[:,:,worst_chromosome_indices[ichromosome]], design_point_bank)
                    end
                    #population[:,:,worst_chromosome_indices[ichromosome]] = remove_duplicates!(population[:,:,worst_chromosome_indices[ichromosome]], design_point_bank)
                    #remove_duplicates!(population[:,:,worst_chromosome_indices[ichromosome]], design_point_bank)
                end

                duprem_time += time() - time0
            else
                time0 = time()
                population[:,:,worst_chromosome_indices[ichromosome]] = generate_chromosome_mixed(nRow, nFactor)
                generation_time += time() - time0
                time0 = time()

                if !allow_replicates
                    if DUPSTYPE == 1
                        tt0 = time()
                        for iRow = 1:nRow
                            bank_index[iRow, worst_chromosome_indices[ichromosome]] = lookup_bank_row(design_point_bank, flip_indices, population[iRow, :, worst_chromosome_indices[ichromosome]])
                        end
                        lookup_time += time() - tt0
                        tt0 = time()
                        remove_duplicates_index_based!(
                            view(population, :, :, worst_chromosome_indices[ichromosome]),
                            design_point_bank, view(bank_index, :, worst_chromosome_indices[ichromosome]) )
                        removal_time += time() - tt0

                    # Phase this out
                    elseif DUPSTYPE == 0
                        #print("here")
                        population[:,:,worst_chromosome_indices[ichromosome]] = remove_duplicates!(population[:,:,worst_chromosome_indices[ichromosome]], design_point_bank)
                    end
                    #population[:,:,worst_chromosome_indices[ichromosome]] = remove_duplicates!(population[:,:,worst_chromosome_indices[ichromosome]], design_point_bank)
                    #remove_duplicates!(population[:,:,worst_chromosome_indices[ichromosome]], design_point_bank)
                end

                duprem_time += time() - time0
            end

            time0 = time()
            fitness[worst_chromosome_indices[ichromosome]] = calc_fitness(population[:,:,worst_chromosome_indices[ichromosome]],
                nModel, k, nFactor, nRow, interactions, intercept, totalFactor, factor_start_cols, slack, penalty, fit_crit)
            fitness_time += time() - time0
        end

        if generation > maxgen
            #print("Maximum generations reached")
            converged = true
        end

        #TODO convergence currently turned off
        if ((minimum(fitness) / maximum(fitness)) > 1.99999) & (minimum(fitness) > 0.0)
            print("iter: ", generation, "\n")
            print("max fitness: ", maximum(fitness), "\n")
            print("min fitness: ", minimum(fitness))
            #print("Population has converged on generation ", generation)
            converged = true
        end

    end

    # # # # Display timing information
    # display(string("============"))
    # display(string("Generation time:   ", generation_time))
    # display(string("fitness time   :   ", fitness_time))
    # display(string("crossover time :   ", crossover_time))
    # display(string("mutation time  :   ", mutation_time))
    # display(string("duprem time    :   ", duprem_time))
    # display(string("lookup time    :   ", lookup_time))
    # display(string("removal time   :   ", removal_time))

    bestdesign::Int64 = argmax(fitness)
    #display(fitness)
    sortedDesign::Array{Int64, 2} = sortslices( population[:,:,bestdesign], dims=1 )

    # Calculate the final fitness value. Should be the same as the fitness value
    # above. this serves as a sanity check
    fitval = calc_fitness(sortedDesign, nModel, k, nFactor, nRow, interactions,
        intercept, totalFactor, factor_start_cols, slack, penalty, fit_crit)

    finaldesign::Design = Design(sortedDesign, nFactor,
        interactions, k, nModel, intercept, nRow)

    return(finaldesign, fitval)

end
