"""
    Design

Structure that stores the output from the genetic algorithm. Variance ratio can
be examined using min_var / max_var or by looking at variance_array.

* `design_matrix` is a 2 dimensional array of integers with the main effects
* `nFactor` is a 1 dimensional array. This first element is the number of 2-level factors and the second element is the number of 3-level factors
* `nInteraction` number of interactions in the models under consideration (k in paper)
* `variance_array` Array storing the variance of each model under consideration.
Can be matched up to `interaction_matrix`
"""

struct Design

    # Design Matrix
    design_matrix::Array{Int64, 2}

    # Number of two and three level factors
    nFactor::Array{Int64, 1}

    # Number of interactions in the model (k)
    nInteraction::Int64

    # Array of the interaction variances/dets for each model
    variance_array::Array{Float64, 1}

    # Min and Max variance across models under consideration
    max_var::Float64
    min_var::Float64

    # Array of the determinant of the information matrix for each model
    det_array::Array{Float64, 1}

    # Matrix of possible interactions (models under consideration)
    interaction_matrix::Array{Int64, 2}

end


"""
Constructor
"""
function Design(chromosome::Array{Int64, 2}, nFactor::Array{Int64, 1},
    interactions::Array{Int64, 2}, nInteraction::Int64, nModel::Int64, intercept::Int64,
    nRow::Int64)

    # Store the design matrix
    design_matrix::Array{Int64, 2} = chromosome
    variance_array::Array{Float64, 1} = zeros(nModel)
    max_var::Float64 = 0.0
    min_var::Float64 = Inf

    # Array of the determinant of the information matrix for each model
    det_array::Array{Float64, 1} = zeros(nModel)

    nPossibleInteraction::Int64 = size(interactions)[1]

    k::Int64 = nInteraction

    factor_start_cols::Array{Int64, 1} = generate_factor_start_cols(sum(nFactor), nFactor,
        intercept, k)

    # Calculate some summary information
    X::Array{Int64, 2} = zeros( nRow, intercept + nFactor[1] + 2*nFactor[2] + k )

    # Storage for results
    detfisher::Float64 = 0.0

    int_settings::Array{Int64, 1} = zeros(k)
    for i = 1:k
        int_settings[i] = i
    end
    # Block indices for lower right part of inverse fisher information
    startIndex::Int64 = intercept + nFactor[1] + 2*nFactor[2] + 1
    endIndex::Int64 = intercept + nFactor[1] + 2*nFactor[2] + k
    detvalue::Float64 = 0.0

    fisher::Array{Float64, 2} = zeros( intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k )
    inv_fisher::Array{Float64, 2} = zeros( intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k)

    # Loop over possible models, calculate fitness for each model
    for iModel = 1:nModel

        # generate the corresponding design matrix for X'X
        generate_X_matrix_mixed!(X, chromosome, iModel, interactions, nRow,
            nFactor, intercept, k, factor_start_cols, int_settings)

        fisher = X' * X
        detfisher = det(fisher)
        det_array[iModel] = detfisher

        # try
        #     inv_fisher(fisher)
        #     if k == 1
        #         if inv_fisher[intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k] > max_var
        #             max_var = inv_fisher[intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k]
        #         end
        #         if inv_fisher[intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k] < min_var
        #             min_var = inv_fisher[intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k]
        #         end
        #         variance_array[iModel] = inv_fisher[intercept + nFactor[1] + 2*nFactor[2] + k,
        #             intercept + nFactor[1] + 2*nFactor[2] + k]
        #     else
        #         detvalue = det(inv_fisher[startIndex:endIndex, startIndex:endIndex])
        #         if detvalue > max_var
        #             max_var = detvalue
        #         end
        #         if detvalue < min_var
        #             min_var = detvalue
        #         end
        #         variance_array[iModel] = detvalue
        #     end
        # catch
        #     max_var = 9999999999.9
        #     min_var = 0.0
        #     variance_array[iModel] = 9999999999.9
        # end

        # calculate variance of interaction terms
        if !(round(detfisher, digits = 20) == 0)
            inv_fisher = inv(fisher)
            if k == 1
                if inv_fisher[intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k] > max_var
                    max_var = inv_fisher[intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k]
                end
                if inv_fisher[intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k] < min_var
                    min_var = inv_fisher[intercept + nFactor[1] + 2*nFactor[2] + k, intercept + nFactor[1] + 2*nFactor[2] + k]
                end
                variance_array[iModel] = inv_fisher[intercept + nFactor[1] + 2*nFactor[2] + k,
                    intercept + nFactor[1] + 2*nFactor[2] + k]
            else
                detvalue = det(inv_fisher[startIndex:endIndex, startIndex:endIndex])
                if detvalue > max_var
                    max_var = detvalue
                end
                if detvalue < min_var
                    min_var = detvalue
                end
                variance_array[iModel] = detvalue
            end

        else
            max_var = 9999999999.9
            min_var = 0.0
            variance_array[iModel] = 9999999999.9
        end

        # Update the model under consideration
        for i = 0:(k-1)
            if int_settings[k-i] < (nPossibleInteraction - i)
                int_settings[k-i] += 1
                if (i > 0)
                    for j = (k-i+1):k
                        int_settings[j] = int_settings[j-1] + 1
                    end
                end
                break
            end
        end

    end

    return( Design(design_matrix, nFactor, nInteraction, variance_array,
        max_var, min_var, det_array, interactions) )

end


# Extend printing function for a design
function Base.show(io::IO, SelDesign::Design)

    # Tell the user what the design was originally constructed for
    print("\nACV design created with ", SelDesign.nFactor, " factors and ", size(SelDesign.design_matrix, 1), " rows.\n")
    print("The design is:\n")
    display(SelDesign.design_matrix)
    display(hcat(SelDesign.variance_array, SelDesign.det_array))
    print("\n")

end
