"""
Functions for working with different models
"""


"""
generate_interaction_matrix_mixed(nInteraction, nFactor)

Returns a matrix that keeps track of what interactions are active. The number of rows in
this matrix corresponds to the number of possible models under consideration. This version
of the function handles both 2 and 3 level factors. This first columns correspond to
the two level factors and the second set of columns corresponds to the 3-level factors.
"""

function generate_interaction_matrix_mixed(nInteraction::Int64, nFactor::Array{Int64, 1})

    # Figure out how many interactions (rows) and factors (cols) are in the model
    totalFactors::Int64 = sum(nFactor)
    interactions::Array{Int64, 2} = zeros(nInteraction, totalFactors)

    # Indices for tracking where we are in the interaction matrix
    iEnd::Int64 = 0
    iStart::Int64 = 0
    nextpos::Int64 = 0

    # Handle the standard interactions (1s not 2s)
    for iFactor = 1:(totalFactors-1)
        iStart = iEnd + 1
        iEnd = iStart + (totalFactors - iFactor) - 1
        #print("Start is ", iStart, " end is: ", iEnd, "\n")
        interactions[ iStart:iEnd, iFactor] .= 1
        # Now add identity matrix to the rest of the columns
        nextpos = iFactor
        for idiag = iStart:iEnd
            nextpos += 1
            interactions[idiag, nextpos] = 1
        end
    end

    # Interactions between 2-level factors, first level of three level factors, and the "2" setting for 3 level factors
    iStart = binomial(totalFactors, 2)
    iEnd = iStart
    # loop is only over the 2 level factors
    for iFactor = 1:(totalFactors - 1)
        iStart = iEnd+1
        if iFactor < nFactor[1]
            iEnd = iEnd + nFactor[2]
        else
            iEnd = iEnd + totalFactors - iFactor
        end
        interactions[iStart:iEnd, iFactor] .= 1
        if iFactor < nFactor[1]
            nextpos = nFactor[1]
        else
            nextpos = iFactor
        end
        for idiag = iStart:iEnd
            nextpos += 1
            interactions[idiag, nextpos] = 2
        end
    end

    # Quadratic - Linear Interactions
    iStart = iEnd
    for iFactor = (nFactor[1] + 1):totalFactors
        iStart = iEnd+1
        iEnd = iStart + (nFactor[2] - iFactor + nFactor[1]) - 1
        interactions[iStart:iEnd, iFactor] .= 2
        nextpos = iFactor
        for idiag = iStart:iEnd
            nextpos += 1
            interactions[idiag, nextpos] = 1
        end
    end

    # Interactions between 2-level factors for 3 level factors
    iStart = iEnd
    for iFactor = (nFactor[1] + 1):totalFactors
        iStart = iEnd+1
        iEnd = iStart + (nFactor[2] - iFactor + nFactor[1]) - 1
        interactions[iStart:iEnd, iFactor] .= 2
        nextpos = iFactor
        for idiag = iStart:iEnd
            nextpos += 1
            interactions[idiag, nextpos] = 2
        end
    end

    return(interactions)
end



"""
generate_factor_start_cols(TODO)

Function to generate an array (factor_start_cols) that keeps track of where
each factor is in the X matrix. This is needed for three level factors since
they use up two columns each.
"""

function generate_factor_start_cols(totalFactor::Int64, nFactor::Array{Int64, 1},
    intercept::Int64, k::Int64)

    factor_start_cols::Array{Int64, 1} = zeros(intercept + totalFactor + k)
    factor_start_cols[1] = 1
    for i = (intercept+1):(nFactor[1]+intercept)
        factor_start_cols[i] = i
    end
    add::Int64 = -1
    for i = (nFactor[1]+intercept+1):(totalFactor+intercept)
        add = add+1
        factor_start_cols[i] = i + add
    end
    for i = 1:(totalFactor+intercept+k-(totalFactor+intercept))
        factor_start_cols[totalFactor+intercept+i] = nFactor[1]+intercept+2*nFactor[2]+i
    end

    return(factor_start_cols)
end


"""
create_design_point_bank(TODO)

Creates a list of all possible design points.
"""
function create_design_point_bank(nFactor::Array{Int64, 1})

    nPoint::Int64 = 2^nFactor[1] * 3^nFactor[2]
    nCol::Int64 = sum(nFactor)
    design_point_bank::Array{Int64, 2} = zeros(nPoint, nCol)

    # Array to store flipping indices
    flip_indices::Array{Int64, 1} = zeros(nCol)

    # Loop over each design point
    flipind::Int64 = nPoint
    cval::Int64 = -999
    for iCol = 1:nCol
        # 2 level factors
        if iCol <= nFactor[1]
            flipind = flipind / 2
            cval = -1
            for iRow = 1:nPoint
                design_point_bank[iRow, iCol] = cval
                if mod(iRow, flipind) == 0
                    cval = -cval
                end
            end
        # 3 level factors
        else
            flipind = flipind / 3
            cval = 0
            for iRow = 1:nPoint
                design_point_bank[iRow, iCol] = cval
                if mod(iRow, flipind) == 0
                    cval = mod(cval + 1, 3)
                end
            end
        end
        # Store the flip index
        flip_indices[iCol] = flipind
    end
    return(design_point_bank, flip_indices)
end


"""
    lookup_bank_row()

    figure out which row of the design point bank the reference point comes from
"""
function lookup_bank_row(design_point_bank::Array{Int64, 2}, flip_indices::Array{Int64, 1},
    design_point::Array{Int64, 1})
    nRow, nCol = size(design_point_bank)
    row_index = 1
    identified = false
    for iCol = 1:nCol
        identified = false
        while !identified
            if design_point_bank[row_index, iCol] != design_point[iCol]
                row_index += flip_indices[iCol]
            else
                identified = true
            end
        end
    end
    return(row_index)
end
