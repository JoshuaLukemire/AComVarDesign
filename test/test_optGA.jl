
"""
test_X_gen(TODO)

tests X fill out. needs to be moved to tests.
"""
function test_X_gen()
nRow = 6
nFactor = [0, 2]
nCol = nFactor[1] + nFactor[2]
totalFactor = sum(nFactor)
factor_start_cols = zeros(totalFactor)
chromosome::Array{Int64, 2} = generate_chromosome_mixed(nRow, nFactor)
intercept = 1
k=1
factor_start_cols::Array{Int64, 1} = zeros(intercept + totalFactor + k)
factor_start_cols[1] = 1
for i = (intercept+1):(nFactor[1]+intercept)
    factor_start_cols[i] = i
end
add = -1
for i = (nFactor[1]+intercept+1):(totalFactor+intercept)
    add = add+1
    factor_start_cols[i] = i + add
end
for i = 1:(totalFactor+intercept+k-(totalFactor+intercept))
    factor_start_cols[totalFactor+intercept+i] = nFactor[1]+intercept+2*nFactor[2]+i
end

# add interaction starting points
display(factor_start_cols)
nInteraction = binomial(totalFactor, 2) + nFactor[1]*nFactor[2] + binomial(nFactor[2], 2) + binomial(nFactor[2], 2)
interactions = generate_interaction_matrix_mixed(nInteraction, nFactor)
X::Array{Int64, 2} = zeros(nRow, intercept + nFactor[1] + 2*nFactor[2] + k)
iModel = 1
generate_X_matrix_mixed!(X, chromosome, iModel, interactions,
    nRow, nFactor, intercept, k, factor_start_cols)
display(X)
end
