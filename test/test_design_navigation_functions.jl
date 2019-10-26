@testset "Interaction matrix generation" begin

    # Test that two factor only interactions come through correctly
    function testint1()
        nFactor = [4, 0]
        intmatrix = generate_interaction_matrix_mixed(binomial(4, 2), nFactor)
        return(all(intmatrix .== [1 1 0 0; 1 0 1 0; 1 0 0 1; 0 1 1 0; 0 1 0 1; 0 0 1 1]))
    end
    @test testint1()

    # Test that three factor only interactions come through correctly
    function testint2()
        nFactor = [0, 3]
        intmatrix = generate_interaction_matrix_mixed(12 , nFactor)
        return(all(intmatrix .== [1 1 0; 1 0 1; 0 1 1; 1 2 0; 1 0 2; 0 1 2;
            2 1 0; 2 0 1; 0 2 1; 2 2 0; 2 0 2; 0 2 2]))
    end
    @test testint2()

end

@testset "X matrix generation from chromosome" begin

    # pg 6, Ghosh and Chowdhury (2017)
    function testcase1()
        nRow = 8
        nFactor = [0, 3]
        #nCol = nFactor[1] + nFactor[2]
        totalFactor = sum(nFactor)
        intercept = 1
        k = 1
        chromosome = [2 2 2;
                      1 1 2;
                      1 2 1;
                      0 0 1;
                      0 2 0;
                      0 0 0;
                      2 1 1;
                      0 1 2]
        factor_start_cols::Array{Int64, 1} = generate_factor_start_cols(totalFactor, nFactor,
            intercept, k)
        nInteraction = binomial(totalFactor, 2) + nFactor[1]*nFactor[2] + 2*binomial(nFactor[2], 2) + binomial(nFactor[2], 2)
        interactions = generate_interaction_matrix_mixed(nInteraction, nFactor)
        #print(interactions)
        X::Array{Int64, 2} = zeros(nRow, intercept + nFactor[1] + 2*nFactor[2] + k)
        # Store the results
        var_results::Array{Float64, 1} = zeros(nInteraction)
        for iModel = 1:nInteraction
            #print(iModel, " ")
            #display(interactions[iModel, :])
            generate_X_matrix_mixed!(X, chromosome, iModel, interactions,
                nRow, nFactor, intercept, k, factor_start_cols)
            #display(X)
            fisher = X' * X
            invfisher = inv(fisher)
            var_results[iModel] = invfisher[8,8]
        end


        return( all(isapprox.(var_results, 0.666, atol=0.001)) )
    end
    @test testcase1()

    function testcase2()
        nRow = 10
        nFactor = [0, 3]
        #nCol = nFactor[1] + nFactor[2]
        totalFactor = sum(nFactor)
        intercept = 1
        k = 1
        chromosome = [2 2 2;
                      1 1 2;
                      1 2 1;
                      0 0 1;
                      0 2 0;
                      1 1 1;
                      0 2 1;
                      2 0 1;
                      2 0 0;
                      0 1 1]
        factor_start_cols::Array{Int64, 1} = generate_factor_start_cols(totalFactor, nFactor,
            intercept, k)
        nInteraction = binomial(totalFactor, 2) + nFactor[1]*nFactor[2] + 2*binomial(nFactor[2], 2) + binomial(nFactor[2], 2)
        interactions = generate_interaction_matrix_mixed(nInteraction, nFactor)
        #print(interactions)
        X::Array{Int64, 2} = zeros(nRow, intercept + nFactor[1] + 2*nFactor[2] + k)
        # Store the results
        var_results::Array{Float64, 1} = zeros(nInteraction)
        for iModel = 1:nInteraction
            #print(iModel, " ")
            generate_X_matrix_mixed!(X, chromosome, iModel, interactions,
                nRow, nFactor, intercept, k, factor_start_cols)
            #display(X)
            fisher = X' * X
            invfisher = inv(fisher)
            var_results[iModel] = invfisher[8,8]
        end

        return( all(isapprox.(var_results, 0.2564, atol=0.001)) )
    end
    @test testcase2()


end;
