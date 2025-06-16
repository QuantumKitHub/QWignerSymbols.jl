using QWignerSymbols, WignerSymbols
smalljlist = 0:(1 // 2):5

# j=1
# for m1 in (-j):j, m2 in (-j):j, m3 in (-j):j
#     @show m1, m2, m3
#     if !isapprox(q_wigner3j(j, j, j, m1, m2, m3, 1.0), wigner3j(j, j, j, m1, m2, m3))
#         @show q_wigner3j(j, j, j, m1, m2, m3, 1.0), wigner3j(j, j, j, m1, m2, m3)
#     end
# end

@testset "q = 1" begin
    q = 1
    for j1 in smalljlist, j2 in smalljlist
        for j3 in abs(j1 - j2):min(5, (j1 + j2))
            for m1 in (-j1):j1, m2 in ((-j2):j2), m3 in (-j3):(j3)
                if !isapprox(q_wigner3j(j1, j2, j3, m1, m2, m3, q),
                             wigner3j(j1, j2, j3, m1, m2, m3); atol=1e-14)
                    @test q_wigner3j(j1, j2, j3, m1, m2, m3, q) ≈
                          wigner3j(j1, j2, j3, m1, m2, m3) atol = 1e-14
                    @show j1, j2, j3, m1, m2, m3
                end
            end
        end
    end

    for j1 in smalljlist, j2 in smalljlist, j3 in smalljlist, j4 in smalljlist,
        j5 in smalljlist, j6 in smalljlist

        @test q_wigner6j(j1, j2, j3, j4, j5, j6, q) ≈ wigner6j(j1, j2, j3, j4, j5, j6) atol = 1e-14
    end
end

@testset "clebschgordan: orthogonality relations" begin
    q = 1.1
    for j1 in smalljlist, j2 in smalljlist
        d1::Int = 2 * j1 + 1
        d2::Int = 2 * j2 + 1
        M = zeros(d1 * d2, d1 * d2)
        ind1 = 1
        for m1 in (-j1):j1, m2 in (-j2):j2
            ind2 = 1
            @inbounds for j3 in abs(j1 - j2):(j1 + j2), m3 in (-j3):j3
                M[ind1, ind2] = q_clebschgordan(j1, m1, j2, m2, j3, m3, q)
                ind2 += 1
            end
            ind1 += 1
        end
        @test M' * M ≈ one(M)
    end
end

@testset "wigner3j: orthogonality relations" begin
    q = 1.1
    for j1 in smalljlist, j2 in smalljlist
        d1::Int = 2 * j1 + 1
        d2::Int = 2 * j2 + 1
        M = zeros(d1 * d2, d1 * d2)
        ind2 = 1
        for m1 in (-j1):j1, m2 in (-j2):j2
            ind1 = 1
            @inbounds for j3 in abs(j1 - j2):(j1 + j2), m3 in (-j3):j3
                d3::Int = 2 * j3 + 1
                M[ind1, ind2] += sqrt(q_number(d3, q)) *
                                 q_wigner3j(j1, j2, j3, m1, m2, m3, q)
                ind1 += 1
            end
            ind2 += 1
        end
        @test M' * M ≈ one(M) # orthogonality relation type 1
        @test M * M' ≈ one(M) # orthogonality relation type 2
    end
end

@testset "wigner6j - wigner3j: consistency relations" begin end

# # @testset "wigner6j: orthogonality relations" begin
# #     for j1 in smalljlist, j2 in smalljlist, j4 in smalljlist
# #         for j5 in
# #             max(abs(j1 - j2 - j4), abs(j1 - j2 + j4), abs(j1 + j2 - j4)):(j1 + j2 + j4)
# #             j6range = max(abs(j2 - j4), abs(j1 - j5)):min((j2 + j4), (j1 + j5))
# #             j3range = max(abs(j1 - j2), abs(j4 - j5)):min((j1 + j2), (j4 + j5))
# #             @test length(j6range) == length(j3range)
# #             M = zeros(Float64, (length(j3range), length(j6range)))
# #             for (k2, j6) in enumerate(j6range)
# #                 for (k1, j3) in enumerate(j3range)
# #                     M[k1, k2] = sqrt(q_number(2 * j3 + 1, q) * q_number(2 * j6 + 1, q)) *
# #                                 q_wigner6j(j1, j2, j3, j4, j5, j6, q)
# #                 end
# #             end
# #             @test M' * M ≈ one(M)
# #         end
# #     end
# # end
