using Test
using Random

@testset "Austin (et.al.) examples" begin
    Random.seed!(137)

    test_cases = [
#       CAS     K     M     N    EPS
        (1,   2000,   20,   4,  1e-12), #pass
        (2,   2000,   20,  10,  1e-8),
        (3,   2000,   20,  10,  1e-11),
        (4,   2000,   20,  10,  1e-12), #pass
        (5,   2000,   10,  10,  2e-04),
        (6,   2000,    3,   3,  1e-12),
        (7,   2000,    2,   2,  1e-12),
        (8,   2000,    4,   4,  1e-12),
        (9,   2000,    2,   2,  1e-12),
        (10,  2000,    2,   3,  1e-12),
        (11,  2000,    3,   2,  1e-12),
        (12,  2000,    4,   4,  1e-12),
        (13,  2000,    3,   4,  1e-12),
        (14,  2000,    4,   3,  1e-12),
        #(15,  6000,   10, 10,  2e-04),  # test case needs work on bounds
        #(16,  6000,    15,  4,  2e-04), # will pass but takes too long
        (17,  4000,    10,   2,  1e-03),
        (18, 10000,    16,   6,  1e-06),
        #(19,  2000,   10,  10,  1e-12), # takes too long
        (20,  2000,    2,  10,  1e-12),
    ]

    for cas in test_cases
        @testset "CAS=$cas[1]" begin
            err = run_case(cas[1]; K=cas[2], M=cas[3], N=cas[4])
            @test err < cas[5]
            @show cas[1], err
        end
    end
end
