using QuantumCumulants
using CatTDVP
using Test



@testset "TDVP Type" begin
    hf = FockSpace(:bosonic_mode)
    a₁ = Destroy(hf, Symbol("a",1), 1)
    a₂ = Destroy(hf, Symbol("a",2), 2)
    a₃ = Destroy(hf, Symbol("a",2), 3)

    #single mode
    H = a₁^2 + (a₁')^2 + a₁ + a₁' + a₁'*a₁
    J = [a₁, a₁^2]
    rates = [1.0, 2.0]
    order = [2]
    sys = @test_nowarn TDVPSystem(H, J, order; rates=rates)
    mat = @test_nowarn TDVPMatrices(sys)
    @test try functions = TDVPFunction(sys)
        true
    catch
        false
    end
    f_dae = make_ODE_problem(sys)
    dimρ = (2*(order[1]+1))^2
    u0 = randn(ComplexF64, 1+dimρ)
    du0 = randn(ComplexF64, 1+dimρ)
    res = zeros(ComplexF64, 1+dimρ)
    @test_nowarn f_dae(res, u0, du0, 0.0, 0.0)
    @test !iszero(res)
    



    
    #two coupled modes
    G₁, G₂ = [4, 4]
    η₁, η₂ = [1.0, 1.0]
    J_hop = 1.0
    H = G₁*(a₁^2 + (a₁')^2) + G₂*(a₂^2 + (a₂')^2) + J_hop*(a₁*a₂' + a₁'*a₂)
    rates = [η₁, η₂]
    J = [a₁^2, a₂^2]

    order = [1, 1]
    sys = @test_nowarn TDVPSystem(H, J, order; rates=rates)
    mat = @test_nowarn TDVPMatrices(sys)
    @test try functions = TDVPFunction(sys)
        true
    catch
        false
    end
    f_dae = make_ODE_problem(sys)
    dimρ = (2*(order[1]+1))^2
    u0 = randn(ComplexF64, 2+dimρ)
    du0 = randn(ComplexF64, 2+dimρ)
    res = zeros(ComplexF64, 2+dimρ)
    @test_nowarn f_dae(res, u0, du0, 0.0, 0.0)
    @test !iszero(res)
end