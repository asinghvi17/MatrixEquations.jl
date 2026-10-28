using Test, MatrixEquations, LinearAlgebra
@testset "default Riccati backend" begin
    A=[0. 1.;-1. -.4];B=reshape([0.,2.],2,1);Q=Matrix{Float64}(I,2,2);R=ones(1,1)
    X=arec(A,B,R,Q)[1]
    @test norm(A'*X+X*A-X*B*(R\B')*X+Q)<1e-12
    @test garec(A,Q,B,R,Q)[1] ≈ X atol=1e-12
    A1=reshape([-2.],1,1);A2=Matrix(Diagonal([-1.,-3.]));G=ones(2,1);C=ones(1,2)
    Z=nare(A1,A2,G,C)[1]
    @test size(Z)==(1,2)
    @test norm(A1*Z+Z*A2-Z*G*Z+C)<1e-12
end
