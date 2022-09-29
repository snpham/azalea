begin
	using Roots
	using DataFrames
	using DifferentialEquations
	using LinearAlgebra
	using Plots
	using PlutoUI

	include("../helpers/propagation.jl")

    # compute masses from GM
    m_earth = gm_earth / G_univ
    m_moon = gm_moon / G_univ
    m_sun = gm_sun / G_univ
    
    ## earth-moon
    # characteristic mass
    m_star_em = m_earth + m_moon
    μ_em = m_moon / m_star_em
    M2_em_nd = μ_em
    M1_em_nd = 1 - μ_em
end


function cr3bp_stm(dx, x, μ, t)

    # transposing for easy element access
    x = x'

    # compute positions
    r1 = sqrt((x[1]+μ)^2 + x[2]^2 + x[3]^2)
    r2 = sqrt((x[1]-1+μ)^2 + x[2]^2 + x[3]^2)

    # compute CR3BP
    dx[1] = x[4]
    dx[2] = x[5]
    dx[3] = x[6]
    dx[4] = 2*x[5] + x[1] - (1-μ)*(x[1]+μ)/r1^3 - (μ*(x[1]-1+μ))/r2^3
    dx[5] = -2*x[4] + x[2] - (1-μ)*x[2]/r1^3 - μ*x[2]/r2^3
    dx[6] = -(1-μ)*x[3]/r1^3 - μ*x[3]/r2^3

    # compute U matrix -> psuedo-potential 2nd derivatives
    Uxx = 1 - (1-μ)/r1^3 - μ/r2^3 + 3*(1-μ)*(x[1]+μ)^2/r1^5 + 3*μ*(x[1]-1+μ)^2/r2^5
    Uxy = 3*(1-μ)*(x[1]+μ)*x[2]/r1^5 + 3*μ*(x[1]-1+μ)*x[2]/r2^5
    Uxz = 3*(1-μ)*(x[1]+μ)*x[3]/r1^5 + 3*μ*(x[1]-1+μ)*x[3]/r2^5
    Uyx = Uxy
    Uyy = 1 - (1-μ)/r1^3 - μ/r2^3 + 3*(1-μ)*x[2]^2/r1^5 + 3*μ*x[2]^2/r2^5
    Uyz = 3*(1-μ)*x[2]*x[3]/r1^5 + 3*μ*x[2]*x[3]/r2^5
    Uzx = Uxz
    Uzy = Uyz
    Uzz = - (1-μ)/r1^3 - μ/r2^3 + 3*(1-μ)*x[3]^2/r1^5 + 3*μ*x[3]^2/r2^5

    # construct A matrix (jacobian of the state vector)
    A = [ 0.0  0.0  0.0  1.0  0.0  0.0;
          0.0  0.0  0.0  0.0  1.0  0.0;
          0.0  0.0  0.0  0.0  0.0  1.0;
          Uxx  Uxy  Uxz  0.0  2.0  0.0; 
          Uyx  Uyy  Uyz -2.0  0.0  0.0;
          Uzx  Uzy  Uzz  0.0  0.0  0.0]

    # get phi matrix (retransposing to normal shape)
    𝚽 = reshape(x[7:42], 6, 6)'

    # compute 𝚽_dot
    𝚽_dot = A * 𝚽
    dx[7:42] = 𝚽_dot'

    # if convert back to matrix for (easy element access)
    # dx = reshape(dx, 7,6)'

end


# A = [1  2  3  4  5  6 ;
#      7  8  9  10 11 12;
#      13 14 15 16 17 18;
#      19 20 21 22 23 24;
#      25 26 27 28 29 30;
#      31 32 33 34 35 36;
#      37 38 39 40 41 42]

if abspath(PROGRAM_FILE) == @__FILE__

## 1 Write a script to simultaneously numerically integrate the state vector and 
# the state transition matrix along a trajectory in the CR3BP. 

    # initial state (CRTBP)
    x_0 = [0.98; 0; 0; 0; 1.7; 0]

    # at time t_0, the initial condition for the STM is the I matrix
    𝚽_0 = Matrix(1.0I, 6, 6)

    # input is a matrix w/ state in first row,
    # STM in rows 2-7
    X0 = vcat(x_0', 𝚽_0)





    # initial condition 2
    x_bar = [0.98, 0, 0, 0, 1.2,0]
    
    # at time t_0, the initial condition for the STM is the I matrix
    𝚽_0 = Matrix(1.0I, 6, 6)

    # input is a matrix w/ state in first row,
    # STM in rows 2-7
    X0 = vcat(x_bar', 𝚽_0)

    # 2 nondimensional units
    tspan = (0.0,2.0)

    # setting up solver
    prob = ODEProblem(cr3bp_stm, X0, tspan, μ_em)
    sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12, maxiters=1e6)
    
    # # 3d plot in earth-moon rotating frame
    # plot(sol[1,:], sol[2,:], sol[3,:], label="trajectory", title=" i.c. #1 (e-m rotating frame)", 
    #     xlabel="x [-]", ylabel="y [-]", zlabel="z [-]", layout = 2, aspect_ratio = 1)
    
    # # adding the moon
    # scatter!([1-mu_earthmoon],[0], [0], color="gray", markersize = 4, label="m2: moon")

    # # 2d plot in earth-moon rotating frame
    # plot!(sol[1,:], sol[2,:], label="trajectory", title="x-y plane view", 
    #     xlabel="x [-]", ylabel="y [-]", subplot=2, aspect_ratio = 1)
    
    # # adding the moon
    # scatter!([1-mu_earthmoon],[0], color="gray", markersize = 4, label="m2: moon", subplot=2)

end

