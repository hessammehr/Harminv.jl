module Harminv

const SMALL = 1e-8
const NPOW = 8
const UNITY_THRESH = 1e-4
const SINGULAR_THRESHOLD = 1e-5

struct HarminvData{C, F} where C, F
    c :: Vector{C}
    n :: Int
    K :: Int
    J :: Int
    nfreqs :: Int
    fmin :: F
    fmax :: F
    z :: Vector{C}
    U0 :: Vector{C}
    U1 :: Vector{C}
    G0 :: Vector{C}
    G0_M :: Vector{C}
    D0 :: Vector{C}
    B :: Vector{C}
    u :: Vector{C}
    amps :: Vector{C}
    errs :: Vector{F}
end

function HarminvData(signal::Vector{C}, fmin, fmax, nf) where C
    @assert nf > 1
    @assert n > 0
    @assert !isempty(signal)
    @assert fmin < fmax

    J = nf
    n = length(signal)
    z = @. exp(-im * 2π * fmin + (fmin-fmax)/(nf-1) * (0:nf-1))

    K = div(n,2) - 1
    U0 = Array{C}(J, J)
    U1 = Array{C}(J, J)
    amps = Array{C}(J)
    errs = Array{Float64}(J)
    G0, G0_M, D0 = generate_U(U0, U0, 0, signal, n, K, J, J, z, z, nothing, nothing, nothing)
    Harminv(signal, n, K, nf, -1, fmin, fmax, z, Matrix{C}(nf,nf), Array{C}(nf,nf), C[],
            C[], C[], C[], C[], C[], C[], Float64[])
    
end

function generate_U(U::AbstractMatrix{C}, U1, p, c::AbstractArray{C},
                    n, K, J, J2, z::AbstractArray{C}, z2::AbstractArray{C},
                    G0 , G0_M, D0) where C
    M = h.K - 1
    @assert h.n >= 2*h.K + p

    if G0 != nothing
        G = G0[:]
        G_M = G0_M[:]
        D = D0[:]
    else
        G = zeros(C, J)
        G_M = zeros(C, J)
        D = zeros(C, J)
    end

    z_m = ones(C, J)
    @. z_inv = 1.0 / z
    @. z_M = z^-M
    if z !== z2
        @. z2_inv = 1 / z2
        z2_m = ones(C, J2)
        @. z2_M = z2^-M
        G2 = zeros(C, J2)
        G2_M = zeros(C, J2)
    else
        # ??
    end

    for m = 1:M
        c1 = c[m + p]
        c2 = c[m + p + M + 1]
        d = float(m)
        d2 = float(M - m - 1)
        if G0 == nothing
            for i = 1:J
                x1 = z_m[i] * c1
                x2 = z_m[i] * c2
                G[i] += x1
                G_M[i] += x2
                D[i] += x1*d + x2*d2*z_M[i] * z_inv[i]
                if m % NPOW == 0
                    z_m[i] = z_inv[i]^m
                else
                    z_m[i] *= z_inv[i]
                end
            end
        end
        if z !== z2
            for i = 1:J2
                G2[i] = z2_m[i] * c1
                G2_M[i] = z2_m[i] * c2
                if m % NPOW == 0
                    z2_m[i] = z2_inv[i]^m
                else
                    z2_m[i] *= z2_inv[i]
                end
            end
        end
    end
    
    if z !== z2
        for i = 1:J
            for j = 1:J2
                if abs(z[i] - z2[j]) < SMALL
                    U[i,j] = D[i]
                else
                    U[i,j] = ( z[i] * G2[j] - z2[j] * G[i] +
                                        z2_M[j] * G_M[i] - z_M[i] * G2_M[j] ) /
                                      (z[i] - z2[j])
                end
            end
        end
        if U1 !== nothing
            for i = 1:J
                for j = 1:J2
                    if abs(z[i] - z2[j]) < SMALL
                        U1[i,j] = z[i] * (D[i] - G[i]) +
                                           z_M[i] * G_M[i]
                    else
                        U1[i,j] = ( z[i] * z2[j] * (G2[j] - G[i]) +
                                             z2_M[j] * z[i] * G_M[i] -
                                             z_M[i] * z2[j] * G2_M[j] ) /
                                           (z[i] - z2[j])
                    end
                end
            end
        end
    else
        for i = 1:J2
            U[i,i] = D[i]
            for j = (i+2):J
                U[i,j] = ( z[i] * G[j] - z[j] * G[i] +
                           z_M[j] * G_M[i] - z_M[i] * G_M[j] ) /
                         (z[i] - z[j])
            end
        end
        if U1 !== nothing
            for i = 1:J2
                U1[i,i] = z[i] * (D[i] - G[i]) + z_M[i] * G_M[i]
                for j = (i+2):J
                    U1[i,j] = ( z[i] * z[j] * (G[j] - G[i]) +
                                z_M[j] * z[i] * G_M[i] -
                                z_M[i] * z[j] * G_M[j] ) /
                              (z[i] - z[j])
                end
            end
        end
        if z === z2
            U .= Symmetric(U, :U)
            if U1 !== nothing
                U1 .= Symmetric(U1, :U)
            end
        end
    end

    return G0, G0_M, D0
end

function symmetric_normalize!(vect)
    symm_norm = sqrt(sum(vect.*vect))
    vect ./= symm_norm
end

function solve_eigenvects(A0)
    v0, mm = eigs(d.U0)
    # symmetric-normalized eigenvectors
    V0 = [symmetric_normalize(V[:,i] for i=1:size(mm)[2]]
    v0, V0
end

function solve_once(d::HarminvData{C,F}) where {C,F}
    one = 1
    zone = one(C)
    zzero = zero(C)
    max_v0 = 0.0

    J = d.J
    v0, V0 = solve_eigenvects(d.U0)
    max_v0 = maximum(v0)
    
    d.nfreqs = J2
    for i = 1:J

    end
end

function solve_again(d::HarminvData)

end

# TODO: actually use ok
function solve_ok_modes!(d::HarminvData, ok::Function, ok_d)
    harminv_solve_once(d)
    cur_nf = d.nfreqs
    while true
        prev_nf = cur_nf
        solve_again(d, ok, ok_d)
        cur_nf = d.nfreqs
        if cur_nf >= prev_nf
            break
        end
    end

end

solve!(d::HarminvData) = solve_ok_modes!(d, nothing, nothing)

frequency_errors!(d::HarminvData)

function u_near_unity(u, n) 
    nlgabsu = n * log(abs(u))
    nlgabsu > log(UNITY_THRESH) && nlgabsu < -log(UNITY_THRESH)
end

function compute_amplitudes(d::HarminvData{C,F}) where {C,F}
    a = Array{C}(h.nfreqs)
    u = Array{C}(h.nfreqs)
    ku = 1

    for k = 1:h.nfreqs
        if u_near_unity(d.u[k], d.n)
            u[ku] = d.u[k]
            ku += 1
        end
    end

    nu = ku - 1
    Uu = Matrix{C}(d.J, nu)
    generate_U(Uu, nothing, 0, d.c, d.n, d.K, d.J, nu, d.z, u)

    for k = 1:d.nfreqs
        asum = zero(C)
        if u_near_unity(d.u[k], d.n)
            for j = 1:d.J
                asum += d.B[k,j] * Uu[j, ku]
            end
            asum /= d.K
            ku += 1
        else
            for j = 1:d.J
                asum += d.B[k, j] * d.G0[j]
            end
        end
        a[k] = asum^2
    end
    return a
end

end