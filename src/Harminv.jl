module Harminv

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

    n = length(signal)
    z = @. exp(-im * 2π * fmin + (fmin-fmax)/(nf-1) * (0:nf-1))

    K = div(n,2) - 1
    Harminv(signal, n, K, nf, -1, fmin, fmax, z, Matrix{C}(nf,nf), Array{C}(nf,nf), C[],
            C[], C[], C[], C[], C[], C[], Float64[])
    
end

function generate_U(u::AbstractArray{C}, u1::AbstractArray{C}, p, c::AbstractArray{C},
                    n, K, J, J2, z::AbstractArray{C}, z2::AbstractArray{C},
                    G = zeros(J), G_M = zeros(J), D = zeros(J)) where C
    M = h.K - 1
    @assert h.n >= 2*h.K + p
    z_inv = Array{C}(J)
    z_m = Array{C}(J)
    z_M = Array{C}(J)
end

function init_z(h :: HarminvData, J :: Integer, z::AbstractArray{Complex64})
    d.J = J
    d.z = z
    d.U0 = 
end

function compute_amplitudes(h :: HarminvData)
    a = Array{Complex128}(h.nfreqs)
    u = Array{Complex128}(h.nfreqs)
    ku = 1

    for k = 1:h.nfreqs
        if u_near_unity(d.u[k], d.n)
            u[ku] = d.u[k]
        end
    end
    nu = ku

    generate_U(d, u)


end

end