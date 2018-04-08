module Harminv

const C = Sys.WORD_SIZE == 64 ? Complex128 : Complex64
const F = Sys.WORD_SIZE == 64 ? Float64 : Float32

struct HarminvData
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

function HarminvData(J, z::AbstractArray{Complex64})
    
end

function generate_U(u :: AbstractArray{C}, u1 :: AbstractArray{C}, p :: Integer)
    M = h.K - 1
    @assert h.n >= 2*h.K + p
    @assert 
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