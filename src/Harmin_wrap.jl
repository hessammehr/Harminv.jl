libpath = "/home/hessammehr/harminv-1.4.1/.libs"

const C = Complex128
const F = Float64

struct HarminvNative
    c::Ptr{C}
    n::Int
    K::Int
    J::Int
    nfreqs::Int
    fmin::F
    fmax::F
    z::Ptr{C}
    U0::Ptr{C}
    U1::Ptr{C}
    G0::Ptr{C}
    G0_M::Ptr{C}
    D0::Ptr{C}
    B::Ptr{C}
    u::Ptr{C}
    amps::Ptr{C}
    errs::Ptr{F}
end

if !(libpath in Base.LD_LOAD_PATH)
    push!(Base.LD_LOAD_PATH, libpath)
end

ccall((:harminv_data_create, "libharminv"), )