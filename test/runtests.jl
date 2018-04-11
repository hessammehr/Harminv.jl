using Base.Test

t = 0:0.1:500
x = @. exp(2π*(-0.001+0.05im)*t) #+ exp(2π*(-0.002+0.03im)*t)
plot(t,real.(x))

import Harminv
reload("Harminv")

d = Harminv.HarminvData(x, 0.0, 0.05, 100)
Harminv.solve_once!(d)
Harminv.solve_again!(d)
Harminv.compute_amplitudes(d)
