# This reporduces the example shown in Fig. 3.1 of "Vandermonde with Arnoldi" and it
# turns out that the Floater-Hormann rational interpolant does OK here too.
using VanderRat
using PyPlot
using BaryRational

xx = [-1.0:0.001:-0.3; 0.3:0.001:1.0;];
yy = [fill(-1.0, 701); fill(1.0, 701)];

d,  H  = polyfitA(xx, yy, 80)
dd, Hd = polyfitD(xx, yy, 80)

s  = [-1.0:0.01:1.0;]
w  = polyvalA(d, H, s)
wd = polyvalD(dd, Hd, s)

plot(s, w,  label="Vandermonde-A")
plot(s, wd, label="Vandermonde-D")

# These turn out to be very close
fh = FHInterp(xx, yy, order=5)
wfh = fh.(s)

plot(s, wfh, label="Floater")

legend()
title("Example 3.1 from 'Vandermonde with Arnoldi'")
