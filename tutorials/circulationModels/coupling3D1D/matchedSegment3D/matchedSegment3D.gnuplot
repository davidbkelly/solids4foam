# B5, the matched segment, against its exact solution: the pressures at the
# pipe's two ends and the reflection at the proximal midpoint, each divided by
# the incident peak. Exact: the incident pulse Zc Qin(t - x/c0), with no
# reflection, delayed by 0.2/c0 = 32.4 ms through the pipe.

run = "< cat postProcessing/circulation1DCoupling/*/coupling.dat | grep -v '^#' | sort -g -k1,1"
c0 = 6.17213399848
P = 1050*c0/(pi*1e-4)*1e-9
s = 1/sqrt(2e4)
t1 = 0.05 + 1/c0
t2 = t1 + 0.2/c0
g(t, t0) = exp(-(t - t0)**2/(2*s*s))

do for [term in "pngcairo pdfcairo"] {
    if (term eq "pngcairo") {
        set terminal pngcairo size 1200,450 font ",11"
        set output "matchedSegment3D.png"
    } else {
        set terminal pdfcairo size 10in,3.75in font ",10"
        set output "matchedSegment3D.pdf"
    }
    set multiplot layout 1,2
    set grid
    set xrange [0.15:0.35]
    set xlabel "t [s]"
    set key top right
    set samples 1000

    set title "Pressure at the ends of the 3D pipe"
    set ylabel "P / incident peak"
    plot run u 1:($4/P) w l lw 2 lc rgb "#1f77b4" t "inlet", \
         run u 1:($6/P) w l lw 2 lc rgb "#d62728" t "outlet", \
         g(x, t1) w l dt 2 lw 2 lc rgb "black" t "exact", \
         g(x, t2) w l dt 2 lw 2 lc rgb "black" notitle

    set title "Proximal midpoint: incident pulse, then any reflection"
    set ylabel "P / incident peak"
    set yrange [-0.01:0.01]
    plot run u 1:($7/P) w l lw 2 lc rgb "#1f77b4" t "run (exact: 0 after the pulse)"
    unset yrange
    unset multiplot
    unset output
}
