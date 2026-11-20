# The rigid segment against its closed form: the flow through the pipe and
# the pressures at its two ends, each divided by the incident pulse's peak.
# The run's inflow is the closed form's times 1e-3: its incident peaks are
# Q = 1e-9 m3/s and P = 0.0206288 Pa, the closed form's 1e-6 and 20.6288.

run   = "postProcessing/circulation1DCoupling/0/coupling.dat"
# the closed form is CSV with two header lines: drop them, commas to spaces
exact = "< grep -v '^[#t]' reference/rigidSegment_pipe.csv | tr ',' ' '"
Qr = 1e-9;  Pr = 0.0206288
Qe = 1e-6;  Pe = 20.6288

do for [term in "pngcairo pdfcairo"] {
    if (term eq "pngcairo") {
        set terminal pngcairo size 1200,450 font ",11"
        set output "rigidSegment3D.png"
    } else {
        set terminal pdfcairo size 10in,3.75in font ",10"
        set output "rigidSegment3D.pdf"
    }
    set multiplot layout 1,2
    set grid
    set xrange [0.15:0.30]
    set xlabel "t [s]"
    set key top right

    set title "Flow through the 3D pipe"
    set ylabel "Q / incident peak"
    plot run u 1:($3/Qr) w l lw 2 lc rgb "#1f77b4" t "pimpleFoam + 1D", \
         exact u 1:($2/Qe) w l dt 2 lw 2 lc rgb "black" t "closed form"

    set title "Pressure at the ends of the pipe"
    set ylabel "P / incident peak"
    plot run u 1:($4/Pr) w l lw 2 lc rgb "#1f77b4" t "inlet", \
         run u 1:($6/Pr) w l lw 2 lc rgb "#d62728" t "outlet", \
         exact u 1:($3/Pe) w l dt 2 lw 2 lc rgb "black" t "closed form", \
         exact u 1:($4/Pe) w l dt 2 lw 2 lc rgb "black" notitle
    unset multiplot
    unset output
}
