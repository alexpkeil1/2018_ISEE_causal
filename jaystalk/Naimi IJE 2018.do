*** Code for Naimi Paper ****************

input a0 z a1 y n
0 0 0  87.29 209271
0 0 1 112.11  93779
0 1 0 119.65  60654
0 1 1 144.84 136293 
1 0 0 105.28 134781
1 0 1 130.18  60789
1 1 0 137.72  93903
1 1 1 162.83 210527
end

expand n

* standard models
gen avga = (a0+a1)/2
reg y avga, nohead cformat(%6.2f)
reg y avga z, nohead cformat(%6.2f)
reg y a0, nohead cformat(%6.2f)
reg y a0 z, nohead cformat(%6.2f)
reg y a1, nohead cformat(%6.2f)
reg y a1 z, nohead cformat(%6.2f)

* g-formula
quietly sum y if a1==0 & z==0 & a0==0
scalar y000 = r(mean)
quietly sum y if a1==0 & z==0 & a0==1	
scalar y001 = r(mean)
quietly sum y if a1==0 & z==1 & a0==0	
scalar y010 = r(mean)
quietly sum y if a1==1 & z==0 & a0==0
scalar y100 = r(mean)
quietly sum y if a1==0 & z==1 & a0==1
scalar y011 = r(mean)
quietly sum y if a1==1 & z==1 & a0==0
scalar y110 = r(mean)
quietly sum y if a1==1 & z==0 & a0==1
scalar y101 = r(mean)
quietly sum y if a1==1 & z==1 & a0==1
scalar y111 = r(mean)
quietly sum z if a0==1
scalar z1 = r(mean)
quietly sum z if a0==0
scalar z0 = r(mean)

scalar Y00 = (y010*z0)+(y000*(1-z0))
scalar Y11 = (y111*z1)+(y101*(1-z1))
scalar Y01 = (y011*z1)+(y001*(1-z1))
scalar Y10 = (y110*z0)+(y100*(1-z0))

disp "Y00 = ", Y00
disp "Y01 = ", Y01
disp "Y10 = ", Y10
disp "Y11 = ", Y11

* IPTW
quietly sum a1
scalar pa1 = r(mean)
quietly sum a1 if z==0
scalar pa1z0 = r(mean)
quietly sum a1 if z==1
scalar pa1z1 = r(mean)

scalar w000 = (0.5 * (1-pa1z0))
scalar w001 = (0.5 * (pa1z0))
scalar w010 = (0.5 * (1-pa1z1))
scalar w011 = (0.5 * (pa1z1))
scalar w100 = (0.5 * (1-pa1z0))
scalar w101 = (0.5 * (pa1z0))
scalar w110 = (0.5 * (1-pa1z1))
scalar w111 = (0.5 * (pa1z1))

scalar sw000 = (0.5*(1-pa1))/(0.5 * (1-pa1z0))
scalar sw001 = (0.5*pa1)/(0.5 * (pa1z0))
scalar sw010 = (0.5*(1-pa1))/(0.5 * (1-pa1z1))
scalar sw011 = (0.5*pa1)/(0.5 * (pa1z1))
scalar sw100 = (0.5*(1-pa1))/(0.5 * (1-pa1z0))
scalar sw101 = (0.5*pa1)/(0.5 * (pa1z0))
scalar sw110 = (0.5*(1-pa1))/(0.5 * (1-pa1z1))
scalar sw111 = (0.5*pa1)/(0.5 * (pa1z1))

gen sw = .
replace sw = (0.5*(1-pa1))/(0.5 * (1-pa1z0)) if a0==0 & z==0 & a1==0
replace sw = (0.5*pa1)/(0.5 * (pa1z0)) if a0==0 & z==0 & a1==1
replace sw = (0.5*(1-pa1))/(0.5 * (1-pa1z1)) if a0==0 & z==1 & a1==0
replace sw = (0.5*pa1)/(0.5 * (pa1z1)) if a0==0 & z==1 & a1==1
replace sw = (0.5*(1-pa1))/(0.5 * (1-pa1z0)) if a0==1 & z==0 & a1==0
replace sw = (0.5*pa1)/(0.5 * (pa1z0)) if a0==1 & z==0 & a1==1
replace sw = (0.5*(1-pa1))/(0.5 * (1-pa1z1)) if a0==1 & z==1 & a1==0
replace sw = (0.5*pa1)/(0.5 * (pa1z1)) if a0==1 & z==1 & a1==1

reg y i.a0##i.a1 [pw=sw], nohead cformat(%6.2f)

table a1 a0 z [pweight = sw], contents(freq )







