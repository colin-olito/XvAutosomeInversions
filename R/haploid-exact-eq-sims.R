for(j in 1:20){
#timescale of simulation
generations = 10000

#parameters
s = 0.01*j
m = 0.01

#initial haplotype frequencies
x1 = 0.25
x2 = 0.25
x3 = 0.25
x4 = 0.25

for(i in 1:generations){
  W = x1 + x2*(1 + s) + x3*(1 + s) + x4*(1 + s)^2
  D.adults = ((1 - m)*(1 + s)^2)*((x1*x4 - x2*x3)*(1 - m)/W + m*x4)/W
  #frequency next generation
  y1 = x1*(1 - m)/W + m - D.adults/2
  y2 = x2*(1 + s)*(1 - m)/W + D.adults/2
  y3 = x3*(1 + s)*(1 - m)/W + D.adults/2
  y4 = x4*((1 - m)*(1 + s)^2)/W - D.adults/2
  x1 = y1
  x2 = y2
  x3 = y3
  x4 = y4
}

#approx.a = m/s
approx.a = ((s + m*(1 + s)) - sqrt((s - m*(1 + s))^2))/(2*s)
freq.a = x1 + x3
freq.b = x1 + x2
LD = (x1*x4 - x2*x3)

print(c(s, approx.a, freq.a, freq.a, LD))

}

