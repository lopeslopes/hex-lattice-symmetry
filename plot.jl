using Plotly
plotlyjs()

f1a = open("lattice1A.dat", "r")
x1a = []
y1a = []
while ! eof(f1a)
    s1 = readline(f1a)
    s2 = split(s1, ";")
    append!(x1a, parse(Float64,s2[1])) 
    append!(y1a, parse(Float64,s2[2]))
end

p1 = scatter(x1a, y1a, label="1A")
gui(p1)
