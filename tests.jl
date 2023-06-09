const sim = Sim(400);
const g = sim.g;
const al1 = g.adjlist[1];
const at1 = g.adjtup[1];

# Copy these functions to repl, otherwise the testrun might hang on the sleep function

# Compile
testrun(sim, 0.1, s = :adjlist)
testrun(sim, 0.1, s = :adjtup)

testrun(sim, 2, s = :adjlist)
testrun(sim, 2, s = :adjlist)

# Now both are much faster
println("Making deepcopies of adjlist and adjtup")
g.adjlist .= deepcopy(g.adjlist)
g.adjtup .= deepcopy(g.adjtup)

testrun(sim, 2, s = :adjlist)
testrun(sim, 2, s = :adjtup)

# Dont understand the interpolation restuls
# Now al is super slow even though it is a const
# Interpolating g.state makes some slower, but @benchmark getEnergyFactor($g.state, $al1) is faster?
println("Benchmark getting the energy factor for al1: ")
display(@benchmark getEnergyFactor(g.state, al1))
println("Benchmark getting the energy factor for at1: ")
display(@benchmark getEnergyFactor(g.state, at1))

# Now it is faster
println("Benchmark getting the energy factor for al1 after deepcopy: ")
display(@benchmark getEnergyFactor(g.state, al1))

println("Benchmark getting the energy factor for at1 after deepcopy: ")
display(@benchmark getEnergyFactor(g.state, at1))
