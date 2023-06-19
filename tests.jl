# NN = 1; 8 neighbors per vertex
const sim = Sim(400);
println("NN = 1; 8 neighbors per vertex")


#Precompile
testrun(sim, :adjtup, 0.1, print = false)
testrun(sim, :adjlist, 0.1, print = false)

println("Benchmark getEnergyFactor")
display(@benchmark getEnergyFactor($sim.g.state, $sim.g.adjtup[1]))     # mean: ~5.287 ns
display(@benchmark getEnergyFactor($sim.g.state, $sim.g.adjlist[1]))    # mean: ~4.003 ns
                                                                            # Rel performance adjtup/adjlist: 0.76
testrun(sim, :adjtup)   # ~23370190 updates
testrun(sim, :adjlist)  # ~25054611 updates
                            # Rel performance adjtup/adjlist: 0.93

println("Relocating")
relocate!(sim)
testrun(sim, :adjtup)   # ~23436213 updates
testrun(sim, :adjlist)  # ~18523010 updates
                            # Rel performance adjtup/adjlist: 1.26

println("Localizing")
localize!(sim)
testrun(sim, :adjtup)   # ~40918586 updates
testrun(sim, :adjlist)  # ~32039635 updates
                            # Rel performance adjtup/adjlist: 1.28

# NN = 3; 48 neighbors per vertex
const sim2 = Sim(400, 3)
println("NN = 3; 48 neighbors per vertex")

println("Benchmark getEnergyFactor")
display(@benchmark getEnergyFactor($sim2.g.state, $sim2.g.adjtup[1]))     # mean: ~19.817 ns
display(@benchmark getEnergyFactor($sim2.g.state, $sim2.g.adjlist[1]))    # mean: ~14.660 ns
                                                                            # Rel performance adjtup/adjlist: 0.74

testrun(sim2, :adjtup)  # ~7582301 updates
testrun(sim2, :adjlist) # ~7331573 updates
                            # Rel performance adjtup/adjlist: 1.03

println("Relocating")
relocate!(sim2)
testrun(sim2, :adjtup)  # ~7369357 updates
testrun(sim2, :adjlist) # ~7179336 updates
                            # Rel performance adjtup/adjlist: 1.03

println("Localizing")
localize!(sim2)
testrun(sim2, :adjtup)  # ~8175987 updates
testrun(sim2, :adjlist) # ~10761614 updates
                            # Rel performance adjtup/adjlist: 0.76

# NN = 10; 440 neighbors per vertex
const sim3 = Sim(400, 10)
println("NN = 10; 440 neighbors per vertex")

println("Benchmark getEnergyFactor")
display(@benchmark getEnergyFactor($sim3.g.state, $sim3.g.adjtup[1]))     # mean: ~233.333 ns
display(@benchmark getEnergyFactor($sim3.g.state, $sim3.g.adjlist[1]))    # mean: ~115.628 ns
                                                                            # Rel performance adjtup/adjlist: 0.50

testrun(sim3, :adjtup)  # ~2732528 updates
testrun(sim3, :adjlist) # ~2975036 updates
                            # Rel performance adjtup/adjlist: 0.92

println("Relocating")
relocate!(sim3)
testrun(sim3, :adjtup)  # ~2854341 updates
testrun(sim3, :adjlist) # ~2845265 updates
                            # Rel performance adjtup/adjlist: 1.00

println("Localizing")
localize!(sim3)
testrun(sim3, :adjtup)  # ~3303150 updates
testrun(sim3, :adjlist) # ~3532045 updates
                            # Rel performance adjtup/adjlist: 0.94
