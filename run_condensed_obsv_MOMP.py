from condensed_test_mito import model as md
from pysb.simulator import ScipyOdeSimulator
import pylab as pl

t = pl.linspace(0, 20000)  # num seconds for simulation
sim_results = ScipyOdeSimulator(md, tspan=t).run()
y_out = sim_results.all



# check y_out, should be observables
print(f"Observables-- y_out: \n{y_out}")

# chck content for 'Bid' observable
## bond num 9
print(f"Check BIDp15_bound_BAX Content::\n{y_out['obBID_p15_bound_BAX']}")
## bond num 10
print(f"Check BIDp15_bound_BAX_pore Content::\n{y_out['obBID_p15_bound_BAX_pore']}")

# plot data
pl.ion()
pl.figure()
pl.plot(t, y_out['obBID_p15_bound_BAX'], label="BIDp15_activated_BAX")
pl.plot(t, y_out['obBID_p15_bound_BAX_pore'], label="formed_BAX_pore")
pl.legend()
pl.xlabel("Time (s)")
pl.ylabel("Molecules/ Cell")
pl.show()


monom = md.monomers
print(f"Monomers of the model::\n{monom}")
print(f"\nList comprehension for model names in monomers:\n{[m.name for m in monom]}")
print(f"\nKeys for Monomers: {monom.keys()}")

print(f"\n\n'PoreFormers' Monomer: {monom['PoreFormers']}")
print(f"'PoreFormers' sites::\n{monom['PoreFormers'].sites}")

print(f"\n\n'BCLs' Monomer: {monom['BCLs']}")
print(f"'BCLs' sites::\n{monom['BCLs'].sites}")

prams = md.parameters
print(f"\nExplore PARAMETERS::\n{prams}")
print(f"\nValue of 'kf' param:: {prams['kf'].value}")


print(f"\n\nModel SPECIES::\n")
for sp in md.species:
    print(f"\n{sp}")