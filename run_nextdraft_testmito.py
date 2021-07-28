from next_draft_test_mito import model

monom = model.monomers
print(f"Monomers of the model::\n{monom}")
print(f"\nList comprehension for model names in monomers:\n{[m.name for m in monom]}")
print(f"\nKeys for Monomers: {monom.keys()}")

print(f"\n\n'PoreFormers' Monomer: {monom['PoreFormers']}")
print(f"'PoreFormers' sites::\n{monom['PoreFormers'].sites}")

print(f"\n\n'BCLs' Monomer: {monom['BCLs']}")
print(f"'BCLs' sites::\n{monom['BCLs'].sites}")

prams = model.parameters
print(f"\nExplore PARAMETERS::\n{prams}")
print(f"\nValue of 'kf' param:: {prams['kf'].value}")