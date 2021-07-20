# monomers - OPTIONS
# anti-apoptotic, GUARDIANS
# BCL-2, BCL-XL, BCL-W, MCL-1, BFL-1/A1


# monomers, - OPTIONS
# PRO-apoptotic (pore-formers)
# BAX, BAK, BOK

# monomers - OPTIONS
# PRO - apoptotic BH3-only proteins (activating others)
# BAD, BID, BIK, BIM, BMF, HRK, NOXA, PUMA, etc

# las reglas
## COMPARTMENT - MOM

# skip BH3 domain for now
# The BH3 domain of activator BH3-
# only proteins binds to the BH3 domain-binding groove in BAX/
# BAK.11,12
# PRO - apoptotic BH3-only proteins - BINDS - BAX &/OR BAK   (figure 2A - in MOM)

# (1)
# BAX, BAK activated (when bound) (binding--transient interaction)
#   ^ gain function as result being bound
#   ^ function is pro-poreforming gain function
# conformational changes that result in BAX/BAK homooligo -- part of activation step
# BAX:(PRO - apoptotic BH3-only proteins) --> (?) BAX-poreforming:(PRO - apoptotic BH3-only proteins)
## --- Binding to tBID activates BAX to insert in the bilayer



# (c) tBID helps BAX more function,
# feed back regulation -- active form provokes inhibition of the activation
# BAX binds BCL-XL, then allows come to MOM, then once in membrane, prevents other BAXs getting actvated



# BAX-poreforming
# oligomerization (lookup/assume how many idividuals) -- grouping together
#
# n * (BAX-poreforming) <--> BAX-pore

# (3) ---- DONE IN NEXT_DRAFT file
# ANTI - apoptotic proteins - BINDS - BAX &/OR BAK
# BAX - held up by anti-apop

# (4) ---- DONE IN NEXT_DRAFT file
# pro-apop = sensitizers
# PRO- BH3-only proteins - bind- ANTI-apop --> inactive ANTI-apop
#   ^^^^^^ holding up anti-apop
# -- inhibiting inhibitors, which is ? true pro-apoptotic --- which one is pushing the whole system


## (5)
# constitutive -- never stops (BAK w/ other things; in MOM, able to bind other proteins; it's ready)
# BAX, BID, BIN -- dynamic, need prep-step
# binding -- balance b/w association & disassociation
# once bound, hard to bring apart
### QUESTION -- happens BAX:tBID in cytosol or membrane, or both
# BAX + tBID <--> BAX:tBID   (recruitment of complexes, assembly, to membrane -- indiv bind then form)
# BAX:tBID (in membrane) --> BAX:tBID (stabilized)
# BAX:tBID (stabilized) <--> BAX + tBID      --- dissociation constant close to zero, hence "stable"
##    ^^ hardly ever happens

# BCL-X binds to BAX,
# BCL-x:BAX (cyto) <--> BCL-x:BAX (mom)   (forward rate is low)
# BAX (cyto) <--> BAX (mom)   (forward rate is high)


# BAX + Ub <--> BAX:Ub (tagging) -- noncovalent
# BAX:Ub --> BAX:Ub (tagged) -- covalent

# mediated -- catalysis -- caspase8 is doing something
## COMPARTMENT - cytosol (C-8 lives there)
# C-8 + BID <--> C-8:BID
# c-8:BID --> C-8 + BID(cleaved - p15 --- called tBID) + BID(cleaved - p7)
# tBID (cytosol) <--> tBID (MOM)



##### read over pseudocode, find what not sound clear
## point out areas that sound ambiguous --- address ambiguity first, then write rules
## heatmap for parameters
###### next citations -- will have data
# rules
# parameters (to use table 1A for the priors)
# use citations for data
# ... each protein individually. see if in reviews for that specific protein; make sure not miss.

# ... likely get up to parameter optimization
# MI -- working on using cell viability data ...
## can make recommendations, this model more reasonable for ambiguous



# proteins in gel, elecric charge, drift to other end of gel. drift slower if heavier
# 15 or 7 band ... p15
# progammming

