# KALE 2018 dance of death review


# monomers - OPTIONS - called "--L"
# anti-apoptotic, GUARDIANS
# BCL-2, BCL-XL, BCL-W, MCL-1, BFL-1/A1


# monomers, - OPTIONS - called "B++"
# PRO-apoptotic (pore-formers)
# BAX, BAK, BOK

# monomers - OPTIONS - called "BH3+"
# PRO - apoptotic BH3-only proteins (activating others)
# BAD, BID, BIK, BIM, BMF, HRK, NOXA, PUMA, etc
####### start FIGURE-2- PARTE A ######
######################################
##### (A1)
## "All binding interactions are reversible and equilibria are governed by local affinities.
## Interactions with the lipid bilayer change the affinities of the interactions and
### therefore have an active role in the functions of the proteins.
## THE GOOD STUFF -- Binding of the activator BH3-only proteins (e.g. BID, BIM) to membranes
## increases their affinity for the pore-formers (e.g. BAX, BAK), which are activated (arrows)
## to permeabilize the mitochondrial outer membrane."
# (A1)(a) pro-apoptotic (BH3-only), the activators, bind to membranes
#     (b) when BH3+ are bound to membrane, they are more likely to be bound to pore-formers (B++) BAX, BAK, BOK
#     (c) when B++ are bound to BH3+, B++ gain a MOM-permeability function, to also make pores
## COMPARTMENT - membrane
#  BH3+ (free)  <-->  BH3+:MOM
#  BH3+:MOM + B++  <-->   BH3+:MOM:B++ (high forward rate, noncovalent?)
#  BH3+:MOM:B++  -->  BH3+:MOM:B++  (covalent bond)
#  BH3+:MOM:B++  -->  BH3:MOM:B++ (pore-forming ability - B++)


##### (A2)
## GOOD STUFF - parte 2 - "The anti-apoptotic proteins (e.g. BCL-XL, BCL-2, MCL-1) inhibit both
## the activator BH3-only proteins and ..." by mutual sequestration
# (A2)(a) anti-apoptotic (--L) bind to BH3+
#     (b) --L and BH3+ are held up from further action
#  --L + BH3+  <--> --L:BH3+  (reversible binding, noncovalent?)
#  --L:BH3+   -->   --L:BH3+  (covalent bond)
#  --L:BH3+ (covalent)  --> --L:BH3+  (both held-up)

##### (A3)
## parte 3 - " and the pore-forming proteins by mutual sequestration (T'd arrows)."
# (A3)(a) anti-apoptotic (--L) bind to B++
#     (b) --L and B++ are held up from further action
#  --L + B++  <--> --L:B++  (reversible binding, noncovalent?)
#  --L:B++   -->   --L:B++  (covalent bond)
#  --L:B++ (covalent)  --> --L:B++  (both held-up)

##### (A4)
## parte 4 - "The sensitizer BH3-only proteins (e.g. BAD, NOXA) bind to and inhibit the anti-apoptotic
## proteins also by mutual sequestration."
# (A4)(a) BH3+ bind to --L
#     (b) --L and BH3+ are held up from further action
#  BH3+ + --L   <-->  BH3+:--L    (reversible binding, noncovalent?)
#  BH3+:--L  -->  BH3+:--L  (covalent bond)
#  BH3+:--L (covalent bound) -->  BH3+:--L  (inactive --L, both held-up)

##### (A5)
##### ---------AMBIGUOUS----------- (figure description; maybe explained in article text)
####### -- Increases everyone's affinity for everyone else?? Increases whose affinity for who?
####### -- Increases local concentrations, of WHO, and by HOW MUCH??
####### -- Reduces diffusion of BCL-2 family proteins by how? hydrophobic balances? direct binding of BCL-2?
##
## parte 5 - "Recruitment of the complexes to the membrane by *constitutive* interactions (e.g. BAK)
## [...] increases the affinities and local concentrations and reduces the diffusion of the BCL-2 family proteins."
# - constitutive - never stops, always ready, no prep step
## "complexes" referring to BH3+:--L and --L:B++?? these from A3 & A4
### QUESTION -- happens BAX:tBID in cytosol or membrane, or both
# (A5)(a) complexes are attracted to the membrane by *constitutive* interactions (BAK and others?)
#     (b) makes everyone more attracted to everyone??



##### (A6)
## parte 6 - "Recruitment of the complexes to the membrane by [...] *dynamic* interactions (e.g. BAX, BID, BIM)
## increases the affinities and local concentrations and reduces the diffusion of the BCL-2 family proteins."
# BAX + tBID <--> BAX:tBID   (recruitment of complexes, assembly, to membrane -- indiv bind then form)
# BAX:tBID (in membrane) --> BAX:tBID (stabilized)
# BAX:tBID (stabilized) <--> BAX + tBID      --- dissociation constant close to zero, hence "stable"
##    ^^ hardly ever happens



##### (A7)  ---------AMBIGUOUS-----------
## parte 7 - "Localization at different intracellular membranes also dictates the binding equilibria between
## each family member."
# (A7) Depending on the respective "intracellular membrane" location (which membrane, or, which location within
#      a membrane?), depending, the two (? - ambiguous) species have varying binding equilibria-tendencies.


##### (A8)
## parte 8 - "The efficiency of inhibition by mutual sequestration of anti-apoptotic proteins depends on both
## affinities and off-rates of the interactions."


##### (A9) ----------AMBIGUOUS---------- in figure, maybe more detail in article-text
##### -----BIIIG-- AMBIGOUS-NESSSS------ components of retrotranslocation
### https://github.com/vw-liane/pysb/blob/master/pysb/examples/noo_proteens/bax_hax.py
########################################
## parte 9 - "Interaction of the BH4 region of the anti-apoptotic proteins with BAX shifts the BAX-membrane
##  binding equilibrium to favor the unbound state (retrotranslocation, not shown.)
# (A9)(a) (alone??) BAX is gravitating towards MOM
#     (b) unidentified --L "interact" via BH4 region with BAX
#                 ^^ Ub also plays role? who knows
#     (c) unidentified --L and BAX (and Ub, etc??) form temporary complex
#                   ***what are the bonds that keep this complex together?
#                   ***and what influences these bonds?
#     (d) BAX is shuttled away from MOM
#
#  general rule:: BAX + MOM + --L <--> BAX:MOM + --L   (high reverse rate, thanks to --L)
#                           ^^not sure how --L "interacts"
#  BAX (free) <--> BAX:MOM  (represents the getting closer, attempting to bond)
#  BAX + --L <--> BAX:--L   (temporary shuttling structure, higher forward rate when near MOM)
#                                                           lower reverse rate when away MOM
####### END FIGURE-2- PARTE A ######
####################################



####### start FIGURE-2- PARTE B ######
######################################
##### (B1)
## parte 1 - "BID is activated by caspase-8 mediated cleavage to cBID (cleaved BID) a protein comprise of two fragments
## BID-P7 and BID-P15 held together by hydrophobic interactions."
## COMPARTMENT - cytosol (C8 lives there)
# mediated catalysis - C8 doing something
# (B1)(a) C8 binds to BID
#     (b) BID is activated
#     (c) BID is split into two pieces, p7, p15   (loosely connected by hydrophobic interactions)
#  C8 + BID  <-->  C8:BID  (reversible binding, noncovalent?)
#  C8:BID  -->  C8 + cBID  (temporary bond chops BID into cBID, which is p15 (tBID) and p7 (unnamed))

##### (B2)
## parte 2 - "Rapid high-affinity binding to membranes dissociates the p7 frament to solution and favors insertion of the
## p15 fragment (tBID; truncated BID) into the membrane."
# (B2)(a) BID attempting to bind to MOM.
#     (b) C8 interrupts BID's binding to MOM, and cleaves BID into cBID (two pieces w/ hydrophobic interactios)
#     (c) Yet, cBID dissociates p7 to 'solution' (cytosol?).
#     (d) And cBID likely dissociates p15 (tBID) to insert into MOM.
#  BID (free)  <-->  BID:MOM  (attempting to bind)
##then same above step from (parte 1)
#  C8 + BID  <-->  C8:BID  (reversible binding, noncovalent?)
#  C8:BID  -->  C8 + cBID  (temporary bond chops BID into cBID, which is p15 (tBID) and p7 (unnamed))
# cBID (with rapid above binding) --> p7 (cyto) + p15 (tBID, in MOM)


##### (B3)   ---  same, w/ more detail as Part A1??
#####  ---------AMBIGUOUS-----------
###### -- is BAX still bound to tBID when forming the pore?
###### -- what influences the forward rate of multiple BAX forming a pore?
###### -- how many BAX are required to form a pore?
##              -+- I think it varies by location on the membrane, i.e.
##              it varies by the chemical/ electrical balance of that part of MOM
##
## parte 3 - "Membrane-bound tBID recruits inactive BAX from the cytosol. Binding to tBID activates BAX to insert in the
## bilayer, oligomerize and permeabilize the mitochondrial outer membrane releasing intermembrane space proteins
## including cytochrome c and SMAC."
########### like general BH3 rule from figure 2-a, part 1
# (B3)(a) tBID is in the MOM
#     (b) tBID attracts inactive BAX from cytosol
#     (c) tBID activates BAX, meaning BAX gets pore-forming abilities
#     (d) Unspecified num BAX go into membrane and form a pore
#  tBID (free)  <-->  tBID:MOM
#  tBID:MOM + BAX  <-->   tBID:MOM:BAX (high forward rate, noncovalent?)
#  tBID:MOM:BAX  -->  tBID:MOM:BAX  (covalent bond)
#  tBID:MOM:BAX  -->  tBID:MOM:BAX (pore-forming ability - BAX)
# n * (tBID:MOM:BAX (poreforming))  <-->  BAX-pore    (unknown forward rate)

####### END FIGURE-2- PARTE B ######
####################################



####### start FIGURE-2- PARTE C ######
######################################
##### (C1)
#####  ---------AMBIGUOUS-----------
#####  -- Active(tBID:BAX), which means the pair of tBID bound to BAX, which gave BAX poreforming ability?
#####       (above is how we will write it)
#####  -- Does this mean BCL-XL takes on tBID, thus unactivating BAX?
#####  -- Does this mean BCL-XL also binds to BAX (as photo shows), or does it only bond to tBID?
#####  -- below article shows either could happen; don't know which is most important
# https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0060147
#####  -- tBID:MOM:BAX attracts proximity of BCL-XL, to self

#####  -- What could BCL-XL:BAX (active shape BAX, inactive abilities) do?

## parte 1 - "Active tBID and BAX can recruit BCL-XL to the membrane resulting in inhibition of both pro and anti-
## apoptotic proteins by mutual sequestration.
## -- BCL-XL prevents tBID from activating BAX and prevents BAX from oligomerizing resulting in the inhibition of MOMP.
## -- BAX bound to BCL-XL is in the active (oligomerization competent) conformation."
# (C1)(a) tBID is bound to BAX (poreforming) in the MOM
#     (b) tBID:MOM:BAX attracts (high affinity for?) BCL-XL to membrane also, or is it really attracting to self?
#     (c) BCL-XL enters membrane and ... maybe binds to tBID
#     (d) BCL-XL enters membrane and ... maybe (also or only) binds to BAX
#  tBID:MOM:BAX (poreforming ability-BAX)  <-->  tBID:MOM:BAX:BCL-XL   (high forward rate?)
#                                                   ^^expresses BCL-XL attraction to tBID:MOM:BAX self, but that
#                                                 binding is not completed as a full three-piece group
# (proximity of) tBID:MOM:BAX:BCL-XL  -->  BCL-XL:tBID
# --and/or--
# (proximity of) tBID:MOM:BAX:BCL-XL  -->  BCL-XL:BAX  (**active shape BAX, but inactive abilities**)

####### END FIGURE-2- PARTE C ######
####################################

# "BCL-XL is a defective BAX"
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2858642/


####### start FIGURE-2- PARTE D-E ######
########################################
##### (D1-E1)
#####  ---------AMBIGUOUS-----------
## parte 1 - "BAD inhibits unbound BCL-XLby mutual sequestration.
##        The affinity of BCL-XL is higher for tBID than for BAX (Tables 1A and 1B) therefre,
##            in the absence of other regulatory interactions or PTMS [<-- SO SKEPTICAL]
##            if BCL-XL is bound to tBID and BAX then high concentrations of BAD will displace [INTERESTING]
##            active BAX (d) and then tBID (e) from BCL-XL resulting in MOMP."
# (D1-E1)(a) BCL-XL attracted more to tBID than to BAX
#        (b) "without other interactions happening" [AS IF],
#               then if BCL-XL:tBID:BAX, "high concentrations of BAD" will steal a dance from BAX
#                       BCL-XL:tBID:BAX  + 8(BAD) <-->  BCL-XL:tBID:BAD   (high forward with lots of BAD around)
#        (c) then also kick off tBID
#                BCL-XL:tBID:BAD  <-->  BCL-XL:BAD + tBID + BAX
#        (d) then thus let BAX form pores :-D :-D
#
##  BCL-XL:tBID:BAX + 8(BAD)  <--> BCL-XL:tBID:BAD  + BAX (free to poreisize)
##  -- diferent or at same time, with below? --
##  BCL-XL:tBID:BAD + BAX  <-->  BCL-XL:BAD + BAX (free to poreisize) + tBID

####### END FIGURE-2- PARTE D-E ######
######################################
