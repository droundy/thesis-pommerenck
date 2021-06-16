## Response to referees

Dear editor,

The responses of referees 1 and 3 could hardly be more divergent.  Our work is impactful precisely because it can tell researchers when to stop looking... or perhaps more relevantly it can tell funding agencies when to stop funding research in this direction.  This is threatening, and Referee 1 doesn't seem interested in knowing whether we should stop looking at porous materials for natural gas and hydrogen storage.

# Referee 1

Referee 1 states that "Phenomenological predictions without any relation to real
materials are not very helpful, if not useless."  The referee fails to grasp
that our theory *does* take into account the gas, and places a fundamental limit
on the capacity of an entire class of materials (inflexible porous materials).
This is far from useless.  The referee suggests as superior an approach
involving screening hundreds of thousands of actual MOFs.  We will point out
that such an approach cannot possibly determine that there isn't a better MOF
out there, while our approach can do so.

TODO: We should find a place to edit the paper to make this more clear, and
state here that we did so.

The referee gives a litany of structural properties that we do not address in
our model: "pore volume, pore size, specific surface area, crystal density,
skeletal density, packing density, different adsorption sites etc."

    We note that some of these (e.g. crystal density, skeletal density)
    manifestly have no impact on the property we are studying, which is
    volumetric storage capacity.  Others of course do impact the volumetric
    storage capacity of a material, but do not relate to an upper bound
    specifically.  TODO: We do add a discussion of the pore volume fraction
    and how it can constrained by our upper bound, to make more clear the
    range of applications of our bound.  TODO: We can put a lower bound on
    pore volume fraction contingent on the maximum density having a given value.
    It will take some explanation, but I think is worthwhile.

The other suggestion of Referee 1 is that we clearly state the capacity that we report, as discussed by Parilla et al.

   The deliverable capacity is defined in the paper as "the density of the gas
   in the material at the storage pressure p_full minus the residual gas that
   remains adsorbed at the lowest pressure p_empty such that sufficient flow is
   maintained to feed the engine."  We believe this is clear. The various
   capacities of Parilla et al. are geared towards experimental measurements and
   account for experimental uncertainty in the volume of material, an issue
   which is not present in a theoretical paper.  We note that the paper Referee 1 brought up as an example does not any more clearly define "usable capacity" as they refer to what we call "deliverable capacity".  They do cite Parilla et al., but only in the context of describing their experimental equipment.

# Referee 3

Referee 3 states that the conclusions are convincing, and makes two suggestions
for improving the paper:

1. The referee suggests that we attempt to quantify the possible effect of a
   flexible MOF on our upper bound.

We have added a section "Flexible two-phase model" discussing flexible MOFs. We
cannot place a formal upper bound on a flexible MOF in general, but we can put
an upper bound on the capacity of a two-phase MOF such as Co(bdp), if we
constrain the free energy difference between the two phases. This suggests that
this type of a system will require a more stable closed phase in order to reach
the DOE targets.

1. The referee asks that we provide some guidelines in terms of working
   conditions to reach the DOE target.  They further suggest looking at a
   temperature difference to see how that could impact the deliverable capacity.

TODO: Add *some* sort of extra discussion on the topic.  Our model cannot
predict a temperature+pressure swing cycle because that would require
information about the entropy of adsorption at the two temperatures, which is
rather trickier, I think.

