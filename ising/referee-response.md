Dear editor,

We appreciate the referees' careful reading of our paper, and feel
that their feedback has led to a significant improvement. We have
addressed all of the first referee's comments by adding two paragraphs
to the introduction, and a section in the paper to 'Address the
general problem of all "flat-histogram" methods'. We also strengthen
the abstract to address the concerns raised by the second referee. We
point out that:

(1) It is relevant to test SAD's convergence properties on the 2D
    Ising where the range of energies is easily known.  This system in
    particular is a kind of worst-case scenario for SAD, and our
    previous paper focused on systems where SAD was more likely to be
    applied in practice.

(2) Our paper is the first of its kind (to the best of our knowledge)
    to compare ``production run'' WL with other recent flat-histogram
    methods.  This is particularly relevant, as it is common in this
    field to discount WL as "not converging to the true answer" (as we
    havepreviously done ourselves, to our chagrin), without
    recognizing that when it is used in practice it is followed by a
    "production run", which has the effect of causing the the
    simulation to converge to the true answer.

Response to first referee
# --------------------------------------------------------------------#

The first referee points out that 'Simulation set-up and analysis
protocols are sound and give confidence in the quality of the
presented data and conclusions.'

The referee also asks that three additions be made to the paper.

(i) "It would be worthwhile to work out this point in the paper
[slower convergence for SAD due to the energy range being an a priori
parameter in WL/MUCA], and, if possible, quantify the loss in
efficiency coming from energy range determination."

We have added a paragraph in the results section that explains the
loss in efficiency SAD suffers by not specifying an energy range a
priori. We explain that the fraction of moves (in the end) that are
outside the range of interesting energies is negligible, and the time
spent before all interesting energies have been identified is not
*much* longer than the time WL takes to find all energies, but that
SAD is more conservative than WL after this point in decreasing
gamma.  We point out that standard WL requires a production run (or
the 1/t variant) in order to prevent *overly* agresssively decreasing
gamma, so SAD and WL are erring in opposite directions.

(ii) "Address the general problem of all "flat histogram" methods that
the performance of these methods is due to hidden barriers almost
always lower than what one would expect from a random walk in energy
(in most cases still exponentially growing). This has been discussed
by Nadler et al (PRE 75 (2007) 1) who also discus a protocol for
optimizing the ensemble which minimizes the residual exponential
slowing down. The authors should point out if and how SAD can be
adapted to such "non-flat-histogram" methods."

We have added a paragraph that addresses the optimized ensemble
approach.  Our understanding is that any of these flat-histogram
methods could be followed by the optimized ensemble method, which
requires as a starting point a set of weights that allow efficient
sampling of the energy range of interest.

(iii) Add a few short remarks on comparison with replica exchange
methods would be helpful.

We added a paragraph very briefly describing how replica exchange has
been adapted to WL, and make the point that the same adaption could
apply to any of these very similar flat-histogram methods.

Response to second referee
# ------------------------------------------------------------------#

We thank the second referee for bringing to our attention two points
that we can better clarify in our paper.

(i) "The work is almost a repetition of what the authors did when they
introduced the algorithm."

Our paper is the first of its kind (to the best of our knowledge) to
compare ``production run'' WL with other flat-histogram methods. The
first (not the second) referee notes rather strongly that this is the
only way that WL should be done; however, no research work has ever
compared this with other flat-histogram methods and certainly not with
the completeness shown here.  This contribution alone makes our paper
valuable, even were SAD not included in the mix.

We believe that this paper helps by bridging the disconnect between
*users* of WL and developers of improved flat-histogram methods. The
majority of developers use pure WL when comparing multiple
flat-histogram methods since there is no direct comparison of WL
followed by a production run with methods such as 1/t-WL or SAMC.

(ii) "It is not clear[ly] understood why the current calculations were
not included in the original article. Although the results presented
by the [authors] are original, they are only a simple extension of the
first work, and do not meet the standards of PRE readers."

The second referee asks an important question here to which there are
two answers. We only became aware of "production run" WL as a result
of the response to our first paper. As mentioned (hinted at) in the
acknowledgement, Johannes Zierenberg emailed us after reading our
paper and suggested that we test WL followed by a production run with
a number of flat-histogram methods.  Due to the lack of literature
detailing its performance vs other flat-histogram methods, the 2D
Ising model seems the logical starting point for such a comparison.

Second, SAD was initially tested on methods where it could benefit
from knowing the temperature range of interest in advance.  It made
sense in the first paper to focus on systems in which the energy range
of interest was *not* easily known as these systems are both more
common and are the kind of systems in which we expect SAD to be
applied.  Nevertheless, it is still valuable to test SAD's convergence
properties on a system where the range of energies is easily known.

We also point out that these simulations are time-consuming.  Some of
those in the current paper ran for more than 6 months.  It is not
reasonable to expect us to complete all possible simulations before
publication of a paper.

In response to the second referee's comments, we have better
highlighted these two important points in the abstract, discussion
paragraph, and conclusion.
