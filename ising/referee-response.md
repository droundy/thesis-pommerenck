Dear editor,

We appreciate the referees' careful reading of our paper, and feel that their
feedback has led to a significant improvement. We have addressed all of the
first referee's comments by added two paragraphs to the introduction and a
section in the paper to 'Address the general problem of all "flat-histogram"
methods'. We also strengthen the abstract to address the concerns raised by the
second referee. We point out that: (i) Our paper is the first of its kind (to
the best of our knowledge) to compare ``production run'' WL with other
flat-histogram methods. (ii) It is important to test SAD's convergence
properties on the 2D Ising where the range of energies is easily known (thus
biasing in favor of WL methods). ...

Response to first referee
# -----------------------------------------------------------------------------#

The first referee points out that 'Simulation set-up and analysis protocols are
sound and give confidence in the quality of the presented data and
conclusions.'

The referee also asks that three additions be made to the paper.

(i) It would be 'worthwhile to work out the [slower convergence for SAD due to
the energy range being an a priori parameter in WL/MUCA] in the paper, and, if
possible, quantify the loss in efficiency coming from energy range
determination.'

We have added a paragraph? subsection in the SAD/results portion of the paper that
address this.

(ii) Address the general problem of all "flat histogram" methods that the
performance of these methods is due to hidden barriers almost always lower than
what one would expect from a random walk in energy (in most cases still
exponentially growing). This has been discussed by Nadler et al (PRE 75 (2007)
026109) who also discus a protocol for optimizing the ensemble which minimizes
the residual exponential slowing down. The authors should point out if and how
SAD can be adapted to such "non-flat-histogram" methods.

We have added a paragraph that 'Address the general problem of all "flat
histogram" methods' in the introduction.

(iii) Add a few short remarks on comparison with replica exchange methods would
be helpful.

We have added a 'few short remarks' on comparison with replica exchange in the
introduction.

Response to second referee
# -----------------------------------------------------------------------------#

We thank the second referee for bringing to our attention two points so that we
can better address them in our paper and make them completely clear. 

(i) The work is almost a repetition of what the authors did when they
introduced the algorithm.

Our paper is the first of its kind (to the best of our knowledge) to compare
``production run'' WL with other flat-histogram methods. The first (not the second) referee
notes rather strongly that this is the only way that WL should be done;
however, no research work has ever compared this with other flat-histogram
methods and certainly not with the completeness shown here. Also, the very fact that
it is compared with pure WL is proof for the first referee's statement that
"this is how WL should be done" since it performs superior to pure WL.

(ii) It is not clear[ly] understood why the current calculations were not
included in the original article. Although the results presented by the
[authors] are original, they are only a simple extension of the first work, and
do not meet the standards of PRE readers.

The second referee an important question here to which there are two answers.
First, we only recently became aware of ``production run'' WL and wanted to
rapidly compare this with other histogram methods due to the lack of literature
detailing its performance vs other flat-histogram methods. Second, SAD was
initially tested on methods where it could benefit from knowing the temperature
range of interest in advance. It is critical to test it's convergence
properties on a different type of system where the range of energies is easily
known (thus biasing in favor of WL methods).

In response to the second referee's comments, we have better highlighted these
important points in the abstract.