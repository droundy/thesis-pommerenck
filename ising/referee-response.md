Dear editor,

We appreciate the referees' careful reading of our paper, and feel that their
feedback has led to a significant improvement. We have addressed all of the
first referee's comments by added two paragraphs to the introduction and a
section in the paper to 'Address the general problem of all "flat-histogram"
methods'. We also strengthen the abstract to address the concerns raised by the
second referee. We point out that: (i) Our paper is the first of its kind (to
the best of our knowledge) to compare ``production run'' WL with other recent
flat-histogram methods. (ii) It is important to test SAD's convergence
properties on the 2D Ising where the range of energies is easily known (thus
biasing in favor of WL methods). ...

Response to first referee
# -----------------------------------------------------------------------------#

The first referee points out that 'Simulation set-up and analysis protocols are
sound and give confidence in the quality of the presented data and
conclusions.'

The referee also asks that three additions be made to the paper.

(i) "It would be worthwhile to work out the [slower convergence for SAD due to
the energy range being an a priori parameter in WL/MUCA] in the paper, and, if
possible, quantify the loss in efficiency coming from energy range
determination."

We have added a paragraph in the discussion section that explains the loss in
efficiency SAD suffers by not specifying an energy range a priori. We are able
to quantify this by identifying the last time SAD finds a new important energy
and then finding the fraction of moves (in the end) that are outside the range
of interesting energies.

(ii) "Address the general problem of all "flat histogram" methods that the
performance of these methods is due to hidden barriers almost always lower than
what one would expect from a random walk in energy (in most cases still
exponentially growing). This has been discussed by Nadler et al (PRE 75 (2007)
1) who also discus a protocol for optimizing the ensemble which minimizes the
residual exponential slowing down. The authors should point out if and how SAD
can be adapted to such "non-flat-histogram" methods."

We have added a paragraph that 'Address the general problem of all "flat
histogram" methods' in the introduction. We discuss how the performance is
impacted due to hidden barriers. We also add that WL has already been adapted
to replica-exchange and point out that this should be possible for SAD as well.

(iii) Add a few short remarks on comparison with replica exchange methods would
be helpful.

We have added a 'few short remarks' on comparison with replica exchange in the
introduction. WL has been formulated as replica-exchange WLMC whereby it adopts
the concept of conformational swapping directly from parallel tempering. Each replica
is allowed to traverse the entire energy range via overlapping energy windows.

While we have performed limited comparison with replica exchange methods (in
our previous work we found the heat capacity to be consistent with computations
performed with replica exchange), SAD could be extended in much the same way
that WL was. In fact, this would make an interesting comparison with replica
exchange statistical temperature Monte Carlo (as they both specify a
temperature range) and should be studied on a temperature range system in a
future work.

Response to second referee
# -----------------------------------------------------------------------------#

We thank the second referee for bringing to our attention two points so that we
can better address them in our paper and make them completely clear. 

(i) "The work is almost a repetition of what the authors did when they
introduced the algorithm."

Our paper is the first of its kind (to the best of our knowledge) to compare
``production run'' WL with other flat-histogram methods. The first (not the
second) referee notes rather strongly that this is the only way that WL should
be done; however, no research work has ever compared this with other
flat-histogram methods and certainly not with the completeness shown here.
Also, the very fact that it is compared with pure WL is proof for the first
referee's statement that "this is how WL should be done" since it performs
superior to pure WL.

We believe that this paper helps by bridging the disconnect between *users* of
WL and developers of improved flat-histogram methods. The majority of
developers use pure WL when comparing multiple flat-histogram methods since
there is no direct comparison of WL followed by a production run with methods
such as 1/t-WL or SAMC.

(ii) "It is not clear[ly] understood why the current calculations were not
included in the original article. Although the results presented by the
[authors] are original, they are only a simple extension of the first work, and
do not meet the standards of PRE readers."

The second referee asks an important question here to which there are two
answers. We only became aware of "production run" WL as a result of the
response to our first paper. As mentioned (hinted at) in the acknowledgement,
Johannes Zierenberg emailed us after reading our paper and suggested that we
test WL followed by a production run with a number of flat-histogram methods.
Due to the lack of literature detailing its performance vs other flat-histogram
methods, the 2D Ising model seems the logical starting point for such a
comparison.

Second, SAD was initially tested on methods where it could benefit from knowing
the temperature range of interest in advance. It is critical to test it's
convergence properties on a different type of system where the range of
energies is easily known (thus biasing in favor of WL methods). Also in the
case of the 128 X 128 2D Ising system, the computation time is quite lengthy
(more than 6 months for some calculations) representing a superior and robust
comparison result.

In response to the second referee's comments, we have better highlighted these
important points in the abstract and conclusion.