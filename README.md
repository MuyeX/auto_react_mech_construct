# mechanism_constructor.py

The mechanism generator file simply grabs a string of reactions, for example: 

r1 = 'A -> B + C'

r2 = 'C -> D'

r3 = 'D -> E'

r4 = 'E -> B + F'

r5 = 'F -> G'

r6 = 'G -> B + H'

r7 = 'H -> I'

reactions = [r1, r2, r3, r4, r5, r6, r7]

mechanism = make_system(reactions)

And then turns them into an ODE system:

print(mechanism)

def kinetic_model(x, init, k1, k2, k3, k4, k5, k6, k7):

    CA,CB,CC,CD,CE,CF,CG,CH,CI = init
    
    dAdt = - k1*CA
    
    dBdt = k1*CA + k4*CE + k6*CG
    
    dCdt = k1*CA - k2*CC
    
    dDdt = k2*CC - k3*CD
    
    dEdt = k3*CD - k4*CE
    
    dFdt = k4*CE - k5*CF
    
    dGdt = k5*CF - k6*CG
    
    dHdt = k6*CG - k7*CH
    
    dIdt = k7*CH
    
    return dAdt,dBdt,dCdt,dDdt,dEdt,dFdt,dGdt,dHdt,dIdt

This is useful because one can easily get the ODE from a mechanism defined reaction mechanism. From here we would need to figure out: how instead of getting the string from the 'print' statement, simply be able to use it (I am not sure if this is possible or not, and if it is possible, if it is difficult to do); and how to iteratively propose the strings that will make up the reaction system.

# case_study.py

The case study file is exactly that, the case study that we will be working on (hopefully this will yield good results, if not, we can always find a different one). From the 'ground-truth' mechanism, I generate some data with some random initial concentrations of the system, some random rate constants and measure the system for a random horizon. You will find that I have only plotted species A, B and I. This is because, in a real experiment, it is likely that we would only be able to measure those species and not the other intermediate. In this case study, we can increase the number of experiments, we can change the initial concentrations we use (we can only vary the concentrations of A, B and I), we can change the rate constants, and we can change the time horizon from which we collect data (this should be more or less reasonable).

# Constraints and rules for string construction

1. A single reaction cannot have more than two species on either side of the arrow. For example: A -> B, A + B -> C, A -> B + C, A + B -> C + D are all okay, but A + B + C -> D would not be okay this is because tri-molecular transformations are extremely rare (anything above that tending to be impossible). That is probably the only constraint.

Rule construction example:
1. Start with the simplest possible mechanism, for the above case, that would be 'A -> B + I'. Get the ODE, minimize the error with respect to the rate constants and calculate Akaike information criterion (AIC).
2. From this starting mechanism, we should also consider the same mechanism but with a reversible step: 'A <-> B + I'. This would have to be written as two separate reactions: 'A -> B + I' and 'B + I -> A'. Again find rate constants and compute AIC.
3. From these two options, we can find the best one based on AIC values. This will be our best mechanism from iteration 1.
4. In iteration 2, we have more flexibility since we are adding one reaction to the mechanism, we could create the following more complex mechanisms: 'A -> C -> B + I' or 'A -> C + B -> I' and their reversible counterparts.
5. We would optimize every possible mechanism and compute the AIC value for each of the mechanisms in iteration 2.
6. The best mechanism from iteration 2 would be compared to the best mechanism in iteration 1. If the mechanism in iteration 1 is better we could terminate (maybe it is a better idea to allow at least two iterations before terminating). If the mechanism in iteration 2 is better, then we would continue to iteration 3 where we would add yet another reaction.
7. This process would continue until we find the smallest feasible reaction mechanism.

I hope that is more or less clear. We can discuss in more detail if it is unclear. For the ESCAPE paper, we can reduce the flexibility of the algorithm (e.g, not considering reversible reactions), to ensure that we are able to obtain the underlying mechanism, if it is needed.

# Packages that turn strings into ODE systems

The one that I used in this example is the Python one that I had talked to you about. You can find it at: https://github.com/EPiCs-group/ODE_fitter where the function that does that has the same name as the one you find in mechanism_constructor.py

I have made some slight changes to the original function to make it work more in line with this work.

The one in Julia, which I have never tried, but could perhaps be better, you can find it at: https://docs.sciml.ai/Catalyst/stable/
