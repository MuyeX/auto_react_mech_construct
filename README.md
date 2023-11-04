# ODE_generator.py

The mechanism generator file simply grabs a string of reactions, for example: 

    r1 = 'A -> B + C'
    
    r2 = 'C -> B + D'
    
    r3 = 'D -> E'
    
    r4 = 'E -> B + F'
    
    reactions = [r1, r2, r3, r4]
    
    mechanism = make_system(reactions)

NOTE: for the above equation we would know that 'A -> 3B + F'.
And then turns them into an ODE system:

    print(mechanism)

    def kinetic_model(x, init, k1, k2, k3, k4):
        
        CA,CB,CC,CD,CE,CF = init
        
        dAdt = - k1*CA
        
        dBdt = k1*CA + k2*CC + k4*CE
        
        dCdt = k1*CA - k2*CC
        
        dDdt = k2*CC - k3*CD
        
        dEdt = k3*CD - k4*CE
        
        dFdt = k4*CE
        
        return dAdt,dBdt,dCdt,dDdt,dEdt,dFdt

This is useful because one can easily get the ODE from a mechanism defined reaction mechanism. From here we would need to figure out: how instead of getting the string from the 'print' statement, simply be able to use it (I am not sure if this is possible or not, and if it is possible, if it is difficult to do); and how to iteratively propose the strings that will make up the reaction system.

# fructose_to_HMF.py

The case study file is exactly that, the case study that we will be working on (hopefully this will yield good results, if not, we can always find a different one). From the 'ground-truth' mechanism, I generate some data with some random initial concentrations of the system, some random rate constants and measure the system for a random time horizon. You will find that I have only plotted species A, B and F. This is because, in a real experiment, it is likely that we would only be able to measure those species and not the other intermediate. In this case study, we can increase the number of experiments, we can change the initial concentrations we use (we can only vary the concentrations of A, B and F), we can change the rate constants, and we can change the time horizon from which we collect data (this should be more or less reasonable).

NOTE: ignore the benzoin_ccondensation case study; this is the more complicated case study, which we will try to tackle after the conference paper.

# Constraints and rules for string construction

1. A single reaction cannot have more than two species on either side of the arrow. For example: A -> B, A + B -> C, A -> B + C, A + B -> C + D are all okay, but A + B + C -> D would not be okay this is because tri-molecular transformations are extremely rare (anything above that tending to be impossible). That is probably the only constraint.
2. Knowing the stoichiometry of the reaction, which in this case is 'A -> 3B + F', we also know that the smallest possible mechanism will be one which has the same number of reactions as the biggest stoichiometric coefficient, in this case 3. 

Rule construction example:
1. Start with the simplest possible mechanism. Get the ODE, minimize the error with respect to the rate constants (to estimate the kinetic parameters) and calculate Akaike information criterion (AIC).
2. From this starting mechanism, we should also consider the same mechanism but with a reversible step.
3. In iteration 2, we have more flexibility since we are adding one intermediate to the mechanism, we could create more complex mechanisms and their reversible counterparts.
4. We would optimize every possible mechanism and compute the AIC value for each of the mechanisms in iteration 2.
5. The best mechanism from iteration 2 would be compared to the best mechanism in iteration 1. If the mechanism in iteration 1 is better we could terminate (maybe it is a better idea to allow at least two iterations before terminating). If the mechanism in iteration 2 is better, then we would continue to iteration 3 where we would add yet another reaction.
6. This process would continue until we find the smallest feasible reaction mechanism.

I hope that is more or less clear. We can discuss in more detail if it is unclear. For the ESCAPE paper, we can reduce the flexibility of the algorithm (e.g, not considering reversible reactions), to ensure that we are able to obtain the underlying mechanism, if it is needed.

# Example for fructose to HMF case study

Given the constraints explained before, the simplest mechanism that is physically possible would be (not considering reversibility, which we can ignore for now):

    r1 = 'A -> B + C'
    
    r2 = 'C -> B + D'
    
    r3 = 'D -> B + F'

This is the only possibility in iteration 1. For iteration 2, we would introduce another intermediate, in this case 'E'. The possible mechanisms with an extra intermediate, keeping in mind above constraints would be the following:

    r1,1 = 'A -> B + C'
    r1,2 = 'C -> B + D'
    r1,3 = 'D -> B + E'
    r1,4 = 'E -> F'
    
    r2,1 = 'A -> B + C'
    r2,2 = 'C -> B + D'
    r2,3 = 'D -> E'
    r2,4 = 'E -> B + F'
        
    r3,1 = 'A -> B + C'
    r3,2 = 'C -> D'
    r3,3 = 'D -> B + E'
    r3,4 = 'E -> B + F'
    
    r4,1 = 'A -> C'
    r4,2 = 'C -> B + D'
    r4,3 = 'D -> B + E'
    r4,4 = 'E -> B + F'

These are the only possibilities (for now) in iteration 2. In the future, we will consider reversibility and potential reactions between intermediates and reactants/products. But let us not complicate the problem further for now.

# Packages that turn strings into ODE systems

The one that I used in this example is the Python one that I had talked to you about. You can find it at: https://github.com/EPiCs-group/ODE_fitter where the function that does that has the same name as the one you find in mechanism_constructor.py

I have made some slight changes to the original function to make it work more in line with this work.

The one in Julia, which I have never tried, but could perhaps be better, you can find it at: https://docs.sciml.ai/Catalyst/stable/
