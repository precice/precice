# Livegraph

This tool repeatedly plots the last given amount of lines from a preCICE conversion or iteration file.

For more information please see the help:
```
livegraph.py --help
```

## Example

When using convergence measures in preCICE, the solvers will generate convergence logs in their working directories.
These have the following naming scheme `precice-XXX-convergence.log` with `XXX` being the participant name.

Example content of `precice-SolverTwo-convergence.log`:
```
Timestep  Iteration  resNorm(0)  resNorm(1)  
1  1  1.0000000000000000  1.0000000000000000  
1  2  0.9989999579094125  0.9989999555832783 
1  3  1.0279806948721781  1.0052154542721008 
1  4  1.0676872486194695  0.9801198387359404 
1  5  1.0344302257268507  0.9680459014699030 
1  6  0.8635571656036130  0.5242130707806104 
1  7  1.1916321820307565  0.6713041685357296 
1  8  1.0319935686713446  0.6783014521877981 
1  9  1.0160402924851608  0.7538397562490677 
1  10  1.0472963750750830  0.2390668233326704  
1  11  1.0439167724591154  0.1990164586784142  
1  12  1.0719680592012939  0.1327153763711783  
1  13  1.2486945738321200  0.0234804736688682  
1  14  1.2970734078552075  0.0067472965414147  
1  15  1.3046996163929816  0.0087649918345370
```

To live-plot the `resNorm(0)` at the 20 last iterations, use the following command:
```
python livegraph.py -x 1 -y 2 -l 20 precice-SolverTwo-convergence.log
```
