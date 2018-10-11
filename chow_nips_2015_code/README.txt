Matlab code for the NIPS paper titled 'Risk-Sensitive
Decision-Making with Conditional Value-at-Risk'

Requirements:
1. Matlab (tested on R2014b)
2. CPlex for Matlab

The script run_gridworld_small.m performs CVaR value iteration on a
small grid-world domain, and produces similar results to the ones
reported in the paper. It takes several minutes to run.

The script run_gridworld_large.m performs CVaR value iteration on
the large grid-world domain reported in the paper. It takes roughly
1 hour to run.

Usage with general MDPs: 
------------------------
- The function VI_CVaR.m performs CVaR value-iteration on a general MDP, specified by the MDP_class.m
interface. 
- For finite MDPs, in which the P and R matrices are known, the straightforward Finite_MDP_class.m may be used. 
- The grid-world MDP class create a finite MDP from a png image of the grid-world. 
- The Chain_MDP_class provides a simple chain MDP, useful for bench-marking risk-sensitive algorithms.
