# Optimal Stopping of Adaptive Clinical Trials
Find the paper at http://snasrol.people.clemson.edu

Diffusion Approx: Implements the diffusion approximation solution to the Bellman's equation descibed in Section 5.3 of the paper. The system of partial equation is solved by adapting a trinomial tree approach.

Gridding: Implements the simulation-based gridding approach which approximates the solution to the Bellman's equation by descretizing the state space and enumerating a large number of future scenarios.

One-step-look-ahead: Implements the one-step look-ahead approximation to the Bellman's equation. The one-step future is projected and assumed to be the last step each time. Then, the optimal decision is derived based on the now given final values.
