# Model-Predictive-control-and-MDP

The given code implements a dynamic programming algorithm for revenue optimization or booking decisions. Here is a breakdown of the code:

1. Initialization: The code initializes variables, loads data from .mat files, and sets parameters such as lambda, mu, sample length, discount factor, and others.

2. Sample Generation: The code generates a sample space by calling the "state_space_generation" function, which creates a set of possible states based on the given parameters.

3. Revenue Generation: The code calculates revenue values for each state in the sample space using the "r_k_horizon_generation_unbounded_bounded" function. These revenue values are stored in the "r_k_horizon_h" variable.

4. Dynamic Programming: The code enters a loop to perform dynamic programming calculations. It starts by initializing the value function ("J") with zeros. Then, for each state in the sample space, it iteratively updates the value function by considering all possible decisions and their associated costs and rewards. The "Exact_ADP_Booking_unbounded_bounded_12" function is likely responsible for updating the value function based on the current state and revenue values.

5. Result Analysis: The code tracks and stores various metrics such as bounded and unbounded reservations, rejection rates, initial state visits, and others.

6. Simulation: The code enters a simulation loop where it runs multiple experiments (num_exper) to simulate revenue generation and booking decisions. It iterates over the sample space, calculates the expected revenue for each state based on the updated value function, and makes decisions accordingly.

7. State Analysis and Transition: The code analyzes the current state, revenue values, and capacity constraints to make the next decision. It considers various factors such as the number of available slots, lambda and mu values, and implements different decision rules based on the conditions.

8. Termination and Output: The simulation loop continues until the end condition is met, and then the results, such as expected revenue and state transitions, are analyzed and stored.

Overall, the code combines dynamic programming techniques with simulation to optimize revenue or booking decisions. It iteratively updates the value function based on the current state and revenue values and uses it to make optimal decisions in a simulated environment.
