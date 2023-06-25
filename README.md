# Model-Predictive-control-and-MDP

This piece of code is related to my PhD thesis when I proposed theorerical approaches for the reservation in resource allocation problem for MDPs. The code is quite heavy to digest and 3 journal articles have  been published based on this procedure which I refer them at the end of this readme file. 


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


# Revenue Optimization with Dynamic Programming

This code implements a revenue optimization algorithm using dynamic programming. It aims to find the optimal booking decisions that maximize revenue in a simulated environment. The code is written in [MATLAB](https://www.mathworks.com/products/matlab.html) and relies on .mat files for data storage.

## Overview

The code consists of several MATLAB scripts and functions that work together to perform revenue optimization. Here's a brief overview of the main components:

1. **Initialization**: The code initializes variables, loads data from .mat files, and sets parameters such as lambda, mu, sample length, discount factor, and others.

2. **Sample Generation**: The code generates a sample space by creating a set of possible states based on the given parameters. The states represent different scenarios or conditions for revenue generation.

3. **Revenue Generation**: The code calculates revenue values for each state in the sample space. It considers various factors like capacity constraints, available slots, and lambda and mu values to estimate the potential revenue.

4. **Dynamic Programming**: The code performs dynamic programming calculations to determine the optimal value function for revenue optimization. It iteratively updates the value function based on the current state and revenue values, considering all possible decisions and their associated costs and rewards.

5. **Simulation**: The code runs simulations to generate revenue and make booking decisions based on the updated value function. It iterates over the sample space, calculates the expected revenue for each state, and determines the optimal decision using decision rules and constraints.

6. **Result Analysis**: The code tracks and stores various metrics such as bounded and unbounded reservations, rejection rates, initial state visits, and others. It provides insights into the performance of the revenue optimization algorithm.

## Usage

To use this code, follow these steps:

1. Ensure you have MATLAB installed on your system.

2. Clone or download the code repository to your local machine.

3. Open MATLAB and navigate to the downloaded code directory.

4. Customize the parameters in the initialization section of the main script according to your specific problem or scenario.

5. Run the main script (e.g., `main_script.m`) to execute the revenue optimization algorithm.

6. Monitor the console output for progress updates and analysis results.

7. After the simulation completes, analyze the generated results to understand the performance and effectiveness of the revenue optimization algorithm.

## Dependencies

The code has the following dependencies:

- MATLAB (version 2016 or higher)

## License

[MIT License](LICENSE)

Feel free to modify and adapt the code according to your needs.

## Contributing

Contributions to this code are welcome. If you find any bugs, issues, or have suggestions for improvements, please submit a pull request or open an issue.

## Contact

For any questions or inquiries, please contact [aliforootani@ieee.org].

## Articles

1. Forootani, Ali, Raffaele Iervolino, and Massimo Tipaldi. "Applying unweighted least‐squares based techniques to stochastic dynamic programming: theory and application." IET Control Theory & Applications 13.15 (2019): 2387-2398.

2. Forootani, Ali, et al. "Allocating resources via price management systems: a dynamic programming-based approach." International Journal of Control 94.8 (2021): 2123-2143.

3. Forootani, Ali, et al. "Approximate dynamic programming for stochastic resource allocation problems." IEEE/CAA Journal of Automatica Sinica 7.4 (2020): 975-990.





