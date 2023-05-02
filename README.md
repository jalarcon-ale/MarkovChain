# MarkovChain
Contact info: alejandrojavier.alarcongonzalez@uantwerpen.be


The Python file Markov_SIR contains the code to generate realizations (simulations) of a discrete-time Markovian SIR model. This model is 
introduced next.

The number of newly infected and removed individuals at time $t+h$ are denoted by $I_{t+h}^{new}$ and $R_{t+h}^{new}$, respectively. 
Their update rule reads as follows


$$
\begin{split}
    I_{t+h}^{new}\sim Binomial(S_t, p_1(t)=&1-\exp(-h\beta I_t)), \\
    %\label{eq:binomial_update}
    R_{t+h}^{new}\sim Binomial(I_t, p_2(t)=&1-\exp(-h\gamma)).
\end{split}
$$


In this way, the number of susceptible, infected and recovered invividuals is given by the next system of recurrence relations 

$$
    \begin{split}
        %\label{eq:recurrence_relations}
        S_{t+h} = {} & S_{t} - I_{t+h}^{new}, \\
        I_{t+h} = {} & I_{t} + I_{t+h}^{new} - R_{t+h}^{new}. \\
        R_{t+h} = {} & R_{t} + R_{t+h}^{new}.
    \end{split}
$$

The script in Markov_SIR takes as an input a range of population sizes varying from $10$ to $100$, and calculates the number of time steps (in hours) 
it takes for the epidemic to die out. We follow the convention of considering that in the beginning there is one infectious individual  and the rest of 
the populations is susceptible to infection. 



%The file functions.py contains the contact matrices used for the calculation of p_star
