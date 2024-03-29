# MarkovChain
Contact info: alejandrojavier.alarcongonzalez@uantwerpen.be


The Python file Markov_SIR contains the code used to generate the data in the pre-print "Computation of Expected Epidemic Duration". The file endofepi.py contains the class being called in Markov_SIR. The files: converging_means.txt, simulated_durations.txt and exact_means.txt contain all the (synthetic) data used in the manuscript. 

The discrete-time Markovian SIR model used in the manuscript is introduced next. The number of newly infected and removed individuals at time $t+h$ are denoted by $I_{t+h}^{new}$ and $R_{t+h}^{new}$, respectively. 
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
