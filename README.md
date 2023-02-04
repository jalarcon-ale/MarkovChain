# MarkovChain

The Python file ChainBinomial contains the code to generate realizations (simulations) of a discrete-time Markovian SIR model. This model is 
introduced next.

The number of newly infected and removed individuals at time $t+h$ are denoted by $I_{t+h}^{new}$ and $R_{t+h}^{new}$, respectively. 
Their update rule reads as follows


\begin{equation*}
\begin{split}
    I_{t+h}^{new}\sim Binomial(S_t, p_1(t)=&1-\exp(-h\beta I_t)), \\
    %\label{eq:binomial_update}
    R_{t+h}^{new}\sim Binomial(I_t, p_2(t)=&1-\exp(-h\gamma)).
\end{split}
\end{equation*}


In this way, the number of susceptible, infected and recovered invividuals is given by the next system of recurrence relations 

\begin{equation*}
    \begin{split}
        \label{eq:recurrence_relations}
        S_{t+h} = {} & S_{t} - I_{t+h}^{new}, \\
        I_{t+h} = {} & I_{t} + I_{t+h}^{new} - R_{t+h}^{new}. \\
        R_{t+h} = {} & R_{t} + R_{t+h}^{new}.
    \end{split}
\end{equation*}

The 



The file functions.py contains the contact matrices used for the calculation of p_star
