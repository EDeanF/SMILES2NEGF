# SMILES2NEGF
compute transmission function within NEGF formalism with SMILES input

In the field of single-molecule electronics, ab-initio calculations of transport characteristics are primarily driven by two techniques: Non-Equilibrium Greens Function (NEGF) methods and quantum master equations.

Density function theory (DFT) is often needed in conjunction with NEGF in order to compute the transport characteristics of a given molecule. However, DFT calculations can be quite slow for large systems and, consequently, does not lend itself to exploration through quick iterations. 

For simple organic molecules, a Huckel or tight-binding model may be sufficient to obtain qualitative understanding of the main features of the transmission functions. For example, different binding sites can be sampled to look for interference features. 

This code is best used in conjunction with DFT calculations. The hope is that by being able to quickly iterate through different models, major features in the DFT calculations can be understood more deeply by non-experts. The quick iteration is aided by the use of SMILES. Commonly used open-source and commerical molecule editors can be used to quickly draw realistic molecular systems or simple graphs which can be exported as SMILES notation.  

Feb 8, 2019: The first version of this code is only applicable to the pi-system. 
