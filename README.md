# Installation Instructions

To use this package, you'll need to import it using the following method:

1. First, download the `DTMC_0.9.1.zip` file and save it to a location on your computer.
2. Then, run the following R code to install and load the package:

```R
# Replace this with the path where you've saved the DTMC_0.9.1.zip file
path <- '...\'

# Construct the full path to the ZIP file
pkgfile.zip <- paste0(path, 'DTMC_0.9.1.zip')

# Install the package
install.packages(pkgfile.zip, repos = NULL, type = "win.binary")

# Load the package
library(DTMC)
```


# Discrete Time Markov Chain

Consider a process that has a value in each time period. Let
*X*<sub>*n*</sub> denote its value in time period n, and we want to
model for successive values
*X*<sub>0</sub>, *X*<sub>1</sub>, *X*<sub>2</sub>.... We can assume that
the current state only depends on the previous state, which means
*X*<sub>*n*</sub> only depends *X*<sub>*n* − 1</sub>. For this kind of
problem, we define it as discrete time Markov chain. In detail:  

$$
\\begin{aligned}
P(X\_{n+1}=j|X\_n=i,X\_{n-1}=i\_{n-1},...,X\_{0}=i\_{0})=P(X\_{n+1}=j|X\_{n}=i)=P\_{ij}
\\end{aligned}
$$

After having some basic understandings of the DTMC, we can do some
further research on it:  
\* Determine whether the transition probability matrix of a Markov chain
is valid.  
\* Calculate n-step transition probability from one step to another
one.  
\* Calculate the probability that the Markov chain has never been in
several states.  
\* Determine whether two states in a Markov chain access, communicate,
or have some other relationships.  
\* Classify states of a Markov chain into “recurrent” and “transient”.  
\* Determine whether the Markov chain is irreducible.  
\* Graphical representation of the relationship between states, applying
knowledge in graph theory.  
\* Calculate long-run proportions of states of the Markov chain. \*
Calculate expected value of the Markov chain for given functions.  
\* Determine whether the Markov chain is time reversible.  
\* Calculate mean time spent in transient states.  
\* Derive the reversed chain of the original Markov chain.

# Numerical Functions

There are in total 16 functions to implement the discrete time Markov
chain, and they are explained as below:

## Check the Validity of a Matrix: is.DTMC()

By definition, the value in matrix P can be denoted by
*P*<sub>*i**j*</sub>, which is called the one-step transition
probability from state i to state j. Obviously,

$$ P\_{ij}≥0,\\ ∑ \\limits\_{j=0}^{p}P\_{ij} = 1,\\ for\\ \\ all\\ \\ i,j ≥0.$$

As above, the package only solves finite and non-negative square matrix,
with each row sum to be 1.

To guarantee the input matrix is valid, before using package to solve
stochastic problems, the user can use function to have a check.

The return is a bool value TRUE or FALSE, while TRUE represents the
matrix is valid, FALSE represents the matrix is not. Moreover, function
has already been added into all the function in who involved an input
*P* matrix or transition matrix inner.

Here’s a quick demo of function :

    P = matrix(c(0.3,0.3,0.4,0.1,0.2,0.7,0,0.2,0.8),byrow = TRUE, nrow = 3)
    is.DTMC(P)
    #> [1] TRUE

## Create a Transition Probability Matrix: runif\_DTMC()

The allows the user to build a random transition matrix *P* with certain
dimension as needed. The user only need to give a row number as needed
*nrow* in matrix.

By default, it will create a 3 by 3 square matrix under the roles in
function if the parameter *nrow* is missing.

Here’s a quick demo of function :

    (P <- runif_DTMC(nrow = 5))
    #>            [,1]       [,2]       [,3]       [,4]       [,5]
    #> [1,] 0.23026677 0.30334623 0.28348955 0.05203037 0.13086709
    #> [2,] 0.14036495 0.29787996 0.24824063 0.09693708 0.21657738
    #> [3,] 0.21206825 0.09836034 0.12633614 0.32085912 0.24237614
    #> [4,] 0.05845882 0.31728815 0.14059078 0.03429052 0.44937173
    #> [5,] 0.11247048 0.38879733 0.02607077 0.38670547 0.08595595
    is.DTMC(P)
    #> [1] TRUE
    rowSums(P)
    #> [1] 1 1 1 1 1
    colSums(P)
    #> [1] 0.7536293 1.4056720 0.8247279 0.8908226 1.1251483

## Chapman-Kolmogorov Equations: CKequ()

When using the discrete time Markov chain, we are interested in the
n-step transition probabilities. Since Markov chain is time homogeneous,
we can represent the probability as

$$
\\begin{aligned}
P\_{ij}^n=P(X\_n=j|X\_0=i). 
\\end{aligned}
$$

In general, we have the equation

$$
\\begin{aligned}
P\_{ij}^{m+n}=\\sum\\limits\_{k}{P}\_{ik}^mP\_{kj}^n,\\ for\\ all\\ m,n\\geq1.
\\end{aligned}
$$

These equations are called Chapman-Kolmogorov equations. By these
equations, we can prove that

$$
\\begin{aligned}
P\_{ij}^n=P^{n}(i,j),\\ for\\ all\\ n\\geq1.
\\end{aligned}
$$

Users can specify the transition probability matrix *P*, starting state
*i*, ending state *j*, and number of steps to calculate the n-step
transition probabilities.

Here is an example of the transition probability from state 2 to state 1
in 2 steps given the probability matrix.

    P <- matrix(c(0.5, 0.4, 0.1, 
                  0.3, 0.4, 0.3, 
                  0.2, 0.3, 0.5), byrow = TRUE, nrow = 3)
    i = 2
    j = 1
    n = 2
    prob <- CKequ(P, i, j, n)
    prob
    #> [1] 0.33

## Probability in Absorbing States: Abprob()

Consider a Markov chain with transition probabilities matrix P. Let A be
a set of states, and suppose users are interested in the probability
that the Markov chain ever enters (or never enters) any of the states in
A by time m. That is,for a given state *i* ∈ *A*, we are interested in
determining
*P*(*X*<sub>*k*</sub>∈𝒜 for some *k*=1,…,*m*∣*X*<sub>0</sub>=*i*) (or
*P*(*X*<sub>*k*</sub>∉𝒜 for any *k*=1,…,*m*∣*X*<sub>0</sub>=*i*))

To determine the preceding probability we will define a Markov chain
{*W*<sub>*n*</sub>:*n*≥0} whose states are the states that are not in A
plus an additional state, which we will call A. Once the
*W*<sub>*n*</sub> enters state A it remains there forever. Any states in
A is absorbing state.

-   **Example**: Consider a Markov chain with states 1, 2, 3, 4, 5, and
    suppose that we want to compute
    *P*(*X*<sub>4</sub> = 2, *X*<sub>3</sub> ≤ 2, *X*<sub>2</sub> ≤ 2, *X*<sub>1</sub> ≤ 2|*X*<sub>0</sub> = 1)
    All we know about the transition probabilities are:
    *P*<sub>11</sub> = 0.3, *P*<sub>12</sub> = 0.3, *P*<sub>21</sub> = 0.1, *P*<sub>22</sub> = 0.2
    We can use to solve this problem.

<!-- -->

    P = matrix(c(0.3, 0.3, 0.4, 0, 0,
                 0.1, 0.2, 0.7, 0, 0,
                 rep(0, 15)),
               nrow = 5, byrow = TRUE)
    # if you don't know absorb state information, you can put any probability such that row sum is 1
    result <- Abprob(P, c(3, 4, 5), 1, 2, 4)
    #> [1] "The Matrix is not valid."
    result
    #> [1] NA

## Check Whether a State is Accessible: is.accessible()

State j is said to be accessible from state i if
*P*<sub>*i**j*</sub><sup>*n*</sup> &gt; 0 for some *n* ≥ 0. Note that
this implies that state j is accessible from state i if and only if,
starting in i, it is possible that the process will ever enter state j.
*i* → *j*

The function give a way to check whether state j is accessible from
state i. This is an example:

    P <- matrix(c(1/2, 0, 1/8, 1/4, 1/8, 0,   0,
                  0,   0, 1,   0,   0,   0,   0,
                  0,   0, 0,   1,   0,   0,   0,
                  1,   0, 0,   0,   0,   0,   0,
                  0,   0, 0,   0,   1/2, 0,   1/2,
                  0,   0, 0,   0,   1/2, 1/2, 0,
                  0,   0, 0,   0,   0,   1/2, 1/2),
                nrow = 7,
                byrow = TRUE)
    is.accessible(P, 1, 2)
    #> [1] FALSE
    is.accessible(P, 1, 3)
    #> [1] TRUE

## Check Whether Two States Communicate: is.communicate()

Two states i and j that are accessible to each other are said to
communicate, and we write *i* ↔ *j* The relation of communication
satisfies the following three properties:

1.  Each state i communicates with itself.
2.  If state i communicates with state j, then state j communicates with
    state i.
3.  If state i communicates with state j, and state j communicates with
    state k, then state i communicates with state k.

The properties of communication are very useful for analyzing DTMC.

There is an example to check communication of DTMC:

    P <- matrix(c(1/2, 0, 1/8, 1/4, 1/8, 0,   0,
                  0,   0, 1,   0,   0,   0,   0,
                  0,   0, 0,   1,   0,   0,   0,
                  1,   0, 0,   0,   0,   0,   0,
                  0,   0, 0,   0,   1/2, 0,   1/2,
                  0,   0, 0,   0,   1/2, 1/2, 0,
                  0,   0, 0,   0,   0,   1/2, 1/2),
                nrow = 7,
                byrow = TRUE)
    is.accessible(P, 1, 2)
    #> [1] FALSE
    is.accessible(P, 2, 1)
    #> [1] TRUE
    is.communicate(P, 1, 2)
    #> [1] FALSE
    is.communicate(P, 2, 1)
    #> [1] FALSE
    is.accessible(P, 1, 3)
    #> [1] TRUE
    is.accessible(P, 3, 1)
    #> [1] TRUE
    is.communicate(P, 1, 3)
    #> [1] TRUE
    is.communicate(P, 3, 1)
    #> [1] TRUE

## Check Whether States Belong to the Same Class: is.class()

We divide the state space into disjoint classes by putting states that
communicate with each other into the same class. Those states in a class
will have some consistent properties.

And, a Markov chain is called irreducible if its state space has only
one class and all the states communicate.

This function can help users check some states of a (block) given
transition matrix whether are i the same class. If the input P is a
whole transition matrix of a DTMC and all state (or state is not
entered), the function will check whether it is irreducible

There is an example:

    P <- matrix(c(1/2, 0, 1/8, 1/4, 1/8, 0,   0,
                  0,   0, 1,   0,   0,   0,   0,
                  0,   0, 0,   1,   0,   0,   0,
                  1,   0, 0,   0,   0,   0,   0,
                  0,   0, 0,   0,   1/2, 0,   1/2,
                  0,   0, 0,   0,   1/2, 1/2, 0,
                  0,   0, 0,   0,   0,   1/2, 1/2),
                nrow = 7,
                byrow = TRUE)
    if (is.class(P,c(1,2,4))) {
      print('state 1,2,4 is in the same class!')
    }
    if (is.class(P,c(1,3,4))) {
      print('state 1,3,4 is in the same class!')
    }
    #> [1] "state 1,3,4 is in the same class!"
    if (is.class(P)) {
      print("The DTMC is irreducible")
    } else {
      print("The DTMC is not irreducible")
    }
    #> [1] "The DTMC is not irreducible"

## Check Whether the Markov Chain is Transient or Recurrent: classify()

For any state i we let *f*<sub>*i*</sub> denote the probability that,
starting in state i, the process will ever reenter state i in finite
steps. State i is said to be recurrent if fi = 1 and transient if fi
&lt; 1.

If state i is recurrent and state j communicates state i, state j is
also recurrent. So all state in the same class are recurrent or are
transient.

The function can check whether a block matrix of the whole transition
matrix with some state is recurrent or transient.

There is an example:

    P <- matrix(c(1/2, 0, 1/8, 1/4, 1/8, 0,   0,
                  0,   0, 1,   0,   0,   0,   0,
                  0,   0, 0,   1,   0,   0,   0,
                  1,   0, 0,   0,   0,   0,   0,
                  0,   0, 0,   0,   1/2, 0,   1/2,
                  0,   0, 0,   0,   1/2, 1/2, 0,
                  0,   0, 0,   0,   0,   1/2, 1/2),
                nrow = 7,
                byrow = TRUE)
    classify(P) # the state in P is not in the same class
    #> [1] "The state of the transition matrix are not in the same class"

    classify(P[c(1, 3, 4), c(1, 3, 4)])
    #> [1] "transient"
    cat("the class of state 1,3,4 is", classify(P[c(1, 3, 4), c(1, 3, 4)]))
    #> the class of state 1,3,4 is transient

## Information of a Class: class\_info()

This function categorizes all states and summarizes all classes, as long
as you enter the transition matrix.

This is an example:

    P <- matrix(c(1/2, 0, 1/8, 1/4, 1/8, 0,   0,
                  0,   0, 1,   0,   0,   0,   0,
                  0,   0, 0,   1,   0,   0,   0,
                  1,   0, 0,   0,   0,   0,   0,
                  0,   0, 0,   0,   1/2, 0,   1/2,
                  0,   0, 0,   0,   1/2, 1/2, 0,
                  0,   0, 0,   0,   0,   1/2, 1/2),
                nrow = 7,
                byrow = TRUE)
    info <- class_info(P)
    info # info$state is the column number of states contained by different classes,
    #> $state
    #> $state[[1]]
    #> [1] 1 3 4
    #> 
    #> $state[[2]]
    #> [1] 2
    #> 
    #> $state[[3]]
    #> [1] 5 6 7
    #> 
    #> 
    #> $prob_martix
    #> $prob_martix[[1]]
    #>      [,1]  [,2] [,3]
    #> [1,]  0.5 0.125 0.25
    #> [2,]  0.0 0.000 1.00
    #> [3,]  1.0 0.000 0.00
    #> 
    #> $prob_martix[[2]]
    #> [1] 0
    #> 
    #> $prob_martix[[3]]
    #>      [,1] [,2] [,3]
    #> [1,]  0.5  0.0  0.5
    #> [2,]  0.5  0.5  0.0
    #> [3,]  0.0  0.5  0.5
    #> 
    #> 
    #> $class
    #> [1] "transient" "transient" "recurrent"
    #> 
    #> $group
    #> [1] 1 2 1 1 3 3 3
         # info$prob_matrix is each block matrix of different classes,
         # info$class is the classfications of each class (recurrent or transient)
         # info$group is the Group of ordinal number of each class

## Derive the Graph: get\_graph()

It is often helpful to lay out the model in the so-called transition
diagram (or transition probability graph), whose vertices are the states
and whose edges are the possible transitions. By recording the numerical
values of *p*<sub>*i**j*</sub> near the corresponding edges, one can
visualize the entire model in a way that can make some of its major
properties readily apparent. We can use some methods of graph theory to
analyze the problem of discrete Markov chains

This function will return a graph object of , so you can use some
functions of to analyze the Discrete time Markov chains.

Now, consider a graph *G* = (*V*, *E*) of a DTMC From the point of view
of graph theory:

-   state *s* ∈ *S*: vertex *v* ∈ *V*

-   transition *P*<sub>*i**j*</sub> &gt; 0: edge (*i*, *j*) ∈ *E*

-   accessible *i* → *j*: there is a path
    *i*, *i*<sub>1</sub>, ..., *i*<sub>*n* − 1</sub>, *j* where
    (*i*, *i*1), (*i*1, *i*2), ..., (*i**n* − 1, *j*) ∈ *E**d**g**e* of
    the graph, and the least possible time to do so is the length of the
    shortest such path.

-   communicate *i* ↔ *j*: i and j are connected.

-   irreducible *i* ↔ *j*, ∀*i*, *j* ∈ *S*: the digraph G is strongly
    connected

This is an example of using distance() in the igraph package to
calculate the shortest transition path length.

    P <- matrix(c(1/2, 0, 1/8, 1/4, 1/8, 0,   0,
                  0,   0, 1,   0,   0,   0,   0,
                  0,   0, 0,   1,   0,   0,   0,
                  1,   0, 0,   0,   0,   0,   0,
                  0,   0, 0,   0,   1/2, 0,   1/2,
                  0,   0, 0,   0,   1/2, 1/2, 0,
                  0,   0, 0,   0,   0,   1/2, 1/2),
                nrow = 7,
                byrow = TRUE)
    g <- get_graph(P)

    library(igraph)
    #> Warning: 程辑包'igraph'是用R版本4.1.2 来建造的
    #> 
    #> 载入程辑包：'igraph'
    #> The following objects are masked from 'package:stats':
    #> 
    #>     decompose, spectrum
    #> The following object is masked from 'package:base':
    #> 
    #>     union
    distances(g, weights = NA) #This is done by calling the distances() function of the igraph package. The (i, j) element of this generated matrix represents a minimum of steps from state i to state j,
    #>      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
    #> [1,]    0    2    1    1    1    2    2
    #> [2,]    2    0    1    2    3    4    4
    #> [3,]    1    1    0    1    2    3    3
    #> [4,]    1    2    1    0    2    3    3
    #> [5,]    1    3    2    2    0    1    1
    #> [6,]    2    4    3    3    1    0    1
    #> [7,]    2    4    3    3    1    1    0

## Graphical Visualization: plot\_graph()

By recording the numerical values of *p*<sub>*i**j*</sub> near the
corresponding edges, one can visualize the entire model in a way that
can make some of its major properties readily apparent.

    P <- matrix(c(1/2, 0, 1/8, 1/4, 1/8, 0,   0,
                  0,   0, 1,   0,   0,   0,   0,
                  0,   0, 0,   1,   0,   0,   0,
                  1,   0, 0,   0,   0,   0,   0,
                  0,   0, 0,   0,   1/2, 0,   1/2,
                  0,   0, 0,   0,   1/2, 1/2, 0,
                  0,   0, 0,   0,   0,   1/2, 1/2),
                nrow = 7,
                byrow = TRUE)
    plot_graph(P)

![unnamed-chunk-15-1](https://github.com/Junxiaolei/DTMC/assets/143832373/0e32aa69-a6f6-46ce-bcfe-7d32aff5d4b1)


    #> $state
    #> $state[[1]]
    #> [1] 1 3 4
    #> 
    #> $state[[2]]
    #> [1] 2
    #> 
    #> $state[[3]]
    #> [1] 5 6 7
    #> 
    #> 
    #> $prob_martix
    #> $prob_martix[[1]]
    #>      [,1]  [,2] [,3]
    #> [1,]  0.5 0.125 0.25
    #> [2,]  0.0 0.000 1.00
    #> [3,]  1.0 0.000 0.00
    #> 
    #> $prob_martix[[2]]
    #> [1] 0
    #> 
    #> $prob_martix[[3]]
    #>      [,1] [,2] [,3]
    #> [1,]  0.5  0.0  0.5
    #> [2,]  0.5  0.5  0.0
    #> [3,]  0.0  0.5  0.5
    #> 
    #> 
    #> $class
    #> [1] "transient" "transient" "recurrent"
    #> 
    #> $group
    #> [1] 1 2 1 1 3 3 3

We can see that this is a very complex Markov chain, but from the image,
we can clearly see the states of different classes, which are
distinguished by different colored circles, and we can also clearly see
the transition relationship.

## Long-run Proportion: LongRunP()

Let *π*<sub>*j*</sub> denote the long-run proportion of time that the
Markov chain is in state j. Consider an irreducible Markov chain. If the
chain is positive recurrent, the long-run proportions are the unique
solution of the equations

$$
\\begin{aligned}
\\pi\_{j}&=\\sum\_{i}\\pi\_{i}P\_{ij},\\ j\\geq1\\\\
\\sum\_{j}\\pi\_j&=1
\\end{aligned}
$$

In the function , we first judge whether the Markov chain is irreducible
and positive recurrent. If so, we can continue the calculation;
otherwise, the Markov chain does not have long-run proportions. When
calculating, there are two cases. If we have known that the Markov chain
is time reversible (a property of Markov chain), then the calculation
will be simplified; otherwise, there will be a more formal way of
calculation.

Users need to specify the Markov chain by giving its transition
probability matrix *P*. They can also choose to specify whether the
Markov chain is known to be time reversible, which is determined by
*sign*. The default value for that parameter is FALSE if they do not
assign the value. As a result, if the Markov chain has long-run
proportions, a vector indicating long-run proportions for state
(1,2,…,n) will be returned, otherwise NA will be returned.

Here is an example to calculate long-run proportions of 3 states given
the transition probability matrix.

    P <- matrix(c(0.3, 0.3, 0.4, 0.1, 0.2, 0.7, 0, 0.2, 0.8), byrow = TRUE, nrow = 3)
    sign <- FALSE
    LRp <- LongRunP(P, sign)
    LRp
    #> [1] 0.02898551 0.20289855 0.76811594

## DTMC’s Expectation: expect\_DTMC()

Let {*X**n*, *n* ≥ 1} be an irreducible Markov chain with stationary
probabilities *π*<sub>*j*</sub>, *j* ≥ 0, and let *r* be a bounded
function on the state space. Then, with probability 1,

$$\\lim\\limits\_{N→∞}\\frac{∑\_{n=1}^{N}r(X\_{n})}{N}=∑\\limits\_{j}r(j)π\_j.$$

By definition, to compute the expectation using formula
$∑\\limits\_{j}r(j)π\_j$, the state function *r* need to be bounded and
the value *r*(*j*) need to be finite.

Function can be used to compute the expectation. There are three
parameters needed: the valid transition matrix *P*; a vector
*state\_value*, which shows the value of different states;a function
object *FUN* as the formula with respect to *state\_value*, which can be
flexibly designed by users to transform *state\_value* to another one if
needed. For example, *FUN* can be *state\_value+1*, *state\_value^2*,
*sin(state\_value)* and so on.

If the final expectation is finite, the value will be returned. While
infinite, the expectation will be shown as “does not exist” rather than
returned as or∞.

Here shows an example and a quick demo of function :

Suppose there are three kinds of weather here: sunny,cloudy and rainy.
If it is sunny today, the probability to be sunny,cloudy and rainy
tomorrow is 0.3,0.3 and 0.4 respectively. If it is cloudy today, the
probability to be sunny,cloudy and rainy tomorrow is 0.1,0.2 and 0.7
respectively. If it is rainy today, the probability of sunny,cloudy and
rainy tomorrow is 0,0.2 and 0.8 respectively. Now there is an ice cream
seller, he can sell 2, 4 and 6 cupboards of ice cream in sunny,cloudy
and rainy respectively. Sold one cupboard of ice cream, he could make
$20 profit. What is his expected daily revenue of in the long run?

    P = matrix(c(0.3, 0.3, 0.4, 
                 0.1, 0.2, 0.7, 
                 0, 0.2, 0.8), byrow = TRUE, nrow = 3)
    state_value = c(2,4,6)
    FUN = function (x) {x*200}
    result <- expect_DTMC(P, state_value, FUN)
    cat("The seller's expected daily revenue in the long run should be ", result, " dollars.")
    #> The seller's expected daily revenue in the long run should be  1095.652  dollars.

## Check Whether the Markov Chain is Time Reversible: isTReversible()

If for a Markov chain of transition probabilities *P*<sub>*i**j*</sub>,
there exists *π*<sub>*i*</sub> ≥ 0 such that

$$
\\begin{aligned}
\\pi\_iP\_{ij}=\\pi\_jP\_{ji},\\ \\sum\_{i}\\pi\_i=1
\\end{aligned}
$$

Then, the chain is time reversible and *π*<sub>*i*</sub>’s are its
stationary probabilities. Users should give the transition probability
matrix of the Markov chain *P*. The program will first judge whether
long-run proportions exist by another function . If exists, the program
will determine whether the Markov chain is time reversible and return
TRUE/FALSE; if not, the function will return FALSE.

Here is an example to judge whether the Markov chain is time reversible
given its transition probability matrix.

    P <- matrix(c(0.3, 0.3, 0.4, 
                  0.1, 0.2, 0.7, 
                  0, 0.2, 0.8), byrow = TRUE, nrow = 3)
    judge <- isTReversible(P)
    judge
    #> [1] FALSE

## Mean Time Spent in Transient States: TranTime()

Consider a finite-state Markov chain with transient states T = {1,2,…,t}
and *P*<sub>*T*</sub> is the corresponding transient transition
probability matrix. Suppose *s*<sub>*i**j*</sub> is the expected number
of times that the Markov chain is in state j given that it starts in
state i. Let *S* = (*s*<sub>*i**j*</sub>)<sub>1 ≤ *i*, *j* ≤ *t*</sub>,
then

$$
\\begin{aligned}
S = (I-P\_T)^{-1}
\\end{aligned}
$$

Users need to input a transient transition probability matrix *P*, the
starting transient state *i*, and the ending transient state *j*. The
function will first check whether the Markov chain has transient states
by another function . If so, the preceding calculation will be performed
and then the mean time from state i to state j will be returned;
otherwise, NA will be returned.

Here is an concrete example: Consider the gambler’s ruin problem with
*p* = 0.4 and *N* = 7. Starting with 3 units, determine (a) the expected
amount of time the gambler has 5 units, (b) the expected amount of time
the gambler has 2 units. To solve the problem, users need to specify the
transient matrix and two states.

    Pt <- matrix(c(0, 0.4, 0, 0, 0, 0, 
                   0.6, 0, 0.4, 0, 0, 0, 
                   0, 0.6, 0, 0.4, 0, 0, 
                   0, 0, 0.6, 0, 0.4, 0, 
                   0, 0, 0, 0.6, 0, 0.4, 
                   0, 0, 0, 0, 0.6, 0), nrow = 6, byrow = TRUE)
    i1 <- 3
    j1 <- 5
    tranT1 <- TranTime(Pt, i1, j1)
    #> [1] "The Matrix is not valid."
    # The expected output is 0.9228
    tranT1
    #> [1] NA

    i2 <- 3
    j2 <- 2
    tranT2 <- TranTime(Pt, i2, j2)
    #> [1] "The Matrix is not valid."
    # The expected output is 2.3677
    tranT2
    #> [1] NA

## Reversed Markov Chain: RevChain()

Consider an irreducible Markov chain with transition probabilities
*P*<sub>*i**j*</sub>. If we can find positive numbers
*π*<sub>*i*</sub>, *i* ≥ 0, summing to one, and a transition probability
matrix *Q* = \[*Q*<sub>*i**j*</sub>\] such that

$$
\\begin{aligned}
\\pi\_iP\_{ij}=\\pi\_jQ\_{ji},
\\end{aligned}
$$

then the *Q*<sub>*i**j*</sub> are the transition probabilities of the
reversed chain and the *π*<sub>*i*</sub> are the stationary
probabilities both for the original and the reversed chain. Users need
to input the transition probability matrix *P*. The function first
compute long-run proportions of the Markov chain. If long-run
proportions exist, transition matrix of the reversed chain will be
returned; otherwise NA will be returned. Here is an example to compute
the transition probability matrix of the reversed chain.

    P <- matrix(c(0.8, 0.2, 0, 
                  0.1, 0.8, 0.1, 
                  0, 0.2, 0.8), nrow = 3, byrow = TRUE)
    revP <- RevChain(P)

# Advantages and Implements of Learned Points in DTMC:

1.Discrete Time Markov Chain (DTMC) is an important theorem in stochastic process, which requires some matrix operations when solving statistical and probability problems. R language has a lot of built-in statistic functions for the Users in the field of statistics, which means using R language to construct statistical package, it is pretty convenient and concise. On the other hand, R language performs well at graphing. This allows us to build functions for building explicit graphs to represent DTMC problems since a good legend can greatly simplify a complex DTMC problem.
1. When making R package by ourselves, we become more familiar with the structure of a R package. When building the graph functions to visualize the Discrete Time Markov Chain (DTMC) problems, we found that there was no good stochastic process package to use directly. Therefore, we decided to build a graph function for DTMC by analyzing and invoking the relevant package about graph theory. After studying and analyzing the source codes of package \code{igraph}, we built two visualization functions named \code{get_graph()} and \code{plot_graph()}, which performs well to help describing a DTMC problem.
1. When building our package \code{DTMC}, we try to use vectorized operations rather than equivalent 
loops as possible. For some problem that is difficult to be vectorized, we only keep the necessary loops and minimize the amount of the codes in loop. This helps the users run \code{DTMC} functions efficiently.
1. DTMC is an important model in stochastic process theorem, which requires a lot of matrix operations in the practical application of statistical and probability problems. 
1. In roxygen, we have included some examples and demos about how to use functions as required. Moreover, in vignette, we use a lot of formulas of random process theory as annotations in vignette to explain the mathematical and statistic logic of our functions.
1. In order to let users have a better experience, we refer to some commonly used functions in R languages when naming functions. For example, the function \code{is.DTMC()}, \code{is.accessible()} has a similar result with \code{is.numeric()}, returned bool value \code{TRUE} or \code{FALSE}. We want users to be familiar with the functions' names and could quickly handle how to use it from past experience when working using DTMC.
1. We also adjusted the code for practical use. In some functions, we added some user-definable content and gave the default values for the general case. For example, in the function \code{expec_DTMC()}, we allow a self-define function as input to convert the value as needed; function {\code{plot_graph} have _..._ as input to adjust the graph if suggested. 
1. In learning R language, we found that logic and specification are important when writing code. When we write DTMC package, we refer to the establishment process of stochastic process theory. In the production of functions, we also strictly follow a certain order. Users can refer to the order in the functions and find out the logical relationship according to the establishment of DTMC theory. At the same time, in order to facilitate the review and modification of the code in the future, we strictly followed the code standards learned in course when writing.
