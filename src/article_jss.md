---
author:
- |
  Daphné Giorgi\
  Sorbonne Université Sarah Kaakaï\
  Le Mans Université Vincent Lemaire\
  Sorbonne Université
bibliography:
- refs.bib
title: "Efficient simulation of individual-based population models: the
  R Package IBMPopSim"
---

# Introduction {#introduction .unnumbered}

In various fields, advances in probability have contributed to the
development of a new mathematical framework for so-called
individual-based stochastic population dynamics, also called stochastic
Individual-Based Models (IBMs).

Stochastic IBMs allow the modeling in continuous time of populations
dynamics structured by age and/or characteristics. In the field of
mathematical biology and ecology, a large community has used this
formalism for the study of the evolution of structured populations (see
e.g.
[@FerTra09; @collet2013rigorous; @BanMel15; @costa2016stochastic; @billiard2016effect; @lavallee2019stochastic; @meleard2019birth; @calvez2020horizontal]),
after the pioneer
works [@FouMel04; @champagnat2006unifying; @tran_2008].

IBMs are also useful in demography and actuarial sciences, for the
modeling of human populations dynamics (see e.g.
[@Ben10; @Bou16; @karoui2021simulating]). They allow the modeling of
heterogeneous and complex population dynamics, which can be used to
compute demographic indicators or simulate the evolution of insurance
portfolios in order to study the basis risk, compute cash flows for
annuity products or pension schemes, or for a fine assessment of
mortality models ([@barrieu2012understanding]). There are other domains
in which stochastic IBMs can be used, for example in epidemiology with
stochastic compartmental models, neurosciences, cyber risk, or
Agent-Based Models (ABMs) in economy and social sciences, which can be
seen as IBMs.\
Many mathematical results have been obtained in the literature cited
above, for quantifying the limit behaviors of IBMs in long time or in
large population. In particular, pathwise representations of IBMs have
been introduced in [@FouMel04] (and extended to age-structured
populations in [@tran_2008]), as measure-valued pure jumps Markov
processes, solutions of SDEs driven by Poisson measures. These pathwise
representations are based on the *thinning* and projection of Poisson
random measures defined on extended spaces. However, the simulation of
large and interacting populations is often referred as computationally
expensive.

The aim of the `R` package `IBMPopSim` is to meet the needs of the
various communities for efficient tools in order to simulate the
evolution of stochastic IBMs. `IBMPopSim` provides a general framework
for the simulation of a wide class of IBMs, where individuals are
characterized by their age and/or a set of characteristics. Different
types of events can be included in the modeling by users, depending on
their needs: births, deaths, entry or exit in/to the population and
changes of characteristics (swap events). Furthermore, the various
events that can happen to individuals in the population can occur at a
non-stationary frequency, depending on the individuals' characteristics
and time, and also including potential interactions between individuals.

We introduce a unified mathematical and simulation framework for this
class of IBMs, generalizing the pathwise representation of IBMs by
thinning of Poisson measures, as well as the associated population
simulation algorithm, based on an acceptance/rejection procedure. In
particular, we provide general sufficient conditions on the event
intensities under which the simulation of a particular model is
possible.

We opted to implement the algorithms of the `IBMPopSim` package using
the `Rcpp` package, a tool facilitating the seamless integration of
high-performance `C++` code into easily callable `R` functions
([@JSSv040i08]). With just a few lines of `C++` code, `IBMPopSim` offers
user-friendly R functions for defining IBMs. Once events and their
associated intensities are specified, an automated procedure creates the
model. This involves integrating the user's source code into the primary
`C++` code using a template mechanism. Subsequently, `Rcpp` is invoked
to compile the model and integrate it into the `R` session. Following
this process, the model becomes callable with varying parameters,
enabling the generation of diverse population evolution scenarios.\
Combined with the design of the simulation algorithms, the package
structure yields very competitive simulation runtimes for IBMs, while
staying user-friendly for `R` users. Several outputs function are also
implemented in `IBMPopSim`. For instance the package allows the
construction and visualization of age pyramids, as well as the
construction of death and exposures table from the censored individual
data, compatible with `R` packages concerned with mortality modelling,
such as [@Rdemography] or [@Rstmomo]. Several examples are provided in
the form of `R` vignettes on the
website <https://daphnegiorgi.github.io/IBMPopSim/>, and in recent works
of [@karoui2021simulating] and [@roget2022positive].

Designed for applications in social sciences, the `R` package
`MicSim` [@Zin14] can be used for continuous time microsimulation. In
continuous-time microsimulation, individual life-courses are usually
specified by sequences of state transitions (events) and the time spans
between these transitions. The state space is usually discrete and
finite, which is no necessarily the case in `IBMPopSim`, where
individuals can have continuous characteristics. But most importantly,
microsimulation does not allow for interactions between individuals.
Indeed, microsimulation produces separately the life courses of all
individuals in the populations, based on the computation of the
distribution functions of the waiting times in the distinct states of
the state space, for each individual ([@Zin14]). This can be slow in
comparison to the simulation by thinning of event times occurring in the
population, which is based on selecting event times among some competing
proposed event times. Finally, `MicSim` simplifies the Mic-Core
microsimulation tool implemented in Java ([@zinn2009mic]). However, the
implementation in `R` of simulation algorithms yields longer simulation
run times than when using `Rcpp`. To the best of our knowledge, there
are no other `R` packages currently available addressing the issue of
IBMs efficient simulation.

In Section [1](#section::IBM){reference-type="ref"
reference="section::IBM"}, we introduce the mathematical framework that
characterizes the class of Stochastic Individual-Based Models (IBMs)
that can be implemented in the `IBMPopSim` package. In particular, a
general pathwise representation of IBMs is presented. The population
dynamics is obtained as the solution of an SDE driven by Poisson
measures, for which we obtain existence and uniqueness results in
Theorem [1](#ThEqZ){reference-type="ref" reference="ThEqZ"}.
Additionally, a succinct overview of the package is provided. In
Section [2](#sec::simulation){reference-type="ref"
reference="sec::simulation"} the two main algorithms for simulating the
population evolution of an IBM across the interval $[0, T]$ are
detailed. In Section [3](#sec::package){reference-type="ref"
reference="sec::package"} we present the main functions of the
`IBMPopSim` package, which allow for the definition of events and their
intensities, the creation of a model, and the simulation of scenarios.
Two examples are detailed in
Sections [4](#SectionInsurancePortofio){reference-type="ref"
reference="SectionInsurancePortofio"}
and [5](#section:ExempleInteraction){reference-type="ref"
reference="section:ExempleInteraction"}, featuring applications
involving an heterogeneous insurance portfolio characterized by entry
and exit events, and an age and size-structured population with
intricate interactions.

# Stochastic Individual-Based Models (IBMs) in IBMPopSim {#section::IBM}

Stochastic Individual-Based Models (IBMs) represent a broad class of
random population dynamics models, allowing the description of
populations evolution on a microscopic scale. Informally, an IBM can be
summarized by the description of the individuals constituting the
population, the various types of events that can occur to these
individuals, along with their respective frequencies. In `IBMPopSim`,
individuals can be characterized by their age and/or a collection of
discrete or continuous characteristics. Moreover, the package enables
users to simulate efficiently populations in which one or more of the
following event types may occur:

-   **Birth event**: addition of an individual of age 0 to the
    population.

-   **Death event**: removal of an individual from the population.

-   **Entry event**: arrival of an individual in the population.

-   **Exit (emigration) event**: exit from the population (other than
    death).

-   **Swap event**: an individual changes characteristics.

Each event type is linked to an associated event kernel, describing how
the population is modified following the occurrence of the event. For
some event types, the event kernel requires explicit specification. This
is the case for entry events when a new individual joins the population.
Then,the model should specify how the age and characteristics of this
new individual are chosen. For instance, the characteristics of a new
individual in the population can be chosen uniformly in the space of all
characteristics, or can depend on the distribution of his parents or
those of the other individuals composing the population.

The last component of an IBM are the event intensities. Informally, an
event intensity is a function $\lambda^e_t(I, Z)$ describing the
frequency at which an event $e$ can occur to an individual $I$ in a
population $Z$ at a time $t$. Given a history of the population
$(\mathcal{F}_t)$, the probability of event $e$ occurring to individual
$I$ during a small interval of time $(t,t+dt]$ is proportional to
$\lambda^e(I,t)$:
$$\mathbb{P}(\text{event } e \text{ occurring to $I$ during } (t,t+dt] | \mathcal{F}_t) \simeq \lambda^e_t(I, Z)dt.$$
The intensity function $\lambda^e$ can include dependency on the
individual's $I$ age and characteristics, the time $t$, or the
population composition $Z$ in the presence of interactions.\

## Brief package overview

Prior to providing a detailed description of an Individual-Based Model
(IBM), we present a simple model of birth and death in an age-structured
"human" population. We assume no interactions between individuals, and
individuals are characterized by their gender, in addition to their age.
In this simple model, all individuals, regardless of gender, can give
birth when their age falls between 15 and 40 years, with a constant
birth rate of 0.05. The death intensity is assumed to follow a
Gompertz-type intensity depending on age. The birth and death
intensities are then given by
$$\lambda^b(t, I) = 0.05 \times \mathbf{1}_{[15,40]}(a(I,t)), \quad \lambda^d(I,t) = \alpha\exp(\beta a(I,t)),$$
with $a(I,t)$ the age of individual $I$ at time $t$. Birth events are
also characterized with a kernel determining the gender of the newborn,
who is male with probability $p_{male}$.

#### Model creation

To implement this model in `IBMPopSim`, it is necessary to individually
define each event type. In this example, the `mk_event_individual`
function is used. The creation of an event involves a few lines of `C++`
instructions defining the intensity and, if applicable, the kernel of
the event. For a more in depth description of the event creation step
and its parameters, we refer to
Section [3.2](#sec::package_events){reference-type="ref"
reference="sec::package_events"}.

The events of this simple model are for example defined through the
following calls. In the `C++` codes, the names `birth_rate`, `p_male`,
`alpha` and `beta` refer to the model parameters defined in the
following list.

In a second step, the model is created by calling the function
`mk_model`. A `C++` source code is automatically created through a
template mechanism based on the events and parameters, subsequently
compiled using the `sourceCpp` function from the `Rcpp` package.

#### Simulation

Once the model is created and compiled, the `popsim` function is called
to simulate the evolution of a population according to this model. To
achieve this, an initial population must be defined. In this example, we
extract a population from a dataset specified in the package (a sample
of $100\,000$ individuals based on the population of England and Wales
in 2014). It is also necessary to set bounds for the events intensities.
In this example, they are obtained by assuming that the maximum age for
an individual is 115 years. The data frame `sim_out$population` contains
the information (birth, death, gender) on individuals who lived in the
population over the period $[0,30]$. Functions of the package allows to
provide aggregated information on the population.

In the remainder of this section, we define rigorously the class of IBMs
that can be simulated in `IBMPopSim`, along with the assumptions that
are required in order for the population to be simulatable. The
representation of age-structured IBMs based on measure-valued processes,
as introduced in [@tran_2008], is generalized to a wider class of
abstract population dynamics. The modeling differs slightly here, since
individuals are "kept in the population" after their death (or exit), by
including the death/exit date as an individual trait.

## Population {#sec::population}

#### Notations

In the remainder of the paper, the filtered probability space is denoted
by $(\Omega,\{\mathcal{F}_t \},{\mathbb{P}})$, under the usual
assumptions. All processes are assumed to be càdlàg and adapted to the
filtration $\{\mathcal{F}_t \}$ (for instance the history of the
population) on a time interval $[0,T]$. For a càdlàg process $X$, we
denote $X_{t^-} := \lim_{\genfrac{}{}{0pt}{2}{s\to t}{s<t}} X_s$.

#### Individuals

An individual is represented by a triplet
$I = (\tau^b, \tau^d, x) \in \mathcal{I}= {\mathbb{R}}\times \bar {\mathbb{R}}\times {\mathcal{X}}$
with:

-   $\tau^b \in {\mathbb{R}}$ the date of birth,

-   $\tau^d \in \bar {\mathbb{R}}$ the death date, with
    $\tau^d = \infty$ if the individual is still alive,

-   a collection $x \in {\mathcal{X}}$ of characteristics where
    ${\mathcal{X}}$ is the space of characteristics.

Note that in IBMs, individuals are usually characterized by their age
$a(t) =t-\tau^b$ instead of their date of birth $\tau^b$. However, using
the latter is actually easier for the simulation, as it remains constant
over time.

#### Population process

The population at a given time $t$ is a random set
$$Z_t=\{ I_k \in \mathcal{I}; \; k= 1,\dots, N_t\},$$ composed of all
individuals (alive or dead) who have lived in the population before time
$t$. As a random set, $Z_t$ can be represented by a random counting
measure on $\mathcal{I}$ , that is an integer-valued measure
$Z: \Omega \times \mathcal{I}\to \bar {\mathbb{N}}$ where for
$A \in \mathcal{I}$, $Z(A)$ is the (random) number of individuals $I$ in
the subset $A$. With this representation: $$\begin{aligned}
\label{eq::popZ}
& Z_t (\mathrm{d}\tau^b, \mathrm{d}\tau^d , \mathrm{d}x)= \sum_{k=1}^{N_t} \delta_{I_k} (\tau^b, \tau^d,x), \nonumber\text{ with }  \int_{\mathcal{I}} f(\tau^b, \tau^d, x)Z_t (\mathrm{d}\tau^b, \mathrm{d}\tau^d , \mathrm{d}x) = \sum_{k=1}^{N_t} f(I_k).
\end{aligned}$$ The number of individuals present in the population
*before time* $t$ is obtained by taking $f\equiv 1$:
$$N_t =  \int_{\mathcal{I}}  Z_t(\mathrm{d}\tau^b, \mathrm{d}\tau^d, \mathrm{d}x) = \sum_{k=1}^{N_t} \boldsymbol{1}_{\mathcal{I}} (I_k).$$
Note that $(N_t)_{t\geq 0}$ is an increasing process since dead/exited
individuals are kept in the population $Z$. The number of alive
individuals in the population at time $t$ is: $$\label{eq:Nta}
N_t^a =  \int_{\mathcal{I}}  \mathsf{1}_{\{\tau^d > t \} }Z_t(\mathrm{d}\tau^b, \mathrm{d}\tau^d, \mathrm{d}x) = \sum_{k=1}^{N_t} \mathsf{1}_{\{\tau^d_k > t \} }.$$
Another example is the number of alive individuals of age over $a$ is
$$N_t([a,+\infty)) :=  \int_{\mathcal{I}}  \boldsymbol{1}_{[a,+\infty)}(t-\tau^b)\mathsf{1}_{]t,\infty]}(\tau^d) Z_t(\mathrm{d}\tau^b, \mathrm{d}\tau^d, \mathrm{d}x) = \sum_{k=1}^{N_t}  \boldsymbol{1}_{\{ t -\tau_k^b \geq a \}}\mathsf{1}_{\{\tau^d_k \geq t \} }.$$

## Events {#sec::events}

The population composition changes at random dates following different
types of events. `IBMPopSim` allows the simulation of IBMs with the
following events types:

-   *A birth* event at time $t$ is the addition of a new individual
    $I'=(t,\infty, X)$ of age $0$ to the population. Their date of birth
    is $\tau^b =t$, and characteristics is $X$, a random variable of
    distribution defined by the birth kernel $k^b(t,I,\mathrm{d}x)$ on
    ${\mathcal{X}}$, depending on $t$ and its parent $I$. The population
    size becomes $N_t = N_{t^-} + 1$, and the population composition
    after the event is $$Z_t  = Z_{t^-} +  \delta_{(t,\infty, X)}.$$

-   An *entry* event at time $t$ is also the addition of an individual
    $I'$ in the population. However, this individual is not of age $0$.
    The date of birth and characteristics of the new individual
    $I'= (\tau^b, \infty, X)$ are random variables of probability
    distribution defined by the entry kernel
    $k^{en}(t, \mathrm{d}s, \mathrm{d}x)$ on
    ${\mathbb{R}}\times {\mathcal{X}}$. The population size becomes
    $N_t = N_{t^-} + 1$, and the population composition after the event
    is: $$Z_t  = Z_{t^-} +  \delta_{(\tau^b, \infty, X)}.$$

-   *A death* or *exit* event of an individual
    $I= (\tau^b,\infty, x)\in Z_{t^-}$ at time $t$ is the modification
    of its death date $\tau^d$ from $+\infty$ to $t$. This event results
    in the simultaneous addition of the individual $(\tau^b,t,x)$ and
    removal of the individual $I$ from the population. The population
    size is not modified, and the population composition after the event
    is $$Z_t  = Z_{t^-} +\delta_{(\tau^b,t,x)}- \delta_{I}.$$

-   *A swap* event (change of characteristics) results in the
    simultaneous addition and removal of an individual. If an individual
    $I= (\tau^b,\infty, x) \in Z_{t^-}$ changes of characteristics at
    time $t$, then it is removed from the population and replaced by
    $I' = (\tau^b,\infty, X)$. The new characteristics $X$ is a random
    variable of distribution $k^s(t, I,\mathrm{d}x)$ on ${\mathcal{X}}$,
    depending on time, the individual's age and previous characteristics
    $x$. In this case, the population size is not modified and the
    population becomes:
    $$Z_t  = Z_{t^-}   +  \delta_{(\tau^b,  \infty, X)} -  \delta_{(\tau^b, \infty, x)}.$$

To summarize, the space of event types is $E = \{ b, en, d, s \}$, and
the jump $\Delta Z_t = Z_t - Z_{t^-}$ (change in the population
composition) generated by an event of type $e \in \{ b, en, d, s \}$ is
denoted by $\phi^e(t, I)$, with:\

:::: center
::: {#TableEvAction}
     Event      Type                           $\phi^e(t, I)$                                              New individual
  ------------ ------ ---------------------------------------------------------------- -------------------------------------------------------
     Birth      $b$                      $\delta_{(t, \infty,  X)}$                          $\tau^b =t, \; X \sim k^b(t,I,\mathrm{d}x)$
     Entry      $en$                  $\delta_{(\tau^b, \infty,  X)}$                   $(\tau^b, X) \sim k^{en}(t,\mathrm{d}s, \mathrm{d}x)$
   Death/Exit   $d$       $\delta_{(\tau^b, t,x)} - \delta_{(\tau^b, \infty, x)}$                           $\tau^d = t$
      Swap      $s$    $\delta_{(\tau^b, \infty , X)} - \delta_{(\tau^b, \infty, x)}$               $X \sim k^s(t,I,\mathrm{d}x)$

  : Events action
:::
::::

::: {#remark::popfinale .remark}
*Remark 1*.

-   At time $T$, the population $Z_T$ contains all individuals who lived
    in the population before $T$, including dead/exited individuals. If
    there are no swap events, or entries,the population state $Z_t$ for
    any time $t\leq T$ can be obtained from $Z_T$. Indeed, if
    $Z_T = \sum_{k=1}^{N_T}  \delta_{I_k}$, then the population at time
    $t\leq T$ is simply composed of the individuals born before $t$:
    $$Z_t  = \sum_{k=1}^{N_T} \boldsymbol{1}_{\{\tau^b_k \leq t \}}   \delta_{I_k}.$$

-   In the presence of entries (open population), a characteristic $x$
    can track the individuals' entry dates. Then, the previous equation
    can be easily modified in order to obtain the population $Z_t$ at
    time $t\leq T$ from $Z_T$.
:::

## Events intensity {#sec::event_intensity}

Once the different event types have been defined in the population
model, the frequency at which each event occur in the population $e$
have to be specified.\
Informally, the intensity $\Lambda^e_t(Z_t)$ at which an event $e$ can
occur is defined by
$$\mathbb P\big( \text{event } e \text { occurs in the population }  Z_t  \in (t,t+\mathrm{d}t] | \mathcal{F}_t \big) \simeq  \Lambda^e_t (Z_t)\mathrm{d}t.$$
For a more formal definition of stochastic intensities, we refer to
[@bremaud1981point] or [@KaaElK20].\
The form of the intensity function $(\Lambda^e_t (Z_t))$ determines the
population simulation algorithm in `IBMPopSim`:

-   When the event intensity does not depend on the population state,
    $$\label{PoissonIntensity}
    (\Lambda^e_t (Z_t))_{t\in [0,T]} = (\mu^e(t))_{t \in [0,T]},$$ with
    $\mu^e$ a deterministic function, the events of type $e$ occur at
    the jump times of an inhomogeneous Poisson process of intensity
    function $(\mu^e(t))_{t \in [0,T]}$. When such an event occurs, the
    individual to whom the event happens to is drawn uniformly among
    alive individuals in the population.\
    In a given model, the set of events $e\in E$ with Poisson
    intensities will be denoted by $\mathcal{P}$.

-   Otherwise, we assume that the global intensity $\Lambda^e_t(Z_t)$ at
    which the events of type $e$ occur in the population can be written
    as the sum of individual intensities $\lambda^e_t(I,Z_t)$:
    $$\begin{aligned}
    \label{eq:GlobalIntensity}
    &\Lambda^e_t (Z_t) = \sum_{k=1}^{N_t} \lambda^e_t ( I_k,Z_t),  \\
    & \nonumber \text{with } \mathbb P\big( \text{event } e \text { occurs to an individual } I \in (t,t+\mathrm{d}t] | \mathcal{F}_t \big) \simeq  \lambda^e_t (I,Z_t)\mathrm{d}t.
    \end{aligned}$$

Obviously, nothing can happen to dead or exited individuals, i.e.
individuals $I= (\tau^b, \tau^d, x)$ with $\tau^d \leq t$. Thus,
individual event intensities are assumed to be null for dead/exited
individuals:
$$\lambda^e_t ( I,Z_t) = 0, \text{ if }\tau^d \leq t, \text{ so that } \Lambda^e_t (Z_t) = \sum_{k=1}^{N_t^a} \lambda^e_t ( I_k,Z_t),$$
with $N^a_t$ the number of alive individuals at time $t$.

The event's individual intensity $\lambda^e_t (I,Z_t)$ can depend on
time (for instance when there is a mortality reduction over time), on
the individual's age $t-\tau^b$ and characteristics, but also on the
population composition $Z_t$. The dependence of $\lambda^e$ on the
population $Z$ models interactions between individuals in the
populations. Hence, two types of individual intensity functions can be
implemented in `IBMPopSim`:

1.  *No interactions:* The intensity function $\lambda^e$ does not
    depend on the population composition. The intensity at which the
    event of type $e$ occur to an individual $I$ only depends on its
    date of birth and characteristics:
    $$\label{eq:intensityNointeraction}
    \lambda^e_t (I,Z_t) = \lambda^e(t, I),$$ where
    $\lambda^e: \mathbb{R}_+ \times \mathcal{I}\to {\mathbb{R}}^+$ is a
    deterministic function. In a given model, we denote by $\mathcal{E}$
    the set of event types with individual intensity
    [\[eq:intensityNointeraction\]](#eq:intensityNointeraction){reference-type="eqref"
    reference="eq:intensityNointeraction"}.

2.  *"Quadratic" interactions:* The intensity at which an event of type
    $e$ occurs to an individual $I$ depends on $I$ and on the population
    composition, through an interaction function $W^e$. The quantity
    $W^e(t, I,J)$ describes the intensity of interactions between two
    alive individuals $I$ and $J$ at time $t$, for instance in the
    presence of competition or cooperation. In this case, we have
    $$\label{eq:intensityInteraction}
    \lambda^e_t(I,Z_t)=\sum_{j=1}^{N_t} W^e(t, I, I_j) = \int_{\mathcal{I}} W^e(t, I, (\tau^b,\tau^d,x)) Z_t (\mathrm{d}\tau^b,\mathrm{d}\tau^d, \mathrm{d}x),$$
    where $W^e(t, I, (\tau^b,\tau^d,x))  = 0$ if the individual
    $J =(\tau^b,\tau^d,x)$ is dead, i.e. $\tau^d \leq t$.\
    In a given model, we denote by $\mathcal{E}_W$ the set of event
    types with individual intensity
    [\[eq:intensityInteraction\]](#eq:intensityInteraction){reference-type="eqref"
    reference="eq:intensityInteraction"}.

To summarize, an individual intensity in IBMPopSim can be written as:
$$\label{IndividualIntensity}
\lambda^e_t(I,Z_t) = \lambda^e(t, I) \mathbf{1}_{\{e \in \mathcal{E}\}} + \biggl( \sum_{j=1}^{N_t} W^e(t, I, I_j) \biggr) \mathbf{1}_{\{e \in \mathcal{E}_W\}}.$$

*Examples*\
*(i)* An example of death intensity without interaction for an
individual $I=(\tau^b, \tau^d, x)$ alive at time $t$ ($t <\tau^d$) is:
$$\lambda^d (t,I) =  \alpha_x \exp(\beta_x a(I,t)), \text{ where }  a(I,t) = t-\tau^b$$
is the age of the individual $I$ at time $t$. In this classical case,
the death rate of an individual $I$ is an exponential (Gompertz)
function of the individual's age, with coefficients depending on the
individual's characteristics $x$.\
*(ii)* In the presence of competition between individuals, the death
intensity of an individual $I$ also depend on other individuals $J$ in
the population. For example, if $I=(\tau^b,\tau^d, x)$, with $x$ its
size, then we can have: $$\label{ex:interaction}
W^d(t,I,J) = (x_J - x)^+ \mathbf{1}_{\{\tau^d_J > t\}}, \quad \forall \; J=(\tau^b_J,\tau^d_J , x_J).$$
This can be interpreted as follows: if the individual $I$ meets randomly
an individual $J$ alive at time $t$, and of bigger size $x_J > x$, then
he can die at the intensity $x_J-x$. If $J$ is smaller than $I$, then he
cannot kill $I$. The bigger is the size $x$ of $I$, the lower is his
death intensity $\lambda^d_t(I,Z_t)$ defined by
$$\lambda^d_t(I,Z_t) = \sum_{\genfrac{}{}{0pt}{2}{J\in Z_t,}{x_J > x}} (x_J -x)\mathbf{1}_{\{\tau^d_J > t\}}.$$
*(iii)* `IBMPopSim` can simulate IBMs that include intensities expressed
as a sum of Poisson intensities and individual intensities of the form
$\Lambda^e(Z_t) =\mu^e_t + \sum_{k=1}^{N_t} \lambda^e(I_k, Z_t)$.\
Other examples are provided in Section
[4](#SectionInsurancePortofio){reference-type="ref"
reference="SectionInsurancePortofio"} and Section
[5](#section:ExempleInteraction){reference-type="ref"
reference="section:ExempleInteraction"}.\
Finally, the global intensity at which an event can occur in the
population is defined by: $$\label{eq:globalEvIntensity}
\Lambda_t(Z_t) = \sum_{e\in \mathcal{P}} \mu^e(t) + \sum_{e \in \mathcal E} \Big(\sum_{k=1}^{N_t} \lambda^e(t, I_k)\Big) + \sum_{e \in \mathcal E_W} \Big(\sum_{k=1}^{N_t}\sum_{j=1}^{N_t} W^e(t, I_k, I_j)\Big).$$
An important point is that for events $e \in \mathcal E$ without
interactions, the global event intensity
$\Lambda^e_t(Z_t) = \sum_{k=1}^{N_t} \lambda^e(t, I_k)$ is "of order"
$N_t^a$ defined in [\[eq:Nta\]](#eq:Nta){reference-type="eqref"
reference="eq:Nta"} (number of alive individuals at time $t$). On the
other hand, for events $e \in \mathcal{E}_W$ with interactions,
$\Lambda^e_t(Z_t) = \sum_{k=1}^{N_t}\sum_{j=1}^{N_t} W^e(t, I_k, I_j)$
is of order $(N_t^a)^2$. Informally, this means that when the population
size increases, events with interaction are more costly to simulate.
Furthermore, the numerous computations of the interaction kernel $W^e$
can also be quite costly. The randomized Algorithm
[\[algo::rzndomized\]](#algo::rzndomized){reference-type="ref"
reference="algo::rzndomized"}, detailed in Section
[2.3](#sec::simulation_algo_randomized){reference-type="ref"
reference="sec::simulation_algo_randomized"}, allows us to overcome
these limitations.

#### Events intensity bounds

The simulation algorithms implemented in `IBMPopSim` are based on an
acceptance/rejection procedure, which requires to specify bounds for the
various events intensities $\Lambda^e_t(Z_t)$. These bounds are defined
differently depending on the expression of the intensity.

::: {#AssumptionIntensityPoisson .assumption}
**Assumption 1**. *For all events $e \in \mathcal{P}$ with Poisson
intensity
[\[PoissonIntensity\]](#PoissonIntensity){reference-type="eqref"
reference="PoissonIntensity"}, the intensity is assumed to be bounded on
$[0,T]$:
$$\forall t \in [0,T], \quad \Lambda^e_t(Z_t) = \mu^e(t) \leq \bar \mu^e.$$*
:::

When $e \in \mathcal{E} \cup \mathcal{E}_W$
($\Lambda^e_t(Z_t) =\sum_{k=1}^{N_t} \lambda^e_t(I_k,Z_t)$), assuming
that $\Lambda^e_t(Z_t)$ is uniformly bounded is too restrictive since
the event intensity depends on the population size. In this case, the
assumption is made on the individual intensity or interaction function
$W^e$, depending on the situation:

::: {#AssumptionIntensity1 .assumption}
**Assumption 2**. *For all event types $e \in \mathcal{E}$, the
associated individual event intensity $\lambda^e$ with no interactions
($\lambda^e$ verifies
[\[eq:intensityNointeraction\]](#eq:intensityNointeraction){reference-type="eqref"
reference="eq:intensityNointeraction"}) is assumed to be uniformly
bounded:
$$\lambda^e(t, I) \leq \bar \lambda^e, \quad \forall \;  t\in [0, T],  \;   I \in \mathcal{I}.$$
In particular, $$\label{EqdefbarLambda}
\forall t \in [0,T], \quad \Lambda^e_t (Z_t) = \sum_{k=1}^{N_t} \lambda^e(t, I) \leq \bar \lambda^e  N_t .$$*
:::

::: {#AssumptionIntensity2 .assumption}
**Assumption 3**. *For all event types $e \in \mathcal{E}_W$, the
associated interaction function $W^e$ is assumed to be uniformly
bounded:
$$W^e(t, I, J) \leq \bar W^e, \quad \forall \; t\in [0,T], \;   I, J \in \mathcal{I}.$$
In particular, $\forall t \in [0,T]$,
$$\lambda^e_t (I,Z_t) = \sum_{j=1}^{N_t} W^e(t, I, I_j)  \leq  \bar W^e N_t, \quad \text{and} \quad \Lambda^e_t (Z_t) \leq \bar W^e (N_t)^2.$$*
:::

Assumptions [1](#AssumptionIntensityPoisson){reference-type="ref"
reference="AssumptionIntensityPoisson"},
[2](#AssumptionIntensity1){reference-type="ref"
reference="AssumptionIntensity1"} and
[3](#AssumptionIntensity2){reference-type="ref"
reference="AssumptionIntensity2"} yield that events in the population
occur with the global event intensity $\Lambda_t(Z_t)$
[\[eq:globalEvIntensity\]](#eq:globalEvIntensity){reference-type="eqref"
reference="eq:globalEvIntensity"}, which is dominated by a polynomial
function in the population size: $$\label{eq:defbarLambda}
\Lambda_t(Z_t) \leq \bar \Lambda(N_t), \quad \text{with }  \bar \Lambda (n) = \sum_{e \in \mathcal{P}} \bar \mu^e + \sum_{e\in \mathcal{E}}\bar \lambda^e  n + \sum_{e \in \mathcal E_W} \bar W^e n^2.$$
This bound is linear in the population size if there are no
interactions, and quadratic if there at least is an event including
interactions. This assumption is the key to the algorithms implemented
in `IBMPopSim`. Before presenting the simulation algorithm, we close
this section with a rigorous definition of an IBM, based on the pathwise
representation of its dynamics a Stochastic Differential Equation (SDE)
driven by Poisson random measures.\

## Pathwise representation

Since the seminal paper of [@FouMel04], it has been shown in many
examples that a stochastic IBM dynamics can be defined rigorously as the
unique solution of an SDE driven by Poisson measures, under reasonable
non explosion conditions. In the following, we introduce a unified
framework for the pathwise representation of the class of stochastic
IBMs introduced above. Some recalls on Poisson random measures are
presented in the Appendix
[6](#section::preliminaries){reference-type="ref"
reference="section::preliminaries"}, and for more details on these
representations on particular examples, we refer to the abundant
literature on the subject.

In the following we consider an individual-based stochastic population
$(Z_t)_{t\in [0,T]}$, keeping the notations introduced in Section
[1.3](#sec::events){reference-type="ref" reference="sec::events"} and
[1.4](#sec::event_intensity){reference-type="ref"
reference="sec::event_intensity"} for the events and their intensities.
In particular, the set of events types that define the population
evolution is denoted by
$\mathcal{P} \cup \mathcal{E} \cup \mathcal{E}_W \subset E$, with
$\mathcal{P}$ the set of events types with Poisson intensity verifying
assumption [1](#AssumptionIntensityPoisson){reference-type="ref"
reference="AssumptionIntensityPoisson"}, $\mathcal{E}$ the set of events
types with individual intensity and no interaction, verifying Assumption
[2](#AssumptionIntensity1){reference-type="ref"
reference="AssumptionIntensity1"}, and finally $\mathcal{E}_W$ the set
of event types with interactions, verifying Assumption
[3](#AssumptionIntensity2){reference-type="ref"
reference="AssumptionIntensity2"}.

#### Non explosion criterion

First, one has to ensure that the number of events occurring in the
population will not explode in finite time, leading to an infinite
simulation time. Assumptions
[2](#AssumptionIntensity1){reference-type="ref"
reference="AssumptionIntensity1"} and
[3](#AssumptionIntensity2){reference-type="ref"
reference="AssumptionIntensity2"} are not sufficient to guarantee the
non explosion of the event number, due to the potential explosion of the
population size in the presence of interactions. An example is the case
when only birth events occur, with an intensity
$\Lambda^b_t(Z_t) = C_b (N_t^a)^2$ ($W^b(t, I,J) =C_b$). Then, the
number of alive individuals $(N_t^a)_{t\geq 0}$ is a well-known pure
birth process of intensity function $g(n) = C_b n^2$ (intensity of
moving from state $n$ to $n+1$). This process explodes in finite time,
since $g$ does not verify the necessary and sufficient non explosion
criterion for pure birth Markov processes:
$\sum_{n=1}^\infty \frac{1}{g(n)} = \infty$ (see e.g. Theorem 2.2 in
[@BanMel15]). There is thus an explosion in finite time of birth events.

This example shows that the important point for non explosion is to
control the population size. We give below a general sufficient
condition on birth and entry event intensities, in order for the
population size to stay finite in finite time. This ensures that the
number of events does not explode in finite time. Informally, the idea
is to control the intensities by a pure birth intensity function
verifying the non-explosion criterion.

::: {#Assumption:nonExplosion .assumption}
**Assumption 4**. *Let $e=b$ or $en$, a birth or entry event type. If
the intensity at which the events of type $e$ occur in the population
are not Poissonian, i.e. $e \in \mathcal{E} \cup \mathcal{E}_W$, then
there exists a function $f^e : {\mathbb{N}}\to (0, +\infty)$, such that
$$\sum_{n=1}^{\infty} \frac{1}{nf^e(n)} = \infty,$$ and for all
individual $I \in \mathcal{I}$ and population measure
$Z = \sum_{k=1}^{n} \delta_{I_k}$ of size $n$,
$$\lambda^e_t (I, Z) \leq f^e(n), \; \forall \; 0\leq t \leq T.$$*
:::

::: remark
*Remark 2*. If $e \in \mathcal{E}$,
$\lambda_t^e(I,Z) = \lambda^e(t,I) \leq \bar{\lambda}^e$ by the
domination Assumption [3](#AssumptionIntensity2){reference-type="ref"
reference="AssumptionIntensity2"}. In this case, Assumption
[4](#Assumption:nonExplosion){reference-type="ref"
reference="Assumption:nonExplosion"} is always verified with
$f^e(n) = \bar{\lambda}^e$.
:::

Assumption [4](#Assumption:nonExplosion){reference-type="ref"
reference="Assumption:nonExplosion"} yields that the global intensity
$\Lambda_t^e(\cdot)$ of event $e$ is bounded by a function $g^e$ only
depending on the population size:
$$\Lambda_t^e (Z) \leq g^e(n) := nf^e(n), \quad \text{with }\sum_{n=1}^{\infty} \frac{1}{g^e(n)} = \infty.$$
If $e\in \mathcal{P}$ has a Poisson intensity, then
$\Lambda_t^e(Z) =\mu^e_t$ always verifies the previous equation with
$g^e(n) = \bar \mu^e$.

####  {#section .unnumbered}

Before introducing the IBM SDE, let us give an idea of the equation
construction. Between two successive events, the population composition
$Z_t$ stays constant, since the population process $(Z_t)_{t \geq 0}$ is
a pure jump process. Furthermore, since each event type is characterized
by an intensity function, the jumps occurring in the population can be
represented by restriction and projection of a Poisson measure defined
on a larger state space. More precisely, we introduce a random Poisson
measure $Q$ on $\mathbb R^+ \times \mathcal{J}\times \mathbb{R}^+$, with
$\mathcal{J}= \mathbb N \times(\mathcal E \cup \mathcal{E}_W)$. $Q$ is
composed of random quadruplets $(\tau, k , e, \theta)$, where $\tau$
represents a potential event time for an individual $I_k$ and event type
$e$. The last variable $\theta$ is used to accept/reject this proposed
event, depending on the event intensity. Hence, the Poisson measure is
restricted to a certain random set and then projected on the space of
interest ${\mathbb{R}}^+ \times \mathcal{J}$. If the event is accepted,
then a jump $\phi^e(\tau,I_k)$ occurs.

The proof of Theorem [1](#ThEqZ){reference-type="ref" reference="ThEqZ"}
is detailed in the Appendix [7.1](#ProofThPathwise){reference-type="ref"
reference="ProofThPathwise"}. Note that Equation
[\[SDE_pop\]](#SDE_pop){reference-type="eqref" reference="SDE_pop"} is
an SDE describing the evolution of the IBM, the intensity of the events
in the right hand side of the equation depending on the population
process $Z$ itself. The main idea of the proof of Theorem
[1](#ThEqZ){reference-type="ref" reference="ThEqZ"} is to use the non
explosion property of Lemma
[2](#lemma:nonExplosionSDE){reference-type="ref"
reference="lemma:nonExplosionSDE"}, and to write the r.h.s of
[\[SDE_pop\]](#SDE_pop){reference-type="eqref" reference="SDE_pop"} as a
sum of simple equations between two successive events, solved by
induction.\
The proof of Lemma [2](#lemma:nonExplosionSDE){reference-type="ref"
reference="lemma:nonExplosionSDE"}, detailed in Appendix
[7.2](#section:proofLemma){reference-type="ref"
reference="section:proofLemma"}, is more technical and rely on pathwise
comparison result, generalizing those obtained in [@KaaElK20]. An
alternative pathwise representation of the population process, inspired
by the randomized Algorithm
[\[algo::rzndomized\]](#algo::rzndomized){reference-type="ref"
reference="algo::rzndomized"}, is given as well in Theorem
[12](#ThEqZrandomized){reference-type="ref"
reference="ThEqZrandomized"}.

::: {#ThEqZ .theo}
**Theorem 1** (Pathwise representation). *Let $T\in {\mathbb{R}}^+$ and
$\mathcal{J}= \mathbb N \times(\mathcal E \cup \mathcal{E}_W)$.\
Let $Q$ be a random Poisson measure on
$\mathbb R^+ \times \mathcal{J}\times \mathbb{R}^+$, of intensity
$\mathrm{d}t \delta_{\mathcal{J}}(\mathrm{d}k,\mathrm{d}e)  \mathbf{1}_{[0,\bar \lambda^e]} (\theta)\mathrm{d}\theta$,
with $\delta_\mathcal{J}$ the counting measure on $\mathcal{J}$.
Finally, let $Q^{\mathcal P}$ be a random Poisson measure on
$\mathbb R^+ \times \mathcal{P}  \times \mathbb{R}^+$, of intensity
$\mathrm{d}t \delta_{\cal P}(\mathrm{d}e)  \mathbf{1}_{[0,\bar \mu^e]} (\theta)\mathrm{d}\theta$,
and $Z_0= \sum_{k=1}^{N_0} \delta_{I_k}$ an initial population.\
Then, under Assumption
[4](#Assumption:nonExplosion){reference-type="ref"
reference="Assumption:nonExplosion"}, there exists a unique
measure-valued population process $Z$, strong solution on the following
SDE driven by the Poisson measure $Q$: $$\begin{aligned}
\label{SDE_pop}
Z_t = Z_0 &+ \int_0^t \int_{ \mathcal{J}\times \mathbb R^+ }\phi^e (s , I_k)  \mathbf{1}_{\{k \leq N_{s^-}\} }\mathbf{1}_{\{\theta \leq \lambda_s^e(I_k, Z_{s^-})\}} Q (\mathrm{d}s ,\mathrm{d}k , \mathrm{d}e, \mathrm{d}\theta )\\
\nonumber &+   \int_0^t \int_{\mathcal{P} \times \mathbb R^+}  \phi^e(s, I_{s^-}) \mathbf{1}_{\{\theta \leq \mu^e(s) \}} Q^{\mathcal{P}} (\mathrm{d}s ,\mathrm{d}e,  \mathrm{d}\theta),  \qquad \forall  0 \leq t \leq T,
\end{aligned}$$ and where $I_{s^-}$ is an individual, chosen uniformly
among alive individuals in the population $Z_{s^-}$.*
:::

::: {#lemma:nonExplosionSDE .lemma}
**Lemma 2**. *Let $Z$ be a solution of
[\[SDE_pop\]](#SDE_pop){reference-type="eqref" reference="SDE_pop"} on
${\mathbb{R}}^+$, with $(T_n)_{n\geq 0}$ its jump times, $T_0 = 0$. If
Assumption [4](#Assumption:nonExplosion){reference-type="ref"
reference="Assumption:nonExplosion"} is satisfied, then
$$\lim_{n \to \infty} T_n = \infty, \quad {\mathbb{P}}\text{-a.s.}$$*
:::

# Population simulation {#sec::simulation}

We now present the main algorithm for simulating the evolution of an IBM
over $[0,T]$.The algorithm implemented in `IBMPopSim` allows the exact
simulation of [\[SDE_pop\]](#SDE_pop){reference-type="eqref"
reference="SDE_pop"}, based on an acceptance/reject algorithm for
simulating random times called *thinning*. The exact simulation of event
times with this acceptance/reject procedure is closely related to the
simulations of inhomogeneous Poisson processes by the so-called thinning
algorithm, often attributed to [@LewShe79]. The simulation methods for
inhomogeneous Poisson processes can be adapted to IBMs, and we introduce
in this section a general algorithm extending those by [@FouMel04] (see
also [@FerTra09], [@Ben10]).

The algorithm is based on exponential "candidate" event times, chosen
with a (constant) intensity which must be greater than the global event
intensity $(\Lambda_t(Z_t))$
[\[eq:GlobalIntensity\]](#eq:GlobalIntensity){reference-type="eqref"
reference="eq:GlobalIntensity"}. Starting from time $t$, once a
candidate event time $t + \bar T_\ell$ has been proposed, a candidate
event type $e$ (birth, death,\...) is chosen with a probability $p^e$
depending on the event intensity bounds $\bar \mu^e$, $\bar \lambda^e$
and $\bar W^e$, as defined in Assumption
[2](#AssumptionIntensity1){reference-type="ref"
reference="AssumptionIntensity1"} and
[3](#AssumptionIntensity2){reference-type="ref"
reference="AssumptionIntensity2"}. An individual $I$ is then drawn from
the population. Finally, it remains to accept or reject the candidate
event with a probability $q^e(t,I,Z_t)$ depending on the true event
intensity. If the candidate event time is accepted, then the event $e$
occurs at time $t + \bar T_\ell$ to the individual $I$. The main idea of
the algorithm implemented can be summarized as follows:

1.  Draw a candidate time $t + \bar T_\ell$ and candidate event type
    $e$.

2.  Draw a uniform variable $\theta \sim \mathcal{U}([0, 1])$ and
    individual $I$.

3.  **If** $\theta \leq q^e(t,I,Z_t)$ **then** event $e$ occur to
    individual $I$, **else** Do nothing and start again from
    $t + \bar T_\ell$.

Before introducing the main algorithms in more details, we recall
briefly the thinning procedure for simulating inhomogeneous Poisson
processes, as well as the links with pathwise representations. Some
recalls on Poisson random measures are presented in Appendix
[6](#section::preliminaries){reference-type="ref"
reference="section::preliminaries"}. For a more general presentation of
thinning of a Poisson random measure, see [@Dev86; @Cin11; @Kal17].

## Thinning of Poisson measure {#sec::thinning}

Let us start with the simulation and pathwise representation of an
inhomogeneous Poisson process on $[0,T]$ with intensity
$(\Lambda(t))_{t\in [0,T]}$. The thinning procedure is based on the
fundamental assumption that $\Lambda(t) \leq  \bar \Lambda$ is bounded
on $[0,T]$. In this case, the inhomogeneous Poisson can be obtained from
an homogeneous Poisson process of intensity $\bar \Lambda$, which can be
simulated easily (see Appendix
[6](#section::preliminaries){reference-type="ref"
reference="section::preliminaries"}).

First, the Poisson process can be extended to a Marked Poisson measure
$\bar Q:= \sum_{\ell \ge 1} \delta_{(\bar T_\ell, \bar \Theta_\ell)}$ on
$(\mathbb{R}^+)^2$, defined as follow:

-   The jump times of $(\bar T_\ell)_{\ell \ge 1}$ of $\bar Q$ are the
    jump times of a Poisson process of intensity $\bar \Lambda$.

-   The marks $(\bar \Theta_\ell)_{\ell \ge 1}$ are *i.i.d.* random
    variables, uniformly distributed on $[0,\bar \Lambda]$.

By Proposition [11](#PropMarkedPoisson){reference-type="ref"
reference="PropMarkedPoisson"}, $\bar{Q}$ is a Poisson random measure
with mean measure
$$\bar \mu(\mathrm{d}t, \mathrm{d}\theta): = \bar \Lambda \mathrm{d}t
    \frac{\mathbf{1}_{[0, \bar \Lambda]}(\theta)}{\bar \Lambda} \mathrm{d}\theta= \mathrm{d}t  \mathbf{1}_{[0, \bar \Lambda]}(\theta) \mathrm{d}\theta.$$
In particular, the average number of atoms
$(\bar T_\ell, \bar \Theta_\ell)$ in $[0,t]\times [0,h]$ is
$${\mathbb{E}}[Q([0,t]\times [0,h])]={\mathbb{E}}[\sum_{\ell} \boldsymbol{1}_{[0,t]\times [0,h]} (\bar T_\ell, \bar \Theta_{\ell})]  = \int_{(\mathbb{R}^+)^2}  \bar \mu(\mathrm{d}t, \mathrm{d}\theta)  = t (\bar \Lambda \wedge h).$$
The thinning is based on the restriction property for Poisson measure:
for a measurable set
$\Delta\subset {\mathbb{R}}^+\times {\mathbb{R}}^+$, the restriction
$Q^\Delta:= \boldsymbol{1}_{\Delta}\bar Q$ of $\bar Q$ to $\Delta$ (by
taking only atoms in $\Delta$) is also a Poisson random measure of mean
measure
$\mu^{\Delta}(\mathrm{d}t, \mathrm{d}\theta)  = \boldsymbol{1}_{\Delta}(t,\theta) \bar \mu(\mathrm{d}t, \mathrm{d}\theta).$

In order to obtain an inhomogeneous Poisson measure of intensity
$(\Lambda(t))$, the "good" choice of $\Delta$ is the hypograph of
$\Lambda$:
$\Delta =\{ (t,\theta) \in [0,T]\times [0,\bar \Lambda] ; \; \theta \leq \Lambda(t)\}$
(see Figure [1](#plot:thinning){reference-type="ref"
reference="plot:thinning"}). Then, $$\begin{aligned}
  &  Q^\Delta = \sum_{\ell \ge 1} \mathbf{1}_{\left\{\bar \Theta_\ell \le \Lambda(\bar T_\ell)\right\}} \delta_{(\bar T_\ell, \bar \Theta_\ell)},
 
\end{aligned}$$ and since $\Lambda(t) \leq \bar \Lambda$, on $[0,T]$:
$$\begin{aligned}
    \mu^{\Delta}(\mathrm{d}t, \mathrm{d}\theta) & = \boldsymbol{1}_{\{ \theta \leq \Lambda(t)\}}  \mathrm{d}t
\mathbf{1}_{[0, \bar \Lambda]}(\theta)\mathrm{d}\theta = \boldsymbol{1}_{\{\theta \leq \Lambda(t)\}} \mathrm{d}t \mathrm{d}\theta.
\end{aligned}$$

<figure id="plot:thinning">

<figcaption>Realization of a Marked Poisson measure <span
class="math inline"><em>Q̄</em></span> on <span
class="math inline">[0,<em>T</em>]</span> with mean measure <span
class="math inline"><em>μ̄</em>(d<em>t</em>,d<em>θ</em>) = d<em>t</em><strong>1</strong><sub>[0,<em>Λ̄</em>]</sub>(<em>θ</em>)d<em>θ</em></span>
(red crosses), and realization of the restriction <span
class="math inline"><em>Q̄</em><sup><em>Δ</em></sup></span> where <span
class="math inline"><em>Δ</em> = {(<em>t</em>,<em>θ</em>) ∈ [0,<em>T</em>] × [0,<em>Λ̄</em>], <em>θ</em> ≤ <em>Λ</em>(<em>t</em>)}</span>
(blue circles). The projection of <span
class="math inline"><em>Q̄</em><sup><em>Δ</em></sup></span> on first
component is an inhomogeneous Poisson process on <span
class="math inline">[0,<em>T</em>]</span> of intensity <span
class="math inline">(<em>Λ</em>(<em>t</em>))</span> and jump times <span
class="math inline">(<em>T</em><sub><em>k</em></sub>)<sub><em>k</em> ≥ 1</sub></span>.</figcaption>
</figure>

Finally, the inhomogeneous Poisson process is obtained by the projection
Proposition [10](#PropProjPoissonMeasure){reference-type="ref"
reference="PropProjPoissonMeasure"}, which states that the jump times of
$Q^\Delta$ are the jump times of an inhomogeneous Poisson process of
intensity $(\Lambda(t))$:

::: {#PropThinning .prop}
**Proposition 3**. *The counting process $N^\Lambda$, projection of
$Q^{\Delta}$ on the time component and defined by,
$$\label{EqThinInhomoeneousPoisson}
N^{\Lambda}_t := Q^{\Delta}( [0,t] \times {\mathbb{R}}^+) = \int_0^t \int_{{\mathbb{R}}^+} \boldsymbol{1}_{\{\theta \leq \Lambda(s)\}} \bar Q(\mathrm{d}s, \mathrm{d}\theta) = \sum_{\ell \geq 1} \mathsf{1}_{\{ \bar T_\ell \leq t \}} \mathsf{1}_{\{\bar \Theta_\ell \leq  \Lambda(\bar T_\ell) \}},  \quad \forall t \in [0,T],$$
is an inhomogeneous Poisson process on $[0,T]$ of intensity function
$(\Lambda(t))_{t\in [0,T]}$. The thinning
Equation [\[EqThinInhomoeneousPoisson\]](#EqThinInhomoeneousPoisson){reference-type="eqref"
reference="EqThinInhomoeneousPoisson"} is a pathwise representation of
$N^\Lambda$ by *restriction and projection* of the Poisson measure $Q$
on $[0,T]$.*
:::

The previous proposition yields a straightforward thinning algorithm to
simulate the jump times $(T_k)_{k \ge 1}$ of an inhomogeneous Poisson
process of intensity $\Lambda(t)$, by selecting jump times $\bar T_\ell$
such that $\bar \Theta_\ell \le \Lambda(\bar T_\ell)$.

#### Multivariate Poisson process {#sec:vector}

This can be extended to the simulation of multivariate inhomogeneous
Poisson processes, which is an important example before tackling the
simulation of an IBM.

Let $(N^j)_{j \in \mathcal{J}}$ be a (inhomogeneous) multivariate
Poisson process indexed by a finite set $\mathcal{J}$, such that
$\forall j \in \mathcal{J}$, the intensity $(\lambda_j(t))_{t\in [0,T]}$
of $N_j$ is bounded on $[0,T]$:
$$\sup_{t\in[0,T]} \lambda_j(t) \le \bar \lambda_j, \text{ and let }
    \bar \Lambda = \sum_{j \in \mathcal{J}} \bar \lambda_j.$$ Recall
that such multivariate counting process can be rewritten as a Poisson
random measure $N= \sum_{k\geq 1} \delta_{(T_k,J_k)}$ on
${\mathbb{R}}^+\times \mathcal{J}$ (see e.g. Sec. 2 of Chapter 6
in [@Cin11]), where $T_k$ is the $k$th jump time of
$\sum_{j\in \mathcal{J}} N^j$ and $J_k$ corresponds to the component of
the the vector which jumps. In particular,
$N^j_t = N([0,t]\times \{j\})$.

Once again the simulation of such process can be obtained from the
simulation of a (homogeneous) multivariate Poisson process of intensity
vector $(\bar{\lambda}_j)_{j \in \mathcal{J}}$, extended into a Poisson
measures by adding marks on ${\mathbb{R}}^+$. Thus, we introduce the
Marked Poisson measure
$\bar Q = \sum \delta_{(\bar T_\ell, \bar J_\ell, \bar \Theta_\ell)}$ on
${\mathbb{R}}^+  \times \mathcal{J}\times {\mathbb{R}}^+$, such that:

-   The jump times $(\bar T_\ell)$ of $\bar Q$ are the jump times of a
    Poisson measure of intensity $\bar \Lambda$.

-   The variables $(\bar J_\ell)$ are *i.i.d.* random variables on
    $\mathcal{J}$, with
    $\displaystyle p_j= {\mathbb{P}}(\bar J_1 = j)= \bar \lambda_j/\bar \Lambda$
    and representing the component of the vector which jumps.

-   The marks $(\bar \Theta_\ell)$ are independent variables with
    $\bar \Theta_\ell$ a uniform random variable on
    $[0,\bar \lambda_{{\bar J}_\ell}]$, $\forall \ell \geq 1$.

By Proposition [11](#PropMarkedPoisson){reference-type="ref"
reference="PropMarkedPoisson"} and
[10](#PropProjPoissonMeasure){reference-type="ref"
reference="PropProjPoissonMeasure"}, each measure
$\bar Q_j (\mathrm{d}t, \mathrm{d}\theta) = \bar Q(\mathrm{d}t, \{j\}, \mathrm{d}\theta) = \sum_{\ell \geq 1} \mathsf{1}_{\{\bar J_\ell=j \}} \delta_{(\bar T_\ell, \bar \Theta_\ell)}$
is a marked Poisson measure of intensity
$$\bar \mu_j ( \mathrm{d}t  ,\mathrm{d}\theta) = \bar{\Lambda}p_j \mathrm{d}t  \frac{\mathsf{1}_{\{\theta \leq \bar \lambda_{j}\}}(\theta)}{\bar \lambda_{j}} \mathrm{d}\theta = \mathrm{d}t  \mathsf{1}_{\{\theta \leq \bar \lambda_{j}\}}(\theta) \mathrm{d}\theta.$$
As a direct application of Proposition
[3](#PropThinning){reference-type="ref" reference="PropThinning"}, the
inhomogeneous multivariate Poisson process is obtained by restriction of
each measures $\bar Q_j$ to
$\Delta_j = \{ (t, \theta) \in [0,T] \times [0,\bar \lambda_j] ;\; \theta \leq \lambda_j(t) \}$
and projection:

::: {#PropThinningVector .prop}
**Proposition 4**. *The multivariate counting process
$(N^j)_{j \in \mathcal{J}}$, defined for all $j \in \mathcal{J}$ and
$t \in [0,T]$ by thinning and projection of $\bar Q$:
$$N^j_t  := \int_0^t \int_{ {\mathbb{R}}^+} \boldsymbol{1}_{\{\theta \leq \lambda_j(s)\}} \bar{Q}(\mathrm{d}s, \{j\}, \mathrm{d}\theta) = \sum_{\ell \geq 1} \mathsf{1}_{\{ \bar T_\ell \leq t \}}\mathsf{1}_{\{\bar J_\ell = j\}} \mathsf{1}_{\{\bar \Theta_\ell \leq  \lambda_j (\bar T_\ell ) \}},$$
is an inhomogeneous Poisson process of intensity vector
$(\lambda_j(t))_{j\in\mathcal{J}}$ on $[0,T]$.*
:::

Proposition [4](#PropThinningVector){reference-type="ref"
reference="PropThinningVector"} yields the following simulation
algorithm for multivariate Poisson processes:\

::: algorithm
Initialization $T_0 \longleftarrow 0$, $\bar T_0 \longleftarrow 0$
:::

::: {#RqAlternateThinningalgo .remark}
*Remark 3*. The acceptance/rejection algorithm
[\[AlgoThinning2\]](#AlgoThinning2){reference-type="ref"
reference="AlgoThinning2"} can be efficient when the functions
$\lambda_j$ are of different order, and thus bounded by different
$\bar \lambda_j$. However, it is important to note that the simulation
of the discrete random variables $(\bar J_\ell)$ can be costly (compared
to a uniform law) when $\mathcal{J}$ is large, for instance when an
individual is drawn from a large population. In this case, an
alternative is to choose the same bound $\bar \lambda_j= \bar \lambda$
for all $j \in \mathcal{J}$. Then the marks
$(\bar J_\ell, \bar \Theta_\ell)$ are i.i.d uniform variables on
$\mathcal{J}\times [0,\bar \lambda]$, faster to simulate.
:::

## Simulation algorithm {#sec::simulation_algo}

Let us now come back to the simulation of the IBM introduced in Section
[1](#section::IBM){reference-type="ref" reference="section::IBM"}. For
ease of notations, we assume that there are no event with Poisson
intensity ($\mathcal{P} =\emptyset$), so that all events that occur are
of type $e \in \mathcal{E} \cup \mathcal{E}_W$, with individual
intensity $\lambda_t^e(I,Z_t)$ depending on the population composition
$Z_t$ ($e \in \mathcal{E}_W$) or not ($e \in \mathcal{E}$), as defined
in
[\[IndividualIntensity\]](#IndividualIntensity){reference-type="eqref"
reference="IndividualIntensity"} and verifying either Assumption
[2](#AssumptionIntensity1){reference-type="ref"
reference="AssumptionIntensity1"} or
[3](#AssumptionIntensity2){reference-type="ref"
reference="AssumptionIntensity2"}. The global
intensity [\[eq:globalEvIntensity\]](#eq:globalEvIntensity){reference-type="eqref"
reference="eq:globalEvIntensity"} at time $t \in [0,T]$ is thus
$$\label{eq:defintensity}
    \Lambda_t(Z_t) = \sum_{e \in \mathcal{E}} \Big( \sum_{k=1}^{N_t} \lambda^e(t, I_k) \Big)
    + \sum_{e \in \mathcal{E}_W} \Big( \sum_{k=1}^{N_t} \sum_{j=1}^{N_t} W^e(t, I_k, I_j) \Big) \leq \bar \Lambda(N_t),$$
with
$\bar \Lambda(n) = \big(\sum_{e \in \mathcal{E}} \bar \lambda^e \big) n + \big( \sum_{e \in \mathcal{E}_W} \bar W^e \big) n^2$.

One of the main difficulty is that the intensity of events is not
deterministic as in the case of inhomogeneous Poisson processes, but a
function $\Lambda_t(Z_t)$ of the population state, bounded by a function
which also depends on the population size. However, the
algorithm [\[AlgoThinning2\]](#AlgoThinning2){reference-type="ref"
reference="AlgoThinning2"} can be adapted to simulate the IBM. The
construction is done by induction, by conditioning on the state of the
population $Z_{T_k}$ at the $k$th event time $T_k$ ($T_0 = 0$).

We first present the construction of the first event at time $T_1$.

#### First event simulation

Before the first event time (on $\{ t< T_1\}$), the population
composition is constant : $Z_t = Z_0 = \{ I_1, \dots, I_{N_0}\}$. For
each type of event $e$ and individual $I_k$, $k \in \{1,\dots N_0\}$, we
denote by $N^{k,e}$ the counting process of intensity
$\lambda_t^e (I_k,Z_t)$, counting the occurrences of the events of type
$e$ happening to the individual $I_k$. Then, the first event $T_1$ is
the first jump time of the multivariate counting vector
$(N^{(k,e)})_{ (k,e) \in \mathcal{J}_0}$, with
$\mathcal{J}_0 = \{1,\dots , N_0\}\times \big(\mathcal{E} \cup \mathcal{E}_W \big)$.

Since the population composition is constant before the first event
time, each counting process $N^{j}$ coincides on $[0,T_1[$ with an
inhomogeneous Poisson process, of intensity $\lambda_t^e (I_k,Z_0)$.
Thus (conditionally to $Z_0$), $T_1$ is also the first jump time of an
inhomogeneous multivariate Poisson process
$N^0 = (N^{0,j})_{j \in \mathcal{J}_0}$ of intensity function
$(\lambda_j)_{j\in \mathcal J_0}$, defined for all
$j = (k,e) \in \mathcal{J}_0$ by:
$$\lambda_j(t) = \lambda^e_t(I_k,Z_0) \le \bar \lambda^e_0 \quad \text{with} \quad \bar \lambda^e_0 = \bar \lambda^e \mathbf{1}_{e \in \mathcal{E}} + \bar W^e N_0 \mathbf{1}_{e \in \mathcal{E}_W},$$
by Assumptions [2](#AssumptionIntensity1){reference-type="ref"
reference="AssumptionIntensity1"}
and [3](#AssumptionIntensity2){reference-type="ref"
reference="AssumptionIntensity2"}. In particular, the jump times of
$N^0$ occur at the intensity
$$\Lambda(t) =\sum_{j \in \mathcal{J}_0} \lambda_j(t)  =\sum_{e \in \mathcal{E} \cup \mathcal{E}_W} \sum_{k=1}^{N_0}  \lambda^e_t(I_k,Z_0) \leq \bar \Lambda(N_0)=N_0 \sum_{e \in \mathcal{E} \cup \mathcal{E}_W} \bar \lambda^e_0.$$
By Proposition [4](#PropThinningVector){reference-type="ref"
reference="PropThinningVector"}, $N^0$ can be obtained by thinning of
the marked Poisson measure
$\bar Q^0 = \sum_{\ell \geq 1} \delta_{(\bar T_\ell , (\bar{K}_\ell, \bar E_\ell), \bar \Theta_\ell)}$
on ${\mathbb{R}}^+\times\mathcal{J}_0 \times {\mathbb{R}}^+$, with:

-   $(\bar T_\ell)_{\ell \in {\mathbb{N}}^*}$ the jump times of a
    Poisson process of rate $\bar \Lambda(N_0)$.

-   $(\bar{K}_\ell, \bar E_\ell)_{\ell \in {\mathbb{N}}^*}$ discrete
    *i.i.d.* random variables on
    $\mathcal{J}_0 = \{1,\dots , N_0\}\times \big(\mathcal{E} \cup \mathcal{E}_W \big)$,
    with $K_\ell$ representing the index of the chosen individual and
    $E_\ell$ the event type for the proposed event, such that:
    $$\mathbb{P}( \bar  K_1 = k, \bar E_1 = e) =  \frac{\bar \lambda^e_0 }{\bar \Lambda(N_0)}
    = \frac{1}{N_0} \frac{\bar \lambda^e_0 N_0}{\bar \Lambda(N_0)},$$
    i.e. $(\bar K_1, \bar E_1)$ are distributed as independent random
    variables where $\bar K_1 \sim \mathcal U(\{1,\dots, N_0\})$ and
    $\bar E_1$ such that $$p_e := \mathbb{P}( \bar E_1 = e)
        = \frac{\bar \lambda^e_0 N_0}{\bar \Lambda(N_0)}.$$

-   $(\bar \Theta_\ell)_{\ell \in {\mathbb{N}}^*}$ are independent
    uniform random variables, with
    $\bar \Theta_\ell \sim \mathcal{U}([0,\bar \lambda^{\bar E_\ell}]).$

Since the first event is the first jump of $N^0$, by Proposition
[4](#PropThinningVector){reference-type="ref"
reference="PropThinningVector"} and Algorithm
[\[AlgoThinning2\]](#AlgoThinning2){reference-type="ref"
reference="AlgoThinning2"}, the first event time $T_1$ is the first jump
time $\bar T_\ell$ of $\bar Q^0$ such that
$\bar \Theta_\ell \leq \lambda^{\bar E_\ell}_{\bar T_\ell}(I_{\bar K_\ell}, Z_0)$.

At $T_1 =\bar T_{\ell}$, the event $\bar E_\ell$ occurs to the
individual $I_{\bar K_\ell} = (\tau^b, \infty, x)$. For instance, if
$\bar E_\ell =d$, a death/exit event occurs, so that
$Z_{T_1} = Z_{0} + \delta_{(\tau^b, T_1, x)} - \delta_{I_{\bar K_\ell}}$
and $N_{T_1} = N_{0}$. If $\bar E_\ell =b$ or $en$, a birth or entry
event occurs, so that $N_{T_1} = N_{0} + 1$, and a new individual
$I_{N_0+1}$ is added to the population, chosen as described in Table
[1](#TableEvAction){reference-type="ref" reference="TableEvAction"}.
Finally, if $\bar E_\ell=s$, a swap event occurs, the population size
stays constant and $I_{\bar K_\ell}$ is replaced by an individual
$I_{\bar K_\ell}'$, chosen as described in Table
[1](#TableEvAction){reference-type="ref" reference="TableEvAction"}.

The steps for simulating the first event in the population can be
iterated in order to simulate the population. At the $k$th step, the
same procedure is repeated to simulate the $k$th event, starting from a
population $Z_{T_{k-1}}$ of size $N_{T_{k-1}}$.\

::: algorithm
Initialization $T_0 \longleftarrow 0$, $\bar T_0 \longleftarrow 0$
:::

The proof of Theorem
[5](#theorem:algoNoInteraction){reference-type="ref"
reference="theorem:algoNoInteraction"} is detailed in the Appendix
[8](#proof:algonointeraction){reference-type="ref"
reference="proof:algonointeraction"}.

::: {#theorem:algoNoInteraction .theo}
**Theorem 5**. *Algorithm
[\[algo:PopNointeraction2\]](#algo:PopNointeraction2){reference-type="ref"
reference="algo:PopNointeraction2"} are exact simulations of Equation
[\[SDE_pop\]](#SDE_pop){reference-type="eqref" reference="SDE_pop"}'s
solution.*
:::

::: {#remark::removaldeadIndividuals .remark}
*Remark 4*. The population $Z_{T_k}$ includes dead/exited individuals
before the event time $T_k$. Thus, $N_{T_k} > N_{T_k}^a$ is greater than
the number of alive individuals at time $T_k$. When a dead individual
$I_{\bar K_l}$ is drawn from the population during the
rejection/acceptance phase of the algorithm, the proposed event
$(\bar T_{\ell}, \bar E_\ell,  I_{\bar K_\ell})$ is automatically
rejected since the event intensity is
$\lambda^{\bar E_\ell}_{T_\ell}(I_{\bar K_\ell}, Z_{T_k}) = 0$ (nothing
can happen to a dead individual). This can slow down the algorithm,
especially when the proportion of dead/exited individuals in the
population increases. However, the computational cost of keeping
dead/exited individuals in the population is much lower than the cost of
removing an individual from the population at each death/exit event,
which is linear in the population size.

Actually, dead/exited individuals are regularly removed from the
population in the `IBMPopSim` algorithm, in order to optimize the
trade-off between having to many dead individuals and removing dead
individuals from the population too often. The frequency at which dead
individuals are "removed from the population" can be chosen by the user,
as an optional argument of the main function `popsim` (see details in
Section [3.4](#simulation){reference-type="ref"
reference="simulation"}).
:::

::: remark
*Remark 5*. In practice, the bounds $\bar \lambda^e$ and $\bar W^e$
should be chosen as sharp as possible. It is easy to see that
conditionally to
$\{\bar E_\ell = e, \bar T_\ell = t, \bar K_\ell = l \}$ the probability
of accepting the event is, depending if there are interactions,
$$\mathbb{P}\big(\bar \Theta_\ell \le \lambda^e_t(I_l,Z_{T_k})| \mathcal{F}_{T_k}\big)
        = \frac{\lambda^e(t, I_l)}{\bar \lambda^e} \mathbf{1}_{e \in \mathcal{E}} +  \frac{\sum_{j=1}^{N_{T_k}} W^e(t, I_l, I_j)}{\bar W^e N_{T_k}} \mathbf{1}_{e \in \mathcal{E}_W}.$$
The sharper the bounds $\bar \lambda^e$ and $\bar W^e$ are, the higher
is the acceptance rate.\
For even sharper bounds, an alternative is to define bounds
$\bar \lambda^e(I_l)$ and $\bar W^e(I_l)$ depending on the individuals'
characteristics. However, the algorithm is modified and the individual
$I_l$ is not chosen uniformly in the population anymore. Due to the
population size, this is way more costly than choosing uniform bounds,
as explained in Remark
[3](#RqAlternateThinningalgo){reference-type="ref"
reference="RqAlternateThinningalgo"}.
:::

## Simulation algorithm with randomization {#sec::simulation_algo_randomized}

Let $e \in \cal E_W$ be an event with interactions. In order to evaluate
the individual intensity
$\lambda^e_t (I,Z_t) = \sum_{j=1}^{N_t} W^e(t, I,I_j)$ one must compute
$W^e(t, I_l, I_j)$ for all individuals in the population. This step can
be computationally costly, especially for large populations. One way to
avoid this summation is to use randomization (see also [@FouMel04] in a
model without age). The randomization consists in replacing the
summation by an evaluation of the interaction function $W^e$ using an
individual $J$ drawn uniformly from the population.

More precisely, if $J \sim \mathcal{U}(\{1, \dots, N_{T_k}\})$ is
independent of $\bar \Theta_\ell$, we have $$\label{eq:randomizedevent}
\mathbb{P}\Big(\bar \Theta_\ell \le \sum_{j=1}^{N_{T_k}} W^e(t, I_l, I_j) | \mathcal{F}_{T_k} \Big)
    = \mathbb{P}\big(\bar \Theta_\ell \le N_{T_k} W^e(t, I_l, I_J) | \mathcal{F}_{T_k}\big).$$
Equivalently, we can write this probability as
$\mathbb{P}\big(\tilde \Theta_\ell \le W^e(t, I_l, I_J) \big)$ where
$\tilde \Theta_\ell  = \frac{\bar \Theta_\ell}{N_{T_k}}\sim \mathcal{U}([0, \bar W^e])$
is independent of $J \sim \mathcal{U}(\{1, \dots, N_{T_k}\})$.

::: remark
*Remark 6*. The efficiency of the randomization procedure increases with
the population homogeneity. If the function $W^e$ varies little
according to the individuals in the population, the randomization
approach is very efficient in practice, especially when the population
is large.
:::

We now present the main algorithm implemented in the `popsim` function
of the `IBMPopSim` package in the case where events arrive with
individual intensities, but also with interactions (using randomization)
and Poisson intensities. In this general case, $\bar \Lambda(n)$ is
defined by
[\[eq:defbarLambda\]](#eq:defbarLambda){reference-type="eqref"
reference="eq:defbarLambda"}.\

::: algorithm
[]{#algo::rzndomized label="algo::rzndomized"} Initialization
$T_0 \longleftarrow 0$, $\bar T_0 \longleftarrow 0$
:::

::: {#prop:algorandomized .prop}
**Proposition 6**. *The population processes $(Z_t)_{t\in [0,T]}$
simulated by the Algorithm
[\[algo:PopNointeraction2\]](#algo:PopNointeraction2){reference-type="ref"
reference="algo:PopNointeraction2"} and
[\[algo:Popinteraction2\]](#algo:Popinteraction2){reference-type="ref"
reference="algo:Popinteraction2"} have the same law.*
:::

::: proof
*Proof.* The only difference between Algorithm
[\[algo:PopNointeraction2\]](#algo:PopNointeraction2){reference-type="ref"
reference="algo:PopNointeraction2"} and
[\[algo:Popinteraction2\]](#algo:Popinteraction2){reference-type="ref"
reference="algo:Popinteraction2"} is in the acceptance/rejection step of
proposed events, in the presence of interactions. In Algorithm
[\[algo:Popinteraction2\]](#algo:Popinteraction2){reference-type="ref"
reference="algo:Popinteraction2"}, a proposed event
$(\bar T_\ell, \bar E_\ell, \bar K_\ell)$, with
$\bar E_l \in \mathcal{E}_W$ an event with interaction, is accepted as a
true event in the population if
$$\bar \Theta_\ell \le W^{\bar E_\ell}(\bar T_\ell, I_{\bar K_\ell}, I_{\bar J_\ell}), \text{ with } (\bar \Theta_\ell, \bar J_\ell) \sim  \mathcal{U}\big([0, \bar W^{\bar E_\ell}] \times \{1, \dots, N_{T_k}\} \big).$$
By [\[eq:randomizedevent\]](#eq:randomizedevent){reference-type="eqref"
reference="eq:randomizedevent"}, the probability of accepting this event
is the same than in Algorithm
[\[algo:PopNointeraction2\]](#algo:PopNointeraction2){reference-type="ref"
reference="algo:PopNointeraction2"}, which achieves the proof. ◻
:::

::: {#theorem:algoInteraction .cor}
**Corollary 7**. *Algorithm
[\[algo:Popinteraction2\]](#algo:Popinteraction2){reference-type="ref"
reference="algo:Popinteraction2"} is an exact simulation of Equation
[\[SDE_pop\]](#SDE_pop){reference-type="eqref" reference="SDE_pop"}'s
solution.*
:::

# Model creation and simulation with IBMPopSim {#sec::package}

The use of the `IBMPopSim` package is mainly done in two steps: a first
model creation followed by the simulation of the population evolution.
The creation of a model is itself based on two steps: the description of
the population $Z_t$, as introduced in
Section [1.2](#sec::population){reference-type="ref"
reference="sec::population"}, and the description of the events types,
along with their associated intensities, as detailed in
Sections [1.3](#sec::events){reference-type="ref"
reference="sec::events"}
and [1.4](#sec::event_intensity){reference-type="ref"
reference="sec::event_intensity"}. A model is compiled by calling the
`mk_model` function, which internally uses a template mechanism to
generate automatically the source code describing the model, which is
subsequently compiled using the `Rcpp` package to produce the object
code.

After the compilation of the model, the simulations are launched by
calling the `popsim` function. This function depends on the previously
compiled model and simulates a random trajectory of the population
evolution based on an initial population and on parameter values, which
can change from a call to another.

In this section, we take a closer look at each component of a model in
`IBMPopSim`. We also refer to the [IBMPopSim
website](https://daphnegiorgi.github.io/IBMPopSim/) and to the
`vignettes` of the package for more details on the package and various
examples of model creation.

## Population {#population}

A population $Z$ is represented by an object of class `population`
containing a data frame where each row corresponds to an individual
$I=(\tau^b, \tau^d, x)$, and which has at least two columns, `birth` and
`death`, corresponding to the birth date $\tau^b$ and death/exit date
$\tau^d$ ($\tau^d$ is set to `NA` for alive individuals). The data frame
can contain more than two columns if individuals are described by
additional characteristics $x= (x_1,\dots x_n)$.

#### Entry and exit events

If entry events can occur in the population, the population shall
contain a characteristic named `entry`. This can be done by setting the
flag `entry=TRUE` in the `population` function, or by calling the
`add_characteristic` function on an existing population. During the
simulation, the date at which an individual enters the population is
automatically recorded in the variable `I.entry`.\
If exit events can occur, the population shall contain a characteristic
named `out`. This can be done by setting the flag `out=TRUE` in the
`population` function, or by calling the `add_characteristic` function.
When an individual `I` exits the population during the simulation,
`I.out` is set to `TRUE` and its exit time is recorded as a "death"
date.

In the example below, individuals are described by their birth and death
dates, as well a Boolean characteristics called male, and the `entry`
characteristic. For instance, the first individual is a female whose age
at $t_0=0$ is $107$ and who was originally in the population.

#### Individual

In the `C++` model which is automatically generated and compiled, an
individual `I` is an object of an internal class containing some
attributes (`birth_date`, `death_date` and the characteristics, here
`male`), and some methods including:

-   `I.age(t)`: a `const` method returning the age of an individual `I`
    at time `t`,

-   `I.set_age(a, t)`: a method to set the age `a` at time `t` of an
    individual `I` (set `birth_date` at `t-a`),

-   `I.is_dead(t)`: a `const` method returning `true` if the individual
    `I` is dead at time `t`.

::: remark
*Remark 7* (Characteristics type). A characteristic $x_i$ must be of
atomic type: logical, integer, double or character. The function
`get_characteristic` allows to easily get characteristics names and
their types from a population data frame. We draw the attention to the
fact that some names for characteristics are forbidden, or reserved to
specific cases : this is the case for `birth, death, entry, out, id`.
:::

## Events {#sec::package_events}

The most important step of the model creation is the events creation.
The call to the function creating an event is of form

where `CLASS` is replaced by the class of the event intensity, described
in Section [1.4](#sec::event_intensity){reference-type="ref"
reference="sec::event_intensity"}, and `type` corresponds to the event
type, described in Section [1.3](#sec::events){reference-type="ref"
reference="sec::events"}.
Tables [2](#tab::intensity_classes){reference-type="ref"
reference="tab::intensity_classes"}
and [3](#tab::event_types){reference-type="ref"
reference="tab::event_types"} summarize the different possible choices
for intensity classes and types of event. The optional argument `name`
gives a name to the event. If not specified, the name of the event is
its type, for instance `death`. However, a name must be specified if the
model is composed of several events with the same type (for instance
when there are multiple death events corresponding to different causes
of death). The other arguments depend on the intensity class and on the
event type.

::: {#tbl-panel layout-ncol=2}

Intensity class         Set               `CLASS`
----------------------- ----------------- -------------------------
Individual              $\mathcal{E}$     `individual`
Interaction             $\mathcal{E}_W$   `interaction`
Poisson                 $\mathcal{P}$     `poisson`
Inhomogeneous Poisson   $\mathcal{P}$     `inhomogeneous_poisson`

: Intensity Classes {#tbl-intensity-classes}

Event type   `TYPE`
------------ ---------
Birth        `birth`
Death        `death`
Entry        `entry`
Exit         `exit`
Swap         `swap`

: Event Types {#tbl-event-types}

:::

The intensity function and the kernel of an event are defined through
arguments of the function `mk_event_CLASS`. These arguments are strings
composed of few lines of code. Since the model is compiled using `Rcpp`,
the code should be written in `C++`. However, thanks to the
functions/variables of the package, even the non-experienced `C++` user
can define a model quite easily. To facilitate the implementation, the
user can also define a list of **model parameters**, which can be used
in the event and intensity definitions. These parameters are stored in a
named list and can be of various types: atomic type, numeric vector or
matrix, predefined function of one variable (`stepfun`, `linfun`,
`gompertz`, `weibull`, `piecewise_x`), piecewise functions of two
variables (`piecewise_xy`). We refer to the `vignette(IBMPopSim_cpp)`
for more details on parameters types and basic `C++` tools. Another
advantage of the model parameters is that their value can be modified
from a simulation to another without changing the model.

### Intensities

In `IBMPopSim`, the intensity of an event can belong to three classes
(see Section [1.4](#sec::event_intensity){reference-type="ref"
reference="sec::event_intensity"}): individual intensities without
interaction between individuals, corresponding to events
$e\in\mathcal{E}$, individual intensities with interaction,
corresponding to events $e\in\mathcal{E}_W$, and Poisson intensities
(homogeneous and inhomogeneous), corresponding to events
$e\in\mathcal{P}$.

#### Event creation with individual intensity

An event $e\in \mathcal{E}$ (see
 [\[eq:intensityNointeraction\]](#eq:intensityNointeraction){reference-type="eqref"
reference="eq:intensityNointeraction"}) has an intensity of the form
$\lambda^e(t, I)$ which depends only on the individual `I` and time.
Events with such intensity are created using the function

The `intensity_code` argument is a character string containing few lines
of `C++` code describing the intensity function $\lambda^e(t, I)$. The
intensity value has to be stored in a variable called `result` and the
available variables for the intensity code are given in
Table [4](#tab::intensity-variables){reference-type="ref"
reference="tab::intensity-variables"}.

For instance, the intensity code below corresponds to an individual
death intensity $\lambda^d(t,I)$ equal to
$d_1(a(I,t)) = \alpha_1 \exp(\beta_1 a(I,t))$ for males and
$d_2(a(I,t)) = \alpha_2 \exp(\beta_2 a(I,t))$ for females, where
$a(I,t)=t-\tau^b$ is the age of the individual $I=(\tau^b, \tau^d,x)$ at
time $t$. In this case, the intensity function depends on the
individuals' age, gender, and on the model parameters
$\alpha = (\alpha_1, \alpha_2)$ and $\beta = (\beta_1, \beta_2)$.

#### Event creation with interaction intensity

An event $e\in \mathcal{E}_W$ is an event which occurs to an individual
at a frequency which is the result of interactions with other members of
the population (see
Equation [\[eq:intensityInteraction\]](#eq:intensityInteraction){reference-type="eqref"
reference="eq:intensityInteraction"}), and which can be written as
$\lambda^e_t(I, Z_t)=\sum_{J\in Z_t} W^e(t, I, J)$ where $W^e(t, I, J)$
is the intensity of the interaction between individual $I$ and
individual $J$.

An event $e\in \mathcal{E}_W$ with such intensity is created by calling
the function

The `interaction_code` argument contains few lines of `C++` code
describing the interaction function $W^e(t, I, J)$. The interaction
function value has to be stored in a variable called `result` and the
available variables for the intensity code are given in
Table [4](#tab::intensity-variables){reference-type="ref"
reference="tab::intensity-variables"}. For example, if we set

the death intensity of an individual `I` is the result of the
competition between individuals, depending on a characteristic named
`size`, as defined in
[\[ex:interaction\]](#ex:interaction){reference-type="eqref"
reference="ex:interaction"}.

The argument `interaction_type`, set by default at `random`, is the
algorithm choice for simulating the model. When `interaction_type=full`,
the simulation follows
Algorithm [\[algo:PopNointeraction2\]](#algo:PopNointeraction2){reference-type="ref"
reference="algo:PopNointeraction2"}, while when
`interaction_type=random` it follows
Algorithm [\[algo:Popinteraction2\]](#algo:Popinteraction2){reference-type="ref"
reference="algo:Popinteraction2"}. In most cases, the `random` algorithm
is much faster than the `full` algorithm, as we illustrate for instance
in Section [5](#section:ExempleInteraction){reference-type="ref"
reference="section:ExempleInteraction"}, where we observe the gain of a
factor of 40 between the two algorithms, on a set of standard
parameters. This allows in particular to explore parameter sets that
give larger population sizes, without reaching computation times that
explode.

::: {#tab::intensity-variables}
  Variable           Description
  ------------------ -------------------------------------------------------------
  Variable           Description
  `I`                Current individual
  `J`                Another individual in the population (only for interaction)
  `t`                Current time
  Model parameters   Depends on the model

  : `C++` variables available for intensity code
:::

#### Events creation with Poisson and Inhomogeneous Poisson intensity

For events $e\in\mathcal{P}$ with an intensity $\mu^e(t)$ which does not
depend on the population, the event intensity is of class
`inhomogeneous_poisson` or `poisson` depending on whether or not the
intensity depends on time (in the second case the intensity is
constant).

For Poisson (constant) intensities the events are created with the
function

The following example creates a death event, where individuals die at a
constant intensity `lambda` (which has to be in the list of model
parameters):

When the intensity $(\mu^e(t))$ depends on time, the event can be
created similarly by using the function

### Event kernel code {#sub_par::event-kernel-code}

When an event occurs, the events kernels $k^e$ specify how the event
modifies the population. The events kernels are defined in the
`kernel_code` parameter of the
`mk_event_CLASS(type = "TYPE", name ="NAME", ...)` function. The
`kernel_code` is `NULL` by default and doesn't have to be specified for
death, exit events and birth events, but mandatory for entry and swap
events. Recall that the `kernel_code` argument is a string composed of a
few lines of `C++` code, characterizing the individual characteristics
following the event.
@tbl-events-variables summarizes the list of available variables that can be used in the `kernel_code`.

-   **Death/Exit event** If the user defines a death event, the death
    date of the current individual `I` is set automatically to the
    current time `t`. Similarly, when an individual `I` exits the
    population,`I.out` is set automatically to `TRUE` and his exit time
    is recorded as a \"death\" date. For these events types, the
    `kernel_code` doesn't have to be specified by the user.

-   **Birth event** The default generated event kernel is that an
    individual `I` gives birth to a new individual `newI` of age 0 at
    the current time `t`, with same characteristics than the parent `I`.
    If no kernel is specified, the default generated `C++` code for a
    birth event is:

    The user can modify the birth kernel, by specify the argument
    `kernel_code` of `mk_event_CLASS`. In this case, the generated code
    is

    where `_KERNEL_CODE_` is replaced by the content of the
    `kernel_code` argument.

-   **Entry event** When an individual `I` enters the population,
    `I.entry` is set automatically as the date at which the individual
    enters the population. When an entry occurs the individual entering
    the population is not of age $0$. In this case, the user must
    specify the `kernel_code` argument indicating how the age and
    characteristics of the new individual are chosen. For instance, the
    code below creates an event of type `entry`, named `ev_example`,
    where individuals enter the population at a Poisson constant
    intensity. When an individual `newI` enters the population at time
    `t`, his age is chosen as a normally distributed random variable,
    with mean 20 and variance 4.

-   **Swap event** The user must specify the `kernel_code` argument
    indicating how the characteristics of an individual are modified
    following a swap.

::: {#tbl-events-variables} 
  Variable                                     Description
  -------------------------------------------- -----------------------------
  Variable                                     Description
  `I`                                          Current individual
  `t`                                          Current time
  `pop`                                        Current population (vector)
  `newI`                                       
  Available only for birth and entry events.   
  Model parameters                             Depends on the model

  : `C++` variables available for events kernel code
:::

::: remark
*Remark 8*. When there are several events of the same type, the user can
identify which events generated a particular event by adding a
characteristic to the population recording the event name/id when it
occurs. See e.g. `vignette(IBMPopSim_human_pop)` for an example with
different causes of death.
:::

## Model creation {#Modelcreation}

[]{#sec::model_creation label="sec::model_creation"}

Once the population, the events, and model parameters are defined, the
IBM model is created using the function `mk_model`.

During this step which can take a few seconds, the model is created and
compiled using the `Rcpp` package. The model structure in `IBMPopSim` is
that the model depends only on the population characteristics' and
parameters names and types, rather than their values. This means that
once the model has been created, various simulations can be done with
different initial populations and different parameters values.

#### Example

Here is an example of model with a population structured by age and
gender, with birth and death events. The death intensity of an
individual of age $a$ is $d(a) = \alpha \exp(\beta a),$ and females
between 15 and 40 can give birth with birth intensity
$b(a) = \bar \lambda^b \mathbf{1}_{[15,40]}.$ The newborn is a male with
probability $p_{male}$.

## Simulation {#simulation}

The simulation of the IBM is based on the algorithms presented in
Sections [2.2](#sec::simulation_algo){reference-type="ref"
reference="sec::simulation_algo"}
and [2.3](#sec::simulation_algo_randomized){reference-type="ref"
reference="sec::simulation_algo_randomized"}. The user has first to
specify bounds for the intensity or interaction functions of each event
type. The random evolution of the population can then be simulated over
a period of time $[0,T]$ by calling the function `popsim`

#### Events bounds

Since the IBM simulation algorithm is based on an acceptance-rejection
method for simulating random times, the user has to specify bounds for
the intensity (or interaction) functions of each event (see Assumptions
[2](#AssumptionIntensity1){reference-type="ref"
reference="AssumptionIntensity1"} and
[3](#AssumptionIntensity2){reference-type="ref"
reference="AssumptionIntensity2"}). These bounds should be stored in a
named vector, where for event $e$, the name corresponding to the event
bound $\bar{\mu}^e$, $\bar{\lambda}^e$ or $\bar{W}^e$ is the event
`name` defined during the event creation step.

In the model example built in the previous section the intensity bound
for birth events is $\bar\lambda_b$. Since the death intensity function
is not bounded, the user will have to specify a maximum age $a_{max}$ in
`popsim` (all individuals above $a_{max}$ die automatically). Then, the
bound for death events is $\bar \lambda_d = \alpha\exp(\beta a_{max}).$
In the example, the death event has been named `my_death_event`. No name
has been specified for the birth event which thus has the default name
`birth`. Then,

Once the model and events bounds have been defined, a random trajectory
of the population can be simulated by calling

#### Optional parameters

If there are no events with intensity of class `interaction`, then the
simulation can be parallelized easily by setting the optional parameter
`multithreading` (`FALSE` by default) to `TRUE`. By default, the number
of threads is the number of concurrent threads supported by the
available hardware implementation. The number of threads can be set
manually with the optional argument `num_threads`. By default, when the
proportion of dead individuals in the population exceeds $10\%$, dead
individuals are removed from the current population used in the
algorithm (see Remark
[4](#remark::removaldeadIndividuals){reference-type="ref"
reference="remark::removaldeadIndividuals"}). The user can modify this
ratio using the optional argument `clean_ratio`, or by removing dead
individuals from the population with a certain frequency, given by the
`clean_step` argument. Finally, the user can also define the seed of the
random number generator stored in the argument `seed`.

#### Outputs and treatment of swap events

The output of the `popsim` function is a list containing three elements:
a data frame `population` containing the output population $Z_T$ (or a
list of populations $(Z_{t_1}, \dots Z_{t_n})$ if `time` is a vector of
times), a numeric vector `logs` of variables related to the simulation
algorithm (including the simulation time and number of proposed/accepted
events), and the list `arguments` of the simulation inputs, including
the initial population, parameters and event bounds used for the
simulation.

When there are no swap events (individuals don't change of
characteristics), the evolution of the population over the period
$[0,T]$ is recorded in a single data frame `sim_out$population` where
each line contains the information of an individual who lived in the
population over the period $[0,T]$ (see Remark
[1](#remark::popfinale){reference-type="ref"
reference="remark::popfinale"}).

When there are swap events (individuals can change of characteristics),
recording the dates of swap events and changes of characteristics
following each swap event and for each individual in the population is a
memory intensive and computationally costly process. To maintain
efficient simulations in the presence of swap events, the argument
`time` of `popsim` can be defined as a vector of dates
$(t_0,\dots, t_n)$. In this case, `popsim` returns in the object
`population` a list of $n$ populations representing the population at
time $t_1,\dots t_n$, simulated from the initial time $t_0$. For
$i=1\dots n$, the $i$th data frame is the population $Z_{t_i}$, i.e.
individuals who lived in the population during the period $[t_0,t_i]$,
with their characteristics at time $t_i$.

It is also possible to isolate the individuals' life course, by adding
an `id` column to the population, which can be done by setting `id=TRUE`
in the population construction, or by calling the `add_characteristic`
function to an existing population, in order to identify each individual
with a unique integer.

Base functions to study the simulation outputs are provided in the
package. For instance, the population age pyramid can computed at a give
time, as well as death and exposure tables. Several illustrations of the
outputs functions are given in the example Sections
[4](#SectionInsurancePortofio){reference-type="ref"
reference="SectionInsurancePortofio"} and
[5](#section:ExempleInteraction){reference-type="ref"
reference="section:ExempleInteraction"}.

# Insurance portfolio {#SectionInsurancePortofio}

This section provides an example of how to use the `IBMPopSim` package
to simulate a heterogeneous life insurance portfolio (see
also `vignette(IBMPopSim_insurance_portfolio)`).

We consider an insurance portfolio consisting of male policyholders, of
age greater than 65. These policyholders are characterized by their age,
assumed to be less than $a_{max} = 110$, and risk class
$x\in \mathcal X =\{1,2\}$.

**Entries in the portfolio** New policyholders enter the population at a
constant Poisson rate $\mu^{en}=\lambda$, which means that on average,
$\lambda$ individuals enter the portfolio each year. A new individual
enters the population at an age a that is uniformly distributed between
65 and 70, and is in risk class 1 with probability $p$.

**Death events** A baseline age and time specific death rate is first
calibrated on "England and Wales (EW)" males mortality historic
data[^1], and projected for 30 years using the Lee-Carter model with the
package `StMoMo` (see [@stmomo]). The forecasted baseline death
intensity is denoted by $d(t,a)$, defined by:
$$\label{eq:insurance-baseline}
d(t,a) = \sum_{k=0}^{29}\mathbf{1}_{\{k\leq t < k+1\}} d_k(a), \quad \forall \; t\in [0,30] \text{ and } a \in [65, a_{max}],$$
with $d_k(a)$ the point estimate of the forecasted mortality rate for
age $a$ and year $k$.\
Individuals in risk class 1 are assumed to have mortality rates that are
20% higher than the baseline mortality (for instance, the risk class
could refer to smokers), while individuals in risk class 2 are assumed
to have mortality rates that are 20% lower than the baseline (non
smokers). The death intensity of an individual $I= (\tau_b, \infty, x)$,
of age $a(I,t) = t - \tau_b$ at time $t$ and in risk class
$x \in \{1, 2\}$ is thus the function $$\label{eq:insurance-deathrates}
\lambda^d(t,I) = \alpha_x d(t,a(I,t)), \quad \alpha_1 = 1.2, \quad \alpha_2 = 0.8.$$
In particular, the death intensity verifies Assumption
[2](#AssumptionIntensity1){reference-type="ref"
reference="AssumptionIntensity1"} since:
$$\label{eq:insurance-bound-deathrates}
\lambda^d(t,I) \leq \bar d : = \alpha_1 \sup_{t \in [0,30]} d(t,a_{max}).$$

**Exits from the portfolio** Individuals exit the portfolio at a
constant (individual) rate $\lambda^{ex}(t,I) = \mu^{i}$ only depending
on their risk class $i\in \{1,2\}$.

## Population {#insurance-population}

We start with an initial population of $30\,000$ males of age 65,
distributed uniformly in each risk class. The population data frame has
thus the two (mandatory) columns `birth` (here the initial time is
$t_0=0$) and `death` (`NA` if alive), and an additional column
`risk_cls` corresponding to the policyholders risk class. Since there
are entry and exit events, the `entry` and `out` flags of the population
constructor are set to `TRUE`.

## Events {#insurance-events}

#### Entry events

The age of the new individual is determined by the `kernel_code`
argument in the `mk_event_poisson` function.

Note that the variables `newI` and `t`, as well as the function
`CUnif()`, are implicitly defined and usable in the `kernel_code`. The
field `risk_cls` comes from the names of characteristics of individuals
in the population. The names `lambda` and `p` are parameter names that
will be specified in the `R` named list `params`.

Here we use a constant $\lambda$ as the event intensity, but we could
also use a rate $\lambda(t)$ that depends on time, using the function
`mk_event_poisson_inhomogeneous`.

#### Death and exit events

The baseline death intensity defined
in [\[eq:insurance-baseline\]](#eq:insurance-baseline){reference-type="eqref"
reference="eq:insurance-baseline"} and obtained with the package
`StMoMo` is stored in the variable `death_male`.

The death and exit intensities are of class `individual` (see Table
[2](#tab::intensity_classes){reference-type="ref"
reference="tab::intensity_classes"} ). Hence, the death and exit events
are created with the `mk_event_individual` function.

## Model creation and simulation {#insurance-simulation}

The model is created from all the previously defined building blocks, by
calling the `mk_model`.

Once the model is compiled, it can be used with different parameters and
run simulations for various scenarios. Similarly, the initial population
(here `pop_df`) can be modified without rerunning the `mk_model`
function. The bounds for entry events is simply the intensity $\lambda$.
For death events, the bound is given by $\bar{d}$ defined
[\[eq:insurance-bound-deathrates\]](#eq:insurance-bound-deathrates){reference-type="eqref"
reference="eq:insurance-bound-deathrates"}, which is stored in the
`death_max` variable.

## Outputs

The data frame `sim_out$population` consists of all individuals present
in the portfolio during the period of $[0, 30]$, including the
individuals in the initial population and those who entered the
portfolio. Each row represents an individual, with their date of birth,
date of death (`NA` if still alive at the end of the simulation), risk
class, and characteristics `entry` and `out`. Recall that if an
individual enters the population at time $t$, his `entry` characteristic
is automatically set up to be equal to $t$. The characteristics `out` is
set to `TRUE` for individuals who left the portfolio due to an exit
event.

In this example, the simulation time over 30 years, starting from an
initial population of 30 000 individuals is of $2\times 10^{-4}$
seconds, for an acceptance rate of proposed event of approximately 25%.
At the end of the simulation, the number of alive individuals is
approximately 430 000.

Initially in the portfolio (at $t=0$), there is the same number of 65
years old policyholders in each risk class. However, policyholders in
the risk class 2 with lower mortality rates leave the portfolio at
higher rate than policyholders in the risk class 1 : $\mu^2 > \mu^1$.
Therefore, the heterogeneous portfolio composition changes with time,
including more and more individuals in risk class 1 with higher
mortality rates, but with variations across age classes. To illustrate
the composition of the total population at the end of the simulation
($t=30$), we present in Figure [2](#fig:insur){reference-type="ref"
reference="fig:insur"}(a) the age pyramid of the final composition of
the portfolio obtained with the `age_pyramid` and `plot` function of the
`pyramid` class.

`IBMPopSim` also allows the fast computation of exact life tables from
truncated and censored individual data (due to entry and exit events),
using the functions `death_table` and `exposure_table`. These function
are particularly efficient, since the computations are made using the
`Rccp` library.

In Figure [2](#fig:insur){reference-type="ref"
reference="fig:insur"}(b), we illustrate the central death rates in the
simulated portfolio at final time. Due to the mortality differential
between risk class 1 and 2, one would expect to observe more individuals
in risk class 2 at higher ages. However, due to exit events, a higher
proportion of individuals in risk class 1 exit the portfolio over time,
resulting in a greater proportion of individuals in risk class 1 at
higher ages than what would be expected in the absence of exit events.
Consequently, the mortality rates in the portfolio are more aligned with
those of risk class 1 at higher ages. This is a simple example of how
composition changes in the portfolio can impact aggregated mortality
rates and potentially compensate or reduce an overall mortality
reduction (see also [@KAAKAI201916]).

<figure id="fig:insur">
<figure>
<img src="images/insur_pyr_group-1.png" />
<figcaption>Portfolio age pyramid at <span
class="math inline"><em>t</em> = 30</span>.</figcaption>
</figure>
<figure>
<img src="images/insur_mortality_rates-1.png" />
<figcaption>Portfolio central death rates at <span
class="math inline"><em>t</em> = 30</span> (black).</figcaption>
</figure>
<figcaption>Information obtained from a simulation of the portfolio
evolving over 30 years with individuals in risk class 1 (blue) and 2
(red).</figcaption>
</figure>

# Population with genetically variable traits {#section:ExempleInteraction}

This section provides an example of how to use the `IBMPopSim` package
to simulate an age-structured population with interactions, based on the
model proposed in Example 1 of [@FerTra09].

In this model, individuals are characterized by their body size at birth
$x_0 \in [0,4]$ and by their physical age $a \in [0,2]$. The body size
of an individual $I=(\tau^b,\infty, x_0)$ at time $t$ is a linear
function of its age $a(I,t) = t-\tau^b$: $$x(t)= x_0 + ga(I,t),$$ where
$g$ is a constant growth rate assumed to be identical for all
individuals.

**Birth events** The birth intensity of each individual
$I=(\tau^b, \infty, x_0)$ depends on a parameter $\alpha > 0$ and on its
initial size, as given by the equation
$$\label{eq::interaction_birth_intensity}
\lambda^b(t,I) = \alpha (4 - x_0) \leq \bar \lambda^b = 4\alpha.$$ Thus,
smaller individuals have a higher birth intensity. When a birth occurs,
the new individual inherit the same birth size $x_0$ as its parent with
high probability $1-p$, or a mutation can occur with probability $p$,
resulting in a birth size given by
$$\label{eq::interaction_birth_kernel}
    x_0' = \min(\max(0, x_0 + G), 4),$$ where $G$ is a Gaussian random
variable with mean 0 and variance $\sigma^2$.

**Death events** Due to competition between individuals, the death
intensity of an individual depends on the size of other individuals in
the population. Bigger individuals have a better chance of survival. If
an individual $I= (\tau^b, \infty, x_0)$ of size $x(t)= x_0 +ga(I,t)$
encounters an individual $J= (\tau^{b}_J, \infty, x_0')$ of size
$x'(t) = x_0'+ ga(J,t)$, then it can die with the intensity
$$W(t, I,J) = U(x(t),x'(t)),$$ where the interaction function $U$ is
defined by $$\label{eq::interaction_death_intensity}
    U(x,y) = \beta \left(1- \frac{1}{1+ c\exp(-4(x-y))}\right) \leq \bar W = \beta.$$
The death intensity of an individual $I$ at time $t$ and in a population
$Z$ is the result of interactions with all individuals in the
population, including itself, and is given by
$$\lambda^d_t(I,Z) = \sum_{J = (\tau^b,\infty, x_0') \in Z}  W (x_0 + g a(I,t), x_0' + g a(J,t)),$$

## Population {#population-1}

We use an initial population of 900 living individuals, all of whom have
the same size and ages uniformly distributed between 0 and 2 years.

## Events {#events}

#### Birth events

The parameters involved in a birth event are the probability of mutation
$p$, the variance of the Gaussian random variable and the coefficient
$\alpha$ of the intensity.

The birth
intensity [\[eq::interaction_birth_intensity\]](#eq::interaction_birth_intensity){reference-type="eqref"
reference="eq::interaction_birth_intensity"} is of class `individual`.
Hence, the event is created by calling the `mk_event_individual`
function. The size of the new individual is given in the kernel
following [\[eq::interaction_birth_kernel\]](#eq::interaction_birth_kernel){reference-type="eqref"
reference="eq::interaction_birth_kernel"}.

#### Death events

The death
intensity [\[eq::interaction_death_intensity\]](#eq::interaction_death_intensity){reference-type="eqref"
reference="eq::interaction_death_intensity"} is of class `interaction`.
Hence, the event is created by calling the `mk_event_interaction`
function. The parameters used for this event are the growth rate $g$,
the amplitude of the interaction function $\beta$, and the strength of
competition $c$.

## Model creation and simulation {#interaction-simulation}

The model is created using the `mk_model` function.

The simulation of one scenario can then be launched with the call of the
`popsim` function, after computing the events bounds
$\bar \lambda^b=4 \alpha$ and $\bar W= \beta$.

Based on the results of a simulation, we can reproduce the numerical
results of [@FerTra09]. In
Figure [3](#fig:interaction){reference-type="ref"
reference="fig:interaction"}(a), we draw a line for each individual in
the population to represent their birth size during their lifetime.

The randomized version of
Algorithm [\[algo:Popinteraction2\]](#algo:Popinteraction2){reference-type="ref"
reference="algo:Popinteraction2"} allows for much faster computation
times than
Algorithm [\[algo:PopNointeraction2\]](#algo:PopNointeraction2){reference-type="ref"
reference="algo:PopNointeraction2"}. This is illustrated in
Figure [3](#fig:interaction){reference-type="ref"
reference="fig:interaction"} (b), where we progressively decrease the
value of the mortality rate parameter $\beta$ and increase the birth
rate parameter $\alpha$. Starting with the values provided
in [@FerTra09], $\alpha=1$ and $\beta=2/300$, resulting in a stationary
population size of approximately $N=360$ individuals for a sample of 50
simulations, we can easily increase the stationary population size to
approximately $N=2600$ individuals with $\alpha=2$ and $\beta=1/300$
[^2]. In the log-scaled figure, we can observe the trend of computation
time as a function of the population size $N$, which is linear for the
randomized algorithm and quadratic for the full one (Algorithm
[\[algo:PopNointeraction2\]](#algo:PopNointeraction2){reference-type="ref"
reference="algo:PopNointeraction2"}). We can also see that the
randomized version of the algorithm is between 17 to 100 times faster
than the full one in this example, taking only 2 seconds in average for
the randomized version versus 211 seconds for Algorithm
[\[algo:PopNointeraction2\]](#algo:PopNointeraction2){reference-type="ref"
reference="algo:PopNointeraction2"} for the biggest population size
($N=2600$) and $T=500$.

<figure id="fig:interaction">
<figure>
<img src="images/inter_output2-1.png" />
<figcaption>Birth size during life time.</figcaption>
</figure>
<figure>
<img src="images/time_pop_size.png" />
<figcaption>Full vs random algorithm.</figcaption>
</figure>
<figcaption>Reproducing the example presented in <span class="citation"
data-cites="FerTra09"></span> and increasing the population size to
observe the difference in computing time between the randomized and full
algorithm.</figcaption>
</figure>

# Appendix {#appendix .unnumbered}

# Recall on Poisson random measures {#section::preliminaries}

We recall below some useful properties of Poisson random measures,
mainly following Chapter 6 of [@Cin11]. We also refer to [@Kal17] for a
more comprehensive presentation of random counting measures.

::: {#DefPoissonRandomMeasure .definition}
**Definition 8** (Poisson Random measures). *Let $\mu$ be a
$\sigma$-finite diffuse measure on a Borel subspace $(E,\mathcal{E})$ of
$({\mathbb{R}}^d, \mathcal{B}({\mathbb{R}}^d))$.. A random counting
measure $Q= \sum_{k\geq 1} \delta_{X_k}$ is a Poisson (counting) random
measure of *mean measure* $\mu$ if*

1.  *$\forall A \in \mathcal{E}$, $Q(A)$ is a Poisson random variable
    with ${\mathbb{E}}[Q(A)]= \mu(A)$.*

2.  *For all disjoints subsets $A_1, \dots , A_n \in \mathcal{E}$,
    $Q(A_1), \dots, Q(A_n)$ are independent Poisson random variables.*
:::

Let us briefly recall here some simple but useful operations on Poisson
measures. In the following, $Q$ is a Poisson measure of mean measure
$\mu$, unless stated otherwise.

::: {#PropRestrictionPoissonMeasure .prop}
**Proposition 9** (Restricted Poisson measure). *If $B \in \mathcal{E}$,
then, the restriction of $Q$ to $B$ defined by
$$Q^B = \boldsymbol{1}_B Q = \sum_{k \ge 1} \mathbf{1}_{B}(X_k) \delta_{X_k}$$
is also a Poisson random measure, of mean measure
$\mu^B = \mu(\cdot \cap B)$.*
:::

::: {#PropProjPoissonMeasure .prop}
**Proposition 10** (Projection of Poisson measure). *If
$E = F_1 \times F_2$ is a product space, then the projection
$$Q_1(\mathrm{d}x) = \int_{F_2} Q(\mathrm{d}x , \mathrm{d}y)$$ is a
Poisson random measure of mean measure
$\mu_1 (\mathrm{d}x ) = \int_{F_2} \mu(\mathrm{d}x, \mathrm{d}y)$.*
:::

#### Link with Poisson processes

Let $Q= \sum_{k\geq 1} \delta_{T_k}$ a Poisson random measure on
$E={\mathbb{R}}^+$ with mean measure
$\mu(\mathrm{d}t) = \Lambda (t) \mathrm{d}t$ absolutely continuous with
respect to the Lebesgue measure
($\mu(A) = \int_A \Lambda(t) \mathrm{d}t$). The counting process
$(N_t)_{t \ge 0}$ defined by $$\label{eq::inhomogeneous_PP}
    N_t = Q([0,t]) = \sum_{k\geq 1} \boldsymbol{1}_{\{T_k \leq t\}}, \quad \forall \; t\geq 0,$$
is an inhomogeneous Poisson process with intensity function (or rate)
$t \mapsto \Lambda(t)$. In particular, when $\Lambda(t) \equiv c$ is a
constant, $N$ is a homogeneous Poisson process with rate $c$. Assuming
that the atoms are ordered $T_1< T_2< \dots$, we recall that the
sequence $(T_{k+1}-T_k)_{k\geq 1}$ is a sequence of *i.i.d.* exponential
variables of parameter $c$.

#### Marked Poisson measures on $E = {\mathbb{R}}^+ \times F$

We are interested in the particular case when $E$ is the product space
${\mathbb{R}}^+ \times F$, with $(F,\mathcal{F})$ a Borel subspace of
${\mathbb{R}}^d$. Then, a random counting measure is defined by a random
set $S =\{ (T_k, \Theta_k ), k \geq 1\}$. The random variables
$T_k\geq 0$ can be considered as time variables, and constitute the jump
times of the random measure, while the variables $\Theta_k \in F$
represent space variables.

We recall in this special case the Theorem VI.3.2 in [@Cin11].

::: {#PropMarkedPoisson .prop}
**Proposition 11** (Marked Poisson measure). *Let $m$ be a
$\sigma$--finite diffuse measure on ${\mathbb{R}}^+$, and $K$ a
transition probability kernel from
$({\mathbb{R}}^+,\mathcal{B}({\mathbb{R}}^+))$ into $(F, \mathcal{F})$.
Assume that the collection $(T_k)_{k \ge 1}$ forms a Poisson process
$(N_t) =(\sum_{k\geq 1} \mathsf{1}_{\{T_k \leq t\}})$ with mean
$m(\mathrm{d}t) =\Lambda(t) \mathrm{d}t$, and that given
$(T_k)_{k \ge 1}$, the variables $\Theta_k$ are conditionally
independent and have the respective distributions $K(T_k, \cdot)$.*

1.  *Then, $\{ (T_k, \Theta_k) ;\; k \ge 1\}$ forms a Poisson random
    measure $Q = \sum_{k\ge 1} \delta_{(T_k, \Theta_k)}$ on
    $({\mathbb{R}}^+ \times F, \mathcal{B}({\mathbb{R}}^+) \otimes \mathcal{F})$,
    called a *Marked point process* , with mean $\mu$ defined by
    $$\mu(\mathrm{d}t, \mathrm{d}y) = \Lambda(t) \mathrm{d}t K(t, \mathrm{d}y).$$*

2.  *Reciprocally let $Q$ be a Poisson random measure of mean measure
    $\mu(\mathrm{d}t, \mathrm{d}y)$, admitting the following
    disintegration with respect to the first coordinate:
    $\mu(\mathrm{d}t , \mathrm{d}y) =\tilde  \Lambda(t) \mathrm{d}t \nu(t, \mathrm{d}y)$,
    with $\nu(t, F)<\infty$. Let
    $K(t, \mathrm{d}y) = \dfrac{\nu(t,\mathrm{d}y) }{\nu(t, F) }$ and
    $\Lambda(t) = \nu(t, F)\tilde  \Lambda(t)$. Then,
    $Q = \sum_{k\ge 1} \delta_{(T_k, \Theta_k)}$ is a marked Poisson
    measure with $(T_k,\Theta_k)_{k\in {\mathbb{N}}^*}$ defined as
    above. In particular, the projection $N= (N_t)_{t\geq0}$ of the
    Poisson measure on the first coordinate,
    $$N_t = Q([0,t] \times F) = \sum_{k\geq 1} \boldsymbol{1}_{[0,t] \times F} (T_k, \Theta_k)  = \sum_{k\geq 1} \boldsymbol{1}_{\{T_k \leq t\}}, \quad \forall \; t \geq 0,$$
    is an inhomogeneous Poisson process of rate
    $\Lambda(t)= \nu(t, F)\tilde  \Lambda(t)$.*
:::

::: remark
*Remark 9*. When the transition probability kernel $K$ does not depend
on the time: $K(t, A) = \nu(A)$ for some probability measure $\nu$, then
the marks $(\Theta_k)_{k \ge 1}$ form an *i.i.d.* sequence with
distribution $\nu$, independent of $(T_k)_{k \ge 1}$.
:::

The preceding proposition thus yields a straight forward iterative
simulation procedure for a Marked Poisson process on $[0,T]\times F$
with mean measure
$\mu(\mathrm{d}t, \mathrm{d}y) = c \mathrm{d}t K(t, \mathrm{d}y)$
($c>0$):\

::: algorithm
Initialization: draw $T_1 \sim \mathcal{E}(c)$ and draw
$Y_1 \sim K(T_1, \mathrm{d}y)$
:::

<figure id="plot:poisson">

<figcaption>Example of Marked Poisson measure on <span
class="math inline">[0,<em>T</em>]</span> with <span
class="math inline"><em>m</em>(d<em>t</em>) = <em>L</em>d<em>t</em></span>
(jump times occur at Poisson arrival times of rate <span
class="math inline"><em>L</em></span>) and with <span
class="math inline">$\nu(\mathrm{d}y) = \frac{1}{L} \mathbf{1}_{[0,
L]}(y) \mathrm{d}y$</span> (marks are drawn uniformly on <span
class="math inline">[0,<em>L</em>]</span>). The mean measure is then
<span
class="math inline"><em>μ</em>(d<em>t</em>,d<em>y</em>) = d<em>t</em><strong>1</strong><sub>[0,<em>L</em>]</sub>(<em>y</em>)d<em>y</em></span>.</figcaption>
</figure>

# Pathwise representation of IBMs {#pathwiserepresentation}

#### Notation reminder

The population's evolution is described by the measure valued process
$(Z_t)_{t\geq 0}$. Several types of events $e$ can occur to individuals
denoted by $I$. In an event of type $e$ occur to the individual $I$ at
time $t$, then the population state $Z_{t^-}$ is modified by
$\phi^e(t,I)$. If $e\in \mathcal{E} \cup \mathcal{E}_W$, then events of
type $e$ occur with an intensity $\sum_{k=1}^{N_t} \lambda_t^e(I,Z_t)$,
with $\lambda_t^e(I,Z_t)$ defined by
[\[IndividualIntensity\]](#IndividualIntensity){reference-type="eqref"
reference="IndividualIntensity"}. If $e \in \mathcal{P}$, then events of
type $e$ occur in the population at a Poisson intensity of $(\mu^e_t)$.

## Proof of Theorem [1](#ThEqZ){reference-type="ref" reference="ThEqZ"} {#ProofThPathwise}

::: proof
*Proof of Theorem [1](#ThEqZ){reference-type="ref" reference="ThEqZ"}.*
For ease of notation, we prove the case when $\mathcal{P} =\emptyset$
(there are no events with Poisson intensity).

**Step 1** The existence of a solution to
[\[SDE_pop\]](#SDE_pop){reference-type="eqref" reference="SDE_pop"} is
obtained by induction. Let $Z^1$ be the unique solution the thinning
equation:
$$Z_t^1 = Z_0 + \int_0^t \int_{ \mathcal{J}\times \mathbb R^+ }\phi^e (s , I_k)  \mathbf{1}_{\{k \leq N_{0}\} }\mathbf{1}_{\{\theta \leq \lambda_s^e(I_k, Z_{0})\}} Q (\mathrm{d}s ,\mathrm{d}k , \mathrm{d}e, \mathrm{d}\theta ), \quad \forall  0 \leq t \leq T.$$
Let $T_1$ be the first jump time of $Z^1$. Since $Z_{s^-}^1 = Z_{0}$ and
$N_{s^-}=N_{0}$ on $[0, T_1]$, $Z^1$ is solution of
[\[SDE_pop\]](#SDE_pop){reference-type="eqref" reference="SDE_pop"} on
$[0,T_1]$.

Let us now assume that [\[SDE_pop\]](#SDE_pop){reference-type="eqref"
reference="SDE_pop"} admits a solution $Z^n$ on $[0,T_n]$, with $T_n$
the $n$--th event time in the population. Let $Z^{n+1}$ be the unique
solution of the thinning equation:
$$Z^{n+1}_t  =  Z_{t\wedge T_n}^n + \int_{t\wedge T_n}^{t} \int_{ \mathcal{J}\times \mathbb R^+ }\phi^e (s , I_k)  \mathbf{1}_{\{\theta \leq \lambda_s^e(I_k, Z_{T_n}^n)\}} \mathbf{1}_{\{k \leq N_{T_n}^n \} }Q (\mathrm{d}s ,\mathrm{d}k , \mathrm{d}e, \mathrm{d}\theta ).$$
First, observe that $Z^{n+1}$ coincides with $Z^n$ on $[0,T_n]$. Let
$T_{n+1}$ be the $(n+1)$--th jump of $Z^{n+1}$. Furthermore,
$Z_{s^-}^{n+1} = Z_{T_n}^n$ and $N_{s^-}^{n+1}=N_{T_n}^{n}$ on
$[T_n, T_{n+1}]$ (nothing happens between two successive event times),
$Z^{n+1}$ verifies for all $t\leq T_{n+1}$: $$\begin{aligned}
Z^{n+1}_t  =  Z_{t\wedge T_n}^n +\int_{t\wedge T_n}^{t} \int_{ \mathcal{J}\times \mathbb R^+ }\phi^e (s , I_k)  \mathbf{1}_{\{\theta \leq \lambda_s^e(I_k, Z_{s^-}^{n+1} )\}} \mathbf{1}_{\{k \leq N_{s^-}^{n+1} \} }Q (\mathrm{d}s ,\mathrm{d}k , \mathrm{d}e, \mathrm{d}\theta ).
\end{aligned}$$ Since, $Z^n$ is a solution of
[\[SDE_pop\]](#SDE_pop){reference-type="eqref" reference="SDE_pop"} on
$[0,T_n]$ coinciding with $Z^{n+1}$, this achieves to prove that
$Z^{n+1}$ is solution of [\[SDE_pop\]](#SDE_pop){reference-type="eqref"
reference="SDE_pop"} on $[0,T_{n+1}]$.\
Finally, let $Z =\lim_{n\to \infty } Z^n$. For all $n\geq 1$, $T_n$ is
the $n$--th event time of $Z$, and $Z$ is solution of
[\[SDE_pop\]](#SDE_pop){reference-type="eqref" reference="SDE_pop"} on
all time intervals $[0,T_n\wedge T]$ by construction.\
By Lemma [2](#lemma:nonExplosionSDE){reference-type="ref"
reference="lemma:nonExplosionSDE"},
$T_n \underset{n\to \infty}{\longrightarrow} \infty$. Thus, by letting
$n\to \infty$ we can conclude that $Z$ is a solution of
[\[SDE_pop\]](#SDE_pop){reference-type="eqref" reference="SDE_pop"} on
$[0,T]$.\
**Step 2** Let $\tilde Z$ be a solution of
[\[SDE_pop\]](#SDE_pop){reference-type="eqref" reference="SDE_pop"}.
Using the same arguments than in Step 1, it is straight forward to show
that $\tilde Z$ coincides with $Z^n$ on $[0,T_n]$, for all $n \geq 1$.
Thus, $\tilde{Z} = Z$, with achieves to prove uniqueness. ◻
:::

## Proof of Lemma [2](#lemma:nonExplosionSDE){reference-type="ref" reference="lemma:nonExplosionSDE"} {#section:proofLemma}

The proof is obtained using pathwise comparison result, generalizing
those obtained in [@KaaElK20].

::: proof
*Proof of Lemma [2](#lemma:nonExplosionSDE){reference-type="ref"
reference="lemma:nonExplosionSDE"}.* Let $Z$ be a solution of
[\[SDE_pop\]](#SDE_pop){reference-type="eqref" reference="SDE_pop"}. For
all $e \in \mathcal{P} \cup \mathcal{E} \cup  \mathcal{E}_W$, let $N^e$
be the process counting the occurrence of events of type $e$ in the
population. $N^e$ is a counting process of $\{\mathcal{F}_t\}$-intensity
$(\Lambda_t^e(Z_{t^-}))$, solution of $$\begin{aligned}
 \label{eq:Ne} & N_t^e = \int_0^t \int_{{\mathbb{N}}\times {\mathbb{R}}^+} \boldsymbol{1}_{\{k\leq N_{s^-}\}}  \boldsymbol{1}_{\{\theta \leq \lambda_s^e(I_k, Z_{s^-})\}} Q(\mathrm{d}s, \mathrm{d}k, \{e\}, \mathrm{d}\theta), & \quad  \textit{ if } e \in \mathcal{E}\cup \mathcal{E}_W, \\
\nonumber& N_t^e = \int_0^t \int_{{\mathbb{R}}^+} \boldsymbol{1}_{\{\theta \leq \mu^e_s \}} Q^{\mathcal{P}}(\mathrm{d}s, \{e\}, \mathrm{d}\theta), & \quad  \textit{ if } e \in \cal P. \\
\end{aligned}$$ By definition, the jump times of the multivariate
counting process
$(N^e)_{e \in \mathcal{P} \cup \mathcal{E}\cup \mathcal{E}_W}$ are the
population event times $(T_n)_{n\geq 0}$. The idea of the proof is to
show that $(N^e)_{e \in \mathcal{P} \cup \mathcal{E}\cup \mathcal{E}_W}$
does not explode in finite time, by pathwise domination with a simpler
multivariate counting process. The first steps are to control the
population size $N_t = N_0 + N^b_t + N^{en}_t$.\
**Step 1** Let $(\bar N^b, \bar N^{en})$ be the 2-dimensional counting
process defined as follows: for $e \in \{b,en\}$, $\bar N^e_0 = 0$ and
$$\begin{aligned}
\label{eq:dominatingprocess}
& \bar N_t^e = \int_0^t \int_{{\mathbb{N}}\times {\mathbb{R}}^+} \boldsymbol{1}_{\{k\leq N_0 + \bar N_{s^-} \}}  \boldsymbol{1}_{\{\theta \leq f^e(N_0 + \bar N_{s^-})\}} Q(\mathrm{d}s, \mathrm{d}k, \{e\}, \mathrm{d}\theta),  \quad  \textit{ if } e \in \mathcal{E}\cup \mathcal{E}_W, \\
& \nonumber  \bar N_t^e = \int_0^t \int_{{\mathbb{R}}^+} \boldsymbol{1}_{\{\theta \leq \bar \mu^e\}} Q^{\mathcal{P}}(\mathrm{d}s, \{e\}, \mathrm{d}\theta)  \quad  \textit{ if } e \in \cal P,
\end{aligned}$$ with $\bar N := \bar N^b + \bar N^{en}$ and $f^e$ the
function introduced in Assumption
[4](#Assumption:nonExplosion){reference-type="ref"
reference="Assumption:nonExplosion"}.\
- If $b,en \in \cal P$, then $\bar N$ is a inhomogeneous Poisson
process.\
- If $b,en \in \mathcal{E}\cup \mathcal{E}_W$, then it is
straightforward to show that conditionally to $N_0$, $\bar N$ is a pure
birth Markov process with birth intensity function
$g(n) = n\big(f^b(N_0+n) + f^{en}(N_0+n)\big)$. In particular, by
Assumption [4](#Assumption:nonExplosion){reference-type="ref"
reference="Assumption:nonExplosion"}, $g$ verifies the standard Feller
condition for pure birth Markov processes (see e.g. [@BanMel15]):
$$\sum_{n=1}^\infty \frac{1}{g(n)}.$$ - Finally, if $b \in \mathcal{E}$
and $en \in \cal P$ (or equivalently if $b \in \cal P$ and
$en \in \mathcal{E}$), then one can show easily that $\bar N$ is a pure
birth Markov process with immigration, of birth intensity function
$g(n)= \bar  \mu^{en} + n f^b(N_0 + n)$ (resp.
$g(n)= \bar  \mu^b+ n f^{en}(N_0 + n)$), also verifying the Feller
condition.

Therefore, there exists a non-exploding solution of
[\[eq:dominatingprocess\]](#eq:dominatingprocess){reference-type="eqref"
reference="eq:dominatingprocess"}, by Proposition 3.3 in [@KaaElK20].\
**Step 2** The second step consists in showing that $(N^b, N^{en})$ is
strongly dominated by $(\bar N^b, \bar N^{en})$, i.e that all jumps of
$(N^b, N^{en})$ are jumps of $(\bar N^b, \bar N^{en})$. Without loss of
generality, we can assume that $f^e:\mathbb{N} \to (0,+\infty)$ is
increasing since $f^e(n)$ can be replaced by
$\sup_{\{m\leq n \} } f^e(m)$.\
Let $e\in \{b, en\}$. If $e \in \mathcal{P}$, then for all $s\in [0,T]$
$$\{\theta \leq \mu_s^e\} \subset \{ \theta \leq \bar \mu^e\},$$ which
yields that all jumps of $N^e$ are jumps of $\bar N^e$.\
If $e \in  \mathcal{E}\cup \mathcal{E}_W$, the proof by induction is
analogous to the proof of Proposition 2.1 in [@KaaElK20]. Let $T_1^e$ be
first jump time of $N^e$, associated with the marks $(K_1^e,\Theta_1^e)$
of $Q$ (or $Q^{\mathcal{P}}$). Then, by Definition of
[\[eq:Ne\]](#eq:Ne){reference-type="eqref"
reference="eq:Ne"},$K_1^e \leq N_0$ and
$\Theta_1^e \leq \lambda_{T_1^e}^e (I_{K_1^e}, Z_0)$.\
By Assumption [4](#Assumption:nonExplosion){reference-type="ref"
reference="Assumption:nonExplosion"}, we have also
$$\Theta_1^e \leq \lambda_{T_1^e}^e (I_{K_1^e}, Z_0) \leq f^e(N_0) \leq f^e(N_0 +\bar N_{T_1^{e,-}}), \quad K_1^e \leq N_0 +  \bar N_{T_1^{e,-}}.$$
Thus, $T_1^e$ is also a jump time of $\bar N^e$. By iterating this
argument, we obtain that all jump times of $N^e$ are jump times of
$\bar N^e$.\
Thus, $(N^b, N^{en})$ does not explode in finite time.\
**Step 3** It remains to show that for $e \notin \{b, en\}$, $N^e$ does
not explode.\
Let $e \neq b, en$. If $e\in \mathcal P$, the proof is the same than in
Step 2. Otherwise, let:
$$h^e_t(n) = \sup_{I \in \mathcal{I},  m \leq n } \lambda^e_t \biggl(I, \sum_{k=1}^{m} \delta_{I_k}\biggr), \quad \forall \; t \in [0,T] \;  n \in {\mathbb{N}}^*.$$
By Assumptions [2](#AssumptionIntensity1){reference-type="ref"
reference="AssumptionIntensity1"} and
[3](#AssumptionIntensity2){reference-type="ref"
reference="AssumptionIntensity2"}, $h^e_t(n) <\infty$, and we can
introduce the non exploding counting process $\bar{N}^e$, defined by the
thinning equation :
$$\bar N_t^e = \int_0^t \int_{{\mathbb{N}}\times {\mathbb{R}}^+} \boldsymbol{1}_{\{k\leq N_0 + \bar N_{s^-} \}}  \boldsymbol{1}_{\{\theta \leq h^e_s (N_0 + \bar N_{s^-})\}} Q(\mathrm{d}s, \mathrm{d}k, \{e\}, \mathrm{d}\theta),$$
with $\bar N_s = \bar N^b_s + \bar N^{en}_s$.\
Finally, by Step 2, for $s\in [0,T]$ the population size
$N_s = N_0 + N^b_s+ N^{en}_s$ is bounded a.s. by $N_0 + \bar N_s$, since
all jumps of $(N^b,N^{en})$ are jumps of $(\bar N^b, \bar N^{en})$.
Thus, for all $s\in [0,T]$,
$$\{k \leq N_{s^-} \}\subset  \{k\leq N_0 + \bar N_{s^-} \}, \text{ and }  \{\theta \leq \lambda_s^e(I_k, Z_{s^-})\} \subset \{\theta \leq h_s^e(N_0 + \bar N_{s^-})\}.$$
This proves that all jumps of $N^e$ are jumps $\bar N^e$, and thus $N^e$
does not explode in finite time. ◻
:::

## Alternative pathwise representation

::: {#ThEqZrandomized .theo}
**Theorem 12**. *Let
$\mathcal{J}_{\mathcal{E}} = \mathbb N \times \mathcal{E}$ and
$\mathcal{J}_W  = \mathbb N \times \mathcal{E}_W$.\
Let $Q^{\mathcal{E}}$ be a random Poisson measure on
$\mathbb R^+ \times \mathcal{J}_\mathcal{E}\times \mathbb{R}^+$, of
intensity
$\mathrm{d}t \delta_{\mathcal{J}_{\mathcal{E}}}(\mathrm{d}k, \mathrm{d}e)  \mathbf{1}_{[0,\bar \lambda^e]} (\theta) \mathrm{d}\theta$,
and $Q^{W}$ a random Poisson measure on
$\mathbb R^+ \times \mathcal{J}_W \times \mathbb{N} \times  \mathbb{R}^+$,
of intensity
$\mathrm{d}t \delta_{\mathcal{J}_{\mathcal{E}}}(\mathrm{d}k,\mathrm{d}e)) \delta_{\mathbb{N}} (\mathrm{d}j) \mathbf{1}_{[0,\bar W^e]} (\theta)\mathrm{d}\theta$.
Finally, let $Q^{\mathcal P}$ be a random Poisson measure on
$\mathbb R^+ \times \mathcal{P}  \times \mathbb{R}^+$, of intensity
$\mathrm{d}t \delta_{\cal P}(\mathrm{d}e)  \mathbf{1}_{[0,\bar \mu^e]} (\theta)\mathrm{d}\theta$.\
There exists a unique measure-valued process $Z$, strong solution on the
following SDE driven by Poisson measure: $$\begin{aligned}
\label{SDE_pop_randomized}
\nonumber Z_t  = Z_0 &  + \int_0^t \int_{\mathcal{J}_\mathcal{E}\times \mathbb R^+ }\phi^e (s , I_k)  \mathbf{1}_{\{k \leq N_{s^-}\} }\mathbf{1}_{\{\theta \leq \lambda_s^e(I_k, Z_{s^-})\}} Q^\mathcal{E}(\mathrm{d}s ,\mathrm{d}k , \mathrm{d}e, \mathrm{d}\theta ) \\
&  + \int_0^t \int_{\mathcal{J}_W  \times {\mathbb{N}}\times  \mathbb R^+ }\phi^e (s , I_k)  \mathbf{1}_{\{k \leq N_{s^-}\} } \mathbf{1}_{\{j \leq N_{s^-}\} }\mathbf{1}_{\{\theta \leq W^e(s, I_k , I_j) \}} Q^W (\mathrm{d}s ,\mathrm{d}k , \mathrm{d}e,  \mathrm{d}j ,\mathrm{d}\theta ),\\
& +  \nonumber    \int_0^t \int_{\mathcal{P} \times \mathbb R^+}  \phi^e(s, I_{s^-}) \mathbf{1}_{\{\theta \leq \mu_s^e \}} Q^{\mathcal{P}} (\mathrm{d}s ,\mathrm{d}e , \mathrm{d}\theta),
\end{aligned}$$ with $I_{s^-}$ an individual taken uniformly in
$Z_{s^-}$.\
Furthermore, the solution of
[\[SDE_pop_randomized\]](#SDE_pop_randomized){reference-type="eqref"
reference="SDE_pop_randomized"} has the same law than the solution of
Equation [\[SDE_pop\]](#SDE_pop){reference-type="eqref"
reference="SDE_pop"}.*
:::

The proof Theorem [12](#ThEqZrandomized){reference-type="ref"
reference="ThEqZrandomized"} follows the same steps than the proof of
Theorem [1](#ThEqZ){reference-type="ref" reference="ThEqZ"}.

# Proof of Theorem [5](#theorem:algoNoInteraction){reference-type="ref" reference="theorem:algoNoInteraction"} {#proof:algonointeraction}

::: proof
*Proof of Theorem [5](#theorem:algoNoInteraction){reference-type="ref"
reference="theorem:algoNoInteraction"}.* For ease of notation, we prove
the case when $\mathcal{P} =\emptyset$ (there are no events with Poisson
intensity).\
Let $Z$ be the population process obtained by Algorithm
[\[algo:PopNointeraction2\]](#algo:PopNointeraction2){reference-type="ref"
reference="algo:PopNointeraction2"}, and $(T_n)_{n\geq 0}$ the sequence
of its jump times ($T_0=0$).\
**Step 1** Let $T_1$ be the first event time in the population, with its
associated marks defining the type $E_1$ of the event and the individual
$I_1$ to which this event occurs. By construction, $(T_1, E_1, I_1)$ is
characterized by the first jump of:
$$Q^0(\mathrm{d}t, \mathrm{d}k , \mathrm{d}e) = \int_{\mathbb R^+} \mathbf{1}_{\{\theta \leq \lambda_{t}^e(I_k,Z_0)\}}\bar Q^0 (\mathrm{d}t ,\mathrm{d}k , \mathrm{d}e, \mathrm{d}\theta ),$$
with $\bar Q^0$ the Poisson measure introduced in the first step of the
algorithm described in Section
[2.2](#sec::simulation_algo){reference-type="ref"
reference="sec::simulation_algo"}.

Since $T_1$ is the first event time, the population composition stays
constant, $Z_t=Z_0$, on $\{t<T_1\}$. In addition, recalling that the
first event has the action $\phi^{E_1}(T_1, I_1)$ (see Table
[1](#TableEvAction){reference-type="ref" reference="TableEvAction"}) on
the population $Z$, we obtain that: $$\begin{aligned}
Z_{t\wedge T_1} & =  Z_0 + \mathbf{1}_{\{t\geq T_1\}} \phi^{E_1} (T_1 , I_1)  \\
 & = Z_0 + \int_0^{t\wedge T_1}  \int_{\mathcal{J}_0} \phi^e (s , I_k)  Q^0 (\mathrm{d}s ,\mathrm{d}k , \mathrm{d}e ) \\
& = Z_0 + \int_0^{t\wedge T_1} \int_{\mathcal{J}_0}  \int_{\mathbb R^+} \phi^e (s , I_k)  \mathbf{1}_{\{\theta \leq \lambda_s^e(I_k,Z_0)\}}\bar Q^0 (\mathrm{d}s ,\mathrm{d}k , \mathrm{d}e, \mathrm{d}\theta ).
\end{aligned}$$ Since $Z_{s^-} = Z_0$ on $\{s \leq T_1\}$, the last
equation can be rewritten as $$\label{eq:DynZT_1}
 Z_{t\wedge T_1}   = Z_0 + \int_0^{t\wedge T_1} \int_{\mathcal{J}_0}  \int_{\mathbb R^+} \phi^e (s , I_k)  \mathbf{1}_{\{\theta \leq \lambda_s^e(I_k,Z_{s^-})\}}\bar Q^0 (\mathrm{d}s ,\mathrm{d}k , \mathrm{d}e, \mathrm{d}\theta ).$$

**Step 2** The population size at the $n$--th event time $T_n$ is
$N_{T_n}$. The $(n+1)$--th event type and the individual to which this
event occur are thus chosen in the set
$$\mathcal{J}_n := \{ 1,\dots, N_{T_n}\} \times (\mathcal{E} \cup \mathcal{E}_W).$$
Conditionally to $\mathcal{F}_{T_n}$, let us first introduce the marked
Poisson measure $\bar Q^n$ on
$[T_n, \infty) \times \mathcal J_n \times \mathbb R^+$, of intensity:
$$\begin{aligned}
\label{eq:barmun}
\bar \mu^n(\mathrm{d}t, \mathrm{d}k, \mathrm{d}e , \mathrm{d}\theta ) & := \mathbf{1}_{\{t > T_n \}}\bar \Lambda (N_{T_n})\mathrm{d}t  \frac{\bar \lambda^e }{\bar \Lambda(N_{T_n})} \delta_{\mathcal J_n}(\mathrm{d}k, \mathrm{d}e) \frac{1}{\bar \lambda^e} \mathbf{1}_{[0,\bar \lambda^e]} (\theta)\mathrm{d}\theta,\\
\nonumber&  = \mathbf{1}_{\{t > T_n \}}\mathrm{d}t  \delta_{\mathcal J_n}(\mathrm{d}k, \mathrm{d}e) \mathbf{1}_{[0,\bar \lambda^e]}(\theta)\mathrm{d}\theta .
\end{aligned}$$ By definition, $\bar Q^n$ has no jumps before $T_n$.\
As for the first event, the triplet $(T_{n+1}, E_{n+1}, I_{n+1})$ is
determined by the first jump of the measure
$Q^n (\mathrm{d}s ,\mathrm{d}k , \mathrm{d}e) := \int_{\mathbb R^+} \mathbf{1}_{\{\theta \leq \lambda_s^e(I_k, Z_{T_n})\}}\bar Q^n (\mathrm{d}s ,\mathrm{d}k , \mathrm{d}e, \mathrm{d}\theta)$,
obtained by thinning of $\bar Q^n$. Finally, since the population
composition is constant on $[T_n, T_{n+1}[$, $Z_t = Z_{T_n}$, the
population on $[0,T_{n+1}]$ is defined by: $$\begin{aligned}
\label{EqZ_Tn}
\nonumber Z_{t\wedge T_{n+1}}  &  = Z_{t \wedge T_n }  + \mathbf{1}_{\{t\geq T_{n+1}\}}\phi^{E_{n+1}}(T_{n+1}, I_{n+1}), \\
& = Z_{t \wedge T_n } + \int_{t \wedge T_n}^{t \wedge T_{n+1}} \int_{\mathcal{J}_n
\times \mathbb R^+} \phi^e (s , I_k)  \mathbf{1}_{\{\theta\leq \lambda_s^e(I_k, Z_{s^-})\}}\bar Q^n (\mathrm{d}s ,\mathrm{d}k , \mathrm{d}e, \mathrm{d}\theta ).
\end{aligned}$$ Applying $n$ times
[\[EqZ_Tn\]](#EqZ_Tn){reference-type="eqref" reference="EqZ_Tn"} yields
that: $$\begin{aligned}
\label{EqZtnRec}
Z_{t\wedge T_{n+1}} = Z_0  + \sum_{l=0}^n \int_{t \wedge T_l}^{t \wedge T_{l+1}}\int_{ \mathcal{J}_l \times  \mathbb R^+} \phi^e (s , I_k)  \mathbf{1}_{\{\theta\leq \lambda_s^e(I_k, \tilde Z_{s^-})\}}\bar Q^l (\mathrm{d}s ,\mathrm{d}k , \mathrm{d}e, \mathrm{d}\theta ).
\end{aligned}$$

**Step 3** Finally, let $\tilde{Z}$ be the solution of
[\[SDE_pop\]](#SDE_pop){reference-type="eqref" reference="SDE_pop"},
with $(\tilde T_n)_{n\geq 0}$ the sequence of its event times. Then, we
can write similarly for all $n\geq 0$: $$\begin{aligned}
 \tilde Z_{t\wedge \tilde T_{n+1}} & = Z_0  + \sum_{l=0}^n \int_{t \wedge \tilde  T_l}^{t \wedge \tilde T_{l+1}} \int_{\mathcal{J}\times \mathbb R^+} \phi^e (s , I_k)  \mathbf{1}_{\{\theta\leq \lambda_s^e(I_k,  \tilde Z_{s^-})\}}\mathbf{1}_{ \{k \leq \tilde N_{s^-} \}}   Q(\mathrm{d}s, \mathrm{d}k , \mathrm{d}e , \mathrm{d}\theta ), \\
& = Z_0  + \sum_{l=0}^n \int_{t \wedge \tilde  T_l}^{t \wedge \tilde T_{l+1}} \int_{\mathcal{J}\times \mathbb R^+} \phi^e (s , I_k)  \mathbf{1}_{\{\theta\leq \lambda_s^e(I_k,  \tilde Z_{s^-})\}}\mathbf{1}_{ \{k \leq \tilde N_{\tilde T_l} \}}   Q(\mathrm{d}s, \mathrm{d}k , \mathrm{d}e , \mathrm{d}\theta ),
\end{aligned}$$ since $\tilde  N_{s^-} = \tilde N_{T_l}$ on
$\tilde [T_l, \tilde T_{l+1}]$.\
For each $l \geq 0$, let
$$\tilde Q^l(\mathrm{d}t, \mathrm{d}k , \mathrm{d}e , \mathrm{d}\theta )  = \mathbf{1}_{\{t >  \tilde T_l\}} \mathbf{1}_{ \{1, \dots , \tilde N_{\tilde T_l} \}}(k)   Q(\mathrm{d}t, \mathrm{d}k , \mathrm{d}e , \mathrm{d}\theta ).$$
By proposition [9](#PropRestrictionPoissonMeasure){reference-type="ref"
reference="PropRestrictionPoissonMeasure"}, $\tilde Q^l$ is,
conditionally to $\mathcal{F}_{T_l}$, a Poisson measure of intensity
$$\mathbf{1}_{\{t >  \tilde T_l\}}   \mathrm{d}t \mathbf{1}_{ \{1, \dots , \tilde N_{\tilde T_l} \}}(k) \delta_{\mathcal{J}}(\mathrm{d}k, \mathrm{d}e) \mathrm{d}\theta.$$
Noticing that
$\mathbf{1}_{ \{1, \dots , \tilde N_{\tilde T_l} \}}(k) \delta_{\mathcal{J}}(\mathrm{d}k, \mathrm{d}e)  =\delta_{\mathcal{J}_l}( \mathrm{d}k , \mathrm{d}e)$,
this shows that $\tilde{Q}^l$ has the conditional intensity $\bar \mu^l$
defined in [\[eq:barmun\]](#eq:barmun){reference-type="eqref"
reference="eq:barmun"} and has thus the same distribution than
$\bar Q^l$. Thus, $Z$ in an exact simulation of
[\[SDE_pop\]](#SDE_pop){reference-type="eqref" reference="SDE_pop"}. ◻
:::

[^1]: source: Human Mortality Database <https://www.mortality.org/>

[^2]: The choices
    $(\alpha, \beta) \in \{(1,2/300),(1, 1/300), (1.5, 1/300), (2, 1/300)\}$
    lead to the stationary population sizes
    $N \in \{360, 900, 1800, 2600\}$. For each set of parameters, we
    generated a new initial population, which was used for a benchmark
    of 50 simulations with both randomized and full algorithm. The
    simulations run on a Intel Core i7-8550U CPU 1.80GHz × 8 processor,
    with 15.3 GiB of RAM, under Debian GNU/Linux 11.
