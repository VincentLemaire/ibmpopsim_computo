
# Repository for Computo submission 
    
## Efficient simulation of individual-based population models

The `R` Package `IBMPopSim` aims to simulate the random evolution of heterogeneous populations using stochastic Individual-Based Models (IBMs). The package enables users to simulate population evolution, in which individuals are characterized by their age and some characteristics, and the population is modified by different types of events, including births/arrivals, death/exit events, or changes of characteristics. The frequency at which an event can occur to an individual can depend on their age and characteristics, but also on the characteristics of other individuals (interactions). Such models have a wide range of applications in fields including  actuarial science, biology, ecology or epidemiology. `IBMPopSim` overcomes the limitations of time-consuming IBMs simulations by implementing new efficient algorithms  based on thinning methods, which are compiled using the `Rcpp` package while providing a user-friendly interface.

To build your document (both in PDF and HTML) you can run the following commands 
```
quarto add computorg/computo-quarto-extension
quarto render article.qmd 
```


