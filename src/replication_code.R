#############################################################
## Replication code for article
##
## Daphné Giorgi, Sarah Kaakaï, Vincent Lemaire
## Efficient simulation of individual-based population models: the R Package IBMPopSim
##
#############################################################

## --------------------------------------------------------------------------------------------------------------------

## Package loading
library(IBMPopSim)
library(StMoMo)
library(reshape2)
library(ggplot2)

## --------------------------------------------------------------------------------------------------------------------
## BRIEF OVERVIEW
## --------------------------------------------------------------------------------------------------------------------


## --------------------------------------------------------------------------------------------------------------------
## Model creation
## --------------------------------------------------------------------------------------------------------------------


birth_event <- mk_event_individual(
  type = "birth",
  intensity_code = "result = birth_rate(I.age(t));",
  kernel_code = "newI.male = CUnif(0,1) < p_male;")

death_event <- mk_event_individual(
  type = "death",
  intensity_code = "result = alpha * exp(beta * I.age(t));")

params <- list(
  "alpha" = 0.008, "beta" = 0.02,
  "p_male" = 0.51,
  "birth_rate" = stepfun(c(15, 40), c(0, 0.05, 0)))

birth_death_model <- mk_model(
  characteristics = c("male" = "bool"),
  events = list(death_event, birth_event),
  parameters = params)

## --------------------------------------------------------------------------------------------------------------------
## Simulation
## --------------------------------------------------------------------------------------------------------------------

a_max <- 115
events_bounds = c(
  "death" = params$alpha * exp(params$beta * a_max),
  "birth" = max(params$birth_rate))

sim_out <- popsim(
  birth_death_model,
  population(EW_pop_14$sample),
  events_bounds,
  parameters = params, age_max = a_max,
  time = 30)

## --------------------------------------------------------------------------------------------------------------------
## PACKAGE DESCRIPTION
## --------------------------------------------------------------------------------------------------------------------


## --------------------------------------------------------------------------------------------------------------------
## Population
## --------------------------------------------------------------------------------------------------------------------

pop_init <- population(EW_pop_14$sample,entry=TRUE)
str(pop_init)

## --------------------------------------------------------------------------------------------------------------------
## Example
## --------------------------------------------------------------------------------------------------------------------

params <- list("p_male"= 0.51,
               "birth_rate" = stepfun(c(15,40),c(0,0.05,0)),
               "death_rate" = gompertz(0.008,0.02))

death_event <- mk_event_individual(type = "death", name= "my_death_event",
                  intensity_code = "result = death_rate(age(I,t));")

birth_event <- mk_event_individual( type = "birth",
                  intensity_code = "if (I.male)
                                        result = 0;
                                    else
                                        result=birth_rate(age(I,t));",
                  kernel_code = "newI.male = CUnif(0, 1) < p_male;")

pop <- population(EW_pop_14$sample)

model <- mk_model(characteristics = get_characteristics(pop),
                  events = list(death_event,birth_event),
                  parameters = params)

a_max <- 120 # maximum age
events_bounds <- c("my_death_event" = params$death_rate(a_max),
                   "birth" = max(params$birth_rate))

sim_out <- popsim(model, pop, events_bounds, params,
                  age_max = a_max, time = 30)

## --------------------------------------------------------------------------------------------------------------------
## INSURANCE PORTFOLIO
## --------------------------------------------------------------------------------------------------------------------

## --------------------------------------------------------------------------------------------------------------------
## Population
## --------------------------------------------------------------------------------------------------------------------

N <- 30000
pop_df <- data.frame("birth" = rep(-65,N), "death" = rep(NA,N),
                     "risk_cls" = rep(1:2,each=N/2))
pop_init <- population(pop_df, entry=TRUE, out=TRUE)

## --------------------------------------------------------------------------------------------------------------------
## Entry event
## --------------------------------------------------------------------------------------------------------------------

entry_params <- list("lambda" = 30000, "p" = 0.5)
entry_event <- mk_event_poisson(
    type = "entry",
    intensity = "lambda",
    kernel_code = "if (CUnif() < p) newI.risk_cls =1;
                   else newI.risk_cls= 2;
                   double a = CUnif(65, 70);
                   newI.set_age(a, t);")

## --------------------------------------------------------------------------------------------------------------------
## Death event
## --------------------------------------------------------------------------------------------------------------------

# StMoMo death rates
EWStMoMoMale <- StMoMoData(EWdata_hmd, series = "male")
LC <- lc()
ages.fit <- 65:100
years.fit <- 1950:2016
LCfitMale <- fit(LC, data = EWStMoMoMale, ages.fit = ages.fit, years.fit = years.fit)
t <- 30
LCforecastMale <- forecast(LCfitMale, h = t)
d_k <- apply(LCforecastMale$rates, 2, function(x) stepfun(66:100, x))
breaks <- 1:29
death_male <- piecewise_xy(breaks,d_k)

death_max <- max(sapply(d_k, function(x) { max(x) }))

death_params <- list("death_male" = death_male, "alpha" = c(1.2, 0.8))
death_event <- mk_event_individual(
    type = "death",
    intensity_code = "result = alpha[I.risk_cls-1] * death_male(t, I.age(t));")

## --------------------------------------------------------------------------------------------------------------------
## Exit event
## --------------------------------------------------------------------------------------------------------------------

exit_params = list("mu" = c(0.001, 0.06))
exit_event <- mk_event_individual(
    type = "exit",
    intensity_code = "result = mu[I.risk_cls-1]; ")

## --------------------------------------------------------------------------------------------------------------------
## Model creation and simulation
## --------------------------------------------------------------------------------------------------------------------

model <- mk_model(
    characteristics = get_characteristics(pop_init),
    events = list(entry_event, death_event, exit_event),
    parameters = c(entry_params, death_params, exit_params))

bounds <- c("entry" = entry_params$lambda,
            "death" = death_max,
            "exit" = max(exit_params$mu))

sim_out <- popsim(
    model = model,
    initial_population = pop_init,
    events_bounds = bounds,
    parameters = c(entry_params, death_params, exit_params),
    time = 30,
    age_max = 110,
    multithreading = TRUE)


## --------------------------------------------------------------------------------------------------------------------
## Outputs
## --------------------------------------------------------------------------------------------------------------------
dim(population_alive(sim_out$population,t = 30))

sim_out$logs[["duration_ns"]]/1e9

sim_out$logs[["effective_events"]] / sim_out$logs[["proposed_events"]]

# Age pyramid (Figure 2(a))

age_grp <- 65:98
pyr = age_pyramid(sim_out$population, time = 30, ages=age_grp)
colnames(pyr)[2]<- "group_name"
pyr$group_name <- as.character(pyr$group_name)
colors <- c("1"="#00AFBB","2"="#FC4E07")
plot(pyr,colors,age_breaks = as.integer(seq(1,length(age_grp)-1,by=2)))

# Portofolio death rates

age_grp <- 65:95
Dx_pop <- death_table(sim_out$population, ages = age_grp, period = 0:30)
Ex_pop <- exposure_table(sim_out$population, ages = age_grp, period = 0:30)
mx_pop <- Dx_pop/Ex_pop

# Figure 2(b)

Dx <- death_table(sim_out$population[sim_out$population$risk_cls==1,],
                  ages = age_grp, period = 0:30)
Ex <- exposure_table(sim_out$population[sim_out$population$risk_cls==1,],
                     ages = age_grp, period = 0:30)


LC <- lc()
LCfitSim1 <- fit(LC, Dxt = Dx , Ext = Ex,ages=age_grp[-length(age_grp)])
time <- 30
age <- age_grp[-length(age_grp)]
estimated <- LCfitSim1$ax + LCfitSim1$bx*LCfitSim1$kt[time]
forecast_1 <- log(1.3*LCforecastMale$rates[1:(length(age_grp)-1),time])
forecast_2 <- log(0.8*LCforecastMale$rates[1:(length(age_grp)-1),time])
portfolio <- log((mx_pop)[,time])


df <- data.frame(age=age, portfolio=portfolio, forecast_1=forecast_1, forecast_1=forecast_2)
colnames(df) <- c("Age", "Portfolio", "Risk class 1", "Risk class 2")
df_melt <- melt(df, id="Age")

ggplot(data=df_melt, aes(x=Age, y=value)) +
  geom_point(aes(color=variable, shape=variable)) +
  geom_line(aes(color=variable, linetype=variable)) +
  xlab("Age") +
  ylab("Log mortality rates") +
  scale_color_manual(values = c("Portfolio" = "black", "Risk class 1" = "blue", "Risk class 2" = "red")) +
  scale_shape_manual(values = c("Portfolio" = 1, "Risk class 1" = NA, "Risk class 2" = NA)) +
  scale_linetype_manual(values = c("Portfolio" = 0, "Risk class 1" = 2, "Risk class 2" = 1)) +
  theme(legend.title = element_blank(), legend.position = c(0.9,0.2), plot.title = element_text(hjust = 0.5))

## --------------------------------------------------------------------------------------------------------------------
## POPULATION WITH GENETICALLY VARIABLE TRAITS
## --------------------------------------------------------------------------------------------------------------------

## --------------------------------------------------------------------------------------------------------------------
## Population
## --------------------------------------------------------------------------------------------------------------------

N <- 900
x0 <- 1.06
agemin <- 0.
agemax <- 2.

pop_df <- data.frame(
  "birth" = -runif(N, agemin, agemax), # Uniform age in [0,2]
  "death" = as.double(NA), # All individuals are alive
  "birth_size" = x0) # All individuals have the same initial birth size x0
pop_init <- population(pop_df)

## --------------------------------------------------------------------------------------------------------------------
## Birth event
## --------------------------------------------------------------------------------------------------------------------

birth_params <- list("p" = 0.03, "sigma" = sqrt(0.01), "alpha" = 1)

birth_event <- mk_event_individual(
  type = "birth",
  intensity_code = "result = alpha*(4 - I.birth_size);",
  kernel_code = "if (CUnif() < p)
                   newI.birth_size = min(max(0.,CNorm(I.birth_size,sigma)),4.);
                 else
                   newI.birth_size = I.birth_size;")

## --------------------------------------------------------------------------------------------------------------------
## Death event
## --------------------------------------------------------------------------------------------------------------------

death_params <- list("g" = 1, "beta" = 2./300., "c" = 1.2)
death_event <- mk_event_interaction(
  type = "death",
  interaction_code = "double x_I = I.birth_size + g * age(I,t);
                      double x_J = J.birth_size + g * age(J,t);
                      result = beta*(1.-1./(1.+c*exp(-4.*(x_I-x_J))));")

## --------------------------------------------------------------------------------------------------------------------
## Model creation and simulation
## --------------------------------------------------------------------------------------------------------------------

model <- mk_model(
    characteristics = get_characteristics(pop_init),
    events = list(birth_event, death_event),
    parameters = c(birth_params, death_params))

sim_out <- popsim(model = model,
    initial_population = pop_init,
    events_bounds = c("birth" = 4 * birth_params$alpha,
                      "death" = death_params$beta),
    parameters = c(birth_params, death_params),
    age_max = 2,
    time = 500)


# Figure 3 (a) 

pop_out <- sim_out$population

ggplot(pop_out) + geom_segment(
  aes(x=birth, xend=death, y=birth_size, yend=birth_size),
  na.rm=TRUE, colour="blue", alpha=0.1) +
  xlab("Time") +
  ylab("Birth size")

# Comparison between randomized and full algorithm 


death_event_full <- mk_event_interaction(type = "death",
                                         interaction_type= "full",
                                         interaction_code = "double x_I = I.birth_size + g * age(I,t);
                      double x_J = J.birth_size + g * age(J,t);
                      result = beta * ( 1.- 1./(1. + c * exp(-4. * (x_I-x_J))));"
)

model_full <- mk_model(characteristics = get_characteristics(pop_init),
                       events = list(birth_event, death_event_full),
                       parameters = c(birth_params, death_params))

sim_out_full <- popsim(model = model_full,
                       initial_population = pop_init,
                       events_bounds =c("birth" = 4 * birth_params$alpha,
                                        "death" = death_params$beta),
                       parameters = c(birth_params, death_params),
                       age_max = 2,
                       time = 500)

sim_out_full$logs["duration_ns"]/sim_out$logs["duration_ns"]




