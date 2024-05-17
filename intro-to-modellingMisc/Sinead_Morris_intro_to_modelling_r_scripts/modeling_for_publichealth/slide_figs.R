library(tidyverse)
library(deSolve)

## Plot settings ---------------------------------------------------------------

mycols3 <- c("grey40", "goldenrod3", "forestgreen")

txtsize <- 18

mytheme <- theme_light() + theme(legend.title = element_blank(),
                                 legend.text = element_text(size = txtsize),
                                 axis.text = element_text(size = txtsize),
                                 axis.title = element_text(size = txtsize))


## ODE equations ---------------------------------------------------------------

sir <- function(t, y, parms) {
  
  # use the names given in y (for the variables) & parms (for the parameters)
  with( as.list(c(y, parms)), {
    
    # Define change in S
    dS <- - beta * S * I
    
    # Define change in I
    dI <- + beta * S * I - gamma * I
    
    # Define change in R
    dR <- + gamma * I
    
    # Keep track of values at each time step
    return(list( c(dS, dI, dR) ))
  })
}


sicr <- function(t, y, parms) {
  
  # use the names given in y (for the variables) & parms (for the parameters)
  with( as.list(c(y, parms)), {
    
    # Define change in S (assume chronic not infectious)
    dS <- - beta * S * I
    
    # Define change in I
    dI <- + beta * S * I - gamma * I
    
    dC <- + c * gamma * I - 0.01 * gamma * C
    
    # Define change in R
    dR <- + (1 - c) * gamma * I + 0.01 * gamma * C
    
    # Keep track of values at each time step
    return(list( c(dS, dI, dC, dR) ))
  })
}


## Fn to simulate ODEs ---------------------------------------------------------

simulate_sir <- function(infectious_period, R0 = 2) {

  parameters <- c(beta = R0 / (N * infectious_period), gamma = 1/infectious_period)
  
  ode(y = c(S = N - 1, I = 1, R = 0), times = 1:150, func = sir, parms = parameters) %>% 
    data.frame() %>% mutate(infectious_period = infectious_period, R0 = R0)
}

simulate_sicr <- function(c, infectious_period = 5.9, R0 = 2) {
  
  parameters <- c(beta = R0 / (N * infectious_period), gamma = 1/infectious_period, c = c)
  
  ode(y = c(S = N - 1, I = 1, C = 0, R = 0), times = 1:150, func = sicr, parms = parameters) %>% 
    data.frame() %>% mutate(infectious_period = infectious_period, R0 = R0, c = c)
}


## Simulate for multiple parameter values --------------------------------------


N     <- 10000
R0vec <- seq(1.4, 2.6, by = 0.1)
c_vec <- seq(0, 1, by = 0.1)
infectious_period_vec <- seq(5.2, 6.4, by = 0.1)

sims   <- lapply(infectious_period_vec, simulate_sir) %>% bind_rows()

simsR0 <- lapply(R0vec, simulate_sir, infectious_period = 5.9) %>% bind_rows()

simsC  <- lapply(c_vec, simulate_sicr) %>% bind_rows()


## Plot models -----------------------------------------------------------------

## SIR - one value for infectious period ----
p_one <- 
  sims %>% filter(infectious_period == 5.9) %>% 
  gather(variable, value, S, I, R) %>% 
  ggplot(aes(x = time, y = value, colour = variable)) +
  geom_line(size = 1.5) +
  annotate(geom = "text", 
           x = c(120, 120, 120), 
           y = c(2550, 550, 8550),
           label = c("Susceptible", "Infectious", "Recovered"), size = 7,
           color = mycols3) +
  scale_colour_manual(guide = "none", values = mycols3, breaks = c("S", "I", "R")) +
  mytheme +
  labs(x = "Days", y = "Number of people")

ggsave(file = "figs/one.png", plot = p_one, dpi = 300, width = 8, height = 5)


## SIR - all values ----
p_all <- sims %>% 
  gather(variable, value, S, I, R) %>% 
  ggplot(aes(x = time, y = value, colour = variable, 
             group = interaction(variable, infectious_period))) +
  geom_line(size = 1.5, alpha = 0.5) +
  annotate(geom = "text", 
           x = c(120, 120, 120), 
           y = c(2550, 550, 8550),
           label = c("Susceptible", "Infectious", "Recovered"), size = 7,
           color = mycols3) +
  scale_colour_manual(guide = "none", values = mycols3, breaks = c("S", "I", "R")) +
  mytheme +
  labs(x = "Days", y = "Number of people")

ggsave(file = "figs/all.png", plot = p_all, dpi = 300, width = 8, height = 5)


## SIR - CIs ----
p_ci <- sims %>% 
  gather(variable, value, S, I, R) %>% 
  group_by(time, variable) %>%
  summarize(low  = quantile(value, probs = 0.025),
            med  = quantile(value, probs = 0.5  ),
            high = quantile(value, probs = 0.975) ) %>%
  ungroup() %>%
  ggplot(aes(x = time, y = med, colour = variable)) +
  geom_line(size = 1.5) +
  geom_ribbon(aes(ymin = low, ymax = high, fill = variable), alpha = 0.5, colour = NA) +
  annotate(geom = "text", 
           x = c(120, 120, 120), 
           y = c(2550, 550, 8550),
           label = c("Susceptible", "Infectious", "Recovered"), size = 7,
           color = mycols3) +
  scale_colour_manual(guide = "none", values = mycols3, breaks = c("S", "I", "R")) +
  scale_fill_manual(  guide = "none", values = mycols3, breaks = c("S", "I", "R")) +
  mytheme +
  labs(x = "Days", y = "Number of people")

ggsave(file = "figs/CIs.png", plot = p_ci, dpi = 300, width = 8, height = 5)


## SIR - CIs for R0 uncertainty ----
p_ciR0 <- simsR0 %>% 
  gather(variable, value, S, I, R) %>% 
  group_by(time, variable) %>%
  summarize(low  = quantile(value, probs = 0.025),
            med  = quantile(value, probs = 0.5  ),
            high = quantile(value, probs = 0.975) ) %>%
  ungroup() %>%
  ggplot(aes(x = time, y = med, colour = variable)) +
  geom_line(size = 1.5) +
  geom_ribbon(aes(ymin = low, ymax = high, fill = variable), alpha = 0.5, colour = NA) +
  scale_colour_manual(guide = "none", values = mycols3, breaks = c("S", "I", "R")) +
  scale_fill_manual(  guide = "none", values = mycols3, breaks = c("S", "I", "R")) +
  mytheme +
  labs(x = "Days", y = "Number of people")

ggsave(file = "figs/CIs_R0.png", plot = p_ciR0, dpi = 300, width = 8, height = 5)



## SICR - one value  ----
p_oneC <- 
  simsC %>% filter(c == 0.5) %>% 
  gather(variable, value, I, C, R) %>% 
  ggplot(aes(x = time, y = value, colour = variable)) +
  geom_line(size = 1.5) +
  annotate(geom = "text", 
           x = c(120, 120, 120), 
           y = c(550, 5250, 2750),
           label = c("Infectious", "Recovered", "Chronic Infected"), size = 7,
           color = c(mycols3[-1], "black")) +
  scale_colour_manual(guide = "none", values = c(mycols3[-1], "black"), breaks = c("I", "R", "C")) +
  mytheme +
  labs(x = "Days", y = "Number of people") +
  ylim(c(0, N))

ggsave(file = "figs/one_sicr.png", plot = p_oneC, dpi = 300, width = 8, height = 5)


