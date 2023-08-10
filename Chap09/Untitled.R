######################################

## functions for responses
f1 <- function(x, s, a, b, c) {s - a*x / (b + c*x)}
f1b <- function(x, s, a, b, c) {s - b*x - a*x^c}
f2 <- function(x, a, b, c) { a / (b + c*x) }

## differential responses
diff.resp <- ggplot() +
  lims(x=c(0,2)) +
  geom_function(
    fun = f1b, 
    args = list(s=2, a=1, b=0, c=2),
    xlim=c(0,2), 
    n=100, linetype=2, color=2
  )  +
  annotate("text", x=1, y=1.25, label="Sp. b") +
  geom_function(
    fun = f1, 
    args = list(s=2.3, a=3, b=.5, c=.5),
    xlim=c(0,2), 
    n=100,
    linetype=1, color=1
  ) +
  annotate("text", x=.3, y=.4, label="Sp. a") +
  labs(y="Growth rate (r)",
       x="Competition, and resource availability") +
  theme_classic() + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position='none')

#############################
### three graphs for two environments

# non-additivity
additive <- 
  ggplot() + lims(x=c(0,2)) +
  geom_function(
    fun = f1, 
    args = list(s=4.3, a=3, b=.5, c=.5),
    xlim=c(0,2), 
    n=100,
    linetype=1, color=1
  )  +
  geom_function(
    fun = f1, 
    args = list(s=1.3, a=3, b=.5, c=.5),
    xlim=c(0,2), 
    n=100,
    linetype=1, color=1
  ) +
  labs(x="Competition", y="Growth rate (r)") +
  annotate("text", x=.5, y=0.5, label="Poor environment") +
  annotate("text", x=.5, y=4, label="Good environment") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position='none')

# subadditivity
sub <- ggplot() + lims(x=c(0,2)) +
  geom_function(
    fun = f2, 
    args = list(a=1, b=1, c=1),
    xlim=c(0,2), 
    n=100
  )   +
  geom_function(
    fun = f2, 
    args = list(a=5, b=1, c=1.5),
    xlim=c(0,2), 
    n=100
  ) + 
  labs(x="Competition", y="Growth rate (r)") +
  annotate("text", x=.5, y=1, label="Poor environment") +
  annotate("text", x=.5, y=2, label="Good environment") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position='none')

# superadditivity
super <- ggplot() + lims(x=c(0,2)) +
  geom_function(
    fun = f2, 
    args = list(a=4, b=1, c=1.5),
    xlim=c(0,2), 
    n=100
  )   +
  annotate("text", x=.75, y=2.5, label="Poor environment") +
  geom_function(
    fun = f2, 
    args = list(a=5, b=1, c=.1),
    xlim=c(0,2), 
    n=100
  ) + 
  labs(x="Competition", y="Growth rate (r)") +
  annotate("text", x=0.75, y=5.25, label="Good environment")  +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position='none')

blank <- ggplot() + lims(x=c(0,2)) +
  theme_minimal() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position='none')

(diff.resp + additive ) / (sub + super)
