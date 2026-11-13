# Surprisal, Shannon entropy, and Hill numbers as means of rarity ----
# A Primer of Ecology with R (2nd ed.) -- demo for the diversity chapter
# Each individual tree carries a surprisal, -ln(p_i).
# Shannon H is the average surprisal per tree, and every Hill number
# is a weighted power mean of the per-tree rarity, 1/p_i.

## Setup ----
library(tidyverse)
library(vegan)   # for the BCI data
data(BCI)        # Barro Colorado Island, 50 one-ha subplots, 225 species
set.seed(42)

## Build two communities ----
# BCI: all 50 ha pooled
bci <- tibble(species = names(BCI), n = colSums(BCI)) |>
  filter(n > 0) |>
  mutate(p = n / sum(n))

# Hypothetical contrast: a geometric series (strong dominance) with
# the same number of individuals as BCI; species expected to have
# fewer than one tree are dropped, so richness is lower
S  <- nrow(bci)
N  <- sum(bci$n)
k  <- 0.05                                # niche pre-emption fraction
geo <- tibble(species = paste0("sp", seq_len(S)),
              p = k * (1 - k)^(seq_len(S) - 1)) |>
  mutate(n = round(N * p / sum(p))) |>
  filter(n > 0) |>
  mutate(p = n / sum(n))

comm <- bind_rows(BCI = bci, `Geometric series` = geo, .id = "community") |>
  mutate(surprisal = -log(p))             # nats per tree

## Hill numbers as weighted power means of rarity ----
# For q != 1, qD = ( sum p * (1/p)^(1-q) )^(1/(1-q)): power mean of order 1-q
# For q == 1, the limit is the weighted geometric mean, exp(H)
hill <- function(p, q) {
  if (q == 1) exp(-sum(p * log(p))) else sum(p^q)^(1 / (1 - q))
}

hill_tab <- comm |>
  summarise(H        = sum(p * surprisal),        # mean surprisal (nats)
            richness = hill(p, 0),                # arithmetic mean of 1/p
            expH     = hill(p, 1),                # geometric mean of 1/p
            invSimp  = hill(p, 2),                # harmonic mean of 1/p
            .by = community)
hill_tab

# Check: vegan agrees with the mean-surprisal definition of H
c(bci = diversity(colSums(BCI), "shannon"),
  by_hand = hill_tab$H[hill_tab$community == "BCI"])

## Figure 1: distribution of surprisal across individual trees ----
# Each species contributes n trees, each with the same surprisal.
# Weighting by n turns the species list into a tree-level distribution.
ggplot(comm, aes(surprisal, weight = n)) +
  geom_histogram(binwidth = 0.5, boundary = 0,
                 fill = "grey70", colour = "white") +
  geom_vline(data = hill_tab, aes(xintercept = H),
             linewidth = 1, colour = "firebrick") +
  geom_text(data = hill_tab,
            aes(x = H, y = Inf,
                label = sprintf("H = %.2f nats\nexp(H) = %.1f spp.", H, expH)),
            hjust = -0.05, vjust = 1.3, colour = "firebrick",
            inherit.aes = FALSE) +
  facet_wrap(~ community, ncol = 1, scales = "free_y") +
  labs(x = expression("Surprisal of one tree, " * -ln~p[i] * " (nats)"),
       y = "Number of trees") +
  theme_bw()

## Figure 2: one tree at a time -- the running mean converges on H ----
# Draw trees at random from BCI and track the surprisal of each.
# Individual trees vary a lot; their running mean settles on H.
n_draws <- 2000
draws <- tibble(tree = seq_len(n_draws),
                species = sample(bci$species, n_draws,
                                 replace = TRUE, prob = bci$p)) |>
  left_join(select(bci, species, p), by = "species") |>
  mutate(surprisal    = -log(p),
         running_mean = cumsum(surprisal) / tree)

H_bci <- hill_tab$H[hill_tab$community == "BCI"]

ggplot(draws, aes(tree)) +
  geom_point(aes(y = surprisal), alpha = 0.15, size = 0.8) +
  geom_line(aes(y = running_mean), colour = "firebrick", linewidth = 0.8) +
  geom_hline(yintercept = H_bci, linetype = 2) +
  scale_x_log10() +
  labs(x = "Trees observed (log scale)",
       y = expression("Surprisal, " * -ln~p[i] * " (nats)")) +
  theme_bw()

## Figure 3: the Hill profile as a continuum of means ----
# As q increases, the average gives rare (high-surprisal) trees less weight.
qs <- seq(0, 3, by = 0.05)
profile <- comm |>
  reframe(q = qs, D = map_dbl(qs, \(q) hill(p, q)), .by = community)

ggplot(profile, aes(q, D, colour = community)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = c(0, 1, 2), linetype = 3) +
  annotate("text", x = c(0, 1, 2), y = Inf, vjust = 1.5, hjust = -0.05,
           label = c("richness\n(arithmetic)", "exp(H)\n(geometric)",
                     "1/Simpson\n(harmonic)"), size = 3) +
  scale_y_log10() +
  labs(x = "Order, q", y = "Effective number of species, qD (log scale)",
       colour = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
