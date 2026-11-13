## ----------------------------------------------------------------------
## Two simple two-stage (Lefkovitch) matrix population models
## One plant, one animal -- both use the same A = [[p11, F12],[p21, p22]]
## convention as the Western toad example in 04-DIDemography.qmd:
##   columns = stage at t, rows = stage at t+1
##   p11 = stage-1 stasis (survive, don't grow up)
##   p21 = stage-1 -> stage-2 (survive AND grow up / mature)
##   p22 = stage-2 stasis (adult survival)
##   F12 = stage-2 -> stage-1 (fecundity contribution of an average adult)
## Both are pre-breeding census models, no regression (g12 = 0), so they
## can be dropped straight into @sec-transition-matrix as a second worked
## example alongside/after the toad, or used standalone as an exercise.
## ----------------------------------------------------------------------


## =========================================================================
## 1. PLANT EXAMPLE: a long-lived perennial herb
##    Stage 1 = vegetative (non-flowering) rosette/juvenile
##    Stage 2 = flowering adult
## -------------------------------------------------------------------------
## Parameter values below are typical/illustrative for this life-history
## type (slow-growing, long-lived perennial: high adult survival, low
## per-capita fecundity, most individuals spend several years vegetative
## before first flowering) rather than a single literature source -- swap
## in real numbers from Silene, Trillium, Lomatium, an orchid, etc. if you
## want to tie it to a specific species. Unlike the garlic mustard example
## later in the chapter, there is no seed bank stage here -- this is the
## "simple" 2-stage case.

# symbolic layout:
#      veg.     adult
# veg  p11      F12
# adu  p21      p22
#
# p11 = probability a vegetative plant survives AND remains vegetative
# p21 = probability a vegetative plant survives AND matures to flowering
# p22 = probability a flowering adult survives (and remains reproductive)
# F12 = expected number of new vegetative recruits (seeds -> seedlings
#       surviving to the census) contributed by an average flowering adult

p11 <- 0.55   # vegetative stasis
p21 <- 0.15   # vegetative -> flowering (maturation)
p22 <- 0.85   # adult survival
F12 <- 0.90   # adult fecundity contribution

(A.plant <- matrix(c(p11, F12,
                     p21, p22), nrow = 2, byrow = TRUE,
                    dimnames = list(c("Vegetative", "Flowering"),
                                    c("Vegetative", "Flowering"))))

# --- eigenanalysis, following the same steps as @sec-eigenanalysis ---
eA.plant   <- eigen(A.plant)
dom.pos    <- which.max(Mod(eA.plant[["values"]]))
(lambda.plant <- Re(eA.plant[["values"]][dom.pos]))            # finite rate of increase
w1.plant   <- Re(eA.plant[["vectors"]][, dom.pos])
(ssd.plant <- w1.plant / sum(w1.plant))                        # stable stage distribution

# reproductive value: dominant left eigenvector (eigen of the transpose)
etA.plant  <- eigen(t(A.plant))
dom.pos2   <- which.max(Mod(etA.plant$values))
v1.plant   <- Re(etA.plant$vectors[, dom.pos2])
(rv.plant  <- v1.plant / v1.plant[1])

cat(sprintf("\nPLANT: lambda = %.3f  |  stable stage dist. (veg, adult) = (%.2f, %.2f)  |  RV (veg, adult) = (%.1f, %.1f)\n",
            lambda.plant, ssd.plant[1], ssd.plant[2], rv.plant[1], rv.plant[2]))

# lambda ~ 1.10: population growing ~10%/yr under these rates
# stable stage dist. ~62% vegetative / 38% flowering -- most individuals
#   are "stuck" pre-reproductive at any given time, typical of slow-maturing
#   perennials
# RV: a flowering adult is worth ~3.6x a vegetative individual, because it
#   is already contributing to reproduction and has a decent chance of
#   surviving to do so again


## =========================================================================
## 2. ANIMAL EXAMPLE: yellow-bellied marmot (Marmota flaviventris)
##    Stage 1 = juvenile (pups through their first year)
##    Stage 2 = breeding adult (>= 1 yr old)
## -------------------------------------------------------------------------
## Survival values are real published rates from the RMBL long-term
## population (Schwartz, Armitage & Van Vuren 1998, J. Zoology, 32-yr
## demography); fecundity uses the reported mean litter size, simplified
## to a single adult stage. See conversation sources below.
##
##   juvenile survival (age 0->1, females): 0.543
##   adult survival    (age 1->2, females): 0.495 (rounded to 0.50 here)
##   mean litter size: 4.1 pups
##
## Real marmot demography is usually built with 3+ stages (juvenile,
## yearling, adult) because females don't start breeding until ~age 2 and
## most 2-yr-olds skip breeding -- collapsing that into one "adult" class
## is the simplification that makes this a 2-stage teaching example rather
## than a publication-ready matrix. F12 below bundles: the probability an
## average adult female breeds in a given year (~0.4, reflecting that many
## adults skip years), litter size (4.1), and a 1:1 sex ratio at birth --
## i.e. exactly the same black-box fecundity trick used for F12 in the
## toad matrix (rho * phi * survival terms all multiplied together).

# symbolic layout:
#       juv.    adult
# juv   p11     F12
# adu   p21     p22
#
# p11 = 0 (no juvenile-juvenile stasis in this simplified model: a juvenile
#          either dies or survives into the adult stage within a year)
# p21 = juvenile survival (survive AND enter the adult stage)
# p22 = adult annual survival
# F12 = rho (annual breeding probability) * phi (litter size) * sex ratio

rho <- 0.40    # proportion of adult females breeding in an average year
phi <- 4.1     # mean litter size (Schwartz et al. 1998)
sr  <- 0.5     # proportion female at birth
rho*phi*sr

p11 <- 0
p21 <- 0.543   # juvenile survival
p22 <- 0.50    # adult survival
F12 <- rho * phi * sr # 0.82

(A.marmot <- matrix(c(p11, F12,
                      p21, p22), nrow = 2, byrow = TRUE,
                     dimnames = list(c("Juvenile", "Adult"),
                                     c("Juvenile", "Adult"))))

eA.marmot  <- eigen(A.marmot)
dom.pos    <- which.max(Mod(eA.marmot[["values"]]))
(lambda.marmot <- Re(eA.marmot[["values"]][dom.pos]))
w1.marmot  <- Re(eA.marmot[["vectors"]][, dom.pos])
(ssd.marmot <- w1.marmot / sum(w1.marmot))

etA.marmot <- eigen(t(A.marmot))
dom.pos2   <- which.max(Mod(etA.marmot$values))
v1.marmot  <- Re(etA.marmot$vectors[, dom.pos2])
(rv.marmot <- v1.marmot / v1.marmot[1])

cat(sprintf("\nMARMOT: lambda = %.3f  |  stable stage dist. (juv, adult) = (%.2f, %.2f)  |  RV (juv, adult) = (%.1f, %.1f)\n",
            lambda.marmot, ssd.marmot[1], ssd.marmot[2], rv.marmot[1], rv.marmot[2]))

# lambda ~ 0.96: population slightly declining under these rates -- this
#   direction (Ro < 1) matches Schwartz et al.'s own Ro estimates of
#   0.67-0.85 for this population, so the simplified 2-stage model
#   reproduces the qualitative finding even though it collapses out the
#   yearling stage
# stable stage dist. ~46% juvenile / 54% adult
# RV: an adult is worth ~1.8x a juvenile -- a smaller gap than the plant,
#   because juvenile survival is already fairly high (0.54) so a juvenile
#   isn't "risking" as much of its future reproduction as the plant's
#   vegetative stage is

## -------------------------------------------------------------------------
## Sources for the marmot parameters:
##  Schwartz, Armitage & Van Vuren 1998. A 32-year demography of
##    yellow-bellied marmots (Marmota flaviventris). J. Zoology 246:337-346.
##    https://sites.lifesci.ucla.edu/eeb-rmbl-marmots/wp-content/uploads/sites/190/2020/03/Schwartz_etal_1998_JZL.pdf
## -------------------------------------------------------------------------
