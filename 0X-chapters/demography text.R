
We use eigenanalysis to solve \@ref(eq:eigsol) and give us the answers. Another way to find $\lambda_1$ is to simply iterate population growth a very large number of times, that is, let $t$ be very large. As $t$ grows, the annual growth rate, $N_{t+1}/N_t$, approaches $\lambda_1$. 



**Integral projection, or How I learned to stop worrying and love the infinite.**

Let's return to our example of garlic mustard. In this biennial species, rosette size is strongly predictive of survival, probability of flowering, and seed production. 

```{r Primula}
pm <- matrix(c(0.084, 0,0,14.676,
0.112, 0, 0, 21.974,
0, 0.226, 0.647, 0.668,
0, 0, 0.267, 0.297), nr=4, byrow=TRUE)
nms <- c("Seed bank", "Seedlings", "Vegetative", "Flowering")
colnames(pm) <- rownames(pm) <- nms

cm <- matrix(c(0.015,0,0,0,10.737,54.469,
0.030,0, 0,0,67.726,343.585,
0,0.076,0.108,0, 0, 0,
0,0,0.241,0.445, 0,0,
0, 0, 0.096, 0.274, 0,0,
0,0,0.008, 0.071,0,0), nr=6, byrow=TRUE)
nms <- c("Seed bank","Seedlings","Small vegetative","Large vegetative",
         "Small flowering","Large flowering")
colnames(cm) <- rownames(cm) <- nms
```

Consider the timber rattlesnake (*Crotalus horridus*). We could model its demography based on three life history stages: eggs, juveniles, and adults. These stages clearly have characteristic traits that cause survival, growtha nd fecundity to differ strongly among stages. However, at least two issues will plague our efforts. First, it is exceedingly hard to tell which females are reproductive and which are not. We are likely to have to rely on a size threshold where juveniles fecundity varies quite a bit among females, with some laying only a few eggs and other laying close to 20. Secondly, the size categories (0-95 cm)(> 95 cm). She found that on average, females produced about 10 eggs every 5 years. She found that, on average, about 60% of the eggs hatched. She marked individuals in a population, and found that each year about 80% of juveniles survived and of those, about 10% of the juveniles reach adulthood. Adults have a 70% annual survival rate. She assumed an even sex ratio. She censused snakes and eggs each year during egg laying season.

One way to visualize a demographic matrix is to make it 3-D, where the
vertical $z$-axis is the magnitude of the transition probability or
fecudnity. Let's show that for the mallwort example, with a slight
change: we set large adult reproduction $a_{1,3} = 2$.
\begin{equation}
\label{Lefk2}
\left(\begin{array}{ccc}
       0   & 0.5 & 2\\
       0.3 & 0  & 0\\
       0   & 0.5 & 0.9
       \end{array} \right)
\end{equation}

```{r,mat, fig=TRUE, echo=FALSE}
A <- matrix(c(0, 0.5, 2,
              0.3, 0, 0,
              0, 0.5, 0.9), nr=3, byrow=TRUE)

levelplot(c(A) ~ c(col(A)) + c(row(A)), ylim = c(3.5, 0.5), col.regions = gray(100:0/100),
          scales=list(alternating=3, at=1:3), ylab="Size at t = 1", xlab = "Size at t")

``` 

Gorchov and Christensen (2002)
\begin{equation}
\mathbf{A}=
  \left(\begin{array}{cccccc}
         0     &  0  & 0 & 0 & 0 & 1.642\\
         0.098 &  0   &  0 & 0 & 0 & 0.437  \\
         0     & 0.342& 0.591& 0.050& 0.095& 0\\
         0 &0.026& 0.295& 0.774& 0.177& 0.194\\
         0& 0& 0& 0.145& 0.596& 0.362\\
         0& 0& 0& 0.016 &0.277 &0.489\\
         \end{array} \right)
\end{equation}


