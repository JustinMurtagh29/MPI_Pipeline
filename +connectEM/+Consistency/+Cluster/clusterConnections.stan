/**
 * We have `N` same-axon same-dendrite spine synapse pairs (see Bartol
 * et al. 2015 eLife). For each spine synapse we've measured the area of
 * the axon-spine interface (ASI) in µm². The values `log10Asi1[i]` and
 * `log10Asi2[i]` contain the log10(ASI area [µm²]) of synapse pair `i`.
 * The ordering within a synapse pair is of no importance.
 *
 * Let's now assume that there exist `K` different types of connections
 * (e.g., induced by the molecular, morphological, electrophysiological
 * state of the pre- and / or postsynaptic neuron). The probability that
 * a connection is of type `k` is given by `theta[k]`. Importantly,
 * spine synapses from the same axon onto the same dendrite are of the
 * same type!
 *
 * Each connection type has its own ASI area distribution. The log10(ASI
 * areas) of connections of type `k` is assumed to follow a normal
 * distribution with mean `mu[k]` and standard deviation `sigma[k]`.
 *
 * The likelihood of spine synapse pair `i` is thus given by:
 *
 *   theta[1] * normal_pdf(log10Asi1[i] | mu[1], sigma[1]
 *            * normal_pdf(log10Asi2[i] | mu[1], sigma[1])
 * + theta[2] * normal_pdf(log10Asi1[i] | mu[2], sigma[2])
 *            * normal_pdf(log10Asi2[i] | mu[2], sigma[2])
 * + ...
 * + theta[K] * normal_pdf(log10Asi1[i] | mu[K], sigma[K])
 *            * normal_pdf(log10Asi2[i] | mu[K], sigma[K])
 *
 * Stan is doing this whole spiel in log-space. The above sum becomes a problem
 * in log-space, so have to use the `log_sum_exp` function. It i) converts the
 * terms of the sum from log to linear space, ii) sums up the likelihood
 * values, and finally iii) converts the total likelihood back to log-space.
 *
 * Some comments about the priors: They were chosen as normal
 * distributions because I don't know any better (clearly negative
 * standard deviations do not make any sense, but this is handled by the
 * lower bound specified in the parameter section). The mean values of
 * the priors were roughly set to the empirically determined values. The
 * standard deviations are chosen large enough not make the priors only
 * weakly informative.
 *
 * Written by
 *   Alessandro Motta <alessandro.motta@brain.mpg.de>
 */
data {
    int<lower = 0> N;
    int<lower = 1> K;
    vector[N] log10Asi1;
    vector[N] log10Asi2;
}

parameters {
    ordered[K] mu;
    real<lower = 0> sigma[K];
    simplex[K] theta;
}

model {
    mu ~ normal(-0.5, 2);
    sigma ~ uniform(0, 4);

    for (n in 1:N) {
        real ps[K];
        for (k in 1:K) {
            ps[k] = log(theta[k])
                  + normal_lpdf(log10Asi1[n] | mu[k], sigma[k])
                  + normal_lpdf(log10Asi2[n] | mu[k], sigma[k]);
        }

        target += log_sum_exp(ps);
    }
}
