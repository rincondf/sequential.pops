stbp_posterior_composite(data = c(1, 2, 3),
                         hypothesis = 2,
                         likelihood_func = function(data, x) {
                           dpois(data, lambda = x)
                         },
                         prior = 0.5,
                         lower_bnd = 0,
                         upper_bnd = Inf)

stbp_posterior_composite(data = c(1, 2, 3), hypothesis = 2,
                         likelihood_func = function(data, x) {dnbinom(data, size = 2, mu = x)},
                         prior = 0.5,
                         lower_bnd = 0,
                         upper_bnd = Inf)

stbp_posterior_simple(data = c(0, 0, 0),
                      hypothesis = 0,
                      likelihood_func = function(data, x)
                        {dpois(data, lambda = x)},
                      prior = 0.5,
                      upper_bnd = Inf)


stbp_simple(data = matrix(rep(0, 30), 10, 3),
            hypothesis = 0,
            likelihood_func= function(data, x)
              {dpois(data, lambda = x)},
            prior = 0.5,
            upper_bnd = Inf,
            lower_criterion = 0,
            upper_criterion = 0.9999)


stbp_simple(data = matrix(rep(0, 90), 30, 3),
            hypothesis = 0,
            likelihood_func= function(data, x)
            {dpois(data, lambda = x)},
            prior = 0.5,
            upper_bnd = Inf,
            lower_criterion = 0,
            upper_criterion = 0.9999)

set.seed(101)
stbp_composite(data = rpois(5, lambda = 3),
               hypothesis = 5,
               likelihood_func = function(data, x)
               {dpois(data, lambda = x)},
               prior = 0.5,
               lower_bnd = 0,
               upper_bnd = Inf,
               lower_criterion = 0.001,
               upper_criterion = 0.999)

H <- c(2, 5, 10, 20, 40, 40, 20, 10, 5, 2)

countP <- matrix(NA, 3, 10)

set.seed(101)
for(i in 1:10){
  countP[, i] <- rpois(3, lambda = (H[i] - 1))
}

stbp_composite(data = countP,
               hypothesis = H,
               likelihood_func = function(data, x)
               {dpois(data, lambda = x)},
               prior = 0.5,
               lower_bnd = 0,
               upper_bnd = Inf,
               lower_criterion = 0.001,
               upper_criterion = 0.999)


