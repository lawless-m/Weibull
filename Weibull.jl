#==

uses JuMP instead of Goal seeker and Solver as described in

https://www.real-statistics.com/distribution-fitting/distribution-fitting-via-maximum-likelihood/fitting-weibull-parameters-mle/

==#

using Roots # Add.pkg("Roots")
using JuMP
using Ipopt

function simplefit(xbar, svar)
    model = Model(solver=IpoptSolver())
    @variable(model, n, start = 1)
    @constraint(model, n >= 1)
    @NLobjective(model, Min, (log(gamma(1+2/n)) - 2log(gamma(1+1/n)) - log(xbar^2+svar^2) + 2log(xbar))^2)
    TT = STDOUT # save original STDOUT stream
    redirect_stdout()
    solve(model)
    redirect_stdout(TT) # restore STDOUT
    beta = getvalue(n)
    alpha = xbar / exp(log(gamma(1+1/beta)))
    alpha, beta
end

function fit(samples)
    xbar = mean(samples)
    svar = std(samples)
    a,b = simplefit(xbar, svar)
    model = Model(solver=IpoptSolver())
    n = length(samples)
    @variable(model, alpha, start = a)
    @variable(model, beta, start = b)
    @constraint(model, alpha >= 0)
    @constraint(model, beta >= 0)
    @NLobjective(model, Max, n * (log(beta) - beta * log(alpha)) + sum((beta - 1) * log(samples[s]) - (samples[s] / alpha) ^ beta for s in 1:length(samples)))

    TT = STDOUT # save original STDOUT stream
    redirect_stdout()
    solve(model)
    redirect_stdout(TT) # restore STDOUT

    alpha = getvalue(alpha)
    beta = getvalue(beta)

    alpha, beta
end
