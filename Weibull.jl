#==

uses JuMP instead of Goal seeker and Solver as described in

https://www.real-statistics.com/distribution-fitting/distribution-fitting-via-maximum-likelihood/fitting-weibull-parameters-mle/

==#
module Weibull

using JuMP
using Ipopt

function simplefit(xbar, svar)
    model = Model(solver=IpoptSolver())
    @variable(model, b, start = 1)
    @constraint(model, b >= 1)
    #@NLobjective(model, Min, (log(gamma(1+2/b)) - 2log(gamma(1+1/b)) - log(xbar^2+svar^2) + 2log(xbar))^2)
    @NLobjective(model, Min, (gamma(1+2/b) - 2gamma(1+1/b) - xbar^2+svar^2 + 2xbar)^2)
    TT = STDOUT # save original STDOUT stream
    redirect_stdout()
    solve(model)
    redirect_stdout(TT) # restore STDOUT
    beta = getvalue(b)
    #alpha = xbar / exp(log(gamma(1+1/beta)))
    alpha = xbar / gamma(1+1/beta)
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

function test()
    a, b = fit([23,19,37,38,40,36,172,48,113,90,54,104,90,54,157,51,77,78,144,34,29,45,16,15,37,218,170,44,121])
    println("alpha $a, beta $b")

end

################
end # module
