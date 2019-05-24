#==

https://stats.stackexchange.com/questions/19866/how-to-fit-a-weibull-distribution-to-input-data-containing-zeroes


implements the Blischke-Scheuer method-of-moments estimation of (a,b)
for the Weibull distribution F(t) = 1 - exp(-(t/a)^b)

==#

using Roots # Add.pkg("Roots")
using JuMP
using Ipopt

gln(n) = log(gamma(n))

function simplefit(xbar, svar)
    model = Model(solver=IpoptSolver())
    @variable(model, n, start = 1)
    @constraint(model, n >= 1)
    @NLobjective(model, Min, (log(gamma(1+2/n)) - 2log(gamma(1+1/n)) - log(xbar*xbar+svar*svar) + 2log(xbar))^2)
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

samples = [23,19,37,38,40,36,172,48,113,90,54,104,90,54,157,51,77,78,144,34,29,45,16,15,37,218,170,44,121]
samples = [509, 660, 386, 753, 811, 613, 848, 725, 315, 872, 487, 512]
a, b = fit(samples)

println("a $a, b $b")
