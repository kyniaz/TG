using Distributions
using Plots
using DataFrames
using CSV
using RCall
R"library(survival)"

#Densidade gompertz, opção de log ou não
function dgompertz(x, a, b, log_opt = false)
    if(log_opt == true)
        return(log(a/b).+ log.(x./b) - a.*(expm1.(x./b)))        
    else
        return(a/b .*exp.(x./b).*exp.(- a.*(expm1.(x./b))))
    end
end

function pgompertz(x, a, b, log_opt = false)
    if(log_opt == true)
        return(log.(1 .- exp.(-a.* expm1.(x./b))))        
    else
        return(1 .- exp.(-a.*expm1.(x./b)))
    end
end

function sgompertz(x, a, b, log_opt = false)
    if(log_opt == true)
        return(-a.*expm1.(x./b))        
    else
        return(exp.(-a.*(expm1.(x./b))))
    end
end

function qgompertz(x, a, b)
    return(b.*log.(1 .- (1 ./a)*log.(1 .- x)))
end

function rgompertz(n, a, b)
    p = rand(Uniform(0, 1), n)
    return(b.*log.(1 .- (1 ./a)*log.(1 .- p)))
end

dados = DataFrame(CSV.File("dados_cancer_mama.csv"))

tempo = dados[!, :tempo]/365
cens = dados[!, :cens]

function log_veros_nomix(tempos, cens, mix, alfa, beta_p)      
    vals = sum(cens .* log(-log(mix))) + 
    sum(dgompertz(tempos, alfa, beta_p, true).* cens) + 
    sum(log(mix).*pgompertz(tempos, alfa, beta_p))
    return(vals)
end

function rprop_beta(ant)
    return (rand(Uniform(0, 1), 1))
end

function ldprop_beta(ant, prop)
    return 1
end

function rprop_gamma(ant)
    return (rand(Gamma(ant, 1), 1))
end

function ldprop_gamma(ant, prop)
    return (logpdf.(Gamma(ant, 1), prop))
end

gamma_prior_a = Gamma(1,1)
gamma_prior_b = Gamma(4,1)

function ldtgt1(t, d, prob, a, b)
    return(log_veros_nomix(t, d, prob, a, b) + logpdf(gamma_prior_a, a) + logpdf(gamma_prior_b, b))
end

function ldtgt2(t, d, prob, a, b)
    return(log_veros_nomix(t, d, prob, a, b) + logpdf(gamma_prior_a, a) + logpdf(gamma_prior_b, b))
end

function ldtgt3(t, d, prob, a, b)
    return(log_veros_nomix(t, d, prob, a, b) + logpdf(gamma_prior_a, a) + logpdf(gamma_prior_b, b))
end

function metropolis_gibbs(delta, prob, a, b, t, d, B = 10^4, start = 0.5) # rprop, ldprop,

    cadeia = ones(B)
    cadeia[1] = start

    if(delta[1] == 1.0)
        for i in 2:B
            prop = rprop_gamma(cadeia[i-1])[1]

            lratio = ldtgt1(t, d, prob, prop, b) - ldtgt1(t, d, prob, cadeia[i-1], b) +
            ldprop_gamma(prop,cadeia[i-1])-
            ldprop_gamma(cadeia[i-1],prop)

            if(log(rand(Uniform(0, 1), 1)[1]) <= lratio) 
                cadeia[i] = prop
            else
                cadeia[i] = cadeia[i-1]
            end
        end
    elseif(delta[1] == 2.0)
        for i in 2:B
        prop = rprop_gamma(cadeia[i-1])[1]

        lratio = ldtgt2(t, d, prob, a, prop) - ldtgt2(t, d, prob, a, cadeia[i-1]) +
          ldprop_gamma(prop,cadeia[i-1]) -
          ldprop_gamma(cadeia[i-1],prop)

        if(log(rand(Uniform(0, 1), 1)[1]) <= lratio) 
            cadeia[i] = prop
        else
            cadeia[i] = cadeia[i-1]
        end
    end
    elseif(delta[1] == 3.0)
        for i in 2:B
        prop = rprop_beta(cadeia[i-1])[1]

        lratio = ldtgt3(t, d, prop, a, b) - ldtgt3(t, d, cadeia[i-1], a, b) +
        ldprop_beta(prop,cadeia[i-1]) -
        ldprop_beta(cadeia[i-1],prop)

        if(log(rand(Uniform(0, 1), 1)[1]) <= lratio) 
            cadeia[i] = prop
        else
            cadeia[i] = cadeia[i-1]
        end
    end
    end
    return(last(cadeia))
end


function metropolis_within_gibbs(tempo, cens, B, BB, start)
    a = ones(BB)
    b = ones(BB)
    prob = ones(BB)

    a[1] = 0.15
    b[1] = 4
    prob[1] = 0.5
    
    for ii in 2:BB
        delta = trunc.(rand(Uniform(1, 4), 1))
        if (delta[1] == 1.0) 
            a[ii] = metropolis_gibbs(delta, prob[ii-1], a[ii-1], b[ii-1], tempo, cens , B, 0.15)[1]
            b[ii] = b[ii-1]
            prob[ii] = prob[ii-1]
        elseif (delta[1] == 2.0) 
            b[ii] = metropolis_gibbs(delta, prob[ii-1], a[ii-1], b[ii-1], tempo, cens , B, 4)[1]
            a[ii] = a[ii-1] 
            prob[ii] = prob[ii-1]
        elseif (delta[1] == 3.0)         
            prob[ii] = metropolis_gibbs(delta, prob[ii-1], a[ii-1], b[ii-1], tempo, cens , B, 0.5)[1]
            a[ii] = a[ii-1] 
            b[ii] = b[ii-1]
        end
    end
    return(Dict("a" => a, "b" => b, "prob" => prob))
end

@time begin
    modelo = metropolis_within_gibbs(tempo, cens, 400, 1000, 0.5)
end

a = modelo["a"]
b = modelo["b"]
prob = modelo["prob"]
vals_dens = [0:0.01:25;]
surv = median(prob).^(pgompertz(vals_dens, median(a), median(b)))

R"plot(survfit(Surv(tempo/365, cens) ~ 1, data = $dados), 
     xlab = 'Anos', 
     ylab = 'Sobrevivencia')"

R"lines($vals_dens, $surv, col = 'red')"

#DIC

dic = -2*log_veros_nomix(tempo, cens, median(prob), median(a), median(b))
print("DIC (modelo sem cura) %.2f:",dic)