
using DataFrames, DataFramesMeta, StatsBase, DynamicPipe, Statistics, Random, Distributions


## This function reflects the changes on population parameters during calibration process
function cont_vect(city::String, environment::String)
    contain = prop_containers |> @> @rsubset(:Site == city, :Environment == environment)
    positivity_c = rand(Binomial(1,((environment == "Household") ? ((city == "Kisumu") ? contain.prop[1]-0.22 : contain.prop[1]+0.3) : (city == "Kisumu") ? contain.prop[1]-0.22 : contain.prop[1]+0.3)),1)[1]
    if positivity_c > 0
        small = rand(Binomial(1,((city == "Kisumu") ? contain.p_s[1] : (environment == "Household") ? contain.p_s[1]+0.15 : contain.p_s[1]+0.1)),1)[1]
        medium = rand(Binomial(1,((city == "Kisumu") ? contain.p_m[1] : (environment == "Household") ? contain.p_m[1]-0.05 : contain.p_m[1]-0.05)),1)[1]
        large = rand(Binomial(1,((city == "Kisumu") ? contain.p_l[1] : (environment == "Household") ? contain.p_l[1]-0.05 : contain.p_l[1]-0.05)),1)[1]
        n_small = (small == 0) ? 0 : round(Integer,rand(truncated(Pareto(contain.shape_s[1],1), upper = 10),1)[1])
        n_medium = (medium == 0) ? 0 : round(Integer,rand(truncated(Pareto(contain.shape_m[1],1), upper = 10),1)[1])
        n_large = (large == 0) ? 0 : round(Integer,rand(truncated(Pareto(contain.shape_l[1],1), upper = 3),1)[1])
        liters = if small+medium+large > 0
            if (medium > 0 || large > 0)
                10
            elseif small > 0
                sum(rand(truncated(Poisson(1),lower=1),sum(n_small)))
            end
        else
            0
        end
        positivity_ae = rand(Binomial(n_small+n_medium+n_large, (city == "Kisumu") ? contain.total_ae[1]-0.1 : contain.total_ae[1]-0.035),1)[1]
        females = if positivity_ae > 0
            params_v = dists |> @> @rsubset(:Site == city && :Environment == environment)
            p_eggs = @>> params_v |> @rsubset(:stage == "eggs")
            p_adult = @>> params_v |> @rsubset(:stage == "adult")
            females_sim = Vector{Integer}()
            for k in 1:positivity_ae
                female_q = Vector{Float64}()
                eggs_n = round.(Integer,rand(LogNormal(p_eggs.par1[1],p_eggs.par2[1]),k))
                for f in 1:length(eggs_n)
                    n_egg_p = cdf(LogNormal(p_eggs.par1[1],p_eggs.par2[1]),eggs_n[f])
                    female_q2 = quantile(LogNormal(p_adult.par1[1],p_adult.par2[1]),n_egg_p)
                    push!(female_q,female_q2)
                end
                push!(females_sim,round(Integer, sum(female_q)))
            end
            sum(females_sim)
        else
            0
        end
        transient_df = DataFrame(Container_s = n_small, Container_ml = n_medium+n_large, liters = liters,
        n_vectors = females)
    else
        transient_df = DataFrame(Container_s = 0, Container_ml = 0, liters = 0, n_vectors = 0)
    end
    return(transient_df)
end

function u_rate(t)
    lf = -1.48e-01*(t-9.16)*(t-37.73)
    return(1/lf)
end

function a_rate(t::Float64)
    a = 2.02e-04*t*(t-13.35)*(40.08-t)^0.5
    return(a)
end

function vector_comp(t::Float64)
    b = 8.49e-04*t*(t-17.05)*(35.83-t)^0.5
    c = 4.91e-04*t*(t-12.22)*(37.46-t)^0.5
    return(b*c)
end

function mosq_prev(a::Float64,b::Float64,Nm::Integer,Nh::Integer,Ih::Integer,um::Float64)
    r = a/Nh
    vm = Nm*um
    prev = (r*b*vm*Ih)/((r*b*Ih+um)*um)
    return(prev)
end

function hrs_str(city::String,str::String)
    times = @rsubset hours_str :Site == city && :spaces == str
    lapse = round(rand(Uniform(times.h_initial[1],times.h_final[1]),1)[1], digits = 1)
    return(lapse)
end

