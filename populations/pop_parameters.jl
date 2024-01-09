using Pkg, DataFrames, DataFramesMeta, Distributions, Statistics, CSV, DynamicPipe, Dates

### Proportion of ages

ages = CSV.read("ages_cities.csv",DataFrame)

pop_counties = ages |>
@> groupby(:County) |>
@combine(:population = sum(:Total))

prop_ages = innerjoin(ages,pop_counties; on = :County) |>
@> @rtransform(:pop_age = :Male + :Female) |> @rtransform(:prop_age = :pop_age/:population) |>
@select(:County, :Age, :prop_age) |>
@rtransform(:activity_age = (if (:Age >= 6 && :Age <= 14) "student" elseif (:Age >= 15 && :Age <= 64) "worker" else "none" end)) |>
groupby([:County, :activity_age]) |>
@combine(:prob_activity = sum(:prop_age))

probs_ages = innerjoin(ages,pop_counties; on = :County) |>
@> @rtransform(:pop_age = :Male + :Female) |> @rtransform(:prop_age = :pop_age/:population) |>
@select(:County, :Age, :prop_age)

# Number of hours
hours_str = CSV.read("hrs_times.csv", DataFrame)


## entomological parameters

prop_str = DataFrame(Site = repeat(["Kisumu","Ukunda"],inner = 2), 
Environment = repeat(["Household","Non-Household"],2), prop = [0.623,0.527,0.538,0.326])

houses_cont = CSV.read("houses_containers.csv",DataFrame)
dist_c = DataFrame(Site = String[], Environment = String[], p_s = Float64[], shape_s = Float64[], p_m = Float64[], 
shape_m = Float64[], p_l = Float64[], shape_l = Float64[])
for r in 1:nrow(prop_str)
    all_n = houses_cont |> @> @rsubset(:Site == prop_str[r,1] && :Environment2 == prop_str[r,2])
    ns = all_n |> @> @rsubset(:n_s > 0) |> _.n_s
    nm = all_n |> @> @rsubset(:n_m > 0) |> _.n_m
    nl = all_n |> @> @rsubset(:n_l > 0) |> _.n_l
    prop_s = mean(all_n.n_s .> 0)
    prop_m = mean(all_n.n_m .> 0)
    prop_l = mean(all_n.n_l .> 0)
    d_s = fit_mle(Pareto, ns)
    d_m = fit_mle(Pareto, nm)
    d_l = fit_mle(Pareto, nl)
    transient_df = DataFrame(Site = prop_str[r,1], Environment = prop_str[r,2], p_s = prop_s, shape_s = shape(d_s), 
    p_m = prop_m, shape_m = shape(d_m), p_l = prop_l, shape_l = shape(d_l))
    dist_c = append!(dist_c,transient_df)
end

cont_prop = houses_cont |> @> groupby([:Site,:Environment2]) |> @combine(:ae_s = mean(:p_ae_s .> 0),
:ae_m = mean(:p_ae_m .> 0), :ae_l = mean(:p_ae_l .> 0),:ae_sm = mean(:p_ae_sm .> 0), 
:ae_sl = mean(:p_ae_sl .> 0), :ae_ml = mean(:p_ae_ml .> 0), :ae_sml = mean(:p_ae_sml .> 0)) |>
@select(:Site,:Environment = :Environment2,:total_ae = :ae_s+:ae_m+:ae_l+:ae_sm+:ae_sl+:ae_ml+:ae_sml) 

prop_containers = innerjoin(prop_str,cont_prop,dist_c,on=[:Site,:Environment])

dists = DataFrame(stage = repeat(["eggs","adult"], inner = 4),
Site = repeat(repeat(["Kisumu","Ukunda"],inner = 2),2), 
Environment = repeat(["Household","Non-Household"],4),
par1 = [3.30,3.40,3.53,3.61,1.356554,2.045432,2.452693,2.588756],
par2 = [0.672,0.645,0.797,0.760,0.4958515,0.4821921,0.6680835,0.9028559])

