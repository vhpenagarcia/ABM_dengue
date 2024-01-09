using DataFrames, DataFramesMeta, Distributions, Statistics, Gadfly, DynamicPipe, StatsBase, CategoricalArrays

include("functions_pops.jl")

### Population

houses_ukunda = convert(Integer,round((77686/7.3)*0.26)) # Creating a population of about 20,000 people
houses_kisumu = convert(Integer,round((397957/4.6)*0.051))


houses = DataFrame(city = ["Ukunda","Kisumu"], houses_n = [houses_ukunda,houses_kisumu])
houses = houses |> @> @rtransform(:markets = round(Integer,:houses_n/30),:churches = round(Integer,:houses_n/50), :social = round(Integer,:houses_n/30)) |> 
@rtransform(:total = :markets+:churches+:social)
population = DataFrame(city = String[],house = String[],subj_cod = String[],age = Integer[],
student = Integer[],worker = Integer[],none = Integer[])

for i in 1:nrow(houses), j in 1:houses.houses_n[i]
    inhabitants = if houses.city[i] == "Ukunda"
        rand(Poisson(7.2990654), 1)
    else houses.city[i] == "Kisumu"
        rand(Poisson(4.6129032), 1)
    end
    cod_h = string(first(houses.city[i],1),j)
    p_ages = probs_ages[probs_ages.County .== houses.city[i],:]
    d_age = DiscreteNonParametric(p_ages.Age,p_ages.prop_age)
    for k in 1:inhabitants[1]
        age = rand(d_age,1)[1]
        occupation = if age <= 3
            "none"
            elseif age > 3 && age <= 5
            (rand(Binomial(1,0.44),1)[1]) == 1 ? "student" : "none"
            elseif age > 5 && age <= 14
            (rand(Binomial(1,0.99),1)[1]) == 1 ? "student" : "none"
            elseif age > 14 && age <= 17
            (rand(Binomial(1,0.71),1)[1]) == 1 ? "student" : (rand(Binomial(1,0.637),1)[1]) == 1 ? "worker" : "none"
            elseif age > 17 && age <= 64
            (rand(Binomial(1,0.637),1)[1]) == 1 ? "worker" : "none"
            elseif age > 64
            "none"
        end
        transient_df = DataFrame(city = houses.city[i], house = cod_h,
        subj_cod = string(cod_h,"_",k), age = age,
        student = (occupation == "student") ? 1 : 0,
        worker = (occupation == "worker") ? 1 : 0,
        none = (occupation == "none") ? 1 : 0)
        append!(population,transient_df)
    end
end


pop_sum = population |> @> groupby(:city) |> @combine(:students = sum(:student), :workers = sum(:worker),
:nones = sum(:none)) |> @rtransform(:pop = :students + :workers + :nones) |> @select(:city, :pop)

### Assign workplace and school

stats_population = population |> @> groupby(:city) |> @combine(:students = sum(:student), :workers = sum(:worker),
:nones = sum(:none)) |> @rtransform(:schools = round(Integer,:students/360), 
:workplaces = round(Integer,:workers/19)) |> @select(:city,:schools,:workplaces)


places = DataFrame(id_place = String[])

for p in 1:nrow(population)
    city_p = population.city[p]
    city_stat = stats_population |> @> @rsubset(:city == city_p)
    place = if population.student[p] == 1
        school = rand(DiscreteUniform(1,city_stat.schools[1]),1)
        string(first(city_p[1],1)[1],school[1],"_s")
    elseif population.worker[p] == 1
        workplace = rand(DiscreteUniform(1,city_stat.workplaces[1]),1)
        string(first(city_p[1],1)[1],workplace[1],"_w")
    else
        population.house[p]
    end
    place_df = DataFrame(; id_place = place)
    places = append!(places,place_df)
end

population = hcat(population,places)

### Infections status at day 1 with probability of 0.0008

p_inf = 0.0008

population = population |>
@> @rtransform(:susc = rand(Binomial(1,1-p_inf),1)[1]) |>
@rtransform(:latent = (:susc == 1) ? 0 : rand(Binomial(1,5/12),1)[1]) |>
@rtransform(:inf = (:susc == 0 && :latent == 0) ? 1 : 0, :recov = 0) |> 
@rtransform(:lat_day = (:latent == 1) ? rand(DiscreteUniform(1,5),1)[1] : 0,
:inf_day = (:inf == 1) ? rand(DiscreteUniform(1,6),1)[1] : 0, :imm = 0,
:day_recov = 0)

population |> @> groupby(:city) |> @combine(:students = sum(:student), :workers = sum(:worker), :nones = sum(:none),
:susc = sum(:susc), :latent = sum(:latent), :infs = sum(:inf), :recovs = sum(:recov))


### Mosquito population

structures = DataFrame(city = String[], type = String[], code = String[])


for s in 1:nrow(houses)
    city_st = houses.city[s]
    for b in 1:houses.houses_n[s]
        cod_house = string(first(city_st,1)[1],b)
        transient_df1 = DataFrame(city = city_st, type = "house", code = cod_house)
        structures = append!(structures,transient_df1)
    end
    places_st = stats_population |> @> @rsubset(:city == city_st)
    for c in 1:places_st.schools[1]
        cod_school = string(first(city_st,1)[1],c,"_s")
        transient_df2 = DataFrame(city = city_st, type = "school", code = cod_school)
        structures = append!(structures,transient_df2)
    end
    for w in 1:places_st.workplaces[1]
        cod_wp = string(first(city_st,1)[1],w,"_w")
        transient_df3 = DataFrame(city = city_st, type = "workplace", code = cod_wp)
        structures = append!(structures,transient_df3)
    end
    str_others = houses |> @> @rsubset(:city == city_st)
    for m in 1:str_others.markets[1]
        cod_market = string(first(city_st,1)[1],m,"_m")
        transient_df5 = DataFrame(city = city_st, type = "market", code = cod_market)
        structures = append!(structures,transient_df5)
    end
    for o in 1:str_others.churches[1]
        cod_hospital = string(first(city_st,1)[1],o,"_c")
        transient_df6 = DataFrame(city = city_st, type = "church", code = cod_hospital)
        structures = append!(structures,transient_df6)
    end
    for s in 1:str_others.social[1]
        cod_social = string(first(city_st,1)[1],s,"_e")
        transient_df7 = DataFrame(city = city_st, type = "social", code = cod_social)
        structures = append!(structures,transient_df7)
    end
end

structures = structures |> @> @rtransform(:hours = hrs_str(:city,:type))

# Number of mosquitoes per structures

mosquitoes_df = DataFrame(Container_s = Integer[], Container_ml = Integer[], liters = Integer[], n_vectors = Integer[])

for s in 1:nrow(structures)
    city_s = structures.city[s]
    type_s = if structures.type[s] == "house"
        "Household"
    else
        "Non-Household"
    end
    container_df = cont_vect(city_s,type_s)
    mosquitoes_df = append!(mosquitoes_df,container_df)
end

vectors = hcat(structures,mosquitoes_df)

vectors |> @> groupby([:city, :type]) |> @combine(:vectors = sum(:n_vectors))

## Mosquito infection ==> look for initial mosquito p_inf

infs_h = population |> @> groupby(:city) |> @combine(:infs = sum(:inf + :latent))

kisumu_prev = vectors |> @> groupby(:city) |> @combine(:n_mosquitoes = sum(:n_vectors)) |> @rsubset(:city == "Kisumu") |>
@rtransform(:a = a_rate(24.7), :u = u_rate(24.7), :b = vector_comp(24.7)) |>
@rtransform(:prevalence = mosq_prev(:a,:b,:n_mosquitoes,pop_sum.pop[2],infs_h.infs[2],:u)) |>
@rtransform(:prob = :prevalence/:n_mosquitoes) |> _.prob[1]

ukunda_prev = vectors |> @> groupby(:city) |> @combine(:n_mosquitoes = sum(:n_vectors)) |> @rsubset(:city == "Ukunda") |>
@rtransform(:a = a_rate(27.56), :u = u_rate(27.56), :b = vector_comp(27.56)) |>
@rtransform(:prevalence = mosq_prev(:a,:b,:n_mosquitoes,pop_sum.pop[1],infs_h.infs[1],:u)) |>
@rtransform(:prob = :prevalence/:n_mosquitoes) |> _.prob[1]

vectors = vectors |> 
@> @rtransform(:inf = if :city == "Kisumu" 
    rand(Binomial(:n_vectors,kisumu_prev),1)[1]
else
    rand(Binomial(:n_vectors,ukunda_prev),1)[1]
end) |>
@rtransform(:exp = (:city == "Kisumu") ? round(Integer,:inf*vector_comp(24.7)) : round(Integer,:inf*vector_comp(27.56))) |>
@rtransform(:susc = :n_vectors-:inf) |> @select(:city,:type,:code,:hours,:Container_s,:Container_ml,:liters,:n_vectors,:susc,:exp,:inf) 

vectors |> @> groupby(:city) |> @combine(:vector = sum(:n_vectors), :sum_susc = sum(:susc), :sum_exp = sum(:exp),
:sum_inf = sum(:inf))

CSV.write("human_pop_20k.csv",population)
CSV.write("vector_pop_20k.csv",vectors)

