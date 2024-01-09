using Pkg, Distributed
addprocs(16)

@everywhere begin
    using Pkg, Distributed, StatsBase, DataFrames, DataFramesMeta, Random, Distributions, Statistics, CSV, Dates
    population = CSV.read("human_pop_20k.csv", DataFrame)
    vectors = CSV.read("vector_pop_20k.csv", DataFrame)
    weather = CSV.read("weather2.csv", DataFrame)
    include("source_info_movimm.jl")
    replicates = 1:200
    days = 1:731
end

infections_def = @distributed (append!) for a in replicates
    infections_transient = DataFrame(City = String[], Day = Integer[], Date = Date[], Inf_h = Integer[], Total_m = Integer[],
    Inf_m = Integer[], I_HH = Integer[], I_s = Integer[], I_w = Integer[], I_m = Integer[],I_c = Integer[], I_e = Integer[],
    susc = Integer[], imm = Integer[])
    population2 = copy(population)
    vectors_2 = copy(vectors)
    for i in days
        update_humans(population2)
        for c in 1:length(cities)
            i_hh = Vector{Integer}()
            i_sc = Vector{Integer}()
            i_wp = Vector{Integer}()
            i_mk = Vector{Integer}()
            i_ch = Vector{Integer}()
            i_et = Vector{Integer}()
            size_h = nrow(population2[population2.city .== cities[c],:])
            day_weather = weather[weather.city .== cities[c],:][i,:]
            vectors2 = vectors_2[vectors_2.city .== cities[c],:]
            temp = day_weather.mean_t
            hr = (!ismissing(day_weather.rh)) ? day_weather.rh : (cities[c] == "Kisumu") ? 65.5 : 77.7
            for r in 1:nrow(vectors2)  ## Updating live and infected mosquitoes
                if vectors2.n_vectors[r] > 0
                    dying_mosquitoes(vectors_2,vectors2[r,:],temp)
                    updated_row = vectors_2[vectors_2.code .== vectors2.code[r],:]
                    vectors2.susc[r] = updated_row.susc[1]
                    vectors2.exp[r] = updated_row.exp[1]
                    vectors2.inf[r] = updated_row.inf[1]
                    vectors2.n_vectors[r] = updated_row.n_vectors[1]
                    emerging_mosquitoes(vectors_2,vectors2[r,:],temp)
                    vectors2.n_vectors[r] = vectors_2[vectors_2.code .== vectors2.code[r],:n_vectors][1]
                    vectors2.susc[r] = vectors_2[vectors_2.code .== vectors2.code[r],:susc][1]                
                    if vectors2.exp[r] > 0 ## How many mosquitoes move to infectious?
                        update_mosquitoes(vectors_2,vectors2[r,:],temp)
                    end
                    vectors2.exp[r] = vectors_2[vectors_2.code .== vectors2.code[r],:exp][1]
                    vectors2.inf[r] = vectors_2[vectors_2.code .== vectors2.code[r],:inf][1]                
                end
            end
            for s in 1:nrow(vectors2) ## Core of loop, where interactions and infections are taking place
                struc = vectors2[s,:]
                cod_s = struc.code
                hours = struc.hours
                if struc.type == "house" || struc.type == "school" || struc.type == "workplace" ## if is in a house, school or workplace
                    inf_h = if struc.type == "house" 
                        sum(population2[population2.house .== cod_s,:inf])
                    else
                        sum(population2[population2.id_place .== cod_s,:inf])
                    end
                    people = if struc.type == "house" 
                        population2[population2.house .== cod_s,:subj_cod]
                    else
                        population2[population2.id_place .== cod_s,:subj_cod]
                    end
                    visit = (struc.type == "house") ? visits(struc,population2) : "0" #Line 67
                    n_people = length(people)
                    if struc.inf > 0 || inf_h > 0 || (struc.type == "house" && visit != "0") #line 68                   
                        bites = bites_n(n_people,hours,struc,temp,visit) #line 69 ## number of total bites
                        if struc.inf > 0  ## if there are infected mosquitoes
                            v_h_inf_bite = rand(Hypergeometric(struc.inf,((struc.n_vectors-struc.inf) < 0 ? 0 : (struc.n_vectors-struc.inf)),sum(bites))) #line 71 ### Are those bites from an infected vector?
                            if v_h_inf_bite > 0                                
                                if !isempty(people)
                                    bitten = if struc.type != "house" || visit == "0" #Line 74 ## who were bitten
                                        sample(people,v_h_inf_bite)
                                    elseif visit != "0"
                                        push!(people,visit)
                                        sample(people,v_h_inf_bite)
                                    end 
                                    for b in bitten
                                        if population2[population2.subj_cod .== b, :susc][1] == 1 && population2[population2.subj_cod .== b, :recov][1] == 0 ## Is this subject susceptible?
                                            population2[population2.subj_cod .== b, :latent] .= 1
                                            if struc.type == "house"
                                                push!(i_hh,1)
                                            elseif struc.type == "school"
                                                push!(i_sc,1)
                                            else
                                                push!(i_wp,1)
                                            end
                                        end                                
                                    end
                                end                            
                            end
                        end
                        visit_inf = population2[population2.subj_cod .== visit,:] # line 90
                        if inf_h > 0 || (nrow(visit_inf) > 1 && visit_inf.inf[1] > 0) #line 90 ## if there are infected humans
                            if struc.susc > 0 ## if there are susceptible mosquitoes
                                humans_str = if struc.type == "house"
                                    population2[population2.house .== cod_s,:] 
                                else
                                    population2[population2.id_place .== cod_s,:]
                                end
                                prob_inf = (inf_h/nrow(humans_str))*(struc.susc/(struc.inf+struc.exp+struc.susc)) ## prob of bite an infected human by susceptible mosquito
                                inf_vector(vectors_2,prob_inf,temp, (struc.type == "house") ? bites[1] : bites,cod_s)
                                if struc.type == "house" && (nrow(visit_inf) > 1 && visit_inf.inf[1] == 1) #line 99
                                    pro_inf_susc = struc.susc/(struc.inf+struc.exp+struc.susc)
                                    inf_vector(vectors_2,pro_inf_susc,temp,bites[2],cod_s)
                                end
                            end          
                        end
                    end
                else  ### What happen in other random places
                    r_n_people = population2[population2.city .== cities[c],:subj_cod]
                    r_occup = sample(r_n_people,rand(DiscreteUniform(10,70)))
                    pop2 = population2[in(r_occup).(population2.subj_cod),:]
                    n_people = nrow(pop2)
                    Inf_h = sum(pop2.inf)
                    if struc.inf > 0 || Inf_h > 0
                        bites = bites_n(n_people,hours,struc,temp,"0") ## how many bites in total
                        if struc.inf > 0
                            v_h_inf_bite = rand(Hypergeometric(struc.inf,((struc.n_vectors-struc.inf) < 0 ? 0 : (struc.n_vectors-struc.inf)),sum(bites))) ### Are those bites from an infected vector?
                            if v_h_inf_bite > 0
                                people = pop2[:,:subj_cod]
                                bitten = sample(people,v_h_inf_bite)
                                for b in bitten
                                    if population2[population2.subj_cod .== b,:susc][1] == 1 && population2[population2.subj_cod .== b, :recov][1] == 0
                                        population2[population2.subj_cod .== b,:latent] .= 1
                                        if struc.type == "market"
                                            push!(i_mk,1)
                                        elseif struc.type == "church"
                                            push!(i_ch,1)
                                        else
                                            push!(i_et,1)
                                        end
                                    end
                                end
                            end
                        end
                        if Inf_h > 0
                            if struc.susc > 0
                                prob_inf = (Inf_h/nrow(pop2))*(struc.susc/(struc.inf+struc.exp+struc.susc))
                                inf_vector(vectors_2,prob_inf,temp,bites,cod_s)                                
                            end                        
                        end
                    end
                end
            end
            tot_inf_h = sum(population2[population2.city .== cities[c],:inf]) ## Collecting city-based daily information
            immunes = sum(population2[population2.city .== cities[c],:imm])
            tot_inf_m = sum(vectors_2[vectors_2.city .== cities[c],:inf])
            tot_m = sum(vectors_2[vectors_2.city .== cities[c],:n_vectors])
            susceptibles = sum(population2[population2.city .== cities[c],:susc].==1 .&& population2[population2.city .== cities[c],:recov].==0)
            transient_df = DataFrame(City = cities[c], Day = days[i], Date = day_weather.date, Inf_h = tot_inf_h, Total_m = tot_m, 
            Inf_m = tot_inf_m, I_HH = sum(i_hh), I_s = sum(i_sc), I_w = sum(i_wp), I_m = sum(i_mk), I_c = sum(i_ch), I_e = sum(i_et),
            susc = susceptibles, imm = immunes)
            append!(infections_transient,transient_df)
        end
    end
    infections_transient = @chain infections_transient begin
        @rtransform(:Week = floor(:Date, Dates.Week))
        groupby([:City,:Week])
        @combine(:Inf_h = sum(:Inf_h), :Total_m = sum(:Total_m),:Inf_m = sum(:Inf_m),:I_HH = sum(:I_HH), :I_s = sum(:I_s),
        :I_w = sum(:I_w), :I_m = sum(:I_m), :I_c = sum(:I_c), :I_e = sum(:I_e), :susc = last(:susc), :imm = last(:imm))
        @rtransform(:rep = a)
    end 
    infections_transient
end

CSV.write("infections_100.csv",infections_def)