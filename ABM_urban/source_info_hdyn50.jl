### Functions needed

function update_humans(population2::AbstractDataFrame)
    for t in 1:nrow(population2)  ## updating human infections and recoveries
        if population2.latent[t] == 1 && population2.lat_day[t] < 5
            population2.susc[t] = 0
            population2.lat_day[t] += 1
        elseif population2.latent[t] == 1 && population2.lat_day[t] == 5
            population2.latent[t] = 0
            population2.inf[t] = 1
            population2.lat_day[t] = 0
        elseif population2.inf[t] == 1 && population2.inf_day[t] < 7
            population2.inf_day[t] += 1 
        elseif population2.inf[t] == 1 && population2.inf_day[t] == 7
            population2.inf[t] = 0
            population2.recov[t] += 1
            population2.imm[t] = 1
            population2.inf_day[t] = 0
        elseif population2.recov[t] == 1 && population2.day_recov[t] < 100
            population2.day_recov[t] += 1
        elseif population2.recov[t] == 1 && population2.day_recov[t] == 100
            population2.day_recov[t] = 0
            population.recov[t] = 0
            population2.susc[t] = 1
        end
    end
end

function u_rate(t::Float64)
    lf = -1.48e-01*(t-9.16)*(t-37.73)
    return(1/lf)
end

function dying_mosquitoes(vectors_2::AbstractDataFrame,vectors2::DataFrameRow,temp::Float64)
    deaths_m = rand(Binomial(vectors2.n_vectors,u_rate(temp)))
    if deaths_m > 0
        p_stages = [vectors2.susc/vectors2.n_vectors,vectors2.exp/vectors2.n_vectors,vectors2.inf/vectors2.n_vectors]
        which_died = wsample(["susc","exp","inf"],p_stages,deaths_m)
        d_susc = sum(which_died .== "susc")
        d_inf = sum(which_died .== "inf")
        d_exp = sum(which_died .== "exp")
        if d_susc > vectors2.susc
            vectors_2[vectors_2.code .== vectors2.code,:susc] .= 0
        else
            vectors_2[vectors_2.code .== vectors2.code,:susc] .-= d_susc
        end
        if d_exp > vectors2.exp
            vectors_2[vectors_2.code .== vectors2.code,:exp] .= 0
        else
            vectors_2[vectors_2.code .== vectors2.code,:exp] .-= d_exp
        end        
        if d_inf > vectors2.inf
            vectors_2[vectors_2.code .== vectors2.code,:inf] .= 0
        else
            vectors_2[vectors_2.code .== vectors2.code,:inf] .-= d_inf
        end        
        new_data = vectors_2[vectors_2.code .== vectors2.code,:]
        new_n = new_data.susc+new_data.exp+new_data.inf    
        vectors_2[vectors_2.code .== vectors2.code,:n_vectors] .= new_n              
    end
end

function a_rate(t::Float64)
    a = (t < 13.35 || t > 40.08) ? 0 : 2.02e-04*t*(t-13.35)*(40.08-t)^0.5
    return(a)
end

function new_mosquitoes(t::Float64, Nm::Int, l::Int)
    if l > 0        
        efd_l = (t<14.58 || t > 34.61) ? 0 : (8.56e-3*t*(t-14.58)*(34.61-t)^0.5)
        efd = rand(Poisson(efd_l))
        mdr = (t<11.36 || t > 39.17) ? 0 : (7.86e-5*t*(t-11.36)*(39.17-t)^0.5)
        p_fed = 1-exp(-a_rate(t)/u_rate(t))
        fed = rand(Binomial(Nm,p_fed))
        dens = round(Integer,efd*fed)
        d = -0.165780731+(t*0.080343300)+((t^2)*-0.001384275)
        f_D = 1/(1+(exp((0.09*(dens/l))-0.55)))*d
        pEA = (t<13.56 || t > 38.29) ? 0 : (-5.99e-3*(t-13.56)*(t-38.29))
        r = mdr*pEA
        Nt_1 = Nm+rand(Binomial(dens,r*f_D))
        return(round(Integer,Nt_1))
    else
        return(Nm)
    end
end

function emerging_mosquitoes(vectors_2::AbstractDataFrame,vectors2::DataFrameRow,temp::Float64)
    new_n = new_mosquitoes(temp,vectors2.n_vectors,vectors2.liters)
    new_s = new_n - (vectors2.inf + vectors2.exp)
    vectors_2[vectors_2.code .== vectors2.code,:n_vectors] .= new_n
    vectors_2[vectors_2.code .== vectors2.code,:susc] .= new_s    
end

function pdr(t::Float64)
    eip = (t < 10.68 || t > 45.9) ? 0 : (6.65e-5*t*(t-10.68)*(45.9-t)^0.5)
    return(eip)
end

function update_mosquitoes(vectors_2::AbstractDataFrame,vectors2::DataFrameRow,temp::Float64)    
    new_infected = rand(Binomial(vectors2.exp,pdr(temp)))
    if new_infected > 0
        vectors_2[vectors_2.code .== vectors2.code,:exp] .-= new_infected
        vectors_2[vectors_2.code .== vectors2.code,:inf] .+= new_infected            
    end
end

function visits(struc::DataFrameRow, population2::AbstractDataFrame)
    p_visit = rand(Binomial(1,0.1))
    if p_visit == 1
        who = sample(population2[population2.city .== struc.city, :subj_cod])
        return(who)
    else
        return("0")
    end        
end

function bites_n(people,hrs,struc,temp,visit)
    n_biting = rand(Binomial(struc.n_vectors,a_rate(temp)))
    p_encounter = 1-exp(-(hrs/24)*people)
    n_bites = rand(Binomial(n_biting,p_encounter))
    if struc.type == "house" 
        visit_bites = if n_biting > n_bites && visit != "0"
            p_visit = 1-exp(-(rand(DiscreteUniform(1,3)/24)))
            rand(Binomial(n_biting-n_bites,p_visit))
        else
            0
        end
        return([n_bites,visit_bites])
    else
        return(n_bites)
    end
end

function vector_comp(t::Float64)
    b = (t<17.05 || t > 35.83) ? 0 : (8.49e-04*t*(t-17.05)*(35.83-t)^0.5)
    c = (t<12.22 || t > 37.46) ? 0 : (4.91e-04*t*(t-12.22)*(37.46-t)^0.5)
    return(b*c)
end

function inf_vector(vectors_2::AbstractDataFrame,prob_inf::Float64,temp::Float64,bites::Int,cod_s::String7)
    h_v_inf_bite = sum(rand(Binomial(bites,prob_inf)))
    bite_to_exp = rand(Binomial(h_v_inf_bite,vector_comp(temp))) ## How many mosquitoes get infected
    if bite_to_exp > 0
        suscept = vectors_2[vectors_2.code .== cod_s,:susc][1]
        if bite_to_exp < suscept
            vectors_2[vectors_2.code .== cod_s,:exp] .+= bite_to_exp
            vectors_2[vectors_2.code .== cod_s,:susc] .-= bite_to_exp
        else
            vectors_2[vectors_2.code .== cod_s,:exp] .+= suscept
            vectors_2[vectors_2.code .== cod_s,:susc] .= 0
        end
    end
end

### information
cities = ["Kisumu","Ukunda"]

## Dynamic model
function chooseurban(vectors_coord::AbstractDataFrame, pop_str::String)
    name_x = "x_"*pop_str
    name_y = "y_"*pop_str
    df = @select vectors_coord $(Between(:city,:inf)) :y = $name_y :x = $name_x
    return(df)
end

function chooseurbanh(population_coord::AbstractDataFrame, pop_str::String)
    name_m = pop_str*"market"
    name_c = pop_str*"church"
    name_p = pop_str*"place"
    df = @select population_coord $(Between(:city,:none)) $(Between(:susc,:day_recov)) :market = $name_m :church = $name_c :id_place = $name_p
    return(df)
end

function assess_intersect(mock_pop_k::AbstractDataFrame, point_x::Float64, point_y::Float64)
    values = (lower_x = point_x-0.15, upper_x=point_x+0.15, lower_y=point_y-0.15, upper_y=point_y+0.15)
    mock_pop = @select(mock_pop_k,:code,:y,:x)
    @rsubset!(mock_pop, (:x > values.lower_x) && (:x < values.upper_x) && (:y > values.lower_y) && (:y < values.upper_y))
    @rtransform!(mock_pop, :dist = sqrt(((point_x-:x)^2)+((point_y-:y)^2))*1000)
    a = 35.591
    b = 1.47541
    @rtransform!(mock_pop, :prob = ((1/(((2π)^(3/2))*b*:dist^2))*exp(-(log(:dist/a)^2)/2b^2)))
    @select!(mock_pop,:neighbors = :code,:prob)
    @rsubset!(mock_pop, isfinite(:prob))
    return(mock_pop)
end

function migration(Nh::Int,hrs::Float64,t::Float64,Nm::Int,l::Int,factor::Float64)
    reps = 1:20
    prod = []
    for i in reps
        sim = Nm/new_mosquitoes(t,Nm,l)
        push!(prod,(sim >= 1))
    end
    vectors = []
    for i in reps
        n_biting = rand(Binomial(Nm,a_rate(t)))
        p_encounter = 1-exp(-(hrs/24)*Nh)
        n_bites = rand(Binomial(n_biting,p_encounter))
        non_bit = n_biting - n_bites
        push!(vectors,(n_biting == 0) ? 0 : (non_bit/n_biting))
    end
    prob = (sum(prod)/20)*(mean(vectors))*factor
    mig = rand(Binomial(Nm,prob))
    return(mig)
end

function migrating(vectors_2::AbstractDataFrame,vectors2::AbstractDataFrame,vector_row::DataFrameRow,migrants::Int)
    inf_mig = rand(Hypergeometric(vector_row.inf,(vector_row.susc + vector_row.inf + vector_row.exp),migrants))
    inf_exp = rand(Hypergeometric(vector_row.exp,(vector_row.susc + vector_row.inf + vector_row.exp),(migrants-inf_mig)))
    pop_mock = assess_intersect(vectors2,vector_row.x,vector_row.y)
    if nrow(pop_mock) > 0 
        if inf_mig > 0 || inf_exp > 0
            vectors_2[vectors_2.code .== vector_row.code,:n_vectors] .-= migrants
            vectors_2[vectors_2.code .== vector_row.code,:susc] .-= (migrants-inf_mig-inf_exp)
            vectors_2[vectors_2.code .== vector_row.code,:inf] .-= inf_mig
            vectors_2[vectors_2.code .== vector_row.code,:exp] .-= inf_exp
            vectors2[vectors2.code .== vector_row.code,:n_vectors] .-= migrants
            vectors2[vectors2.code .== vector_row.code,:susc] .-= (migrants-inf_mig-inf_exp)
            vectors2[vectors2.code .== vector_row.code,:inf] .-= inf_mig
            vectors2[vectors2.code .== vector_row.code,:exp] .-= inf_exp
        else
            vectors_2[vectors_2.code .== vector_row.code,[:n_vectors,:susc]] .-= migrants
            vectors2[vectors2.code .== vector_row.code,[:n_vectors,:susc]] .-= migrants
        end
        @select!(pop_mock[in(wsample(pop_mock.neighbors,pop_mock.prob,migrants)).(pop_mock.neighbors),:],:neighbors)
        dest_inf = sample(pop_mock.neighbors,inf_mig)
        dest_exp = sample(pop_mock.neighbors,inf_exp)
        @rtransform!(pop_mock, :status = (:neighbors in dest_inf) ? "inf" : (:neighbors in dest_exp) ? "exp" : "susc")
        for r in nrow(pop_mock)
            if pop_mock.status == "inf"
                vectors_2[vectors_2.code .== pop_mock.neighbors[r],[:n_vectors,:inf]] .+= 1
                vectors2[vectors2.code .== pop_mock.neighbors[r],[:n_vectors,:inf]] .+= 1
            elseif pop_mock.status == "exp"
                vectors_2[vectors_2.code .== pop_mock.neighbors[r],[:n_vectors,:exp]] .+= 1
                vectors2[vectors2.code .== pop_mock.neighbors[r],[:n_vectors,:exp]] .+= 1
            else
                vectors_2[vectors_2.code .== pop_mock.neighbors[r],[:n_vectors,:susc]] .+= 1
                vectors2[vectors2.code .== pop_mock.neighbors[r],[:n_vectors,:susc]] .+= 1
            end
        end
    end
end

function human_nh(population2::AbstractDataFrame, vectors2::AbstractDataFrame, struc::DataFrameRow, city::String, distan::Float64, perc::Integer)
    people = sample(population2[population2.city .== city,:subj_cod],200)
    houses = String[]
    for p in people
        house = split(p,"_")[1]
        push!(houses,house)
    end
    distns = Float64[]
    for h in houses
        value_y = vectors2[vectors2.code .== h,:y][1]
        value_x = vectors2[vectors2.code .== h,:x][1]
        distance = sqrt(((struc.x-value_x)^2)+((struc.y-value_y)^2))
        push!(distns,distance)
    end
    index = findall(>(distan), distns)
    new_people = people[index]
    new_dists = distns[index]
    probs = 1 .- (new_dists./6)
    d = (perc == 100) ? DiscreteUniform(10,70) : (perc == 50) ? DiscreteUniform(5,35) : (perc == 20) ? DiscreteUniform(2,14) : DiscreteUniform(10,70)
    assistants = wsample(new_people,probs,rand(d))
    return(assistants)
end
