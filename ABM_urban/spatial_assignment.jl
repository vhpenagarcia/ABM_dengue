# Kwale (Ukunda): urban Households: 37,005; area: 63 km2
# Kisumu: urban households: 131,732; area: 139 km2

using StatsBase, Distributions, DataFrames, DataFramesMeta, CSV

vectors = CSV.read("vector_pop_20k.csv", DataFrame)

## Area according to counties: not considered

## Area according to city
area_ukunda = (25.5/15948)*3415
area_kisumu = (2085.4/1155574)*5205

dims_u = sqrt(area_ukunda)
dims_k = sqrt(area_kisumu)

mock_pop_u = DataFrame(x=rand(Uniform(1.25,1.25+(dims_u-0.3)),3415), y=rand(Uniform(1.25,1.25+(dims_u-0.3)),3415))
mock_pop_k = DataFrame(x=rand(Uniform(1.25,1.25+(dims_k-0.3)),5205), y=rand(Uniform(1.25,1.25+(dims_k-0.3)),5205))

@rtransform!(mock_pop_u, :y2 = :y*rand(Normal(1,0.15)), :x2 = :x*rand(Normal(1,0.15)))
@rtransform!(mock_pop_k, :y2 = :y*rand(Normal(1,0.15)), :x2 = :x*rand(Normal(1,0.15)))

### Randomly assigned type of structure: Scattered
#Kisumu
precode_k = 1:5205
mock_pop_k = hcat(mock_pop_k,precode_k)
NH_index = sample(mock_pop_k.x1,793, replace = false)
mock_pop_k_as = @chain mock_pop_k begin
    @rtransform(:env = (issubset(:x1,NH_index)) ? "Non-Household" : "Household")
    @select(:y_coord = :y2, x_coord = :x2, :x1, :env)
end

# Ukunda
precode_u = 1:3415
mock_pop_u = hcat(mock_pop_u,precode_u)
NH_index_u = sample(mock_pop_u.x1,648, replace = false)
mock_pop_u_as = @chain mock_pop_u begin
    @rtransform(:env = (issubset(:x1,NH_index_u)) ? "Non-Household" : "Household")
    @select(:y_coord = :y2, x_coord = :x2, :x1, :env)
end

## Centered

# Kisumu
med_x = mean(mock_pop_k_as.x_coord)
med_y = mean(mock_pop_k_as.y_coord)

mid_val = findall(mock_pop_k_as.x_coord .> med_x-0.025 .&& mock_pop_k_as.x_coord .< med_x+0.025 .&& mock_pop_k_as.y_coord .> med_y-0.025 .&& mock_pop_k_as.y_coord .< med_y+0.025)
middle_x2 = mock_pop_k_as[mid_val,:x_coord][1]
middle_y2 = mock_pop_k_as[mid_val,:y_coord][1]

mock_pop_k2 = @chain mock_pop_k begin
    @rtransform(:x2_dif = (:x2-middle_x2)^2, :y2_dif = (:y2-middle_y2)^2)
    @rtransform(:diff = :x2_dif+:y2_dif)
end
ranks = ordinalrank(mock_pop_k2.diff)

mock_pop_k2 = @chain hcat(mock_pop_k,ranks, makeunique = true) begin
    @rtransform(:prob_type = (:x1_1 <= 793) ? 0.95 : 0.05)
    @rtransform(:type = (rand(Binomial(1,:prob_type)) == 1) ? "Non-Household" : "Household")
    @select(:x2,:y2,:type)
end

@by mock_pop_k2 :type :n = length(:y2)
mock_pop_k2 = @chain hcat(mock_pop_k,ranks, makeunique = true) begin
    @rtransform(:prob_type = (:x1_1 <= 714) ? 0.95 : 0.05)
    @select(:y2,:x2,:x1,:prob_type)
end

NH_index_kcenter = sample(mock_pop_k2.x1,Weights(mock_pop_k2.prob_type),793, replace=false)

mock_pop_k_as2 = @chain mock_pop_k_as begin
    @rtransform(:env_center = (issubset(:x1,NH_index_kcenter)) ? "Non-Household" : "Household")
    @select(:y_coord, :x_coord, :x1, :env_rand = :env, :env_center)
end


## Ukunda
med_x = mean(mock_pop_u_as.x_coord)
med_y = mean(mock_pop_u_as.y_coord)

mid_val_u = findall(mock_pop_u_as.x_coord .> med_x-0.025 .&& mock_pop_u_as.x_coord .< med_x+0.025 .&& mock_pop_u_as.y_coord .> med_y-0.025 .&& mock_pop_u_as.y_coord .< med_y+0.025)
middle_x2_u = mock_pop_u_as[mid_val_u,:x_coord][1]
middle_y2_u = mock_pop_u_as[mid_val_u,:y_coord][1]

mock_pop_u2 = @chain mock_pop_u begin
    @rtransform(:x2_dif = (:x2-middle_x2_u)^2, :y2_dif = (:y2-middle_y2_u)^2)
    @rtransform(:diff = :x2_dif+:y2_dif)
end
ranks_u = ordinalrank(mock_pop_u2.diff)

mock_pop_u2 = @chain hcat(mock_pop_u,ranks_u, makeunique = true) begin
    @rtransform(:prob_type = (:x1_1 <= 583) ? 0.95 : 0.05)
    @select(:y2,:x2,:x1,:prob_type)
end

NH_index_ucenter = sample(mock_pop_u2.x1,Weights(mock_pop_u2.prob_type),648, replace=false)

mock_pop_u_as2 = @chain mock_pop_u_as begin
    @rtransform(:env_center = (issubset(:x1,NH_index_ucenter)) ? "Non-Household" : "Household")
    @select(:y_coord, :x_coord, :x1, :env_rand = :env, :env_center)
end


## many centers (3): Clustered

# Kisumu
x_multiplek = sample(mock_pop_k.x2,3)
y_multiplek = sample(mock_pop_k.y2,3)

mock_pop_k3 = @chain mock_pop_k begin
    @rtransform(:x1_dif = (:x2-x_multiplek[1])^2, :x2_dif = (:x2-x_multiplek[2])^2, :x3_dif = (:x2-x_multiplek[3])^2,
    :y1_dif = (:y2-y_multiplek[1])^2, :y2_dif = (:y2-y_multiplek[2])^2, :y3_dif = (:y2-y_multiplek[3])^2)
    @rtransform(:diff1 = :x1_dif+:y1_dif, :diff2 = :x2_dif+:y2_dif, :diff3 = :x3_dif+:y3_dif)
    @select(:x2,:y2,:x1,:diff1,:diff2,:diff3)
end
ranks1 = ordinalrank(mock_pop_k3.diff1)
ranks2 = ordinalrank(mock_pop_k3.diff2)
ranks3 = ordinalrank(mock_pop_k3.diff3)

mock_pop_k3 = @chain hcat(mock_pop_k,ranks1,ranks2,ranks3,makeunique=true) begin
    @rtransform(:prob_type = (:x1_1 <= 264 || :x1_2 <= 264 || :x1_3 <= 264) ? 0.95 : 0.05)
    @select(:x2,:y2,:x1,:prob_type)
end

NH_index_kscater = sample(mock_pop_k3.x1,Weights(mock_pop_k3.prob_type),793, replace=false)

mock_pop_k_as3 = @chain mock_pop_k_as2 begin
    @rtransform(:env_scater = (issubset(:x1,NH_index_kscater)) ? "Non-Household" : "Household")
    @select(:y_coord, :x_coord, :env_rand, :env_center, :env_scater)
end


# Ukunda
x_multipleu = sample(mock_pop_u.x2,3)
y_multipleu = sample(mock_pop_u.y2,3)

mock_pop_u3 = @chain mock_pop_u begin
    @rtransform(:x1_dif = (:x2-x_multipleu[1])^2, :x2_dif = (:x2-x_multipleu[2])^2, :x3_dif = (:x2-x_multipleu[3])^2,
    :y1_dif = (:y2-y_multipleu[1])^2, :y2_dif = (:y2-y_multipleu[2])^2, :y3_dif = (:y2-y_multipleu[3])^2)
    @rtransform(:diff1 = :x1_dif+:y1_dif, :diff2 = :x2_dif+:y2_dif, :diff3 = :x3_dif+:y3_dif)
    @select(:x2,:y2,:x1,:diff1,:diff2,:diff3)
end
ranks1 = ordinalrank(mock_pop_u3.diff1)
ranks2 = ordinalrank(mock_pop_u3.diff2)
ranks3 = ordinalrank(mock_pop_u3.diff3)

mock_pop_u3 = @chain hcat(mock_pop_u,ranks1,ranks2,ranks3,makeunique=true) begin
    @rtransform(:prob_type = (:x1_1 <= 216 || :x1_2 <= 216 || :x1_3 <= 216) ? 0.95 : 0.05)
    @select(:x2,:y2,:x1,:prob_type)
end

NH_index_uscater = sample(mock_pop_u3.x1,Weights(mock_pop_u3.prob_type),648, replace=false)

mock_pop_u_as3 = @chain mock_pop_u_as2 begin
    @rtransform(:env_scater = (issubset(:x1,NH_index_uscater)) ? "Non-Household" : "Household")
    @select(:y_coord, :x_coord, :env_rand, :env_center, :env_scater)
end


@rtransform!(mock_pop_k_as3, :city = "Kisumu")
@rtransform!(mock_pop_u_as3, :city = "Ukunda")

## catenating data frames
coordinates_abm = @chain vcat(mock_pop_k_as3, mock_pop_u_as3) begin
    @select(:city, :y_coord, :x_coord, :env_rand, :env_cent = :env_center, :env_scat = :env_scater)
end

## Creating the final file with coordinates

function coord_assignment(vectors::AbstractDataFrame, coordinates::AbstractDataFrame)
    cities = ["Kisumu","Ukunda"]
    conformations = ["env_rand","env_cent","env_scat"]
    transient_df = DataFrame(city = String[], type = String[], code = String[], y_rand = Float64[], x_rand = Float64[],
    y_cent = Float64[], x_cent = Float64[], y_scat = Float64[], x_scat = Float64[])
    for c in cities
        houses = @chain vectors begin
            @rsubset(:city == c && :type == "house")
            @select(:city,:type,:code)
        end
        non_houses = @chain vectors begin
            @rsubset(:city == c && :type != "house")
            @select(:city,:type,:code)
        end
        for r in conformations
            name_x = "x"*last(r,5)
            name_y = "y"*last(r,5)
            house_coords = @chain coordinates begin
                @select(:city,:y_coord,:x_coord,cols(Symbol(r)))
                @rsubset(:city == c && cols(Symbol(r)) == "Household")
                @select(cols(Symbol(name_y)) = :y_coord, cols(Symbol(name_x)) = :x_coord)
            end
            houses = hcat(houses,house_coords)
            non_house_coords = @chain coordinates begin
                @select(:city,:y_coord,:x_coord,cols(Symbol(r)))
                @rsubset(:city == c && cols(Symbol(r)) == "Non-Household")
                @select(cols(Symbol(name_y)) = :y_coord, cols(Symbol(name_x)) = :x_coord)
            end
            non_houses = hcat(non_houses,non_house_coords)
        end
        catenated_df = vcat(houses,non_houses)
        append!(transient_df,catenated_df)
    end
    vectors2 = innerjoin(vectors,transient_df, on=[:city,:type,:code])
    return(vectors2)
end

vectors_coord = coord_assignment(vectors,coordinates_abm)

CSV.write("vectors_coordC.csv", vectors_coord)

#### Human populations
population = CSV.read("human_pop_20k.csv", DataFrame)

### Functions to assign locations either the closest, or at least x number of meters.

### Closest distances
function local_assignations(vector::AbstractDataFrame, population::AbstractDataFrame)
    cities = ["Ukunda","Kisumu"]
    urban = ["rand","cent","scat"]
    NHs = ["market","church"]
    randmarket = String[]
    randchurch = String[]
    centmarket = String[]
    centchurch = String[]
    scatmarket = String[]
    scatchurch = String[]
    for c in cities
        houses = unique(population[population.city .== c,:house])
        for h in houses, u in urban, n in NHs
            yname = "y_"*u
            xname = "x_"*u
            y_house = vector[vector.code .== h,Symbol(yname)][1]
            x_house = vector[vector.code .== h,Symbol(xname)][1]
            inhab = nrow(population[population.house .== h,:])
            vector_loc = @select(vector[vector.city .== c .&& vector.type .== n,:],:code,:y = $yname,:x = $xname)
            @rtransform!(vector_loc,:dist = sqrt(((x_house-:x)^2)+((y_house-:y)^2)))
            sort!(vector_loc,:dist)
            lower = first(vector_loc[:,:code],inhab)
            if u == "rand"
                if n == "market"
                    randmarket = vcat(randmarket,lower)
                elseif n == "church"
                    randchurch = vcat(randchurch,lower)
                end
            elseif u == "cent"
                if n == "market"
                    centmarket = vcat(centmarket,lower)
                elseif n == "church"
                    centchurch = vcat(centchurch,lower)
                end
            elseif u == "scat"
                if n == "market"
                    scatmarket = vcat(scatmarket,lower)
                elseif n == "church"
                    scatchurch = vcat(scatchurch,lower)
                end
            end
        end        
    end
    population2 = copy(population)
    population2[:,:randmarket] = randmarket
    population2[:,:randchurch] = randchurch
    population2[:,:centmarket] = centmarket
    population2[:,:centchurch] = centchurch
    population2[:,:scatmarket] = scatmarket
    population2[:,:scatchurch] = scatchurch
    return(population2)
end


function occup_assignations(vector::AbstractDataFrame, population::AbstractDataFrame)
    urban = ["rand","cent","scat"]
    capacity = DataFrame(places = vector[vector.type .== "workplace" .|| vector.type .== "school",:code], number = 0)
    randplace = String[]
    centplace = String[]
    scatplace = String[]
    for p in 1:nrow(population), u in urban
        yname = "y_"*u
        xname = "x_"*u
        house = population.house[p]
        y_house = vector[vector.code .== house,Symbol(yname)][1]
        x_house = vector[vector.code .== house,Symbol(xname)][1]
        place = if last(population.id_place[p],1) == "w"
            locations = vector[vector.city .== population.city[p] .&& vector.type .== "workplace",[:code,Symbol(yname),Symbol(xname)]]
            @rtransform!(locations,:dist = sqrt(((x_house-$xname)^2)+((y_house-$yname)^2)))
            avaloc = capacity[last.(capacity.places,1) .== "w" .&& capacity.number .< 19, :places]
            locations[in(avaloc).(locations.code),:]
            sort!(locations,:dist)
            sample(first(locations[:,:code],3))
        elseif last(population.id_place[p],1) == "s"
            locations = vector[vector.city .== population.city[p] .&& vector.type .== "school",[:code,Symbol(yname),Symbol(xname)]]
            @rtransform!(locations,:dist = sqrt(((x_house-$xname)^2)+((y_house-$yname)^2)))
            avaloc = capacity[last.(capacity.places,1) .== "s" .&& capacity.number .< 360, :places]
            locations[in(avaloc).(locations.code),:]
            sort!(locations,:dist)
            sample(first(locations[:,:code],3))
        else
            population.id_place[p]
        end
        if u == "rand"
            randplace = vcat(randplace,place)
        elseif u == "cent"
            centplace = vcat(centplace,place)
        elseif u == "scat"
            scatplace = vcat(scatplace,place)
        end
        capacity[capacity.places .== place,:number] .+= 1
    end
    population2 = copy(population)
    population2[:,:randplace] = randplace
    population2[:,:centplace] = centplace
    population2[:,:scatplace] = scatplace
    return(population2)
end


### at least x km distances
function local_assignations(vector::AbstractDataFrame, population::AbstractDataFrame)
    cities = ["Ukunda","Kisumu"]
    urban = ["rand","cent","scat"]
    NHs = ["market","church"]
    randmarket = String[]
    randchurch = String[]
    centmarket = String[]
    centchurch = String[]
    scatmarket = String[]
    scatchurch = String[]
    for c in cities
        houses = unique(population[population.city .== c,:house])
        for h in houses, u in urban, n in NHs
            yname = "y_"*u
            xname = "x_"*u
            y_house = vector[vector.code .== h,Symbol(yname)][1]
            x_house = vector[vector.code .== h,Symbol(xname)][1]
            inhab = nrow(population[population.house .== h,:])
            vector_loc = @select(vector[vector.city .== c .&& vector.type .== n,:],:code,:y = $yname,:x = $xname)
            @rtransform!(vector_loc,:dist = sqrt(((x_house-:x)^2)+((y_house-:y)^2)))
            sort!(vector_loc,:dist)
            @rsubset!(vector_loc,:dist > 0.5) ## for 500 meters: ((u == "cent" && n == "church") ? 0.7 : 0.8)
            lower = if nrow(vector_loc) >= inhab
                first(vector_loc[:,:code],inhab)
            else
                add = sample(vector_loc[:,:code],(inhab - nrow(vector_loc)))
                vcat(vector_loc[:,:code],add)
            end
            if u == "rand"
                if n == "market"
                    randmarket = vcat(randmarket,lower)
                elseif n == "church"
                    randchurch = vcat(randchurch,lower)
                end
            elseif u == "cent"
                if n == "market"
                    centmarket = vcat(centmarket,lower)
                elseif n == "church"
                    centchurch = vcat(centchurch,lower)
                end
            elseif u == "scat"
                if n == "market"
                    scatmarket = vcat(scatmarket,lower)
                elseif n == "church"
                    scatchurch = vcat(scatchurch,lower)
                end
            end
        end        
    end
    population2 = copy(population)
    population2[:,:randmarket] = randmarket
    population2[:,:randchurch] = randchurch
    population2[:,:centmarket] = centmarket
    population2[:,:centchurch] = centchurch
    population2[:,:scatmarket] = scatmarket
    population2[:,:scatchurch] = scatchurch
    return(population2)
end

function occup_assignations(vector::AbstractDataFrame, population::AbstractDataFrame)
    urban = ["rand","cent","scat"]
    capacity = DataFrame(places = vector[vector.type .== "workplace" .|| vector.type .== "school",:code], number = 0)
    randplace = String[]
    centplace = String[]
    scatplace = String[]
    for p in 1:nrow(population), u in urban
        yname = "y_"*u
        xname = "x_"*u
        house = population.house[p]
        y_house = vector[vector.code .== house,Symbol(yname)][1]
        x_house = vector[vector.code .== house,Symbol(xname)][1]
        place = if last(population.id_place[p],1) == "w"
            locations = vector[vector.city .== population.city[p] .&& vector.type .== "workplace",[:code,Symbol(yname),Symbol(xname)]]
            @rtransform!(locations,:dist = sqrt(((x_house-$xname)^2)+((y_house-$yname)^2)))
            avaloc = capacity[last.(capacity.places,1) .== "w" .&& capacity.number .< 19, :places]
            locations[in(avaloc).(locations.code),:]
            @rsubset!(locations,:dist > 0.5) ## for 500 meters: ((u == "cent" && n == "church") ? 0.7 : 0.8)
            sort!(locations,:dist)
            sample(first(locations[:,:code],4))
        elseif last(population.id_place[p],1) == "s"
            locations = vector[vector.city .== population.city[p] .&& vector.type .== "school",[:code,Symbol(yname),Symbol(xname)]]
            @rtransform!(locations,:dist = sqrt(((x_house-$xname)^2)+((y_house-$yname)^2)))
            avaloc = capacity[last.(capacity.places,1) .== "s" .&& capacity.number .< 360, :places]
            locations[in(avaloc).(locations.code),:]
            @rsubset!(locations,:dist > 0.5) ## for 500 meters: ((u == "cent" && n == "church") ? 0.7 : 0.8)
            sort!(locations,:dist)
            sample(first(locations[:,:code],4))
        else
            population.id_place[p]
        end
        if u == "rand"
            randplace = vcat(randplace,place)
        elseif u == "cent"
            centplace = vcat(centplace,place)
        elseif u == "scat"
            scatplace = vcat(scatplace,place)
        end
        capacity[capacity.places .== place,:number] .+= 1
    end
    population2 = copy(population)
    population2[:,:randplace] = randplace
    population2[:,:centplace] = centplace
    population2[:,:scatplace] = scatplace
    return(population2)
end

population_coordC_o5 = local_assignations(vectors_coord,population)
population_coordC_05k = occup_assignations(vectors_coord,population_coordC_o5)

CSV.write("population_coordC_05k.csv",population_coordC_05k)
