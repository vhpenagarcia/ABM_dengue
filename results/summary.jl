
using Pkg, CSV, DataFrames, DataFramesMeta, MCMCChains, StatsBase

## With Interquartile Range (IQR)
function summary_table(infections::AbstractDataFrame,Day::Integer ,Environment::String,Percentage::Integer,Container::String)
    temp_df1 = @chain infections begin
        @rtransform(:I_NH = :I_s+:I_w+:I_m+:I_e+:I_c)
        @rtransform(:total = :I_HH+:I_NH)
        groupby([:City,:Week])
        @combine(:Total = median(:total), :Household = median(:I_HH), :Non_Household = median(:I_NH), 
        :total_lower95 = quantile(:total,0.25), :total_upper95 = quantile(:total,0.75),
        :HH_lower95 = quantile(:I_HH,0.25), :HH_upper95 = quantile(:I_HH,0.75),
        :NH_lower95 = quantile(:I_NH,0.25), :NH_upper95 = quantile(:I_NH,0.75))
        groupby(:City)
        @combine(:Cases = sum(:Total), :caseslow95 = sum(:total_lower95), :casesup95 = sum(:total_upper95),
        :Household = sum(:Household), :HH_low95 = sum(:HH_lower95), :HH_up95 = sum(:HH_upper95),
        :Non_Household = sum(:Non_Household), :NH_low95 = sum(:NH_lower95), :NH_up95 = sum(:NH_upper95))
        @rtransform(:Environment = Environment, :Container = Container, :Percentage = Percentage, :Day = Day)
        @select(:City,:Environment,:Container,:Percentage,:Day,:Cases,:caseslow95,:casesup95,:Household,:HH_low95,:HH_up95,
        :Non_Household,:NH_low95,:NH_up95)    
    end
    temp_df2 = @chain infections begin
        @rtransform(:inc_cases = :I_HH+:I_s+:I_w+:I_m+:I_c+:I_e, :I_NH = :I_s+:I_w+:I_m+:I_c+:I_e)
        groupby([:City,:Week])
        @combine(:Cases = median(:inc_cases), :Cases_lower95 = quantile(:inc_cases,0.25), :Cases_upper95 = quantile(:inc_cases,0.75),
        :I_HH = median(:I_HH), :I_HH_lower95 = quantile(:I_HH,0.25), :I_HH_upper95 = quantile(:I_HH,0.75),
        :I_NH = median(:I_NH), :I_NH_lower95 = quantile(:I_NH,0.25), :I_NH_upper95 = quantile(:I_NH,0.75))
        @rtransform(:cases_inc = (:Cases > 0) ? 1 : 0, :caselow_inc = (:Cases_lower95 > 0) ? 1 : 0, :caseup_inc = (:Cases_upper95 > 0) ? 1 : 0,
        :I_HH_inc = (:I_HH > 0) ? 1 : 0, :I_HHlow_inc = (:I_HH_lower95 > 0) ? 1 : 0, :I_HHup_inc = (:I_HH_upper95 > 0) ? 1 : 0,
        :I_NH_inc = (:I_NH > 0) ? 1 : 0, :I_NHlow_inc = (:I_NH_lower95 > 0) ? 1 : 0, :I_NHup_inc = (:I_NH_upper95 > 0) ? 1 : 0)
        groupby(:City)
        @combine(:pweeks_cases = mean(:cases_inc .> 0), :pwk_low95 = mean(:caselow_inc .> 0), :pwk_up95 = mean(:caseup_inc .> 0),
        :pweeks_HHcases = mean(:I_HH_inc .> 0), :pwk_HHlow95 = mean(:I_HHlow_inc .> 0), :pwk_HHup95 = mean(:I_HHup_inc .> 0),
        :pweeks_NHcases = mean(:I_NH_inc .> 0), :pwk_NHlow95 = mean(:I_NHlow_inc .> 0), :pwk_NHup95 = mean(:I_NHup_inc .> 0),
        :weeks_cases = mean(:cases_inc .> 0)*105, :wk_low95 = mean(:caselow_inc .> 0)*105, :wk_up95 = mean(:caseup_inc .> 0)*105,
        :weeks_HHcases = mean(:I_HH_inc .> 0)*105, :wk_HHlow95 = mean(:I_HHlow_inc .> 0)*105, :wk_HHup95 = mean(:I_HHup_inc .> 0)*105,
        :weeks_NHcases = mean(:I_NH_inc .> 0)*105, :wk_NHlow95 = mean(:I_NHlow_inc .> 0)*105, :wk_NHup95 = mean(:I_NHup_inc .> 0)*105)
        @rtransform(:Environment = Environment, :Container = Container, :Percentage = Percentage, :Day = Day)
        @select(:City,:Environment,:Container,:Percentage,:Day,:pweeks_cases,:pwk_low95,:pwk_up95,:weeks_cases,:wk_low95,:wk_up95,
        :pweeks_HHcases,:pwk_HHlow95,:pwk_HHup95,:weeks_HHcases,:wk_HHlow95,:wk_HHup95,
        :pweeks_NHcases,:pwk_NHlow95,:pwk_NHup95,:weeks_NHcases,:wk_NHlow95,:wk_NHup95)     
    end
    df_all = innerjoin(temp_df1,temp_df2, on = [:City,:Environment,:Container,:Percentage,:Day])
    return(df_all)
end


function create_summary(control::AbstractDataFrame)
    summary = DataFrame(City = String[], Environment = String[], Container = String[], Percentage = Integer[], 
    Day = Integer[],Cases = Float16[], caseslow95 = Float16[],casesup95 = Float16[], Household = Float16[], 
    HH_low95 = Float16[],HH_up95 = Float16[],Non_Household = Float16[], NH_low95 = Float16[], NH_up95 = Float16[],
    pweeks_cases = Float64[], pwk_low95 = Float64[], pwk_up95 = Float64[], weeks_cases = Float16[], wk_low95 = Float16[], 
    wk_up95 = Float16[],pweeks_HHcases = Float64[], pwk_HHlow95 = Float64[], pwk_HHup95 = Float64[], 
    weeks_HHcases = Float16[], wk_HHlow95 = Float16[],wk_HHup95 = Float16[], pweeks_NHcases = Float64[], 
    pwk_NHlow95 = Float64[], pwk_NHup95 = Float64[], weeks_NHcases = Float16[], wk_NHlow95 = Float16[], 
    wk_NHup95 = Float16[])
    for r in 1:nrow(control)
        code = convert(String,control[r,:Run_code])
        environment = convert(String,control[r,:Environment])
        container = convert(String,control[r,:Container_type])
        percentage = control[r,:Percentage]
        day = control[r,:Day]
        name1 = join(["infections_out_",code,"_1.csv"])
        name2 = join(["infections_out_",code,"_2.csv"])
        file = CSV.read(name1, DataFrame)
        file2 = @rtransform(CSV.read(name2, DataFrame), :rep = :rep+200)
        append!(file,file2)
        transient_df = summary_table(file,day,environment,percentage,container)
        append!(summary,transient_df)
    end
    return(summary)
end


## Reactive table
control_table = @rsubset(CSV.read("Control_results_ABM.csv", DataFrame), :Container_type != "None")

summary_results = create_summary(control_table)
infections_100 = CSV.read("infections_100.csv",DataFrame)
cero_summary = summary_table(infections_100,0,"None",0,"None")
append!(summary_results,cero_summary)

CSV.write("summary_reactive_iqr.csv",summary_results)

## Preventive table
control_prev = CSV.read("Control_preventive_ABM.csv", DataFrame)

preventive_results = create_summary(control_prev)

CSV.write("summary_preventive_iqr.csv",preventive_results)

