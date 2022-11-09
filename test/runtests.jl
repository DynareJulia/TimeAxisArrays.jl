using ExtendedDates
using TimeAxisArrays
using Test


ta1 = TimeAxisArray(rand(4,3), Undated(11):Undated(14), [:a, :b, :c])
ta2 = TimeAxisArray(rand(4,3), period(Week, 1935, 2):period(Week, 1935, 5), [:a, :b, :c])

@show ta2[:b] 
@show ta2[WeekSE(1935, 3)]
@show ta2[:b, :c]
@show ta2[WeekSE(1935, 3)..WeekSE(1935, 4)]
@show ta2[WeekSE(1935, 3)..WeekSE(1935, 4), [:b, :c]]
@show ta2[WeekSE(1935, 3), :c]

@test ta2.b.data.data ≈ ta2[:b].data
@test ta2[:b].data ≈ ta2.data[:,:b].data
@test ta2[WeekSE(1935, 3)].data ≈ reshape(ta2.data[2, :], 1, 3)
@test ta2[:b, :c].data ≈ ta2.data[:, [2, 3]]
@test ta2[WeekSE(1935, 3)..WeekSE(1935, 4)].data ≈ ta2.data[[2, 3], :]
@test ta2[WeekSE(1935, 3)..WeekSE(1935, 4), [:b, :c]].data ≈ ta2.data[[2, 3], [2,3]]
@test ta2[WeekSE(1935, 3), :c].data.data[1,1] == ta2.data[2, 3]

@test all(log.(ta1).data ≈ log.(ta1.data))
@test ta1 + ta1 == 2*ta1
