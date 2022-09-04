using ExtendedDates
using TimeAxisArrays
using Test


ta1 = TimeAxisArray(rand(4,3), Undated(11):Undated(14), [:a, :b, :c])
ta2 = TimeAxisArray(rand(4,3), period(Week, 1935, 2):period(Week, 1935, 5), [:a, :b, :c])

@show ta2[:b] 
@show ta2[period(Week, 1935, 3)]
@show ta2[:b, :c]
@show ta2[period(Week, 1935, 3)..period(Week, 1935, 4)]
@show ta2[period(Week, 1935, 3)..period(Week, 1935, 4), :b, :c]
@show ta2[period(Week, 1935, 3), :c]
