"""
	mutable struct LBIData <: TreeNodeData

Data used to compute the Local Branching Index.
"""
mutable struct LBIData <: TreeNodeData
	tau::Union{Missing, Float64}
	message_down::Float64
	message_up::Float64
	lbi::Float64
	date
	alive::Bool
end
function LBIData(; tau=0.,
				message_down=0.,
				message_up=0.,
				LBI=0.,
				date=0.,
				alive=true)
	return LBIData(tau, message_down, message_up, LBI, date, alive)
end
