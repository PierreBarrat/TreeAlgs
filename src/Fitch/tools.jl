_alphabet(s::LongSequence{T}) where T = BioSequences.symbols(T())
function get_variable_positions(t::Tree, seqkey)
	L = length(TreeTools.recursive_get(first(t.lleaves)[2].data.dat, seqkey))
	keep = _alphabet(TreeTools.recursive_get(first(t.lleaves)[2].data.dat, seqkey))
	#
	variable = Int64[]
	fixed = collect(1:L)
	state = Array{Any,1}(missing, L)
	# 
	for n in values(t.lleaves)
		seq = TreeTools.recursive_get(n.data.dat, seqkey)
	    todel = Int64[]
	    for i in fixed
	    	if ismissing(state[i]) && !isgap(seq[i]) && in(seq[i], keep) 
	    	# If state not initialized for i, initialize
	    		state[i] = seq[i]
	    	elseif !ismissing(state[i]) && seq[i] != state[i] && !isgap(seq[i]) &&in(seq[i], keep) 
	    	# If state was initialized and changed, this is a variable column
	    		push!(variable, i)
	    		push!(todel, i)
	    	end
	    end	
	    # Deleting variable positions from fixed
	    for i in todel
	    	deleteat!(fixed, findfirst(==(i), fixed))
	    end
	end
	return sort(variable)
end

_in(x) = y -> in(y, x)
function intersect!(aFs::FitchState, fstates::Vararg{FitchState})
	for fs in fstates
		filter!(_in(fs.state), aFs.state)
	end
end
function union!(aFs::FitchState, fs::FitchState)
	for a in Iterators.filter(!_in(aFs.state), fs.state)
		push!(aFs.state,a)
	end
end
