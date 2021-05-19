"""
	CharTable(dat, labels=collect(1:size(dat,2)))


Representation of a binary alignment with sequence labels.
```
	CharTable{T}
dat::Array{Bool,2}
labels::Array{T,1}
```
"""
struct CharTable{T}
	dat::Array{Bool,2}
	labels::Array{T,1}
	CharTable(dat, labels::Array{T,1}) where T = (length(labels) != size(dat, 1) ? error("Number of sequences ($(size(dat,1))) and labels ($(length(labels))) differ") : new{T}(dat, labels)) #(unique(dat, dims=2), labels))
end

function CharTable(aln)
	dat = Array{Bool,2}(aln)
	labels = collect(1:size(dat,2))
	return CharTable(dat, labels)
end

