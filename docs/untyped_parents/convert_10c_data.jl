using DataFrame, CSV

f = open("SNP_def29a.in")
lines = readlines(f)

x = DataFrame(Locus = String[], Chromosome = Int[], 
	Basepairs = Int[], Allele1 = Int[], Allele2 = Int[])

for i in 2:length(lines) #start from 2 to skip header row
	cur = split(lines[i], ',') #turn current line to vector
	
	#get rid of white spaces
	for j in 1:length(cur)
		cur[j] = strip(cur[j])
	end

	#delete last column because we don't want it
	cur = cur[1:5]

	#add row to matrix
	push!(x, [cur[1], parse(Int64, cur[2]), parse(Int64, cur[3]),
		parse(Int64, cur[4]), parse(Int64, cur[5])])
end

CSV.write("SNP_def10c_converted.txt", x)