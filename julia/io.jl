
using DelimitedFiles

function io_saveArray(fileName, variable)
    io = open(fileName, "w")
    writedlm(io, variable, '\n')
    close(io)
end

function io_loadArray(fileName)
	io = open(fileName, "r")
	x = readdlm(io, '\n')
	close(io)
	return x
end

