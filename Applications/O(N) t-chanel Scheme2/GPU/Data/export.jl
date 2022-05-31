using DelimitedFiles
datadir=@__DIR__
writedlm(joinpath(datadir,"Im.dat"),Array(solselfGPU.u[end][:,:,1]))
writedlm(joinpath(datadir,"Re.dat"),Array(solselfGPU.u[end][:,:,2]))
