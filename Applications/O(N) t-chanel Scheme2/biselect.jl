function massmin(inputfun,inimin,inimax;N_kernels=3)
    minidex=inimin
    maxidex=inimax
    list=SharedArray(range(inimin,inimax,N_kernels+2)[2:end-1]|>collect)
    masslist=pmap(inputfun,list)
    for i in 1:N_kernels
        if masslist[i]>0.0
            maxidex=list[i]
            if i>1.5
                minidex=list[i-1]
            end
            break
        elseif i==N_kernels
            minidex=list[i]
        end
    end
    println("massmid=",masslist[ceil(Int64,N_kernels/2)],"\n"," ","minval=",minidex,"maxval=",maxidex)
    return (minidex,maxidex)
end


function getmassmin(inputfun, inimin, inimax; N_kernels = 3, N_iters = 5)
    minidex = inimin
    maxidex = inimax
    for i = 1:N_iters
        (minidex, maxidex) =
            massmin(inputfun, minidex, maxidex, N_kernels = N_kernels)
    end
    return [minidex, maxidex]
end
