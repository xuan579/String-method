
 # just need to add stringenergy and stringforces to complete a string method,if necessary, we need to add the type of configurations when computing energy and forces
"""
`string`:

`stringstring(yA::Array{Float64,2},yB::Array{Float64,2};M::Int64=20,dt::Float64=0.02,Nt::Int64=200,dir=1,redis=1,tol::Float64=1e-6,tau=1,o::Int64=0)` : the string method with two fixed minimas yA and yB.

    input :

    yA,yB : minimas.
    M : numberof images on string(not including yA and yB).
    dt : the time step.
    Nt : the upper bound of the times of string evoluiton.
    dir : choose the evolution direction, if dir=0, -∇V; if dir=1, -∇V normal to string.
    redis : choose to do redistribution or not, redis=1, do redis.
    tau : choose the way to get tangent; tau=0,central difference; tau=1,upwind scheme;
    tol : the tolerance.
    o : decide the output

    output:
    o = 0, p : MEP
    o = 1, Ap: all the string during the process.

 `string(phi::Array{Float64,3}; dt::Float64 = 0.02, Nt::Int64 = 200, dir = 1, redis=1, tol::Float64 = 1e-6, tau=1, o::Int64=0)`: the string method from the string phi

    input :
    phi : the original string.

"""
function string(yA::Array{Float64,2},yB::Array{Float64,2};M::Int64=20,dt::Float64=0.02,Nt::Int64=200,dir=1,redis=1,tol::Float64=1e-6,tau=1,o::Int64=0)

    dir in [0 1] ? nothing :  error(" dir must equals 0 or 1")
    redis in [0 1] ? nothing :  error(" redis must equals 0 or 1")
     tau in [0 1] ? nothing :  error(" tau must equals 0,1,2")

# the orginal string by linear interpolation
    phi = zeros(length(yA[:,1]),length(yA[1,:]),M+2)
    
    for k = 1:M+2
        phi[:,:,k] = (M+2-k)/(M+1).*yA .+ (k-1)/(M+1).*yB
    end
    
    out = string(phi;dt=dt,Nt=Nt,dir=dir,redis=redis,tol=tol,tau=tau,o=o)
    return out
end



function string(phi::Array{Float64,3}; dt::Float64 = 0.02, Nt::Int64 = 200, dir = 1, redis=1, tol::Float64 = 1e-6, tau=1, o::Int64=0)
    dir in [0 1] ? nothing :  error(" dir must equals 0 or 1")
    redis in [0 1] ? nothing :  error(" redis must equals 0 or 1")
    tau in [0 1] ? nothing :  error(" tau must equals 0,1,2")

    M = length(phi[1,1,:]) - 2
    d = length(phi[:,1,1])
    N = length(phi[1,:,1])

    println("dimension :", d, ",  number of atomes in a configuration : ", N, ",  images : ", M, ",  time step : ", dt)

    dir == 0 ? println(" evolution direction : -∇V ") : nothing
    dir == 1 ? println(" evolution direction : (-∇V)^⟂ ") : nothing

    redis == 0 ? println(" redistribution or not : no ") : nothing
    redis == 1 ? println(" redistribution or not : yes ") : nothing

    tau == 0 ? println("choose central difference to get tangent") : nothing
    tau == 1 ? println("choose upwind scheme to get tangent") : nothing

    Ap = zeros(d, N, M+2, Nt)
    Err = zeros(Nt,1)


    k = 1    #k-th string
    Ap[:,:,:,k] = phi
    println(" time for stringforces in each step : ")

    @time begin f = stringforces(phi) end

    tan = tangent(phi,tau)
    f0 = stringf2f0(f,tan)

    println("---------------------------------------------")
    println(" time step |  ||f0||_∞  |  position_of_∞  ")
    println("---------------------------------------------")

    F0min,indmin = f2F(f0)
    println(" ", k, "  ", F0min, "  ", indmin)

    while k<Nt
        dir == 0 ? phi += dt*f : nothing
        dir == 1 ? phi += dt*f0 : nothing

        if redis == 1
            _,s = distance(phi)
            phi = Redistribution(phi,3,s)
        else nothing
        end

#         if redis == 1
#             l,s = distance(phi)
#             lmax,_ = findmax(l)
#             lmin,_ = findmin(l)
#             ll = (lmax + lmin)/2
#             if (lmax-lmin)/ll > 0.2
#                 phi = Redistribution(phi,3,s)
#                 println("redistribution")
#             else nothing
#             end

#         else nothing
#         end

        k+=1 #k-th string
        Ap[:,:,:,k] = phi
        f = stringforces(phi)
        tan = tangent(phi,tau)
        f0 = stringf2f0(f,tan)

        F0min,indmin = f2F(f0)
        println(" ", k, "  ", F0min, "  ", indmin)
        Err[k] = copy(F0min)

        F0min<tol ? break : nothing

#         if k>Nt/5
#             if Err[k]>Err[k-1]>Err[k-2]>Err[k-3]
#                 error("the residual starts to increase ")
#                 p = copy(Ap[:,:,:,k-3])
#                 break
#             end
#         end
        
        if k>Nt/2
            if abs(Err[k]-Err[k-1])/abs(Err[k-1]) < 1e-4 && abs(Err[k-1]-Err[k-2])/abs(Err[k-2]) < 1e-4 && abs(Err[k-2]-Err[k-3])/abs(Err[k-3]) < 1e-4 
                println("the residual does not reach the tolerance")
                p = copy(Ap[:,:,:,k])
                break
            end
        end
    end

     p = copy(Ap[:,:,:,k])
    if o==1
        return Ap
    end
    if o==0
        return p
    end
end
