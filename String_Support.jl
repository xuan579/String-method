

# just need to add energy and forces to complete a string method
"""
`stringenergy(phi::Array{Float64,3})`: the energy of all images on the string `phi`.
"""
function stringenergy(phi::Array{Float64,3})
    N = length(phi[1,1,:]) # the number of images
    E = zeros(N)
    for i = 1:N
        # set_positions!(at, phi[:,:,i])
        E[i] = energy(phi[:,:,i])
    end
    return E
end

"""
`stringforces(phi::Array{Float64,3})`: the forces of all images on the string `phi`.
"""
function stringforces(phi::Array{Float64,3})
    M = length(phi[1,1,:]) - 2 # the number of images (not including yA and yB)
    d = length(phi[:,1,1])
    N = length(phi[1,:,1])

    f = zero(phi)   # -\nabla V

    for i = 2:M+1
        f[:,:,i] = forces(phi[:,:,i])
    end

    return f
end


"""
    `stringf2f0`(f::Array{Float64,3}, tan::Array{Float64,3})` :  get the compoment of the forces `f` normal to string phi.
"""
function stringf2f0(f::Array{Float64,3}, tan::Array{Float64,3})
    M = length(f[1,1,:]) - 2 # the number of images (not including yA and yB)

    f0 = zero(f)   # -nabla V 关于string的法向分量

    for i = 1:M+2
    f0[:,:,i] = f[:,:,i] - sum(f[:,:,i].*tan[:,:,i])*tan[:,:,i]   #计算法向分量
    end

    return f0
end

"""
`tangent(phi::Array{Float64,3},tau::Int64)` : obtain the tangent by finite difference
"""
function tangent(phi::Array{Float64,3},tau::Int64)
    M = size(phi,3) - 2

    tan = zero(phi)
    tan[:,:,1] = (phi[:,:,2] - phi[:,:,1])/norm(phi[:,:,2] - phi[:,:,1])
    tan[:,:,M+2] = (phi[:,:,M+2] - phi[:,:,M+1])/norm(phi[:,:,M+2] - phi[:,:,M+1])
    
    if tau == 0
        for i = 2:M+1
            tan[:,:,i] = (phi[:,:,i+1] - phi[:,:,i-1])/norm(phi[:,:,i+1] - phi[:,:,i-1])
        end
    end
    if tau == 1
        E = stringenergy(phi)

        for i = 2:M+1
            if E[i-1]<E[i]<E[i+1]
                tan[:,:,i] = (phi[:,:,i+1] - phi[:,:,i])/norm(phi[:,:,i+1] - phi[:,:,i])
            else
                if E[i-1]>E[i]>E[i+1]
                    tan[:,:,i] = (phi[:,:,i] - phi[:,:,i-1])/norm(phi[:,:,i] - phi[:,:,i-1])
                else
                    tan[:,:,i] = (phi[:,:,i+1] - phi[:,:,i-1])/norm(phi[:,:,i+1] - phi[:,:,i-1])
                end
            end
        end
    end
    return tan
end


"""
`distance(phi::Array{Float64,3})` :
    input :
        phi : a string.
    output :
        l : 相邻点的距离
        s : 第1个点到第n个点的总距离
"""
function distance(phi::Array{Float64,3})
    M = length(phi[1,1,:]) # the number of images(including two endpoints)
#     N = length(phi[1,:,1]) # the number of atoms in one configuration
#     d = length(phi[:,1,1]) # the dims of one atom

    #计算第i个点右侧的弧长，存于l,前i个点的弧长，存于s
    l = zeros(M-1)
    s = zeros(M)
    for i = 1:M-1
        l[i] = norm(phi[:,:,i] - phi[:,:,i+1])
        s[i+1] = s[i] + l[i]
    end
    return l,s
end

"""
`Redistribution(phi::Array{Float64,3},kk::Int64,s)` : redistribution by kk-th polynomial, when kk = 3, it means cubic spline.
"""
function Redistribution(phi::Array{Float64,3},kk::Int64,s)
    Phi = zero(phi)
    M = length(phi[1,1,:]) # the number of images(including two endpoints)
    N = length(phi[1,:,1]) # the number of atoms in one configuration
    d = length(phi[:,1,1]) # the dims of one atom

#     #计算第i个点右侧的弧长，存于l,前i个点的弧长，存于s
#     l = zeros(M-1)
#     s = zeros(M)
#     for i = 1:M-1
#         l[i] = norm(phi[:,:,i] - phi[:,:,i+1])
#         s[i+1] = s[i] + l[i]
#     end

        alpha = s ./ s[M]    #弧长参数
        for i = 1:d
            for j = 1:N
                spl = Spline1D(alpha,phi[i,j,:],k=kk)
                Alpha = collect(0:1/(M-1):1)
                Phi[i,j,:] = spl(Alpha)
            end
        end
    return Phi
end


function f2F(f::Array{Float64,3})
    M = size(f,3) - 2
    F0 = zeros(M+2)
    for j = 1:M+2
        F0[j],_ = findmax(abs.(f[:,:,j]))
    end
    F0min,indmin = findmax(F0)
    return F0min,indmin
end

"""
 `getEf(pp::Array{Float64,3},dpp::Array{Float64,3})` : numerical calcute the line integral S(φ)= ∫_0^1 |∇ E(ϕ(s))|⋅|̇ϕ'(s)| ds.
    Here, pp is discrete ϕ and dpp is discrete ϕ'.
"""
function getEf(pp::Array{Float64,3},dpp::Array{Float64,3})
    N = size(pp,3)
    F1 = zeros(N)
    F2 = zeros(N)
    
    f1 = stringforces(pp)
    
    for k = 1:N
        F1[k] = norm(f1[:,:,k])
        F2[k] = norm(dpp[:,:,k])
    end
    F = F1.*F2
    
    E = ( 0.5*(F[1]+F[end])+sum(F[2:end-1]) )/(N-1)
end

"""
`getEb(pp::Array{Float64,3})` : The energy barierr of the string pp.
"""
function getEb(pp::Array{Float64,3})
    N = size(pp,3)
    
    E = stringenergy(pp)
    Emax,_ = findmax(E)
    Emin,_ = findmin(E)
    
    Eb = Emax-Emin
end

"""
`getdph(ph::Array{Float64,3})` : get the upwind scheme first order difference.
"""
function getdph(ph::Array{Float64,3})
    Eh =stringenergy(ph)
    _,sad = findmax(Eh)
    M = length(ph[1,1,:]) - 2
    h = 1/(M+1)
    
    dph = zero(ph)    
    for k = 1:sad-1
        dph[:,:,k] = (ph[:,:,k+1]-ph[:,:,k])*(M+1)
    end
    dph[:,:,sad] = (ph[:,:,sad+1]-ph[:,:,sad-1])*(M+1)/2
    for k = sad+1:M+2
        dph[:,:,k] = (ph[:,:,k]-ph[:,:,k-1])*(M+1)
    end
    return dph
end

"""
`addimage(p::Array{Float64,3},n::Int64)` : add n images between every points on curve p by linear interpolation.
"""
function addimage(p::Array{Float64,3},n::Int64)

    M = length(p[1,1,:]) - 2
    d = size(p,1)
    N = size(p,2)
    
    pp = zeros(d,N,1+(n+1)*(M+1))
    
    pp[:,:,1] = p[:,:,1]
    for k = 1:M+1, t = 1:n+1
        pp[:,:,1+(k-1)*(n+1)+t] = (1- t/(n+1))*p[:,:,k] + (t/(n+1))*p[:,:,k+1]
    end

    return pp
end
