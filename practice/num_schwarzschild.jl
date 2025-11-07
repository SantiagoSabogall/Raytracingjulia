using Dierckx, DelimitedFiles

mutable struct BlackHole
    data_N::DataFrame
    N::Spline1D
    data_dNdr::DataFrame
    dNdr::Spline1D
    a::Float64
    EH::Float64
    ISCOco::Float64
    ISCOcounter::Float64

    function BlackHole()
        data_N = readdlm("ruta N")
        N = Spline1D(data_N[:,1], data_N[:,2])
        data_dNdr = readdlm("ruta der")
        dNdr = Spline1D(data_dNdr[:,1], data_dNdr[:,2])
        a = 0
        EH = 2
        ISCOco = 6
        ISCOcounter = 6

        new(data_N, N, data_dNdr, dNdr, a, EH, ISCOco, ISCOcounter)
    end

end

function metric(b::BlackHole, x::Vector{<:Real})

        g_tt = - b.N(x[2])
        g_rr = 1/b.N(x[2])
        g_thth = x[1]^2
        g_phph = (x[1]*sin(x[3]))^2
        g_tph = 0.
        
        return [g_tt, g_rr, g_thth, g_phph, g_tph]

end

function inverse_metric(b::BlackHole, x::Vector{<:Real})

        gtt = - 1/b.N(x[2])
        grr = b.N(x[2])
        gthth = 1/x[2]^2
        gphph = 1/(x[2]*sin(x[3]))^2
        gtph = 0.
        
        return [gtt, grr, gthth, gphph, gtph]

end

function dr_inverse_metric(b::BlackHole, x::Vector{<:Real})

        drgtt =  b.dNdr(x[2])/(b.N(x[2])^2)
        drgrr = b.dNdr(x[2])
        drgthth = -2/x[2]^3
        drgphph = -2/(x[2]^3*sin(x[3])^2)
        drgtph = 0.
        
        return [drgtt, drgrr, drgthth, drgphph, drgtph]

end

function geodesics(b::BlackHole, q::Vector{<:Real}, lmbda::Float64)

        gtt, grr, gthth, gphph, gtph = b.inverse_metric(q[1:2])
        drgtt, drgrr, drgthth, drgphph, drgtph = self.dr_inverse_metric(q[1:5])
        
        # Geodesics differential equations 
        dtdlmbda = gtt*q[5]
        drdlmbda = grr*q[6]
        dthdlmbda = gthth*q[7]
        dphidlmbda = gphph*q[8]
        
        dk_tdlmbda = 0.
        dk_rdlmbda = - (drgtt*q[5]^2)/2 - (drgrr*q[6]^2)/2 \
                     - (drgthth*q[7]^2)/2 - (drgphph*q[8]^2)/2
        dk_thdlmbda = (cos(q[3])/sin(q[3])^3)*(q[8]/q[2])^2
        dk_phidlmbda = 0.
        
        return [dtdlmbda, drdlmbda, dthdlmbda, dphidlmbda, 
                dk_tdlmbda, dk_rdlmbda, dk_thdlmbda, dk_phidlmbda]

end

if abspath(PROGRAM_FILE) == @__FILE__
    println()
    println("THIS IS A MODULE DEFINING ONLY A PART OF THE COMPLETE CODE.")
    println("YOU NEED TO RUN THE main.jl FILE TO GENERATE THE IMAGE")
    println()
end