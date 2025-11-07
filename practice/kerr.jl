
mutable struct BlackHole
a::Float64
EH::Float64
Z1::Float64
Z2::Float64
ISCOco::Float64
ISCOcounter::Float64

  function BlackHole(a::Float64)
    EH = 1 + sqrt(1 - a^2)
    Z1 = 1 + (1 - a^2)^(1/3)*((1 + a)^(1/3) + (1 - a)^(1/3))
    Z2 = sqrt(3*a^2 + Z1^2)
    ISCOco = 3 + Z2 - sqrt((3 - Z1)*(3 + Z1 + 2*Z2))
    ISCOcounter = 3 + Z2 + sqrt((3 - Z1)*(3 + Z1 + 2*Z2))

    new(a, EH, Z1, Z2, ISCOco, ISCOcounter)
  end

end


function Omega(b::BlackHole, r::Real, corotating::Bool=true)
  if corotating
      return 1/(r^(3/2) + b.a)
  else
      return -1/(r^(3/2) - b.a)

  end

end

function metric(b::BlackHole, x::Vector{<:Real})

        r2 = x[2]*x[2]
        a2 = b.a*b.a
        sin_theta2 = sin(x[3])^2
        Delta = r2 - 2*x[2] + a2
        Sigma = r2 + a2*cos(x[3])^2

        g_tt = -(1 - 2*x[2]/Sigma)
        g_rr = Sigma/Delta
        g_thth = Sigma
        g_phph = (r2 + a2 + 2*a2*x[2]*sin_theta2/Sigma)*sin_theta2
        g_tph = -2*b.a*x[2]*sin_theta2/Sigma

        return [g_tt, g_rr, g_thth, g_phph, g_tph]

end

function inverse_metric(b::BlackHole, x::Vector{<:Real})

        r2 = x[2]*x[2]
        a2 = b.a*b.a
        sin_theta2 = sin(x[3])^2
        Delta = r2 - 2*x[2] + a2
        Sigma = r2 + a2*cos(x[3])^2
        A = (r2 + a2)^2 - Delta*a2*sin_theta2

        gtt = - A/(Delta*Sigma)
        grr = Delta/Sigma
        gthth = 1/Sigma
        gphph = (Delta - a2*sin_theta2)/(Delta*Sigma*sin_theta2)
        gtph = - 2*b.a*x[2]/(Delta*Sigma)

        return [gtt, grr, gthth, gphph, gtph]
end

function geodesics(b::BlackHole, q::Vector{<:Real}, lmbda::Float64)

        r2 = q[2]*q[2]
        a2 = b.a*b.a
        sin_th = sin(q[3])
        cos_th = cos(q[3])
        sin_th2 = sin_th*sin_th
        cos_th2 = cos_th*cos_th
        Sigma = r2 + a2*cos_th2
        Sigma2 = Sigma*Sigma
        Delta = r2 - 2*q[2] + a2

        W = -q[5]*(r2 + a2) - b.a*q[8]
        partXi = r2 + (q[8] + b.a*q[5])^2 + a2*(1 + q[5]*q[5])*cos_th2 + q[8]*q[8]*cos_th2/sin_th2
        Xi = W^2 - Delta*partXi

        dXidE = 2*W*(r2 + a2) + 2*b.a*Delta*(q[8] + b.a*q[5]*sin_th2)
        dXidL = -2*b.a*W - 2*b.a*q[5]*Delta - 2*q[8]*Delta/sin_th2

        dXidr = -4*q[2]*q[5]*W - 2*(q[2] - 1)*partXi - 2*q[2]*Delta

        dAdr = (q[2] - 1)/Sigma - (q[2]*Delta)/Sigma2
        dBdr = -q[2]/Sigma^2
        dCdr = dXidr/(2*Delta*Sigma) - (Xi*(q[2]-1))/(Sigma*Delta*Delta) - q[2]*Xi/(Delta*Sigma2)

        auxth = a2*cos_th*sin_th

        dAdth = Delta*auxth/Sigma2
        dBdth = auxth/Sigma2
        dCdth = ((1+q[8]^2)*auxth + q[8]*q[8] * cos_th/(sin_th2*sin_th) )/Sigma + (Xi/(Delta*Sigma2))*auxth


        dtdlmbda = dXidE/(2*Delta*Sigma)
        drdlmbda = (Delta/Sigma)*q[6]
        dthdlmbda = q[7]/Sigma
        dphidlmbda = - dXidL/(2*Delta*Sigma)

        dk_tdlmbda = 0.
        dk_rdlmbda = -dAdr*q[6]*q[6] - dBdr*q[7]*q[7] + dCdr
        dk_thdlmbda = -dAdth*q[6]*q[6] - dBdth*q[7]*q[7] + dCdth
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