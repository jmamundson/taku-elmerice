BEGIN{
    aB = 4000.0   # height of gaussian peak [m] 
    bB = -30000.0 # location of peak of gaussian [m]
    cB = 20000.0  # width of gaussian [m]

    aS = 300.0    # height of sill [m]
    bS = 40000.0  # location of sill [m]
    cS = 3000.0   # width of sill [m]

    zmin = -500.0 # min fjord elevation
    slope = -1.0/40.0
    z0 = 1200.0

    maxz = 1000.0
    minz = 0.0
}
{
    x=$3
    z=$4
    bedrock = aB*exp(-((x-bB)^2.0)/(2.0*cB*cB))
    sill = aS*exp(-((x-bS)^2.0)/(2.0*cS*cS))
    bedelevation = bedrock + sill + zmin
    surfaceelevation = z0 + slope*x
    relz = z/(maxz - minz)
    meshnodeelevation = bedelevation + relz*(surfaceelevation - bedelevation)
    print $1, $2, $3, meshnodeelevation, $5
    #print x, bedelevation, surfaceelevation, relz, meshnodeelevation
    #print x,meshnodeelevation
}
