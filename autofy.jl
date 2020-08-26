

using LinearAlgebra
using SymPy

using DynamicalSystems
using DelimitedFiles

using OrdinaryDiffEq
using Printf

function csetup(Ndim,Npar)
    #simplest: done.
    string=@sprintf "NDIM=   %d\nNPAR = %d\nIPS =   1, IRS =   0, ILP =   1\n#ICP is changed in the auto file\nICP =  [3,2,5,6,7]\nNTST=  55, NCOL=   4, IAD =   3, ISP =   2, ISW = 1, IPLT= 0, NBC= 0, NINT= 0\nNMX=   250, NPR=  25, MXBF=  0, IID =   2, ITMX= 10, ITNW= 10, NWTN= 3, JAC= 0\nEPSL= 1e-08, EPSU = 1e-08, EPSS =1e-06\nDS  =   5e-03, DSMIN= 1e-03 , DSMAX= 1e-02, IADS=   1\nTHL =  {11: 0.0}, THU =  {}\nUZSTOP = {1: [0.01,0.2], 2:[0.01,0.99]}\n" Ndim Npar

    open("c.xxx", "w") do f
        write(f, string)
    end
end

function fsetup(Fstr,Ustr,pars,Ustart,pstart)
    #simplest: done.
    funcstring=@sprintf "SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)\n\n              IMPLICIT NONE\n              INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC\n             DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)\n              DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)\n              DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)\n\n"#%s\nEND SUBROUTINE FUNC\n" Fstr
    open("xxx.f90", "w") do f
        write(f, funcstring)
    end
    #declarations
    for i in 1:length(Ustr)
        funcstring=@sprintf "DOUBLE PRECISION %s\n" Ustr[i]
        open("xxx.f90", "a+") do f
            write(f, funcstring)
        end
    end
    for i in 1:length(pars)
        funcstring=@sprintf "DOUBLE PRECISION %s\n" pars[i]
        open("xxx.f90", "a+") do f
            write(f, funcstring)
        end
    end
    for i in 1:length(Ustr)
        funcstring=@sprintf "%s=U(%d)\n" Ustr[i] i
        open("xxx.f90", "a+") do f
            write(f, funcstring)
        end
    end
    for i in 1:length(Fstr)
        funcstring=@sprintf "F(%d)=%s\n" i Fstr[i]
        open("xxx.f90", "a+") do f
            write(f, funcstring)
        end
    end
    funcstring="\nEND SUBROUTINE FUNC\n"
    open("xxx.f90", "a+") do f
        write(f, funcstring)
    end

    stpstring=@sprintf "\n\nSUBROUTINE STPNT(NDIM,U,PAR,T)\n\nIMPLICIT NONE\nINTEGER, INTENT(IN) :: NDIM\nDOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)\nDOUBLE PRECISION, INTENT(IN) :: T\n\n"#%s\nEND SUBROUTINE STPNT\n" Ustr

    open("xxx.f90", "a+") do f
        write(f, stpstring)
    end
    for i in 1:length(Ustr)
        stpstring=@sprintf "\nU(%d)=%f\n" i Ustart[i]
        open("xxx.f90", "a+") do f
            write(f, stpstring)
        end
    end

    for i in 1:length(pars)
        stpstring=@sprintf "\nPAR(%d)=%f\n" i pstart[i]
        open("xxx.f90", "a+") do f
            write(f, stpstring)
        end
    end

    stpstring="\nEND SUBROUTINE STPNT\n"
    open("xxx.f90", "a+") do f
        write(f, stpstring)
    end

    bcndstring="SUBROUTINE BCND\nEND SUBROUTINE BCND\nSUBROUTINE ICND\nEND SUBROUTINE ICND\nSUBROUTINE FOPT\nEND SUBROUTINE FOPT\n\nSUBROUTINE PVLS(NDIM,U,PAR)\nIMPLICIT NONE\nINTEGER, INTENT(IN) :: NDIM\nDOUBLE PRECISION, INTENT(IN) :: U(NDIM)\nDOUBLE PRECISION, INTENT(INOUT) :: PAR(*)\n\nDOUBLE PRECISION, EXTERNAL :: GETP\n\n\n\nEND SUBROUTINE PVLS\n"
    open("xxx.f90", "a+") do f
        write(f, bcndstring)
    end
end

function autosetup(ICP)
    #for i in 1:length(ICP)
        autputstrings=@sprintf "start1=run(e='xxx',c='xxx',ICP=[%d])" ICP
        open("xxx.auto", "w") do f
            write(f, autputstrings)
        end
    #end
end

function parseinauto(type,ICPs)

        autputstrings=@sprintf "start2=run(start1('%s'),c='xxx',ICP=[%d])" type ICPs
        open("xxx.auto", "a+") do f
            write(f, autputstrings)
        end

        #For plotting.. should i use all ICPs? Why not.
        #@sprintf "%s" "%f" works amazingly
        dummy=""
        #dummy=repeat('%f ', length(ICPs))

#WIP: the number of "%f" above does not equal the start2['K'] number yet. Nasil?
        #outcomm=@sprintf "filehandle.write('%s\n' % (start2['K'][i],start2['r'][i],start2['STA'][i]))" dummy
        autputstrings=""
        #autputstrings=@sprintf "datname='dat/%s.dat'\nif(len(start2)>0):\n\twith open(datname, 'w') as filehandle:\n\t
        #        for i in range(1,len(start2['K'])):\n\tfilehandle.write('%f %f %f\n' % (start2['K'][i],start2['r'][i],start2['STA'][i]))" type
        open("xxx.auto", "a+") do f
            write(f, autputstrings)
        end
end



#not inuse
function autofy(Fs,xs,Parsym,Parsstart)
    #W I P
    for i in 1:length(xs)
        xstrings[i]=@sprintf "U(%d)" i
    end
    for i in 1:length(Fs)
        Fstrings[i]=@sprintf "F(%d)" i
    end

    Fssub=symbols("Fssub")
    Fssub=Fs
    for i in 1:length(Parsym)
        Fssub=Fssub.subs(Parsym[i],Parsstart[i])
    end
    fxpts=sympy.solve(Fssub)

    for i in 1:length(fxpts)

    end

end

function autorun()
    #there is no identifiable renaming of the auto files and so on...

    @printf("Initial run...\n")
    run(`auto xxx.auto`)

    @printf("Continuation of...\n")
    continuationof=readline()

    parseinauto(continuationof,1)

end

function symjlsetup(xs,contents,vecname)
    for i in 1:length(xs)-1
        open("dummy.jl", "a+") do f
            #write(f, "function Jac(U,PAR)\n")
            outi=@sprintf "%s," xs[i]
            write(f, outi)
        end
    end
    open("dummy.jl", "a+") do f
        #write(f, "function Jac(U,PAR)\n")
        outi=@sprintf """%s=symbols( " """ xs[length(xs)]
        write(f, outi)
    end

    for i in 1:length(xs)-1
        open("dummy.jl", "a+") do f
            #write(f, "function Jac(U,PAR)\n")
            outi=@sprintf "%s," xs[i]
            write(f, outi)
        end
    end
    open("dummy.jl", "a+") do f
        outi=@sprintf """%s")\n""" xs[length(xs)]
        write(f, outi)
    end
    if length(contents)>1
        #vecname="Fsym"
        for i in 1:length(contents)
            open("dummy.jl", "a+") do f
                #write(f, "function Jac(U,PAR)\n")
                outi=@sprintf "%s=%s\n" xs[i] contents[i]
                write(f, outi)
            end
        end
    end

    open("dummy.jl", "a+") do f
        outi=@sprintf "%s=[%s " vecname xs[1]
        write(f, outi)
    end
    for i in 2:length(xs)-1
        open("dummy.jl", "a+") do f
            outi=@sprintf "%s " xs[i]
            write(f, outi)
        end
    end
    open("dummy.jl", "a+") do f
        outi=@sprintf "%s]\n" xs[length(xs)]
        write(f, outi)
    end
end

function jljacsetup(Fs,xs,pars)#pstart
    names=Vector{String}()
    for i in 1:length(Fs)
        push!(names,@sprintf "F%d" i)
    end
    symjlsetup(xs, "","xsym")
    symjlsetup(pars,"","par")
    symjlsetup(names,Fs,"Fsym")

    #print(names)
    #print(Fs)
    include("dummy.jl")
    Jac=Array{String}(undef,length(names),length(xs))

    for i in 1:length(Fsym)
        for j in 1:length(xsym)
            #needs to become symbolic first...
            if !isnothing(sympy.rcode(sympy.diff(Fsym[i],xsym[j])))
                Jac[i,j]=sympy.rcode(sympy.diff(Fsym[i],xsym[j]))
                #print(i)
            end
        end
    end

    #next: parse into xxx.jl, with parameters eingesetzt
    #begin
    open("xxx.jl", "a+") do f
        write(f, "function Jac(U,PAR, t)\n")
    end

    for i in 1:length(xs)
        autputstrings=@sprintf "\t%s=U[%d] \n" xs[i] i
        open("xxx.jl", "a+") do f
            write(f, autputstrings)
        end
    end

    for i in 1:length(pars)
        autputstrings=@sprintf "\t%s=PAR[%d] \n" pars[i] i
        open("xxx.jl", "a+") do f
            write(f, autputstrings)
        end
    end

    open("xxx.jl", "a+") do f
        write(f, "\nJac=[")
    end

    for i in 1:length(Fsym)
        for j in 1:length(xsym)

            open("xxx.jl", "a+") do f
                write(f,@sprintf "%s\t" Jac[i,j])
            end
        end
        open("xxx.jl", "a+") do f
            write(f,";\n")
        end
    end
    #end
    open("xxx.jl", "a+") do f
        write(f, "]\nreturn Jac\nend\n")
    end
end

function jlsetup(Fs,xs,pars)
    open("xxx.jl", "w") do f
        write(f, "\n")
    end
    open("xxx.jl", "a+") do f
        write(f, "function Funcs(U,PAR, t)\nF=zeros(length(U));\n")
    end

    for i in 1:length(xs)
        autputstrings=@sprintf "\t %s=U[%d] \n" xs[i] i
        open("xxx.jl", "a+") do f
            write(f, autputstrings)
        end
    end

    for i in 1:length(pars)
        autputstrings=@sprintf "\t%s=PAR[%d] \n" pars[i] i
        open("xxx.jl", "a+") do f
            write(f, autputstrings)
        end
    end

    for i in 1:length(Fs)
        autputstrings=@sprintf "\tF[%d]=%s \n" i Fs[i]
        open("xxx.jl", "a+") do f
            write(f, autputstrings)
        end
    end

    open("xxx.jl", "a+") do f
        write(f, "return SVector{length(F)}(F)\nend\n\n\n")
    end



end

function autosiz()
    run(`rm -f dummy.jl`)
    #core::
        xs=["xA" "mu"]
        pars=["xB" "K" "r"]
        Fs=["xA*(xA*xA-mu)+K*(r*xA+(1-r)*xB)" "xB*(xB*xB-mu)+K*(r*xB+(1-r)*xA)"]

    #@printf("les variables s'il vous plait... quitter avec q\n")
    #xs=Vector{String}()
    #while(readline()!="q")
    #    push!(xs,readline())
    #    @printf("OK\n")
    #end

    #@printf("les parametres s'il vous plait... quitter avec q\n")
    #pars=Vector{String}()
    #while(readline()!="q")
    #    push!(pars,readline())
    #    @printf("OK\n")
    #end


    #@printf("les functions s'il vous plait... quitter avec q\n")
    #Fs=Vector{String}()
    #while(readline()!="q")
    #    push!(Fs,readline())
    #    @printf("OK\n")
    #end


    #now we would need to run this a bit to retrieve Ustart from pstart...


    csetup(length(xs),length(pars))
    run(`cat c.xxx`)
    Ustart=rand(length(xs))

    pstart=rand(length(pars))

    #pstart[3]=0.5

    #if we trust the solvability:
    iwantthat=0
    if iwantthat==1
        Fsym=Vector{String}(undef,length(Fs))
        for i in 1:length(Fs)
            Fsym[i]=Fs[i]
            for j in 1:length(pars)
                #wronk
                Fsym[i]=@sprintf "%s.subs(%s,%f)" Fsym[i] pars[j] pstart[j]
            end

        end
        #print(Fsym)
        sol=(sympy.solve(Symbol.(Fsym)))
    end
        #print(values(sol[2]))
    #in other cases get sol as a return value from simulations.
    jlsetup(Fs,xs,pars)
    jljacsetup(Fs,xs,pars)
    include("xxx.jl")

    dsl=ContinuousDynamicalSystem(Funcs,Ustart,pstart,Jac)
    traj=trajectory(dsl,900,dt=0.4)

    fsetup(Fs,xs,pars,traj[end],pstart)

    autosetup(3)

    run(`ls`)
    autorun()

end
