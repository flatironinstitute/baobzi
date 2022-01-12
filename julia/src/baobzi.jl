module baobzi

mutable struct baobzi_struct
end

function init(fin, dim, order, center, half_length, tol::Float64)
    c_fin = @cfunction($fin, Cdouble, (Ptr{Cdouble},))
    output_ptr = ccall(
        (:baobzi_init, :libbaobzi),
        Ptr{baobzi_struct},
        (Ptr{Cvoid}, Cushort, Cushort, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble,),
        c_fin, dim, order, center, half_length, tol
    )

    return output_ptr
end

function eval(func, x)
    return ccall(
        (:baobzi_eval, :libbaobzi),
        Cdouble,
        (Ptr{Cvoid}, Ptr{Cdouble},),
        func, x
    )
end

end
