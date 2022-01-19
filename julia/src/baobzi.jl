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

function eval(func::Ptr{baobzi_struct}, x)
    return ccall(
        (:baobzi_eval, :libbaobzi),
        Cdouble,
        (Ptr{Cvoid}, Ptr{Cdouble},),
        func, x
    )
end

function free(func::Ptr{baobzi_struct})
    ccall(
        (:baobzi_free, :libbaobzi),
        Ptr{Cvoid},
        (Ptr{Cvoid},),
        func
    )
end

function save(func::Ptr{baobzi_struct}, filename::String)
    ccall(
        (:baobzi_save, :libbaobzi),
        Ptr{Cvoid},
        (Ptr{Cvoid}, Cstring),
        func, filename
    )
end

function restore(filename::String)
    output_ptr = ccall(
        (:baobzi_restore, :libbaobzi),
        Ptr{baobzi_struct},
        (Cstring,),
        filename
    )

    return output_ptr
end

end
